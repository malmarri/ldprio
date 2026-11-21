mod ld;
mod plink;
mod report;
mod select;

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

use clap::{Parser, ValueEnum};

use crate::ld::{build_ld_graph_window, find_backfill_candidates, Window};
use crate::plink::{preview_list, resolve_samples, scan_bed, write_bed_subset, Bim, Fam, QcFilters};
use crate::report::{fmt_int, write_report, ReportData};
use crate::select::{backfill_select, greedy_select, PriorityPool, Weighting};

const POOL_WARN_FRACTION: f64 = 0.70;
const MIN_LD_SAMPLES: usize = 30;

#[derive(Copy, Clone, Debug, ValueEnum)]
enum WeightingArg {
    Adaptive,
    Inverse,
    Fair,
}

impl From<WeightingArg> for Weighting {
    fn from(w: WeightingArg) -> Self {
        match w {
            WeightingArg::Adaptive => Weighting::Adaptive,
            WeightingArg::Inverse => Weighting::Inverse,
            WeightingArg::Fair => Weighting::Fair,
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "ldprio",
    version,
    about = "LD pruning that preserves SNPs covered by nominated low-coverage samples (native Rust, no PLINK required).",
    after_help = "Sample lists: one per line, 'IID' or 'FID<tab>IID'.\nOnly autosomes (chromosomes 1-22) are pruned and written out."
)]
struct Args {
    /// Input PLINK 1 prefix (.bed/.bim/.fam)
    #[arg(long)]
    bfile: String,

    /// Output prefix
    #[arg(long)]
    out: String,

    /// Samples whose covered SNPs should be favoured
    #[arg(long = "priority-samples")]
    priority_samples: PathBuf,

    /// Samples LD is estimated from (high-quality diploids; at least 30 recommended)
    #[arg(long = "ld-samples")]
    ld_samples: PathBuf,

    /// Window in variants (PLINK-style; counted among surviving SNPs)
    #[arg(long, default_value_t = 200)]
    window: usize,

    /// Window as physical distance in kb, from .bim positions (replaces --window;
    /// LD-independence is reached in a single pass, no cleanup needed)
    #[arg(long = "window-kb", conflicts_with = "window")]
    window_kb: Option<f64>,

    /// Unused (kept for CLI compatibility with v1.x)
    #[arg(long, default_value_t = 25)]
    step: usize,

    /// r^2 threshold
    #[arg(long, default_value_t = 0.4)]
    r2: f64,

    /// Minimum minor allele frequency in --ld-samples (e.g. 0.01)
    #[arg(long)]
    maf: Option<f64>,

    /// Maximum missing genotype rate in --ld-samples (e.g. 0.05)
    #[arg(long)]
    geno: Option<f64>,

    /// Minimum one-sided exact excess-heterozygosity p-value in --ld-samples (e.g. 1e-6, or 'auto')
    #[arg(long = "excess-het")]
    excess_het: Option<String>,

    /// How blocks contested between nominated samples are settled
    #[arg(long, value_enum, default_value_t = WeightingArg::Adaptive)]
    weighting: WeightingArg,

    /// Among otherwise equal SNPs, keep transversions over transitions (C<->T, G<->A)
    #[arg(long = "transv", alias = "prefer-transversions")]
    prefer_transversions: bool,

    /// Also write a self-contained HTML report (<out>.report.html)
    #[arg(long)]
    report: bool,

    /// Also write the pruned PLINK fileset (<out>.{bed,bim,fam})
    #[arg(long = "make-bed")]
    make_bed: bool,

    /// Cap on cleanup passes
    #[arg(long = "max-cleanup", default_value_t = 10)]
    max_cleanup: usize,

    /// Number of threads (0 = all available CPU cores)
    #[arg(long, short = 't', default_value_t = 0)]
    threads: usize,
}

fn emit(messages: &mut Vec<String>, text: &str) {
    eprintln!("{}", text);
    let collapsed = text.split_whitespace().collect::<Vec<_>>().join(" ");
    messages.push(collapsed);
}

fn main() {
    let args = Args::parse();
    let start = Instant::now();

    if let Err(err) = run(args, start) {
        eprintln!("{}", err);
        process::exit(1);
    }
}

fn run(args: Args, start: Instant) -> Result<(), String> {
    if !(args.r2 > 0.0 && args.r2 <= 1.0) {
        return Err(format!("ERROR: --r2 must be in (0, 1], got {}", args.r2));
    }
    if let Some(maf) = args.maf {
        if !(0.0..=0.5).contains(&maf) {
            return Err(format!("ERROR: --maf must be in [0, 0.5], got {}", maf));
        }
    }
    if let Some(geno) = args.geno {
        if !(0.0..=1.0).contains(&geno) {
            return Err(format!("ERROR: --geno must be in [0, 1], got {}", geno));
        }
    }
    let excess_het_mode: Option<Result<f64, ()>> = match args.excess_het.as_deref() {
        None => None,
        Some(s) if s.eq_ignore_ascii_case("auto") => Some(Err(())), // Err(()) = auto Bonferroni
        Some(s) => {
            let p: f64 = s.parse().map_err(|_| {
                format!(
                    "ERROR: --excess-het must be a p-value in (0, 1] or 'auto', got '{}'",
                    s
                )
            })?;
            if !(p > 0.0 && p <= 1.0) {
                return Err(format!(
                    "ERROR: --excess-het must be in (0, 1] or 'auto', got {}",
                    p
                ));
            }
            Some(Ok(p))
        }
    };
    if args.window == 0 {
        return Err(format!(
            "ERROR: --window must be a positive integer, got {}",
            args.window
        ));
    }
    if let Some(kb) = args.window_kb {
        if !(kb.is_finite() && kb > 0.0) {
            return Err(format!("ERROR: --window-kb must be a positive number, got {}", kb));
        }
    }
    let window = match args.window_kb {
        Some(kb) => Window::Bp((kb * 1000.0).round() as u64),
        None => Window::Variants(args.window),
    };
    if args.max_cleanup == 0 {
        return Err(format!(
            "ERROR: --max-cleanup must be a positive integer, got {}",
            args.max_cleanup
        ));
    }

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .map_err(|e| format!("ERROR: failed to configure thread pool: {}", e))?;
    }

    for ext in [".bed", ".bim", ".fam"] {
        let p = format!("{}{}", args.bfile, ext);
        if !Path::new(&p).exists() {
            return Err(format!("ERROR: {} not found", p));
        }
    }

    let abs_bfile = fs::canonicalize(format!("{}.bed", args.bfile)).ok();
    let out_path = Path::new(&args.out);
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).map_err(|e| {
                format!("ERROR: failed to create output dir {}: {}", parent.display(), e)
            })?;
        }
    }
    let abs_out_bed = fs::canonicalize(format!("{}.bed", args.out)).ok();
    if abs_bfile.is_some() && abs_bfile == abs_out_bed {
        return Err("ERROR: --out prefix cannot match --bfile prefix".to_string());
    }

    let version = env!("CARGO_PKG_VERSION").to_string();
    let weighting: Weighting = args.weighting.into();
    let mut messages: Vec<String> = Vec::new();

    // ---- 1. Read .fam & .bim ------------------------------------------------
    let fam = Fam::read(Path::new(&format!("{}.fam", args.bfile)))?;
    let bim = Bim::read(Path::new(&format!("{}.bim", args.bfile)))?;
    let n_auto = bim.variants.len();

    // A physical window needs real, sorted positions on every autosomal SNP
    if let Window::Bp(_) = window {
        let n_nopos = bim.variants.iter().filter(|v| v.pos == 0).count();
        if n_nopos > 0 {
            return Err(format!(
                "ERROR: --window-kb needs base-pair positions, but {} autosomal SNP(s) have position 0 or none in the .bim",
                fmt_int(n_nopos)
            ));
        }
        let mut last: [Option<(u64, &str)>; 23] = [None; 23];
        for v in &bim.variants {
            let c = v.chrom as usize;
            if let Some((p, prev)) = last[c] {
                if v.pos < p {
                    return Err(format!(
                        "ERROR: --window-kb needs the .bim sorted by position within each chromosome; {} (chr {}, {}) comes after {} ({})\n  Sort first, e.g. plink --bfile <in> --make-bed --out <sorted>",
                        v.id, v.chrom, v.pos, prev, p
                    ));
                }
            }
            last[c] = Some((v.pos, v.id.as_str()));
        }
    }
    let n_fam = fam.samples.len();

    eprintln!("ldprio {}", version);
    eprintln!(
        "[1/5] input: {} autosomal SNPs ({} non-autosomal dropped), {} samples",
        fmt_int(n_auto),
        fmt_int(bim.n_input - n_auto),
        fmt_int(n_fam)
    );

    // ---- 2. Resolve samples & stream .bed (with optional LD-sample QC) ------
    let prio_indices = resolve_samples(&args.priority_samples, &fam, "--priority-samples")?;
    let ld_indices = resolve_samples(&args.ld_samples, &fam, "--ld-samples")?;

    let (min_excess_het_p, min_excess_het_fis) = match excess_het_mode {
        None => (None, None),
        Some(Ok(p)) => (Some(p), None),
        Some(Err(())) => {
            let p = if ld_indices.len() < 100 { 1e-4 } else { 1e-6 };
            (Some(p), Some(0.20))
        }
    };

    let qc = QcFilters {
        max_geno: args.geno,
        min_maf: args.maf,
        min_excess_het_p,
        min_excess_het_fis,
    };
    let (variants, prio_calls, ld_matrix, qc_stats) = scan_bed(
        Path::new(&format!("{}.bed", args.bfile)),
        n_fam,
        bim.n_input,
        &bim.variants,
        &prio_indices,
        &ld_indices,
        qc,
    )?;
    let n_snps = variants.len();

    let qc_summary = if qc.is_active() {
        let mut parts = Vec::new();
        if let Some(g) = args.geno {
            parts.push(format!("--geno {} [{}]", g, fmt_int(qc_stats.dropped_geno)));
        }
        if let Some(m) = args.maf {
            parts.push(format!("--maf {} [{}]", m, fmt_int(qc_stats.dropped_maf)));
        }
        if let Some(p) = min_excess_het_p {
            let tag = if let Some(ex) = min_excess_het_fis {
                format!("--excess-het auto (p<{:e}, -F_IS>{})", p, ex)
            } else {
                format!("--excess-het {:e}", p)
            };
            parts.push(format!("{} [{}]", tag, fmt_int(qc_stats.dropped_excess_het)));
        }
        let summary = parts.join(", ");
        eprintln!(
            "      LD-sample QC: {} SNPs passed ({} dropped: {})",
            fmt_int(n_snps),
            fmt_int(qc_stats.total_dropped()),
            summary
        );
        Some(summary)
    } else {
        None
    };

    if n_snps == 0 {
        return Err("ERROR: no autosomal SNPs remain after LD-sample QC filters".to_string());
    }

    // Disambiguate priority sample labels if duplicate IIDs exist across families
    let mut prio_iid_counts: HashMap<&str, usize> = HashMap::new();
    for &idx in &prio_indices {
        *prio_iid_counts.entry(fam.samples[idx].1.as_str()).or_insert(0) += 1;
    }
    let prio_labels: Vec<String> = prio_indices
        .iter()
        .map(|&idx| {
            let (fid, iid) = &fam.samples[idx];
            if prio_iid_counts[iid.as_str()] > 1 {
                format!("{}/{}", fid, iid)
            } else {
                iid.clone()
            }
        })
        .collect();

    let empty: Vec<String> = prio_calls
        .coverage
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c == 0)
        .map(|(i, _)| prio_labels[i].clone())
        .collect();
    if !empty.is_empty() {
        emit(
            &mut messages,
            &format!(
                "      WARNING: {} priority sample(s) have no calls in this panel: {}",
                empty.len(),
                preview_list(&empty, 5)
            ),
        );
    }

    let priority = PriorityPool::build(&prio_calls, n_snps, weighting);
    let pool_size = priority.pool_size();
    if pool_size == 0 {
        return Err("ERROR: priority samples cover no SNPs in this panel".to_string());
    }
    let pool_frac = (pool_size as f64) / (n_snps as f64);
    eprintln!(
        "[2/5] priority pool: {} SNPs covered by the nominated samples ({:.1}% of panel, {} sample(s), {} weighting)",
        fmt_int(pool_size),
        100.0 * pool_frac,
        prio_indices.len(),
        weighting.as_str()
    );
    if pool_frac > POOL_WARN_FRACTION {
        emit(
            &mut messages,
            &format!(
                "      WARNING: priority pool is {:.0}% of the panel. Most blocks are now contested between nominated samples, so each keeps a smaller share of its sites -- consider nominating only the samples that actually need it.",
                100.0 * pool_frac
            ),
        );
    }

    // ---- 3. LD sample checks & monomorphic check ----------------------------
    eprintln!(
        "[3/5] estimating LD from {} samples",
        fmt_int(ld_indices.len())
    );
    let ld_set: HashSet<usize> = ld_indices.iter().copied().collect();
    let mut overlap: Vec<String> = prio_indices
        .iter()
        .filter(|idx| ld_set.contains(idx))
        .map(|&idx| format!("{}/{}", fam.samples[idx].0, fam.samples[idx].1))
        .collect();
    overlap.sort();
    if !overlap.is_empty() {
        emit(
            &mut messages,
            &format!(
                "      WARNING: {} sample(s) appear in both --priority-samples and --ld-samples: {}\n      If these are pseudo-haploid ancient samples this is exactly the LD-distortion case --ld-samples exists to avoid.",
                overlap.len(),
                preview_list(&overlap, 5)
            ),
        );
    }
    if ld_indices.len() < MIN_LD_SAMPLES {
        emit(
            &mut messages,
            &format!(
                "      WARNING: LD is estimated from only {} sample(s). With fewer than {}, r2 estimates are noisy: unrelated SNPs can pass the threshold by chance and linked SNPs can be missed, so the pruned panel may be smaller or less independent than intended. Continuing anyway; consider a larger --ld-samples set.",
                ld_indices.len(),
                MIN_LD_SAMPLES
            ),
        );
    }

    let n_uneval_prio = (0..n_snps)
        .filter(|&k| priority.in_pool(k) && ld_matrix.is_constant[k])
        .count();
    if n_uneval_prio > 0 {
        let frac = 100.0 * (n_uneval_prio as f64) / (pool_size as f64);
        emit(
            &mut messages,
            &format!(
                "      note: {}/{} ({:.1}%) priority-pool SNPs are monomorphic (or zero-variance) in --ld-samples and were never evaluated for LD (always kept; may still be correlated with each other among the samples that matter for them, just invisible to the moderns-only estimate)",
                fmt_int(n_uneval_prio),
                fmt_int(pool_size),
                frac
            ),
        );
    }

    // ---- 4. Initial LD graph + priority-first greedy selection --------------
    let n_tv_input = if args.prefer_transversions {
        variants.iter().filter(|v| v.is_tv).count()
    } else {
        0
    };
    if args.prefer_transversions && n_tv_input == 0 {
        emit(
            &mut messages,
            "      WARNING: --transv found no A/C/G/T transversion SNPs in the .bim (alleles missing or coded 0?); the option has no effect.",
        );
    }
    let prefer_tv = args.prefer_transversions && n_tv_input > 0;

    let all_auto: Vec<usize> = (0..n_snps).collect();
    let (adj, n_pairs) = build_ld_graph_window(&ld_matrix, &all_auto, &variants, window, args.r2);
    let mut current = greedy_select(&all_auto, &adj, &priority, &variants, prefer_tv);
    let n_greedy = current.len();
    let n_prio_after_greedy = current.iter().filter(|&&k| priority.in_pool(k)).count();

    eprintln!(
        "[4/5] LD pairs with r2>={}: {} | greedy kept {} SNPs | priority retained {}/{} ({:.1}%)",
        args.r2,
        fmt_int(n_pairs),
        fmt_int(n_greedy),
        fmt_int(n_prio_after_greedy),
        fmt_int(pool_size),
        100.0 * (n_prio_after_greedy as f64) / (pool_size as f64)
    );

    // ---- Cleanup passes to convergence (all in-memory) ----------------------
    let mut converged = false;
    let mut passes: Vec<(usize, usize, usize, usize)> = Vec::new();
    let mut total_cleanup_removed = 0usize;
    for pass in 1..=args.max_cleanup {
        let (adj_i, n_pairs_i) =
            build_ld_graph_window(&ld_matrix, &current, &variants, window, args.r2);
        let (next, n_removed) = if n_pairs_i == 0 {
            (current, 0)
        } else {
            let kept_i = greedy_select(&current, &adj_i, &priority, &variants, prefer_tv);
            let rem = current.len() - kept_i.len();
            (kept_i, rem)
        };
        current = next;
        total_cleanup_removed += n_removed;
        let n_keep = current.len();
        passes.push((pass, n_pairs_i, n_removed, n_keep));
        let pass_label = if matches!(window, Window::Bp(_)) && pass == 1 {
            "verification".to_string()
        } else {
            format!("cleanup pass {}", pass)
        };
        eprintln!(
            "      {}: {} residual pairs, removed {} -> {} SNPs",
            pass_label,
            fmt_int(n_pairs_i),
            fmt_int(n_removed),
            fmt_int(n_keep)
        );
        if n_pairs_i == 0 {
            converged = true;
            break;
        }
    }

    // ---- Backfill recovery --------------------------------------------------
    // If cleanup removed SNPs, previously dropped SNPs whose only kept LD partner
    // was removed during cleanup may now have no LD with any kept SNP within
    // `--window`. Recover them (priority SNPs first), preserving 0 residual LD.
    let mut backfill: Option<(usize, usize, usize)> = None;
    if let (true, true, Window::Variants(win_n)) = (converged, total_cleanup_removed > 0, window) {
        let candidates =
            find_backfill_candidates(&ld_matrix, &current, &variants, win_n, args.r2);
        if !candidates.is_empty() {
            let (recovered, n_added, n_added_prio) = backfill_select(
                &current,
                &candidates,
                &ld_matrix,
                &priority,
                &variants,
                win_n,
                args.r2,
                prefer_tv,
            );
            if n_added > 0 {
                current = recovered;
                backfill = Some((n_added, n_added_prio, current.len()));
                eprintln!(
                    "      backfill: recovered {} SNPs ({} priority) -> {} SNPs",
                    fmt_int(n_added),
                    fmt_int(n_added_prio),
                    fmt_int(current.len())
                );
            }
        }
    }

    if !converged {
        emit(
            &mut messages,
            &format!(
                "WARNING: not converged after {} cleanup passes; residual LD may remain. Raise --max-cleanup.",
                args.max_cleanup
            ),
        );
    }

    // ---- 5. Outputs ---------------------------------------------------------
    let snplist_path = format!("{}.snplist", args.out);
    {
        let f = File::create(&snplist_path)
            .map_err(|e| format!("ERROR: failed to create {}: {}", snplist_path, e))?;
        let mut w = BufWriter::new(f);
        for &k in &current {
            writeln!(w, "{}", variants[k].id)
                .map_err(|e| format!("ERROR: writing {}: {}", snplist_path, e))?;
        }
    }

    let n_prio_final = current.iter().filter(|&&k| priority.in_pool(k)).count();
    eprintln!(
        "[5/5] final: {} SNPs ({}) | priority retained {}/{} ({:.1}%)",
        fmt_int(current.len()),
        if converged {
            "LD-independent, verified"
        } else {
            "NOT verified"
        },
        fmt_int(n_prio_final),
        fmt_int(pool_size),
        100.0 * (n_prio_final as f64) / (pool_size as f64)
    );
    let n_tv_final = if prefer_tv {
        let n_tv = current.iter().filter(|&&k| variants[k].is_tv).count();
        eprintln!(
            "      transversions: {}/{} kept SNPs ({:.1}%)",
            fmt_int(n_tv),
            fmt_int(current.len()),
            100.0 * (n_tv as f64) / (current.len() as f64)
        );
        Some(n_tv)
    } else {
        None
    };
    eprintln!("      wrote {}", snplist_path);

    if args.make_bed {
        write_bed_subset(&args.bfile, &args.out, n_fam, &variants, &current)?;
        eprintln!("      wrote {}.{{bed,bim,fam}}", args.out);
    }

    if args.report {
        let kept_calls = prio_calls.count_calls_on(&current);
        let sample_rows: Vec<(String, usize, usize)> = prio_labels
            .into_iter()
            .zip(prio_calls.coverage.into_iter())
            .zip(kept_calls.into_iter())
            .map(|((lab, av), kp)| (lab, av, kp))
            .collect();
        let cmd = std::env::args().collect::<Vec<_>>().join(" ");
        let report_path = format!("{}.report.html", args.out);
        write_report(
            &report_path,
            &ReportData {
                version,
                out: args.out.clone(),
                command: cmd,
                runtime_secs: start.elapsed().as_secs_f64(),
                weighting: weighting.as_str().to_string(),
                r2: args.r2,
                window: window.describe(),
                single_pass: matches!(window, Window::Bp(_)),
                n_input: bim.n_input,
                n_auto,
                n_snps,
                qc_dropped: qc_stats.total_dropped(),
                qc_summary,
                n_fam,
                n_ld: ld_indices.len(),
                pool: pool_size,
                prio_final: n_prio_final,
                greedy: n_greedy,
                passes,
                backfill,
                final_n: current.len(),
                converged,
                samples: sample_rows,
                tv_share: n_tv_final.map(|n_tv| (n_tv, current.len())),
                messages,
            },
        )?;
        eprintln!("      wrote {}", report_path);
    }

    if !converged {
        process::exit(1);
    }
    Ok(())
}
