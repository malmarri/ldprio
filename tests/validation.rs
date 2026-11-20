use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use ldprio::ld::{build_ld_graph, LdMatrix};
use ldprio::plink::{hwe_exact_p, scan_bed, Bim, Fam, QcFilters};

/// Simple deterministic PRNG (xorshift64*) so tests are 100% reproducible with zero deps.
struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

/// Write a PLINK 1 SNP-major fileset from a `(n_snps x n_samples)` matrix of codes:
/// 0 = hom A1 (dosage 0), 1 = missing, 2 = het (dosage 1), 3 = hom A2 (dosage 2).
fn write_plink_fileset(
    prefix: &Path,
    fam_rows: &[(&str, &str)],
    bim_rows: &[(&str, &str, usize, &str, &str)], // (chrom, id, bp, a1, a2)
    geno: &[Vec<u8>],
) {
    let n_fam = fam_rows.len();
    let rec = (n_fam + 3) / 4;

    // .fam
    {
        let mut w = BufWriter::new(File::create(prefix.with_extension("fam")).unwrap());
        for (fid, iid) in fam_rows {
            writeln!(w, "{}\t{}\t0\t0\t0\t-9", fid, iid).unwrap();
        }
    }
    // .bim
    {
        let mut w = BufWriter::new(File::create(prefix.with_extension("bim")).unwrap());
        for (chrom, id, bp, a1, a2) in bim_rows {
            writeln!(w, "{}\t{}\t0\t{}\t{}\t{}", chrom, id, bp, a1, a2).unwrap();
        }
    }
    // .bed
    {
        let mut w = BufWriter::new(File::create(prefix.with_extension("bed")).unwrap());
        w.write_all(&[0x6c, 0x1b, 0x01]).unwrap();
        let mut buf = vec![0u8; rec];
        for snp in geno {
            buf.fill(0);
            for (s, &code) in snp.iter().enumerate() {
                buf[s / 4] |= (code & 3) << ((s % 4) * 2);
            }
            w.write_all(&buf).unwrap();
        }
    }
}

fn ldprio_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_ldprio"))
}

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("ldprio_val_{}_{}", name, std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// Naive independent f64 Pearson r^2 over jointly non-missing dosages.
fn naive_r2(a: &[u8], b: &[u8]) -> Option<f64> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (&ca, &cb) in a.iter().zip(b.iter()) {
        if ca == 1 || cb == 1 {
            continue;
        }
        let xa = match ca {
            0 => 0.0,
            2 => 1.0,
            3 => 2.0,
            _ => unreachable!(),
        };
        let xb = match cb {
            0 => 0.0,
            2 => 1.0,
            3 => 2.0,
            _ => unreachable!(),
        };
        xs.push(xa);
        ys.push(xb);
    }
    let n = xs.len() as f64;
    if n < 2.0 {
        return None;
    }
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut vx = 0.0;
    let mut vy = 0.0;
    let mut cov = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mx;
        let dy = y - my;
        vx += dx * dx;
        vy += dy * dy;
        cov += dx * dy;
    }
    if vx <= 0.0 || vy <= 0.0 {
        return None;
    }
    Some((cov * cov) / (vx * vy))
}

#[test]
fn test_bitwise_r2_exact_against_naive_with_heavy_missingness_and_edge_cases() {
    let dir = temp_dir("r2_exact");
    let mut rng = Rng::new(42);

    // 137 samples (non-multiple of 64 to stress tail word masking!)
    let n_samples = 137;
    let fam_rows: Vec<(String, String)> = (0..n_samples)
        .map(|i| ("F".to_string(), format!("S{}", i)))
        .collect();
    let fam_refs: Vec<(&str, &str)> = fam_rows
        .iter()
        .map(|(f, i)| (f.as_str(), i.as_str()))
        .collect();

    let n_snps = 60;
    let mut geno: Vec<Vec<u8>> = Vec::with_capacity(n_snps);

    // Special edge-case SNPs first:
    geno.push(vec![0; n_samples]); // 0: monomorphic hom A1
    geno.push(vec![3; n_samples]); // 1: monomorphic hom A2
    geno.push(vec![2; n_samples]); // 2: 100% het (MAF=0.5, variance=0!)
    geno.push(vec![1; n_samples]); // 3: all missing
    let mut single_call = vec![1; n_samples];
    single_call[77] = 2;
    geno.push(single_call); // 4: only 1 called sample

    // 5: base polymorphic SNP with 25% missingness
    let base: Vec<u8> = (0..n_samples)
        .map(|_| {
            let r = rng.next_f64();
            if r < 0.25 {
                1
            } else if r < 0.50 {
                0
            } else if r < 0.75 {
                2
            } else {
                3
            }
        })
        .collect();
    // 6: exact positive clone (r^2 = 1.0)
    let pos_clone = base.clone();
    // 7: exact flipped clone (r = -1.0, r^2 = 1.0)
    let neg_clone: Vec<u8> = base
        .iter()
        .map(|&c| match c {
            0 => 3,
            3 => 0,
            other => other,
        })
        .collect();
    geno.push(base);
    geno.push(pos_clone);
    geno.push(neg_clone);

    // Remaining SNPs: random genotypes with 0%..45% missingness + some correlated copies
    while geno.len() < n_snps {
        let miss_rate = rng.next_f64() * 0.45;
        let copy_from = if rng.next_f64() < 0.35 && geno.len() > 8 {
            Some(5 + (rng.next_u64() as usize) % (geno.len() - 5))
        } else {
            None
        };
        let snp: Vec<u8> = (0..n_samples)
            .map(|s| {
                if rng.next_f64() < miss_rate {
                    return 1;
                }
                if let Some(src) = copy_from {
                    if rng.next_f64() < 0.85 && geno[src][s] != 1 {
                        return geno[src][s];
                    }
                }
                match rng.next_u64() % 3 {
                    0 => 0,
                    1 => 2,
                    _ => 3,
                }
            })
            .collect();
        geno.push(snp);
    }

    let bim_ids: Vec<String> = (0..n_snps).map(|k| format!("rs{}", k)).collect();
    let bim_rows: Vec<(&str, &str, usize, &str, &str)> = bim_ids
        .iter()
        .enumerate()
        .map(|(k, id)| ("1", id.as_str(), (k + 1) * 1000, "A", "C"))
        .collect();

    let prefix = dir.join("panel");
    write_plink_fileset(&prefix, &fam_refs, &bim_rows, &geno);

    let bim = Bim::read(&prefix.with_extension("bim")).unwrap();
    let all_samples: Vec<usize> = (0..n_samples).collect();
    let (_, _, ld_mat, _): (_, _, LdMatrix, _) = scan_bed(
        &prefix.with_extension("bed"),
        n_samples,
        bim.n_input,
        &bim.variants,
        &all_samples[..5],
        &all_samples,
        QcFilters::default(),
    )
    .unwrap();

    // Check constant-SNP flags
    assert!(ld_mat.is_constant[0], "monomorphic hom A1 must be constant");
    assert!(ld_mat.is_constant[1], "monomorphic hom A2 must be constant");
    assert!(ld_mat.is_constant[2], "100% het must be constant (var=0)");
    assert!(ld_mat.is_constant[3], "all-missing must be constant");
    assert!(ld_mat.is_constant[4], "single-call must be constant");

    // Check perfect LD positive & negative clones at threshold = 1.0
    assert!(ld_mat.pair_r2_ge(5, 6, 1.0));
    assert!(ld_mat.pair_r2_ge(5, 7, 1.0));

    // Exhaustively compare ALL pairs against naive f64 Pearson r^2 across multiple thresholds
    for &thresh in &[0.1, 0.25, 0.4, 0.7, 0.99, 1.0] {
        for a in 0..n_snps {
            for b in (a + 1)..n_snps {
                let got = ld_mat.pair_r2_ge(a, b, thresh);
                let expected_r2 = naive_r2(&geno[a], &geno[b]);
                let expected = expected_r2.map(|r2| r2 >= thresh - 1e-12).unwrap_or(false);
                let strictly_below = expected_r2.map(|r2| r2 < thresh - 1e-12).unwrap_or(true);
                if got {
                    assert!(
                        expected,
                        "false positive at ({}, {}) thresh={}: naive={:?}",
                        a, b, thresh, expected_r2
                    );
                } else {
                    assert!(
                        strictly_below || (expected_r2.unwrap() - thresh).abs() < 1e-12,
                        "false negative at ({}, {}) thresh={}: naive={:?}",
                        a, b, thresh, expected_r2
                    );
                }
            }
        }
    }

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn test_multipass_cleanup_and_non_autosomal_interleaving_and_transversions() {
    let dir = temp_dir("multipass_edge");
    let mut rng = Rng::new(1337);

    // 60 LD diploids + 2 priority samples (with duplicate IID across different FIDs!)
    // plus 1 zero-call priority sample
    let mut fam_owned: Vec<(String, String)> = (0..60)
        .map(|i| ("LD".to_string(), format!("M{}", i)))
        .collect();
    fam_owned.push(("FAM1".to_string(), "DUP".to_string()));
    fam_owned.push(("FAM2".to_string(), "DUP".to_string()));
    fam_owned.push(("FAM3".to_string(), "EMPTY".to_string()));
    let n_fam = fam_owned.len();
    let fam_refs: Vec<(&str, &str)> = fam_owned
        .iter()
        .map(|(f, i)| (f.as_str(), i.as_str()))
        .collect();

    // Build a chain on chr1 where window=5 misses long-range LD at distance 7..12
    // until intermediate SNPs are dropped (forces multiple cleanup passes!).
    // Also interleave non-autosomal SNPs (X, MT, chrUn, 0) at start, middle, and end,
    // and use zero-padded "01" and "chr02".
    let mut bim_owned: Vec<(String, String, usize, String, String)> = Vec::new();
    let mut geno: Vec<Vec<u8>> = Vec::new();

    let mut push_snp = |chrom: &str, id: &str, bp: usize, a1: &str, a2: &str, row: Vec<u8>| {
        bim_owned.push((
            chrom.to_string(),
            id.to_string(),
            bp,
            a1.to_string(),
            a2.to_string(),
        ));
        geno.push(row);
    };

    // Leading non-autosomal SNPs
    push_snp("X", "x_start", 100, "A", "T", vec![0; n_fam]);
    push_snp("0", "zero_start", 200, "A", "T", vec![0; n_fam]);

    // Base genotype pattern A (used at positions 0 and 12 — outside window=6 initially,
    // inside window=6 after intermediate block SNPs are pruned!)
    let pat_a: Vec<u8> = (0..60)
        .map(|_| match rng.next_u64() % 3 {
            0 => 0,
            1 => 2,
            _ => 3,
        })
        .collect();
    // Base genotype pattern B (independent of A, used for filler blocks)
    let pat_b: Vec<u8> = (0..60)
        .map(|_| match rng.next_u64() % 3 {
            0 => 0,
            1 => 2,
            _ => 3,
        })
        .collect();
    let pat_c: Vec<u8> = (0..60)
        .map(|_| match rng.next_u64() % 3 {
            0 => 0,
            1 => 2,
            _ => 3,
        })
        .collect();

    // Helper to extend 60 LD genotypes with 3 priority sample calls
    let make_row = |ld_pat: &[u8], p1: u8, p2: u8| {
        let mut r = ld_pat.to_vec();
        r.push(p1);
        r.push(p2);
        r.push(1); // EMPTY always missing
        r
    };

    // SNP 0 (chr01): pattern A, non-priority
    push_snp("01", "a_nonprio", 1000, "A", "G", make_row(&pat_a, 1, 1));
    // SNPs 1..6: 6 copies of pattern B (will collapse to 1 SNP in pass 1!)
    for i in 0..6 {
        let (a1, a2) = if i == 3 { ("A", "C") } else { ("C", "T") }; // only i==3 is transversion!
        push_snp(
            "01",
            &format!("b_block_{}", i),
            2000 + i * 100,
            a1,
            a2,
            make_row(&pat_b, 1, 1),
        );
    }
    // Interleaved non-autosomal SNP right in the middle!
    push_snp("MT", "mt_mid", 5000, "A", "C", make_row(&pat_a, 0, 0));

    // SNPs 7..12: 6 copies of pattern C (will collapse to 1 SNP in pass 1!)
    for i in 0..6 {
        push_snp(
            "01",
            &format!("c_block_{}", i),
            6000 + i * 100,
            "G",
            "A",
            make_row(&pat_c, 1, 1),
        );
    }
    // SNP 13: pattern A again, covered by priority sample FAM1/DUP!
    // Distance from SNP 0 is 13 (> window=6), so pass 1 keeps BOTH a_nonprio and a_prio.
    // After pass 1 collapses blocks B and C to 1 SNP each, distance becomes 3 (< window=6),
    // so cleanup pass 1 MUST catch (a_nonprio, a_prio) and KEEP a_prio (priority-aware)!
    push_snp("01", "a_prio", 9000, "A", "G", make_row(&pat_a, 0, 1));

    // Chr2: tied block with FAM2/DUP covering both a transition and a transversion
    push_snp("chr02", "c2_ti", 1000, "C", "T", make_row(&pat_b, 1, 0));
    push_snp("chr02", "c2_tv", 2000, "A", "T", make_row(&pat_b, 1, 0));

    // Trailing non-autosomal SNP
    push_snp("chrUn_gl000220", "un_end", 100, "A", "C", make_row(&pat_a, 0, 0));

    let bim_refs: Vec<(&str, &str, usize, &str, &str)> = bim_owned
        .iter()
        .map(|(c, id, bp, a1, a2)| (c.as_str(), id.as_str(), *bp, a1.as_str(), a2.as_str()))
        .collect();
    let prefix = dir.join("multi");
    write_plink_fileset(&prefix, &fam_refs, &bim_refs, &geno);

    // Write sample lists
    let ld_list = dir.join("ld.txt");
    fs::write(
        &ld_list,
        (0..60).map(|i| format!("M{}\n", i)).collect::<String>(),
    )
    .unwrap();
    let prio_list = dir.join("prio.txt");
    fs::write(&prio_list, "FAM1\tDUP\nFAM2\tDUP\nFAM3\tEMPTY\n").unwrap();

    // 1. Ambiguous IID without FID must fail
    let bad_prio = dir.join("bad_prio.txt");
    fs::write(&bad_prio, "DUP\n").unwrap();
    let bad_out = Command::new(ldprio_bin())
        .args([
            "--bfile",
            prefix.to_str().unwrap(),
            "--priority-samples",
            bad_prio.to_str().unwrap(),
            "--ld-samples",
            ld_list.to_str().unwrap(),
            "--out",
            dir.join("bad").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!bad_out.status.success(), "ambiguous IID should fail");

    // 2. Run with window=6, --prefer-transversions, --make-bed, --report
    let out_prefix = dir.join("pruned");
    let res = Command::new(ldprio_bin())
        .args([
            "--bfile",
            prefix.to_str().unwrap(),
            "--priority-samples",
            prio_list.to_str().unwrap(),
            "--ld-samples",
            ld_list.to_str().unwrap(),
            "--window",
            "6",
            "--r2",
            "0.4",
            "--prefer-transversions",
            "--make-bed",
            "--report",
            "--out",
            out_prefix.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(res.status.success(), "ldprio failed:\n{}", stderr);

    let kept: Vec<String> = fs::read_to_string(out_prefix.with_extension("snplist"))
        .unwrap()
        .lines()
        .map(|s| s.to_string())
        .collect();

    // Cleanup must have resolved (a_nonprio, a_prio) in favour of a_prio!
    assert!(kept.contains(&"a_prio".to_string()), "a_prio must survive cleanup: {:?}", kept);
    assert!(!kept.contains(&"a_nonprio".to_string()), "a_nonprio must be dropped in cleanup: {:?}", kept);
    // Transversion preference must pick b_block_3 in block B and c2_tv on chr2!
    assert!(kept.contains(&"b_block_3".to_string()), "b_block_3 (TV) must win block B: {:?}", kept);
    assert!(kept.contains(&"c2_tv".to_string()), "c2_tv (TV) must win tied chr2 block: {:?}", kept);
    assert!(!kept.contains(&"c2_ti".to_string()), "c2_ti must be dropped: {:?}", kept);
    // Non-autosomal SNPs must be excluded
    assert!(!kept.iter().any(|s| s == "x_start" || s == "mt_mid" || s == "un_end" || s == "zero_start"));

    // Verify HTML report disambiguates FAM1/DUP and FAM2/DUP
    let html = fs::read_to_string(format!("{}.report.html", out_prefix.display())).unwrap();
    assert!(html.contains("FAM1/DUP"), "report must disambiguate FAM1/DUP");
    assert!(html.contains("FAM2/DUP"), "report must disambiguate FAM2/DUP");
    assert!(html.contains("EMPTY"), "report must include EMPTY sample");

    // 3. Round-trip verification: re-scan the written --make-bed output and verify 0 LD pairs!
    let out_fam = Fam::read(&out_prefix.with_extension("fam")).unwrap();
    let out_bim = Bim::read(&out_prefix.with_extension("bim")).unwrap();
    let ld_idx: Vec<usize> = (0..60).collect();
    let (_, _, out_ld, _) = scan_bed(
        &out_prefix.with_extension("bed"),
        out_fam.samples.len(),
        out_bim.n_input,
        &out_bim.variants,
        &[60, 61, 62],
        &ld_idx,
        QcFilters::default(),
    )
    .unwrap();
    let all_out: Vec<usize> = (0..out_bim.variants.len()).collect();
    let (_, residual_pairs) = build_ld_graph(&out_ld, &all_out, &out_bim.variants, 6, 0.4);
    assert_eq!(residual_pairs, 0, "pruned .bed output must have 0 residual LD pairs");

    // 4. --max-cleanup 1 must fail (exit non-zero) because cleanup pass 1 removes a_nonprio
    // and never gets a pass 2 to verify 0 pairs
    let out_no_conv = dir.join("noconv");
    let res_nc = Command::new(ldprio_bin())
        .args([
            "--bfile",
            prefix.to_str().unwrap(),
            "--priority-samples",
            prio_list.to_str().unwrap(),
            "--ld-samples",
            ld_list.to_str().unwrap(),
            "--window",
            "6",
            "--r2",
            "0.4",
            "--max-cleanup",
            "1",
            "--out",
            out_no_conv.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!res_nc.status.success(), "unconverged run must exit non-zero");

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn test_large_scale_100k_snps_stability_and_thread_determinism() {
    let dir = temp_dir("scale_100k");
    let mut rng = Rng::new(20261120);

    // 100,000 SNPs across 22 chromosomes, 200 LD samples + 50 priority samples (>40 -> histogram!)
    let n_ld = 200;
    let n_prio = 50;
    let n_fam = n_ld + n_prio;
    let n_snps = 100_000;

    let mut fam_owned = Vec::with_capacity(n_fam);
    for i in 0..n_ld {
        fam_owned.push(("MOD".to_string(), format!("M{}", i)));
    }
    for i in 0..n_prio {
        fam_owned.push(("ANC".to_string(), format!("A{}", i)));
    }
    let fam_refs: Vec<(&str, &str)> = fam_owned
        .iter()
        .map(|(f, i)| (f.as_str(), i.as_str()))
        .collect();

    // Generate Markov-chain decaying LD across 22 chromosomes + 5% missingness in LD samples
    // and 1%..15% varying coverage in ancient priority samples.
    let rec = (n_fam + 3) / 4;
    let prefix = dir.join("big100k");

    // Write .fam
    {
        let mut w = BufWriter::new(File::create(prefix.with_extension("fam")).unwrap());
        for (fid, iid) in &fam_refs {
            writeln!(w, "{}\t{}\t0\t0\t0\t-9", fid, iid).unwrap();
        }
    }

    // Stream-write .bim and .bed directly (low memory)
    let mut bim_w = BufWriter::with_capacity(1 << 20, File::create(prefix.with_extension("bim")).unwrap());
    let mut bed_w = BufWriter::with_capacity(1 << 20, File::create(prefix.with_extension("bed")).unwrap());
    bed_w.write_all(&[0x6c, 0x1b, 0x01]).unwrap();

    // Per-ancient coverage probabilities (0.5% to 15%)
    let anc_cov: Vec<f64> = (0..n_prio)
        .map(|i| 0.005 + 0.145 * (i as f64) / ((n_prio - 1) as f64))
        .collect();

    let mut cur_hap1 = vec![0u8; n_ld];
    let mut cur_hap2 = vec![0u8; n_ld];
    let mut row_buf = vec![0u8; rec];
    let alleles = [("A", "G"), ("C", "T"), ("A", "C"), ("G", "T")];

    let snps_per_chr = n_snps / 22;
    for k in 0..n_snps {
        let chrom = (k / snps_per_chr).min(21) + 1;
        let is_chr_start = k % snps_per_chr == 0 && k / snps_per_chr < 22;
        let (a1, a2) = alleles[(rng.next_u64() as usize) % 4];
        writeln!(bim_w, "{}\tsnp_{}\t0\t{}\t{}\t{}", chrom, k, (k % snps_per_chr + 1) * 100, a1, a2).unwrap();

        row_buf.fill(0);
        // Evolve haplotypes with recombination prob 0.08 per SNP (strong local LD blocks)
        let recomb = is_chr_start;
        for s in 0..n_ld {
            if recomb || rng.next_f64() < 0.08 {
                cur_hap1[s] = (rng.next_u64() & 1) as u8;
            }
            if recomb || rng.next_f64() < 0.08 {
                cur_hap2[s] = (rng.next_u64() & 1) as u8;
            }
            let code = if rng.next_f64() < 0.03 {
                1u8 // 3% missing in LD samples
            } else {
                match cur_hap1[s] + cur_hap2[s] {
                    0 => 0,
                    1 => 2,
                    _ => 3,
                }
            };
            row_buf[s / 4] |= code << ((s % 4) * 2);
        }
        // Ancient pseudo-haploid calls
        for a in 0..n_prio {
            let s = n_ld + a;
            let code = if rng.next_f64() < anc_cov[a] {
                if (rng.next_u64() & 1) == 0 { 0 } else { 3 }
            } else {
                1
            };
            row_buf[s / 4] |= code << ((s % 4) * 2);
        }
        bed_w.write_all(&row_buf).unwrap();
    }
    bim_w.flush().unwrap();
    bed_w.flush().unwrap();

    let ld_list = dir.join("ld.txt");
    fs::write(&ld_list, (0..n_ld).map(|i| format!("M{}\n", i)).collect::<String>()).unwrap();
    let prio_list = dir.join("prio.txt");
    fs::write(&prio_list, (0..n_prio).map(|i| format!("A{}\n", i)).collect::<String>()).unwrap();

    // Run with --threads 1 and with default (all cores) using default --weighting adaptive,
    // verify bit-for-bit identical output!
    let out_t1 = dir.join("out_t1");
    let res_t1 = Command::new(ldprio_bin())
        .args([
            "--bfile", prefix.to_str().unwrap(),
            "--priority-samples", prio_list.to_str().unwrap(),
            "--ld-samples", ld_list.to_str().unwrap(),
            "--window", "200", "--r2", "0.4",
            "--weighting", "adaptive",
            "--threads", "1",
            "--out", out_t1.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(res_t1.status.success(), "1-thread run failed:\n{}", String::from_utf8_lossy(&res_t1.stderr));

    let out_mt = dir.join("out_mt");
    let t0 = std::time::Instant::now();
    let res_mt = Command::new(ldprio_bin())
        .args([
            "--bfile", prefix.to_str().unwrap(),
            "--priority-samples", prio_list.to_str().unwrap(),
            "--ld-samples", ld_list.to_str().unwrap(),
            "--window", "200", "--r2", "0.4",
            "--weighting", "adaptive",
            "--make-bed", "--report",
            "--out", out_mt.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    let elapsed = t0.elapsed();
    let stderr_mt = String::from_utf8_lossy(&res_mt.stderr);
    assert!(res_mt.status.success(), "MT run failed:\n{}", stderr_mt);
    eprintln!("[100k benchmark] wall time: {:.2}s\n{}", elapsed.as_secs_f64(), stderr_mt);

    // 1. Thread determinism: 1-thread and multi-thread snplists MUST be byte-for-byte identical
    let s1 = fs::read(out_t1.with_extension("snplist")).unwrap();
    let sm = fs::read(out_mt.with_extension("snplist")).unwrap();
    assert_eq!(s1, sm, "single-threaded and multi-threaded outputs must be identical");

    // 2. HTML histogram report rendered cleanly (>40 priority samples)
    let html = fs::read_to_string(format!("{}.report.html", out_mt.display())).unwrap();
    assert!(html.contains("class=\"hist\""), "expected histogram for 50 priority samples");
    assert!(html.contains("15 samples with the fewest SNPs kept"));

    // 3. Independent round-trip verification of the written --make-bed output: 0 residual LD pairs
    let out_fam = Fam::read(&out_mt.with_extension("fam")).unwrap();
    let out_bim = Bim::read(&out_mt.with_extension("bim")).unwrap();
    let ld_idx: Vec<usize> = (0..n_ld).collect();
    let prio_idx: Vec<usize> = (n_ld..n_fam).collect();
    let (_, _, out_ld, _) = scan_bed(
        &out_mt.with_extension("bed"),
        out_fam.samples.len(),
        out_bim.n_input,
        &out_bim.variants,
        &prio_idx,
        &ld_idx,
        QcFilters::default(),
    )
    .unwrap();
    let all_out: Vec<usize> = (0..out_bim.variants.len()).collect();
    let (_, residual) = build_ld_graph(&out_ld, &all_out, &out_bim.variants, 200, 0.4);
    assert_eq!(residual, 0, "100k pruned output must have 0 residual LD pairs");

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn test_ld_sample_qc_filters_geno_maf_hwe() {
    // 1. Unit-check Wigginton et al. (2005) exact HWE p-values against known reference values
    // Perfect HWE: (50 het, 25 hom1, 25 hom2) -> p == 1.0
    assert!((hwe_exact_p(50, 25, 25) - 1.0).abs() < 1e-12);
    // Monomorphic: (0 het, 100 hom1, 0 hom2) -> p == 1.0
    assert!((hwe_exact_p(0, 100, 0) - 1.0).abs() < 1e-12);
    // Extreme excess homozygosity: (0 het, 50 hom1, 50 hom2) -> p << 1e-20
    assert!(hwe_exact_p(0, 50, 50) < 1e-20);
    // Extreme excess heterozygosity: (100 het, 0 hom1, 0 hom2) -> p << 1e-20
    assert!(hwe_exact_p(100, 0, 0) < 1e-20);
    // Classic Wigginton (2005) example: (hets=57, hom1=14, hom2=50) -> p ~ 0.842
    let p_wig = hwe_exact_p(57, 14, 50);
    assert!(p_wig > 0.80 && p_wig < 0.90, "unexpected HWE p: {}", p_wig);

    // 2. End-to-end CLI test of --geno, --maf, --hwe (calculated ONLY on --ld-samples!)
    let dir = temp_dir("qc_filters");
    let n_ld = 100usize;
    let n_prio = 2usize;
    let mut fam_owned = Vec::new();
    for i in 0..n_ld {
        fam_owned.push(("LD".to_string(), format!("M{}", i)));
    }
    for i in 0..n_prio {
        fam_owned.push(("ANC".to_string(), format!("A{}", i)));
    }
    let fam_refs: Vec<(&str, &str)> = fam_owned
        .iter()
        .map(|(f, i)| (f.as_str(), i.as_str()))
        .collect();

    // 5 SNPs on chr1:
    // snp_good: 0% miss, MAF=0.5, exact HWE (25 hom1, 50 het, 25 hom2) -> PASS
    let mut snp_good = Vec::with_capacity(n_ld + n_prio);
    snp_good.extend(std::iter::repeat(0).take(25));
    snp_good.extend(std::iter::repeat(2).take(50));
    snp_good.extend(std::iter::repeat(3).take(25));
    snp_good.extend([0, 3]); // priority samples called

    // snp_high_miss: 20% missing in LD samples (> --geno 0.10) -> FAIL geno
    let mut snp_high_miss = snp_good.clone();
    for i in 0..20 {
        snp_high_miss[i] = 1;
    }

    // snp_low_maf: 1 het out of 100 LD samples (MAF = 0.005 < --maf 0.01) -> FAIL maf
    let mut snp_low_maf = vec![0u8; n_ld + n_prio];
    snp_low_maf[0] = 2;
    snp_low_maf[n_ld] = 3; // even though ancient has alt allele, LD-sample MAF is 0.005!

    // snp_bad_hwe: 50 hom1, 0 het, 50 hom2 (p < 1e-20 < --hwe 1e-6) -> FAIL hwe
    let mut snp_bad_hwe = Vec::with_capacity(n_ld + n_prio);
    snp_bad_hwe.extend(std::iter::repeat(0).take(50));
    snp_bad_hwe.extend(std::iter::repeat(3).take(50));
    snp_bad_hwe.extend([0, 3]);

    // snp_good2: 5% missing (<= --geno 0.10), MAF=0.25, HWE ok -> PASS
    let mut snp_good2 = Vec::with_capacity(n_ld + n_prio);
    snp_good2.extend(std::iter::repeat(1).take(5));
    snp_good2.extend(std::iter::repeat(0).take(54));
    snp_good2.extend(std::iter::repeat(2).take(35));
    snp_good2.extend(std::iter::repeat(3).take(6));
    snp_good2.extend([3, 1]);

    let bim_rows = [
        ("1", "snp_good", 1000, "A", "C"),
        ("1", "snp_high_miss", 2000, "A", "G"),
        ("1", "snp_low_maf", 3000, "C", "T"),
        ("1", "snp_bad_hwe", 4000, "G", "T"),
        ("2", "snp_good2", 5000, "A", "T"),
    ];
    let geno = vec![snp_good, snp_high_miss, snp_low_maf, snp_bad_hwe, snp_good2];

    let prefix = dir.join("qc_panel");
    write_plink_fileset(&prefix, &fam_refs, &bim_rows, &geno);
    let ld_list = dir.join("ld.txt");
    fs::write(&ld_list, (0..n_ld).map(|i| format!("M{}\n", i)).collect::<String>()).unwrap();
    let prio_list = dir.join("prio.txt");
    fs::write(&prio_list, "A0\nA1\n").unwrap();

    let out_prefix = dir.join("qc_out");
    let res = Command::new(ldprio_bin())
        .args([
            "--bfile", prefix.to_str().unwrap(),
            "--priority-samples", prio_list.to_str().unwrap(),
            "--ld-samples", ld_list.to_str().unwrap(),
            "--geno", "0.10",
            "--maf", "0.01",
            "--hwe", "1e-6",
            "--make-bed", "--report",
            "--out", out_prefix.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(res.status.success(), "QC run failed:\n{}", stderr);

    let kept: Vec<String> = fs::read_to_string(out_prefix.with_extension("snplist"))
        .unwrap()
        .lines()
        .map(|s| s.to_string())
        .collect();
    assert_eq!(kept, vec!["snp_good", "snp_good2"], "stderr:\n{}", stderr);

    let html = fs::read_to_string(format!("{}.report.html", out_prefix.display())).unwrap();
    assert!(html.contains("After LD-sample QC"));
    assert!(html.contains("3 dropped"));

    let _ = fs::remove_dir_all(&dir);
}
