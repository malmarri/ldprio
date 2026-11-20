use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use crate::ld::LdMatrix;

#[derive(Debug, Clone)]
pub struct Fam {
    pub samples: Vec<(String, String)>, // (fid, iid)
    pub by_iid: HashMap<String, Vec<usize>>,
    pub by_pair: HashMap<(String, String), usize>,
}

impl Fam {
    pub fn read(path: &Path) -> Result<Self, String> {
        let file = File::open(path)
            .map_err(|e| format!("ERROR: failed to open {}: {}", path.display(), e))?;
        let reader = BufReader::new(file);
        let mut samples = Vec::new();
        let mut by_iid: HashMap<String, Vec<usize>> = HashMap::new();
        let mut by_pair: HashMap<(String, String), usize> = HashMap::new();

        for line in reader.lines() {
            let line = line.map_err(|e| format!("ERROR: reading {}: {}", path.display(), e))?;
            let mut it = line.split_whitespace();
            let (Some(fid), Some(iid)) = (it.next(), it.next()) else {
                continue;
            };
            let idx = samples.len();
            let fid_s = fid.to_string();
            let iid_s = iid.to_string();
            by_iid.entry(iid_s.clone()).or_default().push(idx);
            by_pair.insert((fid_s.clone(), iid_s.clone()), idx);
            samples.push((fid_s, iid_s));
        }
        Ok(Self {
            samples,
            by_iid,
            by_pair,
        })
    }
}

/// Accept 'IID' or 'FID IID' lines; return validated, order-preserving deduplicated
/// sample indices into `fam.samples`.
pub fn resolve_samples(list_path: &Path, fam: &Fam, label: &str) -> Result<Vec<usize>, String> {
    if !list_path.exists() {
        return Err(format!(
            "ERROR: {} file not found: {}",
            label,
            list_path.display()
        ));
    }
    let file = File::open(list_path)
        .map_err(|e| format!("ERROR: failed to open {}: {}", list_path.display(), e))?;
    let reader = BufReader::new(file);

    let mut keep = Vec::new();
    let mut seen = HashSet::new();
    let mut missing = Vec::new();
    let mut ambiguous = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("ERROR: reading {}: {}", list_path.display(), e))?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        if fields.len() >= 2 {
            let key = (fields[0].to_string(), fields[1].to_string());
            if let Some(&idx) = fam.by_pair.get(&key) {
                if seen.insert(idx) {
                    keep.push(idx);
                }
            } else {
                missing.push(format!("{} {}", fields[0], fields[1]));
            }
        } else {
            let iid = fields[0];
            match fam.by_iid.get(iid) {
                Some(matches) if matches.len() == 1 => {
                    let idx = matches[0];
                    if seen.insert(idx) {
                        keep.push(idx);
                    }
                }
                Some(matches) if matches.len() > 1 => {
                    ambiguous.push(iid.to_string());
                }
                _ => {
                    missing.push(iid.to_string());
                }
            }
        }
    }

    if !ambiguous.is_empty() {
        let preview = preview_list(&ambiguous, 5);
        return Err(format!(
            "ERROR: {} sample(s) in {} match more than one FID in the .fam (duplicate IID across families): {}\n  Specify 'FID<tab>IID' for these samples instead.",
            ambiguous.len(),
            label,
            preview
        ));
    }
    if !missing.is_empty() {
        let preview = preview_list(&missing, 5);
        return Err(format!(
            "ERROR: {} sample(s) in {} not found in the .fam: {}",
            missing.len(),
            label,
            preview
        ));
    }
    if keep.is_empty() {
        return Err(format!(
            "ERROR: {} ({}) is empty",
            label,
            list_path.display()
        ));
    }
    Ok(keep)
}

pub fn preview_list(items: &[String], max_items: usize) -> String {
    let head: Vec<&str> = items.iter().take(max_items).map(|s| s.as_str()).collect();
    if items.len() > max_items {
        format!("{}...", head.join(", "))
    } else {
        head.join(", ")
    }
}

#[derive(Debug, Clone)]
pub struct Variant {
    pub orig_idx: usize,
    pub chrom: u8,
    pub id: String,
    pub is_tv: bool,
    pub raw_line: String,
}

#[derive(Debug)]
pub struct Bim {
    pub n_input: usize,
    pub variants: Vec<Variant>, // autosomes only (1..=22)
}

fn parse_autosome(chrom: &str) -> Option<u8> {
    let lower = chrom.to_ascii_lowercase();
    let stripped = lower.strip_prefix("chr").unwrap_or(&lower);
    let trimmed = stripped.trim_start_matches('0');
    let num: u8 = trimmed.parse().ok()?;
    if (1..=22).contains(&num) {
        Some(num)
    } else {
        None
    }
}

fn is_transversion(a: &str, b: &str) -> bool {
    if a.len() != 1 || b.len() != 1 {
        return false;
    }
    let ab = a.as_bytes()[0].to_ascii_uppercase();
    let bb = b.as_bytes()[0].to_ascii_uppercase();
    if ab == bb {
        return false;
    }
    let is_acgt = |c: u8| matches!(c, b'A' | b'C' | b'G' | b'T');
    if !is_acgt(ab) || !is_acgt(bb) {
        return false;
    }
    !matches!((ab, bb), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'))
}

impl Bim {
    pub fn read(path: &Path) -> Result<Self, String> {
        let file = File::open(path)
            .map_err(|e| format!("ERROR: failed to open {}: {}", path.display(), e))?;
        let reader = BufReader::new(file);

        let mut n_input = 0usize;
        let mut variants = Vec::new();
        let mut seen = HashSet::new();
        let mut dupes = Vec::new();
        let mut dupe_set = HashSet::new();

        for line in reader.lines() {
            let line = line.map_err(|e| format!("ERROR: reading {}: {}", path.display(), e))?;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            let orig_idx = n_input;
            n_input += 1;

            let Some(chrom) = parse_autosome(fields[0]) else {
                continue;
            };
            if fields.len() < 2 {
                continue;
            }
            let id = fields[1].to_string();
            if !seen.insert(id.clone()) && dupe_set.insert(id.clone()) {
                dupes.push(id.clone());
            }
            let tv = if fields.len() >= 6 {
                is_transversion(fields[4], fields[5])
            } else {
                false
            };
            variants.push(Variant {
                orig_idx,
                chrom,
                id,
                is_tv: tv,
                raw_line: line,
            });
        }

        if variants.is_empty() {
            return Err(
                "ERROR: input panel has no autosomal SNPs (chromosomes 1-22)".to_string(),
            );
        }
        if !dupes.is_empty() {
            let preview = preview_list(&dupes, 5);
            return Err(format!(
                "ERROR: {} duplicate variant ID(s) in {}: {}\n  Variant IDs must be unique (rename unnamed '.' SNPs to chr:pos, e.g. plink --set-missing-var-ids @:#).",
                dupes.len(),
                path.display(),
                preview
            ));
        }
        Ok(Self { n_input, variants })
    }
}

/// Normalized excess heterozygosity (`-F_IS = (H_obs - H_exp) / H_exp`).
#[inline]
pub fn het_excess(obs_hets: usize, obs_hom1: usize, obs_hom2: usize) -> f64 {
    let n = obs_hets + obs_hom1 + obs_hom2;
    if n == 0 {
        return 0.0;
    }
    let nf = n as f64;
    let p = (obs_hets + 2 * obs_hom2) as f64 / (2.0 * nf);
    let h_exp = 2.0 * p * (1.0 - p);
    if h_exp <= 0.0 {
        return 0.0;
    }
    let h_obs = (obs_hets as f64) / nf;
    (h_obs - h_exp) / h_exp
}

/// Exact one-sided excess-heterozygosity p-value: P(H >= obs_hets | HWE),
/// computed via the Wigginton, Cutler & Abecasis (2005) recurrence (matching
/// `bcftools` ExcHet and `plink --hwe <p> keep-fewhet`).
/// Homozygote excess (e.g. Wahlund effect in multi-population panels) yields
/// large p-values (~0.5..1.0) and is never filtered.
pub fn excess_het_p(obs_hets: usize, obs_hom1: usize, obs_hom2: usize) -> f64 {
    let obs_homr = obs_hom1.min(obs_hom2);
    let obs_homc = obs_hom1.max(obs_hom2);
    let rare_copies = 2 * obs_homr + obs_hets;
    let genotypes = obs_hets + obs_homc + obs_homr;
    if genotypes == 0 || rare_copies == 0 {
        return 1.0;
    }

    let mut het_probs = vec![0.0f64; rare_copies + 1];
    let mut mid = ((rare_copies as u64) * (2 * (genotypes as u64) - (rare_copies as u64))
        / (2 * (genotypes as u64))) as usize;
    if (rare_copies & 1) != (mid & 1) {
        mid += 1;
    }
    if mid > rare_copies {
        mid -= 2;
    }

    het_probs[mid] = 1.0;
    let mut sum = 1.0f64;

    // Recurse downward from mid
    let mut curr_hets = mid;
    let mut curr_homr = (rare_copies - mid) / 2;
    let mut curr_homc = genotypes - curr_hets - curr_homr;
    while curr_hets >= 2 {
        let next = het_probs[curr_hets] * (curr_hets as f64) * ((curr_hets - 1) as f64)
            / (4.0 * ((curr_homr + 1) as f64) * ((curr_homc + 1) as f64));
        het_probs[curr_hets - 2] = next;
        sum += next;
        curr_hets -= 2;
        curr_homr += 1;
        curr_homc += 1;
    }

    // Recurse upward from mid
    curr_hets = mid;
    curr_homr = (rare_copies - mid) / 2;
    curr_homc = genotypes - curr_hets - curr_homr;
    while curr_hets + 2 <= rare_copies {
        let next = het_probs[curr_hets] * 4.0 * (curr_homr as f64) * (curr_homc as f64)
            / (((curr_hets + 2) as f64) * ((curr_hets + 1) as f64));
        het_probs[curr_hets + 2] = next;
        sum += next;
        curr_hets += 2;
        curr_homr -= 1;
        curr_homc -= 1;
    }

    // Upper tail: P(H >= obs_hets)
    let p_upper: f64 = het_probs[obs_hets..].iter().sum();
    (p_upper / sum).min(1.0)
}

#[derive(Debug, Clone, Copy, Default)]
pub struct QcFilters {
    pub max_geno: Option<f64>,
    pub min_maf: Option<f64>,
    pub min_excess_het_p: Option<f64>,
    /// Optional effect-size floor on `-F_IS = (H_obs - H_exp)/H_exp` (used by `--excess-het auto`).
    pub min_excess_het_fis: Option<f64>,
}

impl QcFilters {
    pub fn is_active(&self) -> bool {
        self.max_geno.is_some() || self.min_maf.is_some() || self.min_excess_het_p.is_some()
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct QcStats {
    pub dropped_geno: usize,
    pub dropped_maf: usize,
    pub dropped_excess_het: usize,
}

impl QcStats {
    pub fn total_dropped(&self) -> usize {
        self.dropped_geno + self.dropped_maf + self.dropped_excess_het
    }
}

/// Compact bitset of priority sample calls per autosomal SNP.
pub struct PriorityCalls {
    pub n_prio: usize,
    pub n_words: usize,
    pub coverage: Vec<usize>,
    /// Flat `n_snps * n_words` bitmask: bit `i` is 1 iff priority sample `i` has a call.
    pub bits: Vec<u64>,
}

impl PriorityCalls {
    pub fn count_calls_on(&self, snp_indices: &[usize]) -> Vec<usize> {
        let mut counts = vec![0usize; self.n_prio];
        for &k in snp_indices {
            let base = k * self.n_words;
            for w in 0..self.n_words {
                let mut word = self.bits[base + w];
                while word != 0 {
                    let tz = word.trailing_zeros() as usize;
                    counts[w * 64 + tz] += 1;
                    word &= word - 1;
                }
            }
        }
        counts
    }
}

/// Single-pass streaming reader over `.bed`:
/// applies optional LD-sample QC filters (`--geno`, `--maf`, `--hwe`) and extracts
/// priority-sample call bitsets and LD-sample 2-bit genotype bitvectors for passing SNPs.
pub fn scan_bed(
    bed_path: &Path,
    n_fam: usize,
    n_input: usize,
    variants: &[Variant],
    prio_indices: &[usize],
    ld_indices: &[usize],
    qc: QcFilters,
) -> Result<(Vec<Variant>, PriorityCalls, LdMatrix, QcStats), String> {
    let mut file = File::open(bed_path)
        .map_err(|e| format!("ERROR: failed to open {}: {}", bed_path.display(), e))?;
    let meta = file
        .metadata()
        .map_err(|e| format!("ERROR: stat {}: {}", bed_path.display(), e))?;

    let bytes_per_snp = (n_fam + 3) / 4;
    let expected_len = 3u64 + (n_input as u64) * (bytes_per_snp as u64);
    if meta.len() != expected_len {
        return Err(format!(
            "ERROR: {} has {} bytes, expected {} (3 + {} variants * {} bytes/variant)",
            bed_path.display(),
            meta.len(),
            expected_len,
            n_input,
            bytes_per_snp
        ));
    }

    let mut magic = [0u8; 3];
    file.read_exact(&mut magic)
        .map_err(|e| format!("ERROR: reading {}: {}", bed_path.display(), e))?;
    if magic != [0x6c, 0x1b, 0x01] {
        return Err(format!(
            "ERROR: {} is not a SNP-major PLINK 1 .bed",
            bed_path.display()
        ));
    }

    let mut reader = BufReader::with_capacity(1 << 20, file);

    let n_snps_in = variants.len();
    let n_prio = prio_indices.len();
    let prio_words = (n_prio + 63) / 64;
    let mut prio_coverage = vec![0usize; n_prio];
    let mut prio_bits: Vec<u64> = Vec::with_capacity(n_snps_in * prio_words);

    // Precompute (byte_idx, shift, word_idx, bit_mask) for priority and LD samples
    let prio_loc: Vec<(usize, u8, usize, u64)> = prio_indices
        .iter()
        .enumerate()
        .map(|(i, &fam_idx)| (fam_idx / 4, ((fam_idx % 4) * 2) as u8, i / 64, 1u64 << (i % 64)))
        .collect();

    let n_ld = ld_indices.len();
    let ld_words = (n_ld + 63) / 64;
    let mut ld_valid: Vec<u64> = Vec::with_capacity(n_snps_in * ld_words);
    let mut ld_het: Vec<u64> = Vec::with_capacity(n_snps_in * ld_words);
    let mut ld_hom2: Vec<u64> = Vec::with_capacity(n_snps_in * ld_words);
    let mut is_constant: Vec<bool> = Vec::with_capacity(n_snps_in);
    let mut kept_variants: Vec<Variant> = Vec::with_capacity(n_snps_in);
    let mut qc_stats = QcStats::default();

    let ld_loc: Vec<(usize, u8, usize, u64)> = ld_indices
        .iter()
        .enumerate()
        .map(|(i, &fam_idx)| (fam_idx / 4, ((fam_idx % 4) * 2) as u8, i / 64, 1u64 << (i % 64)))
        .collect();

    let mut row = vec![0u8; bytes_per_snp];
    let mut cur_orig = 0usize;
    let mut row_valid = vec![0u64; ld_words];
    let mut row_het = vec![0u64; ld_words];
    let mut row_hom2 = vec![0u64; ld_words];
    let mut row_prio = vec![0u64; prio_words];

    for var in variants {
        while cur_orig < var.orig_idx {
            reader
                .read_exact(&mut row)
                .map_err(|e| format!("ERROR: reading {}: {}", bed_path.display(), e))?;
            cur_orig += 1;
        }
        reader
            .read_exact(&mut row)
            .map_err(|e| format!("ERROR: reading {}: {}", bed_path.display(), e))?;
        cur_orig += 1;

        // Decode LD genotypes first so QC filters can reject before touching priority counts
        row_valid.fill(0);
        row_het.fill(0);
        row_hom2.fill(0);
        let mut n_val = 0u64;
        let mut n_het = 0u64;
        let mut n_hom2 = 0u64;
        for &(bidx, shift, w, mask) in &ld_loc {
            match (row[bidx] >> shift) & 3 {
                0 => {
                    row_valid[w] |= mask;
                    n_val += 1;
                }
                1 => {}
                2 => {
                    row_valid[w] |= mask;
                    row_het[w] |= mask;
                    n_val += 1;
                    n_het += 1;
                }
                3 => {
                    row_valid[w] |= mask;
                    row_hom2[w] |= mask;
                    n_val += 1;
                    n_hom2 += 1;
                }
                _ => unreachable!(),
            }
        }

        // Apply LD-sample QC filters (geno -> maf -> hwe, matching PLINK 1.9 order)
        if let Some(max_geno) = qc.max_geno {
            let miss_rate = ((n_ld as u64 - n_val) as f64) / (n_ld as f64);
            if miss_rate > max_geno + 1e-12 {
                qc_stats.dropped_geno += 1;
                continue;
            }
        }
        if let Some(min_maf) = qc.min_maf {
            let maf = if n_val == 0 {
                0.0
            } else {
                let p = ((n_het + 2 * n_hom2) as f64) / (2.0 * (n_val as f64));
                p.min(1.0 - p)
            };
            if maf < min_maf - 1e-12 {
                qc_stats.dropped_maf += 1;
                continue;
            }
        }
        if let Some(min_p) = qc.min_excess_het_p {
            let n_hom1 = (n_val - n_het - n_hom2) as usize;
            let passes_fis_floor = match qc.min_excess_het_fis {
                None => true,
                Some(min_ex) => het_excess(n_het as usize, n_hom1, n_hom2 as usize) > min_ex,
            };
            if passes_fis_floor && excess_het_p(n_het as usize, n_hom1, n_hom2 as usize) < min_p {
                qc_stats.dropped_excess_het += 1;
                continue;
            }
        }

        // Passed QC: record variant, LD bitvectors, and priority calls
        kept_variants.push(var.clone());
        ld_valid.extend_from_slice(&row_valid);
        ld_het.extend_from_slice(&row_het);
        ld_hom2.extend_from_slice(&row_hom2);

        let sum_x = n_het + 2 * n_hom2;
        let sum_x2 = n_het + 4 * n_hom2;
        is_constant.push(n_val == 0 || n_val * sum_x2 == sum_x * sum_x);

        row_prio.fill(0);
        for (i, &(bidx, shift, w, mask)) in prio_loc.iter().enumerate() {
            let code = (row[bidx] >> shift) & 3;
            if code != 1 {
                prio_coverage[i] += 1;
                row_prio[w] |= mask;
            }
        }
        prio_bits.extend_from_slice(&row_prio);
    }

    let prio_calls = PriorityCalls {
        n_prio,
        n_words: prio_words,
        coverage: prio_coverage,
        bits: prio_bits,
    };
    let ld_matrix = LdMatrix {
        n_words: ld_words,
        valid: ld_valid,
        het: ld_het,
        hom2: ld_hom2,
        is_constant,
    };
    Ok((kept_variants, prio_calls, ld_matrix, qc_stats))
}

/// Write pruned `.bed`, `.bim`, `.fam` directly by copying raw `.bed` SNP byte slices.
pub fn write_bed_subset(
    bfile_prefix: &str,
    out_prefix: &str,
    n_fam: usize,
    variants: &[Variant],
    kept_auto_indices: &[usize],
) -> Result<(), String> {
    // 1. Copy .fam verbatim
    let in_fam = format!("{}.fam", bfile_prefix);
    let out_fam = format!("{}.fam", out_prefix);
    fs::copy(&in_fam, &out_fam)
        .map_err(|e| format!("ERROR: failed to copy {} to {}: {}", in_fam, out_fam, e))?;

    // 2. Write kept .bim lines
    let out_bim = format!("{}.bim", out_prefix);
    {
        let f = File::create(&out_bim)
            .map_err(|e| format!("ERROR: failed to create {}: {}", out_bim, e))?;
        let mut w = BufWriter::new(f);
        for &k in kept_auto_indices {
            writeln!(w, "{}", variants[k].raw_line)
                .map_err(|e| format!("ERROR: writing {}: {}", out_bim, e))?;
        }
    }

    // 3. Copy raw .bed byte slices for kept variants
    let in_bed = format!("{}.bed", bfile_prefix);
    let out_bed = format!("{}.bed", out_prefix);
    let mut reader = File::open(&in_bed)
        .map_err(|e| format!("ERROR: failed to open {}: {}", in_bed, e))?;
    let out_file = File::create(&out_bed)
        .map_err(|e| format!("ERROR: failed to create {}: {}", out_bed, e))?;
    let mut writer = BufWriter::with_capacity(1 << 20, out_file);

    writer
        .write_all(&[0x6c, 0x1b, 0x01])
        .map_err(|e| format!("ERROR: writing {}: {}", out_bed, e))?;

    let bytes_per_snp = (n_fam + 3) / 4;
    let mut row = vec![0u8; bytes_per_snp];
    let mut cur_orig: Option<usize> = None;

    for &k in kept_auto_indices {
        let target = variants[k].orig_idx;
        let need_seek = match cur_orig {
            Some(c) if c == target => false,
            _ => true,
        };
        if need_seek {
            let offset = 3u64 + (target as u64) * (bytes_per_snp as u64);
            reader
                .seek(SeekFrom::Start(offset))
                .map_err(|e| format!("ERROR: seeking {}: {}", in_bed, e))?;
        }
        reader
            .read_exact(&mut row)
            .map_err(|e| format!("ERROR: reading {}: {}", in_bed, e))?;
        cur_orig = Some(target + 1);
        writer
            .write_all(&row)
            .map_err(|e| format!("ERROR: writing {}: {}", out_bed, e))?;
    }
    writer
        .flush()
        .map_err(|e| format!("ERROR: flushing {}: {}", out_bed, e))?;

    Ok(())
}
