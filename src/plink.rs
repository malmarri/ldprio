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
/// extracts priority-sample call bitsets and LD-sample 2-bit genotype bitvectors.
pub fn scan_bed(
    bed_path: &Path,
    n_fam: usize,
    n_input: usize,
    variants: &[Variant],
    prio_indices: &[usize],
    ld_indices: &[usize],
) -> Result<(PriorityCalls, LdMatrix), String> {
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

    let n_snps = variants.len();
    let n_prio = prio_indices.len();
    let prio_words = (n_prio + 63) / 64;
    let mut prio_coverage = vec![0usize; n_prio];
    let mut prio_bits = vec![0u64; n_snps * prio_words];

    // Precompute (byte_idx, shift, word_idx, bit_mask) for priority and LD samples
    let prio_loc: Vec<(usize, u8, usize, u64)> = prio_indices
        .iter()
        .enumerate()
        .map(|(i, &fam_idx)| (fam_idx / 4, ((fam_idx % 4) * 2) as u8, i / 64, 1u64 << (i % 64)))
        .collect();

    let n_ld = ld_indices.len();
    let ld_words = (n_ld + 63) / 64;
    let mut ld_valid = vec![0u64; n_snps * ld_words];
    let mut ld_het = vec![0u64; n_snps * ld_words];
    let mut ld_hom2 = vec![0u64; n_snps * ld_words];
    let mut is_constant = vec![false; n_snps];

    let ld_loc: Vec<(usize, u8, usize, u64)> = ld_indices
        .iter()
        .enumerate()
        .map(|(i, &fam_idx)| (fam_idx / 4, ((fam_idx % 4) * 2) as u8, i / 64, 1u64 << (i % 64)))
        .collect();

    let mut row = vec![0u8; bytes_per_snp];
    let mut cur_orig = 0usize;

    for (k, var) in variants.iter().enumerate() {
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

        // Priority calls
        let p_base = k * prio_words;
        for (i, &(bidx, shift, w, mask)) in prio_loc.iter().enumerate() {
            let code = (row[bidx] >> shift) & 3;
            if code != 1 {
                prio_coverage[i] += 1;
                prio_bits[p_base + w] |= mask;
            }
        }

        // LD genotypes
        let l_base = k * ld_words;
        let mut n_val = 0u64;
        let mut n_het = 0u64;
        let mut n_hom2 = 0u64;
        for &(bidx, shift, w, mask) in &ld_loc {
            match (row[bidx] >> shift) & 3 {
                0 => {
                    // hom A1 (dosage 0)
                    ld_valid[l_base + w] |= mask;
                    n_val += 1;
                }
                1 => {
                    // missing
                }
                2 => {
                    // het (dosage 1)
                    ld_valid[l_base + w] |= mask;
                    ld_het[l_base + w] |= mask;
                    n_val += 1;
                    n_het += 1;
                }
                3 => {
                    // hom A2 (dosage 2)
                    ld_valid[l_base + w] |= mask;
                    ld_hom2[l_base + w] |= mask;
                    n_val += 1;
                    n_hom2 += 1;
                }
                _ => unreachable!(),
            }
        }
        let sum_x = n_het + 2 * n_hom2;
        let sum_x2 = n_het + 4 * n_hom2;
        let constant = n_val == 0 || n_val * sum_x2 == sum_x * sum_x;
        is_constant[k] = constant;
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
    Ok((prio_calls, ld_matrix))
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
