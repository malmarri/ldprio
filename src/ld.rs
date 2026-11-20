use rayon::prelude::*;

use crate::plink::Variant;

pub struct LdMatrix {
    pub n_words: usize,
    pub valid: Vec<u64>,
    pub het: Vec<u64>,
    pub hom2: Vec<u64>,
    pub is_constant: Vec<bool>,
}

impl LdMatrix {
    /// Exact unphased genotype Pearson r^2 >= threshold over jointly non-missing LD samples.
    #[inline]
    pub fn pair_r2_ge(&self, a: usize, b: usize, r2_thresh: f64) -> bool {
        if self.is_constant[a] || self.is_constant[b] {
            return false;
        }
        let nw = self.n_words;
        let oa = a * nw;
        let ob = b * nw;

        let va = &self.valid[oa..oa + nw];
        let ha = &self.het[oa..oa + nw];
        let ta = &self.hom2[oa..oa + nw];

        let vb = &self.valid[ob..ob + nw];
        let hb = &self.het[ob..ob + nw];
        let tb = &self.hom2[ob..ob + nw];

        let mut n = 0u64;
        let mut s_ha = 0u64;
        let mut s_ta = 0u64;
        let mut s_hb = 0u64;
        let mut s_tb = 0u64;
        let mut s_hh = 0u64;
        let mut s_ht = 0u64;
        let mut s_th = 0u64;
        let mut s_tt = 0u64;

        for w in 0..nw {
            let v = va[w] & vb[w];
            if v == 0 {
                continue;
            }
            n += v.count_ones() as u64;
            let h1 = ha[w] & v;
            let t1 = ta[w] & v;
            let h2 = hb[w] & v;
            let t2 = tb[w] & v;

            s_ha += h1.count_ones() as u64;
            s_ta += t1.count_ones() as u64;
            s_hb += h2.count_ones() as u64;
            s_tb += t2.count_ones() as u64;

            s_hh += (h1 & h2).count_ones() as u64;
            s_ht += (h1 & t2).count_ones() as u64;
            s_th += (t1 & h2).count_ones() as u64;
            s_tt += (t1 & t2).count_ones() as u64;
        }

        if n == 0 {
            return false;
        }
        let nf = n as f64;
        let sx = (s_ha + 2 * s_ta) as f64;
        let sx2 = (s_ha + 4 * s_ta) as f64;
        let sy = (s_hb + 2 * s_tb) as f64;
        let sy2 = (s_hb + 4 * s_tb) as f64;
        let sxy = (s_hh + 2 * (s_ht + s_th) + 4 * s_tt) as f64;

        let var_x = nf * sx2 - sx * sx;
        let var_y = nf * sy2 - sy * sy;
        if var_x <= 0.0 || var_y <= 0.0 {
            return false;
        }
        let num = nf * sxy - sx * sy;
        num * num >= r2_thresh * var_x * var_y
    }
}

/// Build the sparse LD conflict graph over `surviving` (autosomal SNP indices)
/// within `window` variants on the same chromosome, parallelized with rayon.
/// Returns `(adj, n_pairs)` where `adj` is indexed by `0..surviving.len()`.
pub fn build_ld_graph(
    ld: &LdMatrix,
    surviving: &[usize],
    variants: &[Variant],
    window: usize,
    r2_thresh: f64,
) -> (Vec<Vec<usize>>, usize) {
    let m = surviving.len();
    if m < 2 || window < 2 {
        return (vec![Vec::new(); m], 0);
    }

    // Group local positions 0..m by chromosome (preserving genomic order)
    let mut by_chrom: [Vec<usize>; 23] = std::array::from_fn(|_| Vec::new());
    for (local_idx, &auto_idx) in surviving.iter().enumerate() {
        by_chrom[variants[auto_idx].chrom as usize].push(local_idx);
    }

    // Flatten into parallel work units: chunks of 256 anchor SNPs per chromosome
    const CHUNK: usize = 256;
    let mut tasks: Vec<(usize, usize, usize)> = Vec::new();
    for c in 1..=22 {
        let len = by_chrom[c].len();
        if len < 2 {
            continue;
        }
        let mut start = 0;
        while start + 1 < len {
            let end = (start + CHUNK).min(len - 1);
            tasks.push((c, start, end));
            start = end;
        }
    }

    let edge_chunks: Vec<Vec<(usize, usize)>> = tasks
        .par_iter()
        .map(|&(c, start, end)| {
            let chr_pos = &by_chrom[c];
            let chr_len = chr_pos.len();
            let mut edges = Vec::new();
            for i in start..end {
                let u = chr_pos[i];
                let snp_u = surviving[u];
                if ld.is_constant[snp_u] {
                    continue;
                }
                let j_end = (i + window).min(chr_len);
                for j in (i + 1)..j_end {
                    let v = chr_pos[j];
                    let snp_v = surviving[v];
                    if ld.pair_r2_ge(snp_u, snp_v, r2_thresh) {
                        edges.push((u, v));
                    }
                }
            }
            edges
        })
        .collect();

    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut n_pairs = 0usize;
    for chunk in edge_chunks {
        n_pairs += chunk.len();
        for (u, v) in chunk {
            adj[u].push(v);
            adj[v].push(u);
        }
    }
    (adj, n_pairs)
}

/// Find all dropped autosomal SNPs that have zero LD partners (>= `r2_thresh`)
/// within `±(window - 1)` variants of their insertion position in `kept` on the
/// same chromosome. Runs in parallel across chromosomes.
pub fn find_backfill_candidates(
    ld: &LdMatrix,
    kept: &[usize],
    variants: &[Variant],
    window: usize,
    r2_thresh: f64,
) -> Vec<usize> {
    let n_snps = variants.len();
    if window < 2 || kept.len() == n_snps {
        return Vec::new();
    }

    let mut is_kept = vec![false; n_snps];
    for &k in kept {
        is_kept[k] = true;
    }

    // Group kept and dropped SNPs by chromosome (both sorted by auto_idx = genomic order)
    let mut kept_by_chr: [Vec<usize>; 23] = std::array::from_fn(|_| Vec::new());
    let mut drop_by_chr: [Vec<usize>; 23] = std::array::from_fn(|_| Vec::new());
    for (k, v) in variants.iter().enumerate() {
        if is_kept[k] {
            kept_by_chr[v.chrom as usize].push(k);
        } else {
            drop_by_chr[v.chrom as usize].push(k);
        }
    }

    let reach = window - 1;
    let per_chr: Vec<Vec<usize>> = (1..=22usize)
        .into_par_iter()
        .map(|c| {
            let kchr = &kept_by_chr[c];
            let dchr = &drop_by_chr[c];
            if dchr.is_empty() {
                return Vec::new();
            }
            if kchr.is_empty() {
                return dchr.clone();
            }
            let mut eligible = Vec::new();
            let mut pos = 0usize; // first index in kchr with kchr[pos] > d
            for &d in dchr {
                while pos < kchr.len() && kchr[pos] < d {
                    pos += 1;
                }
                if ld.is_constant[d] {
                    eligible.push(d);
                    continue;
                }
                let lo = pos.saturating_sub(reach);
                let hi = (pos + reach).min(kchr.len());
                let mut conflict = false;
                for &k in &kchr[lo..hi] {
                    if ld.pair_r2_ge(d, k, r2_thresh) {
                        conflict = true;
                        break;
                    }
                }
                if !conflict {
                    eligible.push(d);
                }
            }
            eligible
        })
        .collect();

    per_chr.into_iter().flatten().collect()
}
