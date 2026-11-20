use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use crate::plink::{PriorityCalls, Variant};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Weighting {
    Inverse,
    Fair,
}

impl Weighting {
    pub fn as_str(self) -> &'static str {
        match self {
            Weighting::Inverse => "inverse",
            Weighting::Fair => "fair",
        }
    }
}

pub enum PriorityPool {
    /// Weight per autosomal SNP (0.0 if not in priority pool).
    Inverse(Vec<f64>),
    /// Covering priority sample indices per autosomal SNP (empty if not in pool).
    Fair(Vec<Vec<u32>>),
}

impl PriorityPool {
    pub fn build(calls: &PriorityCalls, n_snps: usize, weighting: Weighting) -> Self {
        match weighting {
            Weighting::Inverse => {
                let inv: Vec<f64> = calls
                    .coverage
                    .iter()
                    .map(|&c| if c > 0 { 1.0 / (c as f64) } else { 0.0 })
                    .collect();
                let mut weights = vec![0.0f64; n_snps];
                for k in 0..n_snps {
                    let base = k * calls.n_words;
                    let mut sum = 0.0f64;
                    for w in 0..calls.n_words {
                        let mut word = calls.bits[base + w];
                        while word != 0 {
                            let tz = word.trailing_zeros() as usize;
                            sum += inv[w * 64 + tz];
                            word &= word - 1;
                        }
                    }
                    weights[k] = sum;
                }
                PriorityPool::Inverse(weights)
            }
            Weighting::Fair => {
                let mut cover = vec![Vec::new(); n_snps];
                for k in 0..n_snps {
                    let base = k * calls.n_words;
                    let mut samples = Vec::new();
                    for w in 0..calls.n_words {
                        let mut word = calls.bits[base + w];
                        while word != 0 {
                            let tz = word.trailing_zeros() as usize;
                            samples.push((w * 64 + tz) as u32);
                            word &= word - 1;
                        }
                    }
                    cover[k] = samples;
                }
                PriorityPool::Fair(cover)
            }
        }
    }

    #[inline]
    pub fn in_pool(&self, auto_idx: usize) -> bool {
        match self {
            PriorityPool::Inverse(w) => w[auto_idx] > 0.0,
            PriorityPool::Fair(c) => !c[auto_idx].is_empty(),
        }
    }

    pub fn pool_size(&self) -> usize {
        match self {
            PriorityPool::Inverse(w) => w.iter().filter(|&&x| x > 0.0).count(),
            PriorityPool::Fair(c) => c.iter().filter(|v| !v.is_empty()).count(),
        }
    }
}

/// Quantize inverse weight to 12 decimal places so exact mathematical ties
/// are not broken by 1-ULP floating-point rounding before `--prefer-transversions`.
#[inline]
fn quantize_weight(w: f64) -> i64 {
    (w * 1e12).round() as i64
}

/// Maximal independent set over `adj` (indexed `0..surviving.len()`):
/// priority-pool SNPs go first (`inverse` or `fair`), then the rest in genomic order
/// (transversions first when `prefer_tv` is enabled). Returns kept autosomal SNP indices
/// in genomic order.
pub fn greedy_select(
    surviving: &[usize],
    adj: &[Vec<usize>],
    priority: &PriorityPool,
    variants: &[Variant],
    prefer_tv: bool,
) -> Vec<usize> {
    let m = surviving.len();
    let mut blocked = vec![false; m];
    let mut kept = vec![false; m];

    let keep = |u: usize, kept: &mut [bool], blocked: &mut [bool]| {
        kept[u] = true;
        for &nb in &adj[u] {
            blocked[nb] = true;
        }
    };

    match priority {
        PriorityPool::Fair(cover) => {
            // Find max sample index
            let mut max_sample = 0usize;
            for &auto_idx in surviving {
                for &s in &cover[auto_idx] {
                    if (s as usize) + 1 > max_sample {
                        max_sample = (s as usize) + 1;
                    }
                }
            }
            let mut queues: Vec<Vec<usize>> = vec![Vec::new(); max_sample];
            for (local_idx, &auto_idx) in surviving.iter().enumerate() {
                for &s in &cover[auto_idx] {
                    queues[s as usize].push(local_idx);
                }
            }
            for q in &mut queues {
                if q.is_empty() {
                    continue;
                }
                // Popped from the end: most-shared first, then transversions, then genomic order
                q.sort_unstable_by(|&a, &b| {
                    let sa = surviving[a];
                    let sb = surviving[b];
                    cover[sa]
                        .len()
                        .cmp(&cover[sb].len())
                        .then_with(|| {
                            let ta = prefer_tv && variants[sa].is_tv;
                            let tb = prefer_tv && variants[sb].is_tv;
                            ta.cmp(&tb)
                        })
                        .then_with(|| b.cmp(&a)) // smaller local_idx goes to end
                });
            }

            let mut n_kept = vec![0usize; max_sample];
            let mut turns: BinaryHeap<Reverse<(usize, usize, usize)>> = BinaryHeap::new();
            for (i, q) in queues.iter().enumerate() {
                if !q.is_empty() {
                    turns.push(Reverse((0, q.len(), i)));
                }
            }

            while let Some(Reverse((n, size, i))) = turns.pop() {
                if n != n_kept[i] {
                    turns.push(Reverse((n_kept[i], size, i)));
                    continue;
                }
                let q = &mut queues[i];
                while let Some(&u) = q.last() {
                    if blocked[u] || kept[u] {
                        q.pop();
                    } else {
                        break;
                    }
                }
                let Some(u) = q.pop() else {
                    continue;
                };
                keep(u, &mut kept, &mut blocked);
                for &j in &cover[surviving[u]] {
                    n_kept[j as usize] += 1;
                }
                turns.push(Reverse((n_kept[i], size, i)));
            }
        }
        PriorityPool::Inverse(weights) => {
            let mut pool: Vec<usize> = (0..m)
                .filter(|&u| weights[surviving[u]] > 0.0)
                .collect();
            pool.sort_by(|&a, &b| {
                let sa = surviving[a];
                let sb = surviving[b];
                let wa = quantize_weight(weights[sa]);
                let wb = quantize_weight(weights[sb]);
                wb.cmp(&wa)
                    .then_with(|| {
                        let ta = prefer_tv && variants[sa].is_tv;
                        let tb = prefer_tv && variants[sb].is_tv;
                        tb.cmp(&ta)
                    })
                    .then_with(|| {
                        // Preserve exact float ordering if quantized values tied only by huge scale
                        weights[sb]
                            .partial_cmp(&weights[sa])
                            .unwrap_or(Ordering::Equal)
                            .then_with(|| a.cmp(&b))
                    })
            });
            for u in pool {
                if !blocked[u] {
                    keep(u, &mut kept, &mut blocked);
                }
            }
        }
    }

    // Non-priority SNPs
    let mut rest: Vec<usize> = (0..m)
        .filter(|&u| !priority.in_pool(surviving[u]))
        .collect();
    if prefer_tv {
        rest.sort_by_key(|&u| !variants[surviving[u]].is_tv);
    }
    for u in rest {
        if !blocked[u] {
            keep(u, &mut kept, &mut blocked);
        }
    }

    (0..m)
        .filter(|&u| kept[u])
        .map(|u| surviving[u])
        .collect()
}
