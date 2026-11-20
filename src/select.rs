use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use crate::plink::{PriorityCalls, Variant};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Weighting {
    Adaptive,
    Inverse,
    Fair,
}

impl Weighting {
    pub fn as_str(self) -> &'static str {
        match self {
            Weighting::Adaptive => "adaptive",
            Weighting::Inverse => "inverse",
            Weighting::Fair => "fair",
        }
    }
}

pub enum PriorityPool {
    /// Coverage-normalized dynamic proportional fairness + GWMIN2 graph awareness.
    Adaptive {
        cover: Vec<Vec<u32>>,
        inv_cov: Vec<f64>,
    },
    /// Static weight per autosomal SNP (0.0 if not in priority pool).
    Inverse(Vec<f64>),
    /// Covering priority sample indices per autosomal SNP (empty if not in pool).
    Fair(Vec<Vec<u32>>),
}

fn build_cover(calls: &PriorityCalls, n_snps: usize) -> Vec<Vec<u32>> {
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
    cover
}

impl PriorityPool {
    pub fn build(calls: &PriorityCalls, n_snps: usize, weighting: Weighting) -> Self {
        match weighting {
            Weighting::Adaptive => {
                let inv_cov: Vec<f64> = calls
                    .coverage
                    .iter()
                    .map(|&c| if c > 0 { 1.0 / (c as f64) } else { 0.0 })
                    .collect();
                PriorityPool::Adaptive {
                    cover: build_cover(calls, n_snps),
                    inv_cov,
                }
            }
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
            Weighting::Fair => PriorityPool::Fair(build_cover(calls, n_snps)),
        }
    }

    #[inline]
    pub fn in_pool(&self, auto_idx: usize) -> bool {
        match self {
            PriorityPool::Adaptive { cover, .. } => !cover[auto_idx].is_empty(),
            PriorityPool::Inverse(w) => w[auto_idx] > 0.0,
            PriorityPool::Fair(c) => !c[auto_idx].is_empty(),
        }
    }

    pub fn pool_size(&self) -> usize {
        match self {
            PriorityPool::Adaptive { cover, .. } => cover.iter().filter(|v| !v.is_empty()).count(),
            PriorityPool::Inverse(w) => w.iter().filter(|&&x| x > 0.0).count(),
            PriorityPool::Fair(c) => c.iter().filter(|v| !v.is_empty()).count(),
        }
    }
}

/// Quantize float score to 12 decimal places so exact mathematical ties
/// are not broken by 1-ULP floating-point rounding before `--prefer-transversions`.
#[inline]
fn quantize_weight(w: f64) -> i64 {
    (w * 1e12).round() as i64
}

#[derive(Clone, Copy)]
struct AdaptiveItem {
    q_ratio: i64,
    is_tv: bool,
    prio_deg: usize,
    idx: usize,
    ratio: f64,
    ver: u32,
}

impl PartialEq for AdaptiveItem {
    fn eq(&self, other: &Self) -> bool {
        self.q_ratio == other.q_ratio
            && self.is_tv == other.is_tv
            && self.prio_deg == other.prio_deg
            && self.idx == other.idx
    }
}
impl Eq for AdaptiveItem {}
impl PartialOrd for AdaptiveItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for AdaptiveItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.q_ratio
            .cmp(&other.q_ratio)
            .then_with(|| self.is_tv.cmp(&other.is_tv))
            .then_with(|| other.prio_deg.cmp(&self.prio_deg)) // lower priority-degree preferred
            .then_with(|| other.idx.cmp(&self.idx)) // earlier genomic order preferred
    }
}

/// Maximal independent set over `adj` (indexed `0..surviving.len()`):
/// priority-pool SNPs go first (`adaptive`, `inverse`, or `fair`), then the rest
/// in genomic order (transversions first when `prefer_tv` is enabled). Returns kept
/// autosomal SNP indices in genomic order.
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
        PriorityPool::Adaptive { cover, inv_cov } => {
            let n_prio = inv_cov.len();
            let mut k_kept = vec![0usize; n_prio];
            let mut version = vec![0u32; m];

            let snp_weight = |u: usize, k_kept: &[usize]| -> f64 {
                cover[surviving[u]]
                    .iter()
                    .map(|&s| inv_cov[s as usize] / ((k_kept[s as usize] + 1) as f64))
                    .sum()
            };

            let eval = |u: usize, k_kept: &[usize], blocked: &[bool]| -> (f64, usize) {
                let w_u = snp_weight(u, k_kept);
                let mut nb_w = 0.0f64;
                let mut deg = 0usize;
                for &v in &adj[u] {
                    if !blocked[v] && !cover[surviving[v]].is_empty() {
                        nb_w += snp_weight(v, k_kept);
                        deg += 1;
                    }
                }
                (w_u / (w_u + nb_w), deg)
            };

            let mut heap = BinaryHeap::new();
            for u in 0..m {
                if !cover[surviving[u]].is_empty() {
                    let (ratio, deg) = eval(u, &k_kept, &blocked);
                    heap.push(AdaptiveItem {
                        q_ratio: quantize_weight(ratio),
                        is_tv: prefer_tv && variants[surviving[u]].is_tv,
                        prio_deg: deg,
                        idx: u,
                        ratio,
                        ver: 0,
                    });
                }
            }

            while let Some(top) = heap.pop() {
                let u = top.idx;
                if blocked[u] || kept[u] || top.ver != version[u] {
                    continue;
                }
                let (true_ratio, true_deg) = eval(u, &k_kept, &blocked);
                let true_q = quantize_weight(true_ratio);
                if (true_ratio - top.ratio).abs() > 1e-12 || true_deg != top.prio_deg {
                    version[u] += 1;
                    heap.push(AdaptiveItem {
                        q_ratio: true_q,
                        is_tv: top.is_tv,
                        prio_deg: true_deg,
                        idx: u,
                        ratio: true_ratio,
                        ver: version[u],
                    });
                    continue;
                }

                kept[u] = true;
                let mut newly_blocked = Vec::new();
                for &nb in &adj[u] {
                    if !blocked[nb] {
                        blocked[nb] = true;
                        newly_blocked.push(nb);
                    }
                }
                for &s in &cover[surviving[u]] {
                    k_kept[s as usize] += 1;
                }
                // Blocking `nb` can increase the GWMIN2 ratio of `nb`'s unblocked priority neighbors
                for &nb in &newly_blocked {
                    for &w in &adj[nb] {
                        if !blocked[w] && !kept[w] && !cover[surviving[w]].is_empty() {
                            version[w] += 1;
                            let (ratio, deg) = eval(w, &k_kept, &blocked);
                            heap.push(AdaptiveItem {
                                q_ratio: quantize_weight(ratio),
                                is_tv: prefer_tv && variants[surviving[w]].is_tv,
                                prio_deg: deg,
                                idx: w,
                                ratio,
                                ver: version[w],
                            });
                        }
                    }
                }
            }
        }
        PriorityPool::Fair(cover) => {
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
                        .then_with(|| b.cmp(&a))
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
