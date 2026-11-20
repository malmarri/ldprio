# ldprio

LD pruning that prioritizes SNPs in low-coverage samples — standalone Rust binary, no PLINK or Python required.

## The problem

LD pruning is commonly performed for PCA and model-based clustering analysis like ADMIXTURE,
which requires relatively independent SNPs. The most used method for this is `plink --indep-pairwise`.
When a block contains several mutually redundant SNPs, PLINK decides which
one survives without looking at which samples have data there — and for a
low-coverage sample, covering a sparse number of sites, the survivor is
usually not one of the few sites it has data for, resulting in an unnecessary loss for ultra-low coverage samples.

## The fix

Do the pruning as a single genome-wide pass, so every LD relationship is
still evaluated exactly once — but choose the SNP within each block toward
sites covered by the low-coverage samples you nominate, and settle
blocks contested between them by `--weighting`.

Well-covered samples are unaffected: they had an equivalent marker either
way. Sparse samples keep substantially more of their real data.

<p align="center">
  <img src="docs/ld-block.svg" alt="One LD block of five redundant SNPs: high coverage samples have calls at all five, the ancient sample only at s3. plink keeps s1, leaving the ancient sample with no SNP from this block; ldprio keeps s3, so the ancient sample keeps its site. High coverage samples keep one SNP either way." width="860">
</p>

## Installation

### Pre-compiled binary (recommended)

Download the binary for your platform from the [Releases](https://github.com/malmarri/ldprio/releases/latest) page — no Rust, PLINK, or Python required:

**Linux (x86_64, static — works on any HPC cluster):**
```bash
curl -sL https://github.com/malmarri/ldprio/releases/latest/download/ldprio-linux-x86_64.tar.gz | tar -xz
chmod +x ldprio
```

**macOS (Apple Silicon):**
```bash
curl -sL https://github.com/malmarri/ldprio/releases/latest/download/ldprio-macos-arm64.tar.gz | tar -xz
chmod +x ldprio
```

*(Intel macOS `ldprio-macos-x86_64.tar.gz` and Linux ARM64 `ldprio-linux-aarch64.tar.gz` are also available on the Releases page.)*

### Build from source

Requires the Rust toolchain (`cargo`):

```bash
cargo build --release
# Binary is written to ./target/release/ldprio
```

Or install onto your `PATH`:

```bash
cargo install --path .
```

## Quick start

```bash
ldprio \
    --bfile panel_qc \
    --priority-samples priority_samples.txt \
    --ld-samples calculate_ld_samples.txt \
    --make-bed \
    --out panel_pruned
```

Writes `panel_pruned.snplist` (and `.bed/.bim/.fam` with `--make-bed`).

Add `--report` to also write `panel_pruned.report.html`, a self-contained page that shows, for each priority sample, how many of its SNPs were kept, along with the panel size through the run and any warnings.

Sample list files are one sample per line, either `IID` or `FID<TAB>IID`.

Only autosomes (chromosomes 1–22) are pruned and written out; X, Y,
mitochondrial and unplaced contigs are dropped.

## Choosing the LD samples

`--ld-samples` is required. It restricts r² estimation to high-quality samples
with reliable diploid calls.

Suggested to use at least 30 samples. With few samples the r² estimates are noisy:
unrelated SNPs pass the threshold by chance. `ldprio` prints a warning below 30
and carries on.

Every sample you list is used, even if the `.fam` file records parents for
it (columns 3–4).

## How it works

1. **Priority pool** — streams the PLINK 1 `.bed` file directly and records every SNP at least one of the `--priority-samples` has a call at.
2. **LD graph** — computes exact unphased genotype Pearson r² on `--ld-samples` across the `--window` variant window using bitwise `POPCNT` parallelized across CPU cores, keeping pairs with r² ≥ `--r2`.
3. **Greedy selection** — maximal independent set over that graph: keep a
   SNP, drop its LD partners. Priority-pool SNPs are visited first, in the
   order set by `--weighting` (see [Choosing `--weighting`](#choosing---weighting)),
   then the rest in genomic order.
4. **Cleanup** — re-evaluates LD on the surviving set in memory and resolves any
   residual conflicts with the same priority-first greedy rule as step 3,
   iterated to convergence.

Step 4 exists because a single pass uses a static window over the
full marker set, so it misses a small number of long-range pairs that only
appear once nearby markers have already been dropped. Convergence typically
takes 1–4 passes and removes <0.5% of markers.

The output is verified LD-independent *within the `--window` variant-count
window used throughout* (the same scope as PLINK's `--indep-pairwise`/`--r2`
guarantee, not a literal unbounded genome-wide claim): the final cleanup
pass finds zero residual pairs within that window. If it does not converge
within `--max-cleanup` passes, `ldprio` exits non-zero.

### Scaling to more priority samples

Nominating more samples grows the priority pool (a SNP joins as soon as any
one nominated sample covers it), so more blocks are contested between
nominated samples and each one keeps a smaller share of its sites.
Nominate only the samples that actually need the help.

### Choosing `--weighting`

When several nominated samples cover different SNPs in the same LD block,
only one SNP survives, so someone loses. `--weighting` decides who:

- **`adaptive`** (default) — combines coverage-normalized dynamic proportional
  fairness ($w_i = \frac{1}{c_i (k_i + 1)}$, where $c_i$ is sample $i$'s total
  called SNPs and $k_i$ is how many SNPs it has kept so far) with LD-graph
  awareness (GWMIN2 ratio $\frac{w(u)}{w(u) + \sum_{v \in N(u)} w(v)}$). This
  protects ultra-sparse samples when coverage varies, splits contested blocks
  fairly when coverage is similar, and avoids "hub" SNPs that would block
  multiple other priority sites on decaying-LD chromosomes.
- **`inverse`** — static weighting: each SNP scores the sum of $1/c_i$ over
  the nominated samples covering it, and higher scores win. Simple and
  favours the sparsest samples, though similar-coverage samples can experience
  winner-take-all behaviour.
- **`fair`** — samples take turns: the one with the fewest SNPs kept so far
  (sparsest first on ties) keeps its next available site, preferring sites
  shared with more nominated samples. Evens out the absolute number of SNPs
  kept per sample.

For most cohorts, the default `adaptive` works best across both uniform and
varying coverage. With one nominated sample, all three modes behave similarly
(with `adaptive` additionally avoiding high-degree LD hub SNPs on decaying LD).

### Preferring transversions

Ancient-DNA damage turns C into T (and G into A on the other strand), so
transition calls (C↔T, G↔A) in ancient samples are less reliable than
transversions. `--prefer-transversions` keeps a transversion over a
transition whenever the two are otherwise equal: between SNPs with the same
priority, and among SNPs no nominated sample covers. It never overrides the
coverage weighting, so the priority samples keep essentially the same
number of calls.

If you will remove transitions later anyway, filter the panel to
transversions *before* running `ldprio` instead. In testing, that kept the
priority samples just as many transversion calls and produced a panel with
far more transversion SNPs than any preference setting could.

### Notes

Filter for rare alleles and monomorphic SNPs beforehand.
SNPs that are monomorphic (or zero-variance) in `--ld-samples` have undefined r² —
`ldprio` cannot evaluate their LD and they always survive, whether or not
they are actually redundant with each other among the samples that matter for
them. This is checked and reported
(`note: N/M priority-pool SNPs are monomorphic (or zero-variance) in --ld-samples...`)
since it disproportionately affects priority-pool SNPs: sites present in ancients
but rare or absent in the modern LD-estimation sample are exactly the case.

## Options

| option | default | description |
| --- | --- | --- |
| `--bfile` | required | input PLINK 1 prefix (`.bed/.bim/.fam`) |
| `--out` | required | output prefix |
| `--priority-samples` | required | samples whose covered SNPs to favour |
| `--ld-samples` | required | samples LD is estimated from (at least 30 recommended) |
| `--window` | 200 | window, in variants |
| `--step` | 25 | unused (kept for CLI compatibility) |
| `--r2` | 0.4 | r² threshold |
| `--maf` | off | minimum minor allele frequency in `--ld-samples` (e.g. `0.01`) |
| `--geno` | off | maximum missing genotype rate in `--ld-samples` (e.g. `0.05`) |
| `--hwe` | off | minimum Hardy-Weinberg equilibrium exact p-value in `--ld-samples` (e.g. `1e-6`) |
| `--weighting` | `adaptive` | how blocks contested between nominated samples are settled: `adaptive`, `inverse`, or `fair` |
| `--prefer-transversions` | off | among otherwise equal SNPs, keep transversions over transitions |
| `--report` | off | also write a self-contained HTML report (`<out>.report.html`) |
| `--make-bed` | off | also write the pruned PLINK fileset (`<out>.{bed,bim,fam}`) |
| `--max-cleanup` | 10 | cap on cleanup passes |
| `-t`, `--threads` | 0 | number of threads (`0` = all available CPU cores) |

## Validation & benchmarks

### Real data: Southern Arc Natufians (1240k)

Evaluated on the full public 1240k dataset from Lazaridis, Alpaslan-Roodenberg et al., *The genetic history of the Southern Arc: a bridge between West Asia and Europe*, Science (2022) ([Harvard Dataverse: `doi:10.7910/DVN/3AR0CD`](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/3AR0CD)) — **1,150,639 autosomal SNPs** across **5,940 individuals**.

LD was estimated from 200 modern high-coverage diploid (`.DG`) samples (`--window 200 --r2 0.4`), nominating all 6 Epipaleolithic Natufian (`ISR_Natufian_EpiP`) individuals as `--priority-samples` (priority pool: 200,048 SNPs, 17.4% of the panel). `ldprio` (`--weighting adaptive`) pruned the entire 1.7 GB dataset to verified zero residual LD in **~5 seconds**, keeping **450,083 SNPs** (vs. **376,252 SNPs** for standard `plink --indep-pairwise 200 25 0.4`):

| Natufian sample | Input calls | Standard PLINK | `ldprio` (`adaptive`) | Raw gain | Share-of-panel gain* |
| --- | ---: | ---: | ---: | ---: | ---: |
| `I1687_d` | 12,883 | 4,352 (33.8%) | **12,189 (94.6%)** | **+180%** (+7,837) | **2.34×** |
| `I1685_d` | 15,262 | 5,270 (34.5%) | **13,630 (89.3%)** | **+159%** (+8,360) | **2.16×** |
| `I0861_d` | 16,613 | 5,692 (34.3%) | **14,252 (85.8%)** | **+150%** (+8,560) | **2.09×** |
| `I1069_published` | 26,699 | 8,946 (33.5%) | **21,029 (78.8%)** | **+135%** (+12,083) | **1.96×** |
| `I1690_published` | 29,423 | 9,829 (33.4%) | **21,839 (74.2%)** | **+122%** (+12,010) | **1.86×** |
| `I1072_d` | 122,713 | 41,854 (34.1%) | **78,850 (64.3%)** | **+88%** (+36,996) | **1.58×** |
| **Mean share** | | **33.9%** | **81.2%** | | |

*\*Share-of-panel gain adjusts for the 1.20× difference in overall pruned panel size (`450,083` vs. `376,252` SNPs).*

Adding `--prefer-transversions` kept the same **81.2%** mean Natufian retention while raising transversions to **25.3%** of the pruned panel (113,452 of 447,947 SNPs).

### Automated tests

Run the built-in Rust test suite (bitwise vs. naive r² accuracy, multi-pass cleanup, edge cases, and a 100,000-SNP multithreading determinism test — no PLINK required):

```bash
cargo test --release
```

To also compare against `plink --indep-pairwise` on the bundled synthetic panels (`test_panel` and `chain_panel`, requires PLINK 1.9 on `PATH`):

```bash
cargo build --release
cd test && ./run_test.sh
```
