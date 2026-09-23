# ldprio

LD pruning that prioritizes the SNPs your low-coverage samples actually cover.

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
sites covered by the low coverage samples you nominate, and settle
blocks contested between them by `--weighting`.

Well-covered samples are unaffected: they had an equivalent marker either
way. Sparse samples keep substantially more of their real data.

<p align="center">
  <img src="docs/ld-block.svg" alt="One LD block of five redundant SNPs: modern samples have calls at all five, the ancient sample only at s3. plink keeps s1, leaving the ancient sample with no SNP from this block; ldprio keeps s3, so the ancient sample keeps its site. Modern samples keep one SNP either way." width="860">
</p>

## Quick start

```bash
ldprio.py \
    --bfile panel_qc \
    --priority-samples low_coverage_samples.txt \
    --ld-samples modern_samples.txt \
    --autosomes-only --make-bed \
    --out panel_pruned
```

Writes `panel_pruned.snplist` (and `.bed/.bim/.fam` with `--make-bed`).

Sample list files are one sample per line, either `IID` or `FID<TAB>IID`.

## Estimate LD from modern samples

Use `--ld-samples` to restrict r² estimation to high-quality modern diploids.
If you omit `--ld-samples` the whole cohort is used and a warning is printed.

Every sample you list is used, even if the `.fam` file records parents for
it (columns 3–4). By default PLINK quietly leaves such samples out of LD
estimates; ldprio turns that off so none of your LD samples are dropped.
Most aDNA and reference panels record no parents (`0 0`), in which case
this changes nothing.

## How it works

1. **Priority pool** — every SNP at least one of the `--priority-samples`
   has a call at.
2. **LD graph** — `plink --r2` on `--ld-samples`, keeping pairs above `--r2`.
3. **Greedy selection** — maximal independent set over that graph: keep a
   SNP, drop its LD partners. Priority-pool SNPs are visited first, in the
   order set by `--weighting` (see [Choosing `--weighting`](#choosing---weighting)),
   then the rest in genomic order.
4. **Cleanup** — re-run `plink --r2` on the surviving set and resolve any
   residual conflicts with the same priority-first greedy rule as step 3,
   iterated to convergence.

Step 4 exists because a single `--r2` pass uses a static window over the
full marker set, so it misses a small number of long-range pairs that only
appear once nearby markers have already been dropped. Convergence typically
takes 3–4 passes and removes <0.5% of markers.

The output is verified LD-independent *within the `--window` variant-count
window used throughout* (the same scope PLINK's own `--indep-pairwise`/`--r2`
guarantee, not a literal unbounded genome-wide claim): the final cleanup
pass finds zero residual pairs within that window. If it does not converge
within `--max-cleanup` passes, `ldprio` exits non-zero rather than silently
handing back a dirty panel — it still writes the (unconverged) output, so
you can inspect it, but a calling pipeline under `set -e` will stop.

### Scaling to more priority samples

Nominating more samples grows the priority pool (a SNP joins as soon as any
one nominated sample covers it), so more blocks are contested between
nominated samples and each one keeps a smaller share of its sites. The
`[2/5] priority pool: X%` line reports where you are; above ~70% you'll see
a warning. Nominate only the samples that actually need the help.

Measured by sweeping K (number of nominated low-coverage samples) on a
2,000-SNP / 200-block synthetic panel with 200 modern + 100 low-coverage
samples. "Mean" is the average share of each nominated sample's own calls
that survive pruning; "fewest" is the fewest SNPs any one nominated sample
keeps. Every run keeps 200 SNPs (one per block).

Every sample at ~4% coverage (~80 calls each):

| K | pool % | plink mean | plink fewest | inverse mean | inverse fewest | fair mean | fair fewest |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 4.1% | 8.5% | 7 | 82.9% | 68 | 82.9% | 68 |
| 5 | 18.4% | 8.5% | 4 | 51.8% | 25 | 50.3% | **39** |
| 10 | 33.5% | 8.9% | 4 | 39.3% | 23 | 38.5% | **31** |
| 20 | 57.1% | 9.7% | 4 | 30.1% | 14 | 29.6% | **23** |
| 50 | 87.3% | 9.7% | 2 | 22.6% | 10 | 21.9% | **15** |
| 100 | 98.0% | 10.3% | 2 | 18.5% | 7 | 18.0% | **11** |

Coverage varying 0.5–20% between samples (5 to 398 calls each), as in a
typical ancient-DNA cohort:

| K | pool % | plink mean | plink fewest | inverse mean | inverse fewest | fair mean | fair fewest |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3.2% | 10.8% | 7 | 92.3% | 60 | 92.3% | 60 |
| 5 | 16.9% | 8.7% | 0 | **70.8%** | **13** | 69.0% | 12 |
| 10 | 41.0% | 8.5% | 0 | **50.9%** | **13** | 48.2% | 10 |
| 20 | 68.3% | 8.8% | 0 | **41.3%** | **11** | 38.9% | 8 |
| 50 | 95.2% | 9.0% | 0 | **31.4%** | 5 | 29.6% | **7** |
| 100 | 99.9% | 9.2% | 0 | **24.8%** | 2 | 23.0% | **4** |

Generated by `test/run_ksweep.sh`; seeded, so these numbers are
reproducible.

### Choosing `--weighting`

When several nominated samples cover different SNPs in the same LD block,
only one SNP survives, so someone loses. `--weighting` decides who:

- **`inverse`** (default) — each SNP scores the sum of 1/(total called SNPs)
  over the nominated samples covering it, and higher scores win. A sample
  with 10 calls contributes 0.1 to each of its sites, one with 1,000 calls
  contributes 0.001, so the sparsest samples win. This maximises the
  average share of each sample's own data that survives. Weak spot: when
  coverage is nearly equal, small differences (say 75 vs 94 calls) still
  rank the samples strictly, so the slightly better-covered one loses most
  contests.
- **`fair`** — samples take turns: the one with the fewest SNPs kept so far
  (sparsest first on ties) keeps its next available site, preferring sites
  shared with more nominated samples. This evens out the absolute number of
  SNPs kept per sample, at a small cost in average retention.

Use `inverse` when coverage varies a lot between the samples you nominate
(the usual ancient-DNA case) and `fair` when they have similar coverage.
With one nominated sample the two are identical.

An earlier version weighted each SNP by the *number* of nominated samples
covering it. That let well-covered samples outvote sparse ones: with
varying coverage the sparsest sample kept 0–3 SNPs at K ≥ 10, barely better
than plain PLINK, so it is no longer offered.

### What it does not do

Finding the *maximum* independent set is NP-hard. This is a greedy
heuristic, as is PLINK's own, so the two land on different valid solutions
that can differ noticeably in size. On clean, disjoint LD blocks both keep
one SNP per block. Where LD decays gradually along the chromosome, as in
real data, `ldprio` kept ~29% more SNPs than `plink --indep-pairwise` on the
test panel below, and PLINK's own check still found no residual LD in
either. So the output is as LD-independent as PLINK's by PLINK's own
criterion, but it is not the same panel: compare samples by their share of
the panel, not raw SNP counts, or part of the apparent gain is just panel
size (`test/run_test.sh` reports both).

SNPs monomorphic in `--ld-samples` have undefined r² and PLINK's `--r2`
simply omits them from its output — `ldprio` cannot evaluate their LD and
they always survive, whether or not they're actually redundant with each
other among the samples that matter for them. This is checked and reported
(`note: N/M priority-pool SNPs are monomorphic in --ld-samples...`) since it
disproportionately affects priority-pool SNPs: sites present in ancients but
rare or absent in the modern LD-estimation sample are exactly the case.

## Options

| option | default | description |
| --- | --- | --- |
| `--bfile` | required | input PLINK prefix |
| `--out` | required | output prefix |
| `--priority-samples` | required | samples whose covered SNPs to favour |
| `--ld-samples` | whole cohort | samples LD is estimated from |
| `--window` | 200 | window, in variants |
| `--step` | 25 | unused (kept for CLI compatibility) — see note below |
| `--r2` | 0.4 | r² threshold |
| `--weighting` | `inverse` | how blocks contested between nominated samples are settled: `inverse` or `fair` (see [Choosing `--weighting`](#choosing---weighting)) |
| `--autosomes-only` | off | restrict to autosomes before pruning |
| `--make-bed` | off | also write the pruned PLINK fileset |
| `--max-cleanup` | 10 | cap on cleanup passes |
| `--keep-intermediates` | off | keep working files (the `.ld` can be ~0.5 GB) |
| `--plink` | `plink` | PLINK 1.9 executable (plink2 is not supported) |

Defaults match the Lazaridis et al. 2016 convention (`200 25 0.4`). `--step`
was only ever consumed by `plink --indep-pairwise`; cleanup no longer uses
it (see [How it works](#how-it-works)), so it has no effect regardless of
value. Kept as a flag rather than removed so existing invocations that pass
it don't break.

## Test

```bash
cd test && ./run_test.sh
```

Runs against two synthetic panels, each with 200 fully called modern
diploids and 5 pseudo-haploid samples at ~4% coverage:

- `test_panel`: 2,000 SNPs in 200 disjoint LD blocks of 10
- `chain_panel`: 3,000 SNPs with LD decaying along each chromosome

Each panel is run with both `--weighting` modes. For each run it checks
that the output is LD-independent, that every
priority sample's *share of the panel* grows (so a gain can't come from
ldprio simply keeping more SNPs), and that modern samples are not penalised.

Expected output with the default `--weighting inverse` (per-panel
summaries; the `fair` runs print tables in the same format):

```
panel size: plink 200 SNPs, ldprio 200 SNPs
sample           plink   ldprio   change   share-of-panel
ancient1             6       43    +617%            7.2x
ancient2            13       30    +131%            2.3x
ancient3             6       39    +550%            6.5x
ancient4             8       63    +688%            7.9x
ancient5            14       36    +157%            2.6x
modern(mean)       200      200    +0.0%

panel size: plink 475 SNPs, ldprio 611 SNPs
sample           plink   ldprio   change   share-of-panel
ancient1            19       66    +247%            2.7x
ancient2            16       80    +400%            3.9x
ancient3            24       65    +171%            2.1x
ancient4            26       97    +273%            2.9x
ancient5            23       74    +222%            2.5x
modern(mean)       475      611   +28.6%
```

On the block panel the difference is entirely in *which* SNP was kept from
each block. On the chain panel ldprio also keeps more SNPs overall; the
share-of-panel column is the gain after allowing for that.

`test/make_test_data.py` and `test/make_chain_data.py` regenerate the
panels; both are seeded, so the numbers above are reproducible.

## Requirements

Python 3.7+ and PLINK 1.9 on `PATH` (or pass `--plink`). No Python
dependencies beyond the standard library.

PLINK 2 is not supported: ldprio relies on PLINK 1.9's pairwise LD report
(`--r2`) and allele-frequency output, which PLINK 2 replaces or drops.

`ldprio.py --version` prints the version; releases are tagged on GitHub
(`v1.1`, …).
