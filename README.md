# ldprio

LD pruning that prioritizes the SNPs your low-coverage samples actually cover.

## The problem

LD pruning is commonly performed for PCA and model-based clustering analysis like ADMIXTURE,
which requires relatively indpendent SNPs. The most used method for this is `plink --indep-pairwise`.
PLINK LD prunning removes correlated SNPs at random, when a block contains several mutually
redundant SNPs it keeps an arbitrary one — and for a low-coverage
sample, covering a sparse number sites, that is
usually not one of the few sites it has data for, resulting in an unnecessary loss for ultra-low coverage samples

## The fix

Do the pruning as a single genome-wide pass, so every LD relationship is
still evaluated exactly once — but choose SNP within each block toward
sites covered by the low coverage samples you nominate.

Well-covered samples are unaffected: they had an equivalent marker either
way. Sparse samples keep substantially more of their real data.

```
                 LD block (all mutually redundant)
                 ┌────┬────┬────┬────┬────┐
     SNPs        │ s1 │ s2 │ s3 │ s4 │ s5 │
     low-cov has │    │    │ ✓  │    │    │
                 └────┴────┴────┴────┴────┘
  plink keeps      ^                          low-cov sample gets nothing
  ldprio keeps               ^                low-cov sample keeps its site
```

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

## How it works

1. **Priority weight** — per-SNP count of non-missing calls among
   `--priority-samples` (0 = not in the pool).
2. **LD graph** — `plink --r2` on `--ld-samples`, keeping pairs above `--r2`.
3. **Greedy selection** — maximal independent set over that graph, visiting
   priority-pool SNPs first (highest weight first, so a SNP covered by more
   of your nominated samples wins a block over one covered by fewer).
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
one nominated sample covers it), and past a certain point the pool
approaches the whole panel — at which point priority-first traversal is
just genomic order again, and the benefit to any individual sample
decreases. Coverage-weighting (a SNP covered by more nominated samples wins a
block over one covered by fewer) keeps the tool useful much further into
this range, but the pool still eventually saturates. The `[2/5] priority pool: X%` line reports
where you are; above ~70% you'll see a warning, since by then most blocks'
winners are governed by genomic order rather than priority. Nominate only
the samples that actually need the help.

Measured by sweeping K (number of nominated low-coverage samples, each at
~4% coverage) on a 2,000-SNP / 200-block synthetic panel, 200 modern + 100
low-coverage samples, baseline (plain `plink --indep-pairwise`) `lowcov001`
= 7 called SNPs out of 200 kept throughout:

| K | pool % of panel | priority retained | `lowcov001` gain vs. plink |
| ---: | ---: | ---: | ---: |
| 1 | 4.1% | 82.9% | 9.7× |
| 5 | 18.4% | 46.3% | 6.1× |
| 10 | 33.5% | 29.0% | 5.1× |
| 20 | 57.1% | 17.5% | 3.1× |
| 50 | 87.3% (warning fires) | 11.5% | 3.1× |
| 100 | 98.0% (warning fires) | 10.2% | 2.6× |

Total SNPs kept stays at 200 throughout (one per block) — you never pay for
this in panel size, prioritisation only changes *which* SNP wins each
block. But the retention rate drops fast as K grows, and even the single
tracked sample's own gain over plain `plink` more than halves between K=1
and K=100 despite that sample's own coverage never changing — it's diluted
by the other nominated samples, not by anything about itself. Generated by
`test/make_ksweep_data.py`; seeded, so these numbers are reproducible.

### What it does not do

Finding the *maximum* independent set is NP-hard. This is a greedy
heuristic, as is PLINK's own, so the two land on different valid solutions
of slightly different size. Do not read a difference in total SNP count
between `ldprio` and `plink --indep-pairwise` as a gain or a loss — it is an
artefact of which heuristic you ran. The meaningful number is the retention
rate *within the priority pool*, which is what the tool actually changes.

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
| `--autosomes-only` | off | restrict to autosomes before pruning |
| `--make-bed` | off | also write the pruned PLINK fileset |
| `--max-cleanup` | 10 | cap on cleanup passes |
| `--keep-intermediates` | off | keep working files (the `.ld` can be ~0.5 GB) |
| `--plink` | `plink` | plink executable |

Defaults match the Lazaridis et al. 2016 convention (`200 25 0.4`). `--step`
was only ever consumed by `plink --indep-pairwise`; cleanup no longer uses
it (see [How it works](#how-it-works)), so it has no effect regardless of
value. Kept as a flag rather than removed so existing invocations that pass
it don't break.

## Test

```bash
cd test && ./run_test.sh
```

Runs against a synthetic panel (2,000 SNPs in 200 LD blocks of 10; 200 fully
called modern diploids; 5 pseudo-haploid samples at ~4% coverage) and checks
that the output is LD-independent, that every priority sample gains, and
that modern samples are not penalised.

Expected output:

```
sample              plink       ldprio     change
ancient1                6           35      +483%
ancient2               13           54      +315%
ancient3                6           44      +633%
ancient4                8           35      +338%
ancient5               14           43      +207%
modern(mean)          200          200      +0.0%
```

Both panels keep 200 SNPs — one per LD block — and both are LD-clean. The
difference is entirely in *which* SNP was kept from each block.

`test/make_test_data.py` regenerates the panel; it is seeded, so the numbers
above are reproducible.

## Requirements

Python 3.6+ and PLINK 1.9 on `PATH` (or pass `--plink`). No Python
dependencies beyond the standard library.
