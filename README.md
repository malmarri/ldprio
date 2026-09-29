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
  <img src="docs/ld-block.svg" alt="One LD block of five redundant SNPs: high coverage samples have calls at all five, the ancient sample only at s3. plink keeps s1, leaving the ancient sample with no SNP from this block; ldprio keeps s3, so the ancient sample keeps its site. High coverage samples keep one SNP either way." width="860">
</p>

## Quick start

```bash
python3 ldprio.py \
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
with reliable diploid calls: modern samples, or high-coverage ancient ones.

Suggested to use at least 30 samples. With few samples the r² estimates are noisy:
unrelated SNPs pass the threshold by chance. ldprio prints a warning below 30
and carries on.

Every sample you list is used, even if the `.fam` file records parents for
it (columns 3–4). By default PLINK quietly leaves such samples out of LD
estimates; ldprio turns that off so none of your LD samples are dropped.


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
within `--max-cleanup` passes, `ldprio` exits non-zero error.

### Scaling to more priority samples

Nominating more samples grows the priority pool (a SNP joins as soon as any
one nominated sample covers it), so more blocks are contested between
nominated samples and each one keeps a smaller share of its sites.
Nominate only the samples that actually need the help.

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

Use `inverse` when coverage varies between the samples you nominate
(the usual ancient-DNA case) and `fair` when they have similar coverage.
With one nominated sample the two are identical.

### Preferring transversions

Ancient-DNA damage turns C into T (and G into A on the other strand), so
transition calls (C↔T, G↔A) in ancient samples are less reliable than
transversions. `--prefer-transversions` keeps a transversion over a
transition whenever the two are otherwise equal: between SNPs with the same
priority, and among SNPs no nominated sample covers. It never overrides the
coverage weighting, so the priority samples keep essentially the same
number of calls.

If you will remove transitions later anyway, filter the panel to
transversions *before* running ldprio instead. In testing, that kept the
priority samples just as many transversion calls and produced a panel with
far more transversion SNPs than any preference setting could.

### Notes

Filter for rare alleles and monomorphic SNPs before hand. 
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
| `--ld-samples` | required | samples LD is estimated from (at least 30 recommended) |
| `--window` | 200 | window, in variants |
| `--step` | 25 | unused (kept for CLI compatibility) — see note below |
| `--r2` | 0.4 | r² threshold |
| `--weighting` | `inverse` | how blocks contested between nominated samples are settled: `inverse` or `fair` (see [Choosing `--weighting`](#choosing---weighting)) |
| `--prefer-transversions` | off | among otherwise equal SNPs, keep transversions over transitions (see [Preferring transversions](#preferring-transversions)) |
| `--report` | off | also write an HTML report of the results (`<out>.report.html`) |
| `--make-bed` | off | also write the pruned PLINK fileset |
| `--max-cleanup` | 10 | cap on cleanup passes |
| `--keep-intermediates` | off | keep working files (the `.ld` can be ~0.5 GB) |
| `--plink` | `plink` | PLINK 1.9 executable (plink2 is not supported) |

`--step` was only ever consumed by `plink --indep-pairwise`; cleanup no longer uses
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
(`--r2`) and allele-frequency output.

`ldprio.py --version` prints the version; releases are tagged on GitHub
(`v1.1`, …).
