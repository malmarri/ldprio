# ldprio

LD pruning that preserves the SNPs your low-coverage samples actually cover — standalone Rust binary, no PLINK or Python required.

Standard `plink --indep-pairwise` picks one surviving SNP per LD block without looking at per-sample missingness, discarding most of the few sites a low-coverage or ancient sample covers. `ldprio` estimates LD from high-quality diploids (`--ld-samples`) and biases which SNP survives each block toward sites covered by your `--priority-samples`, keeping the output strictly LD-independent.

<p align="center">
  <img src="docs/ld-block.svg" alt="One LD block of five redundant SNPs: high coverage samples have calls at all five, the ancient sample only at s3. plink keeps s1, leaving the ancient sample with no SNP from this block; ldprio keeps s3, so the ancient sample keeps its site. High coverage samples keep one SNP either way." width="860">
</p>

## Install

Download a pre-compiled binary from [Releases](https://github.com/malmarri/ldprio/releases/latest) (no dependencies required):

```bash
# Linux (x86_64, static — works on any HPC cluster)
curl -sL https://github.com/malmarri/ldprio/releases/latest/download/ldprio-linux-x86_64.tar.gz | tar -xz && chmod +x ldprio

# macOS (Apple Silicon)
curl -sL https://github.com/malmarri/ldprio/releases/latest/download/ldprio-macos-arm64.tar.gz | tar -xz && chmod +x ldprio
```

*(Intel macOS and Linux ARM64 archives are also on the Releases page, or build from source with `cargo build --release`.)*

## Quick start

```bash
./ldprio \
    --bfile dataset \
    --priority-samples priority_samples.txt \
    --ld-samples ld_samples.txt \
    --maf 0.01 \
    --make-bed --report \
    --out dataset_pruned
```

- Writes `dataset_pruned.snplist`, plus `.bed/.bim/.fam` (`--make-bed`) and a self-contained HTML summary (`--report`).
- Sample lists: one sample per line (`IID` or `FID<TAB>IID`).
- Only autosomes (chromosomes 1–22) are pruned and written out.

## How it works

1. **Priority pool:** identifies SNPs covered by `--priority-samples`. Nominate only the sparse samples that need help.
2. **LD-sample QC (optional):** filters autosomes by `--geno`, `--maf`, and `--excess-het` on `--ld-samples`.
3. **LD graph & greedy selection:** computes exact pairwise Pearson r² on `--ld-samples` (use ≥30 high-quality diploids; avoid pseudo-haploids) across `--window` variants, keeps priority-pool SNPs first according to `--weighting`, then fills the rest in genomic order.
4. **Cleanup:** iterates LD checks in memory until zero residual pairs above `--r2` remain within `--window`.

### Key flags

- **`--weighting`** (how blocks contested between priority samples are settled):
  - `adaptive` *(default)* — combines coverage-normalized fairness with LD-graph awareness; works well for both varying and uniform coverage.
  - `inverse` — static $1 / (\text{called SNPs})$ weighting; favours the sparsest samples.
  - `fair` — samples take turns; evens out total SNPs kept per sample.
- **`--maf` / `--geno` / `--excess-het`** — optional QC calculated on `--ld-samples` before pruning. `--excess-het auto` removes excess-heterozygosity artifacts (`p < 1e-4` for N < 100 else `1e-6`, and `-F_IS > 0.20`) while preserving Wahlund homozygote excess from population structure.
- **`--prefer-transversions`** — breaks ties in favour of transversions over ancient-DNA-damaged transitions (C↔T, G↔A).

## Options

| Option | Default | Description |
| --- | --- | --- |
| `--bfile` | required | Input PLINK 1 prefix (`.bed/.bim/.fam`) |
| `--out` | required | Output prefix |
| `--priority-samples` | required | Samples whose covered SNPs to favour |
| `--ld-samples` | required | High-quality diploid samples for LD & QC (≥30 recommended) |
| `--window` | `200` | Window size in variants |
| `--r2` | `0.4` | r² threshold |
| `--maf` | off | Minimum MAF in `--ld-samples` (e.g. `0.01`) |
| `--geno` | off | Maximum missing rate in `--ld-samples` (e.g. `0.05`) |
| `--excess-het` | off | One-sided excess-het p-value in `--ld-samples` (e.g. `1e-6` or `auto`) |
| `--weighting` | `adaptive` | `adaptive`, `inverse`, or `fair` |
| `--prefer-transversions` | off | Prefer transversions on ties |
| `--make-bed` | off | Write pruned `.bed/.bim/.fam` |
| `--report` | off | Write self-contained HTML report (`<out>.report.html`) |
| `--max-cleanup` | `10` | Maximum cleanup passes |
| `-t`, `--threads` | `0` | Threads (`0` = all cores) |

## Benchmark: Southern Arc Natufians (1240k)

Evaluated on the full public 1240k dataset from Lazaridis, Alpaslan-Roodenberg et al., *Science* (2022) ([Harvard Dataverse: `doi:10.7910/DVN/3AR0CD`](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/3AR0CD)) — **1,150,639 autosomal SNPs × 5,940 individuals**, using 200 modern `.DG` diploids for LD (`--window 200 --r2 0.4`) and the 6 Epipaleolithic Natufian (`ISR_Natufian_EpiP`) samples as `--priority-samples`. `ldprio` pruned the 1.7 GB dataset to verified zero residual LD in **~5 seconds** (450,083 SNPs vs. 376,252 for `plink --indep-pairwise 200 25 0.4`):

| Natufian sample | Input calls | Standard PLINK | `ldprio` (`adaptive`) | Raw gain | Share-of-panel gain* |
| --- | ---: | ---: | ---: | ---: | ---: |
| `I1687_d` | 12,883 | 4,352 (33.8%) | **12,189 (94.6%)** | **+180%** | **2.34×** |
| `I1685_d` | 15,262 | 5,270 (34.5%) | **13,630 (89.3%)** | **+159%** | **2.16×** |
| `I0861_d` | 16,613 | 5,692 (34.3%) | **14,252 (85.8%)** | **+150%** | **2.09×** |
| `I1069_published` | 26,699 | 8,946 (33.5%) | **21,029 (78.8%)** | **+135%** | **1.96×** |
| `I1690_published` | 29,423 | 9,829 (33.4%) | **21,839 (74.2%)** | **+122%** | **1.86×** |
| `I1072_d` | 122,713 | 41,854 (34.1%) | **78,850 (64.3%)** | **+88%** | **1.58×** |
| **Mean share** | | **33.9%** | **81.2%** | | |

*\*Adjusted for the 1.20× difference in overall pruned panel size (`450,083` vs. `376,252` SNPs).*

## Tests

```bash
cargo test --release
```
