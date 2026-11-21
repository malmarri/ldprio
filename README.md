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
- Only autosomes (chromosomes 1–22) are pruned and written out; for other species set the autosome count with `--autosome-num N`, as in PLINK.

## How it works

1. **Priority pool:** identifies SNPs covered by `--priority-samples`. Nominate the samples your analysis is about (see [Choosing priority samples](#choosing-priority-samples)).
2. **LD-sample QC (optional):** filters autosomes by `--geno`, `--maf`, and `--excess-het` on `--ld-samples`.
3. **LD graph & greedy selection:** computes exact pairwise Pearson r² on `--ld-samples` (use ≥30 high-quality diploids; avoid pseudo-haploids) across the window (`--window` variants, or `--window-kb`), keeps priority-pool SNPs first according to `--weighting`, then fills the rest in genomic order.
4. **Cleanup:** with `--window` (counted in variants, as in PLINK), pruning widens the physical span of the window, so LD checks are repeated in memory until zero residual pairs above `--r2` remain, and SNPs left without a kept LD partner are recovered. With `--window-kb` the window is a fixed distance, so the first pass is already LD-independent and is simply verified.

### Key flags

- **`--weighting`** (how blocks contested between priority samples are settled):
  - `adaptive` *(default)* — combines coverage-normalized fairness with LD-graph awareness; works well for both varying and uniform coverage.
  - `inverse` — static $1 / (\text{called SNPs})$ weighting; favours the sparsest samples.
  - `fair` — samples take turns; evens out total SNPs kept per sample.
- **`--maf` / `--geno` / `--excess-het`** — optional QC calculated on `--ld-samples` before pruning. `--excess-het auto` removes excess-heterozygosity artifacts (`p < 1e-4` for N < 100 else `1e-6`, and `-F_IS > 0.20`) while preserving Wahlund homozygote excess from population structure.
- **`--window` / `--window-kb`** — LD window as a variant count (default `200`, PLINK-compatible) or as a physical distance in kb from `.bim` positions (e.g. `--window-kb 500`; needs a position-sorted `.bim`). A kb window means the same distance on any array and needs no cleanup passes; compare against `plink --indep-pairwise 500kb 1 0.4`.
- **`--transv`** — breaks ties in favour of transversions over ancient-DNA-damaged transitions (C↔T, G↔A).

## Options

| Option | Default | Description |
| --- | --- | --- |
| `--bfile` | required | Input PLINK 1 prefix (`.bed/.bim/.fam`) |
| `--out` | required | Output prefix |
| `--priority-samples` | required | Samples whose covered SNPs to favour |
| `--ld-samples` | required | High-quality diploid samples for LD & QC (≥30 recommended) |
| `--window` | `200` | Window size in variants |
| `--window-kb` | off | Window as physical distance in kb (replaces `--window`) |
| `--r2` | `0.4` | r² threshold |
| `--maf` | off | Minimum MAF in `--ld-samples` (e.g. `0.01`) |
| `--geno` | off | Maximum missing rate in `--ld-samples` (e.g. `0.05`) |
| `--excess-het` | off | One-sided excess-het p-value in `--ld-samples` (e.g. `1e-6` or `auto`) |
| `--weighting` | `adaptive` | `adaptive`, `inverse`, or `fair` |
| `--transv` | off | Prefer transversions on ties |
| `--make-bed` | off | Write pruned `.bed/.bim/.fam` |
| `--report` | off | Write self-contained HTML report (`<out>.report.html`) |
| `--max-cleanup` | `10` | Maximum cleanup passes (`--window` only) |
| `--autosome-num` | `22` | Number of autosomes (non-human data; chromosomes `1..N` are pruned) |
| `-t`, `--threads` | `0` | Threads (`0` = all cores) |

## Benchmark: Southern Arc Natufians (1240k)

Evaluated on the full public 1240k dataset from Lazaridis, Alpaslan-Roodenberg et al., *Science* (2022) ([Harvard Dataverse: `doi:10.7910/DVN/3AR0CD`](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/3AR0CD)) — **1,150,639 autosomal SNPs × 5,940 individuals** — with the 6 Epipaleolithic Natufian (`ISR_Natufian_EpiP`) samples as `--priority-samples`.

- **LD samples:** all 300 modern diploids — every `.DG` sample with an `S_`, `A_` or `B_` IID prefix. The six archaic/ancient `.DG` genomes (Altai, Vindija, Denisova, Ust'-Ishim, Loschbour, Stuttgart) are excluded; all other samples in the dataset are pseudo-haploid.
  ```bash
  awk '$2 ~ /^[SAB]_.*\.DG$/ {print $1"\t"$2}' SouthernArc_Public.fam > ld_samples.txt
  ```
- **Same QC for both tools:** `--maf 0.01 --geno 0.05` on the LD samples (1,024,064 SNPs pass, identically in `ldprio` and PLINK), then `--window 200 --r2 0.4` (PLINK: `--indep-pairwise 200 25 0.4`).

`ldprio` pruned the 1.7 GB dataset in **~8 seconds** (Apple M2 Pro) to 372,301 SNPs vs. 345,784 for PLINK. Re-running `plink --indep-pairwise 200 25 0.4` on the `ldprio` output prunes **0 SNPs**.

| Natufian sample | Input calls | Standard PLINK | `ldprio` (`adaptive`) | Raw gain | Share-of-panel gain* |
| --- | ---: | ---: | ---: | ---: | ---: |
| `I1687_d` | 12,883 | 4,046 (31.4%) | **11,028 (85.6%)** | **+173%** | **2.53×** |
| `I1685_d` | 15,262 | 4,855 (31.8%) | **12,249 (80.3%)** | **+152%** | **2.34×** |
| `I0861_d` | 16,613 | 5,256 (31.6%) | **12,663 (76.2%)** | **+141%** | **2.24×** |
| `I1069_published` | 26,699 | 8,155 (30.5%) | **18,727 (70.1%)** | **+130%** | **2.13×** |
| `I1690_published` | 29,423 | 9,096 (30.9%) | **19,246 (65.4%)** | **+112%** | **1.97×** |
| `I1072_d` | 122,713 | 38,596 (31.5%) | **68,781 (56.1%)** | **+78%** | **1.66×** |
| **Mean share** | | **31.3%** | **72.3%** | | |

*\*Adjusted for the 1.08× difference in overall pruned panel size (`372,301` vs. `345,784` SNPs).*

## Choosing priority samples

Nominated samples compete for the same LD blocks, and with the default `--weighting adaptive` the sample with fewer calls tends to win a contested block. Each sample you add therefore reduces the gain of the others, and adding samples **sparser** than the ones you care about costs the most.

To measure this, we reran the [Natufian benchmark](#benchmark-southern-arc-natufians-1240k) above with the same settings, keeping its six Natufians (12k–123k calls) nominated throughout. We then also nominated more low-coverage ancient samples from the same dataset (1,937 samples with 5k–300k calls), adding them either in random order or sparsest first. The table shows the Natufians' gain: the number of their SNPs kept by `ldprio` relative to PLINK, averaged over the six.

| Extra samples nominated | Natufian gain (random extras) | Natufian gain (sparsest extras) |
| ---: | ---: | ---: |
| 0 | 2.31× | 2.31× |
| 10 | 2.09× | 2.00× |
| 25 | 1.93× | 1.75× |
| 50 | 1.81× | 1.53× |
| 100 | 1.67× | 1.36× |
| all 1,937 | 1.40× | 1.40× |

- **Nominate the samples your analysis is about**, not every low-coverage sample in the panel.
- **Adding sparser samples costs the most.** If your focal samples are not the sparsest, every much sparser sample you add takes contested blocks from them.
- **Large nominations are still safe.** With all 1,943 samples nominated, nominated samples kept 1.31× more of their SNPs than with PLINK on average, and samples never nominated still gained (1.14–1.28× on average across all runs, never below PLINK).

*Random-order results varied by at most ±0.1× across 5 replicate orders.*

## Performance

Runtime and peak memory on the full Natufian benchmark dataset (1,150,639 autosomal SNPs × 5,940 individuals) on an Apple M2 Pro laptop (16 GB RAM), with the [benchmark](#benchmark-southern-arc-natufians-1240k) settings:

- **Priority samples:** the 6 Natufians (`ISR_Natufian_EpiP`: I0861_d, I1069_published, I1072_d, I1685_d, I1687_d, I1690_published).
- **LD samples:** the 300 modern `.DG` diploids (IIDs starting `S_`, `A_` or `B_`).
- **Settings:** `--maf 0.01 --geno 0.05`, default window (200 variants, r² 0.4), `adaptive` weighting.

| Threads | Time | Peak memory |
| ---: | ---: | ---: |
| 1 (`-t 1`) | 18.9 s | 686 MB |
| 8 (`-t 8`) | 8.3 s | 683 MB |

Both runs give identical output.

## Tests

```bash
cargo test --release
```

## License

MIT — see [LICENSE](LICENSE).
