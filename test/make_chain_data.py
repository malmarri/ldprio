#!/usr/bin/env python3
"""
Generate the decaying-LD test panel for ldprio (chain_panel.ped/.map).

Complements make_test_data.py, whose panel is 200 disjoint cliques of
near-identical SNPs. Real LD decays with distance instead: each SNP here
copies the previous SNP's genotype except for a random 12% of samples, so
r2 falls off along the chromosome and "blocks" overlap. On this structure
ldprio and plink --indep-pairwise keep different numbers of SNPs, which is
why run_test.sh compares samples by share of the panel, not raw counts.

Same sample IDs as make_test_data.py, so ld_samples.txt and
priority_samples.txt serve both panels. Convert with:
  plink --file chain_panel --make-bed --out chain_panel
"""

import random

random.seed(20260923)

N_CHROM = 10
SNPS_PER_CHROM = 300
N_MODERN = 200
N_ANCIENT = 5
REDRAW_RATE = 0.12
ANCIENT_COVERAGE = 0.04
ALLELES = ("A", "G")

N_SNPS = N_CHROM * SNPS_PER_CHROM


def draw_genotype(freq):
    return sum(1 for _ in range(2) if random.random() < freq)


modern_rows = [[None] * N_SNPS for _ in range(N_MODERN)]
for col in range(N_SNPS):
    freq = random.uniform(0.15, 0.85)
    new_chrom = col % SNPS_PER_CHROM == 0
    for row in modern_rows:
        if new_chrom or random.random() < REDRAW_RATE:
            row[col] = draw_genotype(freq)
        else:
            row[col] = row[col - 1]

ancient_rows = []
for _ in range(N_ANCIENT):
    ancient_rows.append([(0 if random.random() < 0.5 else 2)
                         if random.random() < ANCIENT_COVERAGE else None
                         for _ in range(N_SNPS)])

with open("chain_panel.map", "w") as fh:
    for col in range(N_SNPS):
        chrom = col // SNPS_PER_CHROM + 1
        bp = (col % SNPS_PER_CHROM + 1) * 10_000
        fh.write(f"{chrom}\tcs{col:05d}\t0\t{bp}\n")


def ped_line(fid, iid, row):
    genos = {None: "0 0", 0: "A A", 1: "A G", 2: "G G"}
    return (f"{fid}\t{iid}\t0\t0\t0\t-9\t"
            + "\t".join(genos[d] for d in row) + "\n")


with open("chain_panel.ped", "w") as fh:
    for i, row in enumerate(modern_rows):
        fh.write(ped_line("MODERN", f"modern{i + 1:03d}", row))
    for i, row in enumerate(ancient_rows):
        fh.write(ped_line("ANCIENT", f"ancient{i + 1}", row))

print(f"wrote chain_panel.ped/.map: {N_SNPS} SNPs on {N_CHROM} chromosomes, "
      f"{N_MODERN} modern + {N_ANCIENT} ancient samples")
