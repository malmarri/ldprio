#!/usr/bin/env python3
"""
Generate the synthetic panel behind the README's "Scaling to more priority
samples" table.

Same block structure and per-sample coverage as make_test_data.py (2,000
SNPs in 200 LD blocks of 10; 200 modern samples LD is estimated from), but
100 low-coverage samples instead of 5, so --priority-samples can be swept
from K=1 to K=100 to show how the priority pool saturates. Writes
ksweep.ped/.map, ld_samples.txt, and one priority_K<N>.txt per K in
[1, 5, 10, 20, 50, 100] (each the first N of the 100 low-coverage samples,
so every larger K's list is a superset of every smaller one's).

run_ksweep.sh converts ksweep.ped/.map to PLINK binary and runs the sweep.
"""

import random

random.seed(20260921)

N_BLOCKS = 200
BLOCK_SIZE = 10
N_MODERN = 200
N_LOWCOV = 100
FLIP_RATE = 0.02
LOWCOV_COVERAGE = 0.04
K_VALUES = [1, 5, 10, 20, 50, 100]

N_SNPS = N_BLOCKS * BLOCK_SIZE
ALLELES = ("A", "G")


def draw_genotype(freq):
    """Diploid genotype as an allele-2 dosage, under Hardy-Weinberg."""
    return sum(1 for _ in range(2) if random.random() < freq)


modern_rows = [[None] * N_SNPS for _ in range(N_MODERN)]
for b in range(N_BLOCKS):
    freq = random.uniform(0.15, 0.85)
    founders = [draw_genotype(freq) for _ in range(N_MODERN)]
    for k in range(BLOCK_SIZE):
        col = b * BLOCK_SIZE + k
        for i, g in enumerate(founders):
            if k == 0:
                modern_rows[i][col] = g
            else:
                modern_rows[i][col] = (draw_genotype(freq)
                                       if random.random() < FLIP_RATE else g)

lowcov_rows = []
for a in range(N_LOWCOV):
    row = [None] * N_SNPS
    for col in range(N_SNPS):
        if random.random() < LOWCOV_COVERAGE:
            row[col] = 0 if random.random() < 0.5 else 2   # pseudo-haploid
    lowcov_rows.append(row)

with open("ksweep.map", "w") as fh:
    for b in range(N_BLOCKS):
        for k in range(BLOCK_SIZE):
            col = b * BLOCK_SIZE + k
            chrom = (b // 20) + 1
            bp = (b % 20) * 500_000 + k * 2_000
            fh.write(f"{chrom}\trs{col:05d}\t0\t{bp}\n")


def ped_line(fid, iid, row):
    calls = []
    for dosage in row:
        if dosage is None:
            calls.append("0 0")
        elif dosage == 0:
            calls.append(f"{ALLELES[0]} {ALLELES[0]}")
        elif dosage == 1:
            calls.append(f"{ALLELES[0]} {ALLELES[1]}")
        else:
            calls.append(f"{ALLELES[1]} {ALLELES[1]}")
    return f"{fid}\t{iid}\t0\t0\t0\t-9\t" + "\t".join(calls) + "\n"


with open("ksweep.ped", "w") as fh:
    for i, row in enumerate(modern_rows):
        fh.write(ped_line("MODERN", f"modern{i + 1:03d}", row))
    for i, row in enumerate(lowcov_rows):
        fh.write(ped_line("LOWCOV", f"lowcov{i + 1:03d}", row))

with open("ld_samples.txt", "w") as fh:
    fh.writelines(f"modern{i + 1:03d}\n" for i in range(N_MODERN))
for K in K_VALUES:
    with open(f"priority_K{K}.txt", "w") as fh:
        fh.writelines(f"lowcov{i + 1:03d}\n" for i in range(K))

print(f"wrote ksweep.ped/.map: {N_SNPS} SNPs, {N_MODERN} modern + "
      f"{N_LOWCOV} low-coverage samples")
print(f"wrote ld_samples.txt, priority_K<N>.txt for K in {K_VALUES}")
