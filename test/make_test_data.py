#!/usr/bin/env python3
"""
Generate the synthetic test panel for ldprio.

Writes test_panel.ped/.map, which `make_test_data.sh` converts to PLINK
binary format. Fully synthetic -- no real genotypes.

Panel design, chosen so the effect of prioritisation is unambiguous:

  * 2,000 SNPs in 200 LD blocks of 10. Within a block the SNPs are copies of
    a founder genotype with a small per-call flip rate, giving r2 well above
    the default 0.4 threshold; between blocks they are independent. So a
    correct pruner keeps roughly one SNP per block and there is a genuine
    choice of *which* one -- exactly the situation prioritisation targets.

  * 200 "modern" samples: diploid, fully called. LD is estimated from these.

  * 5 "ancient" samples: ~4% called, and every call homozygous, mimicking
    pseudo-haploid ancient DNA. Their covered sites are sparse and
    essentially random, so under unprioritised pruning they mostly lose the
    coin flip within each block.

With 10 SNPs per block and ~4% coverage, an ancient sample covers a site in
only a minority of blocks, and without prioritisation it keeps roughly 4% of
them -- versus close to the full covered set when prioritised.
"""

import random

random.seed(20240921)

N_BLOCKS = 200
BLOCK_SIZE = 10
N_MODERN = 200
N_ANCIENT = 5
FLIP_RATE = 0.02         # within-block decorrelation; keeps r2 >> 0.4
ANCIENT_COVERAGE = 0.04  # fraction of sites with a call

N_SNPS = N_BLOCKS * BLOCK_SIZE
ALLELES = ("A", "G")


def draw_genotype(freq):
    """Diploid genotype as an allele-2 dosage, under Hardy-Weinberg."""
    return sum(1 for _ in range(2) if random.random() < freq)


# genotype matrix: rows = samples, cols = SNPs (dosages 0/1/2, None = missing)
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
                # copy the founder, flipping occasionally to soften r2 to <1
                modern_rows[i][col] = (draw_genotype(freq)
                                       if random.random() < FLIP_RATE else g)

# ancient samples: sparse, homozygous-only calls drawn from the same
# block structure so their sites are real sites, just thinly covered
ancient_rows = []
for a in range(N_ANCIENT):
    row = [None] * N_SNPS
    for col in range(N_SNPS):
        if random.random() < ANCIENT_COVERAGE:
            row[col] = 0 if random.random() < 0.5 else 2   # pseudo-haploid
    ancient_rows.append(row)

with open("test_panel.map", "w") as fh:
    for b in range(N_BLOCKS):
        for k in range(BLOCK_SIZE):
            col = b * BLOCK_SIZE + k
            chrom = (b // 20) + 1                 # spread across chr1-10
            bp = (b % 20) * 500_000 + k * 2_000   # blocks 500kb apart
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


with open("test_panel.ped", "w") as fh:
    for i, row in enumerate(modern_rows):
        fh.write(ped_line("MODERN", f"modern{i + 1:03d}", row))
    for i, row in enumerate(ancient_rows):
        fh.write(ped_line("ANCIENT", f"ancient{i + 1}", row))

with open("ld_samples.txt", "w") as fh:
    fh.writelines(f"modern{i + 1:03d}\n" for i in range(N_MODERN))

with open("priority_samples.txt", "w") as fh:
    fh.writelines(f"ancient{i + 1}\n" for i in range(N_ANCIENT))

print(f"wrote test_panel.ped/.map: {N_SNPS} SNPs "
      f"({N_BLOCKS} blocks x {BLOCK_SIZE}), "
      f"{N_MODERN} modern + {N_ANCIENT} ancient samples")
print("wrote ld_samples.txt, priority_samples.txt")
