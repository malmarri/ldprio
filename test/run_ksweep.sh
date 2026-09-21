#!/bin/bash
# run_ksweep.sh -- reproduces the README's "Scaling to more priority
# samples" table: baseline (plain plink --indep-pairwise on moderns) vs.
# ldprio at K = 1, 5, 10, 20, 50, 100 nominated low-coverage samples.
#
# Requires plink 1.9 on PATH. Run from this directory.
set -euo pipefail
cd "$(dirname "$0")"

PLINK=${PLINK:-plink}
command -v "$PLINK" >/dev/null || { echo "FAIL: plink not on PATH"; exit 1; }

rm -rf out_ksweep && mkdir -p out_ksweep
(cd out_ksweep && python3 ../make_ksweep_data.py)

$PLINK --file out_ksweep/ksweep --make-bed --out out_ksweep/ksweep >/dev/null 2>&1

echo "=== baseline: plain plink --indep-pairwise (LD from moderns) ==="
awk '{print "MODERN\t"$1}' out_ksweep/ld_samples.txt > out_ksweep/ld.keep
$PLINK --bfile out_ksweep/ksweep --keep out_ksweep/ld.keep \
       --indep-pairwise 200 25 0.4 --allow-no-sex \
       --out out_ksweep/baseline >/dev/null 2>&1
$PLINK --bfile out_ksweep/ksweep --extract out_ksweep/baseline.prune.in \
       --make-bed --allow-no-sex --out out_ksweep/baseline_panel >/dev/null 2>&1
$PLINK --bfile out_ksweep/baseline_panel --missing --allow-no-sex \
       --out out_ksweep/baseline_panel_miss >/dev/null 2>&1
baseline_called=$(awk '$2=="lowcov001"{print $5-$4}' out_ksweep/baseline_panel_miss.imiss)
echo "    kept $(wc -l < out_ksweep/baseline.prune.in) SNPs, lowcov001 called: $baseline_called"

echo
printf "%-6s %10s %10s %10s %14s\n" K pool% retained kept lowcov001_gain
for K in 1 5 10 20 50 100; do
    python3 ../ldprio.py \
        --bfile out_ksweep/ksweep \
        --priority-samples out_ksweep/priority_K${K}.txt \
        --ld-samples out_ksweep/ld_samples.txt \
        --make-bed --out out_ksweep/K${K} --plink "$PLINK" \
        > out_ksweep/K${K}_run.log 2>&1 || true

    pool_line=$(grep "priority pool:" out_ksweep/K${K}_run.log)
    pool_pct=$(echo "$pool_line" | grep -oE '\([0-9.]+%' | tr -d '(%')
    retained=$(grep "^\[5/5\]" out_ksweep/K${K}_run.log | grep -oE '\([0-9.]+%' | tr -d '(%')
    kept=$(wc -l < out_ksweep/K${K}.snplist)

    $PLINK --bfile out_ksweep/K${K} --missing --allow-no-sex \
           --out out_ksweep/K${K}_miss >/dev/null 2>&1
    ldprio_called=$(awk '$2=="lowcov001"{print $5-$4}' out_ksweep/K${K}_miss.imiss)
    gain=$(awk -v b="$baseline_called" -v c="$ldprio_called" 'BEGIN{printf "%.1fx", c/b}')

    printf "%-6s %9s%% %9s%% %10s %14s\n" "$K" "$pool_pct" "$retained" "$kept" "$gain"
done
