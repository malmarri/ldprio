#!/bin/bash
# run_test.sh -- end-to-end test for ldprio on the synthetic panel.
#
# Compares ldprio against the baseline everyone actually uses (plain
# `plink --indep-pairwise`, LD estimated on the same modern samples) and
# checks three things:
#   1. ldprio's output is genuinely LD-independent (plink finds nothing left
#      to prune)
#   2. the ancient samples retain substantially more SNPs than under the
#      baseline
#   3. the modern samples are not meaningfully worse off
#
# Requires plink 1.9 on PATH. Run from this directory.
set -euo pipefail
cd "$(dirname "$0")"

PLINK=${PLINK:-plink}
WINDOW=200 STEP=25 R2=0.4
command -v "$PLINK" >/dev/null || { echo "FAIL: plink not on PATH"; exit 1; }

rm -rf out && mkdir -p out

echo "=== baseline: plain plink --indep-pairwise (LD from moderns) ==="
awk '{print "MODERN\t"$1}' ld_samples.txt > out/ld.keep
$PLINK --bfile test_panel --keep out/ld.keep \
       --indep-pairwise $WINDOW $STEP $R2 \
       --allow-no-sex --out out/baseline >/dev/null 2>&1
$PLINK --bfile test_panel --extract out/baseline.prune.in \
       --make-bed --allow-no-sex --out out/baseline_panel >/dev/null 2>&1
echo "    kept $(wc -l < out/baseline.prune.in) SNPs"

echo "=== ldprio ==="
python3 ../ldprio.py \
    --bfile test_panel \
    --priority-samples priority_samples.txt \
    --ld-samples ld_samples.txt \
    --window $WINDOW --step $STEP --r2 $R2 \
    --make-bed --out out/ldprio --plink "$PLINK"

# ---- check 1: output really is LD-independent ------------------------------
$PLINK --bfile out/ldprio --keep out/ld.keep \
       --indep-pairwise $WINDOW $STEP $R2 \
       --allow-no-sex --out out/verify >/dev/null 2>&1
residual=$(wc -l < out/verify.prune.out)
echo
echo "CHECK 1: residual LD after ldprio = $residual SNPs"
[ "$residual" -eq 0 ] || { echo "FAIL: output is not LD-independent"; exit 1; }

# ---- checks 2 and 3: per-sample SNP counts ---------------------------------
for p in baseline_panel ldprio; do
    $PLINK --bfile out/$p --missing --allow-no-sex --out out/${p}_miss >/dev/null 2>&1
done

echo
printf "%-12s %12s %12s %10s\n" sample baseline ldprio change
fail=0
while read -r iid; do
    b=$(awk -v s="$iid" '$2==s{print $5-$4}' out/baseline_panel_miss.imiss)
    l=$(awk -v s="$iid" '$2==s{print $5-$4}' out/ldprio_miss.imiss)
    printf "%-12s %12s %12s %9s%%\n" "$iid" "$b" "$l" \
           "$(awk -v a="$b" -v c="$l" 'BEGIN{printf "%+.0f", (c-a)*100/a}')"
    [ "$l" -gt "$b" ] || fail=1
done < priority_samples.txt

mb=$(awk '$1=="MODERN"{s+=$5-$4; n++} END{printf "%d", s/n}' out/baseline_panel_miss.imiss)
ml=$(awk '$1=="MODERN"{s+=$5-$4; n++} END{printf "%d", s/n}' out/ldprio_miss.imiss)
printf "%-12s %12s %12s %9s%%\n" "modern(mean)" "$mb" "$ml" \
       "$(awk -v a="$mb" -v c="$ml" 'BEGIN{printf "%+.1f", (c-a)*100/a}')"

echo
echo "CHECK 2: every priority sample gained SNPs"
[ "$fail" -eq 0 ] || { echo "FAIL: a priority sample did not gain"; exit 1; }

echo "CHECK 3: modern samples within 10% of baseline"
awk -v a="$mb" -v c="$ml" 'BEGIN{exit !((c-a)/a > -0.10)}' \
    || { echo "FAIL: modern samples lost >10%"; exit 1; }

echo
echo "ALL CHECKS PASSED"
