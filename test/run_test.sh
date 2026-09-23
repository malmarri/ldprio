#!/bin/bash
# run_test.sh -- end-to-end test for ldprio on two synthetic panels:
#   test_panel   200 disjoint LD blocks of 10 (both tools keep 1 per block)
#   chain_panel  LD decaying along each chromosome (tools keep different
#                numbers of SNPs, so per-sample gains are also reported
#                relative to panel size)
#
# Compares ldprio (both --weighting modes) against plain
# `plink --indep-pairwise` (LD estimated on the same modern samples) and
# checks, per panel and mode:
#   1. ldprio's output is LD-independent (plink finds nothing left to prune)
#   2. every priority sample's share of the panel grows -- i.e. the gain is
#      not just a side effect of ldprio keeping more SNPs overall
#   3. modern samples keep at least 90% of the baseline's SNP count
#
# Requires plink 1.9 on PATH. Run from this directory.
set -euo pipefail
cd "$(dirname "$0")"

PLINK=${PLINK:-plink}
WINDOW=200 STEP=25 R2=0.4
command -v "$PLINK" >/dev/null || { echo "FAIL: plink not on PATH"; exit 1; }

rm -rf out && mkdir -p out
awk '{print "MODERN\t"$1}' ld_samples.txt > out/ld.keep

called() {  # called <imiss> <iid>
    awk -v s="$2" '$2==s{print $5-$4}' "$1"
}

test_panel() {
    local panel=$1 weighting=$2 o=out/$1.$2
    echo "############ $panel, --weighting $weighting ############"

    $PLINK --bfile "$panel" --keep out/ld.keep \
           --indep-pairwise $WINDOW $STEP $R2 \
           --allow-no-sex --out "$o.baseline" >/dev/null 2>&1
    $PLINK --bfile "$panel" --extract "$o.baseline.prune.in" \
           --make-bed --allow-no-sex --out "$o.baseline_panel" >/dev/null 2>&1

    python3 ../ldprio.py \
        --bfile "$panel" \
        --priority-samples priority_samples.txt \
        --ld-samples ld_samples.txt \
        --window $WINDOW --r2 $R2 --weighting "$weighting" \
        --make-bed --out "$o.ldprio" --plink "$PLINK"

    $PLINK --bfile "$o.ldprio" --keep out/ld.keep \
           --indep-pairwise $WINDOW $STEP $R2 \
           --allow-no-sex --out "$o.verify" >/dev/null 2>&1
    local residual nb nl
    residual=$(wc -l < "$o.verify.prune.out" | tr -d ' ')
    nb=$(wc -l < "$o.baseline.prune.in" | tr -d ' ')
    nl=$(wc -l < "$o.ldprio.snplist" | tr -d ' ')

    for p in baseline_panel ldprio; do
        $PLINK --bfile "$o.$p" --missing --allow-no-sex \
               --out "$o.${p}_miss" >/dev/null 2>&1
    done

    echo
    echo "panel size: plink $nb SNPs, ldprio $nl SNPs"
    printf "%-13s %8s %8s %8s %16s\n" sample plink ldprio change "share-of-panel"
    local fail=0 b l iid
    while read -r iid; do
        b=$(called "$o.baseline_panel_miss.imiss" "$iid")
        l=$(called "$o.ldprio_miss.imiss" "$iid")
        printf "%-13s %8s %8s %7s%% %15s\n" "$iid" "$b" "$l" \
            "$(awk -v a="$b" -v c="$l" 'BEGIN{printf "%+.0f", (c-a)*100/a}')" \
            "$(awk -v a="$b" -v c="$l" -v na="$nb" -v nc="$nl" \
                  'BEGIN{printf "%.1fx", (c/nc)/(a/na)}')"
        awk -v a="$b" -v c="$l" -v na="$nb" -v nc="$nl" \
            'BEGIN{exit !(c/nc > a/na)}' || fail=1
    done < priority_samples.txt

    local mb ml
    mb=$(awk '$1=="MODERN"{s+=$5-$4; n++} END{printf "%d", s/n}' "$o.baseline_panel_miss.imiss")
    ml=$(awk '$1=="MODERN"{s+=$5-$4; n++} END{printf "%d", s/n}' "$o.ldprio_miss.imiss")
    printf "%-13s %8s %8s %7s%%\n" "modern(mean)" "$mb" "$ml" \
        "$(awk -v a="$mb" -v c="$ml" 'BEGIN{printf "%+.1f", (c-a)*100/a}')"

    echo
    echo "CHECK 1: residual LD after ldprio = $residual SNPs"
    [ "$residual" -eq 0 ] || { echo "FAIL: output is not LD-independent"; exit 1; }
    echo "CHECK 2: every priority sample's share of the panel grew"
    [ "$fail" -eq 0 ] || { echo "FAIL: a priority sample did not gain"; exit 1; }
    echo "CHECK 3: modern samples keep >= 90% of baseline"
    awk -v a="$mb" -v c="$ml" 'BEGIN{exit !((c-a)/a > -0.10)}' \
        || { echo "FAIL: modern samples lost >10%"; exit 1; }
    echo
}

for panel in test_panel chain_panel; do
    for weighting in inverse fair; do
        test_panel "$panel" "$weighting"
    done
done
echo "ALL CHECKS PASSED"
