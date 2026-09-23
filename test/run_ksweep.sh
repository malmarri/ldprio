#!/bin/bash
# run_ksweep.sh -- reproduces the README's "Scaling to more priority
# samples" tables: plain plink --indep-pairwise (LD from moderns) vs. ldprio
# with each --weighting at K = 1, 5, 10, 20, 50, 100 nominated low-coverage
# samples, with every
# sample at ~4% coverage ("uniform") and with coverage varying 0.5-20%
# ("varying"). Per K it reports, over the K nominated samples, the mean
# share of each sample's own calls that survive and the fewest SNPs any one
# of them keeps.
#
# Requires plink 1.9 on PATH. Run from this directory.
set -euo pipefail
cd "$(dirname "$0")"

PLINK=${PLINK:-plink}
command -v "$PLINK" >/dev/null || { echo "FAIL: plink not on PATH"; exit 1; }

# summarise <full.imiss> <pruned.imiss> <priority list> -> "mean% fewest"
summarise() {
    awk 'FNR==1{f++} f==1&&FNR>1{full[$2]=$5-$4} f==2&&FNR>1{kept[$2]=$5-$4}
         f==3{r=kept[$1]/full[$1]; s+=r; n++; if(n==1||kept[$1]<m)m=kept[$1]}
         END{printf "%.1f%% %d", 100*s/n, m}' "$1" "$2" "$3"
}

for mode in uniform varying; do
    d=out_ksweep/$mode
    rm -rf "$d" && mkdir -p "$d"
    (cd "$d" && python3 ../../make_ksweep_data.py "$mode" >/dev/null)
    $PLINK --file "$d/ksweep" --make-bed --out "$d/ksweep" >/dev/null 2>&1
    $PLINK --bfile "$d/ksweep" --missing --allow-no-sex --out "$d/full" >/dev/null 2>&1

    awk '{print "MODERN\t"$1}' "$d/ld_samples.txt" > "$d/ld.keep"
    $PLINK --bfile "$d/ksweep" --keep "$d/ld.keep" \
           --indep-pairwise 200 25 0.4 --allow-no-sex \
           --out "$d/baseline" >/dev/null 2>&1
    $PLINK --bfile "$d/ksweep" --extract "$d/baseline.prune.in" --missing \
           --allow-no-sex --out "$d/baseline_miss" >/dev/null 2>&1

    echo "=== $mode coverage ==="
    printf "%-5s %7s %14s %14s %14s\n" K pool% plink inverse fair
    printf "%-5s %7s %14s %14s %14s\n" "" "" "mean fewest" "mean fewest" "mean fewest"
    for K in 1 5 10 20 50 100; do
        read -r bmean bmin <<< "$(summarise "$d/full.imiss" "$d/baseline_miss.imiss" "$d/priority_K${K}.txt")"
        cols=()
        for weighting in inverse fair; do
            o=$d/K${K}.$weighting
            python3 ../ldprio.py \
                --bfile "$d/ksweep" \
                --priority-samples "$d/priority_K${K}.txt" \
                --ld-samples "$d/ld_samples.txt" --weighting "$weighting" \
                --out "$o" --plink "$PLINK" > "$o.run.log" 2>&1
            $PLINK --bfile "$d/ksweep" --extract "$o.snplist" --missing \
                   --allow-no-sex --out "$o.miss" >/dev/null 2>&1
            read -r lmean lmin <<< "$(summarise "$d/full.imiss" "$o.miss.imiss" "$d/priority_K${K}.txt")"
            cols+=("$lmean $lmin")
        done
        pool=$(grep -oE '\([0-9.]+% of panel' "$d/K${K}.inverse.run.log" | grep -oE '[0-9.]+%')
        printf "%-5s %7s %14s %14s %14s\n" "$K" "$pool" "$bmean $bmin" \
               "${cols[0]}" "${cols[1]}"
    done
    echo
done
