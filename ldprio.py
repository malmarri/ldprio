#!/usr/bin/env python3
"""
ldprio -- LD pruning that protects the SNPs your low-coverage samples
actually cover.

THE PROBLEM
-----------
`plink --indep-pairwise` decides which SNP survives each LD block without
looking at per-sample missingness. When a block contains several mutually
redundant SNPs it keeps an arbitrary one -- which, for a low-coverage ancient
sample covering a sparse and effectively random subset of sites, is usually
not one of the few sites that sample has data for. The sample's already thin
data then gets thinned further for no analytical gain: the discarded SNP was
interchangeable with the kept one for every well-covered sample, but not for
the sparse one.

THE FIX
-------
Do the pruning as one genome-wide pass (so every LD relationship is still
evaluated exactly once), but bias *which* SNP wins each block toward sites
covered by the samples you nominate. Everyone else is unaffected -- they had
an equivalent marker either way. Blocks contested between nominated samples
are settled by --weighting:
  inverse (default)  each SNP scores the sum of 1/(total called SNPs) over
                     the nominated samples covering it, so the sparsest
                     samples win -- best when coverage varies between them.
  fair               nominated samples take turns: the one with the fewest
                     SNPs kept so far (sparsest first on ties) keeps its
                     next available site, preferring sites shared with more
                     nominated samples -- evens out SNPs kept per sample,
                     best when coverage is similar.

LD SOURCE
---------
LD is estimated only from the samples given by --ld-samples, which should be
high-quality samples with reliable diploid calls -- modern, or high-coverage
ancient. Pseudo-haploid ancient calls (a single random
read reported as a homozygote) systematically distort r2 downward, so
including them hides real LD and leaves the panel under-pruned. If you omit
--ld-samples the whole cohort is used and a warning is printed. If a
priority sample also appears in --ld-samples, a warning is printed -- that's
the exact mistake the above guards against.

METHOD
------
  1. priority pool := SNPs with a non-missing call in at least one
     --priority-samples sample, ranked per --weighting (see THE FIX)
  2. `plink --r2` on --ld-samples  ->  LD-conflict graph (pairs above --r2).
     SNPs monomorphic in --ld-samples never appear in this graph (undefined
     r2) and so always survive uncontested -- reported separately since this
     disproportionately affects priority-pool SNPs (rare-in-moderns sites
     are exactly what low-coverage ancient samples are likely to carry).
  3. greedy maximal independent set over that graph, visiting priority-pool
     SNPs first (in --weighting order, genomic order breaking ties), then
     the rest in genomic order (keep a SNP, block its LD neighbours)
  4. cleanup passes: re-run `plink --r2` on the surviving set (its window is
     re-anchored each pass, catching long-range pairs the first pass's
     static window missed) and resolve any residual conflicts with the same
     priority-first greedy rule -- never plink's own drop heuristic, so a
     priority SNP is never sacrificed for a non-priority one during cleanup.
     Iterated to convergence (typically 3-4 passes, each removing <0.5%).

The result is verified LD-independent within the --window variant-count
window used throughout (the same scope plink's own --indep-pairwise/--r2
guarantee, not a literal unbounded genome-wide claim): the final cleanup
pass finds zero residual pairs within that window. If it does not converge
within --max-cleanup passes, ldprio exits non-zero rather than silently
handing back a dirty panel.

EXAMPLE
-------
  ldprio.py \\
      --bfile data/panel_qc \\
      --priority-samples priority_samples.txt \\
      --ld-samples calculate_ld_samples.txt \\
      --make-bed \\
      --out data/panel_pruned

Sample list files: one sample per line, either "IID" or "FID<tab>IID".
Only autosomes (chromosomes 1-22) are pruned and written out.
"""

import argparse
import heapq
import operator
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

__version__ = "1.2"

COMMON_PLINK_FLAGS = ["--allow-no-sex", "--allow-extra-chr"]
POOL_WARN_FRACTION = 0.70


def run(cmd, label):
    """Run a command, surfacing stderr/stdout only if it fails."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"ERROR: {label} failed (exit {proc.returncode})\n"
                 f"  command: {' '.join(cmd)}\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def read_fam(path):
    """Return {iid: [(fid, iid), ...]} and the set of all (fid, iid) pairs.
    Skips blank/short lines rather than crashing on them (trailing newlines,
    hand-edited files)."""
    by_iid = defaultdict(list)
    fam_pairs = set()
    n_rows = 0
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 2:
                continue
            pair = (f[0], f[1])
            by_iid[f[1]].append(pair)
            fam_pairs.add(pair)
            n_rows += 1
    return by_iid, fam_pairs, n_rows


def resolve_samples(list_path, fam_by_iid, fam_pairs, label):
    """Accept 'IID' or 'FID IID' lines; return validated, order-preserving
    deduplicated (fid, iid) pairs. Every line is checked against the .fam --
    a typo'd FID/IID pair is reported, not silently dropped (plink's own
    --keep silently ignores non-matching lines as long as some match, which
    would otherwise hide a priority/ld sample being missed entirely)."""
    if not os.path.exists(list_path):
        sys.exit(f"ERROR: {label} file not found: {list_path}")

    keep, missing, ambiguous = [], [], []
    with open(list_path) as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            if len(f) >= 2:
                pair = (f[0], f[1])
                if pair in fam_pairs:
                    keep.append(pair)
                else:
                    missing.append(f"{f[0]} {f[1]}")
            else:
                matches = fam_by_iid.get(f[0], [])
                if len(matches) == 1:
                    keep.append(matches[0])
                elif len(matches) > 1:
                    ambiguous.append(f[0])
                else:
                    missing.append(f[0])

    if ambiguous:
        preview = ", ".join(ambiguous[:5]) + ("..." if len(ambiguous) > 5 else "")
        sys.exit(f"ERROR: {len(ambiguous)} sample(s) in {label} match more than "
                 f"one FID in the .fam (duplicate IID across families): "
                 f"{preview}\n  Specify 'FID<tab>IID' for these samples instead.")
    if missing:
        preview = ", ".join(missing[:5]) + ("..." if len(missing) > 5 else "")
        sys.exit(f"ERROR: {len(missing)} sample(s) in {label} not found in the "
                 f".fam: {preview}")
    if not keep:
        sys.exit(f"ERROR: {label} ({list_path}) is empty")
    return list(dict.fromkeys(keep))


def write_keep(keep, path):
    with open(path, "w") as fh:
        fh.writelines(f"{fid}\t{iid}\n" for fid, iid in keep)


def read_bim_variants(bim_path):
    """Ordered variant IDs from a .bim, skipping blank lines. Refuses
    duplicate (including repeated '.') variant IDs -- the LD graph is keyed
    by variant ID, so a duplicate silently merges two unrelated SNPs' LD
    neighbourhoods."""
    order = []
    seen = set()
    dupes = set()
    with open(bim_path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 2:
                continue
            vid = sys.intern(f[1])
            if vid in seen:
                dupes.add(vid)
            seen.add(vid)
            order.append(vid)
    if dupes:
        preview = ", ".join(list(dupes)[:5]) + ("..." if len(dupes) > 5 else "")
        sys.exit(f"ERROR: {len(dupes)} duplicate variant ID(s) in {bim_path}: "
                 f"{preview}\n  Variant IDs must be unique (rename unnamed '.' "
                 f"SNPs to chr:pos, e.g. plink --set-missing-var-ids @:#).")
    return order


def read_ld_graph(ld_path):
    """Parse a plink --r2 .ld file into an adjacency dict + pair count."""
    if not os.path.exists(ld_path):
        sys.exit(f"ERROR: expected LD output not found: {ld_path}\n"
                 f"  (plink returned success but produced no --r2 output -- "
                 f"check that --ld-samples lists enough samples with "
                 f"genotype calls in this panel)")
    adj = defaultdict(list)
    n_pairs = 0
    with open(ld_path) as fh:
        header = fh.readline()
        if not header:
            return adj, 0
        for line in fh:
            f = line.split()
            if len(f) < 6:
                continue
            a, b = sys.intern(f[2]), sys.intern(f[5])
            adj[a].append(b)
            adj[b].append(a)
            n_pairs += 1
    return adj, n_pairs


def compute_priority(plink, base, priority_keep_path, variants, work,
                     weighting):
    """Priority pool from the --priority-samples' calls, keyed by SNP and
    covering only SNPs at least one of them has called. For "inverse" the
    value is sum of 1/(sample's total called SNPs) over the samples calling
    it; for "fair" it is the list of those samples' indices."""
    prefix = os.path.join(work, "priority")
    run([plink, "--bfile", base, "--keep", priority_keep_path, "--make-bed",
         *COMMON_PLINK_FLAGS, "--out", prefix], "priority subset")

    with open(prefix + ".fam") as fh:
        iids = [ln.split()[1] for ln in fh if ln.strip()]
    n = len(iids)
    rec = (n + 3) // 4
    with open(prefix + ".bed", "rb") as fh:
        raw = fh.read()
    if raw[:3] != b"\x6c\x1b\x01":
        sys.exit(f"ERROR: {prefix}.bed is not a SNP-major PLINK 1 .bed")
    data = raw[3:]
    if len(data) != rec * len(variants):
        sys.exit(f"ERROR: {prefix}.bed has {len(data)} genotype bytes, "
                 f"expected {rec * len(variants)}")

    # .bed packs 4 samples per byte, 2 bits each, low bits first; 0b01 = missing.
    # Padding slots past the last sample are excluded.
    def called_slots(p, b):
        return [j for j in range(4)
                if 4 * p + j < n and (b >> (2 * j)) & 3 != 1]

    coverage = [0] * n
    for p in range(rec):
        for b, count in Counter(data[p::rec]).items():
            for j in called_slots(p, b):
                coverage[4 * p + j] += count

    empty = [iids[i] for i in range(n) if coverage[i] == 0]
    if empty:
        preview = ", ".join(empty[:5]) + ("..." if len(empty) > 5 else "")
        print(f"      WARNING: {len(empty)} priority sample(s) have no calls "
              f"in this panel: {preview}", file=sys.stderr)

    if weighting == "fair":
        slots = [[[4 * p + j for j in called_slots(p, b)] for b in range(256)]
                 for p in range(rec)]
        cover = {}
        for k, snp in enumerate(variants):
            off = k * rec
            samples = [i for p in range(rec) for i in slots[p][data[off + p]]]
            if samples:
                cover[snp] = samples
        return cover

    inv = [1.0 / c if c else 0.0 for c in coverage]
    weights = [0.0] * len(variants)
    for p in range(rec):
        table = [sum(inv[4 * p + j] for j in called_slots(p, b))
                 for b in range(256)]
        weights = list(map(operator.add, weights,
                           map(table.__getitem__, data[p::rec])))
    return {snp: wt for snp, wt in zip(variants, weights) if wt > 0}


def greedy_select(order, adj, priority, weighting):
    """Maximal independent set over the LD graph in `adj`: keep a SNP, block
    its LD neighbours. Priority-pool SNPs go first -- highest weight first
    ("inverse"), or by turns to the sample with the fewest SNPs kept so far
    ("fair") -- then the rest in `order`'s order. Returns kept SNPs in
    `order`'s order."""
    blocked, kept_set = set(), set()

    def keep(snp):
        kept_set.add(snp)
        blocked.update(adj.get(snp, ()))

    if weighting == "fair":
        index = {snp: k for k, snp in enumerate(order)}
        queues = defaultdict(list)
        for snp in order:
            for i in priority.get(snp, ()):
                queues[i].append(snp)
        for q in queues.values():
            # popped from the end: most-shared first, then genomic order
            q.sort(key=lambda snp: (len(priority[snp]), -index[snp]))
        n_kept = defaultdict(int)
        turns = [(0, len(q), i) for i, q in queues.items()]
        heapq.heapify(turns)
        while turns:
            n, size, i = heapq.heappop(turns)
            if n != n_kept[i]:
                heapq.heappush(turns, (n_kept[i], size, i))
                continue
            q = queues[i]
            while q and (q[-1] in blocked or q[-1] in kept_set):
                q.pop()
            if not q:
                continue
            snp = q.pop()
            keep(snp)
            for j in priority[snp]:
                n_kept[j] += 1
            heapq.heappush(turns, (n_kept[i], size, i))
    else:
        pool = sorted((snp for snp in order if snp in priority),
                      key=lambda snp: -priority[snp])
        for snp in pool:
            if snp not in blocked:
                keep(snp)

    for snp in order:
        if snp not in priority and snp not in blocked:
            keep(snp)
    return [s for s in order if s in kept_set]


def r2_and_select(plink, base, snplist_path, ld_keep, window, r2,
                   priority, weighting, work, tag):
    """One priority-aware LD-resolution pass: --r2 on the current surviving
    set, then greedy-select over whatever conflicts it finds. Returns
    (new_snplist_path, n_pairs_found, n_removed)."""
    w = lambda name: os.path.join(work, name)
    ld_path = w(f"{tag}.ld")
    if os.path.exists(ld_path):
        os.remove(ld_path)  # never let a stale .ld from an earlier pass/run survive a failed write
    run([plink, "--bfile", base, "--extract", snplist_path, *ld_keep,
         "--make-founders", "--r2",
         "--ld-window", str(window),
         "--ld-window-kb", "999999",
         "--ld-window-r2", str(r2),
         *COMMON_PLINK_FLAGS, "--out", w(f"{tag}")], f"{tag} --r2")

    adj, n_pairs = read_ld_graph(ld_path)
    order = [ln.strip() for ln in open(snplist_path) if ln.strip()]
    if n_pairs == 0:
        return snplist_path, 0, 0

    kept = greedy_select(order, adj, priority, weighting)
    n_removed = len(order) - len(kept)
    out_path = w(f"{tag}.snplist")
    with open(out_path, "w") as fh:
        fh.writelines(s + "\n" for s in kept)
    return out_path, n_pairs, n_removed


def main():
    ap = argparse.ArgumentParser(
        description="LD pruning that preserves SNPs covered by nominated "
                    "low-coverage samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Sample lists: one per line, 'IID' or 'FID<tab>IID'.")
    ap.add_argument("--bfile", required=True, help="input PLINK prefix")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--priority-samples", required=True,
                    help="samples whose covered SNPs should be favoured")
    ap.add_argument("--ld-samples",
                    help="samples LD is estimated from (use high-quality "
                         "samples with reliable diploid calls, modern or "
                         "high-coverage ancient; default: whole cohort, "
                         "with a warning)")
    ap.add_argument("--window", type=int, default=200,
                    help="window in variants (default: 200)")
    ap.add_argument("--step", type=int, default=25,
                    help="unused -- cleanup now re-runs --r2 + priority-"
                         "aware selection instead of plink --indep-pairwise, "
                         "which was the only consumer of --step. Kept for "
                         "CLI compatibility only; has no effect.")
    ap.add_argument("--r2", type=float, default=0.4,
                    help="r^2 threshold (default: 0.4)")
    ap.add_argument("--weighting", choices=("inverse", "fair"),
                    default="inverse",
                    help="how blocks contested between nominated samples are "
                         "settled: 'inverse' favours the sparsest samples "
                         "(best when coverage varies between them); 'fair' "
                         "evens out SNPs kept per sample (best when coverage "
                         "is similar). Default: inverse")
    ap.add_argument("--make-bed", action="store_true",
                    help="also write the pruned PLINK fileset")
    ap.add_argument("--max-cleanup", type=int, default=10,
                    help="cap on cleanup passes (default: 10)")
    ap.add_argument("--keep-intermediates", action="store_true",
                    help="keep working files (the .ld file can be ~0.5GB)")
    ap.add_argument("--plink", default="plink",
                    help="PLINK 1.9 executable (plink2 is not supported)")
    ap.add_argument("--version", action="version",
                    version=f"ldprio {__version__}")
    args = ap.parse_args()

    if not (0 < args.r2 <= 1):
        sys.exit(f"ERROR: --r2 must be in (0, 1], got {args.r2}")
    if args.window <= 0:
        sys.exit(f"ERROR: --window must be a positive integer, got {args.window}")
    if args.max_cleanup <= 0:
        sys.exit(f"ERROR: --max-cleanup must be a positive integer, got {args.max_cleanup}")

    if shutil.which(args.plink) is None:
        sys.exit(f"ERROR: '{args.plink}' not found on PATH")
    for ext in (".bed", ".bim", ".fam"):
        if not os.path.exists(args.bfile + ext):
            sys.exit(f"ERROR: {args.bfile}{ext} not found")
    if os.path.abspath(args.out) == os.path.abspath(args.bfile):
        sys.exit("ERROR: --out prefix cannot match --bfile prefix "
                 "(plink refuses identical input/output filesets)")

    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    work = tempfile.mkdtemp(prefix=os.path.basename(args.out) + ".",
                            suffix=".tmp", dir=out_dir)
    try:
        converged = prune(args, work)
    finally:
        if args.keep_intermediates:
            print(f"      intermediates kept in {work}/", file=sys.stderr)
        else:
            shutil.rmtree(work, ignore_errors=True)
    if not converged:
        sys.exit(1)


def prune(args, work):
    """Run the pipeline in `work`; returns whether cleanup converged."""
    w = lambda name: os.path.join(work, name)

    fam_by_iid, fam_pairs, n_fam = read_fam(args.bfile + ".fam")
    base = args.bfile

    autosomes = {str(c) for c in range(1, 23)}
    n_input = n_auto = 0
    with open(args.bfile + ".bim") as fh:
        for ln in fh:
            f = ln.split()
            if f:
                n_input += 1
                chrom = f[0].lower()
                n_auto += (chrom[3:] if chrom.startswith("chr") else chrom) in autosomes
    if n_auto == 0:
        sys.exit("ERROR: input panel has no autosomal SNPs (chromosomes 1-22)")
    run([args.plink, "--bfile", base, "--autosome", "--make-bed",
         *COMMON_PLINK_FLAGS, "--out", w("auto")], "autosome subset")
    base = w("auto")

    variants = read_bim_variants(base + ".bim")
    n_snps = len(variants)
    if n_snps == 0:
        sys.exit("ERROR: input panel has no autosomal SNPs (chromosomes 1-22)")
    print(f"ldprio {__version__}", file=sys.stderr)
    print(f"[1/5] input: {n_snps:,} autosomal SNPs "
          f"({n_input - n_snps:,} non-autosomal dropped), {n_fam:,} samples",
          file=sys.stderr)

    # ---- priority pool -------------------------------------------------------
    priority_samples = resolve_samples(args.priority_samples, fam_by_iid,
                                       fam_pairs, "--priority-samples")
    write_keep(priority_samples, w("priority.keep"))
    priority = compute_priority(args.plink, base, w("priority.keep"),
                                variants, work, args.weighting)
    if not priority:
        sys.exit("ERROR: priority samples cover no SNPs in this panel")
    pool_frac = len(priority) / n_snps
    print(f"[2/5] priority pool: {len(priority):,} SNPs covered by the "
          f"nominated samples ({pool_frac:.1%} of panel, "
          f"{len(priority_samples)} sample(s), {args.weighting} "
          f"weighting)", file=sys.stderr)
    if pool_frac > POOL_WARN_FRACTION:
        print(f"      WARNING: priority pool is {pool_frac:.0%} of the panel. "
              f"Most blocks are now contested between nominated samples, so "
              f"each keeps a smaller share of its sites -- consider "
              f"nominating only the samples that actually need it.",
              file=sys.stderr)

    # ---- LD samples ------------------------------------------------------------
    ld_keep = []
    if args.ld_samples:
        ld_samples = resolve_samples(args.ld_samples, fam_by_iid, fam_pairs,
                                     "--ld-samples")
        write_keep(ld_samples, w("ld.keep"))
        ld_keep = ["--keep", w("ld.keep")]
        print(f"[3/5] estimating LD from {len(ld_samples):,} samples",
              file=sys.stderr)
        overlap = set(priority_samples) & set(ld_samples)
        if overlap:
            preview = ", ".join(f"{fid}/{iid}" for fid, iid in sorted(overlap)[:5]) + \
                      ("..." if len(overlap) > 5 else "")
            print(f"      WARNING: {len(overlap)} sample(s) appear in both "
                  f"--priority-samples and --ld-samples: {preview}\n"
                  f"      If these are pseudo-haploid ancient samples this is "
                  f"exactly the LD-distortion case --ld-samples exists to "
                  f"avoid (see LD SOURCE in --help).", file=sys.stderr)
    else:
        print("[3/5] WARNING: --ld-samples not given, estimating LD from the "
              "whole cohort. If it contains pseudo-haploid ancient samples "
              "their calls will deflate r2 and leave the panel under-pruned.",
              file=sys.stderr)

    # SNPs monomorphic in --ld-samples have undefined r2 and plink's --r2
    # simply omits them from its output -- they never enter the LD graph and
    # so always survive, regardless of whether they're actually redundant
    # among the samples that matter for them. Checked directly via --freq
    # (NOT inferred from absence in the LD-pairs graph: a SNP can be fully
    # polymorphic and evaluated but simply have no partner above the r2
    # threshold nearby, which is the normal, correct state for most SNPs in
    # any real panel and must not be confused with "never evaluated").
    run([args.plink, "--bfile", base, *ld_keep, "--make-founders", "--freq",
         *COMMON_PLINK_FLAGS, "--out", w("modern_freq")],
        "--freq (monomorphic check)")
    monomorphic = set()
    with open(w("modern_freq.frq")) as fh:
        fh.readline()
        for line in fh:
            f = line.split()
            if len(f) >= 5 and f[4] == "0":
                monomorphic.add(f[1])
    unevaluated_priority = [s for s in priority if s in monomorphic]
    if unevaluated_priority:
        frac = len(unevaluated_priority) / len(priority)
        print(f"      note: {len(unevaluated_priority):,}/{len(priority):,} "
              f"({frac:.1%}) priority-pool SNPs are monomorphic in --ld-samples "
              f"and were never evaluated for LD (always kept; may still be "
              f"correlated with each other among the samples that matter for "
              f"them, just invisible to the moderns-only estimate)",
              file=sys.stderr)

    # ---- initial LD graph + greedy prioritised selection --------------------
    run([args.plink, "--bfile", base, *ld_keep, "--make-founders", "--r2",
         "--ld-window", str(args.window),
         "--ld-window-kb", "999999",
         "--ld-window-r2", str(args.r2),
         *COMMON_PLINK_FLAGS, "--out", w("ld")], "plink --r2")

    adj, n_pairs = read_ld_graph(w("ld.ld"))

    kept = greedy_select(variants, adj, priority, args.weighting)
    current_snplist = w("greedy.snplist")
    with open(current_snplist, "w") as fh:
        fh.writelines(s + "\n" for s in kept)
    n_prio_after_greedy = sum(1 for s in kept if s in priority)
    print(f"[4/5] LD pairs above r2>{args.r2}: {n_pairs:,} | greedy kept "
          f"{len(kept):,} SNPs | priority retained "
          f"{n_prio_after_greedy:,}/{len(priority):,} "
          f"({n_prio_after_greedy / len(priority):.1%})", file=sys.stderr)

    # ---- cleanup passes to convergence, still priority-aware -----------------
    # Re-runs --r2 on the surviving set each pass (its window is re-anchored
    # against the shrunken set, catching long-range pairs the first pass's
    # static window missed) and resolves conflicts with the same
    # priority-first rule as the initial pass -- never plink's own
    # --indep-pairwise heuristic, which has no notion of priority and could
    # drop a priority SNP in favour of a non-priority one.
    converged = False
    for i in range(1, args.max_cleanup + 1):
        current_snplist, n_pairs_i, n_removed = r2_and_select(
            args.plink, base, current_snplist, ld_keep, args.window,
            args.r2, priority, args.weighting, work, f"chk{i}")
        n_keep = sum(1 for _ in open(current_snplist))
        print(f"      cleanup pass {i}: {n_pairs_i:,} residual pairs, "
              f"removed {n_removed:,} -> {n_keep:,} SNPs", file=sys.stderr)
        if n_pairs_i == 0:
            converged = True
            break

    if not converged:
        print(f"WARNING: not converged after {args.max_cleanup} cleanup "
              f"passes; residual LD may remain. Raise --max-cleanup.",
              file=sys.stderr)

    # ---- outputs -----------------------------------------------------------
    final = [ln.strip() for ln in open(current_snplist) if ln.strip()]
    with open(args.out + ".snplist", "w") as fh:
        fh.writelines(s + "\n" for s in final)

    final_set = set(final)
    n_prio_final = sum(1 for s in final_set if s in priority)
    print(f"[5/5] final: {len(final):,} SNPs "
          f"({'LD-independent, verified' if converged else 'NOT verified'}) | "
          f"priority retained {n_prio_final:,}/{len(priority):,} "
          f"({n_prio_final / len(priority):.1%})", file=sys.stderr)
    print(f"      wrote {args.out}.snplist", file=sys.stderr)

    if args.make_bed:
        run([args.plink, "--bfile", base, "--extract", args.out + ".snplist",
             "--make-bed", *COMMON_PLINK_FLAGS, "--out", args.out],
            "final --make-bed")
        print(f"      wrote {args.out}.{{bed,bim,fam}}", file=sys.stderr)

    return converged


if __name__ == "__main__":
    main()
