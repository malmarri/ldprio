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
including them hides real LD and leaves the panel under-pruned. --ld-samples
is required, so the cohort's pseudo-haploid samples are never used by
default. If a priority sample also appears in --ld-samples, a warning is
printed -- that's the exact mistake the above guards against.

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
import datetime
import heapq
import html
import operator
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict

__version__ = "1.5"

COMMON_PLINK_FLAGS = ["--allow-no-sex", "--allow-extra-chr"]
POOL_WARN_FRACTION = 0.70
MIN_LD_SAMPLES = 30
START = time.time()
MESSAGES = []


def emit(text):
    """Print a warning/note to stderr and keep it for the HTML report."""
    print(text, file=sys.stderr)
    MESSAGES.append(" ".join(text.split()))


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


def read_priority_bed(prefix, n_variants):
    """Sample IDs and the packed genotype bytes of a SNP-major .bed subset:
    (iids, n_samples, bytes_per_snp, data)."""
    with open(prefix + ".fam") as fh:
        iids = [ln.split()[1] for ln in fh if ln.strip()]
    n = len(iids)
    rec = (n + 3) // 4
    with open(prefix + ".bed", "rb") as fh:
        raw = fh.read()
    if raw[:3] != b"\x6c\x1b\x01":
        sys.exit(f"ERROR: {prefix}.bed is not a SNP-major PLINK 1 .bed")
    data = raw[3:]
    if len(data) != rec * n_variants:
        sys.exit(f"ERROR: {prefix}.bed has {len(data)} genotype bytes, "
                 f"expected {rec * n_variants}")
    return iids, n, rec, data


def _slots(n, p, b):
    """Samples with a call in byte b of a SNP record's byte p. The .bed packs
    4 samples per byte, 2 bits each, low bits first; 0b01 = missing; padding
    slots past the last sample are excluded."""
    return [j for j in range(4) if 4 * p + j < n and (b >> (2 * j)) & 3 != 1]


def per_sample_calls(n, rec, data, idx=None):
    """Called SNPs per sample, over every SNP or only the SNP records `idx`."""
    counts = [0] * n
    for p in range(rec):
        col = data[p::rec]
        for b, c in Counter(col if idx is None else (col[k] for k in idx)).items():
            for j in _slots(n, p, b):
                counts[4 * p + j] += c
    return counts


REPORT_CSS = """
:root{--bg:#fff;--card:#f6f8fa;--text:#1f2328;--muted:#59636e;--line:#d1d9e0;
--track:#e6eaee;--accent:#0969da;--good:#1a7f37;--bad:#cf222e;--warn:#9a6700}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#0d1117;
--card:#151b23;--text:#e6edf3;--muted:#9198a1;--line:#3d444d;--track:#262c36;
--accent:#4493f8;--good:#3fb950;--bad:#f85149;--warn:#d29922}}
:root[data-theme="dark"]{--bg:#0d1117;--card:#151b23;--text:#e6edf3;--muted:#9198a1;
--line:#3d444d;--track:#262c36;--accent:#4493f8;--good:#3fb950;--bad:#f85149;--warn:#d29922}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);
font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:920px;margin:0 auto;padding:28px 16px 48px}
h1{font-size:24px;margin:0 0 4px}
h2{font-size:17px;margin:34px 0 4px}
.sub,.hint{color:var(--muted);font-size:13px;margin:0 0 12px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin-top:18px}
.tile{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.tile .k{color:var(--muted);font-size:12px;text-transform:uppercase;letter-spacing:.06em}
.tile .v{font-size:26px;font-weight:600;font-variant-numeric:tabular-nums}
.tile .s{color:var(--muted);font-size:13px}
.good{color:var(--good)}.bad{color:var(--bad)}
.row{display:grid;grid-template-columns:minmax(80px,150px) 1fr minmax(120px,auto);
gap:10px;align-items:center;padding:3px 0}
.row.wide{grid-template-columns:minmax(120px,250px) 1fr minmax(70px,auto)}
.lab{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:13px}
.track{height:14px;background:var(--track);border-radius:7px;overflow:hidden}
.fill{height:100%;background:var(--accent);border-radius:7px;min-width:2px}
.fill.zero{min-width:0}
.num{font-size:13px;text-align:right;color:var(--muted);font-variant-numeric:tabular-nums}
.hist{display:flex;gap:6px;align-items:flex-end;height:170px;margin:12px 0 4px}
.bin{flex:1;display:flex;flex-direction:column;justify-content:flex-end;align-items:center;height:100%}
.hb{width:100%;background:var(--accent);border-radius:4px 4px 0 0;min-height:2px}
.hn{font-size:12px;color:var(--muted);margin-bottom:2px;font-variant-numeric:tabular-nums}
.hl{font-size:11px;color:var(--muted);margin-top:4px}
table{border-collapse:collapse;width:100%;font-size:13px;margin-top:8px}
th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:500}td.r,th.r{text-align:right;font-variant-numeric:tabular-nums}
dl{display:grid;grid-template-columns:max-content 1fr;gap:4px 18px;margin:8px 0;font-size:14px}
dt{color:var(--muted)}dd{margin:0;word-break:break-word}
code{font:12px ui-monospace,SFMono-Regular,Menlo,monospace;word-break:break-all}
ul{padding-left:20px;margin:8px 0}li{margin:5px 0;font-size:14px}
li.w{color:var(--warn)}
footer{color:var(--muted);font-size:12px;margin-top:40px}
"""


def write_report(path, r):
    """Write a self-contained HTML report (no scripts, no external files)."""
    esc = html.escape
    n = lambda x: f"{x:,}"
    pct = lambda a, b: f"{100 * a / b:.1f}%" if b else "n/a"

    samples = sorted(r["samples"], key=lambda t: (t[1], t[0]))
    usable = [t for t in samples if t[1] > 0]
    mean_share = (sum(k / a for _, a, k in usable) / len(usable)) if usable else 0
    ok = r["converged"]

    def tile(k, v, sub, cls=""):
        return (f'<div class="tile"><div class="k">{k}</div>'
                f'<div class="v {cls}">{v}</div><div class="s">{sub}</div></div>')

    tiles = "".join([
        tile("Final panel", n(r["final"]),
             f'of {n(r["n_snps"])} autosomal SNPs'),
        tile("Priority samples", n(len(samples)), f'{esc(r["weighting"])} weighting'),
        tile("Average share kept", f"{100 * mean_share:.0f}%",
             "of each priority sample's own SNPs"),
        tile("LD check", "Passed" if ok else "Not verified",
             f'no pair above r&sup2; {r["r2"]} within {r["window"]} variants'
             if ok else "raise --max-cleanup", "good" if ok else "bad"),
    ])

    if len(samples) <= 40:
        rows = []
        for iid, avail, kept in samples:
            share = kept / avail if avail else 0
            cls = "fill zero" if kept == 0 else "fill"
            rows.append(
                f'<div class="row"><div class="lab" title="{esc(iid)}">{esc(iid)}</div>'
                f'<div class="track"><div class="{cls}" style="width:{100 * share:.1f}%"></div></div>'
                f'<div class="num">{n(kept)} of {n(avail)} ({pct(kept, avail)})</div></div>')
        per_sample = "".join(rows)
        per_sample_hint = ('<p class="hint">Bar = share of the sample\'s own SNPs '
                           '(those it has calls at in the input panel) that '
                           'survived pruning.</p>')
    else:
        bins = [0] * 10
        for _, a, k in usable:
            bins[min(int(10 * k / a), 9)] += 1
        top = max(bins) or 1
        hist = "".join(
            f'<div class="bin"><div class="hn">{c}</div>'
            f'<div class="hb" style="height:{100 * c / top * 0.78:.1f}%"></div>'
            f'<div class="hl">{10 * i}&ndash;{10 * i + 10}%</div></div>'
            for i, c in enumerate(bins))
        worst = sorted(samples, key=lambda t: (t[2], t[0]))[:15]
        trs = "".join(
            f'<tr><td>{esc(i)}</td><td class="r">{n(a)}</td><td class="r">{n(k)}</td>'
            f'<td class="r">{pct(k, a)}</td></tr>' for i, a, k in worst)
        per_sample_hint = ('<p class="hint">Number of priority samples by share '
                           'of their own SNPs (those they have calls at in the '
                           'input panel) that survived pruning.</p>')
        per_sample = (
            f'<div class="hist">{hist}</div>'
            '<h3 style="font-size:14px;margin:22px 0 0">15 samples with the fewest SNPs kept</h3>'
            '<table><tr><th>Sample</th><th class="r">SNPs with calls</th>'
            f'<th class="r">SNPs kept</th><th class="r">Share kept</th></tr>{trs}</table>')

    steps = [("Autosomal SNPs in input", r["n_snps"], "")]
    steps.append(("After first LD pass", r["greedy"], ""))
    for i, pairs, removed, keep in r["passes"]:
        steps.append((f"Cleanup pass {i}", keep,
                      f'{n(pairs)} residual pairs, removed {n(removed)}'))
    steps.append(("Final panel", r["final"], ""))
    funnel = "".join(
        f'<div class="row wide"><div class="lab">{esc(lab)}</div>'
        f'<div class="track"><div class="fill" style="width:{100 * v / r["n_snps"]:.1f}%"></div></div>'
        f'<div class="num">{n(v)}</div></div>'
        + (f'<div class="hint" style="margin:-2px 0 6px 0">{extra}</div>' if extra else "")
        for lab, v, extra in steps)

    settings = [
        ("Weighting", r["weighting"]),
        ("LD samples", n(r["n_ld"])),
        ("r&sup2; threshold / window", f'{r["r2"]} / {r["window"]} variants'),
        ("Priority pool", f'{n(r["pool"])} SNPs ({pct(r["pool"], r["n_snps"])} of panel)'),
        ("Pool SNPs retained", f'{n(r["prio_final"])} ({pct(r["prio_final"], r["pool"])})'),
        ("Non-autosomal SNPs dropped", n(r["n_input"] - r["n_snps"])),
        ("Samples in input", n(r["n_fam"])),
    ]
    if r["tv_share"] is not None:
        settings.append(("Transversions in output", pct(*r["tv_share"])))
    settings.append(("Run time", f'{r["runtime"]:.1f} s' if r["runtime"] < 10
                     else f'{r["runtime"]:.0f} s'))
    dl = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in settings)

    def msg(m):
        if m.startswith("WARNING:"):
            return f'<li class="w">{esc(m[8:].strip())}</li>'
        return f'<li>{esc(m[5:].strip() if m.startswith("note:") else m)}</li>'

    msgs = "".join(msg(m) for m in r["messages"]) \
        or '<li class="good">No warnings.</li>'

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ldprio report</title><style>{REPORT_CSS}</style></head><body><main>
<h1>ldprio report</h1>
<p class="sub">v{esc(r["version"])} &middot; {esc(r["when"])} &middot; <code>{esc(r["out"])}</code></p>
<div class="tiles">{tiles}</div>
<h2>What each priority sample kept</h2>
{per_sample_hint}{per_sample}
<h2>Panel size through the run</h2>
<p class="hint">Bars are relative to the autosomal SNPs in the input.</p>
{funnel}
<h2>Run details</h2>
<dl>{dl}</dl>
<h2>Warnings and notes</h2>
<ul>{msgs}</ul>
<h2>Command</h2>
<p><code>{esc(r["command"])}</code></p>
<footer>Generated by ldprio {esc(r["version"])}. Self-contained file: no scripts or external resources.</footer>
</main></body></html>
"""
    with open(path, "w") as fh:
        fh.write(doc)


def compute_priority(plink, base, priority_keep_path, variants, work,
                     weighting):
    """Priority pool from the --priority-samples' calls, keyed by SNP and
    covering only SNPs at least one of them has called. For "inverse" the
    value is sum of 1/(sample's total called SNPs) over the samples calling
    it; for "fair" it is the list of those samples' indices."""
    prefix = os.path.join(work, "priority")
    run([plink, "--bfile", base, "--keep", priority_keep_path, "--make-bed",
         *COMMON_PLINK_FLAGS, "--out", prefix], "priority subset")

    iids, n, rec, data = read_priority_bed(prefix, len(variants))

    def called_slots(p, b):
        return _slots(n, p, b)

    coverage = [0] * n
    for p in range(rec):
        for b, count in Counter(data[p::rec]).items():
            for j in called_slots(p, b):
                coverage[4 * p + j] += count

    empty = [iids[i] for i in range(n) if coverage[i] == 0]
    if empty:
        preview = ", ".join(empty[:5]) + ("..." if len(empty) > 5 else "")
        emit(f"      WARNING: {len(empty)} priority sample(s) have no calls "
              f"in this panel: {preview}")

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


TRANSITIONS = ({"A", "G"}, {"C", "T"})


def read_transversions(bim_path):
    """IDs of biallelic A/C/G/T SNPs whose alleles are not a transition
    (A<->G, C<->T) pair; indels and unknown alleles are never included."""
    tv = set()
    with open(bim_path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 6:
                continue
            a, b = f[4].upper(), f[5].upper()
            if (len(a) == len(b) == 1 and a != b and a in "ACGT" and b in "ACGT"
                    and {a, b} not in TRANSITIONS):
                tv.add(sys.intern(f[1]))
    return tv


def greedy_select(order, adj, priority, weighting, transversions=frozenset()):
    """Maximal independent set over the LD graph in `adj`: keep a SNP, block
    its LD neighbours. Priority-pool SNPs go first -- highest weight first
    ("inverse"), or by turns to the sample with the fewest SNPs kept so far
    ("fair") -- then the rest in `order`'s order. SNPs in `transversions`
    win ties within the pool and go first among the rest. Returns kept SNPs
    in `order`'s order."""
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
            # popped from the end: most-shared first, then transversions,
            # then genomic order
            q.sort(key=lambda snp: (len(priority[snp]), snp in transversions,
                                    -index[snp]))
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
                      key=lambda snp: (-priority[snp], snp not in transversions))
        for snp in pool:
            if snp not in blocked:
                keep(snp)

    rest = [snp for snp in order if snp not in priority]
    if transversions:
        rest.sort(key=lambda snp: snp not in transversions)
    for snp in rest:
        if snp not in blocked:
            keep(snp)
    return [s for s in order if s in kept_set]


def r2_and_select(plink, base, snplist_path, ld_keep, window, r2,
                   priority, weighting, transversions, work, tag):
    """One priority-aware LD-resolution pass: --r2 on the current surviving
    set, then greedy-select over whatever conflicts it finds. Returns
    (new_snplist_path, n_pairs_found, n_removed)."""
    w = lambda name: os.path.join(work, name)
    with open(snplist_path) as fh:
        if sum(1 for ln in fh if ln.strip()) < 2:
            return snplist_path, 0, 0
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

    kept = greedy_select(order, adj, priority, weighting, transversions)
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
    ap.add_argument("--ld-samples", required=True,
                    help="samples LD is estimated from (use high-quality "
                         "samples with reliable diploid calls, modern or "
                         "high-coverage ancient; at least 30)")
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
    ap.add_argument("--prefer-transversions", action="store_true",
                    help="among otherwise equal SNPs, keep transversions "
                         "over transitions (C<->T, G<->A), which ancient-DNA "
                         "damage can mimic. Changes which SNP is kept only "
                         "where the choice barely matters for the priority "
                         "samples")
    ap.add_argument("--report", action="store_true",
                    help="also write an HTML report of the results "
                         "(<out>.report.html)")
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
    ver = subprocess.run([args.plink, "--version"], capture_output=True,
                         text=True).stdout
    if not ver.startswith("PLINK v1."):
        sys.exit(f"ERROR: ldprio needs PLINK 1.9, but '{args.plink}' reports "
                 f"'{ver.strip()[:60]}'. PLINK 2 is not supported.")
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
        emit(f"      WARNING: priority pool is {pool_frac:.0%} of the panel. "
              f"Most blocks are now contested between nominated samples, so "
              f"each keeps a smaller share of its sites -- consider "
              f"nominating only the samples that actually need it.")

    # ---- LD samples ------------------------------------------------------------
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
        emit(f"      WARNING: {len(overlap)} sample(s) appear in both "
              f"--priority-samples and --ld-samples: {preview}\n"
              f"      If these are pseudo-haploid ancient samples this is "
              f"exactly the LD-distortion case --ld-samples exists to "
              f"avoid (see LD SOURCE in --help).")
    if len(ld_samples) < MIN_LD_SAMPLES:
        emit(f"      WARNING: LD is estimated from only {len(ld_samples)} "
              f"sample(s). With fewer than {MIN_LD_SAMPLES}, r2 estimates are "
              f"noisy: unrelated SNPs can pass the threshold by chance and "
              f"linked SNPs can be missed, so the pruned panel may be smaller "
              f"or less independent than intended. Continuing anyway; "
              f"consider a larger --ld-samples set.")

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
            if len(f) >= 5 and f[4] in ("0", "NA"):
                monomorphic.add(f[1])
    unevaluated_priority = [s for s in priority if s in monomorphic]
    if unevaluated_priority:
        frac = len(unevaluated_priority) / len(priority)
        emit(f"      note: {len(unevaluated_priority):,}/{len(priority):,} "
              f"({frac:.1%}) priority-pool SNPs are monomorphic in --ld-samples "
              f"and were never evaluated for LD (always kept; may still be "
              f"correlated with each other among the samples that matter for "
              f"them, just invisible to the moderns-only estimate)")

    # ---- initial LD graph + greedy prioritised selection --------------------
    if n_snps >= 2:
        run([args.plink, "--bfile", base, *ld_keep, "--make-founders", "--r2",
             "--ld-window", str(args.window),
             "--ld-window-kb", "999999",
             "--ld-window-r2", str(args.r2),
             *COMMON_PLINK_FLAGS, "--out", w("ld")], "plink --r2")
        adj, n_pairs = read_ld_graph(w("ld.ld"))
    else:
        adj, n_pairs = defaultdict(list), 0

    transversions = (read_transversions(base + ".bim")
                     if args.prefer_transversions else frozenset())
    if args.prefer_transversions and not transversions:
        emit("      WARNING: --prefer-transversions found no A/C/G/T "
              "transversion SNPs in the .bim (alleles missing or coded 0?); "
              "the option has no effect.")
    kept = greedy_select(variants, adj, priority, args.weighting,
                         transversions)
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
    passes = []
    for i in range(1, args.max_cleanup + 1):
        current_snplist, n_pairs_i, n_removed = r2_and_select(
            args.plink, base, current_snplist, ld_keep, args.window,
            args.r2, priority, args.weighting, transversions, work,
            f"chk{i}")
        n_keep = sum(1 for _ in open(current_snplist))
        passes.append((i, n_pairs_i, n_removed, n_keep))
        print(f"      cleanup pass {i}: {n_pairs_i:,} residual pairs, "
              f"removed {n_removed:,} -> {n_keep:,} SNPs", file=sys.stderr)
        if n_pairs_i == 0:
            converged = True
            break

    if not converged:
        emit(f"WARNING: not converged after {args.max_cleanup} cleanup "
              f"passes; residual LD may remain. Raise --max-cleanup.")

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
    if transversions:
        n_tv = sum(1 for s in final if s in transversions)
        print(f"      transversions: {n_tv:,}/{len(final):,} kept SNPs "
              f"({n_tv / len(final):.1%})", file=sys.stderr)
    print(f"      wrote {args.out}.snplist", file=sys.stderr)

    if args.make_bed:
        run([args.plink, "--bfile", base, "--extract", args.out + ".snplist",
             "--make-bed", *COMMON_PLINK_FLAGS, "--out", args.out],
            "final --make-bed")
        print(f"      wrote {args.out}.{{bed,bim,fam}}", file=sys.stderr)

    if args.report:
        iids, n_prio, rec, data = read_priority_bed(w("priority"), n_snps)
        avail = per_sample_calls(n_prio, rec, data)
        kept_idx = [k for k, snp in enumerate(variants) if snp in final_set]
        kept_calls = per_sample_calls(n_prio, rec, data, kept_idx)
        n_tv = sum(1 for s in final if s in transversions)
        write_report(args.out + ".report.html", dict(
            version=__version__, out=args.out, command=" ".join(sys.argv),
            when=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            runtime=time.time() - START, weighting=args.weighting,
            r2=args.r2, window=args.window, n_input=n_input, n_snps=n_snps,
            n_fam=n_fam, n_ld=len(ld_samples), pool=len(priority),
            prio_final=n_prio_final, greedy=len(kept), passes=passes,
            final=len(final), converged=converged,
            samples=list(zip(iids, avail, kept_calls)),
            tv_share=(n_tv, len(final)) if transversions else None,
            messages=list(MESSAGES)))
        print(f"      wrote {args.out}.report.html", file=sys.stderr)

    return converged


if __name__ == "__main__":
    main()
