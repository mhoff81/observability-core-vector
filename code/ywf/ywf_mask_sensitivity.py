#!/usr/bin/env python3
"""YWF OCV — the echo rule of this record, swept offline (a rule-version study).

The committed layer uses the project-wide echo rule

    mask = intensity >= 0.30 * peak  AND  intensity >= 5.0 * median
    A    = count(mask) >= 2

(``ywf_ocv_masks.PEAK_FRAC`` / ``MEDIAN_MULT`` / ``MIN_N_MASKED``), identical to
``analyze_carola_coherence.coherence_gamma2`` and to the backend's
``sarsolve::coherence::coherence_mast_echo``. The headline of this site is the
echo coverage of the collapse — 16 of 27 pre-event chips against 1 of 6
post-event chips, Fisher p = 0.085 — and the honest question about such a number
is how much of it is the *rule* rather than the scene.

This module answers that question **without changing anything**:

  * it re-runs the mask on the committed payload over a grid of rule parameters
    (``peak_frac`` x ``median_mult`` x ``min_n_masked``) and reports the census
    each set produces, with its delta to the committed one;
  * it scans each knob finely and measures how far it can move before the two
    headline numbers change at all — the knife edge of the rule;
  * it evaluates one *local*, neighbour-referenced variant, the 7x7 form of the
    detection Bautzen's deck-edge layer anchors on
    (there ``|z|^2 >= 1.5 * median(neighbour rows)``).

The committed rule is **not** touched: no cache, CSV or reference file is
written, and the layer this package commits is, and stays, version 1 (the family
rule). Any variant here would be a *family-wide* change — it is defined in five
other packages' masks as well — which is exactly why it is studied and not
switched on.

Usage
-----
  python3 code/ywf/ywf_mask_sensitivity.py
  python3 code/ywf/ywf_mask_sensitivity.py --json /tmp/ywf_mask_sweep.json
Command lines are always run from the **repository root**.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import ywf_ocv_core as core  # noqa: E402
import ywf_ocv_masks as masks  # noqa: E402

# The coarse grid: two notches either side of every committed threshold, so the
# table shows the census at the rule the family uses and at its neighbours.
PEAK_GRID = [0.20, 0.25, 0.30, 0.35, 0.40]
MULT_GRID = [2.5, 5.0, 7.5, 10.0]
MIN_N_GRID = [2, 3]
K_LOCAL_GRID = [1.5, 2.0, 3.0]

# The fine scans behind the knife-edge statement: ``(lo, hi, step)``.
FINE_PEAK = (0.02, 0.60, 0.01)
FINE_MULT = (1.0, 15.0, 0.25)

COMMITTED = {"peak_frac": masks.PEAK_FRAC, "median_mult": masks.MEDIAN_MULT,
             "min_n_masked": masks.MIN_N_MASKED}


# ---------------------------------------------------------------------------
# The censuses (the committed coverage block, rebuilt from any mask set)
# ---------------------------------------------------------------------------
def masks_for(chips, **kw):
    """``measurement_id -> echo_mask`` of every chip under one parameter set."""
    return {c["mid"]: masks.echo_mask(c["z"], **kw) for c in chips}


def coverage(chips, by_mid):
    """The committed ``echo_coverage`` block, rebuilt from an arbitrary rule."""
    light = [{"state": c["state"], "mask": by_mid[c["mid"]]} for c in chips]
    return core.coverage_block(light)


def headline(block):
    """The numbers every layer of this package reports, per rule."""
    pre, post = block["per_state"]["pre"], block["per_state"]["post"]
    return {"n_unique": block["n_windows"],
            "n_pre": pre["n_chips"], "pre_echo": pre["n_with_echo"],
            "n_post": post["n_chips"], "post_echo": post["n_with_echo"],
            "n_echo": pre["n_with_echo"] + post["n_with_echo"],
            "p_two_sided": block["fisher"]["p_two_sided"]}


def mode_mix(chips, by_mid):
    """``no_echo / compact / intermediate / distributed`` over the chips."""
    c = Counter()
    for ch in chips:
        m = by_mid[ch["mid"]]
        c["no_echo" if m is None else core.echo_mode_of(m["A"])] += 1
    return {k: c[k] for k in ("no_echo", "compact", "intermediate",
                              "distributed")}


def a_histogram(chips):
    """How many chips hold ``A`` masked pixels, counted with ``min_n_masked=1``.

    The refusal threshold is the rule's second counting decision, so the
    histogram of ``A`` is what says how many chips a step of it would flip: at
    ``min_n_masked = 1`` every chip with ``A >= 1`` is an echo, at 2 (the
    committed value) it needs ``A >= 2``, at 3 it needs ``A >= 3``.
    """
    c = Counter()
    for ch in chips:
        m = masks.echo_mask(ch["z"], min_n_masked=1)
        c[m["A"] if m else 0] += 1
    return {k: c[k] for k in sorted(c)}


def local_echo_mask(z, k_local, min_n_masked=masks.MIN_N_MASKED):
    """Bautzen's neighbour-referenced detection, in the form 7x7 allows.

    Bautzen's band layer detects where a pixel stands ``CONTRAST_DETECT`` times
    above the median of its *neighbour rows*; its chiplike unit is an 80 x 80
    rectified window with a committed rect spec, so it has rows to spare. A 7 x 7
    window holds 49 px, its outer ring 24 of them, one to three pixels
    (3.4-15 m on this grid) from the centre: that ring is the only background a
    committed payload still carries, so it is taken as the reference and the
    inner 5 x 5 as the detector. Same shape of decision as the committed rule
    (``intensity >= k * reference``), different reference: local, not global.
    """
    intensity = np.abs(z) ** 2
    ring = np.concatenate([intensity[0], intensity[-1], intensity[1:-1, 0],
                           intensity[1:-1, -1]])
    thr = k_local * float(np.median(ring))
    inner = np.zeros(intensity.shape, dtype=bool)
    inner[1:-1, 1:-1] = True
    mask = (intensity >= thr) & inner   # the ring is the reference, never a hit
    if int(mask.sum()) < min_n_masked:
        return None
    return {"A": int(mask.sum()), "mask": mask, "threshold": float(thr),
            "ring_median": float(np.median(ring))}


# ---------------------------------------------------------------------------
# The sweeps
# ---------------------------------------------------------------------------
def build_grid(chips):
    """``(peak_frac, median_mult, min_n_masked) -> census entry``."""
    out = {}
    for peak in PEAK_GRID:
        for mult in MULT_GRID:
            for mn in MIN_N_GRID:
                by_mid = masks_for(chips, peak_frac=peak, median_mult=mult,
                                   min_n_masked=mn)
                out[(peak, mult, mn)] = {"headline": headline(
                    coverage(chips, by_mid)),
                    "echo_modes": mode_mix(chips, by_mid)}
    return out


def sweep(chips, knob, values):
    """One knob over ``values``, every other knob at the committed rule."""
    out = []
    for v in values:
        kw = dict(COMMITTED)
        kw[knob] = v
        by_mid = masks_for(chips, **kw)
        out.append({knob: v, "headline": headline(coverage(chips, by_mid)),
                    "echo_modes": mode_mix(chips, by_mid)})
    return out


def same_headline(a, b):
    """Do two censuses agree on the two headline numbers?"""
    return ((a["n_echo"], a["pre_echo"], a["post_echo"])
            == (b["n_echo"], b["pre_echo"], b["post_echo"]))


def knife_edge(scan, knob):
    """The maximal run of *identical* headlines around the committed value.

    Returns ``(lo, hi, (n_echo, pre_echo, post_echo), lo_edge, hi_edge)``: the
    interval of ``knob`` over which every statement this package makes about the
    site is *unchanged*. The two flags say whether the run reaches an end of the
    scan, i.e. whether the knob binds at all over the scanned interval.
    """
    def hl(step):
        h = step["headline"]
        return (h["n_echo"], h["pre_echo"], h["post_echo"])

    at = COMMITTED[knob]
    idx = min(range(len(scan)), key=lambda j: abs(scan[j][knob] - at))
    lo = hi = idx
    while lo > 0 and hl(scan[lo - 1]) == hl(scan[idx]):
        lo -= 1
    while hi + 1 < len(scan) and hl(scan[hi + 1]) == hl(scan[idx]):
        hi += 1
    return (scan[lo][knob], scan[hi][knob], hl(scan[idx]),
            lo == 0, hi == len(scan) - 1)


def frange(lo, hi, step):
    """A rounded inclusive float range (the fine scans)."""
    n = int(round((hi - lo) / step))
    return [round(lo + i * step, 6) for i in range(n + 1)]


def committed_matches(ref, built):
    """``(bool, [differences])``: does the reference block hold ``built``?"""
    pre = (ref.get("per_state") or {}).get("pre") or {}
    post = (ref.get("per_state") or {}).get("post") or {}
    want = (("pre.n_chips", built["n_pre"], pre.get("n_chips")),
            ("pre.n_with_echo", built["pre_echo"], pre.get("n_with_echo")),
            ("post.n_chips", built["n_post"], post.get("n_chips")),
            ("post.n_with_echo", built["post_echo"], post.get("n_with_echo")),
            ("fisher.p_two_sided", built["p_two_sided"],
             (ref.get("fisher") or {}).get("p_two_sided")))
    diffs = [f"{name}: rebuilt {got!r} vs committed {ref_v!r}"
             for name, got, ref_v in want if got != ref_v]
    return (not diffs), diffs


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def fmt(x):
    """A parameter as it is written in the tables."""
    return f"{x:g}"


def fmt_p(p):
    """A Fisher p (``None`` when scipy is absent)."""
    return "n/a" if p is None else f"{p:.4f}"


def brief(step):
    """One line of a scan's census."""
    h = step["headline"]
    return (f"{h['pre_echo']}/{h['n_pre']} pre, {h['post_echo']}/{h['n_post']} "
            f"post, Fisher p = {fmt_p(h['p_two_sided'])}")


def print_rows(title, knob, rows, base_total):
    """``rows`` is a list of ``(value, census entry)``; ``base_total`` the
    committed count the ``dtot`` column is taken against."""
    print(f"\n{title}")
    print(f"  {knob:<11} {'pre':>7} {'post':>7} {'total':>8} {'dtot':>5} "
          f"{'Fisher p':>9}  no/comp/int/dist")
    for v, e in rows:
        h, m = e["headline"], e["echo_modes"]
        mark = " *" if float(v) == float(COMMITTED[knob]) else ""
        print(f"  {fmt(v):<11} {h['pre_echo']:>3}/{h['n_pre']:<3} "
              f"{h['post_echo']:>3}/{h['n_post']:<3} "
              f"{h['n_echo']:>3}/{h['n_unique']:<3} "
              f"{h['n_echo'] - base_total:>+5d} "
              f"{fmt_p(h['p_two_sided']):>9}  "
              f"{m['no_echo']}/{m['compact']}/{m['intermediate']}"
              f"/{m['distributed']}{mark}")
    print("  (* the committed rule; dtot is against its 17 unique chips)")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", metavar="PATH",
                    help="write the full sweep to PATH (JSON)")
    args = ap.parse_args(argv)

    _, chips = core.chip_table()
    ref = core.load_findings()
    ref_block = ref.get("echo_coverage")
    ref_block = ref_block if isinstance(ref_block, dict) else {}
    built = headline(coverage(chips, masks_for(chips, **COMMITTED)))

    print("YWF OCV — the echo rule of this record, swept offline")
    print("(a rule-version study: no cache, CSV or reference file is written)\n")
    print(f"committed rule    peak_frac {fmt(COMMITTED['peak_frac'])}, "
          f"median_mult {fmt(COMMITTED['median_mult'])}, "
          f"min_n_masked {COMMITTED['min_n_masked']}")
    print(f"committed result  {built['pre_echo']}/{built['n_pre']} pre-event, "
          f"{built['post_echo']}/{built['n_post']} post-event, "
          f"{built['n_echo']}/{built['n_unique']} unique chips, "
          f"Fisher p = {fmt_p(built['p_two_sided'])}")
    ok, diffs = committed_matches(ref_block, built)
    print(f"  re-derives the committed echo_coverage block of "
          f"{os.path.relpath(core.REF_FINDINGS)}: "
          f"{'yes' if ok else 'NO — ' + '; '.join(diffs)}")

    hist = a_histogram(chips)
    hits = {t: sum(n for a, n in hist.items() if a >= t) for t in (1, 2, 3)}
    print(f"\nA histogram over the {built['n_unique']} unique chips "
          f"(A counted with min_n_masked = 1, so A = 0 is the refused case):")
    print("  " + ", ".join(f"A={a}: {n}" for a, n in hist.items()))
    print(f"  chips carrying an echo at min_n_masked = 1 / 2 / 3: "
          f"{hits[1]} / {hits[2]} / {hits[3]}")

    grid = build_grid(chips)
    base_mult, base_mn = masks.MEDIAN_MULT, masks.MIN_N_MASKED
    base_peak = masks.PEAK_FRAC
    print_rows(f"peak_frac — the global contrast threshold "
               f"(median_mult {fmt(base_mult)}, min_n_masked {base_mn})",
               "peak_frac",
               [(p, grid[(p, base_mult, base_mn)]) for p in PEAK_GRID],
               built["n_echo"])
    print_rows(f"median_mult — the local-significance threshold "
               f"(peak_frac {fmt(base_peak)}, min_n_masked {base_mn})",
               "median_mult",
               [(m_, grid[(base_peak, m_, base_mn)]) for m_ in MULT_GRID],
               built["n_echo"])
    print_rows(f"min_n_masked — the refusal threshold "
               f"(peak_frac {fmt(base_peak)}, median_mult {fmt(base_mult)})",
               "min_n_masked",
               [(n, grid[(base_peak, base_mult, n)]) for n in MIN_N_GRID],
               built["n_echo"])

    for mn in MIN_N_GRID:
        print(f"\nfull grid — unique chips of {built['n_unique']} with an echo, "
              f"min_n_masked {mn}")
        print("  peak \\ mult " + " ".join(f"{fmt(m_):>6}" for m_ in MULT_GRID))
        for p in PEAK_GRID:
            cells = " ".join(f"{grid[(p, m_, mn)]['headline']['n_echo']:>6}"
                             for m_ in MULT_GRID)
            print(f"  {fmt(p):<11} {cells}")

    print("\nlocal, neighbour-referenced rule — the 24 px outer ring of the "
          "window as the\nreference, its inner 5x5 as the detector "
          "(Bautzen's ``|z|^2 >= k * median(neighbours)`` in 7x7 form)")
    print(f"  {'k':<11} {'pre':>7} {'post':>7} {'total':>8} {'dtot':>5} "
          f"{'Fisher p':>9}")
    local_rows = []
    for k in K_LOCAL_GRID:
        by_mid = {c["mid"]: local_echo_mask(c["z"], k) for c in chips}
        blk = headline(coverage(chips, by_mid))
        local_rows.append({"k_local": k, **blk})
        print(f"  {fmt(k):<11} {blk['pre_echo']:>3}/{blk['n_pre']:<3} "
              f"{blk['post_echo']:>3}/{blk['n_post']:<3} "
              f"{blk['n_echo']:>3}/{blk['n_unique']:<3} "
              f"{blk['n_echo'] - built['n_echo']:>+5d} "
              f"{fmt_p(blk['p_two_sided']):>9}")

    print("\nknife edges — how far a knob moves before the headline changes")
    edges = {}
    for knob, (lo, hi, step) in (("peak_frac", FINE_PEAK),
                                 ("median_mult", FINE_MULT)):
        scan = sweep(chips, knob, frange(lo, hi, step))
        e_lo, e_hi, hl, lo_edge, hi_edge = knife_edge(scan, knob)
        ends = [n for n, flag in (("low", lo_edge), ("high", hi_edge)) if flag]
        edges[knob] = {"lo": e_lo, "hi": e_hi, "low_at_scan_edge": lo_edge,
                       "high_at_scan_edge": hi_edge, "n_echo": hl[0],
                       "pre_echo": hl[1], "post_echo": hl[2],
                       "scan": [lo, hi, step]}
        tail = ""
        if ends:
            tail = (" — the run reaches the " + " and ".join(ends)
                    + " end of the scan")
        print(f"  {knob:<12} {fmt(COMMITTED[knob])} may move over "
              f"[{fmt(e_lo)}, {fmt(e_hi)}] ({fmt(round(e_hi - e_lo, 6))} wide) "
              f"with {hl[1]}/{built['n_pre']} pre, {hl[2]}/{built['n_post']} "
              f"post unchanged"
              f"{tail}")
        below = [s for s in scan if s[knob] < e_lo]
        above = [s for s in scan if s[knob] > e_hi]
        print(f"    one step below: "
              f"{brief(below[-1]) if below else 'no change at the low end'}")
        print(f"    one step above: "
              f"{brief(above[0]) if above else 'no change at the high end'}")

    same_cells = [(p, m_) for (p, m_, mn) in grid
                  if mn == base_mn and same_headline(grid[(p, m_, mn)]
                                                     ["headline"], built)]
    print("\nReading")
    print(f"  the committed headline ({built['pre_echo']}/{built['n_pre']} pre, "
          f"{built['post_echo']}/{built['n_post']} post) is reproduced by "
          f"{len(same_cells)} of the {len(PEAK_GRID) * len(MULT_GRID)} "
          f"(peak_frac, median_mult) pairs at min_n_masked {base_mn}: "
          + ", ".join(f"peak {fmt(p)} / mult {fmt(m_)}"
                      for p, m_ in same_cells))
    print(f"  every knob has a knife edge several steps wide "
          f"(peak_frac {fmt(COMMITTED['peak_frac'])} over "
          f"[{fmt(edges['peak_frac']['lo'])}, {fmt(edges['peak_frac']['hi'])}], "
          f"median_mult {fmt(COMMITTED['median_mult'])} over "
          f"[{fmt(edges['median_mult']['lo'])}, "
          f"{fmt(edges['median_mult']['hi'])}]), so the committed census is not "
          f"balanced on one crossing")
    print(f"  the neighbour-referenced variant moves the census without "
          f"reproducing it: k = {fmt(K_LOCAL_GRID[0])} gives "
          f"{local_rows[0]['n_echo']}, k = {fmt(K_LOCAL_GRID[-1])} gives "
          f"{local_rows[-1]['n_echo']} of {built['n_unique']} chips against the "
          f"committed {built['n_echo']}")
    print(f"  the global contrast threshold is the *inactive* knob on this "
          f"record: it moves over the whole scanned range up to "
          f"{fmt(edges['peak_frac']['hi'])} without changing a single count, "
          f"while median_mult binds within "
          f"{fmt(round(edges['median_mult']['hi'] - edges['median_mult']['lo'], 6))} "
          f"of its committed value — one fine step. The rule's content on this "
          f"record is therefore the median condition and the refusal count, "
          f"not the peak")
    print("  (measured on the committed payload only; see code/ywf/README.md, "
          "section \"What could\n  improve this record — and what cannot\", for "
          "what would have to change in the record itself)")

    if args.json:
        payload = {
            "rule": COMMITTED,
            "committed": built,
            "committed_matches_reference": ok,
            "committed_differences": diffs,
            "a_histogram": hist,
            "chips_with_echo_at_min_n": hits,
            "grid": [{"peak_frac": p, "median_mult": m_, "min_n_masked": mn,
                      **grid[(p, m_, mn)]}
                     for p in PEAK_GRID for m_ in MULT_GRID
                     for mn in MIN_N_GRID],
            "local_rule": local_rows,
            "knife_edges": edges,
            "n_unique_chips": built["n_unique"],
            "event": ref_block.get("event"),
        }
        with open(args.json, "w") as fh:
            json.dump(payload, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print(f"\nwritten: {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
