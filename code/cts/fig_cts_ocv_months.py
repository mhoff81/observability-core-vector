#!/usr/bin/env python3
"""CTS OCV-Paper — the echo-mask dimensions from the first images to the month
after the collapse.

Question: what do the six dimensions of the **echo-mask vector**
``x = [gamma2, P, D, A, F, S]`` do along the record — from the *first Sentinel-1
images available for this footprint* (``F0``, the 2015 record start) to the
**first month after** the collapse (``2021-06-24``)? The A–F figures of
``fig_cts_ocv_paper.py`` pool the two states (``pre-collapse (healthy)`` vs.
``post-collapse``); this companion keeps the time axis and shows both the
long-baseline contrast the site's request asks for and the four one-month windows
that straddle the event.

Why ``F0`` is 2015, not 2014
----------------------------
Sentinel-1 data for **this** footprint do not exist in 2014: Sentinel-1A was
launched on 2014-04-03 and the suitable orbit (ASC rel-48 IW3) has its **first
acquisition on 2015-09-21** (``data/cts/cts_acquisition_census.json``; the
Step-0 census already queried from 2014-10-01 and found nothing usable earlier).
``F0`` is therefore the record-start epoch that carries the first images that
actually exist, and the report states that rather than inventing a 2014 window.

The five windows (end exclusive), in time order:

    F0    [2015-09-21, 2015-11-01)     the first images available (record start)
    M-3   [2021-03-24, 2021-04-24)     pre-collapse (healthy)
    M-2   [2021-04-24, 2021-05-24)     pre-collapse (healthy)
    M-1   [2021-05-24, 2021-06-24)     pre-collapse (healthy)
    M+1   [2021-06-24, 2021-07-24)     post-collapse (the month after the collapse)

Approach: this is a **standalone companion** of ``fig_cts_ocv_paper.py`` (like
``../carola/fig_carola_ocv_months.py`` and ``../morandi/fig_morandi_ocv_months.py``
are companions of their paper scripts): it reads only the committed channel table
``data/cts/cts_ocv_channels.csv``, recomputes per-window aggregates with
``cts_ocv_stats``, draws one figure and writes a small JSON + English report. It
adds **no** number to the A–F pin contract and touches none of it.

Figure (2 x 3, 183 mm wide):
  a  gamma2 per window   (per-chip dots, window median + IQR, the collapse cut)
  b  P       per window  (the mask-persistence channel)
  c  S       per window  (||centroid - peak||)
  d  D       per window  (mask density A / bbox)
  e  the (P, D) plane    (one dot per chip, coloured by window)
  f  the 6D fingerprint  (each window's six-dim median relative to F0)

The figure is descriptive and carries the site's standing caveats: both sides live
on the *same* structure, so season, weather, pipeline generation and the bare-sand
post state move with the collapse date as well, and the long ``F0 -> M+1`` baseline
conflates seven years of drift with the event. The report states that rather than a
mechanism.

Usage:
  python3 code/cts/fig_cts_ocv_months.py
  python3 code/cts/fig_cts_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import cts_ocv_core as core   # noqa: E402
import cts_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_cts_ocv_months.json")
OUT_MD = os.path.join(FIGDIR, "fig_cts_ocv_months.md")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

# The five windows (end exclusive), in time order. F0 opens the record where the
# first Sentinel-1 images of this footprint exist (2015-09); the collapse date
# itself opens M+1, so the pre side is the three M-* windows (plus the F0 baseline)
# and the post side is M+1 alone.
WINDOWS = [
    ("F0", "2015-09-21", "2015-11-01"),
    ("M-3", "2021-03-24", "2021-04-24"),
    ("M-2", "2021-04-24", "2021-05-24"),
    ("M-1", "2021-05-24", "2021-06-24"),
    ("M+1", "2021-06-24", "2021-07-24"),
]
WINDOW_LABELS = [w[0] for w in WINDOWS]
WINDOW_STATE = {"F0": core.PRE, "M-3": core.PRE, "M-2": core.PRE,
                "M-1": core.PRE, "M+1": core.POST}
WINDOW_RANGE = {w[0]: (w[1], w[2]) for w in WINDOWS}
# The record-start "first images" epoch gets its own colour; the three pre months
# are shades of the pre colour and M+1 the post colour.
WINDOW_COLORS = {"F0": "#55A868", "M-3": "#2E4E7E", "M-2": "#4C72B0",
                 "M-1": "#9DB8DC", "M+1": "#C44E52"}
# The panels (a)-(d): the four echo-mask dimensions in the order of the question.


def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the binaries stay out of git)."""
    out_dir = out_dir or FIGDIR
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png)}")
    return png


def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left")


def jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def window_of(day):
    """The window label of an acquisition day, or ``None`` if outside them."""
    for label, start, end in WINDOWS:
        if start <= day < end:
            return label
    return None


# ---------------------------------------------------------------------------
# Rows and the per-window blocks
# ---------------------------------------------------------------------------
def month_rows(csv_rows):
    """The chips of the five windows as this figure analyses them."""
    out = []
    for r in csv_rows:
        d = core.row_brief(r)
        d["window"] = window_of(d["day"])
        if d["window"] is None:
            continue
        d["state"] = core.state_of(d["day"])
        out.append(d)
    return out


def window_block(rows):
    """Per window: the count, the days, and ``st._stats`` of the six dimensions."""
    out = {}
    for label, start, end in WINDOWS:
        sel = [r for r in rows if r["window"] == label]
        out[label] = {
            "label": label,
            "window": [start, end],
            "state": WINDOW_STATE[label],
            "n": len(sel),
            "days": sorted({r["day"] for r in sel}),
            "date_range": ([min(r["day"] for r in sel), max(r["day"] for r in sel)]
                           if sel else None),
            "channels": {k: st._stats([r.get(k) for r in sel])
                         for k in FINGERPRINT_DIMS},
        }
    return out


def delta_blocks(rows, blocks):
    """M+1 against the baselines — per dimension, three comparisons.

    ``vs_F0`` compares the first post month with the record-start **first-images**
    epoch (the long ~7-year baseline the request asks for); ``vs_M-1`` with the
    last pre month (the two windows the collapse date separates); ``vs_pre_pool``
    with the three pre windows pooled. Positive Cliff's delta => M+1 larger.
    """
    f0 = [r for r in rows if r["window"] == "F0"]
    pre_pool = [r for r in rows if r["window"] in PRE_WINDOWS]
    m1 = [r for r in rows if r["window"] == "M-1"]
    post = [r for r in rows if r["window"] == "M+1"]
    out = {}
    for k in FINGERPRINT_DIMS:
        out[k] = {
            "vs_F0": st.delta_block([r[k] for r in f0], [r[k] for r in post],
                                    "M+1 vs F0 (first images)"),
            "vs_M-1": st.delta_block([r[k] for r in m1], [r[k] for r in post],
                                     "M+1 vs M-1"),
            "vs_pre_pool": st.delta_block([r[k] for r in pre_pool],
                                          [r[k] for r in post],
                                          "M+1 vs pooled pre"),
        }
    return out


def base_median(blocks):
    """Per-dimension mean of the baseline-window medians (the fingerprint base)."""
    base = {}
    for k in FINGERPRINT_DIMS:
        meds = [blocks[w]["channels"][k]["median"] for w in BASELINE_WINDOWS
                if blocks[w]["channels"][k] is not None]
        base[k] = float(np.mean(meds)) if meds else None
    return base


def fingerprint(blocks, base):
    """``{window: [median / baseline median per dim]}`` — the 6D panel data."""
    out = {}
    for label, _s, _e in WINDOWS:
        row = []
        for k in FINGERPRINT_DIMS:
            med = blocks[label]["channels"][k]
            b = base.get(k)
            row.append(med["median"] / b if (med is not None and b not in (None, 0.0))
                       else None)
        out[label] = row
    return out


def compute():
    csv_rows = core.load_csv()
    rows = month_rows(csv_rows)
    t0 = time.time()
    blocks = window_block(rows)
    res = {
        "rows": rows,
        "blocks": blocks,
        "deltas": delta_blocks(rows, blocks),
        "baseline": base_median(blocks),
    }
    res["fingerprint"] = fingerprint(blocks, res["baseline"])
    res["n_rows"] = len(rows)
    res["partition"] = {"rows": len(rows),
                        "per_window": {w[0]: blocks[w[0]]["n"] for w in WINDOWS}}
    res["runtime_s"] = time.time() - t0
    return res

STRIP_DIMS = ["gamma2", "P", "S", "D"]
# The 6D fingerprint runs over the package's own feature order.
FINGERPRINT_DIMS = list(core.FEATURES)
# The fingerprint baseline is the first-images epoch (the record start), which is
# the reference the site's request asks the month-after window to be read against.
BASELINE_WINDOWS = ["F0"]
# The three pre-event months pooled (the local, adjacent baseline).
PRE_WINDOWS = ["M-3", "M-2", "M-1"]
# Where the collapse cut and the long record gap are drawn on the strip axis.
COLLAPSE_CUT_X = WINDOW_LABELS.index("M+1") - 0.5        # between M-1 and M+1
GAP_X = WINDOW_LABELS.index("M-3") - 0.5                 # between F0 and M-3


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def window_strip(ax, rows, key, seed0=300, span=0.20, fs=6.2):
    """One jittered strip per window with the window median and IQR."""
    for i, label in enumerate(WINDOW_LABELS):
        vals = [r[key] for r in rows
                if r["window"] == label and r.get(key) is not None]
        if not vals:
            continue
        x = i + jitter(len(vals), span, seed0 + i)
        ax.plot(x, vals, "o", ms=1.7, alpha=0.32, color=WINDOW_COLORS[label],
                markeredgewidth=0)
        q25, q50, q75 = np.percentile(vals, [25, 50, 75])
        ax.plot([i - 0.30, i + 0.30], [q50, q50], "-", lw=1.9,
                color=WINDOW_COLORS[label])
        ax.plot([i, i], [q25, q75], "-", lw=1.0, color=WINDOW_COLORS[label],
                alpha=0.85)
        ax.text(i, 0.99, f"n={len(vals)}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=fs, color=WINDOW_COLORS[label])
    ax.margins(y=0.10)
    ax.set_xticks(range(len(WINDOW_LABELS)))
    ax.set_xticklabels(WINDOW_LABELS)
    ax.axvline(GAP_X, color="#BBBBBB", ls=":", lw=0.9)
    ax.axvline(COLLAPSE_CUT_X, color="#888888", ls="--", lw=0.9)


def fig_months(res, out_dir=None):
    rows = res["rows"]
    blocks = res["blocks"]
    fp, deltas = res["fingerprint"], res["deltas"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 155 * MM))

    # (a)-(d) the four echo-mask dimensions, one strip panel each
    letters = ["a", "b", "c", "d"]
    for letter, k in zip(letters, STRIP_DIMS):
        i = STRIP_DIMS.index(k)
        ax = axes[i // 3][i % 3]
        window_strip(ax, rows, k, seed0=300 + 7 * i)
        b = deltas[k]["vs_F0"]
        title = core.CHANNEL_AXIS.get(k, k)
        if "delta" in b:
            title += f"\nΔ(M+1−F0) = {b['delta']:+.2f}"
        ax.set_title(f"({letter}) {title}", loc="left", fontsize=8.0)
    # one collapse label, along the dashed cut of the first strip panel
    axes[0][0].text(COLLAPSE_CUT_X + 0.02, 0.5, "collapse 2021-06-24",
                    rotation=90, transform=axes[0][0].get_xaxis_transform(),
                    ha="left", va="center", fontsize=6.0, color="#555555")
    axes[0][0].text(GAP_X + 0.02, 0.5, "7-yr gap", rotation=90,
                    transform=axes[0][0].get_xaxis_transform(),
                    ha="left", va="center", fontsize=6.0, color="#999999")

    # (e) the (P, D) plane, one dot per chip coloured by window
    ax = axes[1][1]
    for label in WINDOW_LABELS:
        sel = [r for r in rows if r["window"] == label
               and r.get("P") is not None and r.get("D") is not None]
        ax.plot([r["P"] for r in sel], [r["D"] for r in sel], "o", ms=1.8,
                alpha=0.35, color=WINDOW_COLORS[label], markeredgewidth=0,
                label=label)
        pb = blocks[label]["channels"]["P"]
        db = blocks[label]["channels"]["D"]
        if pb and db:
            ax.plot(pb["median"], db["median"], "o", ms=6.0, mfc="none",
                    mec=WINDOW_COLORS[label], mew=1.2)
    ax.set_xlabel("P  (mask persistence)")
    ax.set_ylabel("D  (mask density A / bbox)")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.legend(frameon=False, ncol=2, markerscale=3)
    panel(ax, "e", "the (P, D) plane per window")

    # (f) the 6D fingerprint of all six dimensions relative to the F0 baseline
    ax = axes[1][2]
    xs = np.arange(len(FINGERPRINT_DIMS))
    for label in WINDOW_LABELS:
        ys = [v if v is not None else np.nan for v in fp[label]]
        ax.plot(xs, ys, "-o", ms=2.8, lw=1.1, color=WINDOW_COLORS[label],
                label=label)
    ax.axhline(1.0, color="#888888", lw=0.8)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(FINGERPRINT_DIMS)
    ax.set_ylabel("median / F0 median")
    ax.legend(frameon=False, ncol=2)
    panel(ax, "f", "the 6D fingerprint (vs. F0)")

    fig.suptitle(f"(CTS) echo-mask dimensions from the first images to the month "
                 f"after the collapse — {res['n_rows']} chips, F0 + M-3..M+1 of "
                 f"2021-06-24", fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    return save_fig(fig, "fig_cts_ocv_months", out_dir)



# ---------------------------------------------------------------------------
# JSON artifact and English report
# ---------------------------------------------------------------------------
def _fmt(x, nd=4):
    if x is None:
        return "–"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return out


def _delta_columns(b):
    """One ``[median0, median1, delta, sds, ci≠0, ci]`` row for a delta block."""
    if "delta" not in b:
        return ["–", "–", "–", "–", "–", "–"]
    ci = b.get("delta_ci95")
    return [_fmt(b["median"][0]), _fmt(b["median"][1]), _fmt(b["delta"]),
            _fmt(b["sds"]), "yes" if b.get("ci_excludes_zero") else "no",
            f"[{_fmt(ci[0])}, {_fmt(ci[1])}]" if ci else "–"]


def report(res):
    blocks, deltas, base = res["blocks"], res["deltas"], res["baseline"]
    md = []
    md.append("# Champlain Towers South, Surfside (FL) — the echo-mask dimensions "
              "from the first images to the month after the collapse")
    md.append("")
    md.append(f"**Site** {core.SITE} · **event** {core.COLLAPSE_DATE} · "
              f"**chips** {res['n_rows']} (windows F0 + M-3..M+1 only)")
    md.append("")
    md.append("```")
    md.append("x     = [gamma2, P, D, A, F, S]      echo-mask vector")
    md.append("mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2")
    md.append("F0    [2015-09-21, 2015-11-01)       pre-collapse (first images available)")
    md.append("M-3   [2021-03-24, 2021-04-24)       pre-collapse (healthy)")
    md.append("M-2   [2021-04-24, 2021-05-24)       pre-collapse (healthy)")
    md.append("M-1   [2021-05-24, 2021-06-24)       pre-collapse (healthy)")
    md.append("M+1   [2021-06-24, 2021-07-24)       post-collapse (month after)")
    md.append("```")
    md.append("")
    md.append("## Why the baseline is 2015 and not 2014")
    md.append("")
    md.append("Sentinel-1 data for **this** footprint do not exist in 2014. "
              "Sentinel-1A was launched on 2014-04-03, and the suitable orbit "
              "(ASC rel-48 IW3) has its **first acquisition on 2015-09-21** "
              "(`data/cts/cts_acquisition_census.json`; the Step-0 census already "
              "queried from 2014-10-01 and found nothing usable earlier — the "
              "descending stratum carries only 6 acquisitions in six years). `F0` "
              "is therefore the record-start epoch that carries the first images "
              "that actually exist, and is read here as the *first-images* "
              "baseline. No 2014 window is invented.")
    md.append("")
    md.append("## The five windows")
    md.append("")
    md.append("`F0` opens the record where the first Sentinel-1 images of this "
              "footprint exist; the collapse date opens `M+1`, so the pre side is "
              "the three `M-*` windows (plus the `F0` baseline) and the post side "
              "is `M+1` alone. This is a companion of `fig_cts_ocv_paper.py` "
              "(figures A–F): it adds one figure and does not touch the A–F pin "
              "contract.")
    md.append("")
    body = []
    for label, s, e in WINDOWS:
        b = blocks[label]
        body.append([label, f"[{s}, {e})", core.STATE_SHORT[b["state"]], b["n"],
                     b["date_range"][0] if b["date_range"] else "–",
                     b["date_range"][1] if b["date_range"] else "–",
                     len(b["days"])])
    md += md_table(["window", "range", "state", "chips", "first day", "last day",
                    "days"], body)
    md.append("")

    md.append("## Per-window medians of the six dimensions")
    md.append("")
    body = []
    for k in FINGERPRINT_DIMS:
        row = [k, _fmt(base.get(k))]
        for label in WINDOW_LABELS:
            c = blocks[label]["channels"][k]
            row.append(_fmt(c["median"]) if c else "–")
        body.append(row)
    md += md_table(["dim", "F0 (first images)", "F0"] + WINDOW_LABELS[1:], body)
    md.append("")
    md.append("`F0 (first images)` = the record-start baseline (the 1.0 line of "
              "panel f).")
    md.append("")
    md.append("## The long baseline — M+1 vs F0 (the requested contrast)")
    md.append("")
    body = []
    for k in FINGERPRINT_DIMS:
        body.append([k] + _delta_columns(deltas[k]["vs_F0"]))
    md += md_table(["dim", "median F0", "median M+1", "Cliff's delta", "SDS",
                    "CI95 excludes 0", "CI95"], body)
    md.append("")
    md.append("## The adjacent window — M+1 vs M-1")
    md.append("")
    body = []
    for k in FINGERPRINT_DIMS:
        body.append([k] + _delta_columns(deltas[k]["vs_M-1"]))
    md += md_table(["dim", "median M-1", "median M+1", "Cliff's delta", "SDS",
                    "CI95 excludes 0", "CI95"], body)
    md.append("")
    md.append("Cliff's delta > 0 means M+1 larger. The same comparison against the "
              "three pre windows pooled is in the JSON "
              "(`deltas[dim].vs_pre_pool`).")
    md.append("")
    md.append("## Reading (descriptive — not a claim of damage)")
    md.append("")
    n_f0 = sum(1 for k in FINGERPRINT_DIMS
               if deltas[k]["vs_F0"].get("ci_excludes_zero"))
    n_m1 = sum(1 for k in FINGERPRINT_DIMS
               if deltas[k]["vs_M-1"].get("ci_excludes_zero"))
    md.append(f"Against the **first-images** epoch `F0` the long baseline clears "
              f"the bootstrap CI in **{n_f0}** of the six dimensions; against the "
              f"adjacent `M-1` it is **{n_m1}** of six. The `F0` baseline carries "
              f"only {blocks['F0']['n']} chips, seven years and a different "
              "scattering state before the event, so its intervals are wide and "
              "the `F0 -> M+1` step conflates long-term drift, season and the "
              "collapse. The month axis therefore does **not** by itself "
              "reproduce a damage signal; it shows the full spread of the record "
              "and where the collapse sits inside it.")
    md.append("")
    md.append("Both sides live on the *same* structure, and the post state is "
              "**bare sand** whose intensity is moisture-driven, so season, "
              "weather, pipeline generation and the cleared lot move with the "
              "collapse date as well; the panels are read as that co-variation, "
              "**never** as damage.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/cts/cts_ocv_channels_csv.py    # build + verify the table")
    md.append("python3 code/cts/fig_cts_ocv_months.py      # this figure")
    md.append("```")
    return "\n".join(md) + "\n"


def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def build_json(res):
    """The result JSON — the per-window blocks and the deltas, not the per-row rows."""
    return {
        "site": core.SITE,
        "event": core.COLLAPSE_DATE,
        "generated_by": "code/cts/fig_cts_ocv_months.py",
        "features": list(core.FEATURES),
        "strip_dims": list(STRIP_DIMS),
        "baseline_windows": list(BASELINE_WINDOWS),
        "windows": [[w[0], w[1], w[2]] for w in WINDOWS],
        "n_rows": res["n_rows"],
        "partition": res["partition"],
        "baseline": res["baseline"],
        "fingerprint": res["fingerprint"],
        "blocks": res["blocks"],
        "deltas": res["deltas"],
        "figures": ["figures/cts/fig_cts_ocv_months.png"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, windows F0 + M-3..M+1 of {core.COLLAPSE_DATE} …")
    res = compute()
    print(f"  {res['n_rows']} chips "
          f"({res['partition']['per_window']}), {res['runtime_s']:.1f} s")

    os.makedirs(args.figdir, exist_ok=True)
    fig_months(res, args.figdir)
    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(build_json(res), fh, indent=1, sort_keys=True, default=json_default)
    print(f"  {os.path.relpath(args.json)}")
    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res))
    print(f"  {os.path.relpath(args.md)}")
    print(f"DONE ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

