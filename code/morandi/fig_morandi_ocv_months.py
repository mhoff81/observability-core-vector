#!/usr/bin/env python3
"""Morandi OCV-Paper — the echo-mask dimensions month by month around the collapse.

Question: what do the six dimensions of the **echo-mask vector**
``x = [gamma2, P, D, A, F, S]`` do in the months *around* the collapse
(``2018-08-14``) — the **three months before** it and the **first month after**
it? The A-E figures of ``fig_morandi_ocv_paper.py`` compare the two pooled states
(``pre-collapse (healthy)`` vs. ``post-collapse``); this companion figure keeps
the time axis instead and shows the four one-month windows that straddle the
event, so the reader sees *when* the mask moved, not only *that* it moved.

The four windows are rolling one-month buckets **anchored at the event** (the
collapse date belongs to the first post window, so M-1 is pure pre and M+1 pure
post):

    M-3   [2018-05-14, 2018-06-14)     pre-collapse (healthy)
    M-2   [2018-06-14, 2018-07-14)     pre-collapse (healthy)
    M-1   [2018-07-14, 2018-08-14)     pre-collapse (healthy)
    M+1   [2018-08-14, 2018-09-14)     post-collapse

Approach: this is a **standalone companion** of ``fig_morandi_ocv_paper.py``: it
reads only the committed channel table ``data/morandi/morandi_ocv_channels.csv``,
recomputes per-window aggregates with ``morandi_ocv_stats``, draws one figure and
writes a small JSON + English report. It adds **no** number to the A-E pin
contract and touches none of it.

Figure (2 x 3, 183 mm wide):
  a  gamma2 per window   (per-chip dots, window median + IQR, the collapse cut)
  b  P       per window  (the mask-persistence collapse)
  c  S       per window  (||centroid - peak||)
  d  D       per window  (mask density A / bbox)
  e  the (P, D) plane    (one dot per chip, coloured by window)
  f  the 6D fingerprint  (each window's six-dim median relative to the pre-baseline)

The figure is descriptive and carries the site's standing caveat: both sides live
on the *same* bridge, so season, weather, traffic and pipeline generation move
with the collapse date as well. The report states that rather than a mechanism.

Usage:
  python3 code/morandi/fig_morandi_ocv_months.py
  python3 code/morandi/fig_morandi_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
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

import morandi_ocv_core as core   # noqa: E402
import morandi_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_morandi_ocv_months.json")
OUT_MD = os.path.join(FIGDIR, "fig_morandi_ocv_months.md")

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

# The four one-month windows anchored at the collapse (end exclusive), in time
# order. The collapse date itself opens M+1, so the pre side is the three M-*
# windows and the post side is M+1 alone.
WINDOWS = [
    ("M-3", "2018-05-14", "2018-06-14"),
    ("M-2", "2018-06-14", "2018-07-14"),
    ("M-1", "2018-07-14", "2018-08-14"),
    ("M+1", "2018-08-14", "2018-09-14"),
]
WINDOW_LABELS = [w[0] for w in WINDOWS]
WINDOW_STATE = {"M-3": core.PRE, "M-2": core.PRE, "M-1": core.PRE, "M+1": core.POST}
# Three shades of the pre colour for the three pre windows, the post colour for
# the window the collapse opens.
WINDOW_COLORS = {"M-3": "#2E4E7E", "M-2": "#4C72B0", "M-1": "#9DB8DC",
                 "M+1": "#C44E52"}
# The panels (a)-(d): the four echo-mask dimensions in the order of the question.
STRIP_DIMS = ["gamma2", "P", "S", "D"]
# The 6D fingerprint runs over the package's own feature order.
FINGERPRINT_DIMS = list(core.FEATURES)
BASELINE_WINDOWS = ["M-3", "M-2", "M-1"]


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
    """The chips of the four windows as this figure analyses them."""
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
    """M+1 against the pre side — per dimension, both baselines.

    ``vs_M-1`` compares the first post month with the last pre month (the two
    windows the collapse date separates); ``vs_pre_pool`` compares it with the
    three pre windows pooled. Positive Cliff's delta => M+1 larger.
    """
    pre_pool = [r for r in rows if WINDOW_STATE[r["window"]] == core.PRE]
    m1 = [r for r in rows if r["window"] == "M-1"]
    post = [r for r in rows if r["window"] == "M+1"]
    out = {}
    for k in FINGERPRINT_DIMS:
        out[k] = {
            "vs_M-1": st.delta_block([r[k] for r in m1], [r[k] for r in post],
                                     "M+1 vs M-1"),
            "vs_pre_pool": st.delta_block([r[k] for r in pre_pool],
                                          [r[k] for r in post],
                                          "M+1 vs pooled pre"),
        }
    return out


def pre_baseline(blocks):
    """Per-dimension mean of the three pre-window medians (the fingerprint base)."""
    base = {}
    for k in FINGERPRINT_DIMS:
        meds = [blocks[w]["channels"][k]["median"] for w in BASELINE_WINDOWS
                if blocks[w]["channels"][k] is not None]
        base[k] = float(np.mean(meds)) if meds else None
    return base


def fingerprint(blocks, base):
    """``{window: [median / pre-baseline median per dim]}`` — the 6D panel data."""
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
        "baseline": pre_baseline(blocks),
    }
    res["fingerprint"] = fingerprint(blocks, res["baseline"])
    res["n_rows"] = len(rows)
    res["partition"] = {"rows": len(rows),
                        "per_window": {w[0]: blocks[w[0]]["n"] for w in WINDOWS}}
    res["runtime_s"] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def window_strip(ax, rows, key, seed0=300, span=0.20, fs=6.2):
    """Four jittered strips (one per window) with the window median and IQR."""
    for i, label in enumerate(WINDOW_LABELS):
        vals = [r[key] for r in rows if r["window"] == label and r.get(key) is not None]
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
    ax.axvline(2.5, color="#888888", ls="--", lw=0.9)


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
        b = deltas[k]["vs_M-1"]
        title = core.CHANNEL_AXIS.get(k, k)
        if "delta" in b:
            title += f"\nΔ(M+1−M−1) = {b['delta']:+.2f}"
        ax.set_title(f"({letter}) {title}", loc="left", fontsize=8.0)
    # one collapse label, along the dashed cut of the first strip panel
    axes[0][0].text(2.52, 0.5, "collapse 2018-08-14", rotation=90,
                    transform=axes[0][0].get_xaxis_transform(), ha="left",
                    va="center", fontsize=6.2, color="#555555")

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

    # (f) the 6D fingerprint of all six dimensions relative to the pre-baseline
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
    ax.set_ylabel("median / pre-baseline median")
    ax.legend(frameon=False, ncol=2)
    panel(ax, "f", "the 6D fingerprint")

    fig.suptitle(f"(Morandi) echo-mask dimensions month by month around the "
                 f"collapse — {res['n_rows']} chips, windows M-3..M+1 of "
                 f"2018-08-14", fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    return save_fig(fig, "fig_morandi_ocv_months", out_dir)


# ---------------------------------------------------------------------------
# English report
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


def report(res):
    blocks, deltas, part = res["blocks"], res["deltas"], res["partition"]
    md = []
    md.append("# Ponte Morandi (Polcevera), Genova — the echo-mask dimensions "
              "month by month around the collapse")
    md.append("")
    md.append(f"**Site** {core.SITE} · **event** {core.COLLAPSE_DATE} · "
              f"**chips** {res['n_rows']} in the four one-month windows "
              "M-3 … M+1")
    md.append("")
    md.append("```")
    for label, start, end in WINDOWS:
        md.append(f"{label:<4} [{start}, {end})   {WINDOW_STATE[label]}")
    md.append("```")
    md.append("")
    md.append("## Windows")
    md.append("")
    rows = []
    for label, start, end in WINDOWS:
        b = blocks[label]
        dr = b["date_range"]
        rows.append([label, f"[{start}, {end})", core.STATE_SHORT[b["state"]],
                     b["n"], f"{dr[0]} … {dr[1]}" if dr else "–"])
    md += md_table(["window", "interval", "state", "n", "days"], rows)
    md.append("")
    md.append("## Median of each dimension per window")
    md.append("")
    rows = []
    for k in FINGERPRINT_DIMS:
        row = [k]
        for label, _s, _e in WINDOWS:
            c = blocks[label]["channels"][k]
            row.append(_fmt(c["median"]) if c else "–")
        rows.append(row)
    md += md_table(["dimension"] + WINDOW_LABELS, rows)
    md.append("")
    md.append("## M+1 against the pre side")
    md.append("")
    rows = []
    for k in FINGERPRINT_DIMS:
        a = deltas[k]["vs_M-1"]
        b = deltas[k]["vs_pre_pool"]
        rows.append([k, _fmt(a.get("delta")),
                     "yes" if a.get("ci_excludes_zero") else "no",
                     _fmt(b.get("delta")),
                     "yes" if b.get("ci_excludes_zero") else "no"])
    md += md_table(["dimension", "δ vs M-1", "CI≠0 vs M-1",
                    "δ vs pooled pre", "CI≠0 vs pooled pre"], rows)
    md.append("")
    n_excl = sum(1 for k in FINGERPRINT_DIMS
                 if deltas[k]["vs_M-1"].get("ci_excludes_zero"))
    md.append("## What the figure shows")
    md.append("")
    md.append(f"The first post month carries **{part['per_window']['M+1']}** "
              f"chips against "
              f"{part['per_window']['M-3'] + part['per_window']['M-2'] + part['per_window']['M-1']} "
              "before the event, so its intervals are wide and this is a "
              "descriptive companion, not a test. Of the six dimensions, "
              f"**{n_excl}** show a Cliff's-delta 95 % CI that excludes zero "
              "against M-1; the (P, D) plane (panel e) and the 6D fingerprint "
              "(panel f) show *when* the mask moved across the four windows.")
    md.append("")
    md.append("Both sides of the event live on the *same* bridge, so season, "
              "weather and pipeline generation move with the collapse date too; "
              "the panels are read as that co-variation, **never** as damage.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/morandi/morandi_ocv_channels_csv.py    # build + verify the table")
    md.append("python3 code/morandi/fig_morandi_ocv_months.py     # this figure + report")
    md.append("```")
    return "\n".join(md) + "\n"


def build_json(res):
    return {
        "site": core.SITE,
        "event": core.COLLAPSE_DATE,
        "generated_by": "code/morandi/fig_morandi_ocv_months.py",
        "states": dict(core.STATE_SHORT),
        "features": list(core.FEATURES),
        "windows": [{"label": l, "start": s, "end": e, "state": WINDOW_STATE[l]}
                    for l, s, e in WINDOWS],
        "n_rows": res["n_rows"],
        "partition": res["partition"],
        "baseline": res["baseline"],
        "fingerprint": res["fingerprint"],
        "blocks": res["blocks"],
        "deltas": res["deltas"],
        "figure": "figures/morandi/fig_morandi_ocv_months.png",
    }


def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, event {core.COLLAPSE_DATE} …")
    res = compute()
    print(f"  {res['n_rows']} chips in M-3 … M+1 "
          f"{res['partition']['per_window']}")

    os.makedirs(args.figdir, exist_ok=True)
    fig_months(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(build_json(res), fh, indent=1, sort_keys=True,
                  default=json_default)
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res))
    print(f"  {os.path.relpath(args.md)}")
    print(f"DONE ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())