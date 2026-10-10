#!/usr/bin/env python3
"""CTS OCV-Paper — the echo-mask dimensions from the first images to the month
after the collapse, over the 24 pre-collapse months.

Question: what do the six dimensions of the **echo-mask vector**
``x = [gamma2, P, D, A, F, S]`` do along the record — from the *first Sentinel-1
images available for this footprint* (``F0``, the 2015 record start) through the
**24 one-month windows before** the collapse to the **first month after** it
(``2021-06-24``)? The A–F figures of ``fig_cts_ocv_paper.py`` pool the two states
(``pre-collapse (healthy)`` vs. ``post-collapse``); this companion keeps the time
axis and shows both the long-baseline contrast the site's request asks for, the
whole two years of monthly pre-collapse movement, and the terminal month.

Why ``F0`` is 2015, not 2014
----------------------------
Sentinel-1 data for **this** footprint do not exist in 2014: Sentinel-1A was
launched on 2014-04-03 and the suitable orbit (ASC rel-48 IW3) has its **first
acquisition on 2015-09-21** (``data/cts/cts_acquisition_census.json``; the
Step-0 census already queried from 2014-10-01 and found nothing usable earlier).
``F0`` is therefore the record-start epoch that carries the first images that
actually exist, and the report states that rather than inventing a 2014 window.

The windows (end exclusive), in time order:

    F0    [2015-09-21, 2015-11-01)     the first images available (record start)
    M-24  [2019-06-24, 2019-07-24)     pre-collapse (healthy)  first record month
    ...
    M-1   [2021-05-24, 2021-06-24)     pre-collapse (healthy)  last
    M+1   [2021-06-24, 2021-07-24)     post-collapse (the month after)

Approach: this is a **standalone companion** of ``fig_cts_ocv_paper.py`` (like
``../carola/fig_carola_ocv_months.py`` and ``../morandi/fig_morandi_ocv_months.py``
are companions of their paper scripts): it reads only the committed channel table
``data/cts/cts_ocv_channels.csv``, recomputes per-month aggregates with
``cts_ocv_stats``, draws one figure and writes a small JSON + English report. It
adds **no** number to the A–F pin contract and touches none of it.

Figure (2 x 3, 183 mm wide):
  a  gamma2 by month   (per-chip dots, monthly median + IQR, the collapse cut)
  b  P       by month  (the mask-persistence channel)
  c  D       by month  (mask density A / bbox)
  d  A       by month  (mask size, masked px)
  e  the (P, D) plane  (one dot per chip, coloured by month)
  f  the 6D fingerprint (each month's six-dim median relative to the F0 baseline,
     drawn as a line per dimension over the F0 + 24-month + M+1 axis)

The figure is descriptive and carries the site's standing caveats: both sides live
on the *same* structure, so season, weather, pipeline generation and the bare-sand
post state move with the collapse date as well, and the long ``F0 -> M+1`` baseline
spans six years of drift, so it conflates that drift with the event. The report
states that rather than a mechanism.

Usage:
  python3 code/cts/fig_cts_ocv_months.py
  python3 code/cts/fig_cts_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
"""
from __future__ import annotations

import argparse
import datetime as dt
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

# F0 opens the record where the first Sentinel-1 images of this footprint exist
# (2015-09). The 24 rolling one-month buckets anchor at the collapse date, and the
# collapse date itself opens M+1, so the pre side is the 24 M-* windows (plus the
# F0 baseline) and the post side is M+1 alone.
F0_WINDOW = ("F0", "2015-09-21", "2015-11-01")
N_MONTHS = 24          # the pre-collapse series: M-24 … M-1


def shift_month(day, k):
    """``YYYY-MM-DD`` shifted by ``k`` whole calendar months (anchor day kept).

    The event day of this site (24) and of every site in this package is <= 28, so
    the anchor day survives every month; the helper is pure integer arithmetic
    (no dateutil).
    """
    y, m, d = (int(x) for x in day.split("-"))
    ny, nm = divmod(y * 12 + (m - 1) + k, 12)
    return f"{ny:04d}-{nm + 1:02d}-{d:02d}"


def make_windows(event, n_pre=N_MONTHS):
    """The ``n_pre`` one-month buckets before the event, then the month it opens.

    ``M-k`` spans ``[event - k months, event - (k-1) months)``, so the collapse
    date itself belongs to ``M+1`` and ``M-1`` is pure pre-collapse.
    """
    out = [(f"M-{k}", shift_month(event, -k), shift_month(event, -(k - 1)))
           for k in range(n_pre, 0, -1)]
    out.append(("M+1", event, shift_month(event, 1)))
    return out


# F0 (first images) + the 24 pre-collapse months + the month the collapse opens.
WINDOWS = [F0_WINDOW] + make_windows(core.COLLAPSE_DATE)
WINDOW_LABELS = [w[0] for w in WINDOWS]
WINDOW_STATE = {w[0]: (core.POST if w[0] == "M+1" else core.PRE)
                for w in WINDOWS}
WINDOW_RANGE = {w[0]: (w[1], w[2]) for w in WINDOWS}


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
    """The chips that fall in any of the windows (F0, then M-24 … M+1)."""
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
    with all ``N_MONTHS`` pre months pooled. Positive Cliff's delta => M+1 larger.
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

STRIP_DIMS = ["gamma2", "P", "D", "A"]
# A couple of dimensions span decades across the record; those panels get a symlog
# y axis (linthresh) so the monthly median line stays readable.
STRIP_SCALES = {"D": 0.01}
# The 6D fingerprint runs over the package's own feature order.
FINGERPRINT_DIMS = list(core.FEATURES)
# The fingerprint baseline is the first-images epoch (the record start), which is
# the reference the site's request asks the month-after window to be read against.
BASELINE_WINDOWS = ["F0"]
# The 24 pre-event months (M-24 .. M-1) pooled: the local, adjacent baseline.
PRE_WINDOWS = [w[0] for w in WINDOWS if w[0] not in ("F0", "M+1")]
# Where the collapse cut and the long record gap are drawn on the month axis.
COLLAPSE_CUT_X = WINDOW_LABELS.index("M+1") - 0.5     # between M-1 and M+1
GAP_X = WINDOW_LABELS.index("M-24") - 0.5             # between F0 and the months


def span_years(first_day, last_day):
    """Years between two ``YYYY-MM-DD`` days (365.25 days per year)."""
    a = dt.date(*map(int, first_day.split("-")))
    b = dt.date(*map(int, last_day.split("-")))
    return (b - a).days / 365.25


# The two spans the panels and the report name, both measured from the record:
# the F0 -> M-24 acquisition gap (drawn as the dotted line) and the whole
# F0 -> M+1 baseline the month-after window is read against.
GAP_YEARS = span_years(F0_WINDOW[1], WINDOWS[1][1])
GAP_LABEL = f"gap {GAP_YEARS:.1f} yr"
BASE_SPAN_YEARS = span_years(F0_WINDOW[1], WINDOWS[-1][2])


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def time_axis(ax):
    """Month-axis ticks (every third window), the collapse cut and the record gap."""
    ticks = [i for i in range(len(WINDOW_LABELS)) if i % 3 == 0]
    if (len(WINDOW_LABELS) - 1) not in ticks:
        ticks.append(len(WINDOW_LABELS) - 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([WINDOWS[i][1][:7] for i in ticks], rotation=90,
                       fontsize=6.0)
    ax.set_xlim(-1.0, len(WINDOW_LABELS) + 1.2)
    ax.axvline(GAP_X, color="#BBBBBB", ls=":", lw=0.9)
    ax.axvline(COLLAPSE_CUT_X, color="#888888", ls="--", lw=0.9)


def month_line(ax, rows, key, fs=5.2, linthresh=None):
    """One month-by-month panel: per-chip dots, the monthly median and its IQR.

    The per-month chip count is written only on every sixth window (and the last),
    so 26 numbers do not collide above a narrow panel. The 2015 -> 2019 record gap
    breaks the median line instead of joining two epochs nine years apart.
    """
    xs = [float(i) for i in range(len(WINDOW_LABELS))]
    med, q25, q75 = [], [], []
    for i, label in enumerate(WINDOW_LABELS):
        vals = [r[key] for r in rows
                if r["window"] == label and r.get(key) is not None]
        if vals:
            ax.plot(i + jitter(len(vals), 0.16, 300 + 7 * i), vals, "o", ms=1.4,
                    alpha=0.28, color="#8A8A8A", markeredgewidth=0)
            lo, mid, hi = np.percentile(vals, [25, 50, 75])
            if i % 6 == 0 or i == len(WINDOW_LABELS) - 1:
                ax.text(i, 0.995, f"{len(vals)}",
                        transform=ax.get_xaxis_transform(), ha="center",
                        va="top", fontsize=fs, color="#777777")
        else:
            lo = mid = hi = np.nan
        med.append(mid)
        q25.append(lo)
        q75.append(hi)
    # break the line across the 2015 -> 2019 gap (between F0 and M-24)
    xs = [xs[0], GAP_X] + xs[1:]
    med = [med[0], np.nan] + med[1:]
    q25 = [q25[0], np.nan] + q25[1:]
    q75 = [q75[0], np.nan] + q75[1:]
    xs = np.asarray(xs, float)
    ax.fill_between(xs, np.asarray(q25, float), np.asarray(q75, float),
                    color="#4C72B0", alpha=0.16, lw=0)
    ax.plot(xs, np.asarray(med, float), "-o", ms=2.4, lw=1.1, color="#2E4E7E")
    if linthresh is not None:
        ax.set_yscale("symlog", linthresh=linthresh)
    time_axis(ax)
    ax.margins(y=0.14)


def fig_months(res, out_dir=None):
    rows = res["rows"]
    blocks = res["blocks"]
    fp, deltas = res["fingerprint"], res["deltas"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 160 * MM))

    # (a)-(d) the four echo-mask dimensions, one month-by-month panel each
    for letter, k in zip(["a", "b", "c", "d"], STRIP_DIMS):
        i = STRIP_DIMS.index(k)
        ax = axes[i // 3][i % 3]
        month_line(ax, rows, k, linthresh=STRIP_SCALES.get(k))
        b = deltas[k]["vs_F0"]
        title = core.CHANNEL_AXIS.get(k, k)
        if "delta" in b:
            title += f"\nΔ(M+1−F0) = {b['delta']:+.2f}"
        ax.set_title(f"({letter}) {title}", loc="left", fontsize=8.0)
    # the two event markers, along the first panel
    axes[0][0].text(len(WINDOW_LABELS) + 0.35, 0.5,
                    f"collapse {core.COLLAPSE_DATE}", rotation=90,
                    transform=axes[0][0].get_xaxis_transform(), ha="center",
                    va="center", fontsize=6.0, color="#555555")
    axes[0][0].text(GAP_X + 0.12, 0.5, GAP_LABEL, rotation=90,
                    transform=axes[0][0].get_xaxis_transform(), ha="center",
                    va="center", fontsize=5.8, color="#999999")

    # (e) the (P, D) trajectory, one dot per chip coloured by month
    ax = axes[1][1]
    cmap = plt.cm.viridis
    for i, label in enumerate(WINDOW_LABELS):
        sel = [r for r in rows if r["window"] == label
               and r.get("P") is not None and r.get("D") is not None]
        if not sel:
            continue
        ax.plot([r["P"] for r in sel], [r["D"] for r in sel], "o", ms=1.8,
                alpha=0.40, color=cmap(i / max(len(WINDOW_LABELS) - 1, 1)),
                markeredgewidth=0)
    px = [blocks[l]["channels"]["P"]["median"] if blocks[l]["channels"]["P"]
          else np.nan for l in WINDOW_LABELS]
    dy = [blocks[l]["channels"]["D"]["median"] if blocks[l]["channels"]["D"]
          else np.nan for l in WINDOW_LABELS]
    ax.plot(px[1:-1], dy[1:-1], "-", lw=1.0, color="#333333", alpha=0.7)
    ax.plot(px, dy, "o", ms=3.2, mfc="none", mec="#333333", mew=0.8, ls="none")
    ax.plot([px[-1]], [dy[-1]], "o", ms=5.5, color="#C44E52")
    ax.set_xlabel("P  (mask persistence)")
    ax.set_ylabel("D  (mask density A / bbox)")
    ax.set_yscale("symlog", linthresh=0.01)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(0, len(WINDOW_LABELS) - 1))
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.02)
    cb.set_ticks([0, len(WINDOW_LABELS) - 1])
    cb.set_ticklabels([WINDOWS[0][1][:7], WINDOWS[-1][1][:7]])
    cb.ax.tick_params(labelsize=5.6)
    cb.ax.set_title("month", fontsize=6.0, pad=3)
    panel(ax, "e", "the (P, D) trajectory per month")

    # (f) the 6D fingerprint: a line per dimension over the F0 + 24-month axis
    ax = axes[1][2]
    xs = np.arange(len(WINDOW_LABELS))
    for j, k in enumerate(FINGERPRINT_DIMS):
        ys = [(fp[l][j] if (fp[l][j] is not None and fp[l][j] > 0) else np.nan)
              for l in WINDOW_LABELS]
        ax.plot(xs, ys, "-", lw=1.1, label=k, color=plt.cm.tab10(j % 10))
    ax.axhline(1.0, color="#888888", lw=0.8)
    time_axis(ax)
    ax.set_yscale("log")
    ax.set_ylabel("median / F0 median")
    ax.legend(frameon=True, framealpha=0.85, edgecolor="none", ncol=3,
              fontsize=5.8, loc="upper right")
    panel(ax, "f", f"the 6D fingerprint (baseline = {BASELINE_WINDOWS[0]})")

    fig.suptitle(f"(CTS) echo-mask dimensions from the first images (F0) through "
                 f"the {N_MONTHS} months before the collapse — {res['n_rows']} "
                 f"chips, M-{N_MONTHS}..M-1 of {core.COLLAPSE_DATE} plus M+1",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.960))
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
    fp = res["fingerprint"]
    md = []
    md.append("# Champlain Towers South, Surfside (FL) — the echo-mask dimensions "
              f"from the first images through the {N_MONTHS} months before the "
              "collapse to the month after it")
    md.append("")
    md.append(f"**Site** {core.SITE} · **event** {core.COLLAPSE_DATE} · "
              f"**chips** {res['n_rows']} (windows F0 + M-{N_MONTHS}..M+1)")
    md.append("")
    md.append("```")
    md.append("x     = [gamma2, P, D, A, F, S]      echo-mask vector")
    md.append("mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2")
    md.append(f"F0    [{F0_WINDOW[1]}, {F0_WINDOW[2]})       "
              "pre-collapse (first images available)")
    md.append(f"M-{N_MONTHS:<3} [{WINDOWS[1][1]}, {WINDOWS[1][2]})       "
              "pre-collapse (healthy)  first record month")
    md.append(f"M-1   [{WINDOWS[-2][1]}, {WINDOWS[-2][2]})       "
              "pre-collapse (healthy)  last")
    md.append(f"M+1   [{WINDOWS[-1][1]}, {WINDOWS[-1][2]})       "
              "post-collapse (month after)")
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
    md.append(f"## F0 and the {N_MONTHS} pre-collapse months")
    md.append("")
    md.append("`F0` opens the record where the first Sentinel-1 images of this "
              f"footprint exist; the {N_MONTHS} rolling one-month buckets anchor "
              "at the collapse date and cover the two years up to the event, and "
              "the collapse date itself opens `M+1`. The pre side is therefore "
              f"the {N_MONTHS} `M-*` windows (plus the `F0` baseline) and the "
              "post side is `M+1` alone. This is a companion of "
              "`fig_cts_ocv_paper.py` (figures A–F): it adds one figure and does "
              "not touch the A–F pin contract.")
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

    md.append("## Median of each dimension per month")
    md.append("")
    body = []
    for label, _s, _e in WINDOWS:
        b = blocks[label]
        row = [label, b["n"]]
        for k in FINGERPRINT_DIMS:
            c = b["channels"][k]
            row.append(_fmt(c["median"]) if c else "–")
        body.append(row)
    md += md_table(["month", "n"] + FINGERPRINT_DIMS, body)
    md.append("")
    md.append(f"`baseline` = the `{BASELINE_WINDOWS[0]}` (first-images) medians — "
              "the 1.0 line of panel f: "
              + ", ".join(f"`{k}` {_fmt(base.get(k))}" for k in FINGERPRINT_DIMS)
              + ".")
    md.append("")
    md.append("## The 6D fingerprint at the cut (month / baseline)")
    md.append("")
    body = []
    for label in ["F0", PRE_WINDOWS[0], PRE_WINDOWS[len(PRE_WINDOWS) // 2],
                  PRE_WINDOWS[-1], "M+1"]:
        body.append([label] + [_fmt(fp[label][j])
                               for j in range(len(FINGERPRINT_DIMS))])
    md += md_table(["month"] + FINGERPRINT_DIMS, body)
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
    md.append("Cliff's delta > 0 means M+1 larger. The same comparison against all "
              f"{N_MONTHS} pre months pooled is in the JSON "
              "(`deltas[dim].vs_pre_pool`).")
    md.append("")
    md.append("## Reading (descriptive — not a claim of damage)")
    md.append("")
    n_f0 = sum(1 for k in FINGERPRINT_DIMS
               if deltas[k]["vs_F0"].get("ci_excludes_zero"))
    n_m1 = sum(1 for k in FINGERPRINT_DIMS
               if deltas[k]["vs_M-1"].get("ci_excludes_zero"))
    pre_chips = sum(blocks[l]["n"] for l in PRE_WINDOWS)
    md.append(f"The {N_MONTHS} pre-collapse months carry **{pre_chips}** chips "
              f"(median {int(np.median([blocks[l]['n'] for l in PRE_WINDOWS]))} per "
              f"month) and the first post month **{blocks['M+1']['n']}**, so the "
              "monthly intervals are wide and this is a descriptive companion, "
              "not a test.")
    md.append("")
    last_pre = ", ".join(f"`{k}` {_fmt(fp[PRE_WINDOWS[-1]][j])}"
                         for j, k in enumerate(FINGERPRINT_DIMS))
    first_post = ", ".join(f"`{k}` {_fmt(fp['M+1'][j])}"
                           for j, k in enumerate(FINGERPRINT_DIMS))
    md.append(f"Relative to the **first-images** epoch `F0` (= 1.0) the **last** "
              f"pre month M-1 sits at {last_pre}; the first post month M+1 sits at "
              f"{first_post}. The pre-collapse record is **not** a flat line — "
              "season, weather, pipeline generation and the bare-sand post state "
              "move the six dimensions month to month well before the event, and "
              f"the `F0 -> M+1` step spans {BASE_SPAN_YEARS:.0f} years and "
              "conflates that drift with the collapse. Against the "
              "**first-images** baseline `F0` the long "
              f"contrast clears the bootstrap CI in **{n_f0}** of the six "
              f"dimensions; against the adjacent `M-1` it is **{n_m1}** of six. "
              "The month axis therefore does **not** by itself reproduce a damage "
              "signal; it shows the full spread of the record and where the "
              "collapse sits inside it.")
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
        "months": N_MONTHS,
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
    print(f"computing — {core.SITE}, windows F0 + M-{N_MONTHS}..M+1 of "
          f"{core.COLLAPSE_DATE} …")
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

