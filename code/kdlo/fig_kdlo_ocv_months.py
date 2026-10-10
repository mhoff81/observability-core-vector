#!/usr/bin/env python3
"""KDLO OCV-Paper — the echo-mask dimensions month by month over the record.

Question: what do the six dimensions of the **echo-mask vector**
``x = [gamma2, P, D, A, F, S]`` do across the whole KDLO record — from the first
acquisition (``2022-06-02``) through the collapse of the mast to the last one
(``2024-09-19``)? The A-E figures of ``fig_kdlo_ocv_paper.py`` pool the two
**epochs** (``pre`` = standing mast vs. ``rebuild`` = during the reconstruction);
this companion keeps the time axis and rolls **calendar months** over the record,
so the reader sees *when* the mask moved, and where the collapse and the
acquisition gap sit inside it.

Why this is not a 24-month pre-collapse series
----------------------------------------------
The KDLO record is 60 acquisitions in two disjoint epochs, not a dense
twelve-day-repeat series: ``pre`` holds 2022-06-02..2022-12-11 (16 acquisitions,
six months) and ``rebuild`` 2023-05-04..2024-09-19 (43 acquisitions, seventeen
months). There is no pre-collapse year to roll, so the axis here is the whole
record by calendar month (28 buckets, four of them empty), and the epochs are
tinted instead of split into ``M-k`` windows.

Two facts the figure must not paper over, and does not:

  * the **collapse month 2022-12 is mixed** — it holds both the last intact
    acquisitions (``pre``, to 2022-12-11) and the single ``2022-12-23``
    acquisition that belongs to neither epoch (the record's ``single`` row);
  * the record has a **four-month acquisition gap** (2023-01..2023-04), so the
    month before the collapse and the first rebuild month are not adjacent.

Approach: this is a **standalone companion** of ``fig_kdlo_ocv_paper.py``: it
reads only the committed channel table ``data/kdlo/kdlo_ocv_channels.csv``,
recomputes per-month aggregates with ``kdlo_ocv_stats``, draws one figure and
writes a small JSON + English report. It adds **no** number to the A-E pin
contract and touches none of it. The CSV carries no ``day`` column; the day comes
from ``acquisition_ts`` via ``core.day_of``.

Figure (2 x 3, 183 mm wide):
  a  gamma2 by month   (per-chip dots, monthly median + IQR, the collapse cut)
  b  P       by month  (the mask-persistence channel)
  c  D       by month  (mask density A / bbox)
  d  A       by month  (mask size, masked px)
  e  the (P, D) plane  (one dot per acquisition, coloured by month)
  f  the 6D fingerprint (each month's six-dim median relative to the 2022-06
     baseline, drawn as a line per dimension over the whole record)

The figure is descriptive and carries the site's standing caveat: the two epochs
are disjoint time windows, so season, weather, the reconstruction itself and the
change of scattering geometry all move with the collapse as well. The report
states that rather than a mechanism.

Usage:
  python3 code/kdlo/fig_kdlo_ocv_months.py
  python3 code/kdlo/fig_kdlo_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
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

import kdlo_ocv_core as core   # noqa: E402
import kdlo_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_kdlo_ocv_months.json")
OUT_MD = os.path.join(FIGDIR, "fig_kdlo_ocv_months.md")

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

# The two epoch labels of the record, in chronological order. A KDLO "state" is a
# time window, not a damage label.
PRE = st.STATE_ORDER[0]            # "pre"     — standing mast
REBUILD = st.STATE_ORDER[1]        # "rebuild" — during the reconstruction
# The collapse pair the project works with: the last intact acquisition and the
# first collapsed one (the record's single "singleton" row).
COLLAPSE_DAY = core.STATE_WINDOWS[PRE][1]          # 2022-12-11, last intact
FIRST_COLLAPSED_DAY = core.SINGLETON_DAY           # 2022-12-23, first collapsed
# The panels (a)-(d): the four echo-mask dimensions in the question's order.
STRIP_DIMS = ["gamma2", "P", "D", "A"]
# D spans decades over the record; its panel gets a symlog y axis (linthresh).
STRIP_SCALES = {"D": 0.01}
# The 6D fingerprint runs over the package's own feature order.
FINGERPRINT_DIMS = list(core.FEATURES)


def month_span(y, m):
    """``(first day, first day of the next month)`` of the calendar month (y, m)."""
    ny, nm = divmod(y * 12 + m, 12)
    return f"{y:04d}-{m:02d}-01", f"{ny:04d}-{nm + 1:02d}-01"


def next_month(y, m):
    return (y + 1, 1) if m == 12 else (y, m + 1)


def record_months(csv_rows):
    """Every calendar month from the first to the last acquisition (gaps kept)."""
    days = sorted(d for d in (core.day_of(r) for r in csv_rows) if d)
    y, m = int(days[0][:4]), int(days[0][5:7])
    last = days[-1][:7]
    out = []
    while f"{y:04d}-{m:02d}" <= last:
        out.append((f"{y:04d}-{m:02d}",) + month_span(y, m))
        y, m = next_month(y, m)
    return out


# The analysed axis: the record's calendar months, in time order. Built from the
# committed CSV at import so that a change of record cannot silently stale it.
WINDOWS = record_months(core.load_csv())
WINDOW_LABELS = [w[0] for w in WINDOWS]
COLLAPSE_MONTH = COLLAPSE_DAY[:7]                  # "2022-12"


def month_frac_x(day):
    """Fractional x of a day inside its calendar-month bin (the exact cut)."""
    lab = day[:7]
    if lab not in WINDOW_LABELS:
        return float(len(WINDOW_LABELS)) - 0.5
    y, m, d = (int(v) for v in day.split("-"))
    first = dt.date(y, m, 1)
    ndays = ((first + dt.timedelta(days=32)).replace(day=1) - first).days
    return WINDOW_LABELS.index(lab) + (d - 1) / ndays


# Where the collapse falls on the calendar-month axis (2022-12-11 -> x ~ 6.32):
# inside the month, not on a bucket boundary as it is for the anchored sites.
COLLAPSE_CUT_X = month_frac_x(COLLAPSE_DAY)
# The fingerprint baseline is the first record month (2022-06): the only month
# whose value does not depend on any later one.
BASELINE_WINDOWS = [WINDOW_LABELS[0]]


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
    """The calendar-month label of an acquisition day, or ``None`` if outside."""
    for label, start, end in WINDOWS:
        if start <= day < end:
            return label
    return None


# ---------------------------------------------------------------------------
# Rows and the per-month blocks
# ---------------------------------------------------------------------------
def month_rows(csv_rows):
    """The acquisitions of the record with the six echo-mask dimensions."""
    out = []
    for r in csv_rows:
        day = core.day_of(r)
        state = core.state_of(day) if day else None
        if day is None or state is None:            # outside the record
            continue
        d = {"day": day, "state": state, "window": window_of(day)}
        for k in FINGERPRINT_DIMS:
            v = r.get(k)
            d[k] = None if v is None else float(v)
        if d["window"] is not None:
            out.append(d)
    return out


def window_block(rows):
    """Per month: the count, the days, the epoch(s) and ``st._stats`` of the 6D."""
    out = {}
    for label, start, end in WINDOWS:
        sel = [r for r in rows if r["window"] == label]
        states = sorted({r["state"] for r in sel})
        out[label] = {
            "label": label,
            "window": [start, end],
            "states": states,
            "state": "+".join(states) if states else None,
            "n": len(sel),
            "days": sorted({r["day"] for r in sel}),
            "date_range": ([min(r["day"] for r in sel), max(r["day"] for r in sel)]
                           if sel else None),
            "channels": {k: st._stats([r.get(k) for r in sel])
                         for k in FINGERPRINT_DIMS},
        }
    return out


def month_classes(blocks):
    """The month labels by class: pure pre, mixed, rebuild, empty (no acquisition)."""
    states = lambda lab: blocks[lab]["states"]
    return {
        "pre": [lab for lab, _s, _e in WINDOWS if states(lab) == [PRE]],
        "mixed": [lab for lab, _s, _e in WINDOWS if len(states(lab)) > 1],
        "rebuild": [lab for lab, _s, _e in WINDOWS if states(lab) == [REBUILD]],
        "empty": [lab for lab, _s, _e in WINDOWS if not states(lab)],
    }


def delta_blocks(rows, classes):
    """The collapse month against its baselines — per dimension, three views.

    ``vs_prev`` compares the collapse month with the last month before it
    (adjacent calendar buckets); ``vs_pre_pool`` with every pure ``pre`` month
    pooled; ``vs_rebuild`` is the paper's epoch contrast (pooled ``rebuild`` vs
    pooled ``pre``) on the same months. Positive Cliff's delta => the second set
    is larger. The collapse month is *mixed* (``pre`` + ``single``) and the record
    has a four-month gap after it — both are properties of the record.
    """
    pre_pool = [r for r in rows if r["state"] == PRE]
    prev = [r for r in rows if r["window"] == classes["pre"][-1]]
    collapse = [r for r in rows if r["window"] == COLLAPSE_MONTH]
    rebuild_pool = [r for r in rows if r["state"] == REBUILD]
    out = {}
    for k in FINGERPRINT_DIMS:
        out[k] = {
            "vs_prev": st.delta_block([r[k] for r in prev],
                                      [r[k] for r in collapse],
                                      f"{COLLAPSE_MONTH} vs {classes['pre'][-1]}"),
            "vs_pre_pool": st.delta_block([r[k] for r in pre_pool],
                                          [r[k] for r in collapse],
                                          f"{COLLAPSE_MONTH} vs pooled pre"),
            "vs_rebuild": st.delta_block([r[k] for r in pre_pool],
                                         [r[k] for r in rebuild_pool],
                                         "pooled rebuild vs pooled pre"),
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
    """``{month: [median / baseline median per dim]}`` — the 6D panel data."""
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
    classes = month_classes(blocks)
    res = {"rows": rows, "blocks": blocks, "classes": classes,
           "deltas": delta_blocks(rows, classes)}
    res["baseline"] = base_median(blocks)
    res["fingerprint"] = fingerprint(blocks, res["baseline"])
    res["n_rows"] = len(rows)
    res["partition"] = {"rows": len(rows),
                        "per_window": {w[0]: blocks[w[0]]["n"] for w in WINDOWS}}
    res["runtime_s"] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def _runs(labels, subset):
    """The contiguous index runs of ``labels`` whose entries lie in ``subset``."""
    runs = []
    for i, lab in enumerate(labels):
        if lab not in subset:
            continue
        if runs and runs[-1][1] == i - 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    return runs


def draw_epochs(ax, classes):
    """Tint the two epochs, grey the acquisition gap, cut at the collapse day."""
    for key, color in (("pre", "#4C72B0"), ("mixed", "#4C72B0"),
                       ("rebuild", "#DD8452")):
        for i, j in _runs(WINDOW_LABELS, set(classes[key])):
            ax.axvspan(i - 0.5, j + 0.5, color=color, alpha=0.055, lw=0)
    for i, j in _runs(WINDOW_LABELS, set(classes["empty"])):
        ax.axvspan(i - 0.5, j + 0.5, color="#888888", alpha=0.12, lw=0)
        ax.text(0.5 * (i + j), 0.5, "no acquisitions", rotation=90,
                transform=ax.get_xaxis_transform(), ha="center", va="center",
                fontsize=5.4, color="#999999")
    ax.axvline(COLLAPSE_CUT_X, color="#888888", ls="--", lw=0.9)


def time_axis(ax):
    """Month-axis ticks (every third month) and the collapse cut."""
    ticks = [i for i in range(len(WINDOW_LABELS)) if i % 3 == 0]
    if (len(WINDOW_LABELS) - 1) not in ticks:
        ticks.append(len(WINDOW_LABELS) - 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([WINDOW_LABELS[i] for i in ticks], rotation=90,
                       fontsize=6.0)
    ax.set_xlim(-0.9, len(WINDOW_LABELS) - 0.1)
    ax.axvline(COLLAPSE_CUT_X, color="#888888", ls="--", lw=0.9)


def month_line(ax, rows, key, classes, fs=5.2, linthresh=None):
    """One month-by-month panel: per-acquisition dots, median and IQR per month.

    The per-month acquisition count is written only on every third month (and the
    last), so 28 numbers do not collide above a narrow panel. Empty months break
    the median line by construction (their median is NaN).
    """
    xs = np.arange(len(WINDOW_LABELS))
    med, q25, q75 = [], [], []
    for i, label in enumerate(WINDOW_LABELS):
        vals = [r[key] for r in rows
                if r["window"] == label and r.get(key) is not None]
        if vals:
            ax.plot(i + jitter(len(vals), 0.16, 300 + 7 * i), vals, "o", ms=1.9,
                    alpha=0.32, color="#8A8A8A", markeredgewidth=0)
            lo, mid, hi = np.percentile(vals, [25, 50, 75])
            if i % 3 == 0 or i == len(WINDOW_LABELS) - 1:
                ax.text(i, 0.995, f"{len(vals)}",
                        transform=ax.get_xaxis_transform(), ha="center",
                        va="top", fontsize=fs, color="#777777")
        else:
            lo = mid = hi = np.nan
        med.append(mid)
        q25.append(lo)
        q75.append(hi)
    ax.fill_between(xs, np.asarray(q25, float), np.asarray(q75, float),
                    color="#4C72B0", alpha=0.16, lw=0)
    ax.plot(xs, np.asarray(med, float), "-o", ms=2.6, lw=1.1, color="#2E4E7E")
    if linthresh is not None:
        ax.set_yscale("symlog", linthresh=linthresh)
    draw_epochs(ax, classes)
    time_axis(ax)
    ax.margins(y=0.14)


def fig_months(res, out_dir=None):
    rows = res["rows"]
    blocks = res["blocks"]
    classes = res["classes"]
    fp, deltas = res["fingerprint"], res["deltas"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 160 * MM))

    # (a)-(d) the four echo-mask dimensions, one month-by-month panel each
    for letter, k in zip(["a", "b", "c", "d"], STRIP_DIMS):
        i = STRIP_DIMS.index(k)
        ax = axes[i // 3][i % 3]
        month_line(ax, rows, k, classes, linthresh=STRIP_SCALES.get(k))
        b = deltas[k]["vs_prev"]
        title = core.CHANNEL_TITLE.get(k, k)
        if "delta" in b:
            title += (f"\nΔ({COLLAPSE_MONTH}−{classes['pre'][-1]}) = "
                      f"{b['delta']:+.2f}")
        ax.set_title(f"({letter}) {title}", loc="left", fontsize=8.0)
        ax.set_ylabel(core.CHANNEL_AXIS.get(k, k), fontsize=7.0)
    # the collapse day, on the exact fractional cut inside the month: the cut is
    # not a bucket boundary here, so the label sits just past the dashed line
    axes[0][0].text(COLLAPSE_CUT_X + 0.22, 0.34, f"collapse\n{COLLAPSE_DAY}",
                    rotation=90, transform=axes[0][0].get_xaxis_transform(),
                    ha="center", va="center", fontsize=5.4, color="#444444")

    # (e) the (P, D) trajectory, one dot per acquisition coloured by month
    ax = axes[1][1]
    cmap = plt.cm.viridis
    for i, label in enumerate(WINDOW_LABELS):
        sel = [r for r in rows if r["window"] == label
               and r.get("P") is not None and r.get("D") is not None]
        if not sel:
            continue
        ax.plot([r["P"] for r in sel], [r["D"] for r in sel], "o", ms=2.4,
                alpha=0.55, color=cmap(i / max(len(WINDOW_LABELS) - 1, 1)),
                markeredgewidth=0)
    px = [blocks[l]["channels"]["P"]["median"] if blocks[l]["channels"]["P"]
          else np.nan for l in WINDOW_LABELS]
    dy = [blocks[l]["channels"]["D"]["median"] if blocks[l]["channels"]["D"]
          else np.nan for l in WINDOW_LABELS]
    ax.plot(px, dy, "-", lw=1.0, color="#333333", alpha=0.7)
    ax.plot(px, dy, "o", ms=3.2, mfc="none", mec="#333333", mew=0.8, ls="none")
    ci = len(classes["pre"])                    # index of the mixed collapse month
    ax.plot([px[ci]], [dy[ci]], "o", ms=5.5, color="#C44E52")
    ax.annotate(COLLAPSE_MONTH, (px[ci], dy[ci]), textcoords="offset points",
                xytext=(7, -7), fontsize=5.6, color="#C44E52")
    ax.set_xlabel("P  (mask persistence)")
    ax.set_ylabel("D  (mask density A / bbox)")
    ax.set_yscale("symlog", linthresh=0.01)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(0, len(WINDOW_LABELS) - 1))
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.02)
    cb.set_ticks([0, len(WINDOW_LABELS) - 1])
    cb.set_ticklabels([WINDOW_LABELS[0], WINDOW_LABELS[-1]])
    cb.ax.tick_params(labelsize=5.6)
    cb.ax.set_title("month", fontsize=6.0, pad=3)
    panel(ax, "e", "the (P, D) trajectory per month")

    # (f) the 6D fingerprint: a line per dimension over the whole record
    ax = axes[1][2]
    xs = np.arange(len(WINDOW_LABELS))
    for j, k in enumerate(FINGERPRINT_DIMS):
        ys = [(fp[l][j] if (fp[l][j] is not None and fp[l][j] > 0) else np.nan)
              for l in WINDOW_LABELS]
        ax.plot(xs, ys, "-", lw=1.1, label=k, color=plt.cm.tab10(j % 10))
    ax.axhline(1.0, color="#888888", lw=0.8)
    time_axis(ax)
    ax.set_yscale("log")
    ax.set_ylabel(f"median / {BASELINE_WINDOWS[0]} median")
    ax.legend(frameon=True, framealpha=0.85, edgecolor="none", ncol=3,
              fontsize=5.8, loc="upper right")
    panel(ax, "f", f"the 6D fingerprint (baseline = {BASELINE_WINDOWS[0]})")

    pre_labs = classes["pre"] + classes["mixed"]
    reb_labs = classes["rebuild"]
    span = lambda labs: f"{labs[0]}..{labs[-1]}" if labs else "–"
    fig.suptitle(f"(KDLO) echo-mask dimensions month by month over the record — "
                 f"{res['n_rows']} acquisitions, "
                 f"{WINDOW_LABELS[0]}..{WINDOW_LABELS[-1]}; pre-collapse "
                 f"{span(pre_labs)} (blue) vs. rebuild {span(reb_labs)} (orange)",
                 fontsize=8.6, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.960))
    return save_fig(fig, "fig_kdlo_ocv_months", out_dir)


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


def report(res):
    blocks, deltas = res["blocks"], res["deltas"]
    classes, base = res["classes"], res["baseline"]
    fp = res["fingerprint"]
    n_pre = sum(1 for r in res["rows"] if r["state"] == PRE)
    n_single = sum(1 for r in res["rows"] if r["state"] == "single")
    n_reb = sum(1 for r in res["rows"] if r["state"] == REBUILD)
    md = []
    md.append("# KDLO media tower (Garden City SD) — the echo-mask dimensions "
              "month by month over the record")
    md.append("")
    md.append(f"**Site** KDLO media tower (Garden City SD) · **record** "
              f"{WINDOW_LABELS[0]} … {WINDOW_LABELS[-1]} · **acquisitions** "
              f"{res['n_rows']}")
    md.append("")
    md.append("```")
    md.append("x     = [gamma2, P, D, A, F, S]      echo-mask vector")
    md.append("mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2")
    for lab, win, note, cnt in (
            (PRE, f"{core.STATE_WINDOWS[PRE][0]}..{core.STATE_WINDOWS[PRE][1]}",
             "standing mast", f"{n_pre} acquisitions"),
            ("single", FIRST_COLLAPSED_DAY, "in neither epoch",
             f"{n_single} acquisition"),
            (REBUILD,
             f"{core.STATE_WINDOWS[REBUILD][0]}..{core.STATE_WINDOWS[REBUILD][1]}",
             "during the rebuild", f"{n_reb} acquisitions")):
        md.append(f"{lab:<7} {win:<23} {note:<19} {cnt}")
    md.append("```")
    md.append("")
    md.append("## Why this is not a 24-month series")
    md.append("")
    md.append("The other site companions roll 24 one-month buckets before their "
              "event because those records carry a dense pre-collapse stretch. "
              "The KDLO record is 60 acquisitions in two **disjoint epochs** "
              f"(``pre`` {core.STATE_WINDOWS[PRE][0]}..{core.STATE_WINDOWS[PRE][1]} "
              f"and ``rebuild`` {core.STATE_WINDOWS[REBUILD][0]}.."
              f"{core.STATE_WINDOWS[REBUILD][1]}), so the axis here is the whole "
              f"record by **calendar month** ({len(WINDOW_LABELS)} buckets, "
              f"{len(classes['empty'])} of them with no acquisition). Two "
              "properties of the record are kept visible rather than smoothed "
              f"away: the collapse month **{COLLAPSE_MONTH} is mixed** (it holds "
              f"both the last ``pre`` acquisitions, to {COLLAPSE_DAY}, and the "
              f"``single`` {FIRST_COLLAPSED_DAY} that belongs to neither epoch), "
              "and the record has a four-month acquisition gap between the two "
              "epochs.")
    md.append("")
    md.append(f"## The record, month by month "
              f"({len(WINDOW_LABELS)} calendar buckets)")
    md.append("")
    body = []
    for label, _s, _e in WINDOWS:
        b = blocks[label]
        body.append([label, b["state"] or "(no acquisition)", b["n"],
                     b["date_range"][0] if b["date_range"] else "–",
                     b["date_range"][1] if b["date_range"] else "–",
                     len(b["days"])])
    md += md_table(["month", "epoch(s)", "acquisitions", "first day", "last day",
                    "days"], body)
    md.append("")
    md.append(f"`pre` months: {', '.join(classes['pre'])}; collapse month: "
              f"{COLLAPSE_MONTH}; empty: {', '.join(classes['empty'])}; "
              f"`rebuild` months: {classes['rebuild'][0]} … "
              f"{classes['rebuild'][-1]}.")
    md.append("")
    md.append("## Median of each dimension per month")
    md.append("")
    body = []
    for label, _s, _e in WINDOWS:
        b = blocks[label]
        row = [label, b["n"], b["state"] or "–"]
        for k in FINGERPRINT_DIMS:
            c = b["channels"][k]
            row.append(_fmt(c["median"]) if c else "–")
        body.append(row)
    md += md_table(["month", "n", "epoch"] + FINGERPRINT_DIMS, body)
    md.append("")
    md.append(f"`baseline` = the `{BASELINE_WINDOWS[0]}` (first record month) "
              "medians — the 1.0 line of panel f: "
              + ", ".join(f"`{k}` {_fmt(base.get(k))}" for k in FINGERPRINT_DIMS)
              + ".")
    md.append("")
    md.append("## The 6D fingerprint at the cut and at the epoch edge "
              "(month / baseline)")
    md.append("")
    body = []
    for label in [BASELINE_WINDOWS[0], classes["pre"][-1], COLLAPSE_MONTH,
                  classes["rebuild"][0], WINDOW_LABELS[-1]]:
        body.append([label] + [_fmt(fp[label][j])
                               for j in range(len(FINGERPRINT_DIMS))])
    md += md_table(["month"] + FINGERPRINT_DIMS, body)
    md.append("")
    md.append("## The collapse month against its baselines")
    md.append("")
    body = []
    for k in FINGERPRINT_DIMS:
        row = [k]
        for key in ("vs_prev", "vs_pre_pool", "vs_rebuild"):
            b = deltas[k][key]
            row.append(_fmt(b.get("delta")))
            row.append("yes" if b.get("ci_excludes_zero") else "no")
        body.append(row)
    md += md_table(["dimension", f"δ {COLLAPSE_MONTH} vs prev", "CI≠0",
                    f"δ {COLLAPSE_MONTH} vs pooled pre", "CI≠0",
                    "δ rebuild vs pre", "CI≠0"], body)
    md.append("")
    md.append(f"`prev` = the last month before the collapse "
              f"(`{classes['pre'][-1]}`); `pooled pre` = the {len(classes['pre'])} "
              "pure `pre` months; `rebuild vs pre` is the paper's epoch contrast. "
              "Cliff's delta > 0 means the right-hand set is larger. Full CIs, "
              "medians and SDS are in the JSON (`deltas[dim].*`).")
    md.append("")
    md.append("## Reading (descriptive — not a claim of damage)")
    md.append("")
    n_prev = sum(1 for k in FINGERPRINT_DIMS
                 if deltas[k]["vs_prev"].get("ci_excludes_zero"))
    n_epoch = sum(1 for k in FINGERPRINT_DIMS
                  if deltas[k]["vs_rebuild"].get("ci_excludes_zero"))
    last_pre = ", ".join(f"`{k}` {_fmt(fp[classes['pre'][-1]][j])}"
                         for j, k in enumerate(FINGERPRINT_DIMS))
    first_reb = ", ".join(f"`{k}` {_fmt(fp[classes['rebuild'][0]][j])}"
                          for j, k in enumerate(FINGERPRINT_DIMS))
    md.append(f"The {len(classes['pre'])} pure `pre` months carry **{n_pre}** "
              "acquisitions (two to three per month) and the collapse month "
              f"**{blocks[COLLAPSE_MONTH]['n']}**, so every monthly interval is "
              "very wide: this is a descriptive companion, not a test.")
    md.append("")
    md.append(f"Relative to the first record month `{BASELINE_WINDOWS[0]}` (= 1.0) "
              f"the last pure `pre` month `{classes['pre'][-1]}` sits at "
              f"{last_pre} and the first `rebuild` month "
              f"`{classes['rebuild'][0]}` at {first_reb}. Against the last month "
              "before the collapse the mixed bucket clears the bootstrap CI in "
              f"**{n_prev}** of the six dimensions — and that bucket is not a "
              "clean post month. The pooled epoch contrast clears it in "
              f"**{n_epoch}** of six. The axis is read as the record's own "
              "co-variation over two disjoint windows, **never** as damage: the "
              "reconstruction, the season, the acquisition geometry and the "
              "pre/post mix move with the collapse date as well.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/kdlo/kdlo_ocv_channels_csv.py       # build + verify the table")
    md.append("python3 code/kdlo/fig_kdlo_ocv_months.py          # this figure")
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
    """The result JSON — the per-month blocks and the deltas, not the rows."""
    return {
        "site": "KDLO media tower (Garden City SD)",
        "record": [WINDOW_LABELS[0], WINDOW_LABELS[-1]],
        "generated_by": "code/kdlo/fig_kdlo_ocv_months.py",
        "states": dict(core.STATE_LABEL),
        "state_windows": {k: list(v) for k, v in core.STATE_WINDOWS.items()},
        "collapse_day": COLLAPSE_DAY,
        "first_collapsed_day": FIRST_COLLAPSED_DAY,
        "singleton_day": core.SINGLETON_DAY,
        "features": list(core.FEATURES),
        "strip_dims": list(STRIP_DIMS),
        "baseline_windows": list(BASELINE_WINDOWS),
        "months": len(WINDOW_LABELS),
        "classes": res["classes"],
        "windows": [[w[0], w[1], w[2]] for w in WINDOWS],
        "n_rows": res["n_rows"],
        "partition": res["partition"],
        "baseline": res["baseline"],
        "fingerprint": res["fingerprint"],
        "blocks": res["blocks"],
        "deltas": res["deltas"],
        "figures": ["figures/kdlo/fig_kdlo_ocv_months.png"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — KDLO media tower, {len(WINDOW_LABELS)} calendar months "
          f"{WINDOW_LABELS[0]}..{WINDOW_LABELS[-1]} …")
    res = compute()
    print(f"  {res['n_rows']} acquisitions "
          f"({res['partition']['per_window']}), {res['runtime_s']:.1f} s")

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
    raise SystemExit(main())
