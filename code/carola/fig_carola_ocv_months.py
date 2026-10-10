#!/usr/bin/env python3
"""Carola OCV-Paper — the echo-mask dimensions over the 24 pre-collapse months.

Question: what do the six dimensions of the **echo-mask vector**
``x = [gamma2, P, D, A, F, S]`` do over the **24 one-month windows before** the
collapse (``2024-09-11``), and where does the first month *after* it sit on that
axis? The A-E figures of ``fig_carola_ocv_paper.py`` compare the two pooled
states (``pre-collapse (healthy)`` vs. ``post-collapse``); this companion keeps
the time axis instead: it rolls one-month buckets from the record start to the
event, so the reader sees *when* the mask moved across the whole two-year record,
not only *that* it moved between two pooled pools.

The windows are rolling one-month buckets **anchored at the event** (the collapse
date itself opens M+1, so M-1 is pure pre and M+1 pure post):

    M-24  [2022-09-11, 2022-10-11)     pre-collapse (healthy)  first record month
    ...
    M-1   [2024-08-11, 2024-09-11)     pre-collapse (healthy)  last
    M+1   [2024-09-11, 2024-10-11)     post-collapse           terminal contrast

Approach: this is a **standalone companion** of ``fig_carola_ocv_paper.py`` (like
``../bautzen/fig_bautzen_ocv_amp_phase.py`` is a companion of the Bautzen paper
script): it reads only the committed channel table
``data/carola/carola_ocv_channels.csv``, recomputes per-month aggregates with
``carola_ocv_stats``, draws one figure and writes a small JSON + English report.
It adds **no** number to the A-E pin contract and touches none of it.

Figure (2 x 3, 183 mm wide):
  a  gamma2 by month   (per-chip dots, monthly median + IQR, the collapse cut)
  b  P       by month  (the mask-persistence channel)
  c  D       by month  (mask density A / bbox)
  d  A       by month  (mask size, masked px)
  e  the (P, D) plane  (one dot per chip, coloured by month)
  f  the 6D fingerprint (each month's six-dim median relative to the M-24 baseline,
     drawn as a line per dimension over the 24-month axis)

The figure is descriptive and carries the site's standing caveat: both sides live
on the *same* bridge, so season, weather, traffic and pipeline generation move
with the collapse date as well. The report states that rather than a mechanism.

Usage:
  python3 code/carola/fig_carola_ocv_months.py
  python3 code/carola/fig_carola_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
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

import carola_ocv_core as core   # noqa: E402
import carola_ocv_stats as st    # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_carola_ocv_months.json")
OUT_MD = os.path.join(FIGDIR, "fig_carola_ocv_months.md")

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

N_MONTHS = 24          # the pre-collapse series: M-24 … M-1


def shift_month(day, k):
    """``YYYY-MM-DD`` shifted by ``k`` whole calendar months (anchor day kept).

    The event days of every site in this package are <= 28, so the anchor day
    survives every month; the helper is pure integer arithmetic (no dateutil).
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


# The pre-collapse series (M-24 .. M-1) plus the first month after the collapse.
WINDOWS = make_windows(core.COLLAPSE_DATE)
WINDOW_LABELS = [w[0] for w in WINDOWS]
# The pre-collapse monthly labels only (no M+1): the pooled pre side.
PRE_LABELS = [w[0] for w in WINDOWS if w[0] != "M+1"]
WINDOW_STATE = {w[0]: (core.PRE if w[0] != "M+1" else core.POST) for w in WINDOWS}
# Where the collapse cut and (for sites with an earlier baseline epoch) the
# record gap fall on the month axis.
COLLAPSE_CUT_X = WINDOW_LABELS.index("M+1") - 0.5
# The panels (a)-(d): the four echo-mask dimensions in the order of the question.
STRIP_DIMS = ["gamma2", "P", "D", "A"]
# A couple of dimensions span decades across a two-year record; those panels get
# a symlog y axis (linthresh) so the monthly median line stays readable.
STRIP_SCALES = {"D": 0.01}
# The 6D fingerprint runs over the package's own feature order.
FINGERPRINT_DIMS = list(core.FEATURES)
# The fingerprint baseline is the first record month (M-24): the only pre window
# whose value does not depend on any later one.
BASELINE_WINDOWS = [PRE_LABELS[0]]


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
    """The chips that fall in any of the windows (record start … M+1)."""
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
    windows the collapse date separates); ``vs_pre_pool`` compares it with all
    ``N_MONTHS`` pre-collapse months pooled. Positive Cliff's delta => M+1 larger.
    """
    pre_pool = [r for r in rows if r["window"] in PRE_LABELS]
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
    """Per-dimension mean of the baseline-window medians (the fingerprint base)."""
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
def time_axis(ax, gap_x=None):
    """Month-axis ticks (every third month), the collapse cut and any record gap."""
    ticks = [i for i in range(len(WINDOW_LABELS)) if i % 3 == 0]
    if (len(WINDOW_LABELS) - 1) not in ticks:
        ticks.append(len(WINDOW_LABELS) - 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([WINDOWS[i][1][:7] for i in ticks], rotation=90,
                       fontsize=6.0)
    ax.set_xlim(-1.0, len(WINDOW_LABELS) + 1.2)
    ax.axvline(COLLAPSE_CUT_X, color="#888888", ls="--", lw=0.9)
    if gap_x is not None:
        ax.axvline(gap_x, color="#BBBBBB", ls=":", lw=0.9)


def month_line(ax, rows, key, fs=5.2, gap_x=None, linthresh=None):
    """One month-by-month panel: per-chip dots, the monthly median and its IQR.

    The per-month chip count is written only on every sixth month (and the last),
    so 24 numbers do not collide above a narrow panel.
    """
    xs = np.arange(len(WINDOW_LABELS))
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
    ax.fill_between(xs, np.asarray(q25, float), np.asarray(q75, float),
                    color="#4C72B0", alpha=0.16, lw=0)
    ax.plot(xs, np.asarray(med, float), "-o", ms=2.4, lw=1.1, color="#2E4E7E")
    if linthresh is not None:
        ax.set_yscale("symlog", linthresh=linthresh)
    time_axis(ax, gap_x)
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
        b = deltas[k]["vs_M-1"]
        title = core.CHANNEL_AXIS.get(k, k)
        if "delta" in b:
            title += f"\nΔ(M+1−M−1) = {b['delta']:+.2f}"
        ax.set_title(f"({letter}) {title}", loc="left", fontsize=8.0)
    # one collapse label, in the clear right margin of the first panel
    axes[0][0].text(len(WINDOW_LABELS) + 0.35, 0.5,
                    f"collapse {core.COLLAPSE_DATE}", rotation=90,
                    transform=axes[0][0].get_xaxis_transform(), ha="center",
                    va="center", fontsize=6.0, color="#555555")

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
    ax.plot(px[:-1], dy[:-1], "-", lw=1.0, color="#333333", alpha=0.7)
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

    # (f) the 6D fingerprint: a line per dimension over the 24-month axis
    ax = axes[1][2]
    xs = np.arange(len(WINDOW_LABELS))
    for j, k in enumerate(FINGERPRINT_DIMS):
        ys = [(fp[l][j] if (fp[l][j] is not None and fp[l][j] > 0) else np.nan)
              for l in WINDOW_LABELS]
        ax.plot(xs, ys, "-", lw=1.1, label=k, color=plt.cm.tab10(j % 10))
    ax.axhline(1.0, color="#888888", lw=0.8)
    time_axis(ax)
    ax.set_yscale("log")
    ax.set_ylabel("median / baseline median")
    ax.legend(frameon=True, framealpha=0.85, edgecolor="none", ncol=3,
              fontsize=5.8, loc="upper right")
    panel(ax, "f", f"the 6D fingerprint (baseline = {BASELINE_WINDOWS[0]})")

    fig.suptitle(f"(Carola) echo-mask dimensions over the {N_MONTHS} months before "
                 f"the collapse — {res['n_rows']} chips, M-{N_MONTHS}..M-1 of "
                 f"{core.COLLAPSE_DATE} plus M+1", fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.960))
    return save_fig(fig, "fig_carola_ocv_months", out_dir)


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
    blocks, deltas, base = res["blocks"], res["deltas"], res["baseline"]
    fp = res["fingerprint"]
    md = []
    md.append(f"# Carola Bruecke, Dresden — the echo-mask dimensions over the "
              f"{N_MONTHS} months before the collapse")
    md.append("")
    md.append(f"**Site** {core.SITE} · **event** {core.COLLAPSE_DATE} · "
              f"**chips** {res['n_rows']} (windows M-{N_MONTHS} … M-1 plus M+1)")
    md.append("")
    md.append("```")
    md.append("x     = [gamma2, P, D, A, F, S]      echo-mask vector")
    md.append("mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2")
    md.append(f"M-{N_MONTHS:<3} [{WINDOWS[0][1]}, {WINDOWS[0][2]})       "
              "pre-collapse (healthy)  first record month")
    md.append(f"M-1   [{WINDOWS[-2][1]}, {WINDOWS[-2][2]})       "
              "pre-collapse (healthy)  last")
    md.append(f"M+1   [{WINDOWS[-1][1]}, {WINDOWS[-1][2]})       post-collapse")
    md.append("```")
    md.append("")
    md.append(f"## The {N_MONTHS} pre-collapse months and the month after")
    md.append("")
    md.append(f"Rolling one-month buckets anchored at the collapse date, from the "
              f"record start (M-{N_MONTHS}) up to the event (M-1); the event itself "
              f"opens M+1, so the pre side is the {N_MONTHS} M-* windows and the "
              "post side is M+1 alone. This is a companion of "
              "`fig_carola_ocv_paper.py` (figures A–E): it adds one figure and "
              "does not touch the A–E pin contract.")
    md.append("")
    body = []
    for label, s, e in WINDOWS:
        b = blocks[label]
        body.append([label, f"[{s}, {e})", core.STATE_SHORT[b["state"]], b["n"],
                     b["date_range"][0] if b["date_range"] else "–",
                     b["date_range"][1] if b["date_range"] else "–",
                     len(b["days"])])
    md += md_table(["month", "range", "state", "chips", "first day", "last day",
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
    md.append(f"`baseline` = the M-{N_MONTHS} monthly medians (the 1.0 line of "
              "panel f): "
              + ", ".join(f"`{k}` {_fmt(base.get(k))}" for k in FINGERPRINT_DIMS)
              + ".")
    md.append("")
    md.append("## The 6D fingerprint at the cut (month / baseline)")
    md.append("")
    body = []
    for label in [PRE_LABELS[0], PRE_LABELS[len(PRE_LABELS) // 2],
                  PRE_LABELS[-1], "M+1"]:
        body.append([label] + [_fmt(fp[label][j])
                               for j in range(len(FINGERPRINT_DIMS))])
    md += md_table(["month"] + FINGERPRINT_DIMS, body)
    md.append("")
    md.append("## What moved across the cut (M+1 vs M-1)")
    md.append("")
    body = []
    for k in FINGERPRINT_DIMS:
        b = deltas[k]["vs_M-1"]
        if "delta" not in b:
            body.append([k, "–", "–", "–", "–", "–", "–"])
            continue
        ci = b.get("delta_ci95")
        body.append([k, _fmt(b["median"][0]), _fmt(b["median"][1]),
                     _fmt(b["delta"]), _fmt(b["sds"]),
                     "yes" if b.get("ci_excludes_zero") else "no",
                     f"[{_fmt(ci[0])}, {_fmt(ci[1])}]" if ci else "–"])
    md += md_table(["dim", "median M-1", "median M+1", "Cliff's delta", "SDS",
                    "CI95 excludes 0", "CI95"], body)
    md.append("")
    md.append(f"Cliff's delta > 0 means M+1 larger than M-1. The same comparison "
              f"against all {N_MONTHS} pre months pooled is in the JSON "
              "(`deltas[dim].vs_pre_pool`).")
    md.append("")
    md.append("## Reading (descriptive — not a claim of damage)")
    md.append("")
    pre_chips = sum(blocks[l]["n"] for l in PRE_LABELS)
    md.append(f"The {N_MONTHS} pre-collapse months carry **{pre_chips}** chips "
              f"(median {int(np.median([blocks[l]['n'] for l in PRE_LABELS]))} per "
              f"month) and the first post month **{blocks['M+1']['n']}**, so the "
              "monthly intervals are wide and this is a descriptive companion, "
              "not a test.")
    md.append("")
    last_pre = ", ".join(f"`{k}` {_fmt(fp[PRE_LABELS[-1]][j])}"
                         for j, k in enumerate(FINGERPRINT_DIMS))
    first_post = ", ".join(f"`{k}` {_fmt(fp['M+1'][j])}"
                           for j, k in enumerate(FINGERPRINT_DIMS))
    n_excl = sum(1 for k in FINGERPRINT_DIMS
                 if deltas[k]["vs_M-1"].get("ci_excludes_zero"))
    md.append(f"Relative to the record-start baseline `{BASELINE_WINDOWS[0]}` "
              f"(= 1.0) the **last** pre month M-1 sits at {last_pre}; the first "
              f"post month M+1 sits at {first_post}. The pre-collapse record is "
              "therefore **not** a flat line — season, weather, traffic and "
              "pipeline generation move the six dimensions month to month well "
              "before the event. Against that spread the single adjacent contrast "
              f"(M+1 vs M-1) clears the bootstrap CI in **{n_excl}** of the six "
              "dimensions; the full per-month trajectory is in the JSON "
              "(`blocks`, `fingerprint`, `deltas`).")
    md.append("")
    md.append("The month axis therefore does **not** reproduce the pooled A–E "
              "headline by itself; it shows how much of the pre-collapse record "
              "already moves on its own and where the collapse sits inside it. "
              "Both sides live on the *same* bridge, so season, weather, traffic "
              "and pipeline generation move with the collapse date as well — the "
              "panels are read as that co-variation, **never** as damage.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/carola/carola_ocv_channels_csv.py       # build + verify the table")
    md.append("python3 code/carola/fig_carola_ocv_months.py          # this figure")
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
        "generated_by": "code/carola/fig_carola_ocv_months.py",
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
        "figures": ["figures/carola/fig_carola_ocv_months.png"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, windows M-{N_MONTHS}..M+1 of "
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

