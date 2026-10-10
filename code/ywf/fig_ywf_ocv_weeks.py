#!/usr/bin/env python3
"""YWF OCV-Weeks — the OCV time series in the weeks before the collapse.

Question: what does the observability core vector ``OCV = [gamma2, P, D]`` of the
Yeongdeok Wind Farm (Changpo Wind Power Complex, Samgye-ri) do **inside** the
pre-collapse state — in the ~26 weeks between the first committed acquisition
(2025-08-06) and the 100 m monopole of Unit 21 falling onto a public road on
2026-02-02? The paper figure (``fig_ywf_ocv_paper.py``) pools that window against
the 6 post-event chips; this companion keeps the time axis and opens the window
itself.

What this figure claims, and what it does not
---------------------------------------------
  * **the subject is the pre window only.** 27 unique chips on 26 acquisition
    days, one per 6-day repeat step (the ascending and the descending pass
    alternate), with four 12-day steps inside the window (12 Aug -> 24 Aug,
    24 Aug -> 5 Sep, 3 Jan -> 15 Jan, 15 Jan -> 27 Jan). The collapse line sits
    at the right edge (``W-0``) and **no** post-event chip is drawn: this is not
    a contrast, and it adds no number to the A-E pin contract of the paper
    script.
  * **the coordinate is not stationary, and the coverage moves with it.** The
    window walks from late summer into deep winter — temperature falls from
    33.2 C to 0.6 C — and the echo itself is not equally present along it. Any
    "before vs. after" statement computed from the pooled window is therefore
    read against the co-variates of panel (f), never as a change at the event.
  * **the OCV values themselves show no monotone drift.** Over the 16
    echo-bearing chips of the window, Spearman against acquisition time gives
    gamma2 +0.135 (p 0.62), P -0.17 (p 0.52) and D -0.25 (p 0.35): the
    non-stationarity of this window lives in the *coverage* and in the weather,
    not in a trend of the channels. The figure reports both numbers instead of
    picking the flattering one.
  * **``P`` is a window-wide reference, not a per-chip quantity.** It is the mean
    pixel frequency inside the own mask over the **17** echo-bearing chips of the
    whole committed record (the cache's own reference count), so the chips of
    this figure are not independent draws of it — the report says so.

How the time axis is built
--------------------------
  x = -(2026-02-02 - day).days / 7      weeks before the collapse, 0 at the event
  the four 12-day steps are drawn as bands and **never** interpolated over
  the month starts are dotted, and every tick label carries the calendar day
  above the weeks before the collapse of that day (``Sep 1`` over ``22``)

Figures:
  a  the echo of every acquisition (``A``, the masked pixels of 49) as a stem,
     grey where the chip has no echo, coloured by echo mode, with the chip
     ``intensity`` on a twin log axis
  b  gamma2 over the 16 echo-bearing chips + its gap-aware 4-week rolling median
  c  P  (the same)
  d  D  (the same)
  e  the three OCV channels normalised to their pre-window median, one log axis
  f  what else moves over these weeks: temperature, wind speed and gust

Cross-checks: every run re-derives the window from the committed channel table
and the committed mask cache (the mask vector of all 17 echo rows, the per-day
timeline of the paper JSON, the coverage table and the echo rates of the
committed reference, the constant pipeline ``coherence`` column) and exits
non-zero if any of them deviates.

Usage:
  python3 code/ywf/fig_ywf_ocv_weeks.py                  # figure + JSON + report
  python3 code/ywf/fig_ywf_ocv_weeks.py --no-figures      # cross-checks only
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import ywf_ocv_core as core     # noqa: E402
import ywf_ocv_masks as masks   # noqa: E402
import ywf_ocv_stats as st      # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_ywf_ocv_weeks.json")
OUT_MD = os.path.join(FIGDIR, "fig_ywf_ocv_weeks.md")

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

EVENT = dt.date.fromisoformat(core.COLLAPSE_DATE)
OCV = list(core.OCV)
FEATURES = list(core.FEATURES)
# The channels with one dot per chip in the strip panels; ``A`` has its own stem
# panel because it is also the coverage flag (0 = no echo).
STRIP_DIMS = ["gamma2", "P", "D"]
TREND_DIMS = ["gamma2", "P", "D", "A"]
# The co-variates of panel (f). The pipeline's own ``coherence`` column is
# deliberately absent from the figure — it is one constant on the whole record —
# and is reported in the tables instead.
TEMP_KEY, WIND_KEYS = "temperature_c", ["wind_speed_ms", "wind_gust_ms"]
CONTROL_KEYS = ["intensity", "coherence", "coherence_gamma2", "wind_speed_ms",
                "temperature_c", "precipitation_mm", "wind_gust_ms",
                "relative_humidity_pct"]
# The columns whose *movement over the window* panel (f) and the trend table
# report: the two wind columns ride on one axis, temperature on the other. The
# pipeline's constant ``coherence`` is in ``CONTROL_KEYS`` only.
COVARIATE_KEYS = ["temperature_c", "wind_speed_ms", "wind_gust_ms", "intensity"]

CADENCE_DAYS = 6          # the repeat step of the two alternating orbits
GAP_DAYS = 12             # a step above the cadence is a gap in the window
ROLL_SPAN_WEEKS = 4.0     # the trailing window of the rolling medians
ROLL_MIN_N = 2            # a rolling median needs at least two chips
FLOAT_TOL = 1e-12

ECHO_MODE_COLORS = {"compact": "#4C72B0", "intermediate": "#55A868",
                    "distributed": "#C44E52"}
NO_ECHO_COLOR = "#B0B0B0"
CHANNEL_COLORS = {"gamma2": "#4C72B0", "P": "#55A868", "D": "#8172B3",
                  "A": "#937860"}
COVARIATE_COLORS = {"temperature_c": "#C44E52", "wind_speed_ms": "#4C72B0",
                    "wind_gust_ms": "#8172B3", "precipitation_mm": "#55A868",
                    "intensity": "#64B5CD"}
EVENT_COLOR = "#C44E52"
GAP_COLOR = "#E8E8E8"
MONTH_COLOR = "#DDDDDD"
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
          "Oct", "Nov", "Dec"]

# ---------------------------------------------------------------------------
# Small helpers (the shape of the family's companion scripts)
# ---------------------------------------------------------------------------
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


def day_date(day):
    """``YYYY-MM-DD`` -> ``datetime.date``."""
    return dt.date.fromisoformat(str(day))


def weeks_before(day, event=EVENT):
    """Weeks between an acquisition day and the collapse (positive = before)."""
    return round((event - day_date(day)).days / 7.0, 3)


def x_of(day):
    """The x coordinate: negative weeks before the collapse, 0 at the event."""
    return -weeks_before(day)


def orbit_dx(orbit):
    """A deterministic small offset for the two orbits of one day."""
    lab = str(orbit or "").upper()
    if lab.startswith("ASC"):
        return 0.055
    if lab.startswith("DESC"):
        return -0.055
    return 0.0


def row_view(r):
    """The columns of one chip this figure uses (``A`` is None without echo)."""
    d = {
        "id": r.get("id"), "day": core.day_of(r),
        "orbit": r.get("orbit_direction"), "status": r.get("status"),
        "traffic": r.get("traffic_load_label"),
        "echo": bool(r["echo_mask_present"]),
        "A": core._num(r.get("A")), "F": core._num(r.get("F")),
        "S": core._num(r.get("S")), "gamma2": core._num(r.get("gamma2_mask")),
        "P": core._num(r.get("P")), "D": core._num(r.get("D")),
        "echo_mode": r.get("echo_mode"),
        "coherence_masked_pixels": core._num(r.get("coherence_masked_pixels")),
    }
    for k in CONTROL_KEYS:
        d[k] = core._num(r.get(k))
    d["pre"] = bool(r["pre_collapse"])
    d["state"] = core.STATE_SHORT[core.state_of(d["day"])]
    d["wb"] = weeks_before(d["day"])
    return d


def values(rows, key):
    """The non-empty values of one column."""
    return [r[key] for r in rows if r.get(key) is not None]


def stat_block(rows, key):
    """``st._stats`` over one column of the given rows."""
    return st._stats(values(rows, key))


def median_of(rows, key):
    v = values(rows, key)
    return float(np.median(v)) if v else None


def sorted_chips(rows):
    """The chips in acquisition order (day, then orbit label)."""
    return sorted(rows, key=lambda d: (d["day"], str(d["orbit"] or "")))


# ---------------------------------------------------------------------------
# Aggregate blocks of the window
# ---------------------------------------------------------------------------
def gap_block(chips):
    """The steps inside the window that are longer than the 6-day cadence."""
    days = sorted({d["day"] for d in chips})
    out = []
    for a, b in zip(days, days[1:]):
        n = (day_date(b) - day_date(a)).days
        if n > CADENCE_DAYS:
            out.append({"after": a, "before": b, "days": n,
                        "x": [x_of(a), x_of(b)],
                        "wb": [weeks_before(a), weeks_before(b)]})
    return out


def iso_week_block(chips):
    """The Monday-anchored ISO-week buckets of the window (medians per channel)."""
    buckets = {}
    for d in chips:
        day = day_date(d["day"])
        monday = (day - dt.timedelta(days=day.weekday())).isoformat()
        b = buckets.setdefault(monday, {"week": monday, "days": [], "chips": 0,
                                        "echo": 0, "wb": weeks_before(d["day"])})
        b["chips"] += 1
        b["echo"] += 1 if d["echo"] else 0
        if d["day"] not in b["days"]:
            b["days"].append(d["day"])
    out = []
    for monday in sorted(buckets):
        b = buckets[monday]
        sub = [d for d in chips if d["day"] in b["days"]]
        row = dict(b)
        row["weekday"] = "Monday"
        for k in STRIP_DIMS + ["A"]:
            s = stat_block(sub, k)
            row[k] = {"n": s["n"], "median": s["median"]} if s else None
        out.append(row)
    return out


def rolling_block(chips, key):
    """Gap-aware trailing median over ``ROLL_SPAN_WEEKS`` (time keeps moving)."""
    echo = [d for d in chips if d.get(key) is not None]
    out = []
    for d in echo:
        day = day_date(d["day"])
        win = [e for e in echo
               if 0 <= (day - day_date(e["day"])).days <= ROLL_SPAN_WEEKS * 7]
        if len(win) >= ROLL_MIN_N:
            out.append({"day": d["day"], "x": x_of(d["day"]), "n": len(win),
                        "median": float(np.median([e[key] for e in win]))})
    return out


def block_arrangements(n, block_min=st.BLOCK_MIN):
    """The number of distinct moving-block arrangements of ``st.block_perm_p``.

    ``block_perm_p`` cuts ``max(block_min, n // 10)``-long blocks and permutes
    them, so a short series admits only ``nb!`` rearrangements and its
    permutation p can take only that many values. The report states it, because
    on 16 points it is 2 (0.5 or 1.0) and not a continuum.
    """
    block = max(block_min, n // 10)
    nb = max(2, n // block)
    if nb * block > n:
        block = max(1, n // nb)
    return {"block": int(block), "n_blocks": int(nb),
            "n_arrangements": int(math.factorial(nb))}


def trend_block(chips, key, x_key="t"):
    """Spearman + moving-block permutation of one channel against time.

    ``x`` is acquisition time **forward** (days since the first acquisition of
    the record), so a positive rho means "rises over the window".
    """
    xs = np.asarray([d[x_key] for d in chips if d.get(key) is not None],
                    dtype=float)
    ys = np.asarray([float(d[key]) for d in chips if d.get(key) is not None],
                    dtype=float)
    if ys.size < 2:
        return {"n": int(ys.size), "rho": None, "p": None, "perm_p": None}
    rho, p = st.spearman(xs, ys)
    return {"n": int(ys.size), "rho": rho, "p": p,
            "perm_p": st.block_perm_p(xs, ys, xs, rho),
            "blocks": block_arrangements(int(ys.size))}


def halves_block(chips):
    """The window split in two halves by acquisition index (coverage drift)."""
    days = sorted({d["day"] for d in chips})
    mid = (len(days) + 1) // 2
    out = []
    for label, sel in (("first half", days[:mid]), ("second half", days[mid:])):
        sub = [d for d in chips if d["day"] in sel]
        out.append({"half": label, "first_day": sel[0], "last_day": sel[-1],
                    "n_chips": len(sub),
                    "n_echo": sum(1 for d in sub if d["echo"]),
                    "echo_rate": (sum(1 for d in sub if d["echo"]) / len(sub)
                                  if sub else None)})
    return out


def month_marks(chips):
    """``[(month label, x of the first acquisition of that month)]``."""
    out = []
    for d in sorted_chips(chips):
        lab = MONTHS[day_date(d["day"]).month - 1]
        if not out or out[-1][0] != lab:
            out.append((lab, x_of(d["day"])))
    return out


def echo_free_run(chips):
    """The longest run of consecutive pre chips without an echo."""
    best = cur = 0
    for d in chips:
        cur = 0 if d["echo"] else cur + 1
        best = max(best, cur)
    return best


def paper_timeline_block(chips):
    """One entry per pre acquisition day: chips, echoes (the paper's own view)."""
    out = []
    for day in sorted({d["day"] for d in chips}):
        sub = [d for d in chips if d["day"] == day]
        out.append({"day": day, "n_chips": len(sub),
                    "n_with_echo": sum(1 for d in sub if d["echo"])})
    return out


# ---------------------------------------------------------------------------
# The recomputation
# ---------------------------------------------------------------------------
def compute(reference_paper=None):
    """Everything this companion shows, from the committed tables only.

    Inputs: the channel table ``data/ywf/ywf_ocv_channels.csv`` (33 unique chips),
    the committed mask cache ``ywf_windows_mask_cache.txt``, the generated
    ``ywf_ocv_channels_meta.json``, the committed paper result
    ``fig_ywf_ocv_paper.json`` and the committed analysis reference
    ``reference/ywf_ocv_findings.json``.
    """
    t0 = time.time()
    meta = core.load_csv_meta()
    rows = [row_view(r) for r in core.load_csv()]
    first_day = min(day_date(d["day"]) for d in rows)
    for d in rows:
        d["t"] = float((day_date(d["day"]) - first_day).days)

    chips = sorted_chips([d for d in rows if d["pre"]])
    post = sorted_chips([d for d in rows if not d["pre"]])
    echo_rows = [d for d in rows if d["echo"]]
    no_echo_rows = [d for d in rows if not d["echo"]]
    echo = [d for d in chips if d["echo"]]
    no_echo = [d for d in chips if not d["echo"]]
    days = sorted({d["day"] for d in chips})
    flag_rows = [dict(d, echo_flag=1.0 if d["echo"] else 0.0) for d in chips]
    pre_median = {k: median_of(echo, k) for k in STRIP_DIMS}
    normalised = {k: [{"day": d["day"], "x": x_of(d["day"]), "value": d[k],
                       "norm": (d[k] / pre_median[k]) if pre_median[k] else None}
                      for d in echo] for k in STRIP_DIMS}

    reference = core.load_findings()
    paper = core.load_reference(reference_paper or
                                os.path.join(core.DATA, "fig_ywf_ocv_paper.json"))
    cache = core.load_cache()
    rolling = {k: rolling_block(chips, k) for k in STRIP_DIMS}
    spread = max(abs(p["median"] / rolling[k][0]["median"] - 1.0)
                 for k in STRIP_DIMS for p in rolling[k])
    return {
        "rows": rows, "chips": chips, "post": post, "echo": echo,
        "echo_rows": echo_rows, "no_echo_rows": no_echo_rows,
        "no_echo": no_echo, "days": days, "flag_rows": flag_rows,
        "first_day": first_day.isoformat(),
        "last_pre_day": days[-1], "first_post_day": post[0]["day"],
        "n_post_days": len({d["day"] for d in post}),
        "event_weekday": EVENT.strftime("%A"),
        "n_pre_chips": len(chips), "n_pre_days": len(days),
        "n_pre_echo": len(echo), "n_pre_no_echo": len(no_echo),
        "n_post_chips": len(post),
        "n_post_echo": sum(1 for d in post if d["echo"]),
        "echo_rate_pre": len(echo) / len(chips),
        "echo_rate_post": sum(1 for d in post if d["echo"]) / len(post),
        "weeks": iso_week_block(chips), "gaps": gap_block(chips),
        "timeline": paper_timeline_block(chips),
        "halves": halves_block(chips),
        "months": month_marks(chips),
        "echo_free_run": echo_free_run(chips),
        "channels": {k: stat_block(echo, k) for k in FEATURES},
        "covariates": {k: stat_block(chips, k) for k in CONTROL_KEYS},
        "weights": {k: stat_block(echo, k) for k in CONTROL_KEYS},
        "trends": {k: trend_block(echo, k) for k in TREND_DIMS},
        "coverage_trend": trend_block(flag_rows, "echo_flag"),
        "covariate_trends": {k: trend_block(chips, k) for k in COVARIATE_KEYS},
        "rolling": rolling, "rolling_spread": spread,
        "normalised": normalised, "pre_median": pre_median,
        "meta": meta, "manifest": core.load_manifest(),
        "reference": reference, "paper": paper,
        "cache": cache,
        "coherence_unique": len({round(d["coherence"], 12)
                                 for d in rows if d["coherence"] is not None}),
        "p_reference_chips": max(r["P_reference_chips"] for r in cache),
        "n_paper_timeline": len(paper.get("timeline", [])),
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# Pins
# ---------------------------------------------------------------------------
# Exact doubles: a hand copy of a number can be off in the last bits, so the
# tolerance covers the values that come out of numpy/scipy. A deviation above
# the tolerance aborts the script.
TOL_KEYS = {"echo_rate", "rho", "p", "perm_p", "median", "mean", "std", "p5",
            "p25", "p75", "p95", "p_two_sided", "norm", "value"}


class Pins:
    """Collects pin comparisons (recomputed vs. committed) and reports them."""

    def __init__(self):
        self.checks = []
        self.fails = []

    def add(self, path, got, want, tol):
        try:
            exact = bool(got == want)
        except Exception:                       # pragma: no cover - numpy arrays
            exact = False
        dev = 0.0
        if not exact:
            numeric = (isinstance(got, (int, float))
                       and isinstance(want, (int, float))
                       and not isinstance(got, bool) and not isinstance(want, bool))
            dev = abs(float(got) - float(want)) if numeric else None
        ok = exact or (dev is not None and dev <= tol)
        self.checks.append({"path": path, "kind": "exact" if exact else "tolerance",
                            "abs_deviation": dev, "tolerance": tol, "ok": ok})
        if not ok:
            self.fails.append({"path": path, "recomputed": repr(got)[:200],
                               "committed": repr(want)[:200], "abs_deviation": dev})

    def cmp(self, prefix, got, want, key=None):
        """Recursive comparison: ``want`` dictates the structure and expectation."""
        tol = FLOAT_TOL if key in TOL_KEYS else 0.0
        if isinstance(want, dict):
            if not isinstance(got, dict):
                self.add(prefix, f"<{type(got).__name__}>", "<dict>", 0.0)
                return
            for k in want:
                if k not in got:
                    self.add(f"{prefix}.{k}", "<missing>", "<present>", 0.0)
                    continue
                self.cmp(f"{prefix}.{k}", got[k], want[k], k)
            return
        if isinstance(want, (list, tuple)):
            glen = len(got) if hasattr(got, "__len__") else -1
            if not isinstance(got, (list, tuple)) or glen != len(want):
                self.add(prefix, f"<len {glen}>", f"<len {len(want)}>", 0.0)
                return
            for i, (g, w) in enumerate(zip(got, want)):
                self.cmp(f"{prefix}[{i}]", g, w, key)
            return
        self.add(prefix, got, want, tol)

    def summary(self):
        devs = [c["abs_deviation"] for c in self.checks
                if c["abs_deviation"] not in (None, 0.0)]
        return {
            "n_checks": len(self.checks),
            "n_exact": sum(1 for c in self.checks if c["kind"] == "exact"),
            "n_within_tolerance": sum(1 for c in self.checks
                                      if c["kind"] == "tolerance"),
            "n_failed": len(self.fails),
            "max_abs_deviation": max(devs) if devs else 0.0,
            "float_tolerance": FLOAT_TOL,
            "failures": self.fails[:40],
        }


def coverage_numbers(res):
    """The coverage/power block, recomputed from the channel table alone."""
    n_pre, k_pre = res["n_pre_chips"], res["n_pre_echo"]
    n_post, k_post = res["n_post_chips"], res["n_post_echo"]
    p = core.fisher_p(k_pre, n_pre - k_pre, k_post, n_post - k_post)
    per_state = {
        "pre": {"n_chips": n_pre, "n_with_echo": k_pre,
                "n_without_echo": n_pre - k_pre, "echo_rate": k_pre / n_pre},
        "post": {"n_chips": n_post, "n_with_echo": k_post,
                 "n_without_echo": n_post - k_post, "echo_rate": k_post / n_post},
    }
    return {
        "per_state": per_state,
        "fisher": {"table_rows": ["pre", "post"], "table_cols": ["echo", "no_echo"],
                   "table": [[k_pre, n_pre - k_pre], [k_post, n_post - k_post]],
                   "p_two_sided": p},
        "n_unique_chips": n_pre + n_post,
    }


def cache_cross_check(res):
    """The committed mask cache against the channel table, channel by channel.

    ``masks.vector_from_row`` re-verifies the mask *rule* on every cache row
    (``A`` against ``|px|``, the peak against ``max_masked``, the floor against
    ``max(0.3 * peak, 5 * median)``) and derives ``A, D, F, S, gamma2`` from the
    committed coordinates; ``P`` comes from the cache's own reference set. All
    six are compared with the CSV.
    """
    by_id = {d["id"]: d for d in res["rows"]}
    rows, missing = [], 0
    for d in res["cache"]:
        r = by_id.get(d["mid"])
        if r is None:
            missing += 1
            continue
        v = masks.vector_from_row(d)
        rows.append({
            "mid": d["mid"], "date": d["date"], "orbit": d["orbit"],
            "A": (v["A"], r["A"]), "D": (v["D"], r["D"]), "F": (v["F"], r["F"]),
            "S": (v["S"], r["S"]), "gamma2": (v["gamma2"], r["gamma2"]),
            "P": (d["P"], r["P"]),
            "P_reference_chips": d["P_reference_chips"],
        })
    return {"n_cache_rows": len(res["cache"]), "n_compared": len(rows),
            "n_missing_in_csv": missing, "rows": rows,
            "P_reference_chips": sorted({r["P_reference_chips"] for r in rows})}


# ---------------------------------------------------------------------------
# The pins: four committed layers fix the same window
# ---------------------------------------------------------------------------
def committed_timeline(blob, state="pre"):
    """The entries of a committed ``timeline`` block, keyed by day."""
    return {e["day"]: e for e in blob["timeline"] if e.get("state") == state}


def pin_all(res):
    """Re-derive every number of this figure and check it against the commits."""
    meta, man = res["meta"], res["manifest"]
    paper, ref = res["paper"], res["reference"]
    cov = coverage_numbers(res)
    by_day = {d["day"]: d for d in res["timeline"]}
    P = Pins()

    # 1) the window, from four committed artifacts ---------------------------
    P.cmp("meta.n_rows", meta["n_rows"], len(res["rows"]))
    P.cmp("paper.n_rows", paper["n_rows"], len(res["rows"]))
    P.cmp("paper.n_unique_chips", paper["n_unique_chips"], meta["n_rows"])
    P.cmp("paper.n_echo_chips", paper["n_echo_chips"], len(res["echo_rows"]))
    P.cmp("meta.n_echo_masks", meta["n_echo_masks"], len(res["echo_rows"]))
    P.cmp("meta.n_rows_without_echo_mask", meta["n_rows_without_echo_mask"],
          len(res["no_echo_rows"]))
    P.cmp("meta.n_pre_collapse", meta["n_pre_collapse"], res["n_pre_chips"])
    P.cmp("meta.n_post_collapse", meta["n_post_collapse"], res["n_post_chips"])
    P.cmp("meta.n_P_set", meta["n_P_set"], len(res["echo_rows"]))
    P.cmp("meta.mask_window", meta["mask_window"], "7x7")
    P.cmp("paper.mask_window", paper["mask_window"], meta["mask_window"])
    P.cmp("ref.window_px", ref["window_px"], 49)
    for key in ("event", "site", "ocv", "features", "state_labels", "mask_rule"):
        P.cmp(f"ref.{key}", ref[key], paper[key])
    P.cmp("window.last_pre_day", res["last_pre_day"] < core.COLLAPSE_DATE, True)
    P.cmp("window.first_post_day", res["first_post_day"] > core.COLLAPSE_DATE, True)
    P.cmp("window.pre_states", sorted({d["state"] for d in res["chips"]}), ["pre"])
    P.cmp("window.post_states", sorted({d["state"] for d in res["post"]}), ["post"])
    P.cmp("window.event_weekday", res["event_weekday"], "Monday")

    # 2) the mask rule, the guard counts and the committed inputs -------------
    P.cmp("meta.mask_rule", meta["mask_rule"],
          {"peak_frac": masks.PEAK_FRAC, "median_mult": masks.MEDIAN_MULT,
           "min_n_masked": masks.MIN_N_MASKED, "median_kind": "np.median"})
    P.cmp("masks.MIN_N_MASKED", masks.MIN_N_MASKED,
          meta["mask_rule"]["min_n_masked"])
    P.cmp("mask.min_A_vs_min_n_masked",
          min(d["A"] for d in res["echo"]) >= masks.MIN_N_MASKED, True)
    P.cmp("meta.n_duplicate_windows", meta["n_duplicate_windows"],
          man["n_deduplicated_same_payload"])
    P.cmp("manifest.segment_resolved", man["segment_resolved"], False)
    for key, path in (("windows", core.WINDOWS_PATH),
                      ("mask_cache", core.CACHE_PATH),
                      ("mask_cache_manifest", core.MANIFEST_PATH),
                      ("measurements", core.MEAS_PATH),
                      ("segments", core.SEGMENTS_PATH)):
        P.cmp(f"meta.inputs.{key}.sha256", meta["inputs"][key]["sha256"],
              core.sha256(path))

    # 2b) the two resolved references (pass geometry, road anchor): their digests
    # are recorded in the analysis-layer reference, so this independent script
    # re-verifies them against the files on disk
    for key, path in (("burst_geometry", core.REF_BURSTS),
                      ("osm_road_extract", core.REF_OSM_ROAD),
                      ("osm_anchor", core.REF_GEOMETRY)):
        P.cmp(f"reference.sources.{key}.sha256", ref["sources"][key]["sha256"],
              core.sha256(path))

    # 3) the echo coverage: the CSV, the paper json and the reference ---------
    P.cmp("coverage.per_state", cov["per_state"], paper["coverage"]["per_state"])
    P.cmp("coverage.fisher.table", cov["fisher"]["table"],
          paper["coverage"]["fisher"]["table"])
    P.cmp("coverage.fisher.p_two_sided", cov["fisher"]["p_two_sided"],
          paper["coverage"]["fisher"]["p_two_sided"])
    P.cmp("coverage.n_windows", cov["n_unique_chips"],
          paper["coverage"]["n_windows"])
    P.cmp("coverage.window_px", paper["coverage"]["window_px"], ref["window_px"])
    P.cmp("ref.echo_coverage.per_state", ref["echo_coverage"]["per_state"],
          cov["per_state"])
    P.cmp("ref.echo_coverage.fisher", ref["echo_coverage"]["fisher"], cov["fisher"])
    P.cmp("ref.echo_coverage.n_unique_chips", ref["echo_coverage"]["n_unique_chips"],
          cov["n_unique_chips"])
    P.cmp("ref.echo_coverage.power", ref["echo_coverage"]["power"], paper["power"])
    P.cmp("coverage.sum_echo", sum(r[0] for r in cov["fisher"]["table"]),
          len(res["echo_rows"]))
    P.cmp("coverage.sum_chips", sum(sum(r) for r in cov["fisher"]["table"]),
          len(res["rows"]))

    # 4) A is not the pipeline's coherence_masked_pixels ----------------------
    P.cmp("export.n_compared", meta["n_coherence_masked_pixels_compared"],
          len(res["echo_rows"]))
    P.cmp("export.n_matched",
          sum(1 for d in res["echo_rows"]
              if d["coherence_masked_pixels"] is not None
              and abs(d["coherence_masked_pixels"] - d["A"]) < FLOAT_TOL),
          meta["n_coherence_masked_pixels_matched"])
    P.cmp("export.n_differing",
          sum(1 for d in res["echo_rows"]
              if d["coherence_masked_pixels"] is not None
              and abs(d["coherence_masked_pixels"] - d["A"]) >= FLOAT_TOL),
          meta["n_coherence_masked_pixels_differing"])

    # 5) the 26 pre days of the two committed timelines ------------------------
    ct, rt = committed_timeline(paper), committed_timeline(ref)
    P.cmp("paper.timeline.n_pre", len(ct), res["n_pre_days"])
    P.cmp("ref.timeline.n_pre", len(rt), res["n_pre_days"])
    P.cmp("paper.timeline.n_all", len(paper["timeline"]),
          res["n_pre_days"] + res["n_post_days"])
    P.cmp("ref.timeline.n_all", len(ref["timeline"]), len(paper["timeline"]))
    for day, e in sorted(ct.items()):
        P.cmp(f"paper.timeline[{day}].n_unique_chips", e["n_unique_chips"],
              by_day[day]["n_chips"])
        P.cmp(f"paper.timeline[{day}].n_with_echo", e["n_with_echo"],
              by_day[day]["n_with_echo"])
    for day, e in sorted(rt.items()):
        P.cmp(f"ref.timeline[{day}].n_unique_chips", e["n_unique_chips"],
              by_day[day]["n_chips"])
        P.cmp(f"ref.timeline[{day}].n_with_echo", e["n_with_echo"],
              by_day[day]["n_with_echo"])

    # 6) the cadence, the gaps and the week buckets ---------------------------
    committed_gaps = gap_block([{"day": d} for d in ct])
    P.cmp("gaps.n", len(res["gaps"]), len(committed_gaps))
    for i, g in enumerate(committed_gaps):
        P.cmp(f"gaps[{i}].after", res["gaps"][i]["after"], g["after"])
        P.cmp(f"gaps[{i}].before", res["gaps"][i]["before"], g["before"])
        P.cmp(f"gaps[{i}].days", res["gaps"][i]["days"], GAP_DAYS)
    P.cmp("gaps.steps_are_cadence_or_gap",
          sorted({d["days"] for d in res["gaps"]}), [GAP_DAYS])
    P.cmp("weeks.n", len(res["weeks"]), 23)
    P.cmp("weeks.sum_chips", sum(w["chips"] for w in res["weeks"]),
          res["n_pre_chips"])
    P.cmp("weeks.sum_echo", sum(w["echo"] for w in res["weeks"]), res["n_pre_echo"])
    P.cmp("weeks.monday_anchored",
          sorted({w["weekday"] for w in res["weeks"]}), ["Monday"])
    P.cmp("weeks.echo_per_day_max", max(w["echo"] for w in res["weeks"]), 1)
    P.cmp("halves.sum_chips", sum(h["n_chips"] for h in res["halves"]),
          res["n_pre_chips"])
    P.cmp("halves.sum_echo", sum(h["n_echo"] for h in res["halves"]),
          res["n_pre_echo"])
    P.cmp("halves.first_last", [res["halves"][0]["first_day"],
                                res["halves"][-1]["last_day"]],
          [res["days"][0], res["days"][-1]])
    P.cmp("echo_free_run", res["echo_free_run"], 3)

    # 7) the cache: the committed mask vector against the CSV -----------------
    cc = cache_cross_check(res)
    P.cmp("cache.n_cache_rows", cc["n_cache_rows"],
          meta["cache_cross_check"]["n_cache_rows"])
    P.cmp("cache.n_compared", cc["n_compared"], len(res["echo_rows"]))
    P.cmp("cache.n_missing_in_csv", cc["n_missing_in_csv"], 0)
    P.cmp("cache.P_reference_chips", cc["P_reference_chips"], [meta["n_P_set"]])
    for row in cc["rows"]:
        tag = f"cache[{row['date']}/{row['orbit']}]"
        for k in ("A", "D", "F", "S", "gamma2", "P"):
            P.add(f"{tag}.{k}", row[k][0], row[k][1], FLOAT_TOL)

    # 8) this figure's own relations ------------------------------------------
    P.cmp("window.n_pre_echo_plus_no_echo",
          res["n_pre_echo"] + res["n_pre_no_echo"], res["n_pre_chips"])
    P.cmp("window.A_defined_on_echo_rows_only",
          sum(1 for d in res["rows"] if d["A"] is not None),
          len(res["echo_rows"]))
    P.cmp("window.echo_flag_matches_echo",
          sorted(d["day"] for d in res["flag_rows"] if d["echo_flag"] == 1.0),
          sorted(d["day"] for d in res["echo"]))
    P.cmp("window.week_count_of_chips", len(res["chips"]), res["n_pre_chips"])
    for k in STRIP_DIMS:
        P.cmp(f"normalised[{k}].median",
              float(np.median([p["norm"] for p in res["normalised"][k]])), 1.0)
        P.cmp(f"normalised[{k}].n", len(res["normalised"][k]), res["n_pre_echo"])
        P.cmp(f"trends[{k}].n", res["trends"][k]["n"], res["n_pre_echo"])
        P.cmp(f"rolling[{k}].last_x", res["rolling"][k][-1]["x"],
              x_of(res["last_pre_day"]))
    P.cmp("trends[A].n", res["trends"]["A"]["n"], res["n_pre_echo"])
    P.cmp("coverage_trend.n", res["coverage_trend"]["n"], res["n_pre_chips"])
    P.cmp("echo_rate_pre", res["echo_rate_pre"],
          cov["per_state"]["pre"]["echo_rate"])
    P.cmp("echo_rate_post", res["echo_rate_post"],
          cov["per_state"]["post"]["echo_rate"])
    P.cmp("weeks.first_wb", res["chips"][0]["wb"], 25.714)
    P.cmp("weeks.last_wb", res["chips"][-1]["wb"], 0.857)
    P.cmp("months.labels", [m[0] for m in res["months"]],
          ["Aug", "Sep", "Oct", "Nov", "Dec", "Jan"])
    P.cmp("coherence.n_unique", res["coherence_unique"], 1)
    P.cmp("coherence.constant", round(res["covariates"]["coherence"]["median"], 9),
          0.998917749)
    P.cmp("P.reference_chips", res["p_reference_chips"], meta["n_P_set"])
    for k in [c for c in core.CONTROLS if c in res["covariates"]]:
        P.cmp(f"controls_by_state.pre[{k}]", res["covariates"][k],
              paper["controls_by_state"]["pre"][k])
    return P


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def axis_ticks(res):
    """``[(x, label)]`` — the month starts of the window plus the collapse.

    The label carries both coordinates at once: the calendar day on the first
    line, the weeks before the collapse of that day on the second. That keeps
    the shared x axis legible at ~1.5 in per panel, where a numeric tick and a
    separate row of month names would collide.
    """
    d0, d1 = day_date(res["first_day"]), day_date(res["last_pre_day"])
    out = []
    month = dt.date(d0.year, d0.month, 1)
    while month <= d1:
        if month > d0:
            x = x_of(month.isoformat())
            out.append((x, f"{MONTHS[month.month - 1]} {month.day}\n"
                           f"{abs(x):.0f}"))
        month = (month + dt.timedelta(days=32)).replace(day=1)
    out.append((0.0, f"{MONTHS[EVENT.month - 1]} {EVENT.day}\n0"))
    return out


def decorate(ax, res, xlabel=False):
    """The shared x axis: calendar grid, the four cadence gaps, the collapse."""
    d0, d1 = day_date(res["first_day"]), day_date(res["last_pre_day"])
    month = dt.date(d0.year, d0.month, 1)
    while month <= d1:
        if month > d0:
            ax.axvline(x_of(month.isoformat()), color=MONTH_COLOR, lw=0.6,
                       ls=(0, (1, 2)), zorder=0)
        month = (month + dt.timedelta(days=32)).replace(day=1)
    for g in res["gaps"]:
        ax.axvspan(g["x"][0], g["x"][1], color=GAP_COLOR, lw=0, zorder=0)
    ax.axvline(0.0, color=EVENT_COLOR, lw=1.0, zorder=1.5)
    ax.set_xlim(-26.6, 0.95)
    ticks = axis_ticks(res)
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks])
    ax.tick_params(axis="x", labelsize=5.8, pad=1.6)
    ax.set_xlabel("weeks before the collapse" if xlabel else "", fontsize=8.0)


def floor_ticks(ax, res, label=None):
    """The chips *without* an echo, as small bars on the floor of the panel."""
    t = ax.get_xaxis_transform()
    xs = [x_of(d["day"]) + orbit_dx(d["orbit"]) for d in res["no_echo"]]
    ax.plot(xs, [0.022] * len(xs), marker="|", ls="none", ms=6.0, mew=0.9,
            color=NO_ECHO_COLOR, transform=t, clip_on=False, label=label)


def median_text(ax, res, key):
    """The pre-window median of one channel: a dashed line, labelled on it."""
    m = res["pre_median"][key]
    ax.axhline(m, color="#555555", lw=0.7, ls=(0, (4, 2)), zorder=1)
    ax.text(-26.4, m, f"pre median {m:.3f}", ha="left", va="bottom",
            fontsize=5.6, color="#555555", zorder=5,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                      pad=0.6))


def trend_text(ax, block, label="rho", fontsize=5.8):
    """The trend of one series, right-aligned in the panel's top corner."""
    if not block or block.get("rho") is None:
        return
    lines = [f"{label} {block['rho']:+.3f}, p {block['p']:.2f}"]
    if block.get("perm_p") is not None:
        lines.append(f"block-perm p {block['perm_p']:.2f}")
    ax.text(0.985, 0.965, "\n".join(lines), transform=ax.transAxes,
            ha="right", va="top", fontsize=fontsize, color="#444444")


def fig_a(ax, res):
    """The echo of every acquisition, with the chip intensity beside it."""
    xs = [x_of(d["day"]) + orbit_dx(d["orbit"]) for d in res["chips"]]
    ys = [(d["A"] or 0.0) for d in res["chips"]]
    cols = [ECHO_MODE_COLORS.get(d["echo_mode"], NO_ECHO_COLOR) if d["echo"]
            else NO_ECHO_COLOR for d in res["chips"]]
    for x, y, c in zip(xs, ys, cols):
        ax.plot([x, x], [0.0, y], "-", color=c, lw=0.9, zorder=2)
        ax.plot([x], [y], "o", ms=2.5, color=c, markeredgewidth=0, zorder=3)
    handles = [plt.Line2D([], [], marker="o", ls="none", ms=3.2,
                          color=ECHO_MODE_COLORS[m],
                          label=f"{m} (A {cut})")
               for m, cut in (("compact", "2"), ("intermediate", "3-5"),
                              ("distributed", "6+"))]
    handles.append(plt.Line2D([], [], marker="o", ls="none", ms=3.2,
                              color=NO_ECHO_COLOR, label="no echo (A = 0)"))
    ax.set_ylim(-0.5, 8.4)
    ax.set_yticks([0, 2, 4, 6])
    ax.set_ylabel(core.CHANNEL_AXIS["A"])
    ax.legend(handles=handles, frameon=False, ncols=1, loc="upper left",
              fontsize=5.5, handletextpad=0.4, borderpad=0.1,
              labelspacing=0.22, handlelength=1.0)
    tw = ax.twinx()
    inten = [d["intensity"] for d in res["chips"]]
    tw.plot(xs, inten, "-", color=COVARIATE_COLORS["intensity"], lw=0.8,
            alpha=0.85, zorder=1)
    tw.set_yscale("log")
    tw.set_ylim(250.0, 1.2e5)
    tw.set_ylabel(core.CHANNEL_AXIS["intensity"] + " (log)", fontsize=7.5,
                  color=COVARIATE_COLORS["intensity"])
    tw.tick_params(axis="y", labelsize=6.2,
                   colors=COVARIATE_COLORS["intensity"])
    tw.spines["top"].set_visible(False)
    tw.grid(False)
    widest = max(res["gaps"], key=lambda g: g["x"][1] - g["x"][0])
    ax.text(float(np.mean(widest["x"])), 0.60, "12-day gap", rotation=90,
            transform=ax.get_xaxis_transform(), ha="center", va="center",
            fontsize=5.8, color="#8C8C8C")
    ax.text(0.0, 0.60, "collapse", rotation=90, transform=ax.get_xaxis_transform(),
            ha="right", va="center", fontsize=6.0, color=EVENT_COLOR)
    panel(ax, "a", "the echo of every acquisition")


def fig_strip(ax, res, key, letter):
    """One OCV channel: a dot per echo-bearing chip + the rolling median."""
    xs = [x_of(d["day"]) + orbit_dx(d["orbit"]) for d in res["echo"]]
    ys = [d[key] for d in res["echo"]]
    col = CHANNEL_COLORS[key]
    ax.plot(xs, ys, "o", ms=3.0, color=col, alpha=0.85, markeredgewidth=0,
            zorder=3, label="one chip")
    roll = res["rolling"][key]
    ax.plot([p["x"] for p in roll], [p["median"] for p in roll], "-",
            color="#333333", lw=1.1, zorder=4,
            label=f"rolling median ({ROLL_SPAN_WEEKS:.0f} weeks)")
    floor_ticks(ax, res, label="no echo")
    median_text(ax, res, key)
    trend_text(ax, res["trends"][key])
    lo = min(min(ys), min(p["median"] for p in roll))
    hi = max(max(ys), max(p["median"] for p in roll))
    ax.set_ylim(lo - 0.18 * (hi - lo), hi + 0.32 * (hi - lo))
    ax.set_ylabel(core.CHANNEL_AXIS[key])
    if key == "P":
        ax.legend(frameon=False, loc="upper left", fontsize=5.4,
                  handletextpad=0.4, borderpad=0.1, labelspacing=0.22,
                  markerscale=1.2)
    panel(ax, letter, f"{key}, {res['n_pre_echo']} echo chips")
    if key == "D":
        ax.set_ylim(bottom=0.0)


def fig_e(ax, res):
    """The three OCV channels normalised to their pre-window median, one log axis."""
    for k in STRIP_DIMS:
        pts = res["normalised"][k]
        ax.plot([p["x"] for p in pts], [p["norm"] for p in pts], "o", ms=2.8,
                color=CHANNEL_COLORS[k], alpha=0.75, markeredgewidth=0,
                label=f"{k} (median {res['pre_median'][k]:.3f})")
        roll = res["rolling"][k]
        ax.plot([p["x"] for p in roll], [p["median"] / res["pre_median"][k]
                                         for p in roll], "-",
                color=CHANNEL_COLORS[k], lw=0.9, alpha=0.9)
    ax.axhline(1.0, color="#555555", lw=0.7, ls=(0, (4, 2)), zorder=1)
    ax.set_yscale("log")
    ax.set_ylabel("value / pre-window median  (log)")
    ax.legend(frameon=False, loc="upper left", fontsize=5.8, ncols=1,
              handletextpad=0.3, borderpad=0.1, labelspacing=0.25)
    floor_ticks(ax, res)
    ax.text(0.985, 0.035, "median = 1 by construction",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8,
            color="#444444")
    panel(ax, "e", "the three OCV channels")


def fig_f(ax, res):
    """What else moves over these weeks: temperature, wind speed and gust."""
    xs = [x_of(d["day"]) + orbit_dx(d["orbit"]) for d in res["chips"]]
    temp = [d[TEMP_KEY] for d in res["chips"]]
    ax.plot(xs, temp, "o-", ms=3.0, lw=0.8, color=COVARIATE_COLORS[TEMP_KEY],
            alpha=0.9, markeredgewidth=0, label="temperature (C)")
    ax.set_ylabel(core.CHANNEL_AXIS[TEMP_KEY], color=COVARIATE_COLORS[TEMP_KEY])
    ax.tick_params(axis="y", colors=COVARIATE_COLORS[TEMP_KEY])
    ax.set_ylim(min(temp) - 4.0, max(temp) + 9.0)
    tw = ax.twinx()
    for k in WIND_KEYS:
        tw.plot(xs, [d[k] for d in res["chips"]], "o-", ms=2.4, lw=0.7,
                color=COVARIATE_COLORS[k], alpha=0.85, markeredgewidth=0,
                label=core.CHANNEL_TITLE[k])
    tw.set_ylim(0.0, max(max(d[k] for d in res["chips"]) for k in WIND_KEYS) * 1.25)
    tw.set_ylabel("wind / gust (m/s)")
    tw.tick_params(axis="y", labelsize=6.2)
    tw.spines["top"].set_visible(False)
    tw.grid(False)
    floor_ticks(ax, res)
    ct = res["covariate_trends"]
    ax.text(0.985, 0.965, "\n".join([
        f"temperature rho {ct[TEMP_KEY]['rho']:+.3f}, p {ct[TEMP_KEY]['p']:.2f}",
        f"wind rho {ct['wind_speed_ms']['rho']:+.3f}, "
        f"gust rho {ct['wind_gust_ms']['rho']:+.3f}",
        f"block-perm p {ct[TEMP_KEY]['perm_p']:.2f} (temperature)"]),
        transform=ax.transAxes, ha="right", va="top", fontsize=5.8,
        color="#444444")
    panel(ax, "f", f"what else moves, {res['n_pre_chips']} chips")


def fig_weeks(res, out_dir=None):
    """The 2x3 grid: coverage, the three OCV channels, and the co-variates."""
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 0.85 * COL_WIDTH))
    bottom = axes[1]
    for row in axes:
        for ax in row:
            decorate(ax, res, xlabel=ax in bottom)
    fig_a(axes[0][0], res)
    fig_strip(axes[0][1], res, "gamma2", "b")
    fig_strip(axes[0][2], res, "P", "c")
    fig_strip(axes[1][0], res, "D", "d")
    fig_e(axes[1][1], res)
    fig_f(axes[1][2], res)
    fig.suptitle(
        f"(YWF) the OCV time series in the weeks before the collapse of "
        f"Unit 21 ({core.COLLAPSE_DATE})\n"
        f"{res['n_pre_chips']} chips on {res['n_pre_days']} acquisition days "
        f"({res['n_pre_echo']} with an echo), {res['first_day']} -> "
        f"{res['last_pre_day']}", fontsize=8.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    return save_fig(fig, "fig_ywf_ocv_weeks", out_dir)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def fmt(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}" if math.isfinite(x) else "n/a"
    return str(x)


def md_table(header, rows, nd=3):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(fmt(c, nd) for c in r) + " |")
    return out


def acquisition_rows(res):
    """The 27 chips of the window, in acquisition order."""
    out = []
    for d in res["chips"]:
        out.append([d["day"], d["wb"], (d["orbit"] or "-")[:4],
                    "yes" if d["echo"] else "no", d["echo_mode"] or "-",
                    d["A"], d["F"], d["S"], d["gamma2"], d["P"], d["D"],
                    d["intensity"], d["temperature_c"], d["wind_speed_ms"],
                    d["wind_gust_ms"], d["status"]])
    return out


def week_rows(res):
    out = []
    for w in res["weeks"]:
        out.append([w["week"], ", ".join(x[5:] for x in w["days"]), w["wb"],
                    w["chips"], w["echo"],
                    (w["A"] or {}).get("median"), (w["gamma2"] or {}).get("median"),
                    (w["P"] or {}).get("median"), (w["D"] or {}).get("median")])
    return out


def report(res, pins):
    meta, man = res["meta"], res["manifest"]
    cov = coverage_numbers(res)
    ch, ct = res["channels"], res["covariate_trends"]
    md = []
    md.append(f"# {core.SITE} — the OCV time series in the weeks before the collapse")
    md.append("")
    md.append(f"**Site** {core.SITE} · **asset** {core.ASSET_NAME} · "
              f"**event** {core.COLLAPSE_DATE} ({res['event_weekday']}) · "
              f"**window** {res['first_day']} -> {res['last_pre_day']} "
              f"({weeks_before(res['first_day']):.1f} -> "
              f"{weeks_before(res['last_pre_day']):.1f} weeks before the collapse) · "
              f"**chips** {res['n_pre_chips']} on {res['n_pre_days']} acquisition "
              f"days, {res['n_pre_echo']} with an echo")
    md.append("")
    md.append("```")
    md.append("OCV  = [gamma2, P, D]              observability core vector")
    md.append("x    = -(2026-02-02 - day) / 7      weeks before the collapse (0 = event)")
    md.append("tick = the 1st of every month, its weeks before the collapse under it")
    md.append("mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2")
    md.append("step = 6 d (the two orbits alternate), 12 d = a gap: drawn, never "
              "interpolated over")
    md.append("post event chips: none in this figure — the window is the subject")
    md.append("```")
    md.append("")
    md.append("## The window")
    md.append("")
    md.append("| fact | value |")
    md.append("| --- | --- |")
    md.append(f"| acquisition days | **{res['n_pre_days']}**, "
              f"{res['first_day']} -> {res['last_pre_day']}, "
              f"a {CADENCE_DAYS}-day step with "
              f"{len(res['gaps'])} gaps of {GAP_DAYS} days |")
    md.append(f"| chips | **{res['n_pre_chips']}** unique 7x7 chips "
              f"({meta['n_windows']} committed windows, {len(res['rows'])} unique "
              f"chips in the whole record) |")
    md.append(f"| echo coverage | **{res['n_pre_echo']} of {res['n_pre_chips']}** = "
              f"{res['echo_rate_pre'] * 100:.1f} % (post-event, for reference: "
              f"{res['n_post_echo']} of {res['n_post_chips']} = "
              f"{res['echo_rate_post'] * 100:.1f} %, Fisher exact p = "
              f"{cov['fisher']['p_two_sided']:.4f}) |")
    md.append(f"| longest echo-free run | **{res['echo_free_run']}** consecutive "
              f"chips without an echo |")
    md.append(f"| ISO weeks | **{len(res['weeks'])}** Monday-anchored buckets, "
              f"{len([w for w in res['weeks'] if w['echo']])} of them hold an echo |")
    md.append(f"| weather | temperature median "
              f"{fmt(res['covariates']['temperature_c']['median'], 1)} C "
              f"(p5..p95 {fmt(res['covariates']['temperature_c']['p5'], 1)}.."
              f"{fmt(res['covariates']['temperature_c']['p95'], 1)} C), wind "
              f"{fmt(res['covariates']['wind_speed_ms']['median'], 2)} m/s, gust "
              f"{fmt(res['covariates']['wind_gust_ms']['median'], 2)} m/s |")
    md.append(f"| the pipeline `coherence` column | "
              f"{fmt(res['covariates']['coherence']['median'], 12)} on all "
              f"{len(res['rows'])} rows (`n_unique` = {res['coherence_unique']}) — "
              f"a constant, so it cannot be a co-variate of the window |")
    md.append(f"| the `P` reference | the **{res['p_reference_chips']}** "
              f"echo-bearing chips of the whole record; `P` is a window-wide "
              f"reference, not a per-chip quantity |")
    md.append("")
    md.append("## The channels over the weeks")
    md.append("")
    md.append("Dots are one chip each (no pooling); the rolling median is the "
              "trailing 4-week median of the echo-bearing chips with the gap "
              "honoured — a chip after a 12-day step sees the chips of the "
              "previous weeks, never an interpolated point.")
    md.append("")
    rows = []
    for k in STRIP_DIMS + ["A"]:
        s = ch[k]
        roll = res["rolling"].get(k)
        tr = res["trends"][k]
        rows.append([k, s["n"], s["median"], s["p5"], s["p95"],
                     (roll[0]["median"] if roll else None),
                     (roll[-1]["median"] if roll else None),
                     tr["rho"], tr["p"], tr["perm_p"]])
    md += md_table(["channel", "n", "median", "p5", "p95", "rolling first",
                    "rolling last", "rho vs. time", "p", "perm p"], rows)
    md.append("")
    md.append(f"`A` (the echo of the window) is the coverage flag: "
              f"{res['n_pre_echo']} of {res['n_pre_chips']} chips carry one "
              f"(fig. a), and the three OCV channels are defined exactly there. "
              f"The rolling medians of `gamma2`, `P` and `D` stay inside "
              f"±{100 * res['rolling_spread']:.0f} % of their first value over "
              f"the whole window (figs. b-d), while the coverage moves with the "
              f"halves of the window and the temperature falls (figs. a, f).")
    md.append("")
    md.append("## The trend tests, and what they can carry")
    md.append("")
    md.append("Spearman rho against **acquisition time forward** (days since "
              "2025-08-06), with the moving-block permutation p of "
              "`ywf_ocv_stats.block_perm_p` beside it.")
    md.append("")
    rows = []
    for k in TREND_DIMS:
        t = res["trends"][k]
        rows.append([f"{k} (echo chips)", t["n"], t["rho"], t["p"], t["perm_p"]])
    t = res["coverage_trend"]
    rows.append(["echo flag (all chips)", t["n"], t["rho"], t["p"], t["perm_p"]])
    for k in COVARIATE_KEYS:
        t = ct[k]
        rows.append([core.CHANNEL_TITLE.get(k, k) + " (all chips)", t["n"],
                     t["rho"], t["p"], t["perm_p"]])
    md += md_table(["series", "n", "rho", "p", "perm p"], rows)
    md.append("")
    tb = res["coverage_trend"]["blocks"]
    md.append(f"**The permutation cannot resolve this window.** "
              f"`block_perm_p` cuts blocks of "
              f"`max({st.BLOCK_MIN}, n // 10)` along acquisition order, so the "
              f"{res['coverage_trend']['n']}-point series of the trend table is "
              f"cut into {tb['n_blocks']} blocks of {tb['block']} chips and the "
              f"null has only {tb['n_arrangements']} distinct arrangements: its p "
              f"can only take {tb['n_arrangements']} values, and the number in "
              f"that column is one of them. The Spearman p beside it is the "
              f"informative one — and on this window what it reports is the "
              f"*absence* of a channel trend.")
    md.append("")
    md.append("## The coverage drift, the halves and the gaps")
    md.append("")
    md += md_table(["half", "days", "chips", "with echo", "echo rate"],
                   [[h["half"], f"{h['first_day']} -> {h['last_day']}",
                     h["n_chips"], h["n_echo"], h["echo_rate"]]
                    for h in res["halves"]])
    md.append("")
    md.append("The two halves are the window cut at its median acquisition index: "
              "the echo is not equally present along the axis, so a *pooled* "
              "pre-event median is a median of the chips the weather happened to "
              "make observable, and the figure says so instead of smoothing it.")
    md.append("")
    md += md_table(["after", "before", "days", "weeks before (after -> before)"],
                   [[g["after"], g["before"], g["days"],
                     f"{g['wb'][0]:.1f} -> {g['wb'][1]:.1f}"]
                    for g in res["gaps"]])
    md.append("")
    md.append(f"Every other step of the window is the {CADENCE_DAYS}-day repeat "
              f"step of the two alternating orbits (ascending and descending "
              f"passes); the four gaps above are drawn as bands in all six panels "
              f"and are never interpolated across.")
    md.append("")
    md.append("## The acquisition table")
    md.append("")
    md += md_table(["day", "W-", "orbit", "echo", "mode", "A", "F", "S",
                    "gamma2", "P", "D", "intensity", "T (C)", "wind", "gust",
                    "status"], acquisition_rows(res))
    md.append("")
    md.append("`A`, `F`, `S`, `gamma2`, `P` and `D` are empty where the chip "
              "carries no echo: they are undefined there, not zero.")
    md.append("")
    md.append(f"## The ISO-week buckets ({len(res['weeks'])} weeks)")
    md.append("")
    md += md_table(["week (Mon)", "days", "W-", "chips", "echo", "A med",
                    "gamma2 med", "P med", "D med"], week_rows(res))
    md.append("")
    md.append("## Reading (a statement about coverage, not about damage)")
    md.append("")
    md.append(f"Inside this window the echo is present on "
              f"**{res['n_pre_echo']} of {res['n_pre_chips']}** chips and the "
              f"longest run without one is **{res['echo_free_run']}** chips; the "
              f"OCV values themselves show no monotone drift "
              f"(gamma2 rho {res['trends']['gamma2']['rho']:+.3f}, p "
              f"{res['trends']['gamma2']['p']:.2f}; P "
              f"{res['trends']['P']['rho']:+.3f}, p {res['trends']['P']['p']:.2f}; "
              f"D {res['trends']['D']['rho']:+.3f}, p "
              f"{res['trends']['D']['p']:.2f}), while the temperature falls from "
              f"late summer to deep winter (rho "
              f"{ct['temperature_c']['rho']:+.3f}, p "
              f"{ct['temperature_c']['p']:.4f}) and the coverage moves between the "
              f"halves. The pooled pre/post contrast of the paper figure rests on "
              f"a window whose coordinate, coverage and weather are all moving: "
              f"that is the reason this companion prints the time axis instead of "
              f"only its median.")
    md.append("")
    md.append(f"Consequence for the headline of "
              f"[`fig_ywf_ocv_paper.md`](fig_ywf_ocv_paper.md): the "
              f"{res['n_pre_echo']}/{res['n_pre_chips']} against "
              f"{res['n_post_echo']}/{res['n_post_chips']} coverage table "
              f"(Fisher p = {cov['fisher']['p_two_sided']:.4f}) must be read with "
              f"the drift above attached — the pre side is not a stationary "
              f"baseline, and the post side holds six chips. The figure claims "
              f"neither a trend before the event nor a change at it.")
    md.append("")
    md.append("## What this figure does not claim")
    md.append("")
    md.append(f"* it draws **no** post-event chip: it is not a contrast, and it "
              f"adds no number to the A-E pin contract of "
              f"`fig_ywf_ocv_paper.py` ({res['n_pre_chips']} of the record's "
              f"{len(res['rows'])} chips are drawn, the "
              f"{res['n_post_chips']} post-event ones are not);")
    md.append(f"* `A` here is the recomputed echo mask of the committed 7x7 "
              f"window — **not** the pipeline's `coherence_masked_pixels`, which "
              f"coincides with `A` on "
              f"{meta['n_coherence_masked_pixels_matched']} of "
              f"{meta['n_coherence_masked_pixels_compared']} echo-bearing rows;")
    md.append(f"* `P` is computed over the "
              f"{res['p_reference_chips']} echo-bearing chips of the whole "
              f"record, so the chips of this window are not independent draws of "
              f"it and its rolling median must be read as a shape statistic, not "
              f"as a per-chip measurement;")
    md.append("* the moving-block permutation of this window cannot resolve "
              "p-values finer than the number of block arrangements, so no "
              "significance claim is made from it.")
    md.append("")
    md.append("## Cross-checks")
    md.append("")
    ps = pins.summary()
    md.append(f"Each run re-derives the window from the committed channel table "
              f"and compares it, field by field, with the committed mask cache "
              f"(all {len(res['cache'])} echo rows, the mask rule re-verified on "
              f"each of them), the generated metadata, the paper result JSON and "
              f"the analysis reference. **{ps['n_checks']} checks "
              f"({ps['n_exact']} exact, {ps['n_within_tolerance']} within "
              f"{ps['float_tolerance']:g}), {ps['n_failed']} deviations.**")
    md.append("")
    md.append("| committed layer | what is pinned |")
    md.append("| --- | --- |")
    md.append(f"| `data/ywf/{os.path.basename(core.CACHE_PATH)}` | the mask "
              f"vector `A, D, F, S, gamma2` re-derived from the committed "
              f"coordinates and `P` from the cache's own reference set, against "
              f"the CSV — {len(res['cache'])} rows x 6 channels, the mask rule "
              f"re-verified row by row |")
    md.append("| `data/ywf/ywf_ocv_channels_meta.json` | the guard counts "
              "(`n_rows`, `n_unique_chips`, `n_echo_masks`, `n_P_set`, "
              "`n_pre_collapse`, `n_post_collapse`), the mask rule, the "
              "`coherence_masked_pixels` cross-check and the sha256 of all five "
              "committed inputs |")
    md.append("| `data/ywf/fig_ywf_ocv_paper.json` | the record counts, the "
              "coverage table and its Fisher test, the per-day timeline of the "
              f"{res['n_pre_days']} pre-event days, the pre-side co-variate "
              "summary, the mask rule |")
    md.append(f"| `data/ywf/reference/{os.path.basename(core.REF_FINDINGS)}` | "
              "the same coverage table, its power caveat and the per-day "
              "timeline — the two committed timelines must agree day by day |")
    md.append("| `data/ywf/reference/ywf_bursts.json` | the CDSE pass geometry of "
              "the record's chips (sub-swath, relative orbit, polarisation, "
              "azimuth time) and `ywf_osm_road.json` / `ywf_geometry.json` (the "
              "road the tower fell on, the estimated grid) — by the sha256 "
              "digests the reference records |")
    md.append("")
    if ps["n_failed"]:
        md.append("| failed check | recomputed | committed |")
        md.append("| --- | --- | --- |")
        for f in ps["failures"][:40]:
            md.append(f"| `{f['path']}` | `{f['recomputed']}` | "
                      f"`{f['committed']}` |")
        md.append("")
    md.append("```bash")
    md.append("python3 code/ywf/ywf_ocv_masks.py --check              # the mask rule")
    md.append("python3 code/ywf/ywf_ocv_channels_csv.py             # build + verify")
    md.append("bash figures/ywf/fig_ywf_ocv_weeks.sh                 # this figure + pins")
    md.append("bash figures/ywf/fig_ywf_ocv_paper.sh                 # the pooled A-E figure")
    md.append("```")
    return "\n".join(md) + "\n"


# ---------------------------------------------------------------------------
# JSON artifact and the entry point
# ---------------------------------------------------------------------------
def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def json_safe(o):
    """Non-finite floats become ``None`` (JSON has no NaN literal)."""
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, np.floating):
        return json_safe(float(o))
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def build_json(res, pins_summary):
    """The result JSON — every recomputed number, without the per-chip rows."""
    meta = res["meta"]
    return {
        "site": core.SITE, "site_short": core.SITE_SHORT,
        "asset_id": core.ASSET_ID, "asset_name": core.ASSET_NAME,
        "event": core.COLLAPSE_DATE, "event_weekday": res["event_weekday"],
        "generated_by": "code/ywf/fig_ywf_ocv_weeks.py",
        "ocv": list(core.OCV), "features": list(core.FEATURES),
        "mask_window": "7x7", "mask_rule": meta["mask_rule"],
        "axis": {"name": "weeks_before_collapse",
                 "definition": "-(event - day).days / 7, 0 at the event",
                 "ticks": [{"label": t[1].replace("\n", " "),
                            "weeks_before": round(abs(t[0]), 3)}
                           for t in axis_ticks(res)]},
        "cadence_days": CADENCE_DAYS, "gap_days": GAP_DAYS,
        "rolling_span_weeks": ROLL_SPAN_WEEKS, "rolling_min_n": ROLL_MIN_N,
        "n_rows": len(res["rows"]),
        "n_pre_chips": res["n_pre_chips"], "n_pre_days": res["n_pre_days"],
        "n_pre_echo": res["n_pre_echo"], "n_pre_no_echo": res["n_pre_no_echo"],
        "n_post_chips": res["n_post_chips"], "n_post_echo": res["n_post_echo"],
        "echo_rate_pre": res["echo_rate_pre"],
        "echo_rate_post": res["echo_rate_post"],
        "n_iso_weeks": len(res["weeks"]),
        "n_echo_weeks": len([w for w in res["weeks"] if w["echo"]]),
        "first_acquisition": res["first_day"],
        "last_pre_acquisition": res["last_pre_day"],
        "first_post_acquisition": res["first_post_day"],
        "longest_echo_free_run": res["echo_free_run"],
        "coherence_unique_values": res["coherence_unique"],
        "p_reference_chips": res["p_reference_chips"],
        "pre_median": res["pre_median"], "rolling_spread": res["rolling_spread"],
        "channels": res["channels"], "covariates": res["covariates"],
        "trends": res["trends"], "coverage_trend": res["coverage_trend"],
        "covariate_trends": res["covariate_trends"],
        "rolling": res["rolling"], "normalised": res["normalised"],
        "acquisitions": [{k: d[k] for k in
                          ("day", "wb", "orbit", "status", "traffic", "echo",
                           "echo_mode", "A", "F", "S", "gamma2", "P", "D",
                           "intensity", "temperature_c", "wind_speed_ms",
                           "wind_gust_ms", "coherence_masked_pixels")}
                         for d in res["chips"]],
        "iso_weeks": res["weeks"], "gaps": res["gaps"],
        "timeline": res["timeline"], "halves": res["halves"],
        "month_marks": res["months"],
        "coverage": coverage_numbers(res), "power": res["paper"]["power"],
        "pins": pins_summary,
        # no wall-clock time: with the fixed seeds two runs must produce the
        # identical JSON, so the runtime stays on stdout only.
        "figures": ["figures/ywf/fig_ywf_ocv_weeks.png"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    ap.add_argument("--no-figures", action="store_true",
                    help="run the cross-checks and write no PNG")
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, the window before "
          f"{core.COLLAPSE_DATE} …")
    res = compute()
    print(f"  {res['n_pre_chips']} pre-event chips on {res['n_pre_days']} "
          f"acquisition days ({res['n_pre_echo']} carry an echo), "
          f"{len(res['weeks'])} ISO weeks, {res['runtime_s']:.1f} s")

    pins = pin_all(res)
    summary = pins.summary()
    print(f"pins: {summary['n_checks']} checks, {summary['n_exact']} exact, "
          f"{summary['n_within_tolerance']} within "
          f"{summary['float_tolerance']:g}, {summary['n_failed']} failed")

    os.makedirs(args.figdir, exist_ok=True)
    if not args.no_figures:
        fig_weeks(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(json_safe(build_json(res, summary)), fh, indent=1,
                  sort_keys=True, default=json_default)
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res, pins))
    print(f"  {os.path.relpath(args.md)}")

    if summary["n_failed"]:
        for f in summary["failures"][:20]:
            print("PIN FAILED:", f["path"], "recomputed=", f["recomputed"],
                  "committed=", f["committed"])
        raise SystemExit(f"PIN VERIFICATION FAILED ({summary['n_failed']})")
    print(f"PIN VERIFICATION OK  ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
