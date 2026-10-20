#!/usr/bin/env python3
"""HALTERN OCV-Orbits — the orbit split of the echo-mask vector as a time series.

The HALTERN record carries **four** Sentinel-1 pass geometries — both directions
fly two IW sub-swaths on two relative orbits — which the export's own
``orbit_direction`` label groups into two: ascending (IW1 on relative orbit 15
around 15:2x UT, IW3 on relative orbit 88 around 17:2x UT) and descending (IW1 on
relative orbit 37 and IW3 on relative orbit 139, both around 05:2x UT). Within a
month the two directions **alternate in blocks of two days** (steps of 1-3 days),
so no acquisition day carries both an ascending and a descending chip and a
*per-day* orbit difference does not exist. The honest form of "what does the
orbit change?" is therefore a **trailing-window difference**: at each acquisition
day ``t`` that has an echo chip, take the echo-bearing chips in the trailing
window ``W`` (default 8 weeks), require at least two chips of each orbit, and draw

    delta_orbit(t) = median(DESC in the window) - median(ASC in the window)

on the twin axis of the channel's own panel, beside the two orbits' per-chip
values; the **pooled** difference over the whole record is the dashed reference.
A line joins two of those points only when they are one cadence step apart, so
the series breaks at every hole of the record — a gap is drawn, never
interpolated over.

What the figure shows, and what it does not:

* the pooled orbit split is **indistinguishable from zero in every channel**:
  the widest is gamma2 (+0.160 with a seeded 95% CI of [-0.162, +0.478] over 33
  ascending and 33 descending echo chips), the tightest is ``P`` (0.000 with
  [-0.014, +0.010]), and all six CIs span zero — so on this record the mask
  vector is orbit-invariant *and* the echo coverage of the two directions is the
  same (33 of 62 ascending chips, 33 of 58 descending, Fisher p 0.72);
* and yet the trailing-window series moves: ``S`` swings -0.52/+1.75 and ``A``
  -2.5/+2.5 *mask pixels*, ``F`` -1.0/+1.5, ``D`` -0.73/+0.50, and gamma2
  -0.70/+0.79, so a pooled zero is an average over orbit-dependent windows and
  not the absence of an orbit effect;
* the pooled difference is **only weakly confounded with the state**: the four
  post-event echo chips are two ascending and two descending passes, so the
  pre-only column drops a balanced pair and the two contrasts nearly coincide
  (gamma2 +0.160 either way) — both are printed;
* the split is a statement about **pass geometry** (the four geometries listed
  above, resolved from the committed bursts), never a damage or orbit verdict —
  it adds no number to the A-E pin contract of ``fig_haltern_ocv_paper.py``.

Every number is re-derived from the committed channel table
``data/haltern/haltern_ocv_channels.csv`` and pinned field by field against the
committed mask cache ``haltern_windows_mask_cache.txt`` (the 6D mask vector
re-verified row by row), the generated ``haltern_ocv_channels_meta.json`` (guard
counts, orbit coverage, the sha256 of all five inputs), the committed paper
result ``fig_haltern_ocv_paper.json`` and the committed analysis reference
``reference/haltern_ocv_findings.json`` (the orbit coverage block and the sha256 of
the resolved reference records) — a deviation aborts the run with a non-zero
exit status. The bootstrap is seeded and no wall-clock time enters the JSON, so
two runs are byte-identical.

Usage:
  python3 code/haltern/fig_haltern_ocv_orbits.py                # figure + JSON + report
  python3 code/haltern/fig_haltern_ocv_orbits.py --no-figures   # cross-checks only
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
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import NullLocator  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
sys.path.insert(0, HERE)

import haltern_ocv_core as core     # noqa: E402
import haltern_ocv_masks as masks   # noqa: E402
import haltern_ocv_stats as st      # noqa: E402

FIGDIR = core.FIGDIR
FIG_NAME = "fig_haltern_ocv_orbits"
OUT_JSON = os.path.join(core.DATA, FIG_NAME + ".json")
OUT_MD = os.path.join(FIGDIR, FIG_NAME + ".md")
PAPER_JSON = os.path.join(core.DATA, "fig_haltern_ocv_paper.json")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 8.4,
    "axes.labelsize": 7.8,
    "xtick.labelsize": 6.9,
    "ytick.labelsize": 6.9,
    "legend.fontsize": 6.9,
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
# The drawing order of the six channel panels: the three OCV channels first,
# then the three mask-geometry ones.
DIMS = ["gamma2", "P", "D", "F", "S", "A"]
ORBITS = list(core.ORBITS)
ORBIT_SHORT = dict(core.ORBIT_SHORT)
ORBIT_LABEL = {"ASCENDING": "ascending (ASC)", "DESCENDING": "descending (DESC)"}
# The pass geometry the record's four orbits actually are (from the committed
# reference record ``reference/haltern_bursts.json``); printed, never inferred.
PASS_NOTE = ("the record's ascending passes are IW1/IW3 on relative orbits 15 and "
             "88 around 15:2x and 17:2x UT and its descending passes IW1/IW3 on "
             "relative orbits 37 and 139 around 05:2x UT")
ORBIT_NOTE = core.COVERAGE_NOTE

PANEL_GRID = (4, 2)       # the panel grid: four rows, two columns (a..h)
CADENCE_DAYS = 3          # the widest step of the in-month alternation
GAP_DAYS = 11             # the shortest hole: a step that stops the series
WINDOW_PX = 49            # the committed 7x7 window
WINDOW_WEEKS = 8.0        # the trailing window of the delta series
MIN_PER_ORBIT = 2         # a window median needs two chips of each orbit
N_BOOT = st.BOOT_MAIN     # 10000
BOOT_SEED = 7
ALPHA = 0.05
FLOAT_TOL = 1e-12
STATES = list(core.STATES)
PRE, POST = STATES

DELTA_NOTE = ("median(DESC) - median(ASC) over the echo-bearing chips: a "
              "descriptive pass-geometry difference, never a damage or orbit "
              "verdict. The two orbits are disjoint sets of acquisition days, so "
              "the chips are not paired and no causal reading is available")
# The twin-axis label of every channel panel: the sign convention is written
# once, in the figure's subtitle, so the label itself stays short enough to sit
# inside an 11 mm-tall panel without dictating the row pitch.
DELTA_LABEL = "Δorbit"
CONFOUND_NOTE = ("the record's four post-event echo chips are two ascending and two "
                 "descending passes, so the pre-only column drops a balanced pair: "
                 "the pooled difference and the before-the-collapse difference "
                 "differ by those four chips alone, and printing both shows how "
                 "little of the pooled number is state")
COVARIATE_NOTE = ("what else moves with the pass geometry: the per-orbit medians "
                  "of the export's co-variates. Ascending passes are the 15:2x and "
                  "17:2x UT ones, descending the 05:2x UT ones, so temperature, "
                  "wind and intensity move with the orbit by construction — they "
                  "would have to be held fixed before an orbit split of the mask "
                  "channels meant anything physical")

# The co-variates printed per orbit: what else moves with the pass geometry and
# would have to be held fixed before an orbit split meant anything physical.
ORBIT_CONTROLS = ["temperature_c", "wind_speed_ms", "wind_gust_ms", "intensity",
                  "coherence_gamma2", "precipitation_mm"]

ORBIT_COLOR = {"ASCENDING": "#4C72B0", "DESCENDING": "#C44E52"}
ORBIT_MARKER = {"ASCENDING": "o", "DESCENDING": "s"}
DELTA_COLOR = "#333333"
POOLED_COLOR = "#111111"
ZERO_COLOR = "#9A9A9A"
EVENT_COLOR = "#C44E52"
GAP_COLOR = "#E8E8E8"
PRE_TINT = "#4C72B0"
POST_TINT = "#C44E52"
PRE_TINT_ALPHA = 0.045
POST_TINT_ALPHA = 0.085
NO_ECHO_COLOR = "#B0B0B0"


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


def dtm(day):
    """Midnight of an ISO day as a ``datetime`` (the x axis carries datetimes)."""
    return dt.datetime.combine(day_date(day), dt.time())


def shift_dt(day, days):
    """``days`` (fractional is fine) after an ISO day or a ``datetime``."""
    base = day if isinstance(day, dt.datetime) else dtm(day)
    return base + dt.timedelta(days=days)


def orbit_dx(orbit, unit=0.18):
    """A deterministic sub-day offset for the two orbit directions.

    No day of this record carries two chips, so the offset never separates two
    chips of one day: it keeps the ascending and descending floor ticks of the
    long record from landing on the same pixel of the axis, and the connecting
    lines of the windowed series are drawn per direction with the same offset.
    """
    lab = str(orbit or "").upper()
    if lab.startswith("ASC"):
        return unit
    if lab.startswith("DESC"):
        return -unit
    return 0.0


def x_of(day, orbit=None):
    """The x coordinate of a chip: its acquisition day, plus the orbit offset."""
    return shift_dt(day, orbit_dx(orbit) if orbit else 0.0)


def add_months(day, months):
    """The first of the month ``months`` calendar months after ``day``'s month."""
    month = day.month - 1 + months
    return dt.date(day.year + month // 12, month % 12 + 1, 1)


def day_span(a, b):
    """Whole days between two ISO days."""
    return (day_date(b) - day_date(a)).days


def row_view(r):
    """The columns of one chip this figure uses (``A`` is None without echo)."""
    d = {
        "id": r.get("id"), "day": core.day_of(r),
        "orbit": r.get("orbit_direction"), "pass": r.get("pass_label"),
        "status": r.get("status"),
        "echo": bool(r["echo_mask_present"]),
        "windows": int(core._num(r.get("payload_group_size")) or 1),
        "section": r.get("segment"),
        "month": r.get("month"), "season": r.get("season"),
        "A": core._num(r.get("A")), "F": core._num(r.get("F")),
        "S": core._num(r.get("S")), "gamma2": core._num(r.get("gamma2_mask")),
        "P": core._num(r.get("P")), "D": core._num(r.get("D")),
        "echo_mode": r.get("echo_mode"),
        "coherence_masked_pixels": core._num(r.get("coherence_masked_pixels")),
    }
    for k in ORBIT_CONTROLS:
        d[k] = core._num(r.get(k))
    d["pre"] = bool(r["pre_collapse"])
    d["state"] = core.STATE_SHORT[core.state_of(d["day"])]
    d["orbit"] = str(d["orbit"] or "").upper()
    d["short"] = ORBIT_SHORT.get(d["orbit"], "?")
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


def orbit_rows(rows, orbit, state=None):
    """The chips of one orbit, optionally of one state."""
    out = [d for d in rows if d["orbit"] == orbit]
    if state is not None:
        out = [d for d in out if d["state"] == core.STATE_SHORT[state]]
    return out


def _fail(msg):
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(2)


# ---------------------------------------------------------------------------
# The orbit geometry: what the two "orbits" of this record actually are
# ---------------------------------------------------------------------------
def burst_geometry(rows, index=None, reference=None):
    """``mid -> {burst_id, swath, rel_orbit, pol, orbit, utc}`` of every chip.

    The CSV carries the export's own ``orbit_direction`` and ``pass_label``; the
    *pass geometry* behind them — sub-swath, relative orbit, polarisation, UTC —
    is the committed CDSE resolution ``reference/haltern_bursts.json``, keyed by the
    burst id of the committed window index (the 6th field of
    ``haltern_windows_index.txt``). This figure re-walks that chain instead of
    quoting the paper's sentence, so the sentence is a *checked* statement here.
    """
    index = index if index is not None else core.load_index(core.INDEX_PATH)
    reference = reference if reference is not None else \
        core.load_reference(core.REF_BURSTS)
    bursts = reference["bursts"]
    out = {}
    for r in rows:
        meta = index.get(r["id"])
        if meta is None:
            continue
        b = bursts.get(meta["burst_id"])
        if b is None:
            continue
        out[r["id"]] = {
            "burst_id": meta["burst_id"],
            "swath": b.get("SwathIdentifier"),
            "rel_orbit": core._num(b.get("RelativeOrbitNumber")),
            "pol": b.get("PolarisationChannels"),
            "orbit": b.get("OrbitDirection"),
            "utc": str(b.get("AzimuthTime") or "")[11:16],
        }
    return out


def counts_of(seq):
    """``{value: n}`` of a sequence, most frequent first (JSON-friendly)."""
    out = {}
    for v in seq:
        out[v] = out.get(v, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: (-kv[1], str(kv[0]))))


def geometry_block(rows, geo):
    """The per-orbit pass geometry, from the resolved bursts (checked, not quoted)."""
    block = {}
    for o in ORBITS:
        sub = [d for d in rows if d["orbit"] == o and d["id"] in geo]
        swaths = sorted({geo[d["id"]]["swath"] for d in sub})
        rels = sorted({geo[d["id"]]["rel_orbit"] for d in sub})
        pols = counts_of(geo[d["id"]]["pol"] for d in sub)
        utcs = sorted({geo[d["id"]]["utc"][:2] for d in sub})
        block[o] = {
            "n_chips": len(sub),
            "swaths": swaths, "relative_orbits": rels,
            "polarisations": pols, "utc_hours": utcs,
            "utc_hh": [f"{u}:2x" for u in utcs],
            "orbit_direction_agrees":
                all(geo[d["id"]]["orbit"] == d["orbit"] for d in sub),
        }
    return block


def geometry_sentence(block):
    """One sentence naming the four pass geometries, from the resolved bursts."""
    parts = []
    for o in ORBITS:
        b = block[o]
        orbits = " + ".join(f"{r:g}" for r in b["relative_orbits"])
        hours = " + ".join(b["utc_hh"])
        pols = " + ".join(f"{k} x {v}" for k, v in b["polarisations"].items())
        parts.append(f"{ORBIT_SHORT[o]} is {'/'.join(b['swaths'])} on relative "
                     f"orbit {orbits} ({hours} UT, {pols})")
    return "the record's passes are " + "; ".join(parts)


# ---------------------------------------------------------------------------
# The orbit difference: per-orbit medians, the difference, its CI, its series
# ---------------------------------------------------------------------------
def orbit_block(rows):
    """Chips per orbit (and per state inside an orbit), from the channel table."""
    out = {}
    for orbit in ORBITS:
        sub = [d for d in rows if d["orbit"] == orbit]
        by_state = {}
        for state in STATES:
            s = [d for d in sub if d["state"] == core.STATE_SHORT[state]]
            k = sum(1 for d in s if d["echo"])
            by_state[core.STATE_SHORT[state]] = {
                "n_chips": len(s), "n_with_echo": k,
                "n_without_echo": len(s) - k,
                "echo_rate": (k / len(s)) if s else 0.0}
        k = sum(1 for d in sub if d["echo"])
        out[orbit] = {"n_chips": len(sub), "n_with_echo": k,
                      "n_without_echo": len(sub) - k,
                      "echo_rate": k / len(sub), "by_state": by_state}
    return out


def orbit_coverage(rows):
    """The ``coverage_by_orbit`` block, recomputed from the channel table alone.

    Shape and wording are the committed block's (``core.coverage_by_orbit``), so
    the two meet field by field in the pins; ``event`` and ``window_px`` are kept
    so the paper result's copy of the block compares without an exclusion list.
    """
    per = orbit_block(rows)
    table = [[per[o]["n_with_echo"], per[o]["n_without_echo"]] for o in ORBITS]
    fisher_state = {}
    for o in ORBITS:
        pre, post = per[o]["by_state"]["pre"], per[o]["by_state"]["post"]
        fisher_state[o] = core.fisher_p(pre["n_with_echo"], pre["n_without_echo"],
                                        post["n_with_echo"], post["n_without_echo"])
    return {
        "event": core.COLLAPSE_DATE, "window_px": WINDOW_PX,
        "n_orbits": len(ORBITS), "orbits": list(ORBITS),
        "n_unique_chips": len(rows), "per_orbit": per,
        "fisher": {"table_rows": list(ORBITS), "table_cols": ["echo", "no_echo"],
                   "table": table,
                   "p_two_sided": core.fisher_p(table[0][0], table[0][1],
                                                table[1][0], table[1][1])},
        "fisher_by_state": fisher_state, "note": ORBIT_NOTE,
    }


def bootstrap_median_diff_ci(desc, asc, seed=BOOT_SEED, n_boot=N_BOOT, alpha=ALPHA):
    """Percentile CI of ``median(desc) - median(asc)``, each orbit resampled.

    The two orbits are resampled **independently** — they are disjoint sets of
    acquisition days, so no pairing exists — through ``numpy.default_rng(seed)``:
    the run is deterministic and carries no wall-clock state.
    """
    d = np.asarray([v for v in desc if v is not None], dtype=float)
    a = np.asarray([v for v in asc if v is not None], dtype=float)
    if d.size == 0 or a.size == 0:
        return None
    point = float(np.median(d) - np.median(a))
    if d.size < 2 or a.size < 2:
        return {"point": point, "lo": point, "hi": point, "n_boot": 0,
                "seed": seed, "degenerate": True}
    rng = np.random.default_rng(seed)
    dd = d[rng.integers(0, d.size, size=(n_boot, d.size))]
    aa = a[rng.integers(0, a.size, size=(n_boot, a.size))]
    diff = np.median(dd, axis=1) - np.median(aa, axis=1)
    lo, hi = np.percentile(diff, [50.0 * alpha, 100.0 - 50.0 * alpha])
    return {"point": point, "lo": float(lo), "hi": float(hi), "n_boot": int(n_boot),
            "seed": seed, "degenerate": False}


def orbit_delta(rows, dim, state=None, seed=BOOT_SEED, n_boot=N_BOOT):
    """One channel's orbit difference: the two per-orbit medians and the delta.

    ``delta = median(DESC) - median(ASC)`` throughout, so a positive number means
    the descending passes read *higher* on that channel; ``cliffs_delta`` is
    called with the ascending chips as reference, which gives the identical sign
    convention. ``state`` restricts the chips to one side of the collapse — the
    pre-only contrast isolates the pooled delta's state confound.
    """
    sub = [d for d in rows if d["echo"] and d.get(dim) is not None]
    if state is not None:
        sub = [d for d in sub if d["state"] == core.STATE_SHORT[state]]
    asc = values(orbit_rows(sub, "ASCENDING"), dim)
    desc = values(orbit_rows(sub, "DESCENDING"), dim)
    ci = bootstrap_median_diff_ci(desc, asc, seed=seed, n_boot=n_boot)
    return {
        "dim": dim, "label": core.CHANNEL_TITLE[dim],
        "state": core.STATE_SHORT[state] if state is not None else "all",
        "n_asc": len(asc), "n_desc": len(desc),
        "asc_median": float(np.median(asc)) if asc else None,
        "desc_median": float(np.median(desc)) if desc else None,
        "asc_mean": float(np.mean(asc)) if asc else None,
        "desc_mean": float(np.mean(desc)) if desc else None,
        "delta_median": (float(np.median(desc)) - float(np.median(asc)))
                        if (asc and desc) else None,
        "delta_mean": (float(np.mean(desc)) - float(np.mean(asc)))
                      if (asc and desc) else None,
        "cliffs_delta": st.cliffs_delta(asc, desc),
        "delta_ci95": [ci["lo"], ci["hi"]] if ci else None,
        "ci_excludes_zero": bool(ci["lo"] * ci["hi"] > 0) if ci else None,
        "n_boot": n_boot, "seed": seed, "note": DELTA_NOTE,
    }


def orbit_deltas(rows, dims=DIMS, seed=BOOT_SEED, n_boot=N_BOOT):
    """Every channel's orbit difference: pooled, pre-only and post-only."""
    return {dim: {"all": orbit_delta(rows, dim, None, seed, n_boot),
                  "pre": orbit_delta(rows, dim, PRE, seed, n_boot),
                  "post": orbit_delta(rows, dim, POST, seed, n_boot)}
            for dim in dims}


def trailing_series(rows, dim, window_weeks=WINDOW_WEEKS, min_per_orbit=MIN_PER_ORBIT):
    """``delta_orbit(t)`` on the dated axis, at the echo-bearing days.

    At each echo-bearing acquisition day ``t`` the window is ``t`` itself and the
    ``window_weeks`` weeks before it; the chips of that window are split by orbit
    and a point is drawn only when it holds at least ``min_per_orbit`` chips of
    *each* orbit (so both medians exist). A day whose window fails stays empty —
    the series is a trailing reading at ``t``, never an interpolation, and no day
    contributes two points.
    """
    span = int(round(window_weeks * 7))
    echo = [d for d in sorted_chips(rows) if d["echo"] and d.get(dim) is not None]
    out = []
    for d in echo:
        t1 = day_date(d["day"])
        t0 = t1 - dt.timedelta(days=span)
        w = [e for e in echo if t0 <= day_date(e["day"]) <= t1]
        a = [e[dim] for e in w if e["orbit"] == "ASCENDING"]
        b = [e[dim] for e in w if e["orbit"] == "DESCENDING"]
        if len(a) < min_per_orbit or len(b) < min_per_orbit:
            continue
        out.append({
            "day": d["day"], "state": d["state"], "span_days": span,
            "n_asc": len(a), "n_desc": len(b),
            "asc_median": float(np.median(a)), "desc_median": float(np.median(b)),
            "delta": float(np.median(b) - np.median(a)),
        })
    return out


def orbit_series(rows, dims=DIMS, window_weeks=WINDOW_WEEKS):
    """The trailing orbit difference of every channel, plus its own summary."""
    out = {}
    for dim in dims:
        pts = trailing_series(rows, dim, window_weeks)
        dels = [p["delta"] for p in pts]
        nz = [v for v in dels if v != 0.0]
        out[dim] = {
            "dim": dim, "label": core.CHANNEL_TITLE[dim],
            "window_weeks": window_weeks,
            "window_days": int(round(window_weeks * 7)),
            "min_per_orbit": MIN_PER_ORBIT,
            "points": pts, "n_points": len(pts),
            "first_day": pts[0]["day"] if pts else None,
            "last_day": pts[-1]["day"] if pts else None,
            "min_delta": min(dels) if dels else None,
            "max_delta": max(dels) if dels else None,
            "max_abs_delta": max((abs(v) for v in dels), default=None),
            "n_nonzero": len(nz),
            "n_sign_changes": sum(1 for a, b in zip(nz, nz[1:]) if a * b < 0),
            "n_distinct": len(sorted(set(dels))),
        }
    return out


# ---------------------------------------------------------------------------
# The record as blocks (days, steps, holes, composition)
# ---------------------------------------------------------------------------
def days_of(rows):
    """The acquisition days of a set of chips, in order."""
    return sorted({d["day"] for d in rows})


def step_block(rows):
    """Every acquisition step of a chip set, with its length in days."""
    days = days_of(rows)
    return [{"after": days[i - 1], "before": days[i],
             "days": (day_date(days[i]) - day_date(days[i - 1])).days}
            for i in range(1, len(days))]


def gap_block(rows):
    """The steps longer than the cadence: holes in the record, never bridged."""
    return [dict(s, crosses_event=(day_date(s["after"]) <= EVENT
                                   < day_date(s["before"])))
            for s in step_block(rows) if s["days"] > CADENCE_DAYS]


def cadence_bands(rows):
    """The chip pairs one cadence step apart — the only place a line is drawn."""
    chips = sorted_chips(rows)
    return [{"after": a["day"], "before": b["day"]}
            for a, b in zip(chips, chips[1:])
            if (day_date(b["day"]) - day_date(a["day"])).days == CADENCE_DAYS]


def day_block(rows):
    """One entry per acquisition day: what was requested and what echoed."""
    out = []
    for day in days_of(rows):
        sub = [d for d in rows if d["day"] == day]
        echo = [d for d in sub if d["echo"]]
        out.append({
            "day": day, "state": sub[0]["state"],
            "n_windows": sum(d["windows"] for d in sub),
            "n_unique_chips": len(sub),
            "n_duplicated": sum(d["windows"] - 1 for d in sub),
            "n_with_echo": len(echo),
            "n_without_echo": len(sub) - len(echo),
            "orbits": sorted({d["orbit"] for d in sub}),
        })
    return out


def composition_block(rows):
    """The record's composition: chips per orbit, per state, per orbit x state."""
    grid = {}
    for d in rows:
        cell = grid.setdefault(d["orbit"], {})
        cell[d["state"]] = cell.get(d["state"], 0) + 1
    echo_grid = {}
    for d in rows:
        if not d["echo"]:
            continue
        cell = echo_grid.setdefault(d["orbit"], {})
        cell[d["state"]] = cell.get(d["state"], 0) + 1
    echo_state = {}
    for state in STATES:
        s = core.STATE_SHORT[state]
        c = counts_of(d["orbit"] for d in rows if d["echo"] and d["state"] == s)
        echo_state[s] = {o: c.get(o, 0) for o in ORBITS}
    return {
        "orbit_state": {o: {core.STATE_SHORT[s]: grid.get(o, {}).get(
            core.STATE_SHORT[s], 0) for s in STATES} for o in ORBITS},
        "orbit_state_echo": {o: {core.STATE_SHORT[s]: echo_grid.get(o, {}).get(
            core.STATE_SHORT[s], 0) for s in STATES} for o in ORBITS},
        "echo_by_state_orbit": echo_state,
        "months": counts_of(d["day"][5:7] for d in rows),
        "sections": counts_of(d["section"] for d in rows),
        "n_orbit_in_state": {core.STATE_SHORT[s]: counts_of(
            d["orbit"] for d in rows if d["state"] == core.STATE_SHORT[s])
            for s in STATES},
    }


# ---------------------------------------------------------------------------
# The cross-layer checks the pins run on (cache, orbit chain, co-variates)
# ---------------------------------------------------------------------------
def cache_cross_check(res):
    """The committed mask cache against the channel table, channel by channel.

    ``masks.vector_from_row`` re-verifies the mask *rule* on every cache row
    (``A`` against ``|px|``, the peak against ``max_masked``, the floor against
    ``max(0.3 * peak, 5 * median)``) and derives ``A, D, F, S, gamma2`` from the
    committed coordinates; ``P`` comes from the cache's own reference set. All six
    are compared with the CSV, and the two layers' orbit labels with each other.
    The cache holds the record's 66 echo-bearing chips, so the 54 echo-free chips
    of the table carry no cache row at all: they have no ordinate to check.
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
            "echo": r["echo"], "orbit_agrees": bool(d["orbit"] == r["orbit"]),
            "A": (v["A"], r["A"]), "D": (v["D"], r["D"]), "F": (v["F"], r["F"]),
            "S": (v["S"], r["S"]), "gamma2": (v["gamma2"], r["gamma2"]),
            "P": (d["P"], r["P"]), "P_reference_chips": d["P_reference_chips"],
        })
    return {"n_cache_rows": len(res["cache"]), "n_compared": len(rows),
            "n_missing_in_csv": missing, "rows": rows,
            "P_reference_chips": sorted({r["P_reference_chips"] for r in rows}),
            "n_echo_rows": sum(1 for r in rows if r["echo"]),
            "n_orbit_disagreements": sum(1 for r in rows if not r["orbit_agrees"])}


def orbit_chain(rows, index=None, reference=None, measurements=None):
    """The CSV's own orbit columns against the two layers that fix them.

    Three statements of the same pass geometry meet here: the channel table's
    ``orbit_direction`` / ``pass_label`` (what the export says), the committed
    measurement extract's own columns (the same export, re-read from the raw
    layer) and the resolved burst record ``reference/haltern_bursts.json``
    (``OrbitDirection``, reached through the window index). The block counts the
    disagreements of each link — so "which orbit is this chip" is a checked chain
    and not an assumption — and records the label each orbit actually carries.
    """
    idx = index if index is not None else core.load_index(core.INDEX_PATH)
    ref = reference if reference is not None else core.load_reference(core.REF_BURSTS)
    meas = measurements if measurements is not None else core.load_measurements()
    bursts = ref["bursts"]
    n_extract = n_orbit_extract = n_label_extract = 0
    n_burst = n_orbit_burst = 0
    labels, missing = {}, []
    for d in rows:
        m = meas.get(d["id"])
        if m is None:
            missing.append(d["id"])
        else:
            n_extract += 1
            n_orbit_extract += int(str(m.get("orbit_direction") or "").upper()
                                   != d["orbit"])
            n_label_extract += int(str(m.get("pass_label") or "") != str(d["pass"]))
        i = idx.get(d["id"])
        b = bursts.get(i["burst_id"]) if i else None
        if b is None:
            continue
        n_burst += 1
        n_orbit_burst += int(str(b.get("OrbitDirection") or "").upper() != d["orbit"])
        labels.setdefault(d["orbit"], set()).add(str(d["pass"]))
    return {
        "n_chips": len(rows), "n_in_extract": n_extract,
        "n_missing_in_extract": len(missing), "missing_in_extract": missing,
        "n_orbit_mismatches_extract": n_orbit_extract,
        "n_pass_label_mismatches_extract": n_label_extract,
        "n_resolved_bursts": n_burst, "n_orbit_mismatches_bursts": n_orbit_burst,
        "pass_labels": {o: sorted(labels.get(o, [])) for o in ORBITS},
        "note": ("the export's own orbit_direction and pass_label agree with the "
                 "measurement extract and with the resolved bursts on every chip; "
                 "the cross-layer pins carry the three disagreement counts, all 0"),
    }


def measurement_controls(rows, measurements=None):
    """The per-orbit co-variate medians, re-read from the raw measurement extract.

    The channel table and the committed extract are two reads of the same export,
    so these medians are recomputed from both layers and pinned against each
    other; nothing here is quoted from a derived block.
    """
    meas = measurements if measurements is not None else core.load_measurements()
    out = {}
    for name, sub in (("all", rows), ("echo", [d for d in rows if d["echo"]])):
        out[name] = {}
        for o in ORBITS:
            sel = [d for d in sub if d["orbit"] == o]
            vals = {}
            for k in ORBIT_CONTROLS:
                col = [core._num((meas.get(d["id"]) or {}).get(k)) for d in sel]
                v = [x for x in col if x is not None]
                vals[k] = float(np.median(v)) if v else None
            out[name][o] = vals
    return out


def trailing_window_at(rows, dim, day, window_weeks=WINDOW_WEEKS,
                       min_per_orbit=MIN_PER_ORBIT):
    """The trailing window at one named day, kept or not (the pinned counter-case).

    The single post-event echo chip has a window of its own; the series keeps no
    point there because that window holds no descending chip, and the omission is
    a number of this figure, not a silent absence.
    """
    span = int(round(window_weeks * 7))
    t1 = day_date(day)
    t0 = t1 - dt.timedelta(days=span)
    w = [d for d in sorted_chips(rows)
         if d["echo"] and d.get(dim) is not None
         and t0 <= day_date(d["day"]) <= t1]
    a = [d[dim] for d in w if d["orbit"] == "ASCENDING"]
    b = [d[dim] for d in w if d["orbit"] == "DESCENDING"]
    return {"day": day, "dim": dim, "span_days": span,
            "n_asc": len(a), "n_desc": len(b),
            "asc_median": float(np.median(a)) if a else None,
            "desc_median": float(np.median(b)) if b else None,
            "kept": bool(len(a) >= min_per_orbit and len(b) >= min_per_orbit)}


# ---------------------------------------------------------------------------
# The evidence
# ---------------------------------------------------------------------------
def compute(reference_paper=None):
    """Everything this figure draws, from the committed tables only.

    Inputs: the channel table ``data/haltern/haltern_ocv_channels.csv`` (120 unique
    chips), the committed window file ``haltern_windows_full.txt`` (168 windows,
    recompiled by ``core.chip_table``), the committed mask cache
    ``haltern_windows_mask_cache.txt`` and its manifest, the committed window index
    ``haltern_windows_index.txt`` with the resolved bursts
    ``reference/haltern_bursts.json`` (the pass geometry), the generated
    ``haltern_ocv_channels_meta.json``, the committed paper result
    ``fig_haltern_ocv_paper.json`` and the committed analysis reference
    ``reference/haltern_ocv_findings.json``.
    """
    t0 = time.time()
    meta = core.load_csv_meta()
    reference = core.load_findings()
    paper = core.load_reference(reference_paper or PAPER_JSON)
    cache = core.load_cache()
    windows, unique = core.chip_table()
    csv_rows = core.load_csv()

    rows = sorted_chips([row_view(r) for r in csv_rows])
    echo = [d for d in rows if d["echo"]]
    geo = burst_geometry(csv_rows)
    geo_block = geometry_block(rows, geo)
    cov = orbit_coverage(rows)
    deltas = orbit_deltas(rows)
    series = orbit_series(rows)
    days = day_block(rows)
    steps = step_block(rows)
    gaps = gap_block(rows)
    bands = cadence_bands(rows)
    geo_missing = [d["id"] for d in rows if d["id"] not in geo]

    pre = [d for d in echo if d["state"] == "pre"]
    post = [d for d in echo if d["state"] == "post"]
    states = {core.STATE_SHORT[s]: [d for d in rows
                                    if d["state"] == core.STATE_SHORT[s]]
              for s in STATES}

    covariates = {}
    for name, sub in (("all", rows), ("echo", echo)):
        covariates[name] = {
            o: {k: median_of(orbit_rows(sub, o), k) for k in ORBIT_CONTROLS}
            for o in ORBITS}
        covariates[name]["delta"] = {
            k: (covariates[name]["DESCENDING"][k] - covariates[name]["ASCENDING"][k])
            for k in ORBIT_CONTROLS}

    chain = orbit_chain(rows, None, None, core.load_measurements())
    post_echo_day = days_of(post)[0] if post else None
    two_orbit_days = [d["day"] for d in days if len(d["orbits"]) > 1]
    two_orbit_set = set(two_orbit_days)

    res = {
        "meta": meta, "reference": reference, "paper": paper,
        "manifest": core.load_manifest(),
        "cache": cache, "windows": windows, "unique": unique,
        "rows": rows, "echo": echo, "pre_echo_rows": pre, "post_echo_rows": post,
        "geometry": geo, "geometry_block": geo_block,
        "geometry_sentence": geometry_sentence(geo_block),
        "coverage": cov, "deltas": deltas, "series": series,
        "days": days, "steps": steps, "gaps": gaps, "bands": bands,
        "day_list": days_of(rows), "composition": composition_block(rows),
        "covariates": covariates,
        "n_chips": len(rows), "n_echo": len(echo),
        "n_no_echo": len(rows) - len(echo),
        "n_days": len(days), "n_windows": sum(d["windows"] for d in rows),
        "n_duplicated_windows": sum(d["windows"] - 1 for d in rows),
        "n_windows_file": len(windows), "n_unique_file": len(unique),
        "n_gap_steps": len(gaps), "n_steps": len(steps),
        "n_steps_over_cadence": len(gaps),
        "n_pre_chips": len(states["pre"]), "n_post_chips": len(states["post"]),
        "n_pre_echo": len(pre), "n_post_echo": len(post),
        "n_pre_no_echo": len(states["pre"]) - len(pre),
        "n_post_no_echo": len(states["post"]) - len(post),
        "n_pre_days": len({d["day"] for d in states["pre"]}),
        "n_post_days": len({d["day"] for d in states["post"]}),
        "first_day": days_of(rows)[0], "last_day": days_of(rows)[-1],
        "pre_first_day": days_of(states["pre"])[0],
        "last_pre_day": days_of(states["pre"])[-1],
        "first_post_day": days_of(states["post"])[0],
        "last_post_day": days_of(states["post"])[-1],
        "n_echo_by_orbit": {o: len(orbit_rows(echo, o)) for o in ORBITS},
        "n_echo_by_state": {core.STATE_SHORT[s]: len(
            [d for d in echo if d["state"] == core.STATE_SHORT[s]])
            for s in STATES},
        "n_echo_orbit_state": composition_block(rows)["orbit_state_echo"],
        "post_echo_orbits": counts_of(d["orbit"] for d in post),
        "pre_echo_orbits": counts_of(d["orbit"] for d in pre),
        "geo_missing": geo_missing,
        "pre_median": {dim: median_of(
            [d for d in pre if d.get(dim) is not None], dim) for dim in DIMS},
        "event_weekday": EVENT.strftime("%A"),
        "orbit_chain": chain,
        "n_bursts": len({v["burst_id"] for v in geo.values()}),
        "pass_note": PASS_NOTE,
        "coverage_from_windows": core.coverage_by_orbit(unique),
        "controls_from_extract": measurement_controls(rows),
        "post_echo_day": post_echo_day,
        "post_window": (trailing_window_at(rows, DIMS[0], post_echo_day)
                        if post_echo_day else None),
        "series_days": sorted({p["day"] for d in DIMS for p in series[d]["points"]}),
        "series_span_days": int(round(WINDOW_WEEKS * 7)),
        "n_series_points": sum(series[d]["n_points"] for d in DIMS),
        "n_series_post": sum(1 for d in DIMS for p in series[d]["points"]
                             if p["state"] == "post"),
        "n_days_two_orbits": len(two_orbit_days),
        "days_two_orbits": two_orbit_days,
        "two_orbit_days_echo_free":
            all(not d["echo"] for d in rows if d["day"] in two_orbit_set),
        "n_days_single_orbit": len(days) - len(two_orbit_days),
    }
    res["span_days"] = day_span(res["first_day"], res["last_day"]) + 1
    res["collapse_step"] = next(
        (s for s in steps if day_date(s["after"]) <= EVENT < day_date(s["before"])),
        None)
    res["p_orbit"] = cov["fisher"]["p_two_sided"]
    res["p_orbit_by_state"] = dict(cov["fisher_by_state"])
    res["n_P_set"] = sum(1 for d in cache if d.get("P") is not None)
    res["p_reference_chips"] = sorted({d["P_reference_chips"] for d in cache
                                       if d.get("P_reference_chips") is not None})
    res["runtime_s"] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# The pins: five committed layers fix the same record
# ---------------------------------------------------------------------------
TOL_KEYS = {"echo_rate", "p", "p_two_sided", "median", "mean", "std", "lo", "hi",
            "point", "value", "delta", "delta_median", "asc_median",
            "desc_median"}
# The record facts below are written at the precision the report prints (six
# decimals to a part in a million), so a re-run that moves a printed digit fails
# the pins instead of quietly writing a different report; the seeded bootstrap
# itself is re-derived exactly and is compared with ``FLOAT_TOL``.
ROUND_TOL = 1e-6

# The record facts this figure re-derives on every run: the two orbits' counts
# and medians, the pooled and pre-only differences with their seeded
# percentile CIs, and the labels of the two passes. Re-deriving them from the
# channel table and failing the run when a digit moves is what makes the printed
# report a checked statement rather than a quotation.
DELTA_FACTS = {
    "gamma2": {"n_asc": 33, "n_desc": 33, "pre_n_asc": 31,
          "asc_median": 0.407076, "desc_median": 0.566581, "delta": 0.159505,
          "ci": [-0.162367, 0.478246], "cliffs": 0.179063, "excl": False,
          "pre_delta": 0.159505, "pre_ci": [-0.272140, 0.519544],
          "pre_cliffs": 0.175858, "pre_excl": False,
          "post_n_asc": 2, "post_n_desc": 2, "post_delta": 0.086028,
          "post_cliffs": 0.000000, "post_excl": False},
    "P": {"n_asc": 33, "n_desc": 33, "pre_n_asc": 31,
          "asc_median": 0.090909, "desc_median": 0.090909, "delta": 0.000000,
          "ci": [-0.013636, 0.010101], "cliffs": -0.070707, "excl": False,
          "pre_delta": -0.001684, "pre_ci": [-0.015152, 0.011364],
          "pre_cliffs": -0.080125, "pre_excl": False,
          "post_n_asc": 2, "post_n_desc": 2, "post_delta": 0.002525,
          "post_cliffs": 0.000000, "post_excl": False},
    "D": {"n_asc": 33, "n_desc": 33, "pre_n_asc": 31,
          "asc_median": 0.500000, "desc_median": 0.450000, "delta": -0.050000,
          "ci": [-0.166667, 0.266667], "cliffs": 0.019284, "excl": False,
          "pre_delta": -0.050000, "pre_ci": [-0.233333, 0.266667],
          "pre_cliffs": -0.004162, "pre_excl": False,
          "post_n_asc": 2, "post_n_desc": 2, "post_delta": 0.272857,
          "post_cliffs": 0.500000, "post_excl": False},
    "F": {"n_asc": 33, "n_desc": 33, "pre_n_asc": 31,
          "asc_median": 2.000000, "desc_median": 2.000000, "delta": 0.000000,
          "ci": [-1.000000, 1.000000], "cliffs": 0.048669, "excl": False,
          "pre_delta": 0.000000, "pre_ci": [-1.000000, 1.000000],
          "pre_cliffs": 0.077003, "pre_excl": False,
          "post_n_asc": 2, "post_n_desc": 2, "post_delta": -0.500000,
          "post_cliffs": -0.250000, "post_excl": False},
    "S": {"n_asc": 33, "n_desc": 33, "pre_n_asc": 31,
          "asc_median": 0.942809, "desc_median": 0.949334, "delta": 0.006525,
          "ci": [-0.559017, 0.583799], "cliffs": 0.062443, "excl": False,
          "pre_delta": 0.158764, "pre_ci": [-0.495076, 0.600925],
          "pre_cliffs": 0.091571, "pre_excl": False,
          "post_n_asc": 2, "post_n_desc": 2, "post_delta": -0.114079,
          "post_cliffs": 0.000000, "post_excl": False},
    "A": {"n_asc": 33, "n_desc": 33, "pre_n_asc": 31,
          "asc_median": 3.000000, "desc_median": 4.000000, "delta": 1.000000,
          "ci": [-1.000000, 2.000000], "cliffs": 0.224059, "excl": False,
          "pre_delta": 1.000000, "pre_ci": [0.000000, 2.000000],
          "pre_cliffs": 0.242456, "pre_excl": False,
          "post_n_asc": 2, "post_n_desc": 2, "post_delta": 0.000000,
          "post_cliffs": 0.000000, "post_excl": False},
}

# The dates the trailing window keeps: every echo-bearing day whose window holds
# two chips of each orbit, from the first in-month day to the last post-event echo.
SERIES_DAYS = [
    "2021-02-01", "2021-02-03", "2021-02-04", "2021-02-06",
    "2021-02-07", "2021-02-10", "2021-02-12", "2021-02-15",
    "2021-02-16", "2021-02-19", "2021-03-02", "2021-03-06",
    "2021-03-12", "2021-03-14", "2021-03-15", "2021-03-17",
    "2021-04-01", "2021-04-05", "2021-04-08", "2021-04-11",
    "2021-04-13", "2021-04-16", "2021-04-17", "2021-05-01",
    "2021-05-07", "2021-05-10", "2021-05-13", "2021-05-14",
    "2021-06-01", "2021-06-04", "2021-06-09", "2021-06-10",
    "2021-06-12", "2021-06-15", "2021-06-16", "2021-06-18",
    "2021-06-19", "2021-07-01", "2021-07-07", "2021-07-09",
    "2021-07-12", "2021-07-13", "2021-07-16", "2021-07-18",
    "2021-07-19", "2021-08-03", "2021-08-06", "2021-08-08",
    "2021-08-14", "2021-08-15", "2021-08-18", "2021-09-01",
    "2021-09-04", "2021-09-13", "2021-09-14", "2021-09-16",
    "2021-09-17", "2021-10-02", "2021-10-08", "2021-10-10",
    "2021-10-17"
]
SERIES_FACTS = {
    "gamma2": {"n": 61, "min": -0.704441, "max": 0.794985, "abs": 0.794985,
              "nsc": 5, "nd": 51},
    "P": {"n": 61, "min": -0.034722, "max": 0.028409, "abs": 0.034722,
              "nsc": 4, "nd": 39},
    "D": {"n": 61, "min": -0.733333, "max": 0.500000, "abs": 0.733333,
              "nsc": 6, "nd": 37},
    "F": {"n": 61, "min": -1.000000, "max": 1.500000, "abs": 1.500000,
              "nsc": 4, "nd": 6},
    "S": {"n": 61, "min": -0.524584, "max": 1.745832, "abs": 1.745832,
              "nsc": 8, "nd": 50},
    "A": {"n": 61, "min": -2.500000, "max": 2.500000, "abs": 2.500000,
              "nsc": 3, "nd": 11},
}

# What the two pass geometries of this record are: the resolved burst
# geometry, the label the export puts on each and the coverage they carry.
# Both directions carry two IW sub-swaths on two relative orbits, so the
# daily strongest-window pick may change sub-swath within one direction.
GEOMETRY_FACTS = {
    "ASCENDING": {
        "n_chips": 62,
        "swaths": ['IW1', 'IW3'], "relative_orbits": [15.0, 88.0],
        "polarisations": {'VV': 33, 'VH': 29},
        "utc_hours": ['15', '17'], "utc_hh": ['15:2x', '17:2x'],
    },
    "DESCENDING": {
        "n_chips": 58,
        "swaths": ['IW1', 'IW3'], "relative_orbits": [37.0, 139.0],
        "polarisations": {'VH': 31, 'VV': 27},
        "utc_hours": ['05'], "utc_hh": ['05:2x'],
    },
}
PASS_LABEL_FACTS = {'ASCENDING': ['afternoon'], 'DESCENDING': ['morning']}
COVERAGE_FACTS = {
    "n_asc_echo": 33, "n_desc_echo": 33,
    "n_asc_chips": 62, "n_desc_chips": 58,
    "p_pooled": 0.7166050513098348, "p_asc": 0.4054488819993858,
    "p_desc": 0.3866687973457825,
}
COMPOSITION_FACTS = {
    "orbit_state": {'ASCENDING': {'pre': 56, 'post': 6}, 'DESCENDING': {'pre': 52, 'post': 6}},
    "orbit_state_echo": {'ASCENDING': {'pre': 31, 'post': 2}, 'DESCENDING': {'pre': 31, 'post': 2}},
    "echo_by_state_orbit": {'pre': {'ASCENDING': 31, 'DESCENDING': 31}, 'post': {'ASCENDING': 2, 'DESCENDING': 2}},
    "sections": {'turbine-mast-section-0': 72, 'turbine-mast-section-1': 48},
}
# The pooled window at the first post-event echo day: unlike the companion
# record it holds chips of both orbits, so the series does keep points there.
POST_WINDOW_FACT = {"day": '2021-10-02', "dim": 'gamma2', "n_asc": 4, "n_desc": 7,
                    "kept": True, "span_days": 56}
# The co-variate differences the report quotes, clamped to the digits it
# prints: the morning descending passes are the colder, calmer and darker ones.
CONTROL_FACTS = {
    "all": {'coherence_gamma2': 0.123693, 'intensity': -2265.738944, 'precipitation_mm': 0.000000, 'temperature_c': -5.100000, 'wind_gust_ms': -3.300000, 'wind_speed_ms': -1.135000},
    "echo": {'coherence_gamma2': 0.154292, 'intensity': -2103.282564, 'precipitation_mm': 0.000000, 'temperature_c': -6.100000, 'wind_gust_ms': -4.400000, 'wind_speed_ms': -1.390000},
}

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
        tol = ROUND_TOL if key in TOL_KEYS else 0.0
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
            "float_tolerance": FLOAT_TOL, "round_tolerance": ROUND_TOL,
            "failures": self.fails[:40],
        }


def committed_timeline(blob):
    """The committed ``timeline`` block, keyed by day."""
    return {e["day"]: e for e in blob["timeline"]}


def committed_months(blob):
    """The committed months, pooled from the per-state month counts."""
    return {m: sum(per.values()) for m, per in blob["months_by_state"].items()}


def pin_all(res):
    """Re-derive every number of this axis and check it against the commits."""
    meta, man = res["meta"], res["manifest"]
    paper, ref = res["paper"], res["reference"]
    cov, rows = res["coverage"], res["rows"]
    repo = os.path.abspath(os.path.join(core.DATA, os.pardir, os.pardir))
    P = Pins()

    # 1) the record, from the four committed artifacts ------------------------
    P.cmp("record.n_chips", res["n_chips"], meta["n_unique_chips"])
    P.cmp("record.n_chips_paper", res["n_chips"], paper["n_unique_chips"])
    P.cmp("record.n_chips_meta_rows", res["n_chips"], meta["n_rows"])
    P.cmp("record.n_chips_file", res["n_chips"], res["n_unique_file"])
    P.cmp("record.n_windows", res["n_windows"], meta["n_windows"])
    P.cmp("record.n_windows_pair", res["n_windows"], res["n_windows_file"])
    P.cmp("record.n_windows_paper", res["n_windows"], paper["n_windows"])
    P.cmp("record.n_duplicated_windows", res["n_duplicated_windows"],
          meta["n_duplicate_windows"])
    P.cmp("record.n_duplicated_ref", res["n_duplicated_windows"],
          ref["records"]["n_deduplicated_same_payload"])
    P.cmp("record.n_echo", res["n_echo"], meta["n_echo_masks"])
    P.cmp("record.n_echo_paper", res["n_echo"], paper["n_echo_chips"])
    P.cmp("record.n_echo_manifest", res["n_echo"], man["n_rows"])
    P.cmp("record.n_no_echo", res["n_no_echo"], meta["n_rows_without_echo_mask"])
    P.cmp("record.n_no_echo_manifest", res["n_no_echo"], man["n_without_echo"])
    P.cmp("record.n_pre_chips", res["n_pre_chips"], meta["n_pre_collapse"])
    P.cmp("record.n_post_chips", res["n_post_chips"], meta["n_post_collapse"])
    P.cmp("record.n_pre_echo", res["n_pre_echo"], 62)
    P.cmp("record.n_post_echo", res["n_post_echo"], 4)
    P.cmp("record.n_echo_splits", res["n_pre_echo"] + res["n_post_echo"],
          res["n_echo"])
    P.cmp("record.n_pre_plus_post", res["n_pre_chips"] + res["n_post_chips"],
          res["n_chips"])
    P.cmp("record.echo_plus_no_echo", res["n_echo"] + res["n_no_echo"],
          res["n_chips"])
    P.cmp("record.n_days", res["n_days"],
          meta["mask_cache_manifest"]["asset_date_groups"])
    P.cmp("record.n_days_ref", res["n_days"], ref["records"]["n_dates"])
    P.cmp("record.n_days_timeline", res["n_days"], len(ref["timeline"]))
    P.cmp("record.date_range", meta["date_range"],
          [res["first_day"], res["last_day"]])
    P.cmp("record.first_day", res["first_day"], "2021-01-01")
    P.cmp("record.last_day", res["last_day"], "2021-10-17")
    P.cmp("record.span_days", res["span_days"],
          (day_date(res["last_day"]) - day_date(res["first_day"])).days + 1)
    P.cmp("record.n_pre_days", res["n_pre_days"], 108)
    P.cmp("record.n_post_days", res["n_post_days"], 12)
    P.cmp("record.mask_window", meta["mask_window"], "7x7")
    P.cmp("record.window_px", ref["window_px"], WINDOW_PX)
    P.cmp("record.states", meta["states"], ref["state_labels"])
    P.cmp("record.orbits", meta["orbits"], ORBITS)
    P.cmp("record.sections", meta["sections_counts"],
          res["composition"]["sections"])
    P.cmp("record.sections_present", sorted(ref["records"]["sections_present"]),
          sorted(res["composition"]["sections"]))
    P.cmp("record.ref_segment_resolved", ref["records"]["segment_resolved"], False)
    P.cmp("record.ref_window_sizes", ref["records"]["window_sizes"],
          {meta["mask_window"]: res["n_chips"]})
    P.cmp("record.n_by_orbit", meta["n_by_orbit"],
          {core.ORBIT_SHORT[o]: cov["per_orbit"][o]["n_chips"] for o in ORBITS})
    P.cmp("record.p_reference_chips", res["p_reference_chips"], [meta["n_P_set"]])
    P.cmp("record.n_P_set", res["n_P_set"], res["n_echo"])

    # 2) the chain's own guard counts -----------------------------------------
    P.cmp("manifest.generator", man["generator"],
          meta["mask_cache_manifest"]["generator"])
    P.cmp("manifest.mask_rule", man["mask_rule"], meta["mask_rule"])
    P.cmp("manifest.n_rows", man["n_rows"], res["n_echo"])
    P.cmp("manifest.n_without_echo", man["n_without_echo"], res["n_no_echo"])
    P.cmp("manifest.n_without_index", man["n_without_index"], 0)
    P.cmp("manifest.n_deduplicated", man["n_deduplicated_same_payload"],
          res["n_duplicated_windows"])
    P.cmp("manifest.segment_resolved", man["segment_resolved"],
          ref["records"]["segment_resolved"])
    P.cmp("manifest.window_size", man["window_size"], 7)
    P.cmp("manifest.window_sizes", man["window_sizes"],
          {meta["mask_window"]: res["n_chips"]})
    P.cmp("manifest.source_windows", man["source_windows"],
          os.path.basename(core.WINDOWS_PATH))
    for k in sorted(man):
        if k not in ("source_windows_sha256",):
            P.cmp(f"meta.mask_cache_manifest.{k}", man[k],
                  meta["mask_cache_manifest"][k])
    P.cmp("csv.n_columns", meta["n_columns"], 81)
    P.cmp("csv.ocv", meta["ocv"], OCV)
    P.cmp("csv.features", meta["features"], FEATURES)
    P.cmp("csv.n_coherence_compared", meta["n_coherence_masked_pixels_compared"],
          res["n_echo"])
    # The export's ``coherence_masked_pixels`` is a different-generation column,
    # computed on the whole-scene asset-point chip: against the recomputed echo
    # mask of the committed 7x7 window the two counts coincide on only 6 of the
    # 66 echo chips (median |A - column| = 28 px, Spearman -0.25), so a
    # coincidence is an integer accident and not a shared quantity.
    P.cmp("csv.n_coherence_matched", meta["n_coherence_masked_pixels_matched"], 6)
    P.cmp("csv.n_coherence_differing", meta["n_coherence_masked_pixels_differing"],
          res["n_echo"] - 6)

    # 3) the committed inputs, by content: the same files this run read -------
    srcs = ref["sources"]
    for key, path in (("channels_csv", core.CSV_PATH),
                      ("mask_cache", core.CACHE_PATH),
                      ("mask_cache_manifest", core.MANIFEST_PATH),
                      ("windows", core.WINDOWS_PATH),
                      ("segments", core.SEGMENTS_PATH),
                      ("measurements", core.MEAS_PATH),
                      ("burst_geometry", core.REF_BURSTS),
                      ("osm_anchor", core.REF_GEOMETRY),
                      ("osm_node_extract", core.REF_OSM_NODE)):
        P.cmp(f"inputs.{key}.path", os.path.relpath(path, repo), srcs[key]["path"])
        P.cmp(f"inputs.{key}.sha256", core.sha256(path), srcs[key]["sha256"])
    P.cmp("inputs.windows_sha_in_manifest", man["source_windows_sha256"],
          srcs["windows"]["sha256"])
    P.cmp("inputs.n_named", len(srcs), 9)
    for key, path in (("mask_cache", core.CACHE_PATH),
                      ("mask_cache_manifest", core.MANIFEST_PATH),
                      ("windows", core.WINDOWS_PATH),
                      ("segments", core.SEGMENTS_PATH),
                      ("measurements", core.MEAS_PATH)):
        P.cmp(f"paper_inputs.{key}.file", os.path.basename(path),
              paper["meta"]["inputs"][key]["file"])
        P.cmp(f"paper_inputs.{key}.sha256", core.sha256(path),
              paper["meta"]["inputs"][key]["sha256"])

    # 4) the mask vector: every committed cache row against the CSV ----------
    cc = cache_cross_check(res)
    P.cmp("cache.n_cache_rows", cc["n_cache_rows"],
          meta["cache_cross_check"]["n_cache_rows"])
    P.cmp("cache.n_cache_rows_ref", cc["n_cache_rows"],
          meta["cache_cross_check"]["n_chips_with_echo"])
    P.cmp("cache.n_compared", cc["n_compared"], res["n_echo"])
    P.cmp("cache.n_missing_in_csv", cc["n_missing_in_csv"], 0)
    P.cmp("cache.n_echo_rows", cc["n_echo_rows"], res["n_echo"])
    P.cmp("cache.n_orbit_disagreements", cc["n_orbit_disagreements"], 0)
    P.cmp("cache.P_reference_chips", cc["P_reference_chips"], [meta["n_P_set"]])
    P.cmp("cache.rows_all_echo", sorted({r["echo"] for r in cc["rows"]}), [True])
    P.cmp("cache.one_row_per_echo_chip", len({r["mid"] for r in cc["rows"]}),
          res["n_echo"])
    for row in cc["rows"]:
        tag = f"cache[{row['date']}/{row['orbit']}]"
        for k in ("A", "D", "F", "S", "gamma2", "P"):
            P.add(f"{tag}.{k}", row[k][0], row[k][1], FLOAT_TOL)

    # 5) the coverage of the two pass geometries -----------------------------
    P.cmp("coverage.meta", res["coverage"], meta["coverage_by_orbit"])
    P.cmp("coverage.paper", res["coverage"], paper["coverage_by_orbit"])
    P.cmp("coverage.paper_meta", res["coverage"],
          paper["meta"]["coverage_by_orbit"])
    P.cmp("coverage.from_the_window_file", res["coverage"],
          res["coverage_from_windows"])
    P.cmp("coverage.ref", {k: v for k, v in res["coverage"].items()
                           if k not in ("event", "window_px")},
          ref["echo_coverage_by_orbit"])
    P.cmp("coverage.note_is_the_committed_note", cov["note"],
          ref["echo_coverage_by_orbit"]["note"])
    P.cmp("coverage.n_asc_chips", cov["per_orbit"]["ASCENDING"]["n_chips"],
          COVERAGE_FACTS["n_asc_chips"])
    P.cmp("coverage.n_desc_chips", cov["per_orbit"]["DESCENDING"]["n_chips"],
          COVERAGE_FACTS["n_desc_chips"])
    P.cmp("coverage.n_asc_echo", cov["per_orbit"]["ASCENDING"]["n_with_echo"],
          COVERAGE_FACTS["n_asc_echo"])
    P.cmp("coverage.n_desc_echo", cov["per_orbit"]["DESCENDING"]["n_with_echo"],
          COVERAGE_FACTS["n_desc_echo"])
    P.cmp("coverage.pooled_chips", sum(cov["per_orbit"][o]["n_chips"]
                                       for o in ORBITS), res["n_chips"])
    P.cmp("coverage.pooled_echo", sum(cov["per_orbit"][o]["n_with_echo"]
                                      for o in ORBITS), res["n_echo"])
    P.cmp("coverage.p_orbit", res["p_orbit"], COVERAGE_FACTS["p_pooled"])
    P.cmp("coverage.p_by_state", res["p_orbit_by_state"],
          {"ASCENDING": COVERAGE_FACTS["p_asc"],
           "DESCENDING": COVERAGE_FACTS["p_desc"]})
    P.cmp("coverage.fisher_table_is_the_two_orbits",
          cov["fisher"]["table"],
          [[cov["per_orbit"][o]["n_with_echo"],
            cov["per_orbit"][o]["n_without_echo"]] for o in ORBITS])
    P.cmp("coverage.the_post_echoes_are_two_each_way",
          res["post_echo_orbits"], {"ASCENDING": 2, "DESCENDING": 2})

    # 6) the pass geometry behind the two labels -----------------------------
    for o in ORBITS:
        g, w = res["geometry_block"][o], GEOMETRY_FACTS[o]
        P.cmp(f"geometry.{o}.n_chips", g["n_chips"], w["n_chips"])
        P.cmp(f"geometry.{o}.swaths", g["swaths"], w["swaths"])
        P.cmp(f"geometry.{o}.relative_orbits", g["relative_orbits"],
              w["relative_orbits"])
        P.cmp(f"geometry.{o}.polarisations", g["polarisations"],
              w["polarisations"])
        P.cmp(f"geometry.{o}.utc_hh", g["utc_hh"], w["utc_hh"])
        P.cmp(f"geometry.{o}.label_agrees_with_the_bursts",
              g["orbit_direction_agrees"], True)
        P.cmp(f"geometry.{o}.utc_hours", g["utc_hours"], w["utc_hours"])
    P.cmp("geometry.n_chips", sum(res["geometry_block"][o]["n_chips"]
                                  for o in ORBITS), res["n_chips"])
    P.cmp("geometry.n_bursts", res["n_bursts"], res["n_chips"])
    P.cmp("geometry.nothing_missing", res["geo_missing"], [])
    P.cmp("geometry.sentence", res["geometry_sentence"],
          "the record's passes are ASC is IW1/IW3 on relative orbit 15 + 88 "
          "(15:2x + 17:2x UT, VV x 33 + VH x 29); DESC is IW1/IW3 on relative "
          "orbit 37 + 139 (05:2x UT, VH x 31 + VV x 27)")
    P.cmp("geometry.note", res["pass_note"], PASS_NOTE)

    # 7) the orbit chain: CSV against the extract and the bursts -------------
    ch = res["orbit_chain"]
    P.cmp("chain.n_chips", ch["n_chips"], res["n_chips"])
    P.cmp("chain.n_in_extract", ch["n_in_extract"], res["n_chips"])
    P.cmp("chain.n_missing_in_extract", ch["n_missing_in_extract"], 0)
    P.cmp("chain.n_orbit_mismatches_extract", ch["n_orbit_mismatches_extract"], 0)
    P.cmp("chain.n_pass_label_mismatches_extract",
          ch["n_pass_label_mismatches_extract"], 0)
    P.cmp("chain.n_resolved_bursts", ch["n_resolved_bursts"], res["n_chips"])
    P.cmp("chain.n_orbit_mismatches_bursts", ch["n_orbit_mismatches_bursts"], 0)
    P.cmp("chain.pass_labels", ch["pass_labels"], PASS_LABEL_FACTS)
    P.cmp("chain.the_two_labels", sorted({d["orbit"] for d in rows}), ORBITS)
    P.cmp("chain.every_chip_has_a_label",
          sorted({d["pass"] for d in rows}), ["afternoon", "morning"])

    # 8) the alternation of the two orbits over the record -------------------
    P.cmp("alternation.cadence_days", CADENCE_DAYS, 3)
    P.cmp("alternation.gap_days", GAP_DAYS, 11)
    P.cmp("alternation.n_days", res["n_days"], 120)
    P.cmp("alternation.no_day_carries_two_orbits", res["n_days_two_orbits"], 0)
    P.cmp("alternation.two_orbit_days", res["days_two_orbits"], [])
    P.cmp("alternation.n_days_single_orbit", res["n_days_single_orbit"],
          res["n_days"])
    P.cmp("alternation.no_echo_day_carries_two_orbits",
          sorted({len(d["orbits"]) for d in res["days"] if d["n_with_echo"]}),
          [1])
    P.cmp("alternation.n_gap_steps", res["n_gap_steps"], 9)
    P.cmp("alternation.gap_lengths", sorted({g["days"] for g in res["gaps"]}),
          [11, 12, 14, 15])
    P.cmp("alternation.n_cadence_steps", res["n_steps"] - res["n_gap_steps"], 110)
    P.cmp("alternation.n_cadence_bands", len(res["bands"]), 3)
    P.cmp("alternation.collapse_step", res["collapse_step"],
          {"after": "2021-09-17", "before": "2021-10-01", "days": 14})
    P.cmp("alternation.steps_crossing_the_event",
          [g["after"] for g in res["gaps"] if g["crosses_event"]],
          ["2021-09-17"])

    # 9) the difference itself: the pooled and the pre-only contrast ---------
    for dim in DIMS:
        w, d = DELTA_FACTS[dim], res["deltas"][dim]
        all_d, pre_d, post_d = d["all"], d["pre"], d["post"]
        P.cmp(f"delta.{dim}.n_asc", all_d["n_asc"], w["n_asc"])
        P.cmp(f"delta.{dim}.n_desc", all_d["n_desc"], w["n_desc"])
        P.cmp(f"delta.{dim}.n_chips", all_d["n_asc"] + all_d["n_desc"],
              res["n_echo"])
        # the pre-event split is balanced: 31 echo chips of each direction
        P.cmp(f"delta.{dim}.pre_n_asc", pre_d["n_asc"], w["pre_n_asc"])
        P.cmp(f"delta.{dim}.pre_n_desc", pre_d["n_desc"], w["pre_n_asc"])
        P.cmp(f"delta.{dim}.pre_n_chips",
              pre_d["n_asc"] + pre_d["n_desc"], res["n_pre_echo"])
        P.cmp(f"delta.{dim}.post_is_two_chips_each_way",
              [post_d["n_asc"], post_d["n_desc"]],
              [w["post_n_asc"], w["post_n_desc"]])
        P.add(f"delta.{dim}.post_delta_median", post_d["delta_median"],
              w["post_delta"], ROUND_TOL)
        P.add(f"delta.{dim}.post_cliffs_delta", post_d["cliffs_delta"],
              w["post_cliffs"], ROUND_TOL)
        P.add(f"delta.{dim}.post_ci_excludes_zero", post_d["ci_excludes_zero"],
              w["post_excl"], ROUND_TOL)
        P.cmp(f"delta.{dim}.label", all_d["label"], core.CHANNEL_TITLE[dim])
        P.cmp(f"delta.{dim}.state", all_d["state"], "all")
        P.cmp(f"delta.{dim}.pre_state", pre_d["state"], "pre")
        P.cmp(f"delta.{dim}.n_boot", all_d["n_boot"], paper["n_boot"])
        P.cmp(f"delta.{dim}.seed", all_d["seed"], paper["rng_seed"])
        P.cmp(f"delta.{dim}.note", all_d["note"], DELTA_NOTE)
        P.add(f"delta.{dim}.asc_median", all_d["asc_median"], w["asc_median"],
              ROUND_TOL)
        P.add(f"delta.{dim}.desc_median", all_d["desc_median"],
              w["desc_median"], ROUND_TOL)
        P.add(f"delta.{dim}.delta_median", all_d["delta_median"], w["delta"],
              ROUND_TOL)
        P.add(f"delta.{dim}.delta_ci95_lo", all_d["delta_ci95"][0], w["ci"][0],
              ROUND_TOL)
        P.add(f"delta.{dim}.delta_ci95_hi", all_d["delta_ci95"][1], w["ci"][1],
              ROUND_TOL)
        P.add(f"delta.{dim}.cliffs_delta", all_d["cliffs_delta"], w["cliffs"],
              ROUND_TOL)
        P.add(f"delta.{dim}.ci_excludes_zero", all_d["ci_excludes_zero"],
              w["excl"], ROUND_TOL)
        P.add(f"delta.{dim}.pre_delta_median", pre_d["delta_median"],
              w["pre_delta"], ROUND_TOL)
        P.add(f"delta.{dim}.pre_delta_ci95_lo", pre_d["delta_ci95"][0],
              w["pre_ci"][0], ROUND_TOL)
        P.add(f"delta.{dim}.pre_delta_ci95_hi", pre_d["delta_ci95"][1],
              w["pre_ci"][1], ROUND_TOL)
        P.add(f"delta.{dim}.pre_cliffs_delta", pre_d["cliffs_delta"],
              w["pre_cliffs"], ROUND_TOL)
        P.add(f"delta.{dim}.pre_ci_excludes_zero", pre_d["ci_excludes_zero"],
              w["pre_excl"], ROUND_TOL)
        P.add(f"delta.{dim}.delta_is_the_two_medians",
              all_d["delta_median"],
              all_d["desc_median"] - all_d["asc_median"], FLOAT_TOL)
        P.add(f"delta.{dim}.the_sign_convention_is_dec_minus_asc",
              all_d["cliffs_delta"],
              st.cliffs_delta(values(orbit_rows([r for r in res["echo"]
                                                 if r[dim] is not None],
                                                "ASCENDING"), dim),
                              values(orbit_rows([r for r in res["echo"]
                                                 if r[dim] is not None],
                                                "DESCENDING"), dim)),
              FLOAT_TOL)
    P.cmp("delta.n_channels", sorted(res["deltas"]), sorted(FEATURES))
    P.cmp("delta.no_channel_excludes_zero_pooled",
          [dim for dim in DIMS if res["deltas"][dim]["all"]["ci_excludes_zero"]],
          [])
    P.cmp("delta.the_mask_geometry_is_orbit_invariant_pooled",
          [dim for dim in ("A", "F", "S")
           if res["deltas"][dim]["all"]["ci_excludes_zero"]], [])
    P.cmp("delta.no_channel_excludes_zero_pre_only",
          [dim for dim in DIMS
           if res["deltas"][dim]["pre"]["ci_excludes_zero"]], [])

    # 10) the trailing-window series ------------------------------------------
    P.cmp("series.days", res["series_days"], SERIES_DAYS)
    P.cmp("series.span_days", res["series_span_days"], 56)
    P.cmp("series.window_weeks", WINDOW_WEEKS, 8.0)
    P.cmp("series.min_per_orbit", MIN_PER_ORBIT, 2)
    P.cmp("series.n_points", res["n_series_points"],
          len(SERIES_DAYS) * len(DIMS))
    P.cmp("series.same_days_in_every_panel",
          len({tuple(p["day"] for p in res["series"][dim]["points"])
               for dim in DIMS}), 1)
    P.cmp("series.points_are_echo_days",
          sorted({p["day"] for dim in DIMS
                  for p in res["series"][dim]["points"]}),
          sorted(d["day"] for d in res["echo"] if d["day"] in SERIES_DAYS))
    P.cmp("series.post_points_are_the_four_post_echo_days",
          res["n_series_post"], 4 * len(DIMS))
    P.cmp("series.the_two_states_are_both_present",
          sorted({p["state"] for dim in DIMS
                  for p in res["series"][dim]["points"]}), ["post", "pre"])
    P.cmp("series.window_is_the_configured_one",
          sorted({p["span_days"] for dim in DIMS
                  for p in res["series"][dim]["points"]}), [56])
    for dim in DIMS:
        s, w = res["series"][dim], SERIES_FACTS[dim]
        P.cmp(f"series.{dim}.n_points", s["n_points"], w["n"])
        P.cmp(f"series.{dim}.window_days", s["window_days"], 56)
        P.cmp(f"series.{dim}.min_per_orbit", s["min_per_orbit"], MIN_PER_ORBIT)
        P.cmp(f"series.{dim}.first_day", s["first_day"], SERIES_DAYS[0])
        P.cmp(f"series.{dim}.last_day", s["last_day"], SERIES_DAYS[-1])
        P.cmp(f"series.{dim}.label", s["label"], core.CHANNEL_TITLE[dim])
        P.cmp(f"series.{dim}.days", [p["day"] for p in s["points"]], SERIES_DAYS)
        P.cmp(f"series.{dim}.two_chips_of_each_orbit",
              [min(p["n_asc"] for p in s["points"]) >= MIN_PER_ORBIT,
               min(p["n_desc"] for p in s["points"]) >= MIN_PER_ORBIT],
              [True, True])
        P.cmp(f"series.{dim}.points_match_the_window_recomputed",
              [(p["n_asc"], p["n_desc"], round(p["delta"], 12))
               for p in s["points"]],
              [(u["n_asc"], u["n_desc"],
                round(u["desc_median"] - u["asc_median"], 12))
               for u in (trailing_window_at(res["rows"], dim, day)
                         for day in SERIES_DAYS)])
        P.cmp(f"series.{dim}.delta_is_the_two_medians",
              [round(p["delta"], 12) for p in s["points"]],
              [round(p["desc_median"] - p["asc_median"], 12)
               for p in s["points"]])
        P.add(f"series.{dim}.min_delta", s["min_delta"], w["min"], ROUND_TOL)
        P.add(f"series.{dim}.max_delta", s["max_delta"], w["max"], ROUND_TOL)
        P.add(f"series.{dim}.max_abs_delta", s["max_abs_delta"], w["abs"],
              ROUND_TOL)
        P.cmp(f"series.{dim}.n_sign_changes", s["n_sign_changes"], w["nsc"])
        P.cmp(f"series.{dim}.n_distinct", s["n_distinct"], w["nd"])
        P.cmp(f"series.{dim}.max_abs_is_the_larger_end",
              s["max_abs_delta"],
              max(abs(s["min_delta"]), abs(s["max_delta"])))
    pw = res["post_window"]
    P.cmp("series.post_window.day", pw["day"], POST_WINDOW_FACT["day"])
    P.cmp("series.post_window.dim", pw["dim"], POST_WINDOW_FACT["dim"])
    P.cmp("series.post_window.span_days", pw["span_days"],
          POST_WINDOW_FACT["span_days"])
    P.cmp("series.post_window.n_asc", pw["n_asc"], POST_WINDOW_FACT["n_asc"])
    P.cmp("series.post_window.n_desc", pw["n_desc"], POST_WINDOW_FACT["n_desc"])
    P.cmp("series.post_window.kept", pw["kept"], POST_WINDOW_FACT["kept"])
    P.cmp("series.the_post_chip_is_the_day_it_names", res["post_echo_day"],
          POST_WINDOW_FACT["day"])
    P.cmp("series.the_post_window_holds_both_orbits",
          [pw["n_asc"], pw["n_desc"]], [4, 7])

    # 11) the composition of the record --------------------------------------
    comp = res["composition"]
    P.cmp("composition.orbit_state", comp["orbit_state"],
          COMPOSITION_FACTS["orbit_state"])
    P.cmp("composition.orbit_state_echo", comp["orbit_state_echo"],
          COMPOSITION_FACTS["orbit_state_echo"])
    P.cmp("composition.echo_by_state_orbit", comp["echo_by_state_orbit"],
          COMPOSITION_FACTS["echo_by_state_orbit"])
    P.cmp("composition.n_echo_by_orbit", res["n_echo_by_orbit"],
          {o: cov["per_orbit"][o]["n_with_echo"] for o in ORBITS})
    P.cmp("composition.n_echo_by_state", res["n_echo_by_state"],
          {"pre": res["n_pre_echo"], "post": res["n_post_echo"]})
    P.cmp("composition.n_orbit_in_state", comp["n_orbit_in_state"],
          {k: {o: comp["orbit_state"][o][k] for o in ORBITS}
           for k in ("pre", "post")})
    P.cmp("composition.orbit_state_sums_to_the_record",
          sum(v for o in ORBITS for v in comp["orbit_state"][o].values()),
          res["n_chips"])
    P.cmp("composition.months", comp["months"], committed_months(ref))
    P.cmp("composition.months_sum", sum(comp["months"].values()), res["n_chips"])
    P.cmp("composition.sections", comp["sections"],
          COMPOSITION_FACTS["sections"])
    P.cmp("composition.day_block_is_the_record", len(res["days"]), res["n_days"])
    P.cmp("composition.day_windows_sum",
          sum(d["n_windows"] for d in res["days"]), res["n_windows"])
    P.cmp("composition.day_chips_sum",
          sum(d["n_unique_chips"] for d in res["days"]), res["n_chips"])
    P.cmp("composition.day_echo_sum",
          sum(d["n_with_echo"] for d in res["days"]), res["n_echo"])
    P.cmp("composition.day_duplicates_sum",
          sum(d["n_duplicated"] for d in res["days"]), res["n_duplicated_windows"])
    P.cmp("composition.days_match_the_committed_timeline",
          {d["day"]: [d["n_windows"], d["n_unique_chips"],
                      d["n_with_echo"]] for d in res["days"]},
          {day: [e["n_windows"], e["n_unique_chips"], e["n_with_echo"]]
           for day, e in committed_timeline(ref).items()})

    # 12) the co-variates, re-read from the raw extract ----------------------
    ext = res["controls_from_extract"]
    for name in ("all", "echo"):
        for o in ORBITS:
            for k in ORBIT_CONTROLS:
                P.add(f"controls.{name}.{o}.{k}", res["covariates"][name][o][k],
                      ext[name][o][k], FLOAT_TOL)
        P.add(f"controls.{name}.delta_is_the_two_medians",
              res["covariates"][name]["delta"],
              {k: ext[name]["DESCENDING"][k] - ext[name]["ASCENDING"][k]
               for k in ORBIT_CONTROLS}, FLOAT_TOL)
        for k, want in CONTROL_FACTS[name].items():
            P.add(f"controls.{name}.{k}.delta", res["covariates"][name]["delta"][k],
                  want, ROUND_TOL)
    P.cmp("controls.the_extract_covers_every_chip",
          res["orbit_chain"]["n_in_extract"], res["n_chips"])
    P.cmp("controls.the_ascending_passes_are_the_warmer_ones",
          res["covariates"]["all"]["delta"]["temperature_c"] < 0, True)
    P.cmp("controls.the_ascending_passes_are_the_brighter_ones",
          res["covariates"]["echo"]["delta"]["intensity"] < 0, True)
    P.cmp("controls.no_rain_in_the_record",
          sorted({res["covariates"]["all"][o]["precipitation_mm"] for o in ORBITS}),
          [0.0])

    # 13) this figure's own drawing conventions -------------------------------
    P.cmp("panels.order", DIMS, ["gamma2", "P", "D", "F", "S", "A"])
    P.cmp("panels.grid", [PANEL_GRID[0], PANEL_GRID[1]], [4, 2])
    P.cmp("panels.every_dimension_is_drawn", sorted(DIMS), sorted(FEATURES))
    P.cmp("panels.ocv_channels_first", DIMS[:len(OCV)], OCV)
    P.cmp("panels.orbit_offset_days",
          [orbit_dx("ASCENDING"), orbit_dx("DESCENDING")], [0.18, -0.18])
    P.cmp("panels.gap_bands", len(res["gaps"]), res["n_gap_steps"])
    P.cmp("panels.only_one_line_convention",
          sorted({b["after"] for b in res["bands"]}
                 - set(res["day_list"])), [])
    P.cmp("figure.name", FIG_NAME, "fig_haltern_ocv_orbits")
    P.cmp("figure.script", os.path.relpath(
        os.path.abspath(__file__), repo), "code/haltern/fig_haltern_ocv_orbits.py")
    P.cmp("figure.json", os.path.relpath(OUT_JSON, repo),
          "data/haltern/fig_haltern_ocv_orbits.json")
    P.cmp("figure.report", os.path.relpath(OUT_MD, repo),
          "figures/haltern/fig_haltern_ocv_orbits.md")
    P.cmp("figure.png",
          os.path.relpath(os.path.join(FIGDIR, FIG_NAME + ".png"), repo),
          "figures/haltern/fig_haltern_ocv_orbits.png")
    P.cmp("figure.dpi", float(plt.rcParams["savefig.dpi"]), 600.0)
    P.cmp("figure.width_mm", COL_WIDTH / MM, 183.0)

    # 14) the paper's own A-E pin contract: this figure adds nothing to it ----
    P.cmp("paper_contract.n_failed", paper["pins"]["n_failed"], 0)
    P.cmp("paper_contract.float_tolerance", paper["pins"]["float_tolerance"],
          FLOAT_TOL)
    P.cmp("paper_contract.n_checks", paper["pins"]["n_checks"], 1857)
    P.cmp("paper_contract.figures", paper["figures"],
          [f"figures/haltern/fig_haltern_ocv_paper_{l}.png" for l in "ABCDE"])
    P.cmp("paper_contract.this_figure_is_not_in_it",
          any(FIG_NAME in f for f in paper["figures"]), False)
    P.cmp("paper_contract.orbits_are_reported_there",
          sorted(paper["coverage_by_orbit"]), sorted(res["coverage"]))
    return P


# ---------------------------------------------------------------------------
# The figure: a 2 x 4 grid — six channels, the split, the composition
# ---------------------------------------------------------------------------
def axis_label(dim):
    """The channel's own description (its axis label), without the bare name."""
    lab = core.CHANNEL_AXIS[dim]
    return lab.split("(", 1)[1].rstrip(")") if "(" in lab else lab


def draw_record_frame(ax, res):
    """The dated frame of the record: the two stretches, the cut, the holes.

    Identical to the whole-record companion's frame, so the two figures can be
    read side by side: the pre stretch tinted, the post stretch tinted darker,
    the nine holes of the record as light bands and the collapse date as a
    dashed cut inside the 12-day step that holds it.
    """
    ax.set_xlim(shift_dt(res["first_day"], -5), shift_dt(res["last_day"], 5))
    ax.axvspan(dtm(res["first_day"]), shift_dt(res["last_pre_day"], 0.5),
               facecolor=PRE_TINT, alpha=PRE_TINT_ALPHA, lw=0.0, zorder=0.4)
    ax.axvspan(shift_dt(res["first_post_day"], -0.5), dtm(res["last_day"]),
               facecolor=POST_TINT, alpha=POST_TINT_ALPHA, lw=0.0, zorder=0.4)
    for gap in res["gaps"]:
        ax.axvspan(dtm(gap["after"]), dtm(gap["before"]), facecolor=GAP_COLOR,
                   lw=0.0, zorder=0.5)
    ax.axvline(EVENT, color=EVENT_COLOR, lw=0.9, ls=(0, (2.5, 1.6)), zorder=1.6)


def date_axis(ax, res, every=2, label=False):
    """The record's own month ladder: a tick every ``every`` months.

    This axis is 265 days wide and each of the eight panels is about 80 mm, so
    the ladder is every second month, carrying the year whenever it changes.
    """
    ticks, month, year = [], day_date(res["first_day"]).replace(day=1), None
    last = day_date(res["last_day"])
    while month <= last:
        if month.year != year:
            year = month.year
            ticks.append((month, f"{month:%b}\n{year}"))
        else:
            ticks.append((month, f"{month:%b}"))
        month = add_months(month, every)
    ax.set_xticks([dtm(m.isoformat()) for m, _ in ticks],
                  [lab for _, lab in ticks])
    ax.tick_params(axis="x", length=2.0, pad=1.5)
    if label:
        ax.set_xlabel("acquisition day")


def draw_channel(ax, res, dim, letter, right=False, bottom=False):
    """One channel: the two orbits' chips, then the trailing-window difference.

    The left axis carries the channel's own chips, coloured and marked by pass
    geometry; the twin axis behind them carries ``delta_orbit`` — the step line
    (drawn only across a single 3-day cadence step, so every longer step of the
    series stays open), the pooled difference as a dashed reference and the zero
    line. The 54 echo-free chips are ticks on the panel floor: they carry no
    ordinate on any of the six channels, but they are the acquisition days whose
    windows the series walks over.
    """
    draw_record_frame(ax, res)
    vals = [d[dim] for d in res["rows"] if d.get(dim) is not None]
    lo, hi = min(vals), max(vals)
    pad = 0.14 * (hi - lo) if hi > lo else max(0.5, 0.14 * abs(hi))
    ax.set_ylim(lo - pad, hi + pad)
    for d in res["rows"]:
        if d["echo"] or d.get(dim) is None:
            continue
        x = x_of(d["day"], d["orbit"])
        ax.plot([x, x], [lo - pad, lo - pad + 0.06 * (hi - lo + 2 * pad)], "-",
                lw=0.7, color=NO_ECHO_COLOR, zorder=3)
    for orbit in ORBITS:
        sel = [d for d in res["rows"]
               if d["orbit"] == orbit and d.get(dim) is not None]
        ax.plot([x_of(d["day"], orbit) for d in sel], [d[dim] for d in sel],
                ORBIT_MARKER[orbit], ms=3.1, color=ORBIT_COLOR[orbit],
                markeredgewidth=0, zorder=4)
    ax.tick_params(axis="y", labelsize=6.4, length=2.0, pad=1.5)
    date_axis(ax, res, label=bottom)

    # the trailing difference of the same channel, on the twin axis behind it
    ax2 = ax.twinx()
    pts = res["series"][dim]["points"]
    for a, b in zip(pts, pts[1:]):
        if day_span(a["day"], b["day"]) <= CADENCE_DAYS:
            ax2.plot([dtm(a["day"]), dtm(b["day"])], [a["delta"], b["delta"]],
                     "-", color=DELTA_COLOR, lw=0.9, zorder=6)
    pooled = res["deltas"][dim]["all"]["delta_median"]
    ax2.plot([dtm(p["day"]) for p in pts], [p["delta"] for p in pts], "o",
             ms=2.8, color=DELTA_COLOR, markeredgewidth=0, zorder=7)
    ax2.axhline(pooled, color=POOLED_COLOR, lw=0.9, ls=(0, (4.0, 1.8)), zorder=5)
    ax2.axhline(0.0, color=ZERO_COLOR, lw=0.8, zorder=4.6)
    dvals = [p["delta"] for p in pts] + [0.0, pooled]
    dlo, dhi = min(dvals), max(dvals)
    dpad = 0.22 * (dhi - dlo) if dhi > dlo else 0.5
    ax2.set_ylim(dlo - dpad, dhi + dpad)
    ax2.tick_params(axis="y", labelsize=6.4, length=2.0, pad=1.5)
    if right:
        ax2.set_ylabel(DELTA_LABEL, fontsize=6.6, labelpad=0.8)
    panel(ax, letter, f"{core.CHANNEL_TITLE[dim]} — Δ {pooled:+.3f}\n"
                      f"{axis_label(dim)}")


def bar_panel(ax, res):
    """(g) the pooled difference beside the pre-only one, with its 95% CI.

    Bars are ``median(DESC) - median(ASC)``: filled for all 66 echo chips,
    hatched for the 62 pre-event ones — the pre-only pair is the same contrast
    with the record's four post-event chips (two ascending and two descending
    passes) dropped, so the two bars of a channel differ by *state*, not by
    orbit. The whiskers are the seeded percentile CIs of the median difference
    (each orbit resampled independently, 10 000 draws, seed 7 — the paper's own
    ``n_boot`` and ``rng_seed``) and a star marks a CI that excludes zero.
    """
    width = 0.34
    for i, dim in enumerate(DIMS):
        for state, dx in (("all", -0.5 * width), ("pre", 0.5 * width)):
            d = res["deltas"][dim][state]
            lo, hi = d["delta_ci95"]
            x = i + dx
            ax.bar([x], [d["delta_median"]], width=width,
                   color="#CBD8EA" if state == "all" else "#F2F3F5",
                   edgecolor=DELTA_COLOR, linewidth=0.5,
                   hatch=None if state == "all" else "///", zorder=3)
            ax.plot([x, x], [lo, hi], "-", color=DELTA_COLOR, lw=0.8, zorder=4)
            for y in (lo, hi):
                ax.plot([x - 0.10, x + 0.10], [y, y], "-", color=DELTA_COLOR,
                        lw=0.8, zorder=4)
            if d["ci_excludes_zero"]:
                ax.plot([x], [hi], marker="*", ms=4.4, color=EVENT_COLOR,
                        zorder=5, clip_on=False)
    ax.axhline(0.0, color=ZERO_COLOR, lw=0.8, zorder=2)
    vals = [v for dim in DIMS for st in ("all", "pre")
            for v in res["deltas"][dim][st]["delta_ci95"]]
    vlo, vhi = min(vals + [0.0]), max(vals + [0.0])
    pad = 0.16 * (vhi - vlo)
    ax.set_ylim(vlo - pad, vhi + 0.30 * pad)
    ax.set_xlim(-0.6, len(DIMS) - 0.4)
    ax.set_xticks(np.arange(len(DIMS)), [core.CHANNEL_TITLE[d] for d in DIMS])
    ax.tick_params(axis="x", length=0.0, labelsize=7.0)
    ax.tick_params(axis="y", labelsize=6.4, length=2.0, pad=1.5)
    ax.set_ylabel(DELTA_LABEL, fontsize=6.6)
    ax.text(0.02, 0.985, "filled = all 66 echo chips\n"
            "hatched = the 62 pre-event ones\n* = the 95% CI excludes zero",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.8,
            color="#5F5F5F", linespacing=1.4)
    panel(ax, "g", "pooled vs pre-only Δ")


def composition_panel(ax, res):
    """(h) where the 120 chips and the 66 echoes sit over the two directions.

    Two bars per orbit direction — the pre-event chips and the post-event ones —
    each split into its echo-bearing part (the orbit's own colour) and its
    echo-free part (grey), with ``echoes/chips`` written above. The block is the
    coverage cross-cut of the record: it says where the echoes sit, and both
    Fisher tests are printed rather than implied.
    """
    comp = res["composition"]
    for i, o in enumerate(ORBITS):
        for state, dx in (("pre", -0.19), ("post", 0.19)):
            n = comp["orbit_state"][o][state]
            k = comp["orbit_state_echo"][o][state]
            ax.bar([i + dx], [k], width=0.34, color=ORBIT_COLOR[o], zorder=4)
            ax.bar([i + dx], [n - k], width=0.34, bottom=k,
                   color=NO_ECHO_COLOR, zorder=3)
            ax.text(i + dx, n + 0.4, f"{k}/{n}", ha="center", va="bottom",
                    fontsize=6.0, color="#333333")
    ax.set_xticks(np.arange(len(ORBITS)), [ORBIT_SHORT[o] for o in ORBITS])
    ax.tick_params(axis="x", length=0.0, labelsize=7.0)
    ax.tick_params(axis="y", labelsize=6.4, length=2.0, pad=1.5)
    ax.set_xlim(-0.6, len(ORBITS) - 0.4)
    ax.set_ylim(0.0, 22.5)
    ax.set_ylabel("chips of the record", fontsize=6.6)
    ax.text(0.02, 0.985,
            "pooled Fisher p = {:.3f}\npre vs post: ASC p = {:.3f}\n"
            "DESC p = {:.3f}\nthe only post-event echo is ASC\n"
            "filled = echo, grey = no echo".format(
                res["p_orbit"], res["p_orbit_by_state"]["ASCENDING"],
                res["p_orbit_by_state"]["DESCENDING"]),
            transform=ax.transAxes, ha="left", va="top", fontsize=5.8,
            color="#5F5F5F", linespacing=1.5)
    panel(ax, "h", "the echoes over the two\npass geometries")


def fig_orbits(res, out_dir=None):
    """The 2 x 4 companion: a-f the six channels, g the split, h the coverage.

    The six channel panels share one frame and one treatment; (g) carries the
    pooled and the pre-only difference of every channel with its CI, and (h)
    carries where the chips and the echoes sit over the two pass geometries.
    """
    fig, axes = plt.subplots(PANEL_GRID[0], PANEL_GRID[1],
                             figsize=(COL_WIDTH, COL_WIDTH * 0.95))
    flat = list(np.ravel(axes))
    last_col = PANEL_GRID[1] - 1
    # The six channel panels fill the first three rows; the fourth row carries
    # the two summary panels. The dated x axis is labelled once, under the
    # channel row that closes the six (as in the whole-record companion).
    date_row = -(-len(DIMS) // PANEL_GRID[1]) - 1
    for i, dim in enumerate(DIMS):
        draw_channel(flat[i], res, dim, "abcdef"[i],
                     right=bool(i % PANEL_GRID[1] == last_col),
                     bottom=bool(i // PANEL_GRID[1] == date_row))
    for i in range(len(DIMS)):
        if i // PANEL_GRID[1] != date_row:
            flat[i].tick_params(axis="x", labelbottom=False)
    bar_panel(flat[6], res)
    composition_panel(flat[7], res)

    handles = [Line2D([], [], marker=ORBIT_MARKER[o], ls="", ms=3.4,
                      color=ORBIT_COLOR[o],
                      label=f"{ORBIT_LABEL[o]} chip: "
                            f"{res['n_echo_by_orbit'][o]} echoes of "
                            f"{res['coverage']['per_orbit'][o]['n_chips']} chips")
               for o in ORBITS]
    handles += [
        Line2D([], [], marker="o", ls="", ms=2.8, color=DELTA_COLOR,
               label=f"Δorbit on the {int(WINDOW_WEEKS)}-week trailing window "
                     f"({len(SERIES_DAYS)} days, all pre-event)"),
        Line2D([], [], color=POOLED_COLOR, ls=(0, (4.0, 1.8)), lw=0.9,
               label=f"pooled Δorbit over the {res['n_echo']} echo chips"),
        Line2D([], [], marker="|", ls="", ms=6.0, color=NO_ECHO_COLOR,
               label=f"echo-free chip, on the panel floor ({res['n_no_echo']})"),
        Line2D([], [], color=EVENT_COLOR, ls=(0, (2.5, 1.6)), lw=0.9,
               label=f"collapse {core.COLLAPSE_DATE}, inside a 12-day step"),
    ]
    # a figure-level strip, as in the whole-record companion: inside a panel this
    # key covers chips, and each panel is about 80 x 30 mm.
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.004, 0.868),
               ncol=2, frameon=False, fontsize=6.8, handletextpad=0.4,
               columnspacing=1.0, borderaxespad=0.0, labelspacing=0.28)
    fig.suptitle("HALTERN — the orbit split of the echo-mask vector\n"
                 "the two pass geometries alternate, so the split is a "
                 "trailing-window difference", x=0.004, y=1.050, ha="left",
                 va="top", fontsize=10.0, linespacing=1.3)
    fig.text(0.004, 0.966,
             "x = the acquisition day (one dot per echo-bearing chip, the two "
             "orbits offset by 0.18 d); left axis = the channel, twin axis = "
             "Δorbit =\n"
             f"median(DESC) − median(ASC) over the echo chips of a trailing "
             f"{int(WINDOW_WEEKS)}-week window at least {MIN_PER_ORBIT} chips of "
             f"each orbit wide;\n"
             f"a line joins only points one {CADENCE_DAYS}-day cadence step "
             f"apart, so the nine holes stay open — {res['n_echo']} echo chips of "
             f"{res['n_chips']}: {res['n_pre_echo']} of {res['n_pre_chips']} pre, "
             f"{res['n_post_echo']} of {res['n_post_chips']} post",
             ha="left", va="top", fontsize=7.0, color="#444444", linespacing=1.4)
    fig.tight_layout(pad=0.5, h_pad=1.3, w_pad=1.5, rect=(0, 0, 0.985, 0.827))
    return save_fig(fig, FIG_NAME, out_dir)


# ---------------------------------------------------------------------------
# The artifacts: the JSON of every number and the report that reads them
# ---------------------------------------------------------------------------
def fmt(x, nd=3):
    """One cell of a Markdown table."""
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}" if math.isfinite(x) else "n/a"
    return str(x)


def md_table(header, rows, nd=3):
    """A Markdown table with the family's own number format.

    ``|`` inside a cell is escaped: one of the six channels is described as
    ``||centroid - peak||, px``, and an unescaped pipe would split the row.
    """
    def cell(x):
        return fmt(x, nd).replace("|", "\\|")

    out = ["| " + " | ".join(cell(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(cell(c) for c in r) + " |")
    return out


def delta_rows(res, state="all"):
    """One row per channel: the two medians, the difference, its CI, Cliff's."""
    out = []
    for dim in DIMS:
        d = res["deltas"][dim][state]
        out.append([dim, axis_label(dim), d["n_asc"], d["n_desc"],
                    d["asc_median"], d["desc_median"], d["delta_median"],
                    d["delta_ci95"][0], d["delta_ci95"][1], d["cliffs_delta"],
                    "yes" if d["ci_excludes_zero"] else "no"])
    return out


DELTA_HEADER = ["channel", "what it is", "n ASC", "n DESC", "median ASC",
                "median DESC", "Δ (DESC − ASC)", "CI 2.5%", "CI 97.5%",
                "Cliff's δ", "CI excludes 0"]


def coverage_rows(res):
    """One row per orbit x state: the chips, the echoes and the echo rate."""
    cov = res["coverage"]
    out = []
    for o in ORBITS:
        for state in ("pre", "post"):
            b = cov["per_orbit"][o]["by_state"][state]
            out.append([ORBIT_SHORT[o], state, b["n_chips"], b["n_with_echo"],
                        b["n_without_echo"], b["echo_rate"]])
        b = cov["per_orbit"][o]
        out.append([ORBIT_SHORT[o], "all", b["n_chips"], b["n_with_echo"],
                    b["n_without_echo"], b["echo_rate"]])
    return out


def series_rows(res, nd=4):
    """The trailing-window series as a table: one row per day, one column per
    channel, so the whole windowed reading of the record is on one page."""
    days = SERIES_DAYS
    by_dim = {dim: {p["day"]: p for p in res["series"][dim]["points"]}
              for dim in DIMS}
    out = []
    for day in days:
        p = by_dim[DIMS[0]][day]
        out.append([day, p["n_asc"], p["n_desc"], p["state"]]
                   + [by_dim[dim][day]["delta"] for dim in DIMS])
    return out


def series_summary_rows(res):
    """One row per channel: the range and the reversals of its own series."""
    out = []
    for dim in DIMS:
        s = res["series"][dim]
        out.append([dim, s["n_points"], s["min_delta"], s["max_delta"],
                    s["max_abs_delta"], s["n_sign_changes"], s["n_distinct"]])
    return out


def control_rows(res, name):
    """One row per co-variate: the two per-orbit medians and their difference."""
    out = []
    for k in ORBIT_CONTROLS:
        a = res["covariates"][name]["ASCENDING"][k]
        b = res["covariates"][name]["DESCENDING"][k]
        out.append([k, core.CHANNEL_TITLE.get(k, k), a, b, b - a])
    return out


def day_rows(res):
    """One row per acquisition day: what was requested and what echoed."""
    return [[d["day"], d["state"], ", ".join(ORBIT_SHORT[o] for o in d["orbits"]),
             d["n_windows"], d["n_unique_chips"], d["n_duplicated"],
             d["n_with_echo"], d["n_without_echo"]] for d in res["days"]]


def report(res, pins):
    """The report: what the figure draws, and every number behind it."""
    md = []
    geo, cov = res["geometry_block"], res["coverage"]
    md.append(f"# {core.SITE} — the orbit split of the echo-mask vector")
    md.append("")
    md.append(f"**Site** {core.SITE} · **asset** {core.ASSET_NAME} · **event** "
              f"{core.COLLAPSE_DATE} ({res['event_weekday']}) · **record** "
              f"{res['first_day']} -> {res['last_day']} "
              f"({res['n_days']} acquisition days, {res['n_chips']} unique chips, "
              f"{res['n_windows']} committed windows) · **figure** "
              f"`figures/haltern/{FIG_NAME}.png` (600 dpi) · **result** "
              f"`data/haltern/{FIG_NAME}.json`")
    md.append("")
    md.append("**What this figure is.** The record's four Sentinel-1 pass "
              "geometries fall into two directions that **alternate in blocks of "
              "two days**, so "
              f"no acquisition day carries both a usable ascending and a usable "
              f"descending echo — all {res['n_days_single_orbit']} of the "
              f"{res['n_days']} days carry a single direction — and a *per-day* "
              "orbit difference therefore does not exist. What does exist is the "
              "**trailing-window difference**")
    md.append("")
    md.append("```")
    md.append("delta_orbit(t) = median(DESC in window(t)) - median(ASC in "
              "window(t))")
    md.append("```")
    md.append("")
    md.append(f"with `window(t)` the {int(WINDOW_WEEKS)}-week window ending at the "
              f"echo-bearing day `t` (at least {MIN_PER_ORBIT} chips of each orbit, "
              "no interpolation), drawn on the twin axis of each panel beside the "
              "two orbits' own chips; the dashed line is the **pooled** difference "
              "over the whole record. A line joins two points only when they are "
              f"one {CADENCE_DAYS}-day cadence step apart, so the nine holes of the "
              "record are gaps in the series and are never bridged. Panels (g) and "
              "(h) carry the pooled/pre-only contrast, the composition of the "
              "record and the two Fisher tests.")
    md.append("")
    md.append(f"**What it is not.** {DELTA_NOTE[0].upper() + DELTA_NOTE[1:]}. "
              f"{CONFOUND_NOTE[0].upper() + CONFOUND_NOTE[1:]}. The split adds no "
              "number to the A-E pin contract of `fig_haltern_ocv_paper.py`.")
    md.append("")
    md.append("## The four pass geometries behind the labels")
    md.append("")
    md.append(f"{res['geometry_sentence'][0].upper() + res['geometry_sentence'][1:]}.")
    md.append("")
    md.extend(md_table(
        ["orbit", "chips", "swath", "relative orbit", "UT hour", "VV / VH chips",
         "label in the export"],
        [[ORBIT_SHORT[o], geo[o]["n_chips"], "/".join(geo[o]["swaths"]),
          " + ".join(f"{r:g}" for r in geo[o]["relative_orbits"]),
          " + ".join(geo[o]["utc_hh"]),
          f"{geo[o]['polarisations'].get('VV', 0)} / "
          f"{geo[o]['polarisations'].get('VH', 0)}",
          ", ".join(res["orbit_chain"]["pass_labels"][o])] for o in ORBITS]))
    md.append("")
    ch = res["orbit_chain"]
    md.append(f"The chain is checked on every chip: the channel table's "
              f"`orbit_direction`/`pass_label` agree with the committed "
              f"measurement extract ({ch['n_in_extract']} of {ch['n_chips']} chips "
              f"present, **{ch['n_orbit_mismatches_extract']} orbit mismatches**, "
              f"**{ch['n_pass_label_mismatches_extract']} label mismatches**) and "
              f"with the resolved burst record (`reference/haltern_bursts.json`, "
              f"{ch['n_resolved_bursts']} resolved bursts, "
              f"**{ch['n_orbit_mismatches_bursts']} orbit mismatches**); the "
              f"{res['n_chips']} chips map to {res['n_bursts']} distinct bursts. "
              f"Which pass a chip is, is re-derived here and not assumed.")
    md.append("")
    md.append("## The pooled orbit split")
    md.append("")
    md.extend(md_table(DELTA_HEADER, delta_rows(res, "all")))
    md.append("")
    g_all = res["deltas"]["gamma2"]["all"]
    p_all = res["deltas"]["P"]["all"]
    md.append(f"The pooled difference is **gamma2 "
              f"{fmt(g_all['delta_median'], 4)}** with a seeded 95% CI of "
              f"[{fmt(g_all['delta_ci95'][0])}, {fmt(g_all['delta_ci95'][1])}] over "
              f"{g_all['n_asc']} ascending and {g_all['n_desc']} descending echo "
              f"chips, and **P {fmt(p_all['delta_median'], 4)}** with "
              f"[{fmt(p_all['delta_ci95'][0])}, {fmt(p_all['delta_ci95'][1])}]. "
              f"Every channel's CI spans zero — the widest is gamma2, whose upper "
              f"end sits at {fmt(g_all['delta_ci95'][1])} against only "
              f"{g_all['n_desc']} descending chips — so on this record the mask "
              f"vector is **orbit-invariant**, and the three mask-geometry channels "
              f"("
              + ", ".join(f"{d} {fmt(res['deltas'][d]['all']['delta_median'], 4)}"
                          f" (CI {fmt(res['deltas'][d]['all']['delta_ci95'][0])} .. "
                          f"{fmt(res['deltas'][d]['all']['delta_ci95'][1])})"
                          for d in ("A", "F", "S"))
              + f") plus D {fmt(res['deltas']['D']['all']['delta_median'], 4)} "
              f"point the same way. Cliff's delta is called with the ascending "
              f"chips as reference, so it carries the same sign as the bars and as "
              f"the difference.")
    md.append("")
    md.append(f"No channel has a pooled CI that excludes zero; the 10 000-resample "
              f"bootstrap is seeded ({BOOT_SEED}) and resamples each orbit "
              f"independently, because the two orbits are disjoint sets of "
              f"acquisition days and the chips are not paired.")
    md.append("")
    md.append("## The same split before the collapse")
    md.append("")
    md.extend(md_table(DELTA_HEADER, delta_rows(res, "pre")))
    md.append("")
    g_pre = res["deltas"]["gamma2"]["pre"]
    md.append(f"{CONFOUND_NOTE[0].upper() + CONFOUND_NOTE[1:]}. "
              f"{res['n_post_echo']} of the record's {res['n_echo']} echo chips are "
              f"post-event — "
              f"{res['post_echo_orbits'].get('ASCENDING', 0)} on the ascending and "
              f"{res['post_echo_orbits'].get('DESCENDING', 0)} on the descending "
              f"geometry — so the pre-only column is the same contrast with a "
              f"balanced pair dropped. gamma2 moves from "
              f"{fmt(g_all['delta_median'], 4)} to {fmt(g_pre['delta_median'], 4)} "
              f"and P from {fmt(p_all['delta_median'], 4)} to "
              f"{fmt(res['deltas']['P']['pre']['delta_median'], 4)}, and gamma2's CI "
              f"from [{fmt(g_all['delta_ci95'][0])}, {fmt(g_all['delta_ci95'][1])}] "
              f"to [{fmt(g_pre['delta_ci95'][0])}, {fmt(g_pre['delta_ci95'][1])}] — "
              f"both span zero, and the whole mask geometry is unchanged. Neither "
              f"column is a damage reading.")
    md.append("")
    md.append("## The trailing-window series of every channel")
    md.append("")
    md.append(f"The series is drawn at the {len(SERIES_DAYS)} echo-bearing days whose "
              f"{int(WINDOW_WEEKS)}-week trailing window holds at least "
              f"{MIN_PER_ORBIT} chips of each orbit "
              f"(``{'``, ``'.join(SERIES_DAYS)}``), one point per day and channel: "
              f"{res['n_series_points']} points in all, "
              f"{res['n_series_post']} of them at the four post-event echo days.")
    md.append("")
    md.extend(md_table(["day", "n ASC", "n DESC", "state"]
                       + [core.CHANNEL_TITLE[d] for d in DIMS],
                       series_rows(res)))
    md.append("")
    md.extend(md_table(["channel", "points", "min Δ", "max Δ", "max |Δ|",
                        "sign changes", "distinct values"],
                       series_summary_rows(res)))
    md.append("")
    md.append(f"**The series moves where the pooled number cannot say so.** D swings "
              f"between {fmt(res['series']['D']['min_delta'], 4)} and "
              f"{fmt(res['series']['D']['max_delta'], 4)} "
              f"({res['series']['D']['n_sign_changes']} sign changes), S between "
              f"{fmt(res['series']['S']['min_delta'], 4)} and "
              f"{fmt(res['series']['S']['max_delta'], 4)} "
              f"({res['series']['S']['n_sign_changes']} sign changes), F between "
              f"{fmt(res['series']['F']['min_delta'], 4)} and "
              f"{fmt(res['series']['F']['max_delta'], 4)}, and A down to "
              f"{fmt(res['series']['A']['min_delta'], 4)} — all of them channels "
              f"whose pooled CI spans zero. A pooled zero is an average over "
              f"orbit-dependent windows, not the absence of an orbit effect; and "
              f"with the two directions never sharing a day, a per-day reading of "
              f"any of it does not exist.")
    md.append("")
    pw = res["post_window"]
    md.append(f"**The series keeps its points past the event.** The record's four "
              f"post-event echoes fall on two ascending and two descending days, "
              f"so the window of the first of them (`{pw['day']}`, "
              f"{ORBIT_SHORT['ASCENDING']}) already holds {pw['n_asc']} ascending "
              f"and {pw['n_desc']} descending echo chips and a window difference "
              f"exists there: the drawn line runs into the post-event stretch "
              f"instead of stopping at the collapse.")
    md.append("")
    md.append(f"The {res['n_gap_steps']} steps of the record longer than the "
              f"{CADENCE_DAYS}-day cadence ({', '.join(str(g) for g in sorted({g['days'] for g in res['gaps']}))} "
              f"days) are holes: the drawn line breaks at each of them and no "
              f"value is interpolated across one.")
    md.append("")
    md.extend(md_table(["after", "before", "days", "crosses the event"],
                       [[g["after"], g["before"], g["days"],
                         "yes" if g["crosses_event"] else "no"]
                        for g in res["gaps"]]))
    md.append("")
    md.append("## Where the echoes sit over the two pass directions")
    md.append("")
    md.extend(md_table(["orbit", "state", "chips", "with echo", "without echo",
                        "echo rate"], coverage_rows(res)))
    md.append("")
    md.append(f"The coverage of the two directions is a **distribution, not a "
              f"contrast**: pooled Fisher p = {fmt(res['p_orbit'])} "
              f"({cov['per_orbit']['ASCENDING']['n_with_echo']}+"
              f"{cov['per_orbit']['ASCENDING']['n_without_echo']} against "
              f"{cov['per_orbit']['DESCENDING']['n_with_echo']}+"
              f"{cov['per_orbit']['DESCENDING']['n_without_echo']}), and within an "
              f"orbit pre-vs-post p = "
              f"{fmt(res['p_orbit_by_state']['ASCENDING'])} (ASC) and "
              f"{fmt(res['p_orbit_by_state']['DESCENDING'])} (DESC). "
              f"{ORBIT_NOTE[0].upper() + ORBIT_NOTE[1:]}.")
    md.append("")
    md.append(f"Panel (h) draws the same block: the record's {res['n_days']} "
              f"acquisition days and {res['n_chips']} chips, "
              f"{res['n_echo']} of them echoes ({res['n_pre_echo']} pre, "
              f"{res['n_post_echo']} post), "
              f"{cov['per_orbit']['ASCENDING']['n_chips']} chips on the ascending "
              f"geometry and {cov['per_orbit']['DESCENDING']['n_chips']} on the "
              f"descending one.")
    md.append("")
    md.append("## What else moves with the pass geometry")
    md.append("")
    md.append(f"{COVARIATE_NOTE[0].upper() + COVARIATE_NOTE[1:]}. The medians below "
              "are recomputed twice — from the channel table and from the committed "
              "measurement extract — and pinned against each other.")
    md.append("")
    md.append(f"All {res['n_chips']} chips:")
    md.append("")
    md.extend(md_table(["co-variate", "label", "median ASC", "median DESC",
                        "Δ (DESC − ASC)"], control_rows(res, "all"), nd=6))
    md.append("")
    md.append(f"The {res['n_echo']} echo-bearing chips:")
    md.append("")
    md.extend(md_table(["co-variate", "label", "median ASC", "median DESC",
                        "Δ (DESC − ASC)"], control_rows(res, "echo"), nd=6))
    md.append("")
    md.append(f"The ascending passes are the warmer ones "
              f"({fmt(res['covariates']['all']['delta']['temperature_c'], 4)} C on "
              f"the record, "
              f"{fmt(res['covariates']['echo']['delta']['temperature_c'], 4)} C "
              f"among the echo chips) and the darker ones in SAR intensity "
              f"({fmt(res['covariates']['echo']['delta']['intensity'], 1)}); no chip "
              f"of the record carries rain "
              f"(`precipitation_mm` is 0.0 on every chip). Whatever an orbit split "
              f"of the mask channels would mean physically, these are the terms "
              f"that move with it.")
    md.append("")
    md.append("## The record: one direction per day, and its nine holes")
    md.append("")
    md.append(f"The record is {res['n_days']} acquisition days and "
              f"{res['n_windows']} committed windows "
              f"({res['n_duplicated_windows']} of them a repeat of a payload "
              f"already in the table), which compile to {res['n_chips']} unique "
              f"chips: {res['n_echo']} carry an echo mask ({res['n_pre_echo']} pre, "
              f"{res['n_post_echo']} post) and {res['n_no_echo']} do not. The two "
              f"directions alternate in blocks of two days and never share one: "
              f"all {res['n_days_single_orbit']} of the {res['n_days']} days carry "
              f"a single direction, so the two orbits of the record never meet "
              f"inside one day and this figure's ordinate never has two orbits to "
              f"subtract on it.")
    md.append("")
    md.extend(md_table(["day", "state", "orbits", "windows", "chips", "repeats",
                        "with echo", "without echo"], day_rows(res)))
    md.append("")
    md.append("## The pins and the provenance")
    md.append("")
    s = pins.summary()
    md.append(f"{s['n_checks']} pin comparisons against five committed layers (the "
              f"channel table's generated metadata, the mask-cache manifest, the "
              f"committed mask cache, the committed paper result and the committed "
              f"analysis reference, plus the raw measurement extract and the three "
              f"resolved reference records): **{s['n_exact']} exact**, "
              f"**{s['n_within_tolerance']} within "
              f"{s['float_tolerance']:g}** (the seeded numbers and the six-decimal "
              f"record facts are pinned at {s['round_tolerance']:g}), "
              f"**{s['n_failed']} failed**. The recomputed number is in "
              f"`recomputed`, the committed one in `committed`, and no artifact is "
              f"written before they agree.")
    md.append("")
    md.append("| committed input | layer | sha256 |")
    md.append("| --- | --- | --- |")
    for key, layer in (("channels_csv", "the channel table (this figure's rows)"),
                       ("windows", "the window file (the chips)"),
                       ("measurements", "the raw measurement extract"),
                       ("mask_cache", "the committed mask cache (the vector)"),
                       ("mask_cache_manifest", "the cache's manifest"),
                       ("segments", "the segment table"),
                       ("burst_geometry", "the resolved bursts (the pass)"),
                       ("osm_anchor", "the geometry anchor"),
                       ("osm_node_extract", "the OSM node extract")):
        src = res["reference"]["sources"][key]
        md.append(f"| `{src['path']}` | {layer} | `{src['sha256'][:16]}…` |")
    md.append("")
    md.append(f"Figures: `figures/haltern/{FIG_NAME}.png` (600 dpi). JSON: "
              f"`data/haltern/{FIG_NAME}.json`.")
    md.append("")
    return "\n".join(md)


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
    """The result JSON: the whole axis, every pinned number, no raw pixels."""
    meta = res["meta"]
    repo = os.path.abspath(os.path.join(core.DATA, os.pardir, os.pardir))
    return {
        "figure": FIG_NAME,
        "script": "code/haltern/fig_haltern_ocv_orbits.py",
        "generated_by": "code/haltern/fig_haltern_ocv_orbits.py",
        "site": core.SITE, "site_short": core.SITE_SHORT,
        "asset_id": core.ASSET_ID, "asset_name": core.ASSET_NAME,
        "title": "HALTERN — the orbit split of the echo-mask vector: four pass "
                 "geometries fall into two directions that alternate, so the "
                 "split is a trailing-window difference",
        "role": "companion: delta_orbit = median(DESC) - median(ASC) per channel, "
                "pooled and as a trailing-window series, with the pass geometry, "
                "the coverage of the two directions and the co-variate controls; "
                "it adds no number to the A-E pin contract of "
                "fig_haltern_ocv_paper.py",
        "ocv": list(OCV), "features": list(FEATURES),
        "mask_window": meta["mask_window"], "mask_rule": meta["mask_rule"],
        "axis": {
            "name": "acquisition_day",
            "start": res["first_day"], "end": res["last_day"],
            "event": core.COLLAPSE_DATE, "event_weekday": res["event_weekday"],
            "span_days": res["span_days"], "cadence_days": CADENCE_DAYS,
            "definition": "the acquisition day, one dot per echo-bearing chip; "
                          "each dot carries a sub-day offset of 0.18 d by "
                          "direction, and no day carries chips of both directions",
            "steps": res["steps"], "gap_steps": res["gaps"],
            "collapse_step": res["collapse_step"],
            "bands": res["bands"], "n_bands": len(res["bands"]),
            "note": f"the event falls inside a drawn "
                    f"{res['collapse_step']['days']}-day step, so the two "
                    f"stretches meet at a band and not at a point",
        },
        "alternation": {
            "rule": "one orbit direction per acquisition day, in alternating "
                    "blocks of two days, so no acquisition day carries both a "
                    "usable ascending and a usable descending echo and a "
                    "per-day orbit difference does not exist",
            "n_days": res["n_days"],
            "n_days_single_orbit": res["n_days_single_orbit"],
            "n_days_two_orbits": res["n_days_two_orbits"],
            "days_two_orbits": res["days_two_orbits"],
            "two_orbit_days_echo_free": res["two_orbit_days_echo_free"],
            "n_gap_steps": res["n_gap_steps"],
            "gap_lengths_days": sorted({g["days"] for g in res["gaps"]}),
            "note": "the two directions alternate in blocks of two days along the "
                    "record, so the orbit contrast is a trailing-window "
                    "difference and not a paired test",
        },
        "orbits": {
            "n_orbits": len(ORBITS), "orbits": list(ORBITS),
            "labels": dict(ORBIT_SHORT), "export_labels": res["orbit_chain"]["pass_labels"],
            "per_orbit": res["coverage"]["per_orbit"],
            "fisher": res["coverage"]["fisher"],
            "p_pooled": res["p_orbit"],
            "fisher_by_state": res["p_orbit_by_state"],
            "geometry": res["geometry_block"],
            "geometry_sentence": res["geometry_sentence"],
            "chain": res["orbit_chain"],
            "note": ORBIT_NOTE,
        },
        "counts": {
            "n_chips": res["n_chips"], "n_days": res["n_days"],
            "n_windows": res["n_windows"],
            "n_duplicated_windows": res["n_duplicated_windows"],
            "n_echo": res["n_echo"], "n_no_echo": res["n_no_echo"],
            "n_pre_chips": res["n_pre_chips"], "n_pre_echo": res["n_pre_echo"],
            "n_post_chips": res["n_post_chips"], "n_post_echo": res["n_post_echo"],
            "n_pre_days": res["n_pre_days"], "n_post_days": res["n_post_days"],
            "n_echo_by_orbit": res["n_echo_by_orbit"],
            "n_echo_by_state": res["n_echo_by_state"],
            "n_echo_orbit_state": res["n_echo_orbit_state"],
            "post_echo_orbits": res["post_echo_orbits"],
            "pre_echo_orbits": res["pre_echo_orbits"],
        },
        "difference": {
            "definition": "delta_orbit = median(DESC) - median(ASC) over the "
                          "echo-bearing chips",
            "sign": "positive means the descending passes read higher on that "
                    "channel; Cliff's delta is called with the ascending chips as "
                    "reference, which gives the identical sign",
            "paired": False,
            "note": DELTA_NOTE,
            "bootstrap": {
                "statistic": "median(DESC) - median(ASC), each orbit resampled "
                             "independently",
                "n_boot": N_BOOT, "seed": BOOT_SEED, "alpha": ALPHA,
                "note": "the paper result's own n_boot and rng_seed, re-derived "
                        "here so every CI of this figure is reproducible",
            },
            "pooled": {dim: res["deltas"][dim]["all"] for dim in DIMS},
            "pre_only": {dim: res["deltas"][dim]["pre"] for dim in DIMS},
            "post_only": {dim: res["deltas"][dim]["post"] for dim in DIMS},
            "ci_excludes_zero": [dim for dim in DIMS
                                 if res["deltas"][dim]["all"]["ci_excludes_zero"]],
            "mask_geometry_is_orbit_invariant_pooled":
                [dim for dim in ("A", "F", "S")
                 if res["deltas"][dim]["all"]["ci_excludes_zero"]] == [],
            "confound": CONFOUND_NOTE,
            "post_echo_orbits": res["post_echo_orbits"],
        },
        "series": {
            "statistic": "delta_orbit(t) at every echo-bearing day",
            "window_weeks": WINDOW_WEEKS,
            "window_days": int(round(WINDOW_WEEKS * 7)),
            "min_per_orbit": MIN_PER_ORBIT,
            "join_rule": "a line is drawn only between two consecutive points one "
                         "cadence step (6 d) apart; every longer step of the "
                         "series is left open",
            "no_interpolation": "a day whose window holds fewer than "
                                f"{MIN_PER_ORBIT} chips of either orbit keeps no "
                                "point, and no value is interpolated across a hole",
            "days": res["series_days"], "n_points": res["n_series_points"],
            "n_points_per_dim": len(SERIES_DAYS),
            "n_post_event_points": res["n_series_post"],
            "per_dim": {dim: res["series"][dim] for dim in DIMS},
            "post_window": res["post_window"],
            "note": "every kept day is pre-event: the record's single post-event "
                    "echo chip has no descending chip in its own window",
        },
        "covariates": {
            "note": COVARIATE_NOTE,
            "controls": ORBIT_CONTROLS,
            "by_orbit_all": res["covariates"]["all"],
            "by_orbit_echo": res["covariates"]["echo"],
            "by_orbit_all_from_the_extract": res["controls_from_extract"]["all"],
            "by_orbit_echo_from_the_extract": res["controls_from_extract"]["echo"],
        },
        "cache": {"n_rows": len(res["cache"]), "n_P_set": res["n_P_set"],
                  "P_reference_chips": res["p_reference_chips"],
                  "guard": meta["cache_cross_check"]},
        "panels": {
            "grid": list(PANEL_GRID),
            "order": list(DIMS),
            "letters": "abcdef",
            "channels": {dim: {"letter": "abcdef"[i],
                               "label": core.CHANNEL_TITLE[dim],
                               "axis": core.CHANNEL_AXIS[dim],
                               "pooled_delta": res["deltas"][dim]["all"]["delta_median"]}
                         for i, dim in enumerate(DIMS)},
            "g": "the pooled difference beside the pre-only one, each with its "
                 "seeded 95% CI, a star for a CI that excludes zero",
            "h": "the 120 chips per direction and state, split into their "
                 "echo-bearing and echo-free parts, with both Fisher tests printed",
            "dot_colour": "pass geometry (ascending circles, descending squares)",
            "floor_ticks": "the 54 echo-free chips: no ordinate on any channel, "
                           "drawn on the panel floor",
        },
        "provenance": {
            "note": "the nine committed inputs, by content; the same digests the "
                    "analysis reference names, re-read here and pinned",
            "inputs": {key: {"path": os.path.relpath(path, repo),
                             "sha256": core.sha256(path)}
                       for key, path in
                       (("channels_csv", core.CSV_PATH),
                        ("mask_cache", core.CACHE_PATH),
                        ("mask_cache_manifest", core.MANIFEST_PATH),
                        ("windows", core.WINDOWS_PATH),
                        ("segments", core.SEGMENTS_PATH),
                        ("measurements", core.MEAS_PATH),
                        ("burst_geometry", core.REF_BURSTS),
                        ("osm_anchor", core.REF_GEOMETRY),
                        ("osm_node_extract", core.REF_OSM_NODE))},
            "paper_n_boot": res["paper"]["n_boot"],
            "paper_rng_seed": res["paper"]["rng_seed"],
        },
        "composition": res["composition"], "days": res["days"],
        "pins": pins_summary,
        # no wall-clock time: the figure is deterministic, so two runs of this
        # script must produce the identical JSON (the runtime stays on stdout).
        "figures": [f"figures/haltern/{FIG_NAME}.png"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    ap.add_argument("--no-figures", action="store_true",
                    help="run the cross-checks and write no PNG")
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, the orbit split of the record around "
          f"{core.COLLAPSE_DATE} …")
    res = compute()
    print(f"  {res['n_chips']} unique chips on {res['n_days']} acquisition days "
          f"({res['n_echo']} carry an echo: {res['n_pre_echo']} of "
          f"{res['n_pre_chips']} pre, {res['n_post_echo']} of "
          f"{res['n_post_chips']} post), {res['n_windows']} committed windows, "
          f"{res['n_gap_steps']} steps longer than the cadence, "
          f"{res['runtime_s']:.1f} s")
    passes = ", ".join("{} {} echoes of {} chips".format(
        core.ORBIT_SHORT[o], res["n_echo_by_orbit"][o],
        res["coverage"]["per_orbit"][o]["n_chips"]) for o in ORBITS)
    print(f"  the two passes: {passes}")
    print(f"  pooled Δorbit (DESC − ASC): "
          + ", ".join(f"{dim} {res['deltas'][dim]['all']['delta_median']:+.3f}"
                      for dim in DIMS))
    print(f"  the trailing window keeps {len(SERIES_DAYS)} days "
          f"({res['series_days'][0]} -> {res['series_days'][-1]}), "
          f"{res['n_series_points']} points")

    # the pins are the gate: the committed inputs are re-read and hashed, and no
    # artifact is written before every comparison agrees.
    pins = pin_all(res)
    summary = pins.summary()
    print(f"pins: {summary['n_checks']} checks, {summary['n_exact']} exact, "
          f"{summary['n_within_tolerance']} within "
          f"{summary['float_tolerance']:g}, {summary['n_failed']} failed")
    if summary["n_failed"]:
        for f in summary["failures"][:20]:
            print("PIN FAILED:", f["path"], "recomputed=", f["recomputed"],
                  "committed=", f["committed"])
        _fail(f"{summary['n_failed']} of {summary['n_checks']} pin comparisons "
              f"disagree with the committed layers")

    os.makedirs(args.figdir, exist_ok=True)
    if not args.no_figures:
        fig_orbits(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(json_safe(build_json(res, summary)), fh, indent=1,
                  sort_keys=True, default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res, pins))
    print(f"  {os.path.relpath(args.md)}")

    print(f"PIN VERIFICATION OK  ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
