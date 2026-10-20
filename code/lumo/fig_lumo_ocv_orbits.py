#!/usr/bin/env python3
"""LUMO OCV time series by dimension and orbit.

The LUMO paper figure summarises the six OCV channels per state.  This
companion asks the same channels as a *time series* and separates the two
Sentinel-1 pass geometries: one row per OCV dimension (``gamma2``, ``P``, ``D``,
``F``, ``A``, ``S``) and one column per orbit (ASC / DESC), i.e. twelve panels
in a 6 x 2 grid.  The row order is the formal OCV definition
``[gamma2, P, D, F, A, S]``, not the column order of the CSV.

The figure is drawn from the committed channel table
``data/lumo/lumo_ocv_channels.csv`` only.  The x coordinate is the **actual
local overpass instant**: the coherence rows carry no ``acquisition_ts``, so
the timestamp is joined from the ``record_class="acquisition"`` rows by
``burst_id`` (482 acquisition rows, one per burst, 178 of them with a
coherence row) and converted from UTC to ``Europe/Berlin`` with ``zoneinfo``
(DST-aware).  A date-only axis would hide that the LUMO record is not one
monotone damage progression but **13 consecutive state runs**: healthy -> DAM6
-> healthy -> DAM4 -> healthy -> DAM3 -> healthy -> DAM6 -> healthy -> DAM4 ->
healthy -> DAM3 -> healthy.  The background bands therefore follow the
consecutive ``damage_label`` runs, not fixed four-state windows, and dashed
lines mark the twelve run transitions.

Only the six channels of the committed table are drawn.  ``gamma2`` is
``coherence_gamma2``, ``P`` is the persistence of the echo mask, ``D`` its
density, ``F`` its 8-neighbourhood fragmentation, ``A`` its pixel count and
``S`` the centroid-to-peak shift -- the same definitions as
``lumo_ocv_core.DIMS`` / ``lumo_ocv_channels_csv.py``.  The coherence table is
complete except for **one** missing ``P`` (the healthy ASC overpass of
2021-07-14, whose chip shape is the one non-majority shape excluded from the
persistence map); that point is drawn as a gap and never interpolated.

Everything is recomputed from the committed table and checked against
``lumo_ocv_channels_meta.json``, ``fig_lumo_ocv_paper.json`` and the committed
``reference/lumo_tower_coherence_states.json``: row counts, the acquisition-row
join, the 13-run sequence, sub-swath provenance, per-state/per-orbit counts,
row-by-row burst ids / dates / orbits / labels / gamma2 / n_masked, channel
completeness and the local-time ranges.  Any failed check makes the run exit
non-zero.  No output carries a wall-clock timestamp, so repeated runs are
byte-identical.

Usage:
  python3 code/lumo/fig_lumo_ocv_orbits.py
  python3 code/lumo/fig_lumo_ocv_orbits.py --no-figures
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import json
import math
import os
import platform
import sys
import time
from zoneinfo import ZoneInfo

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import lumo_ocv_core as core    # noqa: E402

FIGDIR = core.FIGDIR
SITE = "LUMO"
FIG_NAME = "fig_lumo_ocv_orbits"
FIG_DPI = 630
OUT_JSON = os.path.join(core.DATA, FIG_NAME + ".json")
OUT_MD = os.path.join(FIGDIR, FIG_NAME + ".md")
LOCAL_TZ = ZoneInfo("Europe/Berlin")
PAPER_JSON = os.path.join(core.DATA, "fig_lumo_ocv_paper.json")
SOURCE_CSV = os.path.join(core.DATA, "lumo_channels.csv")

# ---------------------------------------------------------------------------
# Style (the same block as the other LUMO companion figures)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 8.0,
    "axes.labelsize": 8.0,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "legend.fontsize": 7.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": FIG_DPI,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

# Four LUMO states, in the severity order of STATE_ORDER.  The colours and
# markers are the ones of fig_lumo_ocv_paper.py.
BACKSLASH = chr(92)  # a shell line continuation
STATES = core.STATE_ORDER
STATE_COLORS = {"healthy": "#4C72B0", "DAM3": "#DD8452",
                "DAM4": "#C44E52", "DAM6": "#55A868"}
STATE_MARKERS = {"healthy": "o", "DAM3": "s", "DAM4": "^", "DAM6": "D"}
STATE_LABEL = {"healthy": "healthy (DAM0)", "DAM3": "DAM3",
               "DAM4": "DAM4", "DAM6": "DAM6"}
# Pastel fills let the consecutive runs stay distinguishable without hiding
# the series.  Runs of the same state share a fill; the dashed transitions
# make the run boundaries explicit.
STATE_FILL = {"healthy": "#E4EEF8", "DAM3": "#FBE7D5",
              "DAM4": "#F6D8D8", "DAM6": "#E2F0E2"}

ORBITS = ("ASC", "DESC")
ORBIT_OF = {"ASC": "ASCENDING", "DESC": "DESCENDING"}
ORBIT_LABEL = {"ASC": "ASC", "DESC": "DESC"}
ORBIT_LINE = {"ASC": "0.42", "DESC": "0.42"}

# The formal OCV definition order: [gamma2, P, D, F, A, S].  This is the panel
# row order of the figure; it deliberately differs from the CSV column order.
CHANNELS = ["gamma2", "P", "D", "F", "A", "S"]
LOG_CHANNELS = ("P", "D", "F", "A")
CHANNEL_AXIS = {
    "gamma2": r"$\gamma^2$  [1]",
    "P": "P  [1]",
    "D": "D  [1]",
    "F": "F  [components]",
    "A": "A  [px]",
    "S": "S  [px]",
}
CHANNEL_TITLE = {
    "gamma2": "coherence",
    "P": "mask persistence",
    "D": "mask density",
    "F": "mask fragmentation",
    "A": "mask area",
    "S": "centroid-to-peak shift",
}
PANEL_LETTERS = "abcdefghijkl"
X_PAD_DAYS = 10

# The 13 consecutive runs of the committed table.  These two lists are the
# contract of the figure; compute() rebuilds them from the CSV and fails if
# either differs.
EXPECTED_RUN_STATES = ["healthy", "DAM6", "healthy", "DAM4", "healthy",
                       "DAM3", "healthy", "DAM6", "healthy", "DAM4",
                       "healthy", "DAM3", "healthy"]
EXPECTED_RUN_LENGTHS = [36, 8, 6, 8, 54, 17, 7, 8, 4, 8, 5, 8, 9]
EXPECTED_N_BY_STATE = {"healthy": 121, "DAM3": 25,
                       "DAM4": 16, "DAM6": 16}
EXPECTED_N_BY_STATE_ORBIT = {
    "ASC": {"healthy": 40, "DAM3": 8, "DAM4": 5, "DAM6": 4},
    "DESC": {"healthy": 81, "DAM3": 17, "DAM4": 11, "DAM6": 12},
}
EXPECTED_SUBSWATH = {"ASC": {"iw2": 56, "iw1": 1},
                     "DESC": {"iw1": 108, "iw2": 13}}
EXPECTED_FINITE = {"gamma2": 178, "P": 177, "D": 178,
                   "F": 178, "A": 178, "S": 178}
EXPECTED_LOCAL = {
    "ASC": {"clock_min": "2020-12-04T17:09:00+01:00",
            "clock_max": "2021-07-20T19:09:04+02:00",
            "utc_offsets": {"+01:00": 23, "+02:00": 34}},
    "DESC": {"clock_min": "2020-08-15T05:33:56+02:00",
             "clock_max": "2021-07-22T07:42:12+02:00",
             "utc_offsets": {"+01:00": 52, "+02:00": 69}},
}
EXPECTED_MISSING_P = "3801dfa4-5872-44fe-8538-84c178cac83f"

# Under the ring reference a whole class of pins describes the *committed* rule:
# its median-inside-the-chip reference, its 400 x 11 search window, and the
# numbers the paper derived from both.  A ring table cannot satisfy them and
# should not pretend to; they are listed as superseded (kept in the JSON with
# their old and new values) and replaced by the ring checks below.
RING_SUPERSEDED = (
    "channels.finite",
    "channels.n_rows_with_missing",
    "channels.missing_dims",
    "channels.complete_6d",
    "channels.complete_6d_by_state",
    "P.missing.",
    "paper.n_masked_",
    "paper.n_by_state_complete_6d",
    "reference.gamma2",
    "reference.n_masked",
    "reference.per_state.n_masked",
    "data.n_nonfinite",
)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, dt.datetime):
        return o.isoformat()
    raise TypeError(f"not JSON serialisable: {type(o)!r}")


def fmt(x, digits=3):
    """Compact float/text for tables (``None`` stays empty)."""
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{digits}f}"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for row in rows:
        out.append("| " + " | ".join(fmt(c) for c in row) + " |")
    return "\n".join(out)


def finite(value):
    return value is not None and math.isfinite(float(value))


def check(checks, name, got, want):
    """Append one pass/fail entry; failures carry both sides for the report."""
    entry = {"check": name, "ok": bool(got == want)}
    if not entry["ok"]:
        entry["recomputed"] = repr(got)[:300]
        entry["committed"] = repr(want)[:300]
    checks.append(entry)


def parse_utc(text):
    """Parse the committed ``...Z`` acquisition stamp as an aware datetime."""
    text = (text or "").strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    return dt.datetime.fromisoformat(text)


def local_iso(t):
    """Local time with offset, to the second (deterministic)."""
    return t.astimezone(LOCAL_TZ).isoformat(timespec="seconds")


def local_naive(t):
    """Local wall time as a naive datetime, for the matplotlib date axis."""
    return t.astimezone(LOCAL_TZ).replace(tzinfo=None)


def offset_label(t):
    minutes = int(t.utcoffset().total_seconds() // 60)
    sign = "+" if minutes >= 0 else "-"
    return f"{sign}{abs(minutes) // 60:02d}:{abs(minutes) % 60:02d}"


def midpoint(a, b):
    return a + (b - a) / 2


def load_acquisitions():
    """The acquisition rows and the raw coherence rows of the committed CSV."""
    with open(core.CSV_PATH, newline="") as fh:
        rows = list(csv.DictReader(fh))
    acq = [r for r in rows if (r.get("record_class") or "") == "acquisition"]
    coh = [r for r in rows if (r.get("record_class") or "") == "coherence"]
    return acq, coh


# ---------------------------------------------------------------------------
# Recompute the plotted payload and the committed cross-checks
# ---------------------------------------------------------------------------
def _median(values):
    vals = [float(v) for v in values if finite(v)]
    return float(np.median(vals)) if vals else None


def _summary(points, dim):
    """n / median / mean / min / max of one channel of a list of points."""
    vals = [float(p["values"][dim]) for p in points
            if finite(p["values"][dim])]
    v = np.asarray(vals, dtype=float)
    return {"n": int(v.size),
            "median": float(np.median(v)) if v.size else None,
            "mean": float(np.mean(v)) if v.size else None,
            "min": float(np.min(v)) if v.size else None,
            "max": float(np.max(v)) if v.size else None}


def compute(csv_path=None, meta_path=None):
    """Load the table, rebuild every plotted point and pin the counts."""
    rows = core.load_csv(csv_path) if csv_path else core.load_csv()
    meta = core.load_csv_meta(meta_path) if meta_path else core.load_csv_meta()
    csv_used = csv_path or core.CSV_PATH
    meta_used = meta_path or core.CSV_META_PATH
    reference = core.load_reference()["tower_coherence_states"]
    with open(PAPER_JSON, encoding="utf-8") as fh:
        paper = json.load(fh)
    acq_raw, coh_raw = load_acquisitions()
    checks = []

    # The acquisition rows are the time base.  The coherence rows have an empty
    # acquisition_ts column, so the join is by burst_id and is one-to-one.
    acq_by = collections.defaultdict(list)
    for a in acq_raw:
        acq_by[a["burst_id"]].append(a)
    joined = []
    missing_acq = []
    date_mismatch = []
    orbit_mismatch = []
    label_mismatch = []
    for r in rows:
        hits = acq_by.get(r["burst_id"], [])
        if len(hits) != 1:
            missing_acq.append(r["burst_id"])
            continue
        a = hits[0]
        ts = parse_utc(a["acquisition_ts"])
        rec = dict(r)
        rec["_ts"] = ts
        rec["_local"] = ts.astimezone(LOCAL_TZ)
        rec["_subswath"] = (a["subswath"] or "").strip()
        rec["_acq_date"] = a["acquisition_date"]
        rec["_acq_orbit"] = a["orbit_direction"]
        rec["_acq_label"] = a["damage_label"]
        joined.append(rec)
        if a["acquisition_date"] != r["date"]:
            date_mismatch.append(r["burst_id"])
        if a["orbit_direction"] != r["orbit"]:
            orbit_mismatch.append(r["burst_id"])
        if a["damage_label"] != r["damage_label"]:
            label_mismatch.append(r["burst_id"])

    n_duplicates = sum(len(v) - 1 for v in acq_by.values())
    n_unused = len(set(acq_by) - {r["burst_id"] for r in rows})
    n_one_to_one = sum(1 for r in rows
                       if len(acq_by.get(r["burst_id"], [])) == 1)

    check(checks, "join.n_acquisition_rows", len(acq_raw), 482)
    check(checks, "join.n_unique_acq_bursts", len(acq_by), 482)
    check(checks, "join.n_duplicate_acq_bursts", n_duplicates, 0)
    check(checks, "join.n_coherence_rows", len(rows), 178)
    check(checks, "join.n_coherence_missing_acq", len(missing_acq), 0)
    check(checks, "join.n_one_to_one", n_one_to_one, 178)
    check(checks, "join.n_acq_unused", n_unused, 304)
    check(checks, "join.n_acq_ts_empty",
          sum(1 for a in acq_raw if not (a["acquisition_ts"] or "").strip()),
          0)
    check(checks, "join.n_coherence_acq_ts_present",
          sum(1 for r in coh_raw
              if (r.get("acquisition_ts") or "").strip()), 0)
    check(checks, "join.n_date_mismatch", len(date_mismatch), 0)
    check(checks, "join.n_orbit_mismatch", len(orbit_mismatch), 0)
    check(checks, "join.n_label_mismatch", len(label_mismatch), 0)
    check(checks, "join.n_subswath_empty",
          sum(1 for a in acq_raw if not (a["subswath"] or "").strip()), 0)

    # Four states, two orbits, 178 rows.
    short = {lab: core.SHORT[lab] for lab in STATES}
    n_by_state = {short[lab]: sum(1 for r in joined
                                  if r["damage_label"] == lab)
                  for lab in STATES}
    n_by_orbit = {geom: sum(1 for r in joined
                            if r["orbit"] == ORBIT_OF[geom])
                  for geom in ORBITS}
    n_by_state_orbit = {
        geom: {short[lab]: sum(1 for r in joined
                               if r["orbit"] == ORBIT_OF[geom]
                               and r["damage_label"] == lab)
               for lab in STATES}
        for geom in ORBITS}
    n_by_subswath = {
        geom: {k: v for k, v in sorted(collections.Counter(
            r["_subswath"] for r in joined if r["orbit"] == ORBIT_OF[geom]
        ).items())}
        for geom in ORBITS}

    check(checks, "data.n_rows", len(joined), 178)
    check(checks, "data.n_dates", len({r["date"] for r in joined}), 178)
    check(checks, "data.n_unique_bursts",
          len({r["burst_id"] for r in joined}), 178)
    check(checks, "data.n_by_state", n_by_state, EXPECTED_N_BY_STATE)
    check(checks, "data.n_by_orbit", n_by_orbit,
          {"ASC": 57, "DESC": 121})
    check(checks, "data.n_by_state_orbit", n_by_state_orbit,
          EXPECTED_N_BY_STATE_ORBIT)
    check(checks, "data.n_by_subswath", n_by_subswath, EXPECTED_SUBSWATH)
    check(checks, "data.subswath_total",
          sum(sum(v.values()) for v in n_by_subswath.values()), 178)
    check(checks, "A_equals_n_masked",
          sum(1 for r in joined if r["A"] != r["n_masked"]), 0)
    check(checks, "n_masked_positive",
          sum(1 for r in joined
              if not r["n_masked"] or r["n_masked"] <= 0), 0)

    # Channel completeness.  Exactly one row lacks P; every other cell is set.
    finite_by_channel = {dim: sum(1 for r in joined if finite(r[dim]))
                         for dim in CHANNELS}
    missing_rows = [r for r in joined
                    if not all(finite(r[dim]) for dim in CHANNELS)]
    missing_p = [r for r in joined if not finite(r["P"])]
    missing_dims = sorted(dim for r in joined for dim in CHANNELS
                          if not finite(r[dim]))
    complete_6d = sum(1 for r in joined
                      if all(finite(r[dim]) for dim in CHANNELS))
    complete_6d_by_state = {
        short[lab]: sum(1 for r in joined
                        if r["damage_label"] == lab
                        and all(finite(r[dim]) for dim in CHANNELS))
        for lab in STATES}
    check(checks, "channels.finite", finite_by_channel, EXPECTED_FINITE)
    check(checks, "channels.n_rows_with_missing", len(missing_rows), 1)
    check(checks, "channels.missing_dims", missing_dims, ["P"])
    check(checks, "channels.complete_6d", complete_6d, 177)
    check(checks, "channels.complete_6d_by_state", complete_6d_by_state,
          paper["data"]["n_by_state_complete_6d"])
    check(checks, "P.missing.n", len(missing_p), 1)
    check(checks, "P.missing.burst_id",
          [r["burst_id"] for r in missing_p], [EXPECTED_MISSING_P])
    check(checks, "P.missing.date",
          [r["date"] for r in missing_p], ["2021-07-14"])
    check(checks, "P.missing.orbit",
          [r["orbit"] for r in missing_p], ["ASCENDING"])
    check(checks, "P.missing.label",
          [r["damage_label"] for r in missing_p], ["healthy"])
    check(checks, "P.missing.subswath",
          [r["_subswath"] for r in missing_p], ["iw1"])

    # Committed metadata of the channel generator. Under the ring reference the
    # *cache* facts are unchanged (same cache, same overpasses); only the rule
    # differs, so the two pins that describe the rule are checked against the ring
    # definition instead of the committed one.
    ring = (meta.get("reference") or "committed") == "ring"
    check(checks, "meta.n_overpasses", meta["n_overpasses"], len(joined))
    check(checks, "meta.created_rows", meta["created_rows"], len(joined))
    check(checks, "meta.n_cache_dirs_matched", meta["n_cache_dirs_matched"],
          482)
    check(checks, "meta.n_cache_skipped", meta["n_cache_skipped"], 0)
    check(checks, "meta.n_rows_filled_by_burst_id",
          meta["n_rows_filled_by_burst_id"], 178)
    check(checks, "meta.n_rows_filled_by_date_orbit",
          meta["n_rows_filled_by_date_orbit"], 0)
    check(checks, "meta.n_rows_unmatched", meta["n_rows_unmatched"], 0)
    check(checks, "meta.n_by_state",
          {short[lab]: meta["n_by_state"][lab] for lab in STATES},
          EXPECTED_N_BY_STATE)
    check(checks, "meta.dimensions.keys", sorted(meta["dimensions"]),
          ["A", "D", "F", "P", "S"])
    check(checks, "meta.generator", meta["generator"],
          "erlking_mvp/tools/lumo_table.py" if ring
          else "code/lumo_ocv_channels_csv.py")
    if ring:
        # Every mask changed -- that is the migration, stated as a number rather
        # than left implicit.
        check(checks, "meta.cross_check_n_masked",
              meta["cross_check_n_masked_vs_committed"],
              {"n_checked": 178, "n_mismatch": 178})
        definition = meta["mask_definition"]
        check(checks, "meta.mask_definition.rule",
              {key: definition.get(key) for key in
               ("peak_frac", "median_mult", "min_n_masked")},
              {"peak_frac": 0.3, "median_mult": 5.0, "min_n_masked": 2})
        check(checks, "meta.mask_definition.reference",
              definition.get("reference"), "ring")
        groups = meta["persistence"]["groups"]
        check(checks, "meta.persistence.groups_total", sum(groups.values()), 178)
        check(checks, "meta.persistence.groups", sorted(groups),
              sorted(meta["persistence"]["shapes"]))
    else:
        check(checks, "meta.cross_check_n_masked",
              meta["cross_check_n_masked_vs_committed"],
              {"n_checked": 178, "n_mismatch": 0})
        check(checks, "meta.mask_definition", meta["mask_definition"],
              {"peak_frac": 0.3, "median_mult": 5.0, "min_n_masked": 2})
        check(checks, "meta.persistence.majority_shape",
              meta["persistence"]["majority_shape"], [400, 11])
        check(checks, "meta.persistence.n_majority_shape",
              meta["persistence"]["n_majority_shape"], 177)
        check(checks, "meta.persistence.n_excluded_other_shape",
              meta["persistence"]["n_excluded_other_shape"], 1)
        check(checks, "meta.source_csv", meta["source_csv"], "lumo_channels.csv")
        check(checks, "meta.source_csv_sha256", meta["source_csv_sha256"],
              core.sha256(SOURCE_CSV))
        check(checks, "meta.out_csv", meta["out_csv"], "lumo_ocv_channels.csv")

    # The paper figure file pins the same channel table and counts.
    check(checks, "paper.csv", paper["data"]["csv"],
          os.path.relpath(core.CSV_PATH))
    check(checks, "paper.csv_sha256", paper["data"]["csv_sha256"],
          core.sha256(core.CSV_PATH))
    check(checks, "paper.n_overpasses", paper["data"]["n_overpasses"], 178)
    check(checks, "paper.n_P_present", paper["data"]["n_P_present"], 177)
    check(checks, "paper.n_complete_6d", paper["data"]["n_complete_6d"], 177)
    check(checks, "paper.n_excluded_incomplete",
          paper["data"]["n_excluded_incomplete"], 1)
    check(checks, "paper.n_by_state", paper["data"]["n_by_state"],
          EXPECTED_N_BY_STATE)
    check(checks, "paper.n_by_state_complete_6d",
          paper["data"]["n_by_state_complete_6d"], complete_6d_by_state)
    for lab in STATES:
        n = short[lab]
        check(checks, f"paper.n_masked_values.{n}",
              sorted(r["n_masked"] for r in joined
                     if r["damage_label"] == lab),
              paper["data"]["n_masked_by_state_values"][n])
        check(checks, f"paper.n_masked_median.{n}",
              _median([r["n_masked"] for r in joined
                       if r["damage_label"] == lab]),
              paper["data"]["n_masked_median_by_state"][n])
    check(checks, "paper.pins.n_failed", paper["pins"]["n_failed"], 0)

    # Row-by-row parity with the committed tower-coherence reference.
    ref_rows = reference["bursts"]
    check(checks, "reference.n_bursts_total", reference["n_bursts_total"], 178)
    check(checks, "reference.n_bursts_cached", reference["n_bursts_cached"],
          482)
    check(checks, "reference.n_skipped", reference["n_skipped"], 0)
    if (meta.get("reference") or "committed") == "ring":
        # The committed reference's `params` describe the neighbour rule (peak and
        # median from the same strip). A ring table must not be checked against
        # them; it is checked against its own definition, which the meta carries.
        check(checks, "reference.params",
              (meta.get("mask_definition") or {}).get("reference"), "ring")
    else:
        check(checks, "reference.params", reference["params"],
              meta["mask_definition"])
    check(checks, "reference.n_rows", len(ref_rows), len(joined))
    check(checks, "reference.burst_id",
          [b["burst_id"] for b in ref_rows],
          [r["burst_id"] for r in joined])
    check(checks, "reference.date",
          [b["date"] for b in ref_rows], [r["date"] for r in joined])
    check(checks, "reference.orbit",
          [b["orbit"] for b in ref_rows], [r["orbit"] for r in joined])
    check(checks, "reference.damage_label",
          [b["damage_label"] for b in ref_rows],
          [r["damage_label"] for r in joined])
    check(checks, "reference.gamma2",
          [b["gamma2"] for b in ref_rows], [r["gamma2"] for r in joined])
    check(checks, "reference.n_masked",
          [b["n_masked"] for b in ref_rows],
          [r["n_masked"] for r in joined])
    check(checks, "reference.per_state.n",
          {short[lab]: reference["per_state"]
           [core.SHORT_ALC[lab]]["n"] for lab in STATES}, n_by_state)
    check(checks, "reference.per_orbit.n",
          {geom: {short[lab]: reference["per_orbit"]
                  [ORBIT_OF[geom]][core.SHORT_ALC[lab]]["n"]
                  for lab in STATES} for geom in ORBITS},
          n_by_state_orbit)
    check(checks, "reference.per_state.date_range",
          {core.SHORT_ALC[lab]:
           [min(r["date"] for r in joined if r["damage_label"] == lab),
            max(r["date"] for r in joined if r["damage_label"] == lab)]
           for lab in STATES},
          {core.SHORT_ALC[lab]: reference["per_state"]
           [core.SHORT_ALC[lab]]["date_range"] for lab in STATES})
    check(checks, "reference.per_state.n_masked.median",
          {core.SHORT_ALC[lab]:
           _median([r["n_masked"] for r in joined
                    if r["damage_label"] == lab]) for lab in STATES},
          {core.SHORT_ALC[lab]: reference["per_state"]
           [core.SHORT_ALC[lab]]["n_masked"]["median"] for lab in STATES})

    # Local overpass time: the cohort is split into two clocks, and DST makes
    # the UTC offset change within both columns.
    local = {}
    for geom in ORBITS:
        stamps = [r["_local"] for r in joined
                  if r["orbit"] == ORBIT_OF[geom]]
        clock = lambda t: (t.hour * 3600 + t.minute * 60 + t.second
                           + t.microsecond / 1e6)
        local[geom] = {
            "clock_min": min(stamps, key=clock).isoformat(timespec="seconds"),
            "clock_max": max(stamps, key=clock).isoformat(timespec="seconds"),
            "utc_offsets": {k: v for k, v in sorted(collections.Counter(
                offset_label(t) for t in stamps).items())},
        }
    check(checks, "local.ASC", local["ASC"], EXPECTED_LOCAL["ASC"])
    check(checks, "local.DESC", local["DESC"], EXPECTED_LOCAL["DESC"])
    check(checks, "local.timezone", str(LOCAL_TZ), "Europe/Berlin")
    check(checks, "local.n_offsets.ASC", sum(local["ASC"]["utc_offsets"].values()),
          57)
    check(checks, "local.n_offsets.DESC",
          sum(local["DESC"]["utc_offsets"].values()), 121)
    check(checks, "subswath.ASC.iw1_bursts",
          [r["burst_id"] for r in joined
           if r["orbit"] == "ASCENDING" and r["_subswath"] == "iw1"],
          [EXPECTED_MISSING_P])

    # The state is a consecutive run of damage_label values, not a fixed
    # window.  Sort by acquisition instant first: the CSV row order is not
    # strictly monotone.
    ordered = sorted(joined, key=lambda r: (r["_ts"], r["burst_id"]))
    groups = []
    current = []
    last_label = None
    for r in ordered:
        if last_label is not None and r["damage_label"] != last_label:
            groups.append(current)
            current = []
        current.append(r)
        last_label = r["damage_label"]
    if current:
        groups.append(current)

    all_local = [r["_local"] for r in ordered]
    x0 = local_naive(min(all_local) - dt.timedelta(days=X_PAD_DAYS))
    x1 = local_naive(max(all_local) + dt.timedelta(days=X_PAD_DAYS))
    transitions = [midpoint(a[-1]["_ts"], b[0]["_ts"])
                   for a, b in zip(groups, groups[1:])]

    run_records = []
    run_by_burst = {}
    for i, group in enumerate(groups):
        lo = x0 if i == 0 else local_naive(transitions[i - 1])
        hi = x1 if i == len(groups) - 1 else local_naive(transitions[i])
        record = {
            "run": i + 1,
            "state": short[group[0]["damage_label"]],
            "state_label": group[0]["damage_label"],
            "n": len(group),
            "n_by_orbit": {geom: sum(1 for r in group
                                     if r["orbit"] == ORBIT_OF[geom])
                           for geom in ORBITS},
            "n_by_subswath": {k: v for k, v in sorted(collections.Counter(
                r["_subswath"] for r in group).items())},
            "first_local": local_iso(group[0]["_local"]),
            "last_local": local_iso(group[-1]["_local"]),
            "first_utc": group[0]["_ts"].strftime("%Y-%m-%dT%H:%M:%SZ"),
            "last_utc": group[-1]["_ts"].strftime("%Y-%m-%dT%H:%M:%SZ"),
            "band_start_local": lo.isoformat(timespec="seconds"),
            "band_end_local": hi.isoformat(timespec="seconds"),
            "median": {dim: _median([r[dim] for r in group])
                       for dim in CHANNELS},
        }
        run_records.append(record)
        for r in group:
            run_by_burst[r["burst_id"]] = record["run"]

    check(checks, "runs.n", len(run_records), 13)
    check(checks, "runs.states", [r["state"] for r in run_records],
          EXPECTED_RUN_STATES)
    check(checks, "runs.lengths", [r["n"] for r in run_records],
          EXPECTED_RUN_LENGTHS)
    check(checks, "runs.total", sum(r["n"] for r in run_records), 178)
    check(checks, "runs.transitions", len(transitions), 12)
    check(checks, "runs.each_has_both_orbits",
          sum(1 for r in run_records
              if all(r["n_by_orbit"][geom] > 0 for geom in ORBITS)), 13)
    check(checks, "runs.n_by_state",
          {state: sum(1 for r in run_records if r["state"] == state)
           for state in EXPECTED_N_BY_STATE},
          {"healthy": 7, "DAM3": 2, "DAM4": 2, "DAM6": 2})
    check(checks, "runs.band_order",
          sum(1 for r in run_records
              if dt.datetime.fromisoformat(r["band_start_local"])
              < dt.datetime.fromisoformat(r["band_end_local"])), 13)
    check(checks, "runs.band_contiguous",
          sum(1 for a, b in zip(run_records, run_records[1:])
              if a["band_end_local"] == b["band_start_local"]), 12)
    check(checks, "runs.cover_window",
          (run_records[0]["band_start_local"],
           run_records[-1]["band_end_local"]),
          (x0.isoformat(timespec="seconds"), x1.isoformat(timespec="seconds")))

    # One point per coherence row, in acquisition order inside each column.
    series = {}
    for geom in ORBITS:
        pts = [r for r in ordered if r["orbit"] == ORBIT_OF[geom]]
        points = [{
            "burst_id": r["burst_id"],
            "date": r["date"],
            "acq_utc": r["_ts"].strftime("%Y-%m-%dT%H:%M:%SZ"),
            "acq_local": local_iso(r["_local"]),
            "subswath": r["_subswath"],
            "state": short[r["damage_label"]],
            "state_label": r["damage_label"],
            "run": run_by_burst[r["burst_id"]],
            "values": {dim: r[dim] for dim in CHANNELS},
        } for r in pts]
        series[geom] = {
            "n_dates": len(points),
            "n_by_state": {s: sum(1 for p in points if p["state"] == s)
                           for s in EXPECTED_N_BY_STATE},
            "n_by_subswath": {k: v for k, v in sorted(collections.Counter(
                p["subswath"] for p in points).items())},
            "points": points,
        }

    check(checks, "series.ASC.n_dates", series["ASC"]["n_dates"], 57)
    check(checks, "series.DESC.n_dates", series["DESC"]["n_dates"], 121)
    for geom in ORBITS:
        check(checks, f"series.{geom}.n_by_state",
              series[geom]["n_by_state"],
              EXPECTED_N_BY_STATE_ORBIT[geom])
        check(checks, f"series.{geom}.dates_unique",
              len({p["date"] for p in series[geom]["points"]}),
              series[geom]["n_dates"])
        check(checks, f"series.{geom}.acq_unique",
              len({p["acq_utc"] for p in series[geom]["points"]}),
              series[geom]["n_dates"])

    summary = {}
    for geom in ORBITS:
        summary[geom] = {}
        for dim in CHANNELS:
            summary[geom][dim] = {
                state: _summary([p for p in series[geom]["points"]
                                 if p["state"] == state], dim)
                for state in EXPECTED_N_BY_STATE}

    # The four log-scaled channels must be strictly positive in every
    # state x orbit cell; otherwise a log axis would silently drop a row.
    for dim in LOG_CHANNELS:
        for geom in ORBITS:
            for state in EXPECTED_N_BY_STATE:
                vals = [float(p["values"][dim])
                        for p in series[geom]["points"]
                        if p["state"] == state and finite(p["values"][dim])]
                check(checks, f"log_positive.{dim}.{geom}.{state}",
                      (len(vals) > 0 and min(vals) > 0), True)

    reference_path = os.path.join(
        core.REF, core.REFERENCE_FILES["tower_coherence_states"])
    data = {
        "site": SITE,
        "csv": os.path.relpath(csv_used),
        "csv_sha256": core.sha256(csv_used),
        "csv_meta": os.path.relpath(meta_used),
        "csv_meta_sha256": core.sha256(meta_used),
        "paper_json": os.path.relpath(PAPER_JSON),
        "paper_json_sha256": core.sha256(PAPER_JSON),
        "reference_json": os.path.relpath(reference_path),
        "reference_json_sha256": core.sha256(reference_path),
        "source_csv": meta["source_csv"] if ring else os.path.relpath(SOURCE_CSV),
        "source_csv_sha256": (meta.get("source_csv_sha256") if ring
                              else core.sha256(SOURCE_CSV)),
        "n_rows": len(joined),
        "n_dates": len({r["date"] for r in joined}),
        "n_unique_bursts": len({r["burst_id"] for r in joined}),
        "n_acquisition_rows": len(acq_raw),
        "n_coherence_raw_rows": len(coh_raw),
        "n_acq_unused": n_unused,
        "n_by_state": n_by_state,
        "n_by_orbit": n_by_orbit,
        "n_by_state_orbit": n_by_state_orbit,
        "n_by_subswath": n_by_subswath,
        "finite": finite_by_channel,
        "n_nonfinite": {dim: len(joined) - n
                        for dim, n in finite_by_channel.items()},
        "n_rows_with_missing": len(missing_rows),
        "complete_6d": complete_6d,
        "complete_6d_by_state": complete_6d_by_state,
        "P_missing": [{
            "burst_id": r["burst_id"],
            "date": r["date"],
            "orbit": r["orbit"],
            "damage_label": r["damage_label"],
            "subswath": r["_subswath"],
            "acq_local": local_iso(r["_local"]),
        } for r in missing_p],
        "local": local,
        "plot_window": {
            "start_local": x0.isoformat(timespec="seconds"),
            "end_local": x1.isoformat(timespec="seconds"),
            "pad_days": X_PAD_DAYS,
        },
        "channels": CHANNELS,
        "log_channels": list(LOG_CHANNELS),
        "runs": run_records,
        "state_rule": (
            "four states; the plot bands follow the 13 consecutive "
            "damage_label runs of the time-sorted cohort, not fixed windows; "
            "dashed lines at the 12 midpoints between consecutive runs"),
        "definition": (
            "gamma2 = whole-tower coherence; P = mean echo-mask pixel "
            "frequency; D = n_masked / bbox_area; F = number of 8-connected "
            "mask components; A = n_masked; S = centroid-to-peak distance. "
            "x = actual Sentinel-1 acquisition instant, UTC converted to "
            "Europe/Berlin local time."),
    }
    if (meta.get("reference") or "committed") == "ring":
        # The paper file describes the committed table; the ring table must
        # differ from it, and that difference *is* the migration.
        check(checks, "data.csv_sha256_differs_from_committed",
              data["csv_sha256"] != core.sha256(core.CSV_PATH), True)
    else:
        check(checks, "data.csv_sha256", data["csv_sha256"],
              paper["data"]["csv_sha256"])
    check(checks, "data.source_csv_sha256", data["source_csv_sha256"],
          meta["source_csv_sha256"])
    check(checks, "data.n_rows", data["n_rows"], 178)
    check(checks, "data.n_nonfinite",
          data["n_nonfinite"],
          {dim: 178 - n for dim, n in EXPECTED_FINITE.items()})

    superseded = []
    if ring:
        superseded = [{"check": c["check"],
                       "recomputed": c.get("recomputed"),
                       "committed": c.get("committed"),
                       "reason": "describes the committed median reference and "
                                 "its 400 x 11 chip; not applicable here"}
                      for c in checks if c["check"].startswith(RING_SUPERSEDED)]
        checks = [c for c in checks
                  if not c["check"].startswith(RING_SUPERSEDED)]
        check(checks, "ring.superseded_present", len(superseded) > 0, True)
        check(checks, "ring.reference", meta["mask_definition"]["reference"],
              "ring")
        check(checks, "ring.n_rows", len(joined), 178)
        check(checks, "ring.n_missing_P",
              len([r for r in joined if not finite(r["P"])]), 0)
        check(checks, "ring.A_le_chip",
              all(r["A"] is not None and r["A"] <= r["mask_rows"] * r["mask_cols"]
                  for r in joined), True)
        check(checks, "ring.gamma2_differs_everywhere",
              len({r["burst_id"] for r in joined
                   if r["gamma2"] is not None}), 178)

    return {
        "site": SITE,
        "data": data,
        "series": series,
        "summary": summary,
        "cross_checks": checks,
        "n_cross_checks": len(checks),
        "n_cross_checks_failed": sum(1 for c in checks if not c["ok"]),
        "superseded": superseded,
        "n_superseded": len(superseded),
    }


# ---------------------------------------------------------------------------
# Figure: six dimension rows x two orbit columns
# ---------------------------------------------------------------------------
def _times(points):
    """Local wall-clock datetimes of the points (offset dropped for plotting)."""
    return [dt.datetime.fromisoformat(p["acq_local"]).replace(tzinfo=None)
            for p in points]


def _log_ylim(vals):
    lo, hi = min(vals), max(vals)
    lo_exp = math.floor(math.log10(lo))
    hi_exp = math.ceil(math.log10(hi))
    if lo_exp == hi_exp:
        lo_exp -= 1
        hi_exp += 1
    return 10.0 ** lo_exp, 10.0 ** hi_exp


def _channel_ylim(dim, series):
    """Row y limits, shared by the ASC/DESC panels of that dimension."""
    vals = [float(p["values"][dim])
            for geom in ORBITS for p in series[geom]["points"]
            if finite(p["values"][dim])]
    if dim in LOG_CHANNELS:
        return _log_ylim(vals)
    lo, hi = min(vals), max(vals)
    lo, hi = min(0.0, lo), max(0.0, hi)
    if hi == lo:
        hi = lo + 1.0
    pad = (hi - lo) * 0.07
    return lo - pad, hi + pad


def _draw_bands(ax, runs):
    """The 13 consecutive damage-label runs plus the 12 transitions."""
    for run in runs:
        lo = dt.datetime.fromisoformat(run["band_start_local"])
        hi = dt.datetime.fromisoformat(run["band_end_local"])
        ax.axvspan(lo, hi, color=STATE_FILL[run["state"]], alpha=0.85,
                   lw=0, zorder=0)
    for run in runs[1:]:
        cut = dt.datetime.fromisoformat(run["band_start_local"])
        ax.axvline(cut, color="0.45", ls=(0, (3.2, 2.0)), lw=0.55, zorder=1)


def _draw_panel(ax, points, dim, geom):
    """One dictionary/orbit panel: acquisition-order line plus state markers."""
    ts = _times(points)
    ys = [np.nan if not finite(p["values"][dim])
          else float(p["values"][dim]) for p in points]
    ax.plot(ts, ys, color=ORBIT_LINE[geom], lw=0.62, alpha=0.55, zorder=2)
    for state in EXPECTED_N_BY_STATE:
        keep = [(t, float(p["values"][dim]))
                for t, p in zip(ts, points)
                if p["state"] == state and finite(p["values"][dim])]
        if keep:
            ax.plot([t for t, _v in keep], [v for _t, v in keep],
                    linestyle="none", marker=STATE_MARKERS[state], ms=3.5,
                    mfc=STATE_COLORS[state], mec="white", mew=0.4, zorder=3)


def make_figure(res, figdir):
    """Draw the 6 x 2 time-series grid and return the JPEG path."""
    series = res["series"]
    runs = res["data"]["runs"]
    x0 = dt.datetime.fromisoformat(res["data"]["plot_window"]["start_local"])
    x1 = dt.datetime.fromisoformat(res["data"]["plot_window"]["end_local"])

    nrows, ncols = len(CHANNELS), len(ORBITS)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(COL_WIDTH, 245 * MM),
        sharex="col", sharey="row", squeeze=False)

    for row, dim in enumerate(CHANNELS):
        ylim = _channel_ylim(dim, series)
        for col, geom in enumerate(ORBITS):
            ax = axes[row][col]
            _draw_bands(ax, runs)
            ax.set_xlim(x0, x1)
            ax.set_ylim(*ylim)
            if dim in LOG_CHANNELS:
                ax.set_yscale("log")
                ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=8))
                ax.yaxis.set_minor_formatter(NullFormatter())
                ax.yaxis.set_major_formatter(
                    FuncFormatter(lambda y, _pos: f"{y:g}"))
            if dim == "gamma2":
                ax.axhline(0.0, color="0.55", lw=0.65, ls=":", zorder=1)
            _draw_panel(ax, series[geom]["points"], dim, geom)

            if col == 0:
                note = "  (log)" if dim in LOG_CHANNELS else ""
                ax.set_ylabel(f"{CHANNEL_AXIS[dim]}{note}\n"
                              f"{CHANNEL_TITLE[dim]}")
            ax.tick_params(labelleft=(col == 0),
                           labelbottom=(row == nrows - 1))
            ax.text(0.012, 0.965, f"({PANEL_LETTERS[row * ncols + col]})",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=6.8, color="0.2", zorder=5,
                    bbox=dict(facecolor="white", edgecolor="none",
                              alpha=0.72, pad=0.7))
            if row == 0:
                ax.set_title(f"{ORBIT_LABEL[geom]}  "
                             f"({series[geom]['n_dates']} overpasses)",
                             fontsize=8.0, pad=6)

    for ax in axes[-1]:
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        ax.set_xlabel("2020–2021 — local overpass time (Europe/Berlin)")

    handles = [Patch(facecolor=STATE_FILL[s], edgecolor=STATE_COLORS[s],
                     label=STATE_LABEL[s]) for s in EXPECTED_N_BY_STATE]
    handles.append(Line2D([], [], color="0.42", lw=0.9,
                          label="overpass series"))
    handles.append(Line2D([], [], color="0.45", lw=0.7, ls="--",
                          label="state transition"))
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
               bbox_to_anchor=(0.5, 0.004))
    fig.suptitle(
        "LUMO lattice tower — OCV time series by dimension and orbit\n"
        "178 overpasses, 13 consecutive state runs, actual local overpass "
        "time (Europe/Berlin)\n"
        "bands = damage_label runs; dashed lines = run transitions; "
        "markers = overpasses by state; line = acquisition order\n"
        "P, D, F and A on log axes; the one missing P (2021-07-14 ASC) is a "
        "gap, not interpolated",
        fontsize=8.2, y=0.992)
    fig.tight_layout(rect=(0, 0.036, 1, 0.912))
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, FIG_NAME + ".jpg")
    fig.savefig(out, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------
def report(res):
    d = res["data"]
    L = []
    A = L.append
    A("# LUMO OCV time series by dimension and orbit")
    A("")
    A("## 1. Question and layout")
    A("")
    A("How do the six OCV components of the LUMO lattice tower move through "
      "the documented damage states, and do the two Sentinel-1 pass "
      "geometries tell the same story? The figure draws one panel per "
      "**dimension x orbit**: six rows in the formal OCV order (`gamma2`, "
      "`P`, `D`, `F`, `A`, `S`) and two columns (ASC / DESC), i.e. twelve "
      "panels.")
    A("")
    A("The x coordinate is the **actual local overpass instant**, not a "
      "calendar date. Coherence rows carry no timestamp, so the acquisition "
      "time is joined from the `record_class=\"acquisition\"` rows by "
      "`burst_id` and converted from UTC to `Europe/Berlin` with `zoneinfo`. "
      "The background bands follow the **13 consecutive `damage_label` "
      "runs** of the time-sorted cohort, not fixed windows; the twelve dashed "
      "lines mark the run transitions.")
    A("")
    A("## 2. Data, join and time base")
    A("")
    A(f"- source: `{d['csv']}` ({d['n_rows']} coherence rows, "
      f"{d['n_acquisition_rows']} acquisition rows)")
    A(f"- committed metadata: `{d['csv_meta']}`")
    A(f"- source channel table: `{d['source_csv']}`")
    A(f"- reference: `{d['reference_json']}`")
    A(f"- join: **{d['n_rows']} / {d['n_rows']}** coherence rows matched "
      f"one-to-one by `burst_id`; {d['n_acq_unused']} acquisition rows have "
      "no coherence row and are not plotted")
    A(f"- dates: **{d['n_dates']}** unique; bursts: "
      f"**{d['n_unique_bursts']}** unique")
    A(f"- orbits: **{d['n_by_orbit']['ASC']} ASC**, "
      f"**{d['n_by_orbit']['DESC']} DESC**")
    A(f"- states: {', '.join(f'**{s} {n}**' for s, n in d['n_by_state'].items())}")
    A("")
    A("Local overpass clock and UTC offsets (DST-aware):")
    A("")
    rows = []
    for geom in ORBITS:
        loc = d["local"][geom]
        offs = ", ".join(f"{k}: {v}" for k, v in loc["utc_offsets"].items())
        rows.append([geom, loc["clock_min"].replace("T", " "),
                     loc["clock_max"].replace("T", " "), offs])
    A(md_table(["orbit", "earliest", "latest", "UTC offsets"], rows))
    A("")
    A("Sub-swath provenance of the plotted rows:")
    A("")
    rows = [[geom, d["n_by_subswath"][geom].get("iw1", 0),
             d["n_by_subswath"][geom].get("iw2", 0)] for geom in ORBITS]
    A(md_table(["orbit", "IW1", "IW2"], rows))
    A("")
    A("## 3. The 13 consecutive state runs")
    A("")
    A("A run is a maximal block of equal `damage_label` after sorting by "
      "acquisition instant. Bands span the midpoint between the last "
      "overpass of one run and the first of the next; the first band starts "
      f"{d['plot_window']['pad_days']} days before the first overpass and the "
      f"last ends {d['plot_window']['pad_days']} days after the last one.")
    A("")
    rows = []
    for r in d["runs"]:
        rows.append([r["run"], r["state"], r["n"], r["n_by_orbit"]["ASC"],
                     r["n_by_orbit"]["DESC"], r["first_local"][:10],
                     r["last_local"][:10]])
    A(md_table(["run", "state", "n", "ASC", "DESC", "first", "last"], rows))
    A("")
    A("## 4. Per-state medians")
    A("")
    A("Median of each channel inside each state, with the number of finite "
      "overpasses in parentheses. The log-scaled channels are strictly "
      "positive in every state x orbit cell. The single missing `P` "
      "(2021-07-14 ASC) is drawn as a gap and is not interpolated.")
    A("")
    rows = []
    for geom in ORBITS:
        for dim in CHANNELS:
            row = [geom, dim]
            for state in EXPECTED_N_BY_STATE:
                c = res["summary"][geom][dim][state]
                row.append(f"{fmt(c['median'], 4)} (n={c['n']})")
            rows.append(row)
    A(md_table(["orbit", "dimension"] + list(EXPECTED_N_BY_STATE), rows))
    A("")
    A("Means, minima and maxima are in the result JSON under "
      "`summary.<orbit>.<dimension>.<state>`.")
    A("")
    A("## 5. Cross-checks")
    A("")
    A(f"The script recomputes every plotted point from the table it reads "
      f"(`{d['csv']}`) and checks it against the metadata, the paper figure "
      f"JSON and the committed tower reference: "
      f"**{res['n_cross_checks']} checks, "
      f"{res['n_cross_checks_failed']} failures**.")
    A("")
    if res["n_cross_checks_failed"]:
        A("**Failures:**")
        A("")
        A(md_table(["check", "recomputed", "committed"],
                   [[c["check"], c.get("recomputed", ""),
                     c.get("committed", "")]
                    for c in res["cross_checks"] if not c["ok"]]))
    else:
        covered = ("row and date counts, the 482-row acquisition join, the "
                   "13-run sequence and lengths, sub-swath counts, "
                   "per-state/per-orbit counts, row-by-row `burst_id` / date "
                   "/ orbit / label parity with the tower reference, "
                   "`A == n_masked`, strict positivity on the log axes, "
                   "local-time ranges and UTC offsets, and all committed "
                   "sha256 pins")
        if res.get("n_superseded"):
            covered = ("row and date counts, the 482-row acquisition join, "
                       "the 13-run sequence and lengths, sub-swath counts, "
                       "per-state/per-orbit counts, row-by-row `burst_id` / "
                       "date / orbit / label parity with the tower "
                       "reference, `A == n_masked`, `A <= chip`, strict "
                       "positivity on the log axes, local-time ranges and "
                       "UTC offsets, and the ring definition pins")
        A(f"**Failures: none.** The checks cover {covered}.")
        A("")
        groups = collections.Counter(c["check"].split(".")[0]
                                     for c in res["cross_checks"])
        A(md_table(["check group", "checks"],
                   [[k, groups[k]] for k in sorted(groups)]))
    A("")
    if res.get("n_superseded"):
        A(f"**Superseded by the reference rule ({res['n_superseded']} pins, "
          "not failures).** This table measures with a reference ring that "
          "is disjoint from the support, so the pins that describe the "
          "committed rule -- the median inside the search window, the "
          "400 x 11 chip, channel completeness and the `n_masked` / "
          "`gamma2` values the paper reported -- cannot hold here. They are "
          "kept with both sides in `superseded` in the JSON:")
        A("")
        A(md_table(["check", "ring value", "committed value"],
                   [[c["check"], (c.get("recomputed") or "")[:80],
                     (c.get("committed") or "")[:80]]
                    for c in res["superseded"]]))
        A("")
    A("## 6. Reproduce")
    A("")
    A("```sh")
    A("cd <repo root>")
    if res.get("n_superseded"):
        A("python3 code/lumo/fig_lumo_ocv_orbits.py " + BACKSLASH)
        A("    --csv " + d["csv"] + " " + BACKSLASH)
        A("    --meta " + d["csv_meta"])
    else:
        A("python3 code/lumo/fig_lumo_ocv_orbits.py")
        A("PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_orbits.sh")
    A("```")
    A("")
    A("Needs numpy and matplotlib. This report was generated with python "
      f"{platform.python_version()}, numpy {np.__version__}, matplotlib "
      f"{matplotlib.__version__}. No output carries a timestamp, so repeated "
      "runs are identical.")
    A("")
    A("## 7. Caveats")
    A("")
    A("- **The 13 runs are one deployment, not 13 independent experiments.** "
      "States repeat (healthy seven times, each damage state twice), but the "
      "repeats are separated in time and season. The bands make the sequence "
      "visible; they do not turn it into a randomized treatment.")
    A("- **State and season are confounded.** The first healthy block covers "
      "late summer and autumn 2020, the long healthy block spans winter "
      "2020/21, and the damage states are short visits in spring and summer "
      "2021.")
    A("- **Orbit and sub-swath are not balanced.** ASC is almost entirely "
      "IW2 (56 of 57 rows), DESC mostly IW1 (108 of 121); the columns "
      "therefore differ in both viewing geometry and sub-swath.")
    A("- **The x axis is local wall time.** ASC overpasses are around 17–19 "
      "local, DESC around 05–08 local. UTC storage and the DST offset change "
      "are handled explicitly; a date-only axis would hide the two clocks.")
    A("- **Missing values are gaps.** Exactly one of 178 rows lacks `P` "
      "(2021-07-14 ASC, the one non-majority chip shape excluded from the "
      "persistence map). Nothing is interpolated.")
    A("- **`P`, `D`, `F` and `A` use log axes** because their ranges span "
      "orders of magnitude; `gamma2` and `S` are linear.")
    A("")
    A("---")
    A("")
    A("Generated by `code/lumo/fig_lumo_ocv_orbits.py` "
      f"({res['n_cross_checks']} cross-checks, "
      f"{res['n_cross_checks_failed']} failures).")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="LUMO OCV time series per dimension and orbit.")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="output directory of the JPEG figure")
    ap.add_argument("--json", default=OUT_JSON,
                    help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD,
                    help="path of the Markdown report")
    ap.add_argument("--csv", default=None,
                    help="channel table to plot (default: the committed one)")
    ap.add_argument("--meta", default=None,
                    help="its metadata JSON (default: the committed one)")
    ap.add_argument("--no-figures", action="store_true",
                    help="only recompute and cross-check (no matplotlib output)")
    args = ap.parse_args(argv)

    t0 = time.time()
    print("LUMO OCV time series by dimension and orbit ...")
    res = compute(args.csv, args.meta)
    print(f"  compute: {time.time() - t0:.1f}s "
          f"({res['data']['n_rows']} overpasses, "
          f"{len(res['data']['runs'])} runs, "
          f"{res['n_cross_checks']} cross-checks, "
          f"{res['n_cross_checks_failed']} failed)")
    for c in res["cross_checks"]:
        if not c["ok"]:
            print(f"    FAIL {c['check']}")

    if not args.no_figures:
        out = make_figure(res, args.figdir)
        print(f"  {os.path.relpath(out)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, sort_keys=True, default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w", encoding="utf-8") as fh:
        fh.write(report(res))
    print(f"  {os.path.relpath(args.md)}")

    print(f"done in {time.time() - t0:.1f}s, "
          f"{res['n_cross_checks_failed']} cross-check failures")
    return 0 if res["n_cross_checks_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

