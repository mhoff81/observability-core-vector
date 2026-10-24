#!/usr/bin/env python3
"""LUMO — physical dynamics against OCV behaviour.

The LUMO paper reports an OCV *behaviour* vector (gamma2, P) per damage state.
This companion asks a different question: do states that occupy similar
positions in an **independently measured physical** plane also occupy similar
positions in that behaviour plane?  Two panels, one state colour scheme:

  A  physical response      x = df_H1 / f_H1,0   y = dzeta_H1 / zeta_H1,0
  B  OCV behaviour          x = z(gamma2)        y = z(P)

Both planes start at the healthy state (0, 0).  In A the origin is the healthy
block of the *same* controlled campaign (the tower as it stood days before the
state), so the arrow is a genuine before/after vector; in B the origin is the
healthy median of the same orbit, and z is measured in healthy MADs of the log
channel.  Panel A has **one** point per state: each state was measured in one
campaign of five healthy plus five damaged recordings.  Panel B has **two**
visits per state (the two SAR state runs), so reproducibility can be assessed on
the OCV side only -- that asymmetry is the main caveat of the figure and is
stated in the report.

Panel B never pools the orbits: colour = state, marker shape = orbit (circle
ASC, square DESC), fill = visit (filled first, open repeat), marker area
proportional to n.  A thin dashed arrow shows the state-level (visit-pooled)
OCV vector, which is the reproducible object; the visit markers show how much
the individual visits scatter.

Primary analysis (descriptive, no p-values): the direction of each vector, its
magnitude, visit-to-visit repeatability (cosine and angle between the two visit
vectors), the state ordering by |v| in both planes, and the cosine between
v_phys and v_OCV per state and orbit.  A matching direction is evidence of
*consistency between two measurement channels* -- it is **not** evidence that
(gamma2, P) measures stiffness or damping.

Data.  Panel A reads ``data/lumo/reference/lumo_state_dynamics.json`` (built and
verified by ``lumo_state_dynamics_json.py``: H1 frequency shift from the Welch
mode consensus, H1 damping from the half-power bandwidth of 30 recordings).
Panel B reads the committed LUMO channel table -- by default the **support**
chip ``data/lumo/support/lumo_ocv_channels.csv``; ``--chip committed`` repeats
everything on ``data/lumo/lumo_ocv_channels.csv`` as a sensitivity, because the
two chips define gamma2 differently (the support re-cut raises gamma2 by an
order of magnitude, see ``data/lumo/support/README.md``) and the OCV direction
turns out to depend on it.  No output carries a wall-clock timestamp, so
repeated runs are byte-identical.

Usage:
  python3 code/lumo/fig_lumo_ocv_physics.py
  python3 code/lumo/fig_lumo_ocv_physics.py --chip committed
  python3 code/lumo/fig_lumo_ocv_physics.py --selftest
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import os
import random
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import lumo_ocv_core as core              # noqa: E402
import lumo_state_dynamics_json as dyn    # noqa: E402

SITE = "LUMO"
FIG_NAME = "fig_lumo_ocv_physics"
FIG_DPI = 630
OUT_JSON = os.path.join(core.DATA, FIG_NAME + ".json")
OUT_MD = os.path.join(core.FIGDIR, FIG_NAME + ".md")

DYNAMICS_JSON = os.path.join(core.REF, "lumo_state_dynamics.json")
CHIPS = {
    "support": os.path.join(core.DATA, "support", "lumo_ocv_channels.csv"),
    "committed": core.CSV_PATH,
}

STATES = ["DAM 3", "DAM 4", "DAM 6"]
ORBITS = ["ASC", "DESC"]
ORBIT_OF = {"ASC": "ASCENDING", "DESC": "DESCENDING"}
STATE_COLORS = {"DAM 3": "#DD8452", "DAM 4": "#C44E52", "DAM 6": "#55A868"}
ORBIT_MARKERS = {"ASC": "o", "DESC": "s"}
N_BOOT = 2000
BOOT_SEED = core.RNG_SEED
LOG_CHANNELS = ("gamma2", "P")

# Mathtext twins of the plain labels: markdown stays grepable ("gamma^2"),
# only the drawing calls read these.
FIG_GAMMA2 = r"$\gamma^2$"
FIG_ZETA = r"$\zeta$"
FIG_DF = r"$\Delta f_{\mathrm{H1}}/f_{\mathrm{H1},0}$"
FIG_DZETA = r"$\Delta\zeta_{\mathrm{H1}}/\zeta_{\mathrm{H1},0}$"
FIG_ZG = "z(" + FIG_GAMMA2 + ")"
FIG_ZP = "z(P)"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 8.4,
    "axes.labelsize": 8.2,
    "xtick.labelsize": 7.2,
    "ytick.labelsize": 7.2,
    "legend.fontsize": 7.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.22,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": FIG_DPI,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
MM = 1.0 / 25.4


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def fmt(x, digits=3):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{digits}f}"
    return str(x)


def pct(x, digits=1):
    if x is None:
        return ""
    return f"{100.0 * x:+.{digits}f} %"


def md_table(header, rows):
    def cell(value):
        return fmt(value).replace("|", r"\|")

    out = ["| " + " | ".join(cell(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for row in rows:
        out.append("| " + " | ".join(cell(c) for c in row) + " |")
    return "\n".join(out)


def check(checks, name, got, want):
    entry = {"check": name, "ok": bool(got == want)}
    if not entry["ok"]:
        entry["recomputed"] = repr(got)[:300]
        entry["committed"] = repr(want)[:300]
    checks.append(entry)
    return entry["ok"]


def norm(v):
    return math.hypot(v[0], v[1])


def angle_deg(v):
    """Direction in the plotted plane, degrees, atan2(y, x) on data units."""
    return math.degrees(math.atan2(v[1], v[0]))


def cosine(a, b):
    na, nb = norm(a), norm(b)
    if na == 0.0 or nb == 0.0:
        return None
    return (a[0] * b[0] + a[1] * b[1]) / (na * nb)


def angle_between(a, b):
    c = cosine(a, b)
    return None if c is None else math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _median(values):
    return float(np.median(values))


def bootstrap_ratio(healthy, state, n_boot=N_BOOT, seed=BOOT_SEED):
    """CI of (median(state) - median(healthy)) / median(healthy).

    Resampling both five-recording blocks independently is the honest statement
    for n = 5 per block: it carries the recording-to-recording spread, not a
    model.
    """
    rng = random.Random(seed)
    nh, ns = len(healthy), len(state)
    out = []
    for _ in range(n_boot):
        bh = [healthy[rng.randrange(nh)] for _ in range(nh)]
        bs = [state[rng.randrange(ns)] for _ in range(ns)]
        mh = _median(bh)
        out.append((_median(bs) - mh) / mh)
    out.sort()
    lo = out[int(0.025 * n_boot)]
    hi = out[int(0.975 * n_boot) - 1]
    return lo, hi


# ---------------------------------------------------------------------------
# Panel A - physical response, from the committed dynamics reference
# ---------------------------------------------------------------------------
def panel_a(doc):
    """Relative frequency and damping change per state, H1 primary."""
    out = []
    for s in doc["states"]:
        modes = {m["name"]: m for m in s["modes"]}
        h1 = modes["first_higher"]
        fund = modes["fundamental"]
        lo, hi = bootstrap_ratio(h1["zeta_healthy_recordings"],
                                 h1["zeta_state_recordings"])
        lo_f, hi_f = bootstrap_ratio(fund["zeta_healthy_recordings"],
                                     fund["zeta_state_recordings"])
        x = h1["df_rel"]
        y = h1["zeta_rel"]
        out.append({
            "state": s["state"],
            "healthy_campaign": s["healthy_campaign"],
            "damaged_campaign": s["damaged_campaign"],
            "healthy_window": s["healthy_window"],
            "state_window": s["state_window"],
            "df_rel": x,
            "dzeta_rel": y,
            "v": [x, y],
            "norm": norm([x, y]),
            "angle_deg": angle_deg([x, y]),
            "df_bin_uncertainty_rel": h1["df_bin_uncertainty_rel"],
            "dzeta_ci95": [lo, hi],
            "zeta_healthy": h1["zeta_healthy"],
            "zeta_state": h1["zeta_state"],
            "n_healthy_recordings": len(h1["zeta_healthy_recordings"]),
            "n_state_recordings": len(h1["zeta_state_recordings"]),
            "sensitivity_fundamental": {
                "df_rel": fund["df_rel"],
                "dzeta_rel": fund["zeta_rel"],
                "v": [fund["df_rel"], fund["zeta_rel"]],
                "norm": norm([fund["df_rel"], fund["zeta_rel"]]),
                "angle_deg": angle_deg([fund["df_rel"], fund["zeta_rel"]]),
                "dzeta_ci95": [lo_f, hi_f],
            },
        })
    return out


# ---------------------------------------------------------------------------
# Panel B - OCV behaviour, from a channel table
# ---------------------------------------------------------------------------
def runs_of(rows):
    """Consecutive ``damage_label`` runs in acquisition order."""
    ordered = sorted(rows, key=lambda r: (r["date"], r["burst_id"]))
    runs = []
    for r in ordered:
        if not runs or r["damage_label"] != runs[-1]["state"]:
            runs.append({"state": r["damage_label"], "rows": []})
        runs[-1]["rows"].append(r)
    for i, run in enumerate(runs, 1):
        run["run"] = i
        run["n"] = len(run["rows"])
        run["first_date"] = run["rows"][0]["date"]
        run["last_date"] = run["rows"][-1]["date"]
    return ordered, runs


def visits_of(runs):
    """The state runs of each damage state, in acquisition order."""
    out = collections.OrderedDict((s, []) for s in STATES)
    for run in runs:
        if run["state"] in out:
            out[run["state"]].append(run)
    return out


def healthy_reference(rows, orbit):
    """Robust centre and scale of the log channels of the healthy rows."""
    ref = [r for r in rows
           if r["orbit"] == ORBIT_OF[orbit] and r["damage_label"] == "healthy"]
    out = {}
    for ch in LOG_CHANNELS:
        vals = [math.log(r[ch]) for r in ref if r[ch] and r[ch] > 0]
        centre = _median(vals)
        scale = _median([abs(v - centre) for v in vals])
        out[ch] = {"n": len(vals), "centre_log": centre, "scale_log_mad": scale}
    return out


def zvector(sel, stats):
    """z of the median log channel against the healthy median, in healthy MADs."""
    out = {}
    for ch in LOG_CHANNELS:
        vals = [math.log(r[ch]) for r in sel if r[ch] and r[ch] > 0]
        st = stats[ch]
        z = (float(np.median(vals)) - st["centre_log"]) / st["scale_log_mad"]
        out[ch] = {"n": len(vals), "z": z,
                   "median_log": float(np.median(vals))}
    v = [out["gamma2"]["z"], out["P"]["z"]]
    return {
        "v": v,
        "norm": norm(v),
        "angle_deg": angle_deg(v),
        "n_dates": len(sel),
        "channels": out,
        "thin": min(out["gamma2"]["n"], out["P"]["n"]) < 3,
    }


def panel_b(rows):
    """State-level and visit-level OCV vectors per orbit, orbits never pooled."""
    ordered, runs = runs_of(rows)
    by_state = visits_of(runs)
    out = {"n_rows": len(ordered), "n_runs": len(runs),
           "runs": [{"run": r["run"], "state": r["state"], "n": r["n"],
                     "first": r["first_date"], "last": r["last_date"]}
                    for r in runs],
           "orbits": {}}
    for orbit in ORBITS:
        stats = healthy_reference(rows, orbit)
        block = {"healthy": stats, "states": {}}
        for state in STATES:
            pooled_rows = [r for r in ordered
                           if r["damage_label"] == state
                           and r["orbit"] == ORBIT_OF[orbit]]
            visits = []
            for idx, run in enumerate(by_state[state], 1):
                sel = [r for r in run["rows"] if r["orbit"] == ORBIT_OF[orbit]]
                entry = zvector(sel, stats)
                entry.update({"visit": idx, "first": run["first_date"],
                              "last": run["last_date"]})
                visits.append(entry)
            entry = zvector(pooled_rows, stats)
            entry.update({"state": state, "n_visits": len(visits),
                          "visits": visits})
            block["states"][state] = entry
        out["orbits"][orbit] = block
    return out


# ---------------------------------------------------------------------------
# Analysis - directions, repeatability, ordering, cross-plane agreement
# ---------------------------------------------------------------------------
def analysis(pa, pb):
    phys = {e["state"]: e for e in pa}
    out = {"repeatability": {}, "ordering": {}, "cross": {}}
    for orbit in ORBITS:
        st_block = pb["orbits"][orbit]["states"]

        rep = {}
        for state in STATES:
            visits = st_block[state]["visits"]
            entry = {"n_visits": len(visits),
                     "visit_n": [v["n_dates"] for v in visits],
                     "visit_v": [v["v"] for v in visits]}
            if len(visits) >= 2:
                a, b = visits[0]["v"], visits[1]["v"]
                entry["cosine"] = cosine(a, b)
                entry["angle_deg"] = angle_between(a, b)
                entry["thin"] = bool(visits[0]["thin"] or visits[1]["thin"])
            else:
                entry["cosine"] = None
                entry["angle_deg"] = None
                entry["thin"] = True
            rep[state] = entry
        out["repeatability"][orbit] = rep

        order = {
            "physical_by_norm": sorted(STATES, key=lambda s: phys[s]["norm"]),
            "ocv_by_norm": sorted(STATES, key=lambda s: st_block[s]["norm"]),
        }
        order["same_order"] = (order["physical_by_norm"]
                               == order["ocv_by_norm"])
        out["ordering"][orbit] = order

        cross = {}
        for state in STATES:
            v_a = phys[state]["v"]
            v_b = st_block[state]["v"]
            cross[state] = {
                "cosine": cosine(v_a, v_b),
                "angle_deg": angle_between(v_a, v_b),
                "norm_physical": phys[state]["norm"],
                "norm_ocv": st_block[state]["norm"],
                "n_ocv_dates": st_block[state]["n_dates"],
            }
        pairs = []
        for i, s1 in enumerate(STATES):
            for s2 in STATES[i + 1:]:
                a_phys = angle_between(phys[s1]["v"], phys[s2]["v"])
                a_ocv = angle_between(st_block[s1]["v"], st_block[s2]["v"])
                pairs.append({"pair": [s1, s2],
                              "angle_physical": a_phys,
                              "angle_ocv": a_ocv,
                              "abs_diff": (None if a_phys is None
                                           or a_ocv is None
                                           else abs(a_phys - a_ocv))})
        cross["pairwise"] = pairs
        out["cross"][orbit] = cross
    return out


# ---------------------------------------------------------------------------
# Compute and cross-check
# ---------------------------------------------------------------------------
EXPECTED_N_BY_STATE = {"healthy": 121, "DAM 3": 25, "DAM 4": 16, "DAM 6": 16}
EXPECTED_HEALTHY_BY_ORBIT = {"ASC": 40, "DESC": 81}
EXPECTED_VISIT_N = {
    "ASC": {"DAM 3": [5, 3], "DAM 4": [3, 2], "DAM 6": [2, 2]},
    "DESC": {"DAM 3": [12, 5], "DAM 4": [5, 6], "DAM 6": [6, 6]},
}
EXPECTED_RUNS = 13


def compute(chip="support", csv_path=None):
    path = csv_path or CHIPS[chip]
    other = "committed" if chip == "support" else "support"
    with open(DYNAMICS_JSON, encoding="utf-8") as fh:
        doc = json.load(fh)
    ref = dyn.verify(doc)
    checks = []
    for entry in ref["checks"]:
        new = {"check": "reference." + entry["check"], "ok": entry["ok"]}
        if not entry["ok"]:
            new["recomputed"] = entry.get("recomputed", "")
            new["committed"] = entry.get("committed", "")
        checks.append(new)

    rows = core.load_csv(path)
    pa = panel_a(doc)
    pb = panel_b(rows)
    an = analysis(pa, pb)
    pb_other = panel_b(core.load_csv(CHIPS[other]))

    n_by_state = {s: sum(1 for r in rows if r["damage_label"] == s)
                  for s in ["healthy"] + STATES}
    n_healthy_orbit = {o: sum(1 for r in rows
                              if r["damage_label"] == "healthy"
                              and r["orbit"] == ORBIT_OF[o]) for o in ORBITS}
    visit_n = {o: {s: [v["n_dates"]
                       for v in pb["orbits"][o]["states"][s]["visits"]]
                   for s in STATES} for o in ORBITS}
    pooled_mismatch = []
    for o in ORBITS:
        for s in STATES:
            block = pb["orbits"][o]["states"][s]
            if block["n_dates"] != sum(v["n_dates"] for v in block["visits"]):
                pooled_mismatch.append((o, s))

    check(checks, "data.n_rows", len(rows), 178)
    check(checks, "data.n_runs", pb["n_runs"], EXPECTED_RUNS)
    check(checks, "data.n_by_state", n_by_state, EXPECTED_N_BY_STATE)
    check(checks, "data.healthy_by_orbit", n_healthy_orbit,
          EXPECTED_HEALTHY_BY_ORBIT)
    check(checks, "data.n_visits_per_state",
          {s: pb["orbits"]["ASC"]["states"][s]["n_visits"] for s in STATES},
          {s: 2 for s in STATES})
    check(checks, "data.visit_n_by_orbit", visit_n, EXPECTED_VISIT_N)
    check(checks, "data.pooled_equals_visit_sum", pooled_mismatch, [])
    check(checks, "panel_a.n_states", len(pa), 3)
    check(checks, "panel_a.primary_mode",
          doc["params"]["primary_mode"], "first_higher")
    check(checks, "panel_b.orbits_separate", sorted(pb["orbits"]), ORBITS)
    check(checks, "panel_b.healthy_scale_positive",
          all(v["scale_log_mad"] > 0
              for o in ORBITS for v in pb["orbits"][o]["healthy"].values()), True)
    check(checks, "panel_b.visits_kept_per_orbit",
          all(pb["orbits"][o]["states"][s]["n_visits"] == 2 for o in ORBITS
              for s in STATES), True)
    differs = sum(1 for o in ORBITS for s in STATES
                  if pb["orbits"][o]["states"][s]["v"]
                  != pb_other["orbits"][o]["states"][s]["v"])
    check(checks, "sensitivity.chip_rows_equal", pb["n_rows"], pb_other["n_rows"])
    check(checks, "sensitivity.n_state_vectors_changed", differs, 6)

    n_ok = sum(1 for c in checks if c["ok"])
    return {
        "site": SITE,
        "figure": FIG_NAME,
        "chip": chip,
        "csv": path,
        "csv_sha256": core.sha256(path),
        "reference": {"path": DYNAMICS_JSON,
                      "n_checks": ref["n_checks"],
                      "n_failed": ref["n_failed"]},
        "panel_a": pa,
        "panel_b": pb,
        "analysis": an,
        "sensitivity": {"chip": other, "csv": CHIPS[other],
                        "panel_b": pb_other},
        "cross_checks": checks,
        "n_cross_checks": len(checks),
        "n_cross_checks_failed": len(checks) - n_ok,
    }


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def _pad(values, frac=0.18, floor=0.6):
    lo, hi = min(values), max(values)
    span = max(hi - lo, floor)
    return lo - frac * span, hi + frac * span


def _draw_panel_a(ax, pa):
    ax.axhline(0.0, color="0.45", lw=0.7, ls=":")
    ax.axvline(0.0, color="0.45", lw=0.7, ls=":")
    for e in pa:
        col = STATE_COLORS[e["state"]]
        x, y = 100.0 * e["df_rel"], 100.0 * e["dzeta_rel"]
        ax.annotate("", xy=(x, y), xytext=(0.0, 0.0),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.25,
                                    shrinkA=0.0, shrinkB=0.0), zorder=4)
        ax.errorbar([x], [y], xerr=[100.0 * e["df_bin_uncertainty_rel"]],
                    yerr=[[y - 100.0 * e["dzeta_ci95"][0]],
                          [100.0 * e["dzeta_ci95"][1] - y]],
                    fmt="none", ecolor=col, elinewidth=0.8, capsize=2.2,
                    alpha=0.85, zorder=3)
        ax.plot([x], [y], marker="o", ms=5.2, mfc=col, mec="white", mew=0.7,
                zorder=5)
    ax.text(0.015, 0.02, "healthy\n(0, 0)", fontsize=6.6, color="0.35",
            ha="left", va="bottom", transform=ax.transAxes)
    ax.set_xlabel(FIG_DF + "  [%]  (signed; loss to the left)")
    ax.set_ylabel(FIG_DZETA + "  [%]")
    x_lo, x_hi = _pad([100.0 * e["df_rel"] for e in pa] + [0.0])
    ax.set_xlim(x_lo, x_hi + 0.42 * (x_hi - x_lo))
    ys = ([100.0 * e["dzeta_rel"] for e in pa]
          + [100.0 * e["dzeta_ci95"][0] for e in pa]
          + [100.0 * e["dzeta_ci95"][1] for e in pa] + [0.0])
    ax.set_ylim(*_pad(ys, frac=0.10))
    ax.set_title("(A)  physical response — H1 mode, one campaign per state",
                 loc="left")
    ax.text(0.015, 0.975,
            "arrows: healthy → state (same campaign)\n"
            "bars: 95 % bootstrap over 5+5 recordings\n"
            "x: ±0.125 Hz mode-grid quantisation",
            fontsize=6.3, color="0.3", ha="left", va="top",
            transform=ax.transAxes,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=1.0))


def _draw_panel_b(ax, pb, chip):
    ax.axhline(0.0, color="0.45", lw=0.7, ls=":")
    ax.axvline(0.0, color="0.45", lw=0.7, ls=":")
    xs, ys = [0.0], [0.0]
    for orbit in ORBITS:
        marker = ORBIT_MARKERS[orbit]
        for state in STATES:
            col = STATE_COLORS[state]
            block = pb["orbits"][orbit]["states"][state]
            for visit in block["visits"]:
                x, y = visit["v"]
                ms = 3.0 + 2.1 * math.sqrt(max(visit["n_dates"], 1))
                filled = visit["visit"] == 1
                ax.plot([x], [y], linestyle="none", marker=marker, ms=ms,
                        mfc=(col if filled else "white"), mec=col,
                        mew=(0.7 if filled else 1.0),
                        alpha=(0.95 if filled else 0.9), zorder=4)
                if visit.get("thin"):
                    ax.plot([x], [y], linestyle="none", marker="o",
                            ms=ms + 4.4, mfc="none", mec="0.55", mew=0.5,
                            zorder=3)
                xs.append(x)
                ys.append(y)
            px, py = block["v"]
            ax.annotate("", xy=(px, py), xytext=(0.0, 0.0),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=0.9,
                                        ls="--", shrinkA=0.0, shrinkB=0.0,
                                        alpha=0.75), zorder=2)
            xs.append(px)
            ys.append(py)
    ax.set_xlabel(FIG_ZG + "   (healthy median + MAD of the log channel)")
    ax.set_ylabel(FIG_ZP)
    x_lo, x_hi = _pad(xs)
    y_lo, y_hi = _pad(ys)
    ax.set_xlim(x_lo, x_hi + 0.30 * (x_hi - x_lo))
    ax.set_ylim(y_lo, y_hi)
    ax.set_title(f"(B)  OCV behaviour — {chip} chip, visits separate, "
                 "orbits not pooled", loc="left")
    ax.text(0.985, 0.03,
            "marker area ∝ n   ·   grey ring: n < 3\n"
            "dashed arrow: state-level (visit-pooled)",
            fontsize=6.3, color="0.3", ha="right", va="bottom",
            transform=ax.transAxes,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=1.0))


def _place_labels(ax, fig, items, offsets, bbox=True):
    """Place annotations without overlap, deterministically.

    ``items`` is a list of ``(text, x, y, colour)``; each label is tried at the
    candidate offsets in order and kept at the first position that (a) does not
    overlap an already placed label and (b) stays inside the axes.  The order is
    fixed, so repeated runs place the same labels in the same spots.
    """
    placed = []
    axbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    for text, x, y, colour in items:
        chosen = None
        for dx, dy in offsets:
            ann = ax.annotate(
                text, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                fontsize=6.2, ha=("left" if dx > 0 else "right"), va="center",
                color=colour, zorder=7,
                bbox=(dict(facecolor="white", edgecolor="none", alpha=0.80,
                           pad=0.7) if bbox else None),
                arrowprops=dict(arrowstyle="-", color=colour, lw=0.45,
                                shrinkA=0.0, shrinkB=2.0))
            fig.canvas.draw()
            bb = ann.get_window_extent(renderer=fig.canvas.get_renderer())
            inside = (bb.x0 >= axbox.x0 and bb.x1 <= axbox.x1
                      and bb.y0 >= axbox.y0 and bb.y1 <= axbox.y1)
            if inside and all(not bb.overlaps(p) for p in placed):
                chosen = bb
                break
            ann.remove()
        if chosen is None:
            ann = ax.annotate(
                text, xy=(x, y), xytext=offsets[0], textcoords="offset points",
                fontsize=6.2, ha="left", va="center", color=colour, zorder=7)
            fig.canvas.draw()
            chosen = ann.get_window_extent(renderer=fig.canvas.get_renderer())
        placed.append(chosen)


def _label_panel_a(ax, pa):
    """Short state tags next to the points (numbers live in the caption/table)."""
    positions = {"DAM 3": (-46, 9, "right"), "DAM 4": (-46, -9, "right"),
                 "DAM 6": (10, 0, "left")}
    for e in pa:
        dx, dy, ha = positions[e["state"]]
        ax.annotate(
            e["state"].replace(" ", ""),
            xy=(100.0 * e["df_rel"], 100.0 * e["dzeta_rel"]),
            xytext=(dx, dy), textcoords="offset points", fontsize=6.8,
            color=STATE_COLORS[e["state"]], ha=ha, va="center", zorder=7,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.6))


def _label_panel_b(ax, pb):
    items = []
    for orbit in ORBITS:
        for state in STATES:
            px, py = pb["orbits"][orbit]["states"][state]["v"]
            items.append((f"{state.replace(' ', '')} {orbit}", px, py,
                          STATE_COLORS[state]))
    offsets = [(6, 8), (6, -10), (-6, 8), (-6, -10),
               (6, 20), (6, -22), (-6, 20), (-6, -22),
               (16, 0), (-16, 0), (6, 32), (6, -34)]
    _place_labels(ax, ax.figure, items, offsets, bbox=False)


def make_figure(res, figdir):
    pa, pb, chip = res["panel_a"], res["panel_b"], res["chip"]
    fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 108 * MM))
    _draw_panel_a(axes[0], pa)
    _draw_panel_b(axes[1], pb, chip)
    fig.tight_layout(rect=(0, 0.055, 1, 0.895))
    _label_panel_a(axes[0], pa)
    _label_panel_b(axes[1], pb)

    handles = [Patch(facecolor=STATE_COLORS[s], edgecolor="none",
                     label=s.replace(" ", "")) for s in STATES]
    handles += [Line2D([], [], color="0.35", lw=1.25, label="physical arrow"),
                Line2D([], [], color="0.35", lw=0.9, ls="--",
                       label="OCV pooled arrow"),
                Line2D([], [], color="0.35", marker="o", ls="none", ms=5,
                       label="ASC visit"),
                Line2D([], [], color="0.35", marker="s", ls="none", ms=5,
                       label="DESC visit"),
                Line2D([], [], color="0.35", marker="o", ls="none", ms=5,
                       mfc="white", label="repeat visit"),
                Line2D([], [], color="0.35", marker="o", ls="none", ms=3.4,
                       label="n = 2"),
                Line2D([], [], color="0.35", marker="o", ls="none", ms=8.2,
                       label="n = 12")]
    fig.legend(handles=handles, loc="lower center", ncol=8, frameon=False,
               bbox_to_anchor=(0.5, 0.002), handletextpad=0.4,
               columnspacing=1.1)
    fund_txt = " / ".join(f"{100.0 * e['sensitivity_fundamental']['dzeta_rel']:+.1f}"
                          for e in pa)
    fig.suptitle(
        "LUMO lattice tower — do states that sit together in the physical "
        "plane also sit together in the OCV plane?\n"
        "A: relative H1 frequency change vs relative H1 damping change, one "
        "campaign per state (5 healthy + 5 damaged recordings)\n"
        "B: OCV behaviour " + FIG_GAMMA2 + " vs P, two SAR visits per state, "
        "z against the healthy median of the same orbit, orbits never pooled\n"
        "fundamental-mode " + FIG_DZETA + " [" + fund_txt + "] % is the "
        "sensitivity, not the axis; marker area ∝ n; the visit-pooled arrow is "
        "the reproducible object",
        fontsize=8.0, y=0.995)
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, FIG_NAME + ".jpg")
    fig.savefig(out, dpi=FIG_DPI, facecolor="white")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------
def report(res):
    pa = {e["state"]: e for e in res["panel_a"]}
    pb = res["panel_b"]
    an = res["analysis"]
    chip = res["chip"]
    other = res["sensitivity"]["chip"]
    L = []
    A = L.append

    A("# LUMO — physical dynamics against OCV behaviour")
    A("")
    A("Two planes, same three states.  Panel A is the **independently measured "
      "physical** response of the tower, panel B the **OCV behaviour** of the "
      "same states.  The question is whether states that sit close together in "
      "A also sit close together in B.")
    A("")
    A("```")
    A("v_phys = [ df_H1 / f_H1,0 ,  dzeta_H1 / zeta_H1,0 ]")
    A("v_OCV  = [ z(gamma2)      ,  z(P)               ]")
    A("```")
    A("")
    A("`f_H1` is the first higher mode of the Welch mode consensus, `zeta_H1` "
      "its half-power bandwidth damping, both measured in the controlled "
      "campaign of the state (five healthy recordings, then five damaged).  "
      "`z` is the median log channel of the state against the healthy median of "
      "the **same orbit**, in healthy MADs of the log channel.  The origin of "
      "both panels is the corresponding healthy state.")
    A("")
    A(f"Inputs: `{os.path.relpath(res['csv'])}` ({chip} chip, "
      f"{pb['n_rows']} overpasses, sha256 `{res['csv_sha256'][:12]}…`) and "
      f"`{os.path.relpath(res['reference']['path'])}` "
      f"({res['reference']['n_checks']} reference checks, "
      f"{res['reference']['n_failed']} failures).")
    A("")

    A("## 1. The asymmetry that shapes the figure")
    A("")
    A("| plane | objects | replicate |")
    A("| --- | --- | --- |")
    A("| A physical | 3 states | **none** — one campaign per state (5 + 5 "
      "recordings) |")
    A(f"| B OCV | 3 states x 2 orbits x 2 visits | 2 SAR visits ({pb['n_runs']} "
      "consecutive state runs in the record) |")
    A("")
    A("So visit-to-visit repeatability can only be assessed on the OCV side.  "
      "Panel A carries a 95 % bootstrap interval over the ten recordings "
      "instead, and its arrows are drawn only for the unambiguous "
      "healthy → state pairing of one campaign.")
    A("")

    A("## 2. Panel A — physical response")
    A("")
    rows = []
    for s in res["panel_a"]:
        rows.append([s["state"], s["healthy_campaign"], s["damaged_campaign"],
                     pct(s["df_rel"], 2), pct(s["dzeta_rel"], 2),
                     f"[{pct(s['dzeta_ci95'][0], 1)}; "
                     f"{pct(s['dzeta_ci95'][1], 1)}]",
                     f"{s['norm']:.3f}", f"{s['angle_deg']:.1f}",
                     pct(s["sensitivity_fundamental"]["dzeta_rel"], 1)])
    A(md_table(["state", "healthy", "damaged", "df/f0", "dzeta/zeta0",
                "95 % CI (damping)", "|v|", "angle [deg]", "fund. dzeta"], rows))
    A("")
    A("`angle` is `atan2(dzeta/zeta0, df/f0)` in the plotted plane; the two "
      "axes are both relative changes, so the direction is scale free.")
    A("")

    A("## 3. Panel B — OCV behaviour per orbit (never pooled)")
    A("")
    rows = []
    for orbit in ORBITS:
        for state in STATES:
            b = pb["orbits"][orbit]["states"][state]
            rows.append([orbit, state, b["n_dates"],
                         f"{b['v'][0]:+.2f}", f"{b['v'][1]:+.2f}",
                         f"{b['norm']:.2f}", f"{b['angle_deg']:.1f}",
                         f"{b['channels']['gamma2']['n']}/"
                         f"{b['channels']['P']['n']}"])
    A(md_table(["orbit", "state", "n", "z(gamma2)", "z(P)", "|v|",
                "angle [deg]", "n gamma2/P"], rows))
    A("")
    A("Healthy reference (median and MAD of the log channel):")
    A("")
    rows = []
    for orbit in ORBITS:
        h = pb["orbits"][orbit]["healthy"]
        rows.append([orbit, h["gamma2"]["n"], f"{h['gamma2']['centre_log']:.4f}",
                     f"{h['gamma2']['scale_log_mad']:.4f}",
                     h["P"]["n"], f"{h['P']['centre_log']:.4f}",
                     f"{h['P']['scale_log_mad']:.4f}"])
    A(md_table(["orbit", "n gamma2", "median log gamma2", "MAD log gamma2",
                "n P", "median log P", "MAD log P"], rows))
    A("")

    A("### Visit-level points")
    A("")
    rows = []
    for orbit in ORBITS:
        for state in STATES:
            for v in pb["orbits"][orbit]["states"][state]["visits"]:
                rows.append([orbit, state, v["visit"], f"{v['first']}",
                             f"{v['last']}", v["n_dates"],
                             f"{v['v'][0]:+.2f}", f"{v['v'][1]:+.2f}",
                             "yes" if v["thin"] else ""])
    A(md_table(["orbit", "state", "visit", "first", "last", "n", "z(gamma2)",
                "z(P)", "n < 3"], rows))
    A("")

    A("## 4. Visit-to-visit repeatability (OCV side only)")
    A("")
    rows = []
    for orbit in ORBITS:
        for state in STATES:
            r = an["repeatability"][orbit][state]
            rows.append([orbit, state, r["visit_n"][0], r["visit_n"][1],
                         f"{r['cosine']:+.2f}",
                         f"{r['angle_deg']:.1f}",
                         "yes" if r["thin"] else ""])
    A(md_table(["orbit", "state", "n visit 1", "n visit 2", "cos(v1, v2)",
                "angle [deg]", "thin"], rows))
    A("")
    good = [(o, s) for o in ORBITS for s in STATES
            if an["repeatability"][o][s]["cosine"] is not None
            and an["repeatability"][o][s]["cosine"] > 0.9]
    bad = [(o, s) for o in ORBITS for s in STATES
           if an["repeatability"][o][s]["cosine"] is not None
           and an["repeatability"][o][s]["cosine"] < 0.5]
    A(f"Cosine above 0.9: {', '.join(f'{s} {o}' for o, s in good) or 'none'}.  "
      f"Cosine below 0.5: {', '.join(f'{s} {o}' for o, s in bad) or 'none'}.")
    A("")
    A("The two visits are months to a year apart, so a low cosine is a statement "
      "about the **record**, not about a laboratory repeat: season, wind "
      "climate, chip shape and sub-swath all change between the visits.  Where "
      "n is 2 or 3 per visit the point is drawn with a grey ring and its cosine "
      "carries almost no information.")
    A("")

    A("## 5. Ordering and cross-plane agreement (descriptive)")
    A("")
    rows = []
    for orbit in ORBITS:
        o = an["ordering"][orbit]
        rows.append([orbit,
                     " < ".join(s.replace(" ", "") for s in o["physical_by_norm"]),
                     " < ".join(s.replace(" ", "") for s in o["ocv_by_norm"]),
                     "yes" if o["same_order"] else "**no**"])
    A(md_table(["orbit", "order by |v_phys|", "order by |v_OCV|",
                "same order"], rows))
    A("")
    A("Per state, the cosine between the physical and the OCV vector:")
    A("")
    rows = []
    for orbit in ORBITS:
        for state in STATES:
            c = an["cross"][orbit][state]
            rows.append([orbit, state, f"{c['cosine']:+.2f}",
                         f"{c['angle_deg']:.1f}",
                         f"{c['norm_physical']:.3f}", f"{c['norm_ocv']:.2f}",
                         c["n_ocv_dates"]])
    A(md_table(["orbit", "state", "cos(v_phys, v_OCV)", "angle [deg]",
                "|v_phys|", "|v_OCV|", "n OCV"], rows))
    A("")
    A("The informative comparison is not this cosine but the **pairwise angle "
      "pattern**, because the three physical vectors are nearly parallel:")
    A("")
    rows = []
    for orbit in ORBITS:
        for p in an["cross"][orbit]["pairwise"]:
            rows.append([orbit,
                         f"{p['pair'][0].replace(' ', '')} vs "
                         f"{p['pair'][1].replace(' ', '')}",
                         f"{p['angle_physical']:.1f}",
                         f"{p['angle_ocv']:.1f}", f"{p['abs_diff']:.1f}"])
    A(md_table(["orbit", "pair", "angle physical [deg]", "angle OCV [deg]",
                "|difference|"], rows))
    A("")
    A("**Answer to the figure question, as far as three states can carry it:** "
      "**no.**  The three physical vectors are nearly parallel (1.0°, 7.5° and "
      "8.5° apart), so panel A orders the states by distance along one "
      "direction; panel B orders them differently (DAM4 closest to healthy in "
      "both orbits), its visit-to-visit repeatability is good only where n is "
      "largest, and its vectors change sign with the chip definition.  What the "
      "data support is a monotone **dose** response in the physical plane and a "
      "state- and orbit-dependent OCV signature — not a shared geometry.")
    A("")
    A("`cos(v_phys, v_OCV)` is reported for completeness.  It is a "
      "**descriptive** statistic on three points per orbit: the physical "
      "vectors point almost the same way, so the cosine is essentially the "
      "direction of the OCV vector, and a small angle between the two planes is "
      "therefore not evidence of a shared mechanism.  A matching direction "
      "would support *consistency between two measurement channels*; it would "
      "**not** show that (gamma2, P) measures stiffness or damping.")
    A("")

    A("## 6. Sensitivity — chip definition and damping mode")
    A("")
    A(f"Panel B is drawn on the **{chip}** chip.  The two chips define gamma2 "
      "differently: the support re-cut measures the echo cells inside the dwell "
      "support against a disjoint ring reference and raises gamma2 by an order "
      "of magnitude (see `data/lumo/support/README.md`).  Repeating the whole "
      f"panel on the **{other}** chip:")
    A("")
    rows = []
    for orbit in ORBITS:
        for state in STATES:
            a = pb["orbits"][orbit]["states"][state]
            b = res["sensitivity"]["panel_b"]["orbits"][orbit]["states"][state]
            rows.append([orbit, state,
                         f"{a['v'][0]:+.2f} / {a['v'][1]:+.2f}",
                         f"{b['v'][0]:+.2f} / {b['v'][1]:+.2f}",
                         f"{angle_between(a['v'], b['v']):.1f}"])
    A(md_table(["orbit", "state", f"z ({chip})", f"z ({other})",
                "angle between [deg]"], rows))
    A("")
    A("The chip definition moves the OCV points by tens of degrees and in "
      "places reverses them, so the OCV plane is **definition dependent** at "
      "this sample size, and any statement that joins the two planes inherits "
      "that dependence.")
    A("")
    rows = []
    for s in res["panel_a"]:
        f = s["sensitivity_fundamental"]
        rows.append([s["state"], pct(s["df_rel"], 2), pct(s["dzeta_rel"], 2),
                     f"{s['angle_deg']:.1f}", pct(f["df_rel"], 2),
                     pct(f["dzeta_rel"], 2), f"{f['angle_deg']:.1f}"])
    A(md_table(["state", "df/f0 (H1)", "dzeta H1", "angle H1",
                "df/f0 (fund.)", "dzeta fund.", "angle fund."], rows))
    A("")
    A("The fundamental mode barely moves in frequency but changes its damping "
      "in the **opposite direction** for DAM4 and DAM6.  Because its frequency "
      "shift is 0.00 % on the 0.25 Hz grid, its direction is ±90°: the damping "
      "axis alone, which is exactly why it is a sensitivity and not the primary "
      "pair.  H1 is the primary axis only because it is the mode that supplies "
      "the frequency shift.")
    A("")

    A("## 7. Cross-checks")
    A("")
    A(f"The script recomputes every plotted number and checks it: "
      f"**{res['n_cross_checks']} checks, "
      f"{res['n_cross_checks_failed']} failures** "
      f"({res['reference']['n_checks']} of them verify the dynamics reference "
      "itself: medians, spreads, relative changes, campaign names and the "
      "agreement with the committed frequency reference).")
    A("")
    if res["n_cross_checks_failed"]:
        A("**Failures:**")
        A("")
        A(md_table(["check", "recomputed", "committed"],
                   [[c["check"], c.get("recomputed", ""),
                     c.get("committed", "")]
                    for c in res["cross_checks"] if not c["ok"]]))
    else:
        groups = collections.Counter(c["check"].split(".")[0]
                                     for c in res["cross_checks"])
        A(md_table(["check group", "checks"],
                   [[k, groups[k]] for k in sorted(groups)]))
    A("")

    A("## 8. Caveats")
    A("")
    A("- **No physical replicate.** Each state was measured once (5 healthy + 5 "
      "damaged recordings in one campaign).  The damping interval is a "
      "bootstrap over those ten recordings, not a repeat of the experiment.")
    A("- **The physical plane is nearly one-dimensional.** DAM3 and DAM4 sit "
      "almost on top of each other (identical 1 Hz H1 shift on the 0.25 Hz "
      "grid, damping differing by 1.5 percentage points), so the plane cannot "
      "separate them; the ordering by |v| rests on DAM6 being about twice as "
      "far out.")
    A("- **The H1 frequency shift is quantised.** The mode consensus resolves "
      "0.25 Hz; the drawn x error is ±0.125 Hz (±0.9 %), and the identical "
      "−7.41 % of DAM3 and DAM4 is the same bin difference.")
    A("- **Thin visit points.** ASC visits carry 2–5 dates each; only the DESC "
      "visits reach n = 5–12.  Points with n < 3 are ringed.")
    A("- **The OCV plane depends on the chip definition** (section 6), and the "
      "two orbits are not comparable: they differ in viewing geometry and "
      "sub-swath, and are therefore never pooled.")
    A("- **Damping mode matters.** The H1 damping falls for all three states, "
      "the fundamental-mode damping rises for DAM4 and DAM6 (section 6).  The "
      "second-higher mode is carried in the reference as context only; its "
      "DAM3 value (a factor of 13) is implausible and shows how unstable that "
      "estimate is.")
    A("- **State and season are confounded**, as in the rest of the LUMO "
      "record: the visits are separated by months to a year.")
    A("")
    A("---")
    A("")
    A("Generated by `code/lumo/fig_lumo_ocv_physics.py` "
      f"({res['n_cross_checks']} cross-checks, "
      f"{res['n_cross_checks_failed']} failures).")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# Selftest and command line
# ---------------------------------------------------------------------------
def selftest():
    """Check the geometry helpers on cases with a known answer."""
    fails = []

    def eq(name, got, want, tol=1e-9):
        ok = (got is None and want is None) or (
            got is not None and want is not None and abs(got - want) <= tol)
        if not ok:
            fails.append(f"{name}: {got!r} != {want!r}")

    eq("norm.3-4-5", norm([3.0, 4.0]), 5.0)
    eq("angle.x-pos", angle_deg([1.0, 0.0]), 0.0)
    eq("angle.y-neg", angle_deg([0.0, -1.0]), -90.0)
    eq("angle.third-quadrant", angle_deg([-1.0, -1.0]), -135.0)
    eq("cos.parallel", cosine([2.0, 0.0], [5.0, 0.0]), 1.0)
    eq("cos.antiparallel", cosine([2.0, 0.0], [-5.0, 0.0]), -1.0)
    eq("cos.orthogonal", cosine([1.0, 0.0], [0.0, 3.0]), 0.0)
    eq("cos.scale-invariant", cosine([0.2, -0.4], [2.0, -4.0]), 1.0)
    eq("angle_between.45", angle_between([1.0, 0.0], [1.0, 1.0]), 45.0)
    eq("cosine.zero-vector", cosine([0.0, 0.0], [1.0, 1.0]), None)
    eq("angle_between.zero-vector", angle_between([0.0, 0.0], [1.0, 1.0]), None)

    # bootstrap_ratio: deterministic for a fixed seed, brackets the estimate.
    h = [0.0040, 0.0041, 0.0042, 0.0039, 0.0043]
    s = [0.0028, 0.0029, 0.0030, 0.0027, 0.0031]
    lo1, hi1 = bootstrap_ratio(h, s)
    lo2, hi2 = bootstrap_ratio(h, s)
    eq("bootstrap.deterministic_lo", lo1, lo2)
    eq("bootstrap.deterministic_hi", hi1, hi2)
    point = (_median(s) - _median(h)) / _median(h)
    if not lo1 < point < hi1:
        fails.append("bootstrap interval does not bracket the estimate")
    if not lo1 < hi1:
        fails.append("bootstrap interval is empty")

    # The run and visit splitters on a hand-built sequence.
    seq = (["healthy"] * 2 + ["DAM 3"] * 3 + ["healthy"] * 1
           + ["DAM 6"] * 2 + ["DAM 3"] * 2)
    rows = [{"date": f"2021-01-{i + 1:02d}", "burst_id": f"b{i:02d}",
             "orbit": "ASCENDING", "damage_label": lab, "gamma2": 0.5, "P": 2.0}
            for i, lab in enumerate(seq)]
    _ordered, runs = runs_of(rows)
    if [r["state"] for r in runs] != ["healthy", "DAM 3", "healthy", "DAM 6",
                                      "DAM 3"]:
        fails.append(f"runs_of sequence wrong: {[r['state'] for r in runs]}")
    if [r["n"] for r in runs] != [2, 3, 1, 2, 2]:
        fails.append(f"runs_of lengths wrong: {[r['n'] for r in runs]}")
    visits = visits_of(runs)
    if [len(visits[st]) for st in STATES] != [2, 0, 1]:
        fails.append(f"visits_of counts wrong: "
                     f"{[len(visits[st]) for st in STATES]}")

    for f in fails:
        print(f"  FAIL {f}")
    print(f"selftest: {'ok' if not fails else str(len(fails)) + ' failures'}")
    return 0 if not fails else 1


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="LUMO physical dynamics against OCV behaviour.")
    ap.add_argument("--chip", default="support", choices=sorted(CHIPS),
                    help="channel table to draw (default: support)")
    ap.add_argument("--csv", default=None,
                    help="explicit channel table (overrides --chip)")
    ap.add_argument("--figdir", default=core.FIGDIR,
                    help="output directory of the JPEG figure")
    ap.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD, help="path of the Markdown report")
    ap.add_argument("--no-figures", action="store_true",
                    help="only recompute and cross-check (no matplotlib output)")
    ap.add_argument("--selftest", action="store_true",
                    help="run the geometry self-test and exit")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest()

    t0 = time.time()
    print("LUMO physical dynamics against OCV behaviour ...")
    res = compute(args.chip, args.csv)
    print(f"  compute: {time.time() - t0:.1f}s "
          f"({res['panel_b']['n_rows']} overpasses on the {res['chip']} chip, "
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
        json.dump(res, fh, indent=1, sort_keys=True)
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
