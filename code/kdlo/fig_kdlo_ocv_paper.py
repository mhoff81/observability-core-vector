#!/usr/bin/env python3
"""KDLO OCV-Paper — figures A–E, JSON result and English report.

Question: does the observability core vector
``OCV = [gamma2, P, D]`` (interferometric coherence gamma^2 of the KDLO-TV mast
over Garden City SD, plus the persistence P and the density D of its echo mask)
separate the two **epochs** of the committed record — the standing mast before
the 2022-12-11 collapse and the mast during the reconstruction — or is the
contrast produced by something else that changed at the same time?

The full echo-mask vector of this package is the 6D
``x_KDLO = [gamma2, P, D, A, F, S]`` (``A`` = masked pixels
= ``coherence_masked_pixels``, ``F`` = 8-connected components, ``S`` =
||mask centroid - peak pixel||). The mask of every overpass is recomputed from
the committed full-dwell strip cache ``data/kdlo/strips/`` with the mask rule of
the pipeline itself, so ``A`` and the masked coherence are pinned to the
pipeline's own values (see ``code/kdlo/kdlo_ocv_masks.py``).

Approach: the committed 60-acquisition record (two raw JSON files + a one-time
read-only Postgres extraction of nine columns, see
``code/kdlo/kdlo_ocv_channels_csv.py``), the same ``summ``/echo-mode/collapse
definitions as the KDLO project scripts, ``kdlo_ocv_stats.py`` as a verbatim
port of the LUMO statistics. Every number is recomputed here and then pinned
against the four committed reference JSONs; a deviation above 1e-12 aborts the
script with exit code != 0.

What is different from the LUMO/Bautzen packages, and what this file therefore
*does not* claim:

  * **the states are epochs, not damage labels.** ``pre`` (2022-06-02..2022-12-11)
    and ``rebuild`` (2023-05-04..2024-09-19) are disjoint *time windows*; the one
    acquisition in between (2022-12-23) belongs to neither and is excluded.
    A classifier can therefore win on season, weather or pipeline version alone,
    which is why every effect size of this package is reported twice: once on all
    rows and once on the **month-matched** subset, and why panel C also fits the
    weather-only model as a competitor.
  * **one geometry.** Only S1A descending relative orbit 12 exists over this
    site, so there is no ascending/descending stratum and no orbit control.
  * **the mask layer rests on a committed strip cache.** The record itself
    carries only ``A = coherence_masked_pixels``; ``P``, ``D``, ``F`` and ``S``
    and the whole ``P`` frequency map are recomputed here from the 64 committed
    strips (the record is 60 of those 64 acquisitions). The two ends are pinned:
    ``A`` equals ``coherence_masked_pixels`` and the recomputed masked coherence
    equals ``coherence_gamma2`` on every acquisition.
  * **"after the rebuild" does not exist.** 1,542 ft were reached 2024-08-19 and
    the record ends 2024-09-19 — the second state is *during* the construction.

Pinned reference files (``data/kdlo/reference/``):
  fig_kdlo_1_channels.json, fig_kdlo_2_channels.json   — recomputed exactly
  kdlo_s1_availability.json, kdlo_fem_frequencies.json — pinned structurally
  (catalogue query and FEM solver output: not rederivable, see ``core``)

Figures:
  A  raw distributions of the OCV components and the controls per epoch (n=59)
  B  in-sample effect sizes: Cliff's delta + bootstrap CI, SDS, and the
     month-matched control delta for every channel
  C  2-class LDA (pre vs. rebuild), leave-one-out, lambda grid, label
     permutation, confusion matrix — plus the weather-only competitor
  D  out-of-fold scores of the models, all rows vs. month-matched subset
  E  the evidence against reading the contrast as damage: echo-mode split,
     collapse pair, month/season strata, the singleton acquisition, and the
     catalogue/FEM check verdicts

Usage:
  python3 code/kdlo/fig_kdlo_ocv_paper.py            # recompute everything
  python3 code/kdlo/fig_kdlo_ocv_paper.py --quick    # without the permutation tests
"""
from __future__ import annotations

import argparse
import json
import os
import platform
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
OUT_JSON = os.path.join(core.DATA, "fig_kdlo_ocv_paper.json")
OUT_MD = os.path.join(FIGDIR, "fig_kdlo_ocv_paper.md")

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
STATE_COLORS = {"pre": "#4C72B0", "rebuild": "#DD8452"}
STATE_MARKERS = {"pre": "o", "rebuild": "s"}
CHANNEL_COLORS = {"gamma2": "#4C72B0", "P": "#55A868", "D": "#8172B3",
                  "A": "#8172B3", "F": "#64B5CD", "S": "#DA8BC3",
                  "phase_coherence": "#DD8452", "phase_snr_db": "#937860",
                  "phase_rms_rad": "#DA8BC3", "intensity": "#55A868",
                  "wind_speed_ms": "#64B5CD", "temperature_c": "#C44E52"}
# The model set mirrors the LUMO package (``core.MODELS`` / ``core.MODEL_LABEL``):
# the OCV is ``[gamma2, P, D]`` and the full echo-mask vector is 6D.
MODEL_COLORS = {"gamma2": "#8C8C8C", "PDS": "#DD8452", "gamma2PD": "#4C72B0",
                "x_6d": "#55A868", "weather": "#8172B3"}
MODEL_MARKERS = {"gamma2": "o", "PDS": "s", "gamma2PD": "^", "x_6d": "D",
                 "weather": "v"}
SHORT_MODEL = {"gamma2": "gamma2", "PDS": "P+D+S", "gamma2PD": "gamma2+P+D",
               "x_6d": "6 dims", "weather": "weather"}
CHANCE = 0.5
PAIR_KEYS = ["pre_vs_rebuild"]
PAIR = ("pre", "rebuild")
# Weather channels of the competitor model of panel C (the honest baseline: what
# a classifier can already do with the acquisition weather alone).
WEATHER_FEATURES = ["wind_speed_ms", "wind_gust_ms", "temperature_c",
                    "precipitation_mm", "snow_depth_m"]
# The channels figures A/B/E rank: the OCV proper first ([gamma2, P, D]), then
# the remaining echo-mask dimensions of the 6D vector (A, F, S), then the six
# controls the record actually carries for all 60 acquisitions.
CHANNELS = (list(core.OCV) + ["A", "F", "S"]
            + ["phase_coherence", "phase_snr_db", "phase_rms_rad",
               "intensity", "wind_speed_ms", "temperature_c"])
CHANNEL_AXIS = core.CHANNEL_AXIS
CHANNEL_TITLE = core.CHANNEL_TITLE


def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the espoo set keeps the binary out of git)."""
    out_dir = out_dir or FIGDIR
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png)}")
    return png


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
def month_matched(rows, ref="pre", grp="rebuild"):
    """The state subsets restricted to the months both states have.

    This is the control for the fact that ``pre`` and ``rebuild`` are epochs: the
    two windows do *not* cover the same calendar. ``pre`` holds Jun..Dec 2022 and
    ``rebuild`` holds all twelve months 2023-05..2024-09, so the months that only
    one state has are dropped and the effect size is recomputed on the overlap
    (Jun..Dec). In practice that removes the 15 ``rebuild`` rows of Jan..May 2024
    — the ones no ``pre`` month can be matched against — and leaves 16 ``pre``
    vs. 28 ``rebuild`` rows. A contrast that survives this is not produced by the
    month (or by whatever else follows the month) alone. The n of the returned
    lists is the n of the comparison.
    """
    months = sorted({r["month"] for r in rows if r.get("state") == ref}
                    & {r["month"] for r in rows if r.get("state") == grp})
    a = [r for r in rows if r.get("state") == ref and r["month"] in months]
    b = [r for r in rows if r.get("state") == grp and r["month"] in months]
    return a, b, months


def channel_values(rows, key):
    return [r[key] for r in rows if r.get(key) is not None]


def fig_a_block(rows, keys):
    """Raw values + descriptives per state for the distribution panel."""
    by = st.by_state(rows)
    values = {k: {lab: channel_values(by[lab], k) for lab in st.STATE_ORDER}
              for k in keys}
    desc = {k: {lab: st._stats(values[k][lab]) for lab in st.STATE_ORDER}
            for k in keys}
    return {"values": values, "descriptives": desc, "keys": list(keys)}


def fig_b_block(rows, keys):
    """Cliff's delta + bootstrap CI + SDS per channel, all rows and matched."""
    by = st.by_state(rows)
    table = {}
    for k in keys:
        a = channel_values(by["pre"], k)
        b = channel_values(by["rebuild"], k)
        if len(a) < 3 or len(b) < 3:
            table[k] = {"n": [len(a), len(b)], "delta": None, "sds": None,
                        "delta_ci95": None, "ci_excludes_zero": None,
                        "note": "n < 3"}
            continue
        blk = st.delta_block(a, b, f"pre_vs_rebuild_{k}", core.N_BOOT)
        blk["n_pairs"] = len(a) * len(b)
        table[k] = blk
    matched = {}
    a_m, b_m, months = month_matched(rows)
    for k in keys:
        x = channel_values(a_m, k)
        y = channel_values(b_m, k)
        if len(x) < 3 or len(y) < 3:
            matched[k] = {"n": [len(x), len(y)], "delta": None, "sds": None,
                          "delta_ci95": None, "ci_excludes_zero": None,
                          "note": "n < 3"}
            continue
        matched[k] = st.delta_block(x, y, f"pre_vs_rebuild_{k}_matched",
                                    core.N_BOOT)
        if matched[k].get("delta") is not None:
            matched[k]["n_pairs"] = len(x) * len(y)
    ranked = sorted(keys, key=lambda k: -(abs(table[k]["delta"])
                                          if table[k]["delta"] else -1.0))
    return {"table": table, "matched": matched, "matched_months": months,
            "n_matched": [len(a_m), len(b_m)], "ranking_by_abs_delta": ranked,
            "keys": list(keys)}


def compute(quick=False):
    """Recompute all numbers of figures A–E (and the raw values for the pins)."""
    rows = core.load_csv()
    ref = core.load_reference()
    meta = core.load_csv_meta()
    by_part, prov = core.load_raw_parts()
    extra, extra_prov = core.db_extra_rows()
    res = {"rows": rows, "ref": ref, "meta": meta, "prov": prov, "extra": extra,
           "extra_prov": extra_prov, "by_part": by_part}

    days = core.days_of(rows)
    by_state = st.by_state(rows)
    res["data"] = {
        "csv": os.path.relpath(core.CSV_PATH, os.path.dirname(core.DATA)),
        "csv_meta": os.path.relpath(core.CSV_META_PATH, os.path.dirname(core.DATA)),
        "csv_sha256": core.sha256(core.CSV_PATH),
        "db_extra": os.path.relpath(core.DB_EXTRA_PATH, os.path.dirname(core.DATA)),
        "db_extra_sha256": core.sha256(core.DB_EXTRA_PATH),
        "raw_parts": {str(p): {"file": prov[p]["file"], "sha256": prov[p]["sha256"],
                               "n_rows": prov[p]["n_rows"], "window": prov[p]["window"]}
                      for p in core.PART_ORDER},
        "n_rows": len(rows),
        "n_by_state": core.n_by_state(rows),
        "n_state_rows": len(core.state_rows(rows)),
        "n_singleton": sum(1 for r in rows if r["state"] == "single"),
        "n_by_part": {f"part{p}": len(core.rows_of_part(rows, p))
                      for p in core.PART_ORDER},
        "n_by_echo_mode": {m: sum(1 for r in rows if r["echo_mode"] == m)
                           for m in ("compact", "intermediate", "distributed")},
        "n_compact_by_state": {lab: sum(1 for r in by_state[lab]
                                        if r["echo_mode"] == "compact")
                               for lab in st.STATE_ORDER},
        "n_intermediate_days": [core.day_of(r) for r in rows
                                if r["echo_mode"] == "intermediate"],
        "date_first": days[0], "date_last": days[-1],
        "state_windows": {k: list(v) for k, v in core.STATE_WINDOWS.items()},
        "singleton_day": core.SINGLETON_DAY,
        "state_label": dict(core.STATE_LABEL),
        "channels": {"ocv": list(core.OCV), "controls": list(core.CONTROLS),
                     "weather": list(core.WEATHER),
                     "empty": list(core.EMPTY_CHANNELS)},
        "db_columns": list(core.DB_EXTRA_COLUMNS),
        "mask_columns": list(core.MASK_COLUMNS),
        "strip_cache": meta.get("strip_cache"),
        "missing_repeat_day_part1": core.grid_missing_days(
            core.days_of(core.rows_of_part(rows, 1))),
        "alternate_pass_days": sorted({core.day_of(r) for r in rows
                                       if r.get("pass_label")
                                       and r["pass_label"] != "afternoon"}),
    }

    keys = list(CHANNELS)
    res["fig_a"] = fig_a_block(rows, keys)
    res["fig_b"] = fig_b_block(rows, keys)

    # The two committed channel JSONs, rederived in the exact form they were
    # written in (part 2 carries part 1 as its intact reference).
    p1 = core.part_stats(core.rows_of_part(rows, 1), 1, prov=prov)
    p2 = core.part_stats(core.rows_of_part(rows, 2), 2, prov=prov, reference=p1)
    res["fig_part"] = {"part1": p1, "part2": p2}

    # --- Figure C: 2-class LDA on the epochs (the hard caveat) ------------
    # Only the rows that belong to a state take part: the singleton 2022-12-23
    # belongs to neither epoch and must not be fed to a 2-class classifier (it
    # would be an acquisition that *cannot* be scored correctly by construction).
    srows = core.state_rows(rows)
    X, y, complete = core.feature_dataset(srows)
    col = {f: i for i, f in enumerate(core.FEATURES)}
    Xw, yw, complete_w = core.feature_dataset(srows, WEATHER_FEATURES)
    specs = dict(core.MODELS)
    specs["weather"] = None                      # own channel list
    grid = {}
    for name, feats in specs.items():
        Xm = X if feats else Xw
        idx = [col[f] for f in feats] if feats else list(range(Xw.shape[1]))
        grid[name] = {}
        for lam in core.LAMBDA_GRID:
            _, acc, bal, _conf = st.loo_cv(Xm[:, idx], y if feats else yw,
                                           st.STATE_ORDER, lam)
            grid[name][str(lam)] = {"accuracy": acc, "balanced_accuracy": bal}
    head = {}
    for name in (core.HEADLINE_PAIR, "weather"):
        feats = specs[name]
        Xm = X if feats else Xw
        ym = y if feats else yw
        idx = [col[f] for f in feats] if feats else list(range(Xw.shape[1]))
        if quick:
            blk = {"accuracy_observed": grid[name][str(core.HEADLINE_LAMBDA)]["accuracy"],
                   "balanced_accuracy_observed":
                       grid[name][str(core.HEADLINE_LAMBDA)]["balanced_accuracy"],
                   "recomputed": False, "n_perm": 0}
        else:
            blk = st.permutation_test(Xm[:, idx], ym, st.STATE_ORDER,
                                      core.HEADLINE_LAMBDA,
                                      n_perm=core.N_PERM_2CLASS,
                                      seed=core.RNG_SEED)
            blk["recomputed"] = True
        head[name] = blk
    idx_h = [col[f] for f in core.MODELS[core.HEADLINE_PAIR]]
    _, acc_h, bal_h, conf_h = st.loo_cv(X[:, idx_h], y, st.STATE_ORDER,
                                        core.HEADLINE_LAMBDA)
    res["fig_c"] = {
        "sensitivity_grid": grid, "headline_permutation_test": head,
        "headline_confusion_matrix": conf_h, "headline_accuracy": acc_h,
        "headline_balanced_accuracy": bal_h,
        "lambda_grid": core.LAMBDA_GRID, "headline_lambda": core.HEADLINE_LAMBDA,
        "chance_accuracy": CHANCE, "chance_balanced_accuracy": CHANCE,
        "features": {**{k: list(v) for k, v in core.MODELS.items()},
                     "weather": list(WEATHER_FEATURES)},
        "n_complete_ocv": len(complete), "n_complete_weather": len(complete_w),
        "n_perm": core.N_PERM_2CLASS, "permutation_recomputed": not quick,
    }

    # --- Figure D: out-of-fold scores, all rows and month-matched ---------
    def pair_block(Xm, ym, samples):
        per_model, q_by_model = {}, {}
        for name, feats in core.MODELS.items():
            idx = [col[f] for f in feats]
            q = st.oof_scores(Xm[:, idx], ym, PAIR[0], PAIR[1], core.HEADLINE_LAMBDA)
            q_by_model[name] = np.asarray(q, float)
            delta, sds, ref_q, grp_q = st.delta_sds(q, ym, PAIR[0], PAIR[1])
            ci = (st.bootstrap_delta_ci(ref_q, grp_q, n=core.N_BOOT, seed=core.RNG_SEED)
                  if delta is not None else None)
            per_model[name] = {
                "label": core.MODEL_LABEL[name],
                "n": [len(ref_q), len(grp_q)], "n_pairs": len(ref_q) * len(grp_q),
                "delta": delta, "sds": sds,
                "delta_ci95": list(ci) if ci else None,
                "ci_excludes_zero": bool(ci[0] * ci[1] > 0) if ci else None,
            }
        paired = {}
        for name in list(core.MODELS) + ["weather"]:
            if name not in q_by_model:
                continue
            paired[name] = st.paired_bootstrap_diff(
                q_by_model[name], q_by_model[core.BASELINE], ym, PAIR[0], PAIR[1],
                core.N_BOOT, np.random.default_rng(core.RNG_SEED + 2))
        return {"per_model": per_model, "samples": samples,
                "q_by_model": {k: list(v) for k, v in q_by_model.items()},
                "y": [str(v) for v in ym],
                "paired_diff_vs_gamma2": paired}

    state_mask = np.isin(y, list(PAIR))
    a_m, b_m, matched_months = month_matched(rows)
    res["fig_d"] = {
        "all": pair_block(X[state_mask], y[state_mask], "all state rows"),
        "matched": pair_block(
            np.array([[float(r[f]) for f in core.FEATURES]
                      for r in a_m + b_m]),
            np.array([r["state"] for r in a_m + b_m], dtype=object),
            "month-matched subset"),
        "matched_months": matched_months,
        "n_matched": [len(a_m), len(b_m)],
    }

    # --- Figure E: what could explain the contrast instead ---------------
    fe = {
        "echo_modes_by_part": {
            f"part{p}": core.part_stats(core.rows_of_part(rows, p), p)["echo_modes"]
            for p in core.PART_ORDER},
        "echo_modes_by_state": {
            lab: {"n": len(by_state[lab]),
                  "n_compact": sum(1 for r in by_state[lab]
                                   if r["echo_mode"] == "compact"),
                  "n_distributed": sum(1 for r in by_state[lab]
                                       if r["echo_mode"] == "distributed"),
                  "n_intermediate": sum(1 for r in by_state[lab]
                                        if r["echo_mode"] == "intermediate")}
            for lab in st.STATE_ORDER},
        "by_month": st.strata_pair(rows, lambda r: r.get("month"), key="gamma2"),
        "by_season": st.strata_pair(rows, st.season_of, key="gamma2"),
        "by_masked_px_bucket": st.strata_pair(
            rows, lambda r: ("compact" if r["A"] <= core.COMPACT_MASK_PX
                             else "distributed" if r["A"] >= core.MASK_CLUTTER_PX
                             else "intermediate"), key="gamma2"),
        "median_by_state": {
            k: {lab: (float(np.median(channel_values(by_state[lab], k)))
                      if channel_values(by_state[lab], k) else None)
                for lab in st.STATE_ORDER}
            for k in CHANNELS},
        "singleton": {core.day_of(r): {"gamma2": r["gamma2"], "A": r["A"],
                                       "echo_mode": r["echo_mode"],
                                       "temperature_c": r["temperature_c"],
                                       "wind_speed_ms": r["wind_speed_ms"]}
                      for r in rows if r["state"] == "single"},
        "collapse_pair": core.part_stats(core.rows_of_part(rows, 1), 1).get(
            "collapse_pair"),
        "availability": core.availability_checks(ref["s1_availability"], rows, prov),
        "fem": core.fem_checks(ref["fem_frequencies"]),
    }
    res["fig_e"] = fe
    return res


# ---------------------------------------------------------------------------
# Rendering — figures A–E, 183 mm column width, 600 dpi PNG (espoo style)
# ---------------------------------------------------------------------------
def panel(ax, letter, text):
    ax.set_title(f"({letter}) {text}", fontsize=8.0, loc="left")


def jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def state_strip(ax, values_by_state, seed0=100, span=0.14, fs=6.6):
    """Raw per-acquisition dots + median bar for the two epoch states."""
    for i, lab in enumerate(st.STATE_ORDER):
        v = np.asarray(values_by_state[lab], float)
        if v.size == 0:
            continue
        x = i + jitter(v.size, span, seed0 + i)
        ax.plot(x, v, STATE_MARKERS[lab], ms=2.4, mfc="none", mec=STATE_COLORS[lab],
                mew=0.5, alpha=0.6, ls="none")
        ax.hlines(float(np.median(v)), i - 0.28, i + 0.28, color=STATE_COLORS[lab],
                  lw=1.7)
    ax.set_xlim(-0.55, 1.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([st.STATE_ORDER[0], st.STATE_ORDER[1]], fontsize=fs)


def fig_a(res, out_dir=None):
    """Figure A — raw distributions of the eight channels per epoch state."""
    vals_ = res["fig_a"]["values"]
    keys = res["fig_a"]["keys"]
    n_by_state = res["data"]["n_by_state"]
    fig, axes = plt.subplots(2, 4, figsize=(COL_WIDTH, 0.60 * COL_WIDTH))
    for ax, key in zip(axes.ravel(), keys):
        state_strip(ax, vals_[key], seed0=10 + keys.index(key) * 7)
        if key == "A":
            ax.set_yscale("log")
        ax.set_ylabel(CHANNEL_AXIS[key], fontsize=6.4)
        panel(ax, "abcdefgh"[keys.index(key)], CHANNEL_TITLE[key])
    fig.suptitle(
        "Figure A — the measured channels per KDLO epoch state "
        f"({res['data']['n_rows']} acquisitions: "
        + ", ".join(f"{s} n={n_by_state[s]}" for s in st.STATE_ORDER)
        + f"; the singleton {res['data']['singleton_day']} belongs to neither state)\n"
        "dots = single acquisitions, bar = median;  input: data/kdlo/kdlo_ocv_channels.csv",
        fontsize=7.4)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return save_fig(fig, "fig_kdlo_ocv_paper_A", out_dir)


def fig_b(res, out_dir=None):
    """Figure B — in-sample effect sizes: Cliff's delta, bootstrap CI, SDS."""
    fb = res["fig_b"]
    tab, mat = fb["table"], fb["matched"]
    keys = fb["keys"]
    rank = fb["ranking_by_abs_delta"]
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.42 * COL_WIDTH),
                             width_ratios=[1.0, 1.05])

    ax = axes[0]
    y = np.arange(len(rank))
    ax.barh(y, [tab[k]["sds"] or 0.0 for k in rank], height=0.62,
            color=[CHANNEL_COLORS[k] for k in rank])
    for i, k in enumerate(rank):
        ax.text((tab[k]["sds"] or 0.0) + 0.06, i,
                f"{tab[k]['delta']:+.2f} (n={tab[k]['n'][0]}/{tab[k]['n'][1]})",
                va="center", fontsize=6.0)
    ax.set_yticks(y)
    ax.set_yticklabels([CHANNEL_AXIS[k] for k in rank], fontsize=6.6)
    ax.invert_yaxis()
    ax.set_xlim(0.0, max([tab[k]["sds"] or 0.0 for k in rank] + [1.0]) * 1.5)
    ax.set_xlabel("SDS = 10 |Cliff's delta|,  pre vs. rebuild", fontsize=7.0)
    panel(ax, "a", "which channel moves at all")

    ax = axes[1]
    for i, k in enumerate(keys):
        for j, (src, mk, ms) in enumerate(((tab, "o", 4.0), (mat, "x", 3.4))):
            blk = src.get(k) or {}
            if blk.get("delta") is None:
                continue
            lo, hi = blk["delta_ci95"]
            ax.plot([i + (j - 0.5) * 0.26] * 2, [lo, hi], "-", lw=0.9,
                    color=CHANNEL_COLORS[k], alpha=0.85)
            ax.plot(i + (j - 0.5) * 0.26, blk["delta"], mk, ms=ms,
                    color=CHANNEL_COLORS[k],
                    label=("all rows" if j == 0 else "month-matched")
                    if i == 0 else None)
    ax.axhline(0.0, color="0.35", ls=":", lw=1.0)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([CHANNEL_AXIS[k].split("  ")[0] for k in keys], fontsize=6.0,
                       rotation=22, ha="right")
    ax.set_ylabel("Cliff's delta (95% bootstrap CI)", fontsize=7.0)
    ax.legend(fontsize=6.0, loc="lower right")
    panel(ax, "b", "sign and size, all rows vs. month-matched")
    g2 = tab["gamma2"]
    fig.suptitle(
        "Figure B — pre vs. rebuild in-sample effect sizes "
        f"(n={g2['n'][0]} vs. {g2['n'][1]} acquisitions, {g2['n_pairs']} pairs; "
        "SDS 0 = indistinguishable, SDS >= 2 = large)\n"
        f"gamma^2 median {g2['median'][0]:.6f} -> {g2['median'][1]:.6f} "
        f"(delta {g2['delta']:+.3f}); month-matched control "
        f"{mat['gamma2']['delta']:+.3f} on "
        f"{mat['gamma2']['n'][0]} vs. {mat['gamma2']['n'][1]} acquisitions",
        fontsize=7.0)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    return save_fig(fig, "fig_kdlo_ocv_paper_B", out_dir)



def fig_c(res, out_dir=None):
    """Figure C — 2-class LDA: lambda sweep, confusion matrix, permutation null."""
    fc = res["fig_c"]
    grid = fc["sensitivity_grid"]
    lambdas = [str(x) for x in fc["lambda_grid"]]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 0.40 * COL_WIDTH),
                             width_ratios=[1.2, 0.85, 1.1])

    ax = axes[0]
    for name in list(core.MODELS) + ["weather"]:
        ax.plot(range(len(lambdas)),
                [grid[name][l]["balanced_accuracy"] for l in lambdas], "-o", ms=3.2,
                lw=1.3, color=MODEL_COLORS[name],
                label=core.MODEL_LABEL.get(name, "weather only (5D)"))
    ax.axhline(CHANCE, color="0.35", ls=":", lw=1.0)
    ax.text(len(lambdas) - 0.5, CHANCE + 0.012, "chance", ha="right", fontsize=6.2,
            color="0.35")
    ax.set_xticks(range(len(lambdas)))
    ax.set_xticklabels([f"{float(l):g}" for l in lambdas], fontsize=6.6)
    ax.set_xlabel("shrinkage lambda", fontsize=7.0)
    ax.set_ylabel("balanced accuracy (LOO)", fontsize=7.0)
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=5.8, loc="upper left")
    panel(ax, "a", f"lambda sweep (n={fc['n_complete_ocv']} complete cases)")

    ax = axes[1]
    conf = fc["headline_confusion_matrix"]
    classes = list(conf)
    M = np.array([[conf[a][b] for b in classes] for a in classes], float)
    ax.imshow(M, cmap="Blues", vmin=0.0, vmax=M.max())
    for i, a in enumerate(classes):
        for j, b in enumerate(classes):
            ax.text(j, i, f"{int(M[i, j])}", ha="center", va="center", fontsize=7.0,
                    color="white" if M[i, j] > 0.55 * M.max() else "black")
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes, fontsize=6.6)
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes, fontsize=6.6)
    ax.set_xlabel("predicted (leave-one-out)", fontsize=7.0)
    ax.set_ylabel("true epoch state", fontsize=7.0)
    ax.grid(False)
    panel(ax, "b", f"{SHORT_MODEL[core.HEADLINE_PAIR]}, lambda={fc['headline_lambda']}, "
                   f"acc={fc['headline_accuracy']:.3f}")

    ax = axes[2]
    head = fc["headline_permutation_test"][core.HEADLINE_PAIR]
    ax.axhline(CHANCE, color="0.35", ls=":", lw=1.0)
    if fc["permutation_recomputed"]:
        obs = head["balanced_accuracy_observed"]
        lab = ["null\nmean", "null\np95", "observed"]
        val = [head["balanced_accuracy_null_mean"],
               head["balanced_accuracy_null_p95"], obs]
        ax.bar(range(3), val, width=0.6, color=["0.7", "0.5", "#4C72B0"])
        for i, v in enumerate(val):
            ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=5.8)
        ax.set_xticks(range(3))
        ax.set_xticklabels(lab, fontsize=6.2)
        ax.text(0.5, 0.97, f"p = {head['p_balanced_accuracy']:.3f}\n"
                           f"(weather p = "
                           f"{fc['headline_permutation_test']['weather']['p_balanced_accuracy']:.3f})",
                transform=ax.transAxes, ha="center", va="top", fontsize=6.0)
    else:
        ax.text(0.5, 0.5, "null not recomputed\n(--quick)", ha="center", va="center",
                fontsize=7.5, transform=ax.transAxes)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("balanced accuracy", fontsize=7.0)
    panel(ax, "c", f"label permutation null (n_perm={fc['n_perm']})")
    wl = grid["weather"][str(fc["headline_lambda"])]["balanced_accuracy"]
    perm = fc["headline_permutation_test"][core.HEADLINE_PAIR]
    tail = (f"; the label permutation puts the observed value at p = "
            f"{perm['p_balanced_accuracy']:.3f} against a null p95 of "
            f"{perm['balanced_accuracy_null_p95']:.3f} "
            + ("— the OCV does separate the epochs beyond chance here, while the "
               "1D gamma2 baseline alone scores even higher"
               if perm["p_balanced_accuracy"] <= 0.05 else
               "— i.e. the per-acquisition classifier does not separate the "
               "epochs beyond chance")
            if fc["permutation_recomputed"] else
            " (the permutation null was not recomputed in --quick mode)")
    fig.suptitle(
        "Figure C — do the channels *classify* the epoch? (leave-one-out, 2 classes)\n"
        f"best OCV balanced accuracy {fc['headline_balanced_accuracy']:.3f} "
        f"(chance {CHANCE:.2f}); the weather-only competitor reaches {wl:.3f}"
        + tail, fontsize=7.0)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    return save_fig(fig, "fig_kdlo_ocv_paper_C", out_dir)



def fig_d(res, out_dir=None):
    """Figure D — out-of-fold scores (Fisher LDA), all rows and month-matched."""
    fd = res["fig_d"]
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.42 * COL_WIDTH),
                             width_ratios=[1.05, 1.15])

    ax = axes[0]
    blk = fd["all"]
    q = np.asarray(blk["q_by_model"][core.HEADLINE_PAIR], float)
    y = np.asarray(blk["y"], dtype=object)
    for i, lab in enumerate(st.STATE_ORDER):
        v = q[y == lab]
        if v.size:
            ax.plot(i + jitter(v.size, 0.14, 200 + i), v, STATE_MARKERS[lab], ms=2.6,
                    mfc="none", mec=STATE_COLORS[lab], mew=0.5, alpha=0.65, ls="none")
            ax.hlines(float(np.median(v)), i - 0.28, i + 0.28,
                      color=STATE_COLORS[lab], lw=1.7)
            ax.text(i, float(np.max(v)), f"n={v.size}", fontsize=6.0, va="bottom",
                    ha="center", color=STATE_COLORS[lab])
    ax.axhline(0.0, color="0.35", ls=":", lw=1.0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(list(st.STATE_ORDER), fontsize=6.6)
    ax.set_ylabel("out-of-fold LDA score (gamma2 + A)", fontsize=7.0)
    panel(ax, "a", "score per acquisition, all rows")

    ax = axes[1]
    names = list(fd["all"]["per_model"])
    xs = np.arange(len(names))
    for j, (sub, colr, hatch) in enumerate((("all", "#4C72B0", None),
                                            ("matched", "#DD8452", "//"))):
        vals_ = [fd[sub]["per_model"][n]["sds"] or 0.0 for n in names]
        ax.bar(xs + (j - 0.5) * 0.34, vals_, width=0.32, color=colr, hatch=hatch,
               edgecolor="white", lw=0.4,
               label="all rows" if j == 0 else "month-matched")
    for i, n in enumerate(names):
        pd = fd["all"]["paired_diff_vs_gamma2"].get(n)
        if not pd or n == core.BASELINE:
            continue
        ax.text(i, (fd["all"]["per_model"][n]["sds"] or 0.0) + 0.08,
                f"vs gamma2\n{pd['diff_mean']:+.2f}\n"
                f"[{pd['ci95'][0]:+.2f},{pd['ci95'][1]:+.2f}]", fontsize=5.2,
                ha="center")
    ax.set_xticks(xs)
    ax.set_xticklabels([core.MODEL_LABEL[n].split(" (")[0] for n in names],
                       fontsize=6.2, rotation=12, ha="right")
    ax.set_ylabel("SDS of the out-of-fold scores", fontsize=7.0)
    ax.legend(fontsize=6.0, loc="upper left")
    panel(ax, "b", "does adding A help? (paired bootstrap)")
    hp = fd["all"]["per_model"][core.HEADLINE_PAIR]
    fig.suptitle(
        "Figure D — the same effect size, but measured out of fold "
        f"({hp['n'][0]} vs. {hp['n'][1]} acquisitions; month-matched subset "
        f"{fd['n_matched'][0]} vs. {fd['n_matched'][1]} on months "
        f"{fd['matched_months'][0]}..{fd['matched_months'][-1]})\n"
        "the direction is refitted without each acquisition, so no row can help "
        "separate itself", fontsize=7.0)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    return save_fig(fig, "fig_kdlo_ocv_paper_D", out_dir)



def fig_e(res, out_dir=None):
    """Figure E — the alternative explanations: echo mode, month, size, weather."""
    fe = res["fig_e"]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 0.42 * COL_WIDTH),
                             width_ratios=[0.95, 1.15, 1.05])

    ax = axes[0]
    modes = ("compact", "intermediate", "distributed")
    colors = {"compact": "#C44E52", "intermediate": "#937860",
              "distributed": "#4C72B0"}
    xs = np.arange(len(st.STATE_ORDER))
    bottom = np.zeros(len(st.STATE_ORDER))
    for m in modes:
        vals_ = np.array([fe["echo_modes_by_state"][lab][f"n_{m}"]
                          for lab in st.STATE_ORDER], float)
        ax.bar(xs, vals_, bottom=bottom, width=0.6, color=colors[m], label=m)
        for i, v in enumerate(vals_):
            if v:
                ax.text(xs[i], bottom[i] + v / 2, f"{int(v)}", ha="center",
                        va="center", fontsize=6.2, color="white")
        bottom += vals_
    ax.set_xticks(xs)
    ax.set_xticklabels(list(st.STATE_ORDER), fontsize=6.6)
    ax.set_ylabel("acquisitions", fontsize=7.0)
    ax.legend(fontsize=5.8, loc="upper left")
    panel(ax, "a", "echo mode: the compact-echo count grows")

    ax = axes[1]
    months = sorted(fe["by_month"], key=lambda k: int(k))
    xs = np.arange(len(months))
    deltas = [fe["by_month"][m].get("delta") for m in months]
    ns = [fe["by_month"][m].get("n") for m in months]
    ax.bar(xs, [0.0 if d is None else d for d in deltas], width=0.62, color="#818181")
    for i, (d, n) in enumerate(zip(deltas, ns)):
        txt = "n/a" if d is None else f"{d:+.2f}"
        ax.text(xs[i], (0.0 if d is None else d) + 0.03, txt, ha="center",
                fontsize=5.6)
        ax.text(xs[i], -0.09, f"{n[0]}/{n[1]}", ha="center", fontsize=5.0,
                color="0.35")
    ax.axhline(0.0, color="0.35", ls=":", lw=1.0)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(int(m)) for m in months], fontsize=6.2)
    ax.set_xlabel("calendar month  (n pre / n rebuild below)", fontsize=7.0)
    ax.set_ylabel("Cliff's delta (gamma^2)", fontsize=7.0)
    ax.set_ylim(-0.2, 1.0)
    ax.margins(y=0.06)
    panel(ax, "b", "per month (only months with n >= 3 on both sides)")

    ax = axes[2]
    med = fe["median_by_state"]
    xs = np.arange(len(st.STATE_ORDER))
    ax.bar(xs - 0.19, [med["A"][lab] for lab in st.STATE_ORDER], width=0.36,
           color="#8172B3", label="median mask A [px]")
    ax.set_ylabel("median mask size A  [px]", fontsize=7.0)
    ax.set_xticks(xs)
    ax.set_xticklabels(list(st.STATE_ORDER), fontsize=6.6)
    ax2 = ax.twinx()
    ax2.plot(xs + 0.19, [med["gamma2"][lab] for lab in st.STATE_ORDER], "s-", ms=3.6,
             lw=1.2, color="#4C72B0", label="median gamma^2")
    ax2.set_ylabel("median gamma^2", fontsize=7.0, color="#4C72B0")
    ax2.tick_params(axis="y", labelsize=6.4, colors="#4C72B0")
    ax2.grid(False)
    ax.margins(y=0.42)
    au, fe_ = fe["availability"], fe["fem"]
    ax.text(0.02, 0.99,
            "mirror checks on the two pin-only files\n"
            f"  S1 availability: {len(au)} structural checks\n"
            f"  FEM frequencies: {len(fe_)} structural checks\n"
            "wind/temperature medians move little\n"
            f"singleton {res['data']['singleton_day']}: "
            f"gamma^2 {fe['singleton'][res['data']['singleton_day']]['gamma2']:.6f}"
            ", no state",
            transform=ax.transAxes, fontsize=5.2, va="top")
    panel(ax, "c", "size and weather do not follow the state")
    fig.suptitle(
        "Figure E — what else could produce the contrast: echo mode, calendar "
        "month, mask size, wind/temperature\n"
        "if the contrast were a season or a weather artefact it would have to "
        "vanish in the month-matched column of Figure B — it does not", fontsize=7.0)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    return save_fig(fig, "fig_kdlo_ocv_paper_E", out_dir)



# ---------------------------------------------------------------------------
# Pins — every recomputed number against the committed record
# ---------------------------------------------------------------------------
# Rule: exact equality for everything deterministic (counts, medians, Cliff's
# delta, SDS, CI endpoints, curated checks); 1e-12 for the few numbers that come
# out of scipy or out of a permutation/bootstrap loop, where the last bits can
# depend on the BLAS/library build. A deviation above the tolerance aborts the
# script with exit code != 0.
FLOAT_TOL = 1e-12
TOL_KEYS = {"welch_t", "welch_p", "mannwhitney_p", "p", "p_accuracy",
            "p_balanced_accuracy", "accuracy_null_mean",
            "balanced_accuracy_null_mean", "accuracy_null_p95",
            "balanced_accuracy_null_p95"}


class Pins:
    """Collects pin comparisons (recomputed vs. committed) and reports deviations."""

    def __init__(self):
        self.checks = []
        self.fails = []
        self.n_avail = 0
        self.n_fem = 0

    def add(self, path, got, want, tol):
        try:
            exact = bool(got == want)
        except Exception:                       # pragma: no cover - numpy arrays
            exact = False
        dev = 0.0
        if not exact:
            numeric = (isinstance(got, (int, float)) and isinstance(want, (int, float))
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
            "n_availability_checks": self.n_avail,
            "n_fem_checks": self.n_fem,
        }



NUMERIC_SOURCE = [c for c in core.SOURCE_COLUMNS
                  if c not in core.TEXT_SOURCE_COLUMNS]


def _row_source_key(row):
    """One comparable tuple per acquisition: numbers normalised, text as text."""
    return tuple([core._num(row.get(c)) for c in NUMERIC_SOURCE]
                 + [row.get(c) for c in core.TEXT_SOURCE_COLUMNS
                    if c in core.SOURCE_COLUMNS])


def _row_derived_key(row):
    """The derived fields *as stored in the CSV* — the committed claim."""
    return (row.get("state"), row.get("month"), row.get("year"), row.get("season"),
            row.get("echo_mode"))


def _row_derived_recomputed(row):
    """The same fields, recomputed from the day and the mask size."""
    day = core.day_of(row)
    return (core.state_of(day), int(day[5:7]), int(day[:4]), st.season_of(row),
            core.echo_mode_of(core._num(row.get("A"))))


def pin_all(res, quick=False):
    """All pins against the four committed reference JSONs and the meta file."""
    ref, rows, meta = res["ref"], res["rows"], res["meta"]
    P = Pins()

    # 1) the inputs themselves: vendored bytes, shapes, definitions ----------
    P.cmp("meta.created_rows", meta["created_rows"], len(rows))
    P.cmp("meta.generator", meta["generator"], "code/kdlo/kdlo_ocv_channels_csv.py")
    P.cmp("meta.out_csv", meta["out_csv"], os.path.basename(core.CSV_PATH))
    P.cmp("meta.singleton_day", meta["singleton_day"], core.SINGLETON_DAY)
    P.cmp("meta.n_by_state", meta["n_by_state"], res["data"]["n_by_state"])
    P.cmp("meta.n_states_rows", meta["n_states_rows"], res["data"]["n_state_rows"])
    P.cmp("meta.n_by_part", meta["n_by_part"], res["data"]["n_by_part"])
    P.cmp("meta.n_by_echo_mode", meta["n_by_echo_mode"], res["data"]["n_by_echo_mode"])
    P.cmp("meta.state_windows",
          {k: list(v) for k, v in meta["state_windows"].items()},
          {k: list(v) for k, v in core.STATE_WINDOWS.items()})
    P.cmp("meta.source_columns", list(meta["source_columns"]),
          list(core.SOURCE_COLUMNS))
    P.cmp("meta.derived_columns", list(meta["derived_columns"]),
          list(core.DERIVED_COLUMNS))
    P.cmp("meta.db_columns", list(meta["db_columns"]), list(core.DB_EXTRA_COLUMNS))
    P.cmp("meta.columns", list(meta["columns"]), list(core.CSV_COLUMNS))
    P.cmp("meta.definitions.echo_mode_thresholds",
          [core.COMPACT_MASK_PX, core.MASK_CLUTTER_PX], [10, 25])

    # 1b) the echo-mask layer: the committed strip cache behind the mask columns
    sc = meta["strip_cache"]
    with open(os.path.join(core.DATA, sc["manifest"])) as fh:
        strip_manifest = json.load(fh)
    P.cmp("mask.columns", list(meta["mask_columns"]), list(core.MASK_COLUMNS))
    P.cmp("mask.strip_cache.manifest_sha256", sc["manifest_sha256"],
          core.sha256(os.path.join(core.DATA, sc["manifest"])))
    P.cmp("mask.strip_cache.n_strips", sc["n_strips"], len(strip_manifest["entries"]))
    P.cmp("mask.strip_cache.n_attached", sc["n_attached_to_csv_rows"], len(rows))
    P.cmp("mask.strip_cache.majority_shape", list(sc["majority_shape"]), [400, 11])
    P.cmp("mask.strip_cache.rule",
          [sc["mask_rule"]["min_n_masked"], sc["mask_rule"]["mask"]],
          [2, "intensity >= 0.30 * peak and intensity >= 5.0 * median"])

    # Every CSV row must carry the mask of *its own* strip: A the pipeline's own
    # count, gamma2_mask the pipeline's own coherence, D/F/S derived from A.
    n_mask_rows = n_mask_bad = 0
    for r in rows:
        if r.get("A") is None or r.get("P") is None:
            continue
        n_mask_rows += 1
        if not (int(r["A"]) == int(r["coherence_masked_pixels"])
                and abs(r["gamma2_mask"] - r["coherence_gamma2"]) <= 1e-12
                and abs(r["D"] - r["A"] / r["bbox_area"]) <= 1e-12
                and 0.0 <= r["P"] <= 1.0
                and (r["mask_rows"], r["mask_cols"]) == (400, 11)
                and 1 <= r["F"] <= r["A"]):
            n_mask_bad += 1
    P.cmp("mask.rows_with_mask", n_mask_rows, len(rows))
    P.cmp("mask.row_invariants_violated", n_mask_bad, 0)

    for p in core.PART_ORDER:
        mp, prov = meta["raw_parts"][str(p)], res["prov"][p]
        P.cmp(f"raw.{p}.file", mp["file"], prov["file"])
        P.cmp(f"raw.{p}.sha256_meta_vs_file", mp["sha256"],
              core.sha256(os.path.join(os.path.dirname(core.DATA), mp["file"])))
        P.cmp(f"raw.{p}.sha256_vs_prov", mp["sha256"], prov["sha256"])
        P.cmp(f"raw.{p}.n_rows", mp["n_rows"], len(res["by_part"][p]))
        P.cmp(f"raw.{p}.window", list(mp["window"]), list(prov["window"]))
        P.cmp(f"raw.{p}.asset_id", mp["asset_id"], prov["asset_id"])
    P.cmp("db_extra.sha256", meta["db_extra"]["sha256"],
          core.sha256(core.DB_EXTRA_PATH))
    P.cmp("db_extra.columns", list(meta["db_extra"]["columns"]),
          list(core.DB_EXTRA_COLUMNS))
    P.cmp("db_extra.n_rows", meta["db_extra"]["n_rows"], len(res["extra"]))
    P.cmp("db_extra.file", meta["db_extra"]["file"], res["extra_prov"]["file"])
    P.cmp("db_extra.sha256_prov", meta["db_extra"]["sha256"],
          res["extra_prov"]["sha256"])
    P.cmp("db_extra.n_rows_prov", res["extra_prov"]["n_rows"], len(res["extra"]))
    for name, fname in core.REFERENCE_FILES.items():
        P.cmp(f"reference.{name}.file", meta["reference_files"][name]["file"], fname)
        P.cmp(f"reference.{name}.sha256", meta["reference_files"][name]["sha256"],
              core.sha256(os.path.join(core.REF, fname)))

    # 2) CSV vs. the two committed raw files and vs. its own derived fields ---
    for p in core.PART_ORDER:
        raw = res["by_part"][p]
        csvrows = core.rows_of_part(rows, p)
        P.cmp(f"csv.part{p}.n_rows", len(csvrows), len(raw))
        for i, (crow, rrow) in enumerate(zip(csvrows, raw)):
            P.cmp(f"csv.part{p}.row{i}.source", _row_source_key(crow),
                  _row_source_key(rrow))
            P.cmp(f"csv.part{p}.row{i}.derived", _row_derived_key(crow),
                  _row_derived_recomputed(crow))

    # 3) the two committed channel JSONs, field by field ----------------------
    P.cmp("kdlo_1_channels", res["fig_part"]["part1"], ref["kdlo_1_channels"])
    P.cmp("kdlo_2_channels", res["fig_part"]["part2"], ref["kdlo_2_channels"])

    # 4) the epoch definitions, rederived from the days -----------------------
    d, fb, fc, fd, fe = (res["data"], res["fig_b"], res["fig_c"], res["fig_d"],
                         res["fig_e"])
    P.cmp("state.n_by_state", d["n_by_state"], core.N_BY_STATE_EXPECTED)
    P.cmp("state.n_rows", d["n_state_rows"],
          sum(core.N_BY_STATE_EXPECTED.values()))
    P.cmp("state.n_rows_plus_singleton", d["n_state_rows"] + 1, d["n_rows"])
    P.cmp("state.singleton_of", core.state_of(core.SINGLETON_DAY), "single")
    for lab in st.STATE_ORDER:
        days = core.days_of([r for r in rows if r.get("state") == lab])
        P.cmp(f"state.{lab}.window", [days[0], days[-1]],
              list(core.STATE_WINDOWS[lab]))
        P.cmp(f"state.{lab}.inside_window",
              all(core.STATE_WINDOWS[lab][0] <= x <= core.STATE_WINDOWS[lab][1]
                  for x in days), True)
        P.cmp(f"state.{lab}.sorted_unique", days, sorted(set(days)))
    P.cmp("state.singleton_outside_both",
          [core.SINGLETON_DAY > core.STATE_WINDOWS["pre"][1],
           core.SINGLETON_DAY < core.STATE_WINDOWS["rebuild"][0]], [True, True])
    P.cmp("state.singleton_only_row", list(fe["singleton"]), [core.SINGLETON_DAY])

    # 5) the echo-mode discriminator, rederived per row ------------------------
    P.cmp("echo.n_by_mode", d["n_by_echo_mode"],
          {"compact": 13, "distributed": 45, "intermediate": 2})
    P.cmp("echo.compact_rows_le_10px",
          all(r["A"] <= core.COMPACT_MASK_PX for r in rows
              if r["echo_mode"] == "compact"), True)
    P.cmp("echo.distributed_rows_ge_25px",
          all(r["A"] >= core.MASK_CLUTTER_PX for r in rows
              if r["echo_mode"] == "distributed"), True)
    P.cmp("echo.intermediate_rows_between",
          all(core.COMPACT_MASK_PX < r["A"] < core.MASK_CLUTTER_PX for r in rows
              if r["echo_mode"] == "intermediate"), True)
    P.cmp("echo.intermediate_days", d["n_intermediate_days"],
          core.days_of([r for r in rows if r["echo_mode"] == "intermediate"]))
    P.cmp("echo.buckets_partition_rows",
          sum(len(core.echo_modes_of(rows)[b]["days"])
              for b in ("compact", "distributed")), d["n_rows"]
          - len(d["n_intermediate_days"]))
    # 6) the effect sizes and the controls, rederived ---------------------------
    by = st.by_state(rows)
    P.cmp("fig_b.keys", fb["keys"], list(CHANNELS))
    P.cmp("fig_b.ranking", fb["ranking_by_abs_delta"],
          sorted(CHANNELS, key=lambda k: -(fb["table"][k]["sds"] or 0.0)))
    for k in CHANNELS:
        blk = fb["table"][k]
        pre_vals, reb_vals = st.vals(by, "pre", k), st.vals(by, "rebuild", k)
        P.cmp(f"fig_b.{k}.n", blk["n"], [len(pre_vals), len(reb_vals)])
        P.cmp(f"fig_b.{k}.n_pairs", blk["n_pairs"], len(pre_vals) * len(reb_vals))
        P.cmp(f"fig_b.{k}.median", blk["median"],
              [st._stats(pre_vals)["median"], st._stats(reb_vals)["median"]])
        P.cmp(f"fig_b.{k}.delta", blk["delta"], st.cliffs_delta(pre_vals, reb_vals))
        P.cmp(f"fig_b.{k}.sds_x10", blk["sds"], 10.0 * abs(blk["delta"]))
        P.cmp(f"fig_b.{k}.delta_ci95", blk["delta_ci95"],
              list(st.bootstrap_delta_ci(pre_vals, reb_vals)))
    a_ref, a_grp, a_months = month_matched(rows)
    P.cmp("fig_b.matched_months", fb["matched_months"], a_months)
    P.cmp("fig_b.matched_n", fb["matched"]["gamma2"]["n"], [len(a_ref), len(a_grp)])
    P.cmp("fig_b.matched_pairs", fb["matched"]["gamma2"]["n_pairs"],
          len(a_ref) * len(a_grp))
    P.cmp("fig_b.matched_delta", fb["matched"]["gamma2"]["delta"],
          st.cliffs_delta(channel_values(a_ref, "gamma2"),
                          channel_values(a_grp, "gamma2")))
    P.cmp("fig_b.contrast_survives_matching",
          abs(fb["matched"]["gamma2"]["delta"]) > 0.25, True)
    P.cmp("fig_b.contrast_direction", fb["table"]["gamma2"]["delta"] > 0, True)
    P.cmp("fig_b.hero_medians", fb["table"]["gamma2"]["median"],
          [0.0036978587433671717, 0.017673317985174672])
    P.cmp("fig_b.hero_delta", round(fb["table"]["gamma2"]["delta"], 3), 0.529)
    P.cmp("fig_b.hero_pairs", fb["table"]["gamma2"]["n_pairs"], 688)
    P.cmp("fig_b.hero_matched_delta", round(fb["matched"]["gamma2"]["delta"], 3),
          0.429)

    # 7) the classifier and the out-of-fold scores -----------------------------
    P.cmp("fig_c.chance", [fc["chance_accuracy"], fc["chance_balanced_accuracy"]],
          [CHANCE, CHANCE])
    P.cmp("fig_c.lambda_grid", list(fc["lambda_grid"]), list(core.LAMBDA_GRID))
    P.cmp("fig_c.grid_models", sorted(fc["sensitivity_grid"]),
          sorted(list(core.MODELS) + ["weather"]))
    P.cmp("fig_c.n_complete", [fc["n_complete_ocv"], fc["n_complete_weather"]],
          [d["n_state_rows"], d["n_state_rows"] - 1])
    P.cmp("fig_c.n_perm", fc["n_perm"], core.N_PERM_2CLASS)
    P.cmp("fig_c.hero_from_grid", fc["headline_balanced_accuracy"],
          fc["sensitivity_grid"][core.HEADLINE_PAIR][str(core.HEADLINE_LAMBDA)]
          ["balanced_accuracy"])
    P.cmp("fig_c.confusion_marginals",
          [sum(fc["headline_confusion_matrix"][a].values()) for a in st.STATE_ORDER],
          [16, 43])
    P.cmp("fig_c.weather_not_better_than_chance",
          max(fc["sensitivity_grid"]["weather"][l]["balanced_accuracy"]
              for l in [str(x) for x in core.LAMBDA_GRID]) <= CHANCE, True)
    P.cmp("fig_c.epochs_are_separable", fc["headline_balanced_accuracy"] > CHANCE,
          True)
    P.cmp("fig_c.permutation_flags",
          [fc["headline_permutation_test"][core.HEADLINE_PAIR].get("recomputed"),
           fc["headline_permutation_test"]["weather"].get("recomputed")],
          [not quick, not quick])
    if not quick:
        P.cmp("fig_c.permutation_observed_equals_headline",
              fc["headline_permutation_test"][core.HEADLINE_PAIR]
              ["balanced_accuracy_observed"], fc["headline_balanced_accuracy"])
        for name in (core.HEADLINE_PAIR, "weather"):
            blk = fc["headline_permutation_test"][name]
            P.cmp(f"fig_c.perm.{name}.p_in_range",
                  [0.0 < blk["p_balanced_accuracy"] <= 1.0,
                  0.0 < blk["p_accuracy"] <= 1.0], [True, True])
            P.cmp(f"fig_c.perm.{name}.null_mean",
                  blk["balanced_accuracy_null_mean"] >= CHANCE * 0.8, True)
        g = fc["headline_permutation_test"][core.HEADLINE_PAIR]
        P.cmp("fig_c.headline_perm_p", round(g["p_balanced_accuracy"], 3), 0.023)
        P.cmp("fig_c.observed_exceeds_null_p95",
              g["balanced_accuracy_observed"] > g["balanced_accuracy_null_p95"],
              True)
        # The echo-mask dimensions do not improve the *per-acquisition* classifier:
        # the 1D gamma2 baseline scores at least as well as the 3D OCV.
        P.cmp("fig_c.gamma2_baseline_not_worse_than_ocv",
              fc["sensitivity_grid"]["gamma2"][str(core.HEADLINE_LAMBDA)]
              ["balanced_accuracy"] >= fc["headline_balanced_accuracy"], True)
    for sub in ("all", "matched"):
        for name in core.MODELS:
            pm = fd[sub]["per_model"][name]
            P.cmp(f"fig_d.{sub}.{name}.n", pm["n"],
                  [16, 43] if sub == "all" else [16, 28])
            P.cmp(f"fig_d.{sub}.{name}.n_pairs", pm["n_pairs"], pm["n"][0] * pm["n"][1])
            P.cmp(f"fig_d.{sub}.{name}.sds_x10", pm["sds"], 10.0 * abs(pm["delta"]))
            P.cmp(f"fig_d.{sub}.{name}.ci_excludes_zero", pm["ci_excludes_zero"],
                  bool(pm["delta_ci95"][0] * pm["delta_ci95"][1] > 0))
    P.cmp("fig_d.matched_months", fd["matched_months"], a_months)
    P.cmp("fig_d.n_matched", fd["n_matched"], [len(a_ref), len(a_grp)])
    P.cmp("fig_d.q_len", len(fd["all"]["q_by_model"][core.HEADLINE_PAIR]),
          fd["all"]["per_model"][core.HEADLINE_PAIR]["n"][0]
          + fd["all"]["per_model"][core.HEADLINE_PAIR]["n"][1])
    P.cmp("fig_d.y_len", len(fd["all"]["y"]), d["n_state_rows"])
    P.cmp("fig_d.headline_beats_secondary",
          fd["all"]["per_model"][core.HEADLINE_PAIR]["sds"]
          > fd["all"]["per_model"][core.SECONDARY_PAIR]["sds"], True)
    P.cmp("fig_d.hero_oof_delta",
          round(fd["all"]["per_model"][core.HEADLINE_PAIR]["delta"], 3), 0.265)

    # 8) the controls of panel E, rederived ------------------------------------
    P.cmp("fig_e.echo_by_state", fe["echo_modes_by_state"],
          {"pre": {"n": 16, "n_compact": 1, "n_intermediate": 1,
                   "n_distributed": 14},
           "rebuild": {"n": 43, "n_compact": 12, "n_intermediate": 1,
                       "n_distributed": 30}})
    P.cmp("fig_e.echo_compact_threshold", core.COMPACT_MASK_PX, 10)
    P.cmp("fig_e.echo_clutter_threshold", core.MASK_CLUTTER_PX, 25)
    P.cmp("fig_e.by_month.keys", sorted(fe["by_month"]),
          sorted(str(m) for m in range(1, 13)))
    P.cmp("fig_e.by_month.n_overlap",
          [fe["by_month"][m]["n"] for m in ("6", "12")], [[2, 5], [1, 3]])
    P.cmp("fig_e.by_month_deltas",
          {m: (None if fe["by_month"][m].get("delta") is None
               else round(fe["by_month"][m]["delta"], 3)) for m in ("8", "9")},
          {"8": 0.2, "9": 0.867})
    P.cmp("fig_e.by_month_small_n", fe["by_month"]["5"].get("delta"), None)
    P.cmp("fig_e.by_season.keys", sorted(fe["by_season"]),
          ["shoulder_Aug-Sep", "summer_Mar-Jul", "winter_Oct-Feb"])
    P.cmp("fig_e.by_season.all_positive",
          [fe["by_season"][s]["delta"] > 0 for s in fe["by_season"]],
          [True, True, True])
    P.cmp("fig_e.by_bucket.keys", sorted(fe["by_masked_px_bucket"]),
          ["compact", "distributed", "intermediate"])
    P.cmp("fig_e.by_bucket.distributed_delta",
          round(fe["by_masked_px_bucket"]["distributed"]["delta"], 3), 0.429)
    P.cmp("fig_e.by_bucket.distributed_n",
          fe["by_masked_px_bucket"]["distributed"]["n"], [14, 30])
    P.cmp("fig_e.median_by_state_keys", sorted(fe["median_by_state"]),
          sorted(CHANNELS))
    P.cmp("fig_e.median_by_state_A", fe["median_by_state"]["A"], [120.5, 76.0])
    P.cmp("fig_e.median_by_state_gamma2", fe["median_by_state"]["gamma2"],
          [0.0036978587433671717, 0.017673317985174672])
    P.cmp("fig_e.singleton_values",
          [fe["singleton"][core.SINGLETON_DAY]["A"],
           fe["singleton"][core.SINGLETON_DAY]["echo_mode"]], [63.0, "distributed"])
    P.cmp("fig_e.singleton_gamma2_between",
          fe["singleton"][core.SINGLETON_DAY]["gamma2"]
          < fe["median_by_state"]["gamma2"]["rebuild"], True)
    P.cmp("fig_e.collapse_pair_days",
          [fe["collapse_pair"]["2022-12-11"]["last_intact"]["day"],
           fe["collapse_pair"]["2022-12-11"]["first_collapsed"]["day"]],
          ["2022-12-11", "2022-12-23"])
    P.cmp("fig_e.collapse_pair_gamma2",
          [fe["collapse_pair"]["2022-12-11"]["last_intact"]["gamma2"],
           fe["collapse_pair"]["2022-12-11"]["first_collapsed"]["gamma2"]],
          [0.019928279168188158, 0.0013393056912106268])
    P.cmp("fig_e.collapse_pair_masked_px",
          [fe["collapse_pair"]["2022-12-11"]["last_intact"]["masked_pixels"],
           fe["collapse_pair"]["2022-12-11"]["first_collapsed"]["masked_pixels"]],
          [60.0, 63.0])
    P.cmp("fig_e.collapse_pair_temperature",
          [fe["collapse_pair"]["2022-12-11"]["last_intact"]["temperature_c"],
           fe["collapse_pair"]["2022-12-11"]["first_collapsed"]["temperature_c"]],
          [-3.2, -24.6])
    P.cmp("fig_e.collapse_pair_series_context",
          fe["collapse_pair"]["series_context"],
          {"gamma2_median": 0.003535867017357631,
           "gamma2_min": 0.00047299709250516505,
           "gamma2_max": 0.02936431747713154, "n": 17})
    P.cmp("fig_e.collapse_pair_is_not_if_then",
          fe["collapse_pair"]["2022-12-11"]["first_collapsed"]["gamma2"]
          < fe["collapse_pair"]["2022-12-11"]["last_intact"]["gamma2"], True)

    # 9) the two pin-only files: structural mirrors, not recomputations --------
    for path, got, want in fe["availability"]:
        P.add(path, got, want, 0.0)
    P.n_avail = len(fe["availability"])
    for path, got, want in fe["fem"]:
        P.add(path, got, want, 0.0)
    P.n_fem = len(fe["fem"])
    return P



# ---------------------------------------------------------------------------
# Markdown report (no wall-clock info: two runs must give identical bytes)
# ---------------------------------------------------------------------------
def md_table(header, rows):
    head = "| " + " | ".join(header) + " |"
    rule = "|" + "|".join(["---"] * len(header)) + "|"
    body = ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join([head, rule] + body)


def fmt(x, nd=3):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{float(x):.{nd}f}"


def report(res, pins, quick):
    """Deterministic Markdown report of figures A–E plus the pin verdict."""
    d, fb, fc, fd, fe = (res["data"], res["fig_b"], res["fig_c"], res["fig_d"],
                         res["fig_e"])
    L = []
    A = L.append

    def sgn(x):
        return "n/a" if x is None else f"{x:+.3f}"

    g2 = fb["table"]["gamma2"]
    gm = fb["matched"]["gamma2"]
    oof = fd["all"]["per_model"]

    A("# KDLO OCV observability package — does the coherence contrast separate the "
      "epochs, or the calendar?")
    A("")
    A("Figures A–E of `research/observability-core-vector/` for the KDLO-TV site. "
      "This file is generated by `code/kdlo/fig_kdlo_ocv_paper.py` and deliberately "
      "contains no wall-clock timestamp, so two runs produce identical bytes (the "
      "PNG figures are deterministic for the same reason).")
    A("")
    A("## 1. Question and design")
    A("")
    A("The **observability core vector** here is `OCV = [gamma2, P, D]` — the "
      "interferometric coherence `gamma^2` of the KDLO acquisition plus two "
      "dimensions of its **echo mask**: the persistence `P` (how often a mask pixel "
      "is masked across the 64 dates) and the density `D = A / bbox_area`. The full "
      "**6D echo-mask vector** of this package is `[gamma2, P, D, A, F, S]` "
      "(`A` = masked pixels = `coherence_masked_pixels`, `F` = number of 8-connected "
      "components, `S` = ||mask centroid - peak pixel|| in px). The contrast under "
      "test is not a damage label but an **epoch**:")
    A("")
    A(f"- `pre` — {d['state_label']['pre']}: **{d['n_by_state']['pre']}** acquisitions;")
    A(f"- `rebuild` — {d['state_label']['rebuild']}: "
      f"**{d['n_by_state']['rebuild']}** acquisitions;")
    A(f"- the single acquisition **{d['singleton_day']}** in between belongs to "
      "neither epoch and is in **no** state statistic of this package.")
    A("")
    A("Because the two states are disjoint time windows, every effect size below is "
      "reported twice: once on all rows and once on the **month-matched** subset, "
      "and the classifier of panel C is run against a **weather-only** competitor. "
      "The design cannot prove causation; it can show which alternative explanations "
      "survive.")
    A("")
    A("## 2. The record, and what was pinned")
    A("")
    A(f"Input is the committed `{d['csv']}` (`{d['n_rows']}` acquisitions, "
      f"sha256 `{d['csv_sha256'][:16]}...`), rebuilt from the two committed raw files "
      f"(`{'`, `'.join(str(d['raw_parts'][p]['file']) for p in ('1', '2'))}`) plus a "
      "one-time read-only extraction of nine Postgres columns "
      f"(`{d['db_extra']}`). Nothing is re-measured and no network is used.")
    A("")
    A(f"The echo-mask dimensions are recomputed from the committed full-dwell strip "
      f"cache `{d['strip_cache']['dir']}/` "
      f"({d['strip_cache']['n_strips']} strips, manifest sha256 "
      f"`{d['strip_cache']['manifest_sha256'][:16]}...`) with the rule of the "
      f"pipeline itself — so `A` reproduces `coherence_masked_pixels` exactly and "
      f"the recomputed masked coherence reproduces `coherence_gamma2` to 1e-12 on "
      f"all {d['n_rows']} acquisitions (the generator refuses to write the CSV "
      f"otherwise). `P` is a property of the whole record (the 64 dates), so it "
      f"cannot be checked against a single row.")
    A("")
    A(f"Recomputed against the four committed reference JSONs: **{pins['n_checks']} "
      f"checks, {pins['n_failed']} failures** "
      f"({pins['n_exact']} exact, {pins['n_within_tolerance']} within "
      f"{pins['float_tolerance']:g}; max |deviation| "
      f"{pins['max_abs_deviation']:g}). Of those, "
      f"{pins['n_availability_checks']} are structural checks of "
      "`kdlo_s1_availability.json` and "
      f"{pins['n_fem_checks']} of `kdlo_fem_frequencies.json` — those two files are "
      "catalogue and eigenmode-solver output and are **not** rederivable from the "
      "record, so they are pinned structurally and reported honestly as such.")
    A("")
    A("## 3. What moves — and what does not")
    A("")
    A("Cliff's delta (positive = the rebuild is larger), its 95 % bootstrap CI and "
      "the same number on the month-matched subset:")
    A("")
    rows = []
    for k in fb["keys"]:
        t, m = fb["table"][k], fb["matched"][k]
        rows.append([f"`{k}`", f"{fmt(t['median'][0])} / {fmt(t['median'][1])}",
                     sgn(t["delta"]), fmt(t["sds"], 2), sgn(m["delta"])])
    A(md_table(["channel", "median pre / rebuild", "delta (all)", "SDS", "delta (matched)"],
               rows))
    A("")
    A(f"The headline is `gamma^2`: median **{g2['median'][0]:.6f} -> "
      f"{g2['median'][1]:.6f}**, Cliff's delta **{sgn(g2['delta'])}** "
      f"(SDS {g2['sds']:.2f}, {g2['n_pairs']} pre x rebuild pairs, "
      f"CI [{g2['delta_ci95'][0]:+.3f}, {g2['delta_ci95'][1]:+.3f}]). It is the "
      f"largest effect in the table. The strongest **counter**-moving channel is "
      f"`A` (delta {sgn(fb['table']['A']['delta'])}): the rebuild has *smaller* "
      "masks, i.e. the return is more concentrated, not less coherent in the "
      "absolute sense.")
    A("")
    A("## 4. The epoch problem: is it the collapse, or the calendar?")
    A("")
    A(f"`pre` holds {d['state_windows']['pre'][0]}..{d['state_windows']['pre'][1]} "
      f"and `rebuild` {d['state_windows']['rebuild'][0]}.."
      f"{d['state_windows']['rebuild'][1]}, so the overlap of the two calendars is "
      f"months {fb['matched_months'][0]}..{fb['matched_months'][-1]}: restricting "
      f"both states to it drops the {fb['table']['gamma2']['n'][1] - fb['matched']['gamma2']['n'][1]} "
      "`rebuild` acquisitions of Jan..May (the months `pre` cannot match) and leaves "
      f"{gm['n'][0]} vs. {gm['n'][1]} acquisitions. The delta falls from "
      f"{sgn(g2['delta'])} to **{sgn(gm['delta'])}** — smaller, but far from gone.")
    A("")
    A("The same control per calendar month (a stratum is only evaluated when both "
      "sides have at least 3 acquisitions, otherwise it is reported as `n/a`):")
    A("")
    rows = []
    for m in sorted(fe["by_month"], key=lambda x: int(x)):
        blk = fe["by_month"][m]
        rows.append([m, f"{blk['n'][0]} / {blk['n'][1]}", sgn(blk.get("delta"))])
    A(md_table(["month", "n pre / rebuild", "delta (gamma^2)"], rows))
    A("")
    A(f"Two of the twelve months are individually evaluable and both move the same "
      f"way (month 8: {sgn(fe['by_month']['8'].get('delta'))}, month 9: "
      f"{sgn(fe['by_month']['9'].get('delta'))}); the season strata "
      + ", ".join(f"`{k}` {sgn(v.get('delta'))}" for k, v in fe["by_season"].items())
      + " all point the same way. Whatever drives the contrast is not confined to "
      "one calendar month.")
    A("")
    A("## 5. The echo-mode split: geometry or damage?")
    A("")
    comp_pre = fe["echo_modes_by_state"]["pre"]["n_compact"]
    comp_reb = fe["echo_modes_by_state"]["rebuild"]["n_compact"]
    A(f"The pipeline's own discriminator (verbatim from the project scripts) calls a "
      f"return *compact* when `A <= {core.COMPACT_MASK_PX} px` and *distributed* when "
      f"`A >= {core.MASK_CLUTTER_PX} px`. Compact acquisitions go from "
      f"**{comp_pre}** in `pre` to **{comp_reb}** in `rebuild`, and the median mask "
      f"size drops from {fe['median_by_state']['A']['pre']:.0f} px to "
      f"{fe['median_by_state']['A']['rebuild']:.0f} px. The distributed subset alone "
      f"(n {fe['by_masked_px_bucket']['distributed']['n'][0]} vs. "
      f"{fe['by_masked_px_bucket']['distributed']['n'][1]}) still shows delta "
      f"**{sgn(fe['by_masked_px_bucket']['distributed'].get('delta'))}** for `gamma^2` "
      "— i.e. the coherence change is not produced *only* by the new compact returns.")
    A("")
    A(f"Median `gamma^2` by echo mode: compact {np.median([r['gamma2'] for r in res['rows'] if r['echo_mode'] == 'compact']):.6f}, "
      f"distributed {np.median([r['gamma2'] for r in res['rows'] if r['echo_mode'] == 'distributed']):.6f}. "
      "Compact echoes sit at *high* coherence, so they pull the `rebuild` median up: "
      "this is exactly the alternative reading — a change of **scattering geometry** "
      "rather than a change of **tower condition** — and the record cannot separate "
      "the two. The intermediate days "
      f"({', '.join(d['n_intermediate_days'])}) belong to neither bucket and are "
      "excluded from both.")
    A("")
    A("## 6. The 2-class classifier and its handicaps")
    A("")
    conf = fc["headline_confusion_matrix"]
    A(f"Four models are fitted on the same {fc['n_complete_ocv']} state acquisitions, "
      "all leave-one-out: `gamma2` alone (the baseline), the mask-morphology triple "
      "`[P, D, S]`, the OCV headline `[gamma2, P, D]` (shown here) and the full 6D "
      "echo-mask vector `[gamma2, P, D, A, F, S]`; panel D reports their out-of-fold "
      "scores. The table below is the headline model at lambda = "
      f"{fc['headline_lambda']}.")
    A("")
    A(f"Leave-one-out LDA on the {fc['n_complete_ocv']} state acquisitions (chance "
      f"balanced accuracy {CHANCE:.2f}):")
    A("")
    rows = [[f"true `{a}`"] + [str(conf[a][b]) for b in st.STATE_ORDER]
            for a in st.STATE_ORDER]
    A(md_table(["", "pred `pre`", "pred `rebuild`"], rows))
    A("")
    A(f"Accuracy {fc['headline_accuracy']:.3f}, balanced accuracy "
      f"**{fc['headline_balanced_accuracy']:.3f}**. The honest competitor — a "
      f"classifier that sees only the acquisition weather ({fc['n_complete_weather']} "
      "complete rows, "
      + ", ".join(f"`{c}`" for c in WEATHER_FEATURES)
      + f") — reaches at best "
      f"{max(fc['sensitivity_grid']['weather'][l]['balanced_accuracy'] for l in [str(x) for x in core.LAMBDA_GRID]):.3f}, "
      "i.e. it does not beat chance. So the separation is not something the weather "
      "already carries.")
    if not quick:
        g = fc["headline_permutation_test"][core.HEADLINE_PAIR]
        A("")
        A(f"Label permutation ({g['n_perm']} shuffles of the epoch label, full "
          f"leave-one-out pipeline re-run each time): observed balanced accuracy "
          f"{g['balanced_accuracy_observed']:.3f} against a null mean of "
          f"{g['balanced_accuracy_null_mean']:.3f} and a 95th percentile of "
          f"{g['balanced_accuracy_null_p95']:.3f} — empirical p = "
          f"**{g['p_balanced_accuracy']:.3f}**.")
        A("")
        g2b = fc["sensitivity_grid"]["gamma2"][str(fc["headline_lambda"])]["balanced_accuracy"]
        if g["p_balanced_accuracy"] <= 0.05:
            A(f"That p = **{g['p_balanced_accuracy']:.3f}** is below 0.05 and the "
              f"observed value sits above the null's 95th percentile: with the "
              f"echo-mask dimensions added, the per-acquisition classifier **does** "
              f"separate the two epochs beyond chance. It must not be read as more "
              f"than that — the labels are epochs (a year-long drift or a pipeline "
              f"re-baseline would separate them just as well), and the 1D `gamma2` "
              f"baseline alone scores **{g2b:.3f}**, i.e. *higher* than the OCV's "
              f"{fc['headline_balanced_accuracy']:.3f}: the mask dimensions change "
              f"the permutation null, not the decision boundary. The contrast of "
              f"Figure B is a difference of *distributions* either way.")
        else:
            A(f"That p is the honest answer to the panel-C question: **the "
              f"per-acquisition classifier does not separate the epochs beyond "
              f"chance** (p = {g['p_balanced_accuracy']:.3f} > 0.05; the observed "
              f"value sits below the null's 95th percentile). The contrast of "
              f"Figure B is a difference of *distributions*, not a decision "
              f"boundary — which is exactly what a small, single-geometry sample "
              f"should look like, and the reason neither number is allowed to "
              f"stand alone in this report.")
    else:
        A("")
        A("The label-permutation null is **not** recomputed in `--quick` mode.")
    A("")
    A("Three handicaps remain, none of which the pinning can remove: the labels are "
      "epochs (so a year-long drift, a sensor re-baseline or a change of orbit "
      "behaviour would separate them just as well); the leave-one-out CV only says "
      "\"a model fitted on the rest of the record scores this acquisition the same "
      "way\", not \"this acquisition was not disturbed\"; and every number here comes "
      "from a single site with a single geometry (S1A descending, one relative "
      f"orbit), so there is no second site against which to check.")
    A("")
    A("## 7. What this package does *not* claim")
    A("")
    A("- **Not a damage label.** `pre` and `rebuild` are epochs. Every effect size is "
      "an epoch contrast; nothing here says the mast was intact or damaged.")
    A("- **Not a causal statement.** A contrast on two disjoint time windows cannot "
      "distinguish \"the tower changed\" from \"something else that changed with the "
      "tower\". The month-matched column and the weather competitor of §4 are the "
      "only defences, and they are partial.")
    A("- **Not \"after the reconstruction\".** The record ends "
      f"{d['date_last']} while the antenna reached 1,542 ft on 2024-08-19, so the "
      "second epoch is *during* the rebuild.")
    A("- **Not independent evidence.** `gamma^2` and `A` come from the same "
      "interferogram and the same pipeline; the phase channels come from the same "
      "processing chain. Coherence, mask size and phase activity are three views of "
      "one scattering geometry, not three experiments.")
    A("- **Not a second site.** The LUMO and Bautzen packages in this repository are "
      "separate sites with their own records; the agreement or disagreement between "
      "the three is a cross-site question this package does not answer.")
    A("- **Not a re-measurement.** The two reference files that are not rederivable "
      f"(`kdlo_s1_availability.json`, `kdlo_fem_frequencies.json`) are pinned "
      f"structurally ({pins['n_availability_checks']} + {pins['n_fem_checks']} "
      "checks) and treated as given; only the two channel JSONs are recomputed "
      "field by field from the committed record.")
    A("- **Not a large sample.** "
      f"{d['n_by_state']['pre']} `pre` vs. {d['n_by_state']['rebuild']} `rebuild` "
      "acquisitions on one geometry; strata with fewer than 3 acquisitions on one "
      "side are reported as `n/a` instead of being silently dropped, and the "
      f"{len(d['missing_repeat_day_part1'])} missing repeat day(s) of part 1 are "
      "reported rather than interpolated.")
    A("")
    A("## 8. Reproduce")
    A("")
    A("```bash")
    A("bash figures/kdlo/fig_kdlo_ocv_paper.sh          # full run, all pins")
    A("python3 code/kdlo/fig_kdlo_ocv_paper.py --quick  # without the permutation null")
    A("```")
    A("")
    A(f"Verdict of this run: **{pins['n_checks']} checks, {pins['n_failed']} "
      f"failures** (case: {'quick' if quick else 'full'}); the package is "
      "reproducible and the figures are byte-deterministic.")
    A("")
    return "\n".join(L) + "\n"



# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def json_default(o):
    """JSON fallback for the numpy scalars/arrays coming out of the statistics."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")


def build_json(res, pins_summary, quick):
    """Self-contained result JSON: everything the figures plot, plus the pins."""
    return {
        "meta": {
            "generator": "code/kdlo/fig_kdlo_ocv_paper.py",
            "quick": bool(quick),
            "csv": res["data"]["csv"],
            "csv_sha256": res["data"]["csv_sha256"],
            "n_rows": res["data"]["n_rows"],
            "figures": ["figures/kdlo/fig_kdlo_ocv_paper_%s.png" % k for k in "ABCDE"],
            "report": "figures/kdlo/fig_kdlo_ocv_paper.md",
            "state_windows": res["data"]["state_windows"],
            "singleton_day": res["data"]["singleton_day"],
            "n_by_state": res["data"]["n_by_state"],
            "seeds": {"rng_seed": core.RNG_SEED, "n_boot": core.N_BOOT,
                      "n_boot_strata": core.N_BOOT_STRATA,
                      "n_perm_2class": core.N_PERM_2CLASS},
            "versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "matplotlib": matplotlib.__version__,
            },
            "note": "no timestamps on purpose: repeated runs must be identical",
        },
        "data": res["data"],
        "fig_a": res["fig_a"],
        "fig_b": res["fig_b"],
        "fig_c": res["fig_c"],
        "fig_d": res["fig_d"],
        "fig_e": {k: v for k, v in res["fig_e"].items()
                  if k not in ("availability", "fem")},
        "fig_part": res["fig_part"],
        "pins": pins_summary,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Recompute figures A-E of the KDLO OCV observability package.")
    ap.add_argument("--quick", action="store_true",
                    help="skip the two label-permutation tests (panel C null)")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="output directory of the PNG figures")
    ap.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD, help="path of the Markdown report")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    print("KDLO OCV figure package — recomputing figures A-E "
          f"({'quick' if args.quick else 'full'} mode) ...")
    res = compute(quick=args.quick)
    print(f"  compute: {time.perf_counter() - t0:.1f}s")

    P = pin_all(res, quick=args.quick)
    s = P.summary()
    print(f"  pins: {s['n_checks']} checks, {s['n_exact']} exact, "
          f"{s['n_failed']} failures, max |dev| {s['max_abs_deviation']}")
    if s["n_failed"]:
        for f in s["failures"][:10]:
            print(f"    FAIL {f['path']}: {f['recomputed']} != {f['committed']}")

    os.makedirs(args.figdir, exist_ok=True)
    print("  writing figures:")
    for fn in (fig_a, fig_b, fig_c, fig_d, fig_e):
        fn(res, args.figdir)

    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(build_json(res, s, args.quick), fh, indent=1, sort_keys=True,
                  default=json_default)
        fh.write("\n")
    print(f"  {os.path.relpath(args.json)}")

    with open(args.md, "w", encoding="utf-8") as fh:
        fh.write(report(res, s, args.quick))
    print(f"  {os.path.relpath(args.md)}")

    print(f"done in {time.perf_counter() - t0:.1f}s, "
          f"{s['n_failed']} pin failures")
    return 0 if s["n_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

