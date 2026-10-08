#!/usr/bin/env python3
"""YWF OCV-Paper — figures A-E, JSON result and English report.

Question: does the observability core vector ``OCV = [gamma2, P, D]`` of the
Yeongdeok Wind Farm (Changpo Wind Power Complex, Samgye-ri) — the
interferometric coherence gamma^2 of the echo of Unit 21, plus the persistence
``P`` and the density ``D`` of that echo mask — separate the two states of the
site, the steel monopile tower **before** its 2026-02-02 collapse
(``pre-collapse (healthy)``) and **after** it (``post-collapse``)?

What this package claims, and what it does not
----------------------------------------------
  * **one small record, two states.** ``pre-collapse (healthy)`` and
    ``post-collapse`` are the two sides of the collapse date of the *same*
    tower, so season (pre Aug-Jan, post Feb-Apr), weather and pipeline
    generation all move with the state. Every contrast is reported beside a
    co-variate control — never as damage.
  * **the post side cannot carry a contrast at all.** The record holds 6 unique
    chips after the event, of which **1** carries an echo; before it, 16 of 27.
    The headline of this figure is therefore an **echo-coverage** statement with
    an explicit power caveat (panel B), not a damage verdict: the two states are
    not separable here, and the figure says so.
  * **the 6D vector exists only where the mask rule holds.** ``A`` (masked
    pixels), ``P``, ``D``, ``F``, ``S`` and the masked ``gamma2`` are undefined on
    a chip without echo, so the channel panels rest on 17 of the 33 chips.
  * **one mask layer.** The committed record is *not* segment-resolved: 24 of its
    57 windows are byte-identical duplicates that differ only in the requested
    mast section (2 vs. 3), so a per-mast-section layer would be dishonest and
    this package carries the echo mask of the whole 7x7 window only (panel E d).
  * **``A`` is not the pipeline's ``coherence_masked_pixels``.** That column was
    computed on the different-generation asset-point chip (75, 133, ... px) while
    this package recomputes the echo mask of the committed 7x7 window (2..7 px);
    the two coincide on **0** of 17 rows, which is reported, not hidden (panel E a).

Reference file (``data/ywf/reference/``)
----------------------------------------
  ywf_ocv_findings.json   the committed analysis-layer reference. This site has
                          **no upstream analysis project** to pin against (the
                          record *is* the site's own export), so the reference is
                          *defined* by this package: ``--write-reference`` writes
                          it, and every later run re-derives it field by field
                          and aborts on the first deviation.

Figures:
  A  raw distributions of the six channels per state (the 17 echo-bearing chips)
  B  the headline: echo coverage per state + its Fisher exact test, and the
     coverage timeline across the collapse
  C  in-sample effect sizes: Cliff's delta + bootstrap CI and SDS, for the mask
     channels and for the co-variate controls beside them
  D  the OCV plane (gamma2, P, D): 16 pre chips against 1 post chip — no
     classifier is fitted, because one echo chip cannot be one
  E  why this contrast is not a damage signal: the export-vs-recomputed mask,
     the echo-mode mix, the season composition and the payload duplication

Usage:
  python3 code/ywf/fig_ywf_ocv_paper.py                    # figures + pins
  python3 code/ywf/fig_ywf_ocv_paper.py --write-reference  # (re)define the
                                                           # committed reference
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
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import ywf_ocv_core as core    # noqa: E402
import ywf_ocv_masks as masks  # noqa: E402
import ywf_ocv_stats as st     # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_ywf_ocv_paper.json")
OUT_MD = os.path.join(FIGDIR, "fig_ywf_ocv_paper.md")
REF_FINDINGS = core.REF_FINDINGS

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

STATE_COLORS = {core.PRE: "#4C72B0", core.POST: "#C44E52"}
NO_ECHO_COLOR = "#B0B0B0"
CHANNEL_COLORS = {"gamma2": "#4C72B0", "P": "#55A868", "D": "#8172B3",
                  "A": "#937860", "F": "#64B5CD", "S": "#DA8BC3",
                  "wind_speed_ms": "#4C72B0", "temperature_c": "#C44E52",
                  "precipitation_mm": "#55A868", "wind_gust_ms": "#8172B3",
                  "brightness_ratio": "#937860", "intensity": "#64B5CD",
                  "coherence": "#DA8BC3", "coherence_gamma2": "#BCBD22",
                  "cumulative_mm": "#7F7F7F",
                  "structural_frequency_hz": "#17BECF"}
# The mask channels and, beside them, the columns that move with the state for
# reasons other than the collapse: the pipeline's own observables and weather.
CHANNELS = list(core.FEATURES) + list(core.CONTROLS)
ECHO_MODES = ["compact", "intermediate", "distributed"]



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


def _jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def day_date(day):
    """``YYYY-MM-DD`` -> ``datetime.date``, for the time axes."""
    return dt.date.fromisoformat(str(day))


# ---------------------------------------------------------------------------
# Row views
# ---------------------------------------------------------------------------
def mask_rows(csv_rows):
    """The rows as this package analyses them (``A`` is ``None`` without echo)."""
    out = []
    for r in csv_rows:
        d = core.row_brief(r)
        # row_brief carries the weather subset; the pipeline's own observables and
        # the remaining controls have to be added for the co-variate panels.
        for c in ("precipitation_mm", "wind_gust_ms", "relative_humidity_pct",
                  "brightness_ratio", "intensity", "coherence",
                  "coherence_gamma2", "cumulative_mm",
                  "structural_frequency_hz"):
            d[c] = core._num(r.get(c))
        d["section_index"] = r.get("section_index")
        d["month"] = core._num(r.get("month"))
        d["day"] = core.day_of(r)
        d["echo_mode"] = r.get("echo_mode")
        d["A_export"] = core._num(r.get("coherence_masked_pixels"))
        d["gamma2_export"] = core._num(r.get("coherence_gamma2"))
        out.append(d)
    return out


def _vals(rows, key):
    return [r[key] for r in rows if r.get(key) is not None]


def delta_pair(rows, key, ref=core.PRE, grp=core.POST,
               n_boot=core.N_BOOT, seed=core.RNG_SEED):
    """``st.delta_block`` for one column between the two states."""
    a = _vals([r for r in rows if r.get("state") == ref], key)
    b = _vals([r for r in rows if r.get("state") == grp], key)
    return st.delta_block(a, b,
                          f"{core.STATE_SHORT[ref]} -> {core.STATE_SHORT[grp]}",
                          n_boot=n_boot, seed=seed)


# ---------------------------------------------------------------------------
# Aggregate blocks (the reference file is built from these)
# ---------------------------------------------------------------------------
def timeline_block(windows, unique):
    """One entry per acquisition day: committed windows, unique chips, echoes."""
    out = []
    for day in sorted({c["date"] for c in windows}):
        w = [c for c in windows if c["date"] == day]
        u = [c for c in unique if c["date"] == day]
        out.append({
            "day": day, "state": core.STATE_SHORT[core.state_of(day)],
            "n_windows": len(w), "n_unique_chips": len(u),
            "n_duplicated": len(w) - len(u),
            "n_with_echo": sum(1 for c in u if c["mask"] is not None),
        })
    return out


def payload_groups(windows):
    """The payload-digest groups of the committed windows (the duplication)."""
    groups = {}
    for c in windows:
        g = groups.setdefault(c["digest"], {"segments": set(), "n_windows": 0,
                                            "days": set(), "with_echo": False})
        g["segments"].add(c["section_index"])
        g["n_windows"] += 1
        g["days"].add(c["date"])
        g["with_echo"] = g["with_echo"] or c["mask"] is not None
    return [{"digest": d, "segments": sorted(v["segments"]),
             "n_windows": v["n_windows"], "days": sorted(v["days"]),
             "with_echo": bool(v["with_echo"])}
            for d, v in sorted(groups.items())]


def echo_mode_block(rows):
    """Echo-mode mix per state over the echo-bearing chips."""
    out = {}
    for s in core.STATES:
        sub = _vals([r for r in rows if r.get("state") == s], "A")
        counts = {m: sum(1 for a in sub if core.echo_mode_of(a) == m)
                  for m in ECHO_MODES}
        out[core.STATE_SHORT[s]] = {
            "n_echo_chips": len(sub), "counts": counts,
            "fractions": {m: (c / len(sub)) if sub else None
                          for m, c in counts.items()},
        }
    return out


def composition_block(rows):
    """Rows per month, per season and per mast section, per state."""
    out = {"months": {}, "seasons": {}, "sections": {}}
    for s in core.STATES:
        sub = [r for r in rows if r.get("state") == s]
        out[core.STATE_SHORT[s]] = {"n": len(sub)}
        for r in sub:
            month = f"{int(r['month']):02d}" if r.get("month") else "?"
            out["months"].setdefault(month, {})
            out["months"][month][core.STATE_SHORT[s]] = \
                out["months"][month].get(core.STATE_SHORT[s], 0) + 1
    for lab, sub in core.by_section(rows).items():
        if not sub:                      # the section table lists 0..4, but the
            continue                     # record only decoded mast section 3
        out["sections"][lab] = {
            "n": len(sub),
            "n_pre": sum(1 for r in sub if r["state"] == core.PRE),
            "n_post": sum(1 for r in sub if r["state"] == core.POST),
            "n_with_echo": sum(1 for r in sub if r.get("A") is not None),
        }
    return out


def post_chips_needed(n_pre, k_pre, alpha=0.05, nmax=400):
    """Post-event chips needed for a significant Fisher result if *none* echoes.

    The echo rate before the event is known from the record; the question this
    answers is how many post-event chips it would take to reject "the coverage
    is unchanged" at ``alpha`` when all of them turn out echo-free. ``None``
    when even ``nmax`` chips would not do it.
    """
    if st.sps is None:                                  # pragma: no cover
        return None
    for n in range(1, nmax + 1):
        p = core.fisher_p(k_pre, n_pre - k_pre, 0, n)
        if p is not None and p < alpha:
            return n
    return None


def compute():
    t0 = time.time()
    rows = mask_rows(core.load_csv())
    windows, unique = core.chip_table()
    coverage = core.coverage_block(unique)
    manifest = core.load_manifest()
    meta = core.load_csv_meta()
    by_state = st.by_state(rows)
    echo_rows = [r for r in rows if r.get("A") is not None]

    channels = {k: delta_pair(rows, k) for k in CHANNELS
                if any(r.get(k) is not None for r in rows)}
    cov = coverage["per_state"]
    n_pre, k_pre = cov["pre"]["n_chips"], cov["pre"]["n_with_echo"]
    power = {
        "n_pre_chips": n_pre, "n_pre_with_echo": k_pre,
        "n_post_chips": cov["post"]["n_chips"],
        "n_post_with_echo": cov["post"]["n_with_echo"],
        "echo_rate_pre": cov["pre"]["echo_rate"],
        "echo_rate_post": cov["post"]["echo_rate"],
        "fisher_p_two_sided": coverage["fisher"]["p_two_sided"],
        "fisher_alpha": 0.05,
        "post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha":
            post_chips_needed(n_pre, k_pre),
        "note": ("1 of 6 post-event chips carries an echo, so there is no "
                 "post-state distribution: the collapse contrast of this site "
                 "is an echo-coverage statement with this caveat attached, "
                 "not a damage verdict."),
    }
    return {
        "rows": rows,
        "echo_rows": echo_rows,
        "windows": windows,
        "unique": unique,
        "coverage": coverage,
        "timeline": timeline_block(windows, unique),
        "payload_groups": payload_groups(windows),
        "echo_modes": echo_mode_block(rows),
        "composition": composition_block(rows),
        "channels": channels,
        "controls_by_state": {core.STATE_SHORT[s]: {
            k: st._stats(_vals(by_state[s], k)) for k in core.CONTROLS}
            for s in core.STATES},
        "controls_delta": {k: delta_pair(rows, k) for k in core.CONTROLS
                           if any(r.get(k) is not None for r in rows)},
        "power": power,
        "meta": meta,
        "manifest": manifest,
        "n_rows": len(rows),
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# The committed reference file
# ---------------------------------------------------------------------------
def json_safe(o):
    """Non-finite floats (``nan``/``inf``) become ``None``.

    A test that cannot be computed — here the Welch/Mann-Whitney pair of a
    control whose two states are nearly identical — is stored as JSON ``null``
    rather than as the non-standard ``NaN`` literal, and compared as such.
    """
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if isinstance(o, np.floating):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def echo_mask_by_state(rows):
    """``{state_short: {channel: st._stats}}`` over the echo-bearing chips."""
    return {core.STATE_SHORT[s]: {
        k: st._stats(_vals([r for r in rows if r["state"] == s], k))
        for k in core.FEATURES} for s in core.STATES}


def findings(res):
    """The committed analysis-layer reference (``data/ywf/reference/``).

    Everything in here follows from the committed tables alone, so the file is a
    *definition* of this package's analysis layer: written once by
    ``--write-reference``, re-derived and pinned field by field on every later
    run.
    """
    man, meta = res["manifest"], res["meta"]
    return {
        "site": core.SITE,
        "asset_id": core.ASSET_ID,
        "asset_name": core.ASSET_NAME,
        "event": core.COLLAPSE_DATE,
        "state_labels": list(core.STATES),
        "state_short": dict(core.STATE_SHORT),
        "window_px": 49,
        "ocv": list(core.OCV),
        "features": list(core.FEATURES),
        "mask_rule": man["mask_rule"],
        "echo_mode_cuts": core.ECHO_MODE_CUTS,
        "records": {
            "n_windows": man["n_rows"] + man["n_deduplicated_same_payload"],
            "n_unique_chips": man["n_rows"],
            "n_dates": len(res["timeline"]),
            "n_deduplicated_same_payload": man["n_deduplicated_same_payload"],
            "segment_resolved": man["segment_resolved"],
            "window_sizes": man["window_sizes"],
            "sections_present": list(res["composition"]["sections"]),
            "asset_date_groups": man["asset_date_groups"],
            "asset_date_groups_with_several_segments":
                man["asset_date_groups_with_several_segments"],
            "asset_date_groups_with_differing_payloads":
                man["asset_date_groups_with_differing_payloads"],
            "asset_date_groups_with_repeated_segment":
                man["asset_date_groups_with_repeated_segment"],
        },
        "echo_coverage": {
            "n_unique_chips": res["coverage"]["n_windows"],
            "per_state": res["coverage"]["per_state"],
            "fisher": res["coverage"]["fisher"],
            "power": res["power"],
        },
        "echo_mask_by_state": echo_mask_by_state(res["rows"]),
        "echo_modes_by_state": res["echo_modes"],
        "pipeline_column": {
            "name": "coherence_masked_pixels",
            "n_compared": meta["n_coherence_masked_pixels_compared"],
            "n_matched": meta["n_coherence_masked_pixels_matched"],
            "note": ("a different generation: the column was computed on the "
                     "asset-point chip of the whole-scene window, while this "
                     "package recomputes the echo mask of the committed 7x7 "
                     "window"),
            "export_by_state": {core.STATE_SHORT[s]: st._stats(
                _vals([r for r in res["rows"] if r["state"] == s], "A_export"))
                for s in core.STATES},
            "recomputed_by_state": {core.STATE_SHORT[s]: st._stats(
                _vals([r for r in res["rows"] if r["state"] == s], "A"))
                for s in core.STATES},
        },
        "controls_by_state": res["controls_by_state"],
        "controls_delta": res["controls_delta"],
        "months_by_state": res["composition"]["months"],
        "sections": res["composition"]["sections"],
        "timeline": res["timeline"],
        "payload_groups": res["payload_groups"],
        "sources": {
            "windows": {"path": "data/ywf/ywf_windows_full.txt",
                        "sha256": core.sha256(core.WINDOWS_PATH)},
            "mask_cache": {"path": "data/ywf/ywf_windows_mask_cache.txt",
                           "sha256": core.sha256(core.CACHE_PATH)},
            "mask_cache_manifest": {
                "path": "data/ywf/ywf_windows_mask_cache.manifest.json",
                "sha256": core.sha256(core.MANIFEST_PATH)},
            "measurements": {"path": "data/ywf/ywf_measurements_full.txt",
                             "sha256": core.sha256(core.MEAS_PATH)},
            "segments": {"path": "data/ywf/ywf_segments.txt",
                         "sha256": core.sha256(core.SEGMENTS_PATH)},
            "channels_csv": {"path": "data/ywf/ywf_ocv_channels.csv",
                             "sha256": core.sha256(core.CSV_PATH)},
        },
        "figures": [f"figures/ywf/fig_ywf_ocv_paper_{L}.png" for L in "ABCDE"],
    }




# ---------------------------------------------------------------------------
# The figures A-E
# ---------------------------------------------------------------------------
def fig_A(res, out_dir=None):
    """Raw distributions of the six channels per state (echo-bearing chips)."""
    rows = res["echo_rows"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 0.62 * COL_WIDTH))
    for i, k in enumerate(core.FEATURES):
        ax = axes[i // 3][i % 3]
        for j, s in enumerate(core.STATES):
            vals = _vals([r for r in rows if r["state"] == s], k)
            if not vals:
                continue
            if len(vals) >= 3:
                bp = ax.boxplot([vals], positions=[j], widths=0.5,
                                showfliers=False, patch_artist=True,
                                medianprops=dict(color="white", lw=1.1))
                bp["boxes"][0].set(facecolor=STATE_COLORS[s], alpha=0.65,
                                   edgecolor=STATE_COLORS[s])
                ax.plot(j + _jitter(len(vals), 0.14, 11 + i), vals, ".",
                        ms=2.3, alpha=0.5, color=STATE_COLORS[s],
                        markeredgewidth=0)
            else:
                ax.plot([j] * len(vals), vals, "*", ms=10.0,
                        color=STATE_COLORS[s], markeredgewidth=0)
            ax.text(j, 0.99, f"n={len(vals)}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=6.4, color=STATE_COLORS[s])
        ax.set_xticks(range(len(core.STATES)))
        ax.set_xticklabels([core.STATE_SHORT[s] for s in core.STATES])
        ax.set_ylabel(core.CHANNEL_TITLE.get(k, k))
        if k in ("A", "D", "S"):
            ax.set_yscale("log")
        panel(ax, "abcdef"[i], core.CHANNEL_AXIS.get(k, k))
    n_pre = len([r for r in rows if r["state"] == core.PRE])
    fig.suptitle(f"(YWF) echo-mask channels per state — the "
                 f"{res['n_rows']} chips of the record, the {len(rows)} of them "
                 f"that carry an echo ({n_pre} pre, {len(rows) - n_pre} post)",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return save_fig(fig, "fig_ywf_ocv_paper_A", out_dir)


def fig_B(res, out_dir=None):
    """The headline: echo coverage per state, and the coverage timeline."""
    cov = res["coverage"]["per_state"]
    fis = res["coverage"]["fisher"]
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.40 * COL_WIDTH),
                             gridspec_kw={"width_ratios": [0.72, 1.6]})

    ax = axes[0]
    for i, s in enumerate(core.STATES):
        b = cov[core.STATE_SHORT[s]]
        ax.bar(i, b["n_with_echo"], 0.55, color=STATE_COLORS[s], alpha=0.85,
               label="carries an echo" if i == 0 else None)
        ax.bar(i, b["n_without_echo"], 0.55, bottom=b["n_with_echo"],
               color=NO_ECHO_COLOR, alpha=0.55,
               label="no echo (A undefined)" if i == 0 else None)
        ax.text(i, b["n_chips"] + 0.5,
                f"{b['n_with_echo']}/{b['n_chips']} = "
                f"{(b['echo_rate'] or 0.0) * 100:.0f} %",
                ha="center", va="bottom", fontsize=6.6,
                color=STATE_COLORS[s])
    p = fis["p_two_sided"]
    ax.set_xticks(range(len(core.STATES)))
    ax.set_xticklabels([core.STATE_SHORT[s] for s in core.STATES])
    ax.set_ylabel("unique chips")
    ax.set_ylim(0, 40)
    ax.legend(frameon=False, loc="upper right", fontsize=6.4)
    panel(ax, "a", "echo coverage per state\n"
                   f"(Fisher exact p = {'n/a' if p is None else f'{p:.3f}'})")

    ax = axes[1]
    for e in res["timeline"]:
        xd = mdates.date2num(day_date(e["day"]))
        col = STATE_COLORS[core.PRE if e["state"] == "pre" else core.POST]
        if e["n_with_echo"]:
            ax.plot(xd, 1.0, "o", ms=3.5 + 3.0 * e["n_with_echo"], color=col,
                    alpha=0.90, markeredgewidth=0)
        n_no = e["n_unique_chips"] - e["n_with_echo"]
        if n_no:
            ax.plot(xd, 0.0, "o", ms=3.5 + 3.0 * n_no, color=NO_ECHO_COLOR,
                    alpha=0.75, markeredgecolor=col, markeredgewidth=0.8)
    ax.axvline(mdates.date2num(day_date(core.COLLAPSE_DATE)), color="#333333",
               lw=1.0, ls="--")
    ax.text(mdates.date2num(day_date(core.COLLAPSE_DATE)), 1.62,
            f" collapse {core.COLLAPSE_DATE}", fontsize=6.4, color="#333333")
    ax.set_yticks([0.0, 1.0])
    ax.set_yticklabels(["no echo", "echo"])
    ax.set_ylim(-0.45, 1.85)
    ax.set_xlim(mdates.date2num(day_date(res["timeline"][0]["day"])) - 8,
                mdates.date2num(day_date(res["timeline"][-1]["day"])) + 8)
    ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=9))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    for lab in ax.get_xticklabels():
        lab.set_rotation(28)
        lab.set_ha("right")
    panel(ax, "b", "the coverage timeline \u2014 marker area = chips of the day")
    fig.suptitle("(YWF) the headline is echo coverage, not a damage verdict",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_ywf_ocv_paper_B", out_dir)


def _delta_panel(ax, blocks, keys, fs=6.0):
    """Cliff's delta + bootstrap CI per column, one row each."""
    keys = [k for k in keys if "delta" in (blocks.get(k) or {})]
    labels = []
    for i, k in enumerate(keys):
        b = blocks.get(k) or {}
        col = CHANNEL_COLORS.get(k, "#333333")
        d = b["delta"]
        ax.plot([d], [i], "o", ms=5.0, color=col)
        if "delta_ci95" in b:
            lo, hi = b["delta_ci95"]
            ax.plot([lo, hi], [i, i], "-", lw=1.4, color=col)
            txt, tcol = (("CI \u2260 0", "#444444") if b.get("ci_excludes_zero")
                         else ("CI \u220b 0", "#AA2222"))
        else:
            txt, tcol = (f"no CI (n_post = {b['n'][1]})", "#AA2222")
        ax.text(d, i - 0.32, txt, fontsize=fs, color=tcol, ha="center",
                va="center")
        labels.append(f"{k}  (SDS {b['sds']:.1f})")
    ax.axvline(0.0, color="#888888", lw=0.9)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.margins(y=0.16)
    ax.set_xlabel("Cliff's delta  (post \u2212 pre)")


def fig_C(res, out_dir=None):
    """In-sample effect sizes: the mask channels and the co-variate controls."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.44 * COL_WIDTH))
    _delta_panel(axes[0], res["channels"], list(core.FEATURES))
    panel(axes[0], "a", "the six mask channels")
    _delta_panel(axes[1], res["controls_delta"], list(core.CONTROLS))
    panel(axes[1], "b", "the co-variate controls")
    fig.suptitle("(YWF) in-sample effect sizes \u2014 the mask channels (16 pre / 1 "
                 "post echo chip) beside the controls that move with the state",
                 fontsize=8.6, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return save_fig(fig, "fig_ywf_ocv_paper_C", out_dir)


def fig_D(res, out_dir=None):
    """The OCV plane (gamma2, P, D) per state - 16 pre chips against 1 post."""
    rows = res["echo_rows"]
    pairs = [("gamma2", "P"), ("gamma2", "D"), ("P", "D")]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 0.34 * COL_WIDTH))
    for i, (xa, ya) in enumerate(pairs):
        ax = axes[i]
        for s in core.STATES:
            sel = [r for r in rows if r["state"] == s
                   and r.get(xa) is not None and r.get(ya) is not None]
            if not sel:
                continue
            pre = s == core.PRE
            ax.plot([r[xa] for r in sel], [r[ya] for r in sel],
                    "o" if pre else "*", ms=3.0 if pre else 10.0,
                    alpha=0.55 if pre else 1.0, color=STATE_COLORS[s],
                    markeredgewidth=0, label=core.STATE_SHORT[s])
        ax.set_xlabel(core.CHANNEL_TITLE.get(xa, xa))
        ax.set_ylabel(core.CHANNEL_TITLE.get(ya, ya))
        panel(ax, "abc"[i], f"{core.CHANNEL_TITLE.get(xa, xa)} vs "
                            f"{core.CHANNEL_TITLE.get(ya, ya)}")
    axes[0].legend(frameon=False, fontsize=6.8)
    axes[2].text(0.03, 0.03, "1 post-event echo chip:\nno classifier is fitted",
                 transform=axes[2].transAxes, fontsize=6.4, color="#AA2222")
    fig.suptitle("(YWF) the OCV plane \u2014 one post-event echo chip against "
                 "sixteen pre-event ones", fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return save_fig(fig, "fig_ywf_ocv_paper_D", out_dir)


def fig_E(res, out_dir=None):
    """Why this contrast is not a damage signal: the four controls."""
    rows = res["rows"]
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, 0.62 * COL_WIDTH))
    width = 0.38

    # (a) the pipeline's own masked-pixel column vs. the recomputed echo A
    ax = axes[0][0]
    sel = [r for r in rows if r.get("A_export") and r.get("A") is not None]
    for s in core.STATES:
        sub = [r for r in sel if r["state"] == s]
        if not sub:
            continue
        ax.plot([r["A_export"] for r in sub], [r["A"] for r in sub], "o", ms=3.2,
                alpha=0.65, color=STATE_COLORS[s], markeredgewidth=0,
                label=core.STATE_SHORT[s])
    if sel:
        lo = min(min(r["A_export"], r["A"]) for r in sel)
        hi = max(max(r["A_export"], r["A"]) for r in sel)
        ax.plot([lo, hi], [lo, hi], "-", color="#888888", lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
    meta = res["meta"]
    ax.set_xlabel("A export = coherence_masked_pixels")
    ax.set_ylabel("A recomputed (7x7 echo mask)")
    ax.legend(frameon=False, loc="lower right", fontsize=6.8)
    ax.text(0.03, 0.03, f"{meta['n_coherence_masked_pixels_matched']} of "
                        f"{meta['n_coherence_masked_pixels_compared']} rows "
                        f"coincide", transform=ax.transAxes, fontsize=6.4,
            color="#AA2222")
    panel(ax, "a", "export vs. recomputed mask")

    # (b) echo-mode mix per state
    ax = axes[0][1]
    modes = list(ECHO_MODES)
    for si, s in enumerate(core.STATES):
        b = res["echo_modes"][core.STATE_SHORT[s]]
        vals = [b["fractions"][m] or 0.0 for m in modes]
        ax.bar([i + (si - 0.5) * width for i in range(len(modes))], vals,
               width=width, color=STATE_COLORS[s], alpha=0.80,
               label=f"{core.STATE_SHORT[s]} (n={b['n_echo_chips']})")
    ax.set_xticks(range(len(modes)))
    ax.set_xticklabels(modes)
    ax.set_ylabel("fraction of echo chips")
    ax.legend(frameon=False, fontsize=6.8)
    panel(ax, "b", "echo-mode mix of the echo chips")

    # (c) the season moves with the state
    ax = axes[1][0]
    months = res["composition"]["months"]
    keys = sorted(months)
    for si, s in enumerate(core.STATES):
        short = core.STATE_SHORT[s]
        vals = [months[m].get(short, 0) for m in keys]
        ax.bar([i + (si - 0.5) * width for i in range(len(keys))], vals,
               width=width, color=STATE_COLORS[s], alpha=0.80, label=short)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys)
    ax.set_xlabel("month of the acquisition")
    ax.set_ylabel("committed windows")
    ax.legend(frameon=False, fontsize=6.8)
    ax.text(0.985, 0.97, "pre Aug-Jan\npost Feb-Apr", transform=ax.transAxes,
            fontsize=6.2, color="#444444", ha="right", va="top")
    panel(ax, "c", "the season co-moves")

    # (d) the duplication: windows per day vs. unique chips per day
    ax = axes[1][1]
    tl = res["timeline"]
    xs = np.arange(len(tl))
    ax.bar(xs, [e["n_windows"] for e in tl], width=0.8, color="#C9C9C9",
           label="committed windows")
    ax.bar(xs, [e["n_duplicated"] for e in tl], width=0.8, color="#8C8C8C",
           label="byte-identical duplicates")
    ax.plot(xs, [e["n_unique_chips"] for e in tl], "k.", ms=3.0,
            label="unique chips")
    n_dup = sum(e["n_duplicated"] for e in tl)
    n_win = sum(e["n_windows"] for e in tl)
    ax.set_ylim(0, 2.6)
    ax.set_xticks(xs[::3])
    ax.set_xticklabels([tl[i]["day"][5:] for i in range(0, len(tl), 3)],
                       rotation=45, ha="right")
    ax.set_ylabel("windows per day")
    ax.legend(frameon=True, fontsize=6.2, loc="upper right", framealpha=0.92,
              borderpad=0.25, labelspacing=0.25)
    panel(ax, "d", f"{n_dup} of {n_win} windows are duplicates")

    fig.suptitle("(YWF) why the coverage contrast is not a damage signal",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return save_fig(fig, "fig_ywf_ocv_paper_E", out_dir)


def fig_paper(res, out_dir=None):
    out = []
    for fn in (fig_A, fig_B, fig_C, fig_D, fig_E):
        out.append(fn(res, out_dir))
    return out



# ---------------------------------------------------------------------------
# Pins
# ---------------------------------------------------------------------------
# Exact doubles: a hand copy of a number can be off in the last bits, so the
# tolerance covers the values that come out of numpy/scipy or a bootstrap loop.
# A deviation above the tolerance aborts the script.
FLOAT_TOL = 1e-12
TOL_KEYS = {"p_two_sided", "echo_rate", "delta", "delta_ci95", "sds", "median",
            "mean", "std", "p5", "p25", "p75", "p95", "welch_t", "welch_p",
            "mannwhitney_u", "mannwhitney_p"}


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


def cache_cross_check(res):
    """The committed mask cache against the committed channel table."""
    by_mid = {r["id"]: r for r in res["rows"]}
    n_rows = n_seen = n_ok = 0
    for d in core.load_cache():
        n_rows += 1
        r = by_mid.get(d["mid"])
        if r is None:
            continue
        n_seen += 1
        ok = (int(d["n_masked"]) == int(r["A"])
              and int(d["n_components"]) == int(r["F"]))
        n_ok += 1 if ok else 0
    return {"n_cache_rows": n_rows, "n_chips_with_echo": n_seen,
            "n_matched": n_ok}


def pin_all(res, ref):
    """Every committed number of this package is re-derived and checked."""
    meta, man = res["meta"], res["manifest"]
    P = Pins()

    # 1) the channel table's own metadata against the recomputed layer --------
    P.cmp("meta.site", meta["site"], core.SITE)
    P.cmp("meta.generator", meta["generator"], "code/ywf/ywf_ocv_channels_csv.py")
    P.cmp("meta.out_csv", meta["out_csv"], "ywf_ocv_channels.csv")
    P.cmp("meta.n_rows", meta["n_rows"], res["n_rows"])
    P.cmp("meta.n_columns", meta["n_columns"], len(core.CSV_COLUMNS))
    P.cmp("meta.event", meta["event"], core.COLLAPSE_DATE)
    P.cmp("meta.ocv", meta["ocv"], list(core.OCV))
    P.cmp("meta.features", meta["features"], list(core.FEATURES))
    P.cmp("meta.mask_columns", meta["mask_columns"], list(core.MASK_COLUMNS))
    P.cmp("meta.mask_window", meta["mask_window"], "7x7")
    P.cmp("meta.mask_rule", meta["mask_rule"], man["mask_rule"])
    P.cmp("meta.echo_mode_cuts", meta["echo_mode_cuts"], core.ECHO_MODE_CUTS)
    P.cmp("meta.summary_by_state", echo_mask_by_state(res["rows"]),
          meta["summary_by_state"])
    P.cmp("meta.coverage", res["coverage"], meta["coverage"])
    P.cmp("meta.states", meta["states"], list(core.STATES))
    P.cmp("meta.sections", meta["sections"], list(core.SECTIONS))
    P.cmp("meta.sections_counts",
          {k: v["n"] for k, v in res["composition"]["sections"].items()},
          meta["sections_counts"])
    P.cmp("meta.date_range",
          [res["timeline"][0]["day"], res["timeline"][-1]["day"]],
          meta["date_range"])
    P.cmp("meta.n_windows", meta["n_windows"], sum(e["n_windows"]
                                                   for e in res["timeline"]))
    P.cmp("meta.n_unique_chips", meta["n_unique_chips"], len(res["unique"]))
    P.cmp("meta.n_duplicate_windows", meta["n_duplicate_windows"],
          man["n_deduplicated_same_payload"])
    P.cmp("meta.n_echo_masks", meta["n_echo_masks"], len(res["echo_rows"]))
    P.cmp("meta.n_rows_without_echo_mask", meta["n_rows_without_echo_mask"],
          res["n_rows"] - len(res["echo_rows"]))
    P.cmp("meta.n_pre_collapse", meta["n_pre_collapse"],
          res["coverage"]["per_state"]["pre"]["n_chips"])
    P.cmp("meta.n_post_collapse", meta["n_post_collapse"],
          res["coverage"]["per_state"]["post"]["n_chips"])
    P.cmp("meta.n_P_set", meta["n_P_set"],
          len(_vals(res["rows"], "P")))
    for key, path in (("windows", core.WINDOWS_PATH),
                      ("mask_cache", core.CACHE_PATH),
                      ("mask_cache_manifest", core.MANIFEST_PATH),
                      ("measurements", core.MEAS_PATH),
                      ("segments", core.SEGMENTS_PATH)):
        P.cmp(f"meta.inputs.{key}.sha256", meta["inputs"][key]["sha256"],
              core.sha256(path))
    P.cmp("manifest.source_windows_sha256", man["source_windows_sha256"],
          core.sha256(core.WINDOWS_PATH))

    # 2) the committed mask cache against the channel table -------------------
    cc = cache_cross_check(res)
    P.cmp("cache_cross_check.n_cache_rows", cc["n_cache_rows"],
          meta["cache_cross_check"]["n_cache_rows"])
    P.cmp("cache_cross_check.n_chips_with_echo", cc["n_chips_with_echo"],
          meta["cache_cross_check"]["n_chips_with_echo"])
    P.cmp("cache_cross_check.n_matched", cc["n_matched"], len(res["echo_rows"]))

    # 3) the committed reference file, field by field ------------------------
    P.cmp("findings", json_safe(findings(res)), ref)

    # 4) this package's own numbers against each other ------------------------
    P.cmp("timeline.n_duplicated_sum",
          sum(e["n_duplicated"] for e in res["timeline"]),
          man["n_deduplicated_same_payload"])
    P.cmp("payload_groups.sum_n_windows",
          sum(g["n_windows"] for g in res["payload_groups"]),
          sum(e["n_windows"] for e in res["timeline"]))
    P.cmp("payload_groups.pairs",
          len([g for g in res["payload_groups"] if g["n_windows"] == 2]),
          man["n_deduplicated_same_payload"])
    P.cmp("payload_groups.segments_of_pairs",
          sorted({tuple(g["segments"]) for g in res["payload_groups"]
                  if g["n_windows"] == 2}),
          [(2, 3)])
    for k in core.FEATURES:
        b = res["channels"].get(k) or {}
        P.cmp(f"channels[{k}].n_pre", b.get("n", [None])[0],
              len(_vals([r for r in res["rows"] if r["state"] == core.PRE], k)))
        P.cmp(f"channels[{k}].n_post", b.get("n", [None, None])[1],
              len(_vals([r for r in res["rows"] if r["state"] == core.POST], k)))
    P.cmp("power.n_post_with_echo",
          res["power"]["n_post_with_echo"],
          res["coverage"]["fisher"]["table"][1][0])
    P.cmp("power.echo_rate_pre",
          res["power"]["echo_rate_pre"],
          res["coverage"]["fisher"]["table"][0][0] /
          res["coverage"]["per_state"]["pre"]["n_chips"])
    P.cmp("power.post_chips_needed",
          res["power"]["post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha"],
          post_chips_needed(res["power"]["n_pre_chips"],
                            res["power"]["n_pre_with_echo"]))
    P.cmp("coverage.table_n_echo_sum",
          sum(row[0] for row in res["coverage"]["fisher"]["table"]),
          len(res["echo_rows"]))
    P.cmp("coverage.table_n_chip_sum",
          sum(sum(row) for row in res["coverage"]["fisher"]["table"]),
          len(res["unique"]))
    return P


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def fmt(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}" if math.isfinite(x) else "n/a"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for r in rows:
        out.append("| " + " | ".join(fmt(c) for c in r) + " |")
    return out


def report(res, pins, wrote_reference):
    man, meta = res["manifest"], res["meta"]
    cov = res["coverage"]["per_state"]
    fis = res["coverage"]["fisher"]
    pw, ch, ct = res["power"], res["channels"], res["controls_delta"]
    sec_labels = [lab for lab, b in res["composition"]["sections"].items()
                  if b["n"]]
    sec_txt = ", ".join(f"`{lab}`" for lab in sec_labels) or "none"
    md = []
    md.append(f"# {core.SITE} — observability core vector, pre/post collapse")
    md.append("")
    md.append(f"**Site** {core.SITE} · **asset** {core.ASSET_NAME} · "
              f"**event** {core.COLLAPSE_DATE} · "
              f"**rows** {res['n_rows']} unique 7x7 chips "
              f"({meta['n_windows']} committed windows on "
              f"{len(res['timeline'])} dates) · "
              f"**states** `{core.STATE_SHORT[core.PRE]}` / "
              f"`{core.STATE_SHORT[core.POST]}`")
    md.append("")
    md.append("```")
    md.append("OCV  = [gamma2, P, D]              observability core vector")
    md.append("x    = [gamma2, P, D, A, F, S]      echo-mask vector")
    md.append("mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2")
    md.append("layer: echo mask (whole 7x7 window) — the record is not "
              "segment-resolved, so there is no per-mast-section layer")
    md.append("```")
    md.append("")
    md.append("## What is recomputed, and what is pinned")
    md.append("")
    md.append("| input | role |")
    md.append("| --- | --- |")
    md.append(f"| `data/ywf/{os.path.basename(core.WINDOWS_PATH)}` | the committed "
              f"window payloads (38 kB, committed in full) |")
    md.append(f"| `data/ywf/{os.path.basename(core.CACHE_PATH)}` | the echo mask of "
              f"every committed chip — the mask layer's committed input |")
    md.append(f"| `data/ywf/{os.path.basename(core.SEGMENTS_PATH)}` | the "
              f"mast-section labels + FEM baselines |")
    md.append(f"| `data/ywf/ywf_ocv_channels.csv` | generated: one row per unique "
              f"chip, {meta['n_columns']} columns |")
    md.append(f"| `data/ywf/reference/{os.path.basename(core.REF_FINDINGS)}` | "
              f"the committed analysis-layer reference this script defines and "
              f"pins |")
    md.append("")
    md.append(f"The record holds **{meta['n_windows']}** committed windows on "
              f"**{len(res['timeline'])}** acquisition dates, **{meta['n_rows']}** of "
              f"them unique chips after dropping "
              f"**{man['n_deduplicated_same_payload']}** byte-identical duplicates; "
              f"`segment_resolved = {man['segment_resolved']}`. All "
              f"{meta['n_rows']} unique chips are the same mast section "
              f"({sec_txt}): the 24 duplicated "
              f"windows differ only in the requested section index (2 vs. 3) and "
              f"carry the identical payload, so **one** mask layer is all the "
              f"record supports (fig. E d).")
    md.append("")
    md.append(f"`A` is *not* the pipeline's own `coherence_masked_pixels`: that "
              f"column was computed on the asset-point chip of the whole-scene "
              f"window while this package recomputes the echo mask of the "
              f"committed 7x7 window — the two are equal on "
              f"**{meta['n_coherence_masked_pixels_matched']}** of "
              f"**{meta['n_coherence_masked_pixels_compared']}** rows compared "
              f"(fig. E a).")
    md.append("")
    md.append("## A — the headline: echo coverage")
    md.append("")
    md.append("The 6D vector exists only where the mask rule holds, and the "
              "collapse contrast this site can carry is the *coverage* of that "
              "echo:")
    md.append("")
    body = []
    for s in core.STATES:
        b = cov[core.STATE_SHORT[s]]
        body.append([core.STATE_SHORT[s], b["n_chips"], b["n_with_echo"],
                     b["n_without_echo"],
                     f"{(b['echo_rate'] or 0.0) * 100:.1f} %"])
    md += md_table(["state", "unique chips", "with echo", "without echo",
                    "echo rate"], body)
    md.append("")
    need = pw["post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha"]
    md.append(f"Fisher exact on the 2x2 table "
              f"`{fis['table_rows']} x {fis['table_cols']}` = {fis['table']}: "
              f"**p = {fis['p_two_sided']:.4f}** — not significant, and it cannot "
              f"be: the post side holds {pw['n_post_chips']} chips of which "
              f"{pw['n_post_with_echo']} still carries an echo. Had all "
              f"{pw['n_post_chips']} been echo-free, the same test would need only "
              f"**{need}** post-event chips at alpha = {pw['fisher_alpha']} — so "
              f"this record's coverage contrast is limited by the *surviving "
              f"echo*, not only by the length of the record.")
    md.append("")

    md.append("## B — the channels (where the mask holds)")
    md.append("")
    body = []
    for k in core.FEATURES:
        b = ch.get(k, {})
        if "delta" not in b:
            continue
        body.append([k, f"{b['n'][0]}/{b['n'][1]}", fmt(b["median"][0], 4),
                     fmt(b["median"][1], 4), fmt(b["delta"]), fmt(b["sds"]),
                     "yes" if b.get("ci_excludes_zero") else
                     ("no CI" if "delta_ci95" not in b else "no")])
    md += md_table(["channel", "n pre/post", "pre median", "post median", "delta",
                    "SDS", "CI95 excludes 0"], body)
    md.append("")
    md.append(f"The post column is **one** chip "
              f"(`{pw['n_post_with_echo']}` of the {pw['n_post_chips']} post-event "
              f"chips carry an echo), so the bootstrap CI is undefined there "
              f"(`st.delta_block` declines to resample a single value) and `SDS` "
              f"documents only the spread of the pre side. Fig. A draws those "
              f"distributions, fig. C their effect sizes beside the co-variate "
              f"controls.")
    md.append("")
    md.append("## C — the OCV plane")
    md.append("")
    md.append("Fig. D plots `gamma2` against `P` and `D` for the 17 echo-bearing "
              "chips. No classifier is fitted and no leave-one-out accuracy is "
              "reported: one post-event echo chip is not a class, and an accuracy "
              "computed from it would be an artefact of the sample size rather "
              "than a property of the site.")
    md.append("")
    md.append("## D — what else moves with the collapse date?")
    md.append("")
    md.append("The pre side is Aug-Jan and the post side Feb-Apr, so season, "
              "weather and acquisition generation move with the state. The "
              "co-variate controls beside the mask channels (fig. C b, fig. E):")
    md.append("")
    body = []
    for k in core.CONTROLS:
        b = ct.get(k, {})
        if "delta" not in b:
            continue
        ci = b.get("delta_ci95") or [None, None]
        body.append([k, f"{b['n'][0]}/{b['n'][1]}", fmt(b["median"][0], 3),
                     fmt(b["median"][1], 3), fmt(b["delta"]),
                     f"[{fmt(ci[0])}, {fmt(ci[1])}]"])
    md += md_table(["control", "n pre/post", "pre median", "post median", "delta",
                    "CI95"], body)
    md.append("")
    miss = []
    for k in core.CONTROLS:
        if "delta" in (ct.get(k) or {}):
            continue
        pre_n = (res["controls_by_state"]["pre"].get(k) or {}).get("n", 0)
        post_n = (res["controls_by_state"]["post"].get(k) or {}).get("n", 0)
        miss.append(f"`{k}` ({pre_n}/{post_n})")
    if miss:
        md.append("Carried in the CSV but not contrasted here, because one side "
                  "holds no value: " + ", ".join(miss) + ".")
        md.append("")
    em = res["echo_modes"]
    md.append(f"Echo-mode mix (on this site's own scale, "
              f"`{core.ECHO_MODE_CUTS['compact']}`, "
              f"`{core.ECHO_MODE_CUTS['intermediate']}`, "
              f"`{core.ECHO_MODE_CUTS['distributed']}` of the 49-px window): "
              f"pre {em['pre']['counts']} over {em['pre']['n_echo_chips']} echo "
              f"chips, post {em['post']['counts']} over "
              f"{em['post']['n_echo_chips']}.")
    md.append("")

    md.append("## Pins")
    md.append("")
    if wrote_reference:
        md.append(f"`--write-reference` wrote "
                  f"`data/ywf/reference/{os.path.basename(core.REF_FINDINGS)}` "
                  f"(reloaded and re-compared bit for bit). This site has no "
                  f"upstream analysis project, so **this package defines** the "
                  f"reference it later pins against.")
    else:
        md.append(f"Every number above was re-derived from the committed tables "
                  f"and compared field by field against "
                  f"`data/ywf/reference/{os.path.basename(core.REF_FINDINGS)}`.")
    md.append("")
    ps = pins.summary()
    md.append(f"{ps['n_checks']} checks ({ps['n_exact']} exact, "
              f"{ps['n_within_tolerance']} within {ps['float_tolerance']:g}), "
              f"**{ps['n_failed']} deviations**.")
    md.append("")
    md.append("| committed layer | what is pinned |")
    md.append("| --- | --- |")
    md.append(f"| `data/ywf/reference/{os.path.basename(core.REF_FINDINGS)}` | "
              f"the echo-coverage table, its Fisher test, the mask dimensions per "
              f"state, the pipeline column, the co-variate controls, the timeline "
              f"and the record's provenance — field by field |")
    md.append("| `data/ywf/ywf_ocv_channels_meta.json` | the table's own "
              "definition and summary blocks against this script's recomputation |")
    md.append(f"| `data/ywf/{os.path.basename(core.CACHE_PATH)}` | each cached "
              f"chip's `n_masked` / `n_components` against the CSV's `A` / `F` "
              f"({meta['cache_cross_check']['n_cache_rows']} rows) |")
    md.append("| the committed inputs | their sha256 digests "
              f"(`{os.path.basename(core.WINDOWS_PATH)}`, cache, manifest, "
              f"measurements, segments) |")
    md.append("")
    md.append("### Reading (not a claim of damage)")
    md.append("")
    md.append(f"Before the collapse **{cov['pre']['n_with_echo']} of "
              f"{cov['pre']['n_chips']}** unique chips carry an echo; after it "
              f"**{cov['post']['n_with_echo']} of {cov['post']['n_chips']}**. "
              f"With {pw['n_post_chips']} post-event acquisitions the site cannot "
              f"reach significance, and the single post-event echo chip lies "
              f"inside the pre-event cloud of `gamma2` / `P` / `D` (fig. D). What "
              f"the record supports is the *observability* statement: the echo "
              f"mask is the quantity that can be observed on this collapsed "
              f"tower, and its coverage is the number to extend before any state "
              f"contrast is attempted.")
    md.append("")
    md.append("```bash")
    md.append("python3 code/ywf/ywf_ocv_masks.py --check                  # the mask rule")
    md.append("python3 code/ywf/ywf_ocv_channels_csv.py                 # build + verify")
    md.append("python3 code/ywf/fig_ywf_ocv_paper.py                    # figures + pins")
    md.append("python3 code/ywf/fig_ywf_ocv_paper.py --write-reference  # re-define")
    md.append("```")
    return "\n".join(md) + "\n"


# ---------------------------------------------------------------------------
# JSON artifact, the reference writer and the entry point
# ---------------------------------------------------------------------------
def json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


def build_json(res, pins_summary, wrote_reference):
    """The result JSON — every recomputed number, without the per-row rows."""
    return {
        "site": core.SITE,
        "site_short": core.SITE_SHORT,
        "asset_id": core.ASSET_ID,
        "asset_name": core.ASSET_NAME,
        "event": core.COLLAPSE_DATE,
        "state_labels": list(core.STATES),
        "ocv": list(core.OCV),
        "features": list(core.FEATURES),
        "mask_window": "7x7",
        "mask_rule": res["manifest"]["mask_rule"],
        "echo_mode_cuts": core.ECHO_MODE_CUTS,
        "rng_seed": core.RNG_SEED,
        "n_boot": core.N_BOOT,
        "generated_by": "code/ywf/fig_ywf_ocv_paper.py",
        "wrote_reference": wrote_reference,
        "reference_file": "data/ywf/reference/ywf_ocv_findings.json",
        "n_rows": res["n_rows"],
        "n_windows": sum(e["n_windows"] for e in res["timeline"]),
        "n_unique_chips": len(res["unique"]),
        "n_echo_chips": len(res["echo_rows"]),
        "meta": res["meta"],
        "coverage": res["coverage"],
        "power": res["power"],
        "echo_modes": res["echo_modes"],
        "composition": res["composition"],
        "channels": res["channels"],
        "controls_by_state": res["controls_by_state"],
        "controls_delta": res["controls_delta"],
        "summary_by_state": echo_mask_by_state(res["rows"]),
        "timeline": res["timeline"],
        "payload_groups": res["payload_groups"],
        "pins": pins_summary,
        # no wall-clock time: with the fixed seeds two runs must produce the
        # identical JSON, so the runtime stays on stdout only.
        "figures": [f"figures/ywf/fig_ywf_ocv_paper_{L}.png" for L in "ABCDE"],
    }


def write_reference(res, path):
    """Write the reference JSON and prove the round trip.

    This site has no upstream analysis project, so the committed reference is
    *defined* here once (``--write-reference``) and re-derived by every later run
    of :func:`pin_all`. The file carries no timestamp and is written with sorted
    keys, so two runs are byte-identical.
    """
    payload = json_safe(findings(res))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True, default=json_default)
        fh.write("\n")
    with open(path) as fh:
        back = json.load(fh)
    pins = Pins()
    pins.cmp("reference.round_trip", payload, back)
    return pins


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    ap.add_argument("--reference", default=REF_FINDINGS)
    ap.add_argument("--write-reference", action="store_true",
                    help="write the reference JSON instead of pinning against it")
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, event {core.COLLAPSE_DATE} …")
    res = compute()
    print(f"  {res['n_rows']} unique chips "
          f"({res['coverage']['per_state']['pre']['n_with_echo']} pre / "
          f"{res['coverage']['per_state']['post']['n_with_echo']} post carry an "
          f"echo), {res['runtime_s']:.1f} s")

    if args.write_reference:
        pins = write_reference(res, args.reference)
        print(f"  reference written: {os.path.relpath(args.reference)}")
    else:
        if not os.path.isfile(args.reference):
            raise SystemExit(f"missing reference file {args.reference}: run this "
                             f"script once with --write-reference")
        pins = pin_all(res, core.load_findings(args.reference))
    summary = pins.summary()
    print(f"pins: {summary['n_checks']} checks, {summary['n_exact']} exact, "
          f"{summary['n_within_tolerance']} within "
          f"{summary['float_tolerance']:g}, {summary['n_failed']} failed")

    os.makedirs(args.figdir, exist_ok=True)
    fig_paper(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(json_safe(build_json(res, summary, args.write_reference)), fh,
                  indent=1, sort_keys=True, default=json_default)
    print(f"  {os.path.relpath(args.json)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.md)), exist_ok=True)
    with open(args.md, "w") as fh:
        fh.write(report(res, pins, args.write_reference))
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

