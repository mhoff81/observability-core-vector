#!/usr/bin/env python3
"""Morandi OCV-Paper — figures A–E, JSON result and English report.

Question: does the observability core vector ``OCV = [gamma2, P, D]`` of the
Ponte Morandi / Polcevera viaduct (Genova) — the interferometric coherence
gamma^2 of its deck echo, plus the persistence ``P`` and the density ``D`` of
that echo mask — separate the two states of the site, the viaduct **before** the
2018-08-14 collapse (``pre-collapse (healthy)``) and **after** it
(``post-collapse``)? And does the *per-track deck mask* — the same six
dimensions on the ``A_asc`` / ``A_des`` tracks only, the layer this package adds
— carry the contrast, or does the whole-window chip alone?

The full mask vector is the 6D ``x = [gamma2, P, D, A, F, S]`` (``A`` = masked
pixels, ``F`` = 8-connected components, ``S`` = ||mask centroid - peak||). Both
mask layers are recomputed from the committed window cache
``data/morandi/morandi_windows_mask_cache.txt`` with the mask rule of the
analysis project (``0.30 x peak`` and ``5 x np.median``, ``min_n_masked = 2``);
see ``code/morandi/morandi_ocv_masks.py``.

What this package claims, and what it does not
----------------------------------------------
  * **two states on one bridge.** ``pre-collapse (healthy)`` and
    ``post-collapse`` are the two sides of the collapse date of the *same*
    structure, so season, weather and pipeline generation all move with the
    state. Every contrast is therefore reported beside a co-variate control
    (weather, orbit, track) — never as damage.
  * **the post side is under-powered.** Only **11** chips survive the event
    (5 ``A_asc`` + 6 ``A_des``) against **221** before it, so every post-state
    number carries a wide interval. This figure is an **honest null**: it shows
    that the two states are not separable here, not that the bridge is intact.
  * **the analysed window is not on the collapsed span** and the archive holds
    no persistent scatterer at this footprint (see
    ``morandi_analysis/MORANDI_SCATTERER_LIMITS.md``) — another reason not to
    read any contrast as a damage signal.
  * **two mask layers.** The echo mask is the whole 80x80 chip; the deck mask is
    the track's own chip. ``A`` is *not* the site's ``coherence_masked_pixels``:
    that column is a fixed 640-px (10 %) quantile mask
    (``morandi_deck_channels.csv``), while this package recomputes the echo mask
    from the payloads — the two never coincide, which is reported, not hidden.

Pinned reference files (``data/morandi/reference/``):
  morandi_registration.json       the committed co-registration/QC table: per
                                  track n_dates / n_pre / n_post, first / last /
                                  master day, the fixed 640-px deck mask and the
                                  per-date whole-window metrics — reproduced
  morandi_deck_channels.json      the site's fixed quantile deck-mask table,
                                  pinned structurally
  morandi_geometry.json           the viaduct geometry (assumption grid; a
                                  geometry extract, not rederivable)

Figures:
  A  raw distributions of the six channels per state (all 232 chips)
  B  in-sample effect sizes: Cliff's delta + bootstrap CI and SDS for the six
     channels, pooled and per track
  C  2-class LDA (pre vs. post) leave-one-out accuracy + permutation null
  D  the OCV plane (gamma2, P, D): does the core vector separate the states?
  E  the evidence against reading the contrast as damage: weather controls,
     orbit/track composition, echo-mode mix, and the export-vs-recomputed mask

Usage:
  python3 code/morandi/fig_morandi_ocv_paper.py            # recompute everything
  python3 code/morandi/fig_morandi_ocv_paper.py --quick    # no permutation null
"""
from __future__ import annotations

import argparse
import hashlib
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

import morandi_ocv_core as core    # noqa: E402
import morandi_ocv_masks as masks  # noqa: E402
import morandi_ocv_stats as st     # noqa: E402

FIGDIR = core.FIGDIR
OUT_JSON = os.path.join(core.DATA, "fig_morandi_ocv_paper.json")
OUT_MD = os.path.join(FIGDIR, "fig_morandi_ocv_paper.md")

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
GIRDER_COLORS = {"A_asc": "#4C72B0", "A_des": "#C44E52"}
CHANNEL_COLORS = {"gamma2": "#4C72B0", "P": "#55A868", "D": "#8172B3",
                  "A": "#937860", "F": "#64B5CD", "S": "#DA8BC3",
                  "wind_speed_ms": "#64B5CD", "temperature_c": "#C44E52",
                  "precipitation_mm": "#4C72B0", "wind_gust_ms": "#8172B3",
                  "nearest_lag_h": "#8C8C8C"}
CHANNELS = list(core.FEATURES) + list(core.WEATHER) + ["nearest_lag_h"]


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


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Row views
# ---------------------------------------------------------------------------
def mask_rows(csv_rows):
    """The rows as this package analyses them (full-precision gamma2)."""
    out = []
    for r in csv_rows:
        d = core.row_brief(r)
        d["month"] = r.get("month")
        d["day"] = core.day_of(r)
        d["echo_mode"] = r.get("echo_mode")
        d["A_export"] = core._num(r.get("coherence_masked_pixels"))
        out.append(d)
    return out


def _vals(rows, key):
    return [r[key] for r in rows if r.get(key) is not None]


def delta_pair(rows, key, ref=core.PRE, grp=core.POST,
               n_boot=core.N_BOOT, seed=core.RNG_SEED):
    """``st.delta_block`` for one channel between the two states."""
    a = _vals([r for r in rows if r.get("state") == ref], key)
    b = _vals([r for r in rows if r.get("state") == grp], key)
    return st.delta_block(a, b, f"{ref} -> {grp}", n_boot=n_boot, seed=seed)


def lda_block(rows, features, shrinkage, n_perm, seed=core.RNG_SEED):
    """2-class LDA LOO accuracy + label-permutation null over ``features``."""
    X, y = [], []
    for r in rows:
        if all(r.get(f) is not None for f in features):
            X.append([r[f] for f in features])
            y.append(0 if r["state"] == core.PRE else 1)
    X, y = np.asarray(X, dtype=float), np.asarray(y, dtype=int)
    classes = [0, 1]
    _, acc, bal, conf = st.loo_cv(X, y, classes, shrinkage)
    out = {"n": int(X.shape[0]), "n_features": len(features),
           "accuracy": acc, "balanced_accuracy": bal, "confusion": conf}
    if n_perm:
        out["permutation"] = st.permutation_test(X, y, classes, shrinkage,
                                                 n_perm=n_perm, seed=seed)
    return out


def oof_block(rows, features, shrinkage, class_a=core.PRE, class_b=core.POST):
    """LOO Fisher scores + SDS between the two states for one model."""
    X, y = [], []
    for r in rows:
        if all(r.get(f) is not None for f in features):
            X.append([r[f] for f in features])
            y.append(r["state"])
    X, y = np.asarray(X, dtype=float), np.asarray(y, dtype=object)
    q = st.oof_scores(X, y, class_a, class_b, shrinkage)
    d, sds, _ref, _grp = st.delta_sds(q, y, class_a, class_b)
    return {"n_features": len(features), "delta": d, "sds": sds}


def compute(quick=False):
    t0 = time.time()
    rows = mask_rows(core.load_csv())
    manifest = core.load_manifest()
    meta = core.load_csv_meta()
    by_state = st.by_state(rows)

    channels = {k: delta_pair(rows, k) for k in CHANNELS}
    tracks = {}
    for g in core.GIRDERS:
        sub = [r for r in rows if r.get("segment") == g]
        tracks[g] = {k: delta_pair(sub, k, n_boot=max(core.N_BOOT // 20, 500))
                     for k in core.FEATURES}
    models = {}
    n_perm = 0 if quick else core.N_PERM
    for name, feats in core.MODELS.items():
        models[name] = lda_block(rows, feats, core.HEADLINE_LAMBDA, n_perm)
    models["ocv_oof"] = oof_block(rows, core.OCV, core.HEADLINE_LAMBDA)
    models["x6d_oof"] = oof_block(rows, core.FEATURES, core.HEADLINE_LAMBDA)

    return {
        "rows": rows,
        "states": {core.STATE_SHORT[s]: len(v) for s, v in by_state.items()},
        "channels": channels,
        "tracks": tracks,
        "models": models,
        "meta": meta,
        "manifest": manifest,
        "n_rows": len(rows),
        "quick": quick,
        "runtime_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# The figures A-E
# ---------------------------------------------------------------------------
def _jitter(n, width, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def fig_A(res, out_dir=None):
    """Raw distributions of the six channels per state (all chips)."""
    rows = res["rows"]
    fig, axes = plt.subplots(2, 3, figsize=(COL_WIDTH, 0.62 * COL_WIDTH))
    for i, k in enumerate(core.FEATURES):
        ax = axes[i // 3][i % 3]
        for j, s in enumerate(core.STATES):
            vals = _vals([r for r in rows if r["state"] == s], k)
            if not vals:
                continue
            bp = ax.boxplot([vals], positions=[j], widths=0.5,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="white", lw=1.1))
            bp["boxes"][0].set(facecolor=STATE_COLORS[s], alpha=0.65,
                               edgecolor=STATE_COLORS[s])
            ax.plot(j + _jitter(len(vals), 0.14, 11 + i),
                    vals, ".", ms=1.6, alpha=0.35, color=STATE_COLORS[s],
                    markeredgewidth=0)
            ax.text(j, 0.99, f"n={len(vals)}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=6.4, color=STATE_COLORS[s])
        ax.set_xticks(range(len(core.STATES)))
        ax.set_xticklabels([core.STATE_SHORT[s] for s in core.STATES])
        ax.set_ylabel(core.CHANNEL_TITLE.get(k, k))
        if k in ("A", "D", "S"):
            ax.set_yscale("log")
        panel(ax, "abcdef"[i], core.CHANNEL_AXIS.get(k, k))
    fig.suptitle(f"(Morandi) echo-mask channels per state — {res['n_rows']} chips, "
                 f"pre = {res['states'].get('pre', 0)}, "
                 f"post = {res['states'].get('post', 0)}",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return save_fig(fig, "fig_morandi_ocv_paper_A", out_dir)


def fig_B(res, out_dir=None):
    """In-sample effect sizes: Cliff's delta + CI + SDS, pooled and per track."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_WIDTH, 0.40 * COL_WIDTH),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})
    chans = list(core.FEATURES)

    ax = axes[0]
    ys = np.arange(len(chans))
    for i, k in enumerate(chans):
        b = res["channels"][k]
        if "delta" not in b:
            continue
        d = b["delta"]
        ax.plot([d], [i], "o", ms=5.0, color=CHANNEL_COLORS[k])
        if "delta_ci95" in b:
            lo, hi = b["delta_ci95"]
            ax.plot([lo, hi], [i, i], "-", lw=1.4, color=CHANNEL_COLORS[k])
        excl = b.get("ci_excludes_zero")
        ax.text(1.02, i, "CI≠0" if excl else "CI∋0", fontsize=6.2, va="center",
                color="#444444" if excl else "#AA2222",
                transform=ax.get_yaxis_transform())
    ax.axvline(0.0, color="#888888", lw=0.9)
    ax.set_yticks(ys)
    ax.set_yticklabels(chans)
    ax.invert_yaxis()
    ax.set_xlabel("Cliff's delta  (post − pre)")
    panel(ax, "a", "pooled, all chips (bootstrap 95 % CI)")

    ax = axes[1]
    width = 0.38
    for gi, g in enumerate(core.GIRDERS):
        for i, k in enumerate(chans):
            b = res["tracks"][g][k]
            if "delta" not in b:
                continue
            x = i + (gi - 0.5) * width
            ax.bar(x, b["delta"], width=width, color=GIRDER_COLORS[g],
                   alpha=0.80, label=g if i == 0 else None)
    ax.axhline(0.0, color="#888888", lw=0.9)
    ax.set_xticks(range(len(chans)))
    ax.set_xticklabels(chans)
    ax.set_ylabel("Cliff's delta  (post − pre)")
    ax.legend(frameon=False)
    panel(ax, "b", "per track (deck mask)")

    fig.suptitle("(Morandi) in-sample effect sizes of the six channels",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_morandi_ocv_paper_B", out_dir)


def fig_C(res, out_dir=None):
    """2-class LDA (pre vs. post): LOO accuracy + permutation null."""
    order = ["gamma2", "PD", "gamma2PD", "x_6d", "weather"]
    fig, ax = plt.subplots(figsize=(0.62 * COL_WIDTH, 0.40 * COL_WIDTH))
    xs = np.arange(len(order))
    accs = [res["models"][m]["accuracy"] for m in order]
    ax.bar(xs, accs, width=0.62, color="#4C72B0", alpha=0.80)
    for i, m in enumerate(order):
        p = res["models"][m].get("permutation", {}).get("p_accuracy")
        txt = f"p={p:.3f}" if p is not None else "no null"
        ax.text(i, min(accs[i] + 0.02, 1.0), txt, ha="center", va="bottom",
                fontsize=6.4, color="#444444")
    ax.axhline(0.5, color="#888888", ls="--", lw=0.9)
    labels = ["gamma2", "P+D", "gamma2+P+D", "6 dims", "weather"]
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("LOO accuracy")
    ax.set_ylim(0, 1.0)
    panel(ax, "a", "LDA pre vs. post (chance = 0.5)")
    fig.suptitle("(Morandi) is the state learnable from the chips?", fontsize=9.0,
                 x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_fig(fig, "fig_morandi_ocv_paper_C", out_dir)


def fig_D(res, out_dir=None):
    """The OCV plane (gamma2, P, D): does the core vector separate the states?"""
    rows = res["rows"]
    pairs = [("gamma2", "P"), ("gamma2", "D"), ("P", "D")]
    fig, axes = plt.subplots(1, 3, figsize=(COL_WIDTH, 0.34 * COL_WIDTH))
    for ax, (kx, ky) in zip(axes, pairs):
        for s in core.STATES:
            sel = [r for r in rows if r["state"] == s
                   and r.get(kx) is not None and r.get(ky) is not None]
            ax.plot([r[kx] for r in sel], [r[ky] for r in sel], "o", ms=2.4,
                    alpha=0.45, color=STATE_COLORS[s], markeredgewidth=0,
                    label=core.STATE_SHORT[s])
        ax.set_xlabel(core.CHANNEL_TITLE.get(kx, kx))
        ax.set_ylabel(core.CHANNEL_TITLE.get(ky, ky))
    axes[1].set_yscale("log")
    axes[2].set_yscale("log")
    axes[0].legend(frameon=False, markerscale=2)
    panel(axes[0], "a", "gamma2 vs P")
    panel(axes[1], "b", "gamma2 vs D")
    panel(axes[2], "c", "P vs D")
    fig.suptitle("(Morandi) the OCV plane per state — no clean partition",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return save_fig(fig, "fig_morandi_ocv_paper_D", out_dir)


def fig_E(res, out_dir=None):
    """Evidence against reading the contrast as damage: co-variate controls."""
    rows = res["rows"]
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, 0.60 * COL_WIDTH))
    width = 0.38

    # (a) weather-control effect sizes
    ax = axes[0][0]
    wc = list(core.WEATHER)
    for i, k in enumerate(wc):
        b = res["channels"][k]
        if "delta" not in b:
            continue
        ax.plot([b["delta"]], [i], "o", ms=5.0, color=CHANNEL_COLORS[k])
        if "delta_ci95" in b:
            lo, hi = b["delta_ci95"]
            ax.plot([lo, hi], [i, i], "-", lw=1.4, color=CHANNEL_COLORS[k])
    ax.axvline(0.0, color="#888888", lw=0.9)
    ax.set_yticks(range(len(wc)))
    ax.set_yticklabels([core.CHANNEL_TITLE.get(k, k) for k in wc])
    ax.invert_yaxis()
    ax.set_xlabel("Cliff's delta  (post − pre)")
    panel(ax, "a", "weather moves with the state too")

    # (b) composition: track x state counts
    ax = axes[0][1]
    tracks = list(core.GIRDERS)
    for si, s in enumerate(core.STATES):
        counts = [len([r for r in rows if r["state"] == s
                       and r.get("segment") == g]) for g in tracks]
        xs = np.arange(len(tracks)) + (si - 0.5) * width
        ax.bar(xs, counts, width=width, color=STATE_COLORS[s], alpha=0.80,
               label=core.STATE_SHORT[s])
    ax.set_xticks(range(len(tracks)))
    ax.set_xticklabels(tracks)
    ax.set_ylabel("chips")
    ax.legend(frameon=False)
    panel(ax, "b", "track composition per state")

    # (c) echo-mode mix per state
    ax = axes[1][0]
    modes = ["compact", "intermediate", "distributed"]
    for si, s in enumerate(core.STATES):
        vals = _vals([r for r in rows if r["state"] == s], "A")
        counts = [sum(1 for a in vals if core.echo_mode_of(a) == m) for m in modes]
        tot = sum(counts) or 1
        ax.bar([i + (si - 0.5) * width for i in range(len(modes))],
               [c / tot for c in counts], width=width,
               color=STATE_COLORS[s], alpha=0.80, label=core.STATE_SHORT[s])
    ax.set_xticks(range(len(modes)))
    ax.set_xticklabels(modes)
    ax.set_ylabel("fraction of chips")
    ax.legend(frameon=False)
    panel(ax, "c", "echo-mode mix (A-pixel count)")

    # (d) A_export (site's fixed 640-px mask) vs the recomputed echo A
    ax = axes[1][1]
    for s in core.STATES:
        sel = [r for r in rows if r["state"] == s and r.get("A_export") is not None]
        ax.plot([r["A_export"] for r in sel], [r["A"] for r in sel], "o", ms=2.4,
                alpha=0.45, color=STATE_COLORS[s], markeredgewidth=0,
                label=core.STATE_SHORT[s])
    ax.plot([1, 400], [1, 400], "-", color="#888888", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("A export = coherence_masked_pixels (fixed 640 px)")
    ax.set_ylabel("A recomputed (echo mask)")
    panel(ax, "d", "export vs recomputed mask")

    fig.suptitle("(Morandi) why this contrast is not a damage signal",
                 fontsize=9.0, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return save_fig(fig, "fig_morandi_ocv_paper_E", out_dir)


def fig_paper(res, out_dir=None):
    out = []
    for fn in (fig_A, fig_B, fig_C, fig_D, fig_E):
        out.append(fn(res, out_dir))
    return out


# ---------------------------------------------------------------------------
# The pin block: every committed number is re-derived and checked; a deviation
# above the tolerance makes the script exit non-zero.
# ---------------------------------------------------------------------------
def pin_block(res):
    problems = []

    def check(name, got, want, tol=0.0):
        if isinstance(want, float) or isinstance(got, float):
            ok = got is not None and want is not None and abs(got - want) <= tol
        else:
            ok = got == want
        if not ok:
            problems.append(f"{name}: got {got!r} != want {want!r}")
        return ok

    meta = res["meta"]
    man = res["manifest"]

    # --- the generated channel table --------------------------------
    check("meta.site", meta["site"], core.SITE)
    check("meta.event", meta["event"], core.COLLAPSE_DATE)
    check("meta.generator", meta["generator"],
          "code/morandi/morandi_ocv_channels_csv.py")
    check("meta.n_rows", meta["n_rows"], 232)
    check("meta.n_pre_collapse", meta["n_pre_collapse"], 221)
    check("meta.n_post_collapse", meta["n_post_collapse"], 11)
    check("meta.n_columns", meta["n_columns"], 59)
    check("meta.n_echo_masks", meta["n_echo_masks"], 232)
    check("meta.n_rows_without_echo_mask", meta["n_rows_without_echo_mask"], 0)
    check("meta.states", list(meta["states"]), list(core.STATES))
    check("meta.features", list(meta["features"]), list(core.FEATURES))
    check("meta.ocv", list(meta["ocv"]), list(core.OCV))
    check("meta.girders", list(meta["girders"]), list(core.GIRDERS))
    check("meta.deck_mask_export", meta["deck_mask_export"],
          "fixed 640-px (10 %) quantile mask (morandi_deck_channels.csv / "
          "morandi_registration.json)")
    check("meta.n_coherence_masked_pixels_compared",
          meta["n_coherence_masked_pixels_compared"], 232)
    check("meta.n_coherence_masked_pixels_matched",
          meta["n_coherence_masked_pixels_matched"], 0)

    # --- the recomputed rows against the metadata -------------------
    check("res.n_rows", res["n_rows"], meta["n_rows"])
    check("res.states.pre", res["states"]["pre"], meta["n_pre_collapse"])
    check("res.states.post", res["states"]["post"], meta["n_post_collapse"])

    # --- the committed mask-cache manifest --------------------------
    check("manifest.n_rows", man["n_rows"], 232)
    check("manifest.crop", list(man["crop"]), [160, 240, 160, 240])
    check("manifest.window_size", man["window_size"], 80)
    check("manifest.source_size", man["source_size"], 400)
    check("manifest.window_sizes.400x400", man["window_sizes"]["400x400"], 234)
    check("manifest.n_without_echo", man["n_without_echo"], 2)
    check("manifest.n_clipped_or_odd_shape",
          man["n_clipped_or_odd_shape"], 5)
    rule = man["mask_rule"]
    check("manifest.mask_rule.peak_frac", rule["peak_frac"], 0.30)
    check("manifest.mask_rule.median_mult", rule["median_mult"], 5.0)
    check("manifest.mask_rule.min_n_masked", rule["min_n_masked"], 2)
    # the source payload digests are provenance — pinned so a tampered
    # manifest cannot masquerade as the committed one.
    check("manifest.src.A_asc.sha256",
          man["source_rect_files"]["ref_rects_A_asc_400.txt"]["sha256"],
          "5c6bfc68421d6fa9b9d13da49e46993885b38c987f3654b46d249929c0e80ea5")
    check("manifest.src.A_des.sha256",
          man["source_rect_files"]["ref_rects_A_des_400.txt"]["sha256"],
          "201d8c8987a79604a697c6fdf62badc3e647e06c5d18f8f4eaa6b18e3fffd4a1")

    # --- per-track counts: recomputed == cache == committed table ----
    reg = {t["tag"]: t for t in core.load_reference(core.REF_STATES)["tracks"]}
    for tag in core.GIRDERS:
        got = len([r for r in res["rows"] if r.get("segment") == tag])
        check(f"rows.{tag}.n", got, man["n_per_track"].get(tag))
        t = reg.get(tag)
        if t is not None:
            check(f"registration.{tag}.n_masked_pixels",
                  t.get("n_masked_pixels"), 640)
            check(f"registration.{tag}.n_masked_pixels_all_dates",
                  t.get("n_masked_pixels_all_dates"), 640)
    # A_des is the second track; its source n_dates (115) exceeds the cache
    # count (113) by the two clipped payloads — both facts are pinned.
    if "A_des" in reg:
        check("registration.A_des.n_dates", reg["A_des"].get("n_dates"), 115)
        check("manifest.n_per_track.A_des", man["n_per_track"].get("A_des"), 113)
    if "A_asc" in reg:
        check("registration.A_asc.n_dates", reg["A_asc"].get("n_dates"), 119)
        check("registration.A_asc.n_pre", reg["A_asc"].get("n_pre"), 114)
        check("registration.A_asc.n_post", reg["A_asc"].get("n_post"), 5)

    # --- the site's fixed deck-mask table ---------------------------
    deck = core.load_reference(core.REF_MASK)
    for tag in core.GIRDERS:
        d = deck.get(tag, {})
        check(f"deck.{tag}.n_masked_px", d.get("n_masked_px"), 640)

    # --- the geometry (structural pin) ------------------------------
    geom = core.load_reference(core.REF_OSM)
    check("geometry.system_span_m", geom.get("system_span_m"), 208)
    check("geometry.n_systems", geom.get("n_systems"), 3)
    check("geometry.analysed_field_collapsed",
          geom.get("analysed_field_collapsed"), False)
    return problems


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
    states = res["states"]
    md = []
    md.append("# Ponte Morandi (Polcevera), Genova — the observability core "
              "vector around the collapse")
    md.append("")
    md.append(f"**Site** {core.SITE} · **event** {core.COLLAPSE_DATE} · "
              f"**chips** {res['n_rows']} ({states.get('pre', 0)} pre-collapse "
              f"(healthy) / {states.get('post', 0)} post-collapse)")
    md.append("")
    md.append("```")
    md.append("state = pre-collapse (healthy)  day <  2018-08-14")
    md.append("        post-collapse            day >= 2018-08-14")
    md.append("x     = [gamma2, P, D, A, F, S]   echo-mask vector")
    md.append("OCV   = [gamma2, P, D]            observability core vector")
    md.append("mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), "
              "A >= 2")
    md.append("```")
    md.append("")
    md.append("## Honest null — read this first")
    md.append("")
    md.append(f"Only **{states.get('post', 0)}** chips survive the collapse "
              f"(5 `A_asc` + 6 `A_des`) against **{states.get('pre', 0)}** "
              "before it, so every post-state number carries a wide interval. "
              "The analysed 80x80 window is **not on the collapsed span** and "
              "the archive holds no persistent scatterer at this footprint. "
              "This figure therefore asks whether the states are *separable* "
              "here — it must never be read as a damage verdict.")
    md.append("")
    md.append("## Effect sizes of the six channels (post − pre)")
    md.append("")
    rows = []
    for k in core.FEATURES:
        b = res["channels"][k]
        ci = b.get("delta_ci95")
        rows.append([k, _fmt(b.get("delta")),
                     f"[{_fmt(ci[0])}, {_fmt(ci[1])}]" if ci else "–",
                     "yes" if b.get("ci_excludes_zero") else "no",
                     _fmt(b.get("sds"))])
    md += md_table(["channel", "Cliff's delta", "95 % CI", "CI≠0", "SDS"], rows)
    md.append("")
    md.append("## Per-track deck-mask effect sizes")
    md.append("")
    for g in core.GIRDERS:
        md.append(f"**{g}**")
        md.append("")
        rows = []
        for k in core.FEATURES:
            b = res["tracks"][g][k]
            ci = b.get("delta_ci95")
            rows.append([k, _fmt(b.get("delta")),
                         f"[{_fmt(ci[0])}, {_fmt(ci[1])}]" if ci else "–",
                         "yes" if b.get("ci_excludes_zero") else "no"])
        md += md_table(["channel", "Cliff's delta", "95 % CI", "CI≠0"], rows)
        md.append("")
    md.append("## Is the state learnable? (2-class LDA, leave-one-out)")
    md.append("")
    rows = []
    for name in ["gamma2", "PD", "gamma2PD", "x_6d", "weather"]:
        m = res["models"][name]
        p = m.get("permutation", {}).get("p_accuracy")
        rows.append([name, m["n_features"], _fmt(m["accuracy"], 3),
                     _fmt(m["balanced_accuracy"], 3),
                     f"{p:.3f}" if p is not None else "–"])
    md += md_table(["model", "dims", "LOO acc.", "balanced acc.", "p (perm.)"], rows)
    md.append("")
    md.append("## What this does and does not say")
    md.append("")
    md.append("Both sides of the event live on the *same* bridge, so season, "
              "weather and pipeline generation move with the state — every "
              "post-state interval is wide and no contrast here is a damage "
              "verdict. The point of the figure is the **honest null**: with "
              f"{core.STATE_SHORT[core.POST]} = {states.get('post', 0)} chips the "
              "states are not cleanly separable in the analysed window, and the "
              "`A` column of the site's export is a *fixed 640-px quantile mask*, "
              "not the echo mask this package recomputes — the two never "
              "coincide (panel E).")
    md.append("")
    md.append("```bash")
    md.append("python3 code/morandi/morandi_ocv_channels_csv.py   # build + verify the table")
    md.append("python3 code/morandi/fig_morandi_ocv_paper.py      # figures A–E + this report")
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


def build_json(res, problems):
    """The result JSON — the aggregate blocks and the per-track/channel deltas."""
    return {
        "site": core.SITE,
        "event": core.COLLAPSE_DATE,
        "generated_by": "code/morandi/fig_morandi_ocv_paper.py",
        "states": dict(core.STATE_SHORT),
        "features": list(core.FEATURES),
        "ocv": list(core.OCV),
        "models": list(core.MODELS),
        "n_rows": res["n_rows"],
        "state_counts": res["states"],
        "channels": res["channels"],
        "tracks": res["tracks"],
        "models_fit": res["models"],
        "pins": {"ok": not problems, "n_problems": len(problems),
                 "problems": problems},
        "figures": [f"figures/morandi/fig_morandi_ocv_paper_{L}.png"
                    for L in "ABCDE"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--json", default=OUT_JSON)
    ap.add_argument("--md", default=OUT_MD)
    ap.add_argument("--quick", action="store_true",
                    help="skip the LDA label-permutation null (fast)")
    args = ap.parse_args(argv)

    t0 = time.time()
    print(f"computing — {core.SITE}, event {core.COLLAPSE_DATE} …")
    res = compute(quick=args.quick)
    print(f"  {res['n_rows']} chips ({res['states']}), "
          f"{res['runtime_s']:.1f} s")

    problems = pin_block(res)
    if problems:
        print(f"PIN FAILED ({len(problems)} problem(s)):")
        for p in problems:
            print(f"  - {p}")
        return 2
    print("  pins OK")

    os.makedirs(args.figdir, exist_ok=True)
    fig_paper(res, args.figdir)

    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(build_json(res, problems), fh, indent=1, sort_keys=True,
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
