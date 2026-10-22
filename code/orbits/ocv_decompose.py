#!/usr/bin/env python3
"""Decompose a site's post-collapse signature into four contributions.

One module, one site at a time. The registry below is the only place a site is
described: its table, its event date, how its acquisition geometry becomes
strata, which targets exist, and which model tree its post group can carry.

  1. the observed area ``A``   (explanatory variable, never a target)
  2. the acquisition geometry  (angle strata, channel strata, or pooled)
  3. temperature and season
  4. what is left              (the adjusted post effect)

Two ladders are reported in parallel:

    A-adjusted   M0 -> M1 (+log A) -> M2 (+T) -> M3 (+T, season)
    A-free       M0 -> M2- (+T)   -> M3- (+T, season)

A post effect that attenuates under A adjustment while staying put in the A-free
ladder **co-moves with the observed area**. That is co-movement, not cause: the
area is itself moved by a collapse, so conditioning on it can remove the change
being measured -- and can also create association where the raw data show none.

Model trees, chosen by the post group
-------------------------------------
A model is not fitted when it is not estimable. The rule is ``n_post >= 3 * k``
with k the number of columns of the design; a cell that fails it is reported as
``refused`` with the reason, never as a converged fit on ten points.

    full      k up to 6   A-adjusted M0..M3, A-free M0..M3-, sensitivities
    reduced   k up to 3   A-adjusted M0..M1, A-free M0..M2- (temperature only)
    minimal   k up to 2   A-adjusted M0..M1, A-free M0

The tree a site *intends* is in the registry; the guard still refuses per cell,
so the report can show a site that is thinner than its own label.

Usage
-----
  python3 code/orbits/ocv_decompose.py --site all
  python3 code/orbits/ocv_decompose.py --site cts --boot 400
  python3 code/orbits/ocv_decompose.py --selftest
  python3 code/orbits/ocv_decompose.py --verify --site carola

Outputs
-------
  data/<site>/support/<site>_ocv_decompose.json
  data/<site>/support/<site>_ocv_decompose.md
  data/decompose_summary.md              (--site all)
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
for path in (HERE, os.path.join(REPO, "code", "carola"),
             os.path.join(REPO, "code", "espoo"),
             os.path.join(REPO, "code", "lumo")):
    if path not in sys.path:
        sys.path.insert(0, path)

# ``carola_ocv_stats`` is the house port of the LUMO statistics layer; the other
# site packages carry byte-identical copies of ``ols_fit``/``f_test``, so one
# import serves every site (see code/carola/README.md).
import carola_ocv_stats as stats_mod  # noqa: E402
import espoo_ocv_stats as espoo  # noqa: E402
from carola_ocv_precollapse import holm, sha256_of  # noqa: E402
from lumo_ocv_stats import cliffs_delta  # noqa: E402

SCHEMA = "ocv-decompose-v2"

#: The registry: the only place a site is described.
##``strata`` is how the acquisition geometry becomes strata -- ``angle`` splits
# the ASC channel at ``angle_split`` degrees (a strict 2-cycle in Carola),
###``channel`` uses ASC and DESC as they are, ``single`` pools everything.
SITES = {
    "carola": {
        "label": "Carolabr\u00fccke deck",
        "event": dt.datetime(2024, 9, 11),
        "strata": "angle", "angle_split": 40.0,
        "targets": ("gamma2", "P", "D", "F", "S"),
        "tree": "full",
    },
    "cts": {
        "label": "CTS deck",
        "event": dt.datetime(2021, 6, 24),
        "strata": "channel",
        "targets": ("gamma2", "P", "D", "F", "S"),
        "tree": "full",
    },
    "haltern": {
        "label": "Haltern canal bridge",
        "event": dt.datetime(2021, 9, 29),
        "strata": "single",
        "targets": ("gamma2", "P", "D", "F", "S"),
        "tree": "reduced",
    },
    "morandi": {
        "label": "Ponte Morandi",
        "event": dt.datetime(2018, 8, 14),
        "strata": "single",
        "targets": ("gamma2", "P", "D", "F", "S"),
        "tree": "minimal",
    },
}

#: Live configuration (set by ``configure``). The analysis functions read these,
##so the module is configured for one site at a time rather than threaded through
###every call -- a deliberate trade for keeping the numbers reproducible.
SITE = "carola"
LABEL = ""
TABLE = ""
OUT_JSON = ""
OUT_MD = ""
EVENT = dt.datetime(2024, 9, 11)
STRATA = ()
TARGETS = ("gamma2", "P", "D", "F", "S")
LOG_TARGETS = ("gamma2", "P", "D", "F")
STRATUM_MODE = "angle"
ANGLE_SPLIT = 40.0
TREE = "full"
#: ``gamma2`` is the **mask-derived echo coherence** (``gamma2_mask``), the
# OCV dimension the orbit figures and the SDS work use -- see the ``"gamma2"``
# entry of every site in ``code/orbits/fig_ocv_orbits.py``. It is *not* the
# pipeline's ``coherence_gamma2``, which is copied verbatim from the committed
# export and measured on a different chip: the two are different quantities
# (r = 0.09 on Carola, medians 0.136 vs 0.003). Morandi has no
# ``coherence_gamma2`` at all, only the mask.
#
# The chips are not 80 x 80 and differ per site, because ``tools/ocv_table.py``
# cuts them from the support geometry (guard 30 m, reference ring 30-200 m):
# Carola 69 x 125, CTS 39 x 145, Haltern 35 x 229, Morandi 51 x 137. That
# builder recomputes ``A, D, F, S, P`` and ``gamma2`` on that chip and keeps
# the committed values as ``*_committed``; its docstring says the two are
# "not comparable with the committed table row by row". A support-chip model
# therefore uses the support columns, and a comparison with the pre-migration
# analysis goes through the ``*_committed`` ones.
COLUMNS = {"gamma2": "gamma2_mask", "P": "P", "D": "D", "F": "F",
           "A": "A", "S": "S", "T": "temperature_c"}
BLOCK = 5
BLOCKS = (3, 5, 7)
BOOT = 1000
N_PERM = 2000
N_PLACEBO = 8
SEED = 7
SMOOTH_DAYS = 30.0
#: ``n_post >= FITS_PER_COLUMN * k`` or the cell is refused.
FITS_PER_COLUMN = 3

#: The specs each tree may use; the guard still refuses per cell.
TREES = {
    "full": {
        "A_adjusted": ("M0", "M1", "M2", "M3"),
        "A_free": ("M0", "M2", "M3"),
        "sensitivity": ("M1_spline", "M1_smoothT", "M3_month"),
    },
    "reduced": {
        "A_adjusted": ("M0", "M1"),
        "A_free": ("M0", "M2"),
        "sensitivity": ("M1_smoothT",),
    },
    "minimal": {
        "A_adjusted": ("M0", "M1"),
        "A_free": ("M0",),
        "sensitivity": ("M1_smoothT",),
    },
}

#: Phrases the report must never produce: the design is co-movement, not cause.
FORBIDDEN = ("area-driven", "area driven", "caused by", "explained by",
             "because of the collapse", "mediated by")


# ---------------------------------------------------------------------------
# Configuration: one site at a time
# ---------------------------------------------------------------------------
def configure(name, table=None):
    """Point the module at one site."""
    global SITE, LABEL, TABLE, OUT_JSON, OUT_MD, EVENT, TARGETS, STRATUM_MODE
    global ANGLE_SPLIT, TREE, GEOMETRY_COVARIATE
    spec = SITES[name]
    SITE = name
    LABEL = spec["label"]
    TABLE = table or os.path.join(REPO, "data", name, "support",
                                  f"{name}_ocv_channels.csv")
    OUT_JSON = os.path.join(REPO, "data", name, "support",
                            f"{name}_ocv_decompose.json")
    OUT_MD = os.path.join(REPO, "data", name, "support",
                          f"{name}_ocv_decompose.md")
    EVENT = spec["event"]
    TARGETS = tuple(spec["targets"])
    STRATUM_MODE = spec["strata"]
    ANGLE_SPLIT = float(spec.get("angle_split", 40.0))
    TREE = spec["tree"]
    GEOMETRY_COVARIATE = bool(spec.get("geometry_covariate"))
    return spec


GEOMETRY_COVARIATE = False


def block_for(n_group):
    """Block length for one group: never more than a third of it, never above 5."""
    return max(2, min(BLOCK, int(n_group // 3)))


def fits(n_post, k):
    """The estimability rule: ``n_post >= 3 * k``."""
    return n_post >= FITS_PER_COLUMN * max(1, int(k))


# ---------------------------------------------------------------------------
# Data layer -- one dict per overpass, nothing imputed
# ---------------------------------------------------------------------------
def _naive(stamp):
    when = dt.datetime.fromisoformat(stamp)
    if when.tzinfo is not None:
        when = when.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return when


def stratum_of(row):
    """One stratum label per overpass, from the site's geometry rule."""
    if STRATUM_MODE == "single":
        return "ALL"
    orbit = (row.get("orbit_direction") or "").upper()
    short = "ASC" if orbit.startswith("ASC") else "DESC"
    if STRATUM_MODE == "channel":
        return short
    if short == "DESC":
        return "DESC"
    angle = row.get("incidence_angle_deg")
    if angle in (None, ""):
        return "ASC"
    return ("ASC low" if float(angle) < ANGLE_SPLIT else "ASC high")


def _number(row, key):
    raw = (row.get(key) or "").strip()
    if raw in ("", "NA", "null", "NaN"):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def load_rows(table=None):
    """Every overpass of the committed table, with the fields the models need."""
    table = table or TABLE
    rows = []
    with open(table, newline="") as handle:
        for raw in csv.DictReader(handle):
            when = _naive(raw["acquisition_ts"])
            values = {dim: _number(raw, COLUMNS[dim]) for dim in TARGETS}
            values["A"] = _number(raw, COLUMNS["A"])
            rows.append({
                "t": when,
                "doy": float(when.timetuple().tm_yday),
                "month": when.month,
                "year": when.year,
                "stratum": stratum_of(raw),
                "post": 0 if (raw.get("state") or "").startswith("pre") else 1,
                "angle": _number(raw, "incidence_angle_deg"),
                "T": _number(raw, COLUMNS["T"]),
                "values": values,
            })
    rows.sort(key=lambda row: row["t"])
    _add_smoothed_temperature(rows)
    return rows


def _add_smoothed_temperature(rows):
    """T_smooth30: 30-day centred mean of the overpass temperatures.

    A proxy for the thermal state of the deck, not a measurement of it -- the
    record has no hourly series, so this is the closest honest stand-in for what
    the spec calls T24/T72.
    """
    half = SMOOTH_DAYS / 2.0
    stamps = [row["t"] for row in rows]
    for row in rows:
        window = [other["T"] for other in rows
                  if other["T"] is not None
                  and abs((other["t"] - row["t"]).total_seconds()) <= half * 86400]
        row["T_smooth30"] = float(np.mean(window)) if window else None
        row["dT30"] = (None if row["T"] is None or not window
                        else float(row["T"] - np.mean(window)))
    return rows


def rows_of(rows, stratum):
    return [row for row in rows if row["stratum"] == stratum]


def stratum_display(rows):
    """Readable labels: an angle stratum carries its median incidence angle."""
    labels = {}
    for name in ordered_strata(rows):
        if STRATUM_MODE == "angle" and name.startswith("ASC "):
            angles = [row["angle"] for row in rows
                      if row["stratum"] == name and row.get("angle") is not None]
            labels[name] = (f"ASC {float(np.median(angles)):.1f}" if angles
                            else name)
        else:
            labels[name] = name
    return labels


def ordered_strata(rows):
    """The strata present, in a stable order (ASC low, ASC high, ASC, DESC, ALL)."""
    order = {"ASC low": 0, "ASC high": 1, "ASC": 2, "DESC": 3, "ALL": 4}
    present = {row["stratum"] for row in rows}
    return sorted(present, key=lambda name: (order.get(name, 9), name))


def counts(rows):
    pre = [row for row in rows if row["post"] == 0]
    post = [row for row in rows if row["post"] == 1]
    return {"pre": len(pre), "post": len(post)}


# ---------------------------------------------------------------------------
# Support and overlap
# ---------------------------------------------------------------------------
def overlap_coefficient(a, b, bins=20):
    """Histogram overlap of two samples, in [0, 1]; 0 = disjoint supports."""
    a = np.asarray([v for v in a if v is not None], float)
    b = np.asarray([v for v in b if v is not None], float)
    if a.size < 2 or b.size < 2:
        return None
    lo = min(a.min(), b.min())
    hi = max(a.max(), b.max())
    if hi <= lo:
        return 1.0
    edges = np.linspace(lo, hi, bins + 1)
    pa, _ = np.histogram(a, bins=edges, density=False)
    pb, _ = np.histogram(b, bins=edges, density=False)
    pa = pa / max(1, pa.sum())
    pb = pb / max(1, pb.sum())
    return float(np.minimum(pa, pb).sum())


def support_block(pre, post, use_log=False):
    """Common support of one variable: p5-p95 of pre, share of post outside."""
    a = np.asarray([v for v in pre if v is not None], float)
    b = np.asarray([v for v in post if v is not None], float)
    if a.size < 5 or b.size < 5:
        return {"available": False, "n": {"pre": int(a.size), "post": int(b.size)}}
    if use_log:
        a = np.log(a)
        b = np.log(b)
    lo, hi = np.percentile(a, [5, 95])
    inside = float(((b >= lo) & (b <= hi)).mean())
    return {
        "available": True,
        "n": {"pre": int(a.size), "post": int(b.size)},
        "pre_p5_p95": [float(lo), float(hi)],
        "pre_range": [float(a.min()), float(a.max())],
        "post_range": [float(b.min()), float(b.max())],
        "post_outside_share": float(1.0 - inside),
        "overlap": overlap_coefficient(a, b),
        "log_scale": bool(use_log),
    }


# ---------------------------------------------------------------------------
# Design matrices, fits, moving-block bootstrap
# ---------------------------------------------------------------------------
def transform(target, value):
    """g(Y): log for the strictly positive skewed dimensions, identity for S."""
    if value is None:
        return None
    if target in LOG_TARGETS:
        return math.log(value) if value > 0 else None
    return value


def rcs_basis(x, df=3):
    """Restricted cubic spline basis (linear + df-1 non-linear columns)."""
    x = np.asarray(x, float)
    knots = np.percentile(x, np.linspace(0, 100, df + 1))
    knots = np.unique(knots)
    if knots.size < 3:
        return [x]
    scale = (knots[-1] - knots[0]) ** 2

    def cube(t):
        return np.maximum(t, 0.0) ** 3

    cols = [x]
    for j in range(knots.size - 2):
        last, penultimate = knots[-1], knots[-2]
        term = (cube(x - knots[j])
                - cube(x - penultimate) * (last - knots[j]) / (last - penultimate)
                + cube(x - last) * (penultimate - knots[j]) / (last - penultimate))
        cols.append(term / scale)
    return cols


def design(rows, spec):
    """X and column names for one model spec."""
    n = len(rows)
    post = np.array([row["post"] for row in rows], float)
    angle = np.array([row["angle"] for row in rows], float)
    centre = float(np.mean(angle)) if n else 0.0
    cols = [np.ones(n)]
    names = ["intercept"]
    if spec.get("angle"):
        cols.append(angle - centre)
        names.append("angle")
    cols.append(post)
    names.append("post")
    if spec.get("area"):
        log_a = np.log(np.array([row["values"]["A"] for row in rows], float))
        basis = rcs_basis(log_a) if spec["area"] == "spline" else [log_a]
        for index, column in enumerate(basis):
            cols.append(column)
            names.append("logA" if len(basis) == 1 else f"logA_s{index + 1}")
    if spec.get("temp"):
        key = spec.get("temp_column", "T")
        cols.append(np.array([row[key] for row in rows], float))
        names.append(key)
    if spec.get("season") == "harmonic":
        doy = np.array([row["doy"] for row in rows], float)
        cols.append(np.sin(2.0 * math.pi * doy / 365.25))
        names.append("season_sin")
        cols.append(np.cos(2.0 * math.pi * doy / 365.25))
        names.append("season_cos")
    elif spec.get("season") == "month":
        months = np.array([row["month"] for row in rows])
        for month in range(2, 13):
            cols.append((months == month).astype(float))
            names.append(f"month_{month:02d}")
    if spec.get("post_angle"):
        cols.append(post * (angle - centre))
        names.append("post_x_angle")
    return np.column_stack(cols), names


def prepare(rows, target, spec):
    """Complete cases for one target and one spec -- never imputed."""
    usable = []
    key = spec.get("temp_column", "T")
    for row in rows:
        if transform(target, row["values"][target]) is None:
            continue
        if spec.get("area"):
            area = row["values"]["A"]
            if area is None or area <= 0:
                continue
        if spec.get("temp") and row.get(key) is None:
            continue
        usable.append(row)
    if len(usable) < 8:
        return None
    y = np.array([transform(target, row["values"][target]) for row in usable], float)
    X, names = design(usable, spec)
    return {"rows": usable, "y": y, "X": X, "names": names,
            "target": target, "spec": dict(spec)}


def block_index(n, block, rng):
    """One moving-block resample of ``range(n)`` -- consecutive blocks, in order."""
    if n <= 1:
        return np.zeros(n, int)
    block = max(1, min(int(block), n))
    starts = rng.integers(0, n - block + 1, int(np.ceil(n / block)))
    return np.concatenate([np.arange(s, s + block) for s in starts])[:n]


def fast_beta(y, X, column):
    """Least-squares coefficient only -- the bootstrap does not need HC3 or R2."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(beta[column])


def bootstrap_beta(prepared, name, block=None, boot=BOOT, seed=SEED):
    """Moving-block bootstrap of one coefficient of a prepared model.

    Blocks are drawn inside each state separately, so no block ever crosses the
    pre/post boundary (and, because a stratum is fitted on its own, no block
    crosses a geometry either). Rows enter time-ordered, so the blocks keep the
    serial dependence.
    """
    post = np.array([row["post"] for row in prepared["rows"]])
    groups = [np.where(post == 0)[0], np.where(post == 1)[0]]
    lengths = [block_for(g.size) if block is None else int(block) for g in groups]
    column = prepared["names"].index(name)
    observed = stats_mod.ols_fit(prepared["y"], prepared["X"], prepared["names"])
    rng = np.random.default_rng(seed)
    draws = np.empty(boot)
    for draw in range(boot):
        index = np.concatenate([g[block_index(g.size, length, rng)]
                                for g, length in zip(groups, lengths)])
        draws[draw] = fast_beta(prepared["y"][index], prepared["X"][index], column)
    share_le = float((draws <= 0.0).mean())
    share_ge = float((draws >= 0.0).mean())
    return {
        "beta": float(observed["beta"][column]),
        "ci": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
        "p_boot": min(1.0, 2.0 * min(share_le, share_ge)),
        "se_boot": float(draws.std(ddof=1)),
        "p_classical": float(observed["p_classical"][column]),
        "n": int(observed["n"]), "k": int(observed["k"]),
        "r2": float(observed["r2"]), "rss": float(observed["rss"]),
        "block": lengths[0] if len(set(lengths)) == 1 else lengths,
        "boot": int(boot), "seed": int(seed),
    }


def bootstrap_post(prepared, block=BLOCK, boot=BOOT, seed=SEED):
    """The post coefficient under the same bootstrap."""
    return bootstrap_beta(prepared, "post", block=block, boot=boot, seed=seed)


def bootstrap_delta(pre_values, post_values, block=BLOCK, n=BOOT, seed=SEED):
    """Block-bootstrap CI of Cliff's delta on the transformed scale."""
    rng = np.random.default_rng(seed)
    a = np.asarray(pre_values, float)
    b = np.asarray(post_values, float)
    if a.size < 2 or b.size < 2:
        return None
    draws = np.empty(n)
    for draw in range(n):
        draws[draw] = cliffs_delta(a[block_index(a.size, block, rng)],
                                   b[block_index(b.size, block, rng)])
    return {"delta": float(cliffs_delta(a, b)),
            "ci": [float(np.percentile(draws, 2.5)),
                   float(np.percentile(draws, 97.5))],
            "block": int(block), "boot": int(n), "seed": int(seed)}


# ---------------------------------------------------------------------------
# Test 1 -- the raw, unadjusted post signal
# ---------------------------------------------------------------------------
def raw_test(rows, target, boot=BOOT, seed=SEED):
    """Median ratio, Cliff's delta and exp(beta_post) from the state-only model."""
    pre, post = [], []
    for row in rows:
        value = transform(target, row["values"][target])
        if value is None:
            continue
        (post if row["post"] else pre).append(value)
    out = {"target": target, "n": {"pre": len(pre), "post": len(post)},
           "transform": "log" if target in LOG_TARGETS else "identity"}
    if len(pre) < 8 or len(post) < 8:
        out["available"] = False
        return out
    median_pre = float(np.median(pre))
    median_post = float(np.median(post))
    out.update({
        "available": True,
        "median_pre": median_pre, "median_post": median_post,
        "median_shift_g": median_post - median_pre,
        "ratio": (float(math.exp(median_post - median_pre))
                  if target in LOG_TARGETS else None),
        "delta": bootstrap_delta(pre, post, n=boot, seed=seed),
        "model": bootstrap_post(prepare(rows, target, {}), boot=boot, seed=seed),
    })
    return out


# ---------------------------------------------------------------------------
# Test 2 -- does A absorb the signature? plus the pre-fitted residual check
# ---------------------------------------------------------------------------
def ladder(rows, target, specs, boot=BOOT, seed=SEED):
    """Fit a sequence of nested specs and return the post coefficient of each.

    A spec whose design needs more post observations than the data can support
    (``n_post < 3 k``) is **refused** and named, not fitted.
    """
    out = {}
    for label, spec in specs.items():
        prepared = prepare(rows, target, spec)
        if prepared is None:
            out[label] = {"available": False, "spec": spec,
                          "reason": "no complete cases"}
            continue
        n_post = int(sum(1 for row in prepared["rows"] if row["post"]))
        k = int(prepared["X"].shape[1])
        if not fits(n_post, k):
            out[label] = {"available": False, "spec": spec, "k": k,
                          "n_post": n_post,
                          "reason": f"refused: n_post {n_post} < "
                                    f"{FITS_PER_COLUMN} x k {k}"}
            continue
        block = bootstrap_post(prepared, block=None, boot=boot, seed=seed)
        block.update({"available": True, "spec": spec,
                      "n_pre": int(sum(1 for r in prepared["rows"] if not r["post"])),
                      "n_post": int(sum(1 for r in prepared["rows"] if r["post"]))})
        out[label] = block
    return out


def attenuation(raw, adjusted):
    """1 - |beta_adjusted| / |beta_raw|; descriptive only, no threshold."""
    if not raw or not adjusted or not raw.get("available") or not adjusted.get("available"):
        return None
    if raw["beta"] == 0:
        return None
    return float(1.0 - abs(adjusted["beta"]) / abs(raw["beta"]))


def residual_test(rows, target, boot=BOOT, seed=SEED):
    """Fit g(Y) ~ log A on the pre data only, then compare the residuals.

    Linear ``f`` only: a spline basis built on the post data would put different
    knots on the same curve, so the comparison would mix two models. The linear
    form is the one the pre-fitted residual check is defined for.
    """
    pre_rows = [row for row in rows if row["post"] == 0]
    post_rows = [row for row in rows if row["post"] == 1]
    spec = {"area": "linear"}
    pre_prep = prepare(pre_rows, target, spec)
    post_prep = prepare(post_rows, target, spec)
    if pre_prep is None or post_prep is None:
        return {"available": False}
    fit = stats_mod.ols_fit(pre_prep["y"], pre_prep["X"], pre_prep["names"])
    beta = np.array(fit["beta"], float)
    resid_pre = pre_prep["y"] - pre_prep["X"] @ beta
    resid_post = post_prep["y"] - post_prep["X"] @ beta
    delta = bootstrap_delta(resid_pre, resid_post, n=boot, seed=seed)
    return {
        "available": True,
        "n": {"pre": int(resid_pre.size), "post": int(resid_post.size)},
        "pre_median_residual": float(np.median(resid_pre)),
        "post_median_residual": float(np.median(resid_post)),
        "delta": delta,
        "pre_sd": float(resid_pre.std(ddof=1)),
    }


# ---------------------------------------------------------------------------
# Test 3 -- geometry: the joint ASC model and the interaction
# ---------------------------------------------------------------------------
def joint_asc(rows, target, spec, terms=("post",), boot=BOOT, seed=SEED):
    """Fit one spec on both ASC tracks, bootstrap the named terms, F-test the
    interaction. With two angle values the interaction is the stratum contrast,
    so this adds no information beyond the stratified fits -- it is reported as
    family B, not as a second headline."""
    asc = [row for row in rows if row["stratum"] != "DESC"]
    prepared = prepare(asc, target, spec)
    if prepared is None:
        return {"available": False}
    out = {"available": True, "n": int(prepared["y"].size),
           "n_pre": int(sum(1 for r in prepared["rows"] if not r["post"])),
           "n_post": int(sum(1 for r in prepared["rows"] if r["post"])),
           "terms": {}}
    for term in terms:
        if term in prepared["names"]:
            out["terms"][term] = bootstrap_beta(prepared, term, boot=boot, seed=seed)
    if "post_angle" in spec:
        restricted_spec = {key: value for key, value in spec.items()
                           if key != "post_angle"}
        restricted = prepare(asc, target, restricted_spec)
        if restricted is not None and restricted["y"].size == prepared["y"].size:
            full_fit = stats_mod.ols_fit(prepared["y"], prepared["X"], prepared["names"])
            rest_fit = stats_mod.ols_fit(restricted["y"], restricted["X"],
                                    restricted["names"])
            out["f_interaction"] = stats_mod.f_test(
                rest_fit["rss"], full_fit["rss"], 1,
                max(1, prepared["y"].size - prepared["X"].shape[1]))
    return out


# ---------------------------------------------------------------------------
# Test 4/5 -- nested F-tests and the block-length ladder
# ---------------------------------------------------------------------------
def nested_test(rows, target, restricted_spec, full_spec):
    """Partial F-test between two nested specs on identical rows."""
    left = prepare(rows, target, restricted_spec)
    right = prepare(rows, target, full_spec)
    if left is None or right is None:
        return {"available": False}
    if left["y"].size != right["y"].size:
        return {"available": False, "reason": "different complete-case sets"}
    lf = stats_mod.ols_fit(left["y"], left["X"], left["names"])
    rf = stats_mod.ols_fit(right["y"], right["X"], right["names"])
    dfn = right["X"].shape[1] - left["X"].shape[1]
    dfd = max(1, right["y"].size - right["X"].shape[1])
    out = stats_mod.f_test(lf["rss"], rf["rss"], dfn, dfd)
    out.update({"available": True, "n": int(right["y"].size),
                "rss_restricted": lf["rss"], "rss_full": rf["rss"]})
    return out


def block_sensitivity(rows, target, spec, boot=BOOT, seed=SEED):
    """The same post estimate at block lengths 3, 5 and 7."""
    prepared = prepare(rows, target, spec)
    if prepared is None:
        return {"available": False}
    out = {"available": True, "by_block": {}}
    for block in BLOCKS:
        fit = bootstrap_post(prepared, block=block, boot=boot, seed=seed)
        out["by_block"][str(block)] = {"beta": fit["beta"], "ci": fit["ci"],
                                       "se_boot": fit["se_boot"]}
    widths = [ci["ci"][1] - ci["ci"][0] for ci in out["by_block"].values()]
    signs = {math.copysign(1.0, ci["beta"]) for ci in out["by_block"].values()}
    crossing = [ci["ci"][0] <= 0.0 <= ci["ci"][1] for ci in out["by_block"].values()]
    out["stable_sign"] = len(signs) == 1
    out["same_verdict"] = len(set(crossing)) == 1
    out["width_ratio"] = (max(widths) / min(widths)) if min(widths) > 0 else None
    return out


# ---------------------------------------------------------------------------
# Test 6 -- serial dependence of the residuals
# ---------------------------------------------------------------------------
def residual_diagnostics(prepared):
    """Lag-1, runs test and ACF of the residuals, in time order."""
    fit = stats_mod.ols_fit(prepared["y"], prepared["X"], prepared["names"])
    resid = prepared["y"] - prepared["X"] @ np.array(fit["beta"], float)
    runs = espoo.runs_test(resid) or {}
    acf = {}
    for lag in range(1, 9):
        if resid.size > lag + 2:
            acf[str(lag)] = espoo.lag1_autocorr(resid) if lag == 1 else None
    for lag in range(1, 9):
        if resid.size > lag + 2:
            a, b = resid[:-lag], resid[lag:]
            if a.std() > 0 and b.std() > 0:
                acf[str(lag)] = float(np.corrcoef(a, b)[0, 1])
    return {"available": True, "n": int(resid.size),
            "lag1": espoo.lag1_autocorr(resid),
            "runs_p": runs.get("p"), "acf": acf,
            "resid": resid}


# ---------------------------------------------------------------------------
# Test 3/4 -- the seasonal variants used by the ladders
# ---------------------------------------------------------------------------
def specs_for_ladders(tree=None):
    """The specs of the two ladders, cut down to what this site's tree allows."""
    everything = {
        "A_adjusted": {
            "M0": {},
            "M1": {"area": "linear"},
            "M2": {"area": "linear", "temp": True},
            "M3": {"area": "linear", "temp": True, "season": "harmonic"},
        },
        "A_free": {
            "M0": {},
            "M2": {"temp": True},
            "M3": {"temp": True, "season": "harmonic"},
        },
        "sensitivity": {
            "M1_spline": {"area": "spline"},
            "M1_smoothT": {"area": "linear", "temp": True,
                           "temp_column": "T_smooth30"},
            "M3_month": {"area": "linear", "temp": True, "season": "month"},
        },
    }
    allowed = TREES[tree or TREE]
    return {group: {label: spec for label, spec in everything[group].items()
                    if label in allowed[group]}
            for group in everything}


def month_matched(rows, target):
    """Same-calendar-month contrast: post versus pre, month by month.

    A seasonal control that does not assume a sinusoid. Only months with at
    least three values on both sides contribute, and the reported ratio is the
    median of the monthly ratios, so no month can dominate by sample size.
    """
    by_month = {}
    for row in rows:
        value = transform(target, row["values"][target])
        if value is None:
            continue
        bucket = by_month.setdefault(row["month"], {"pre": [], "post": []})
        bucket["post" if row["post"] else "pre"].append(value)
    ratios, used = [], []
    for month in sorted(by_month):
        pre = by_month[month]["pre"]
        post = by_month[month]["post"]
        if len(pre) < 3 or len(post) < 3:
            continue
        shift = float(np.median(post) - np.median(pre))
        ratios.append(math.exp(shift) if target in LOG_TARGETS else shift)
        used.append({"month": month, "n": [len(pre), len(post)],
                     "shift_g": shift})
    if not ratios:
        return {"available": False}
    return {"available": True, "months": used,
            "median_shift_g": float(np.median([u["shift_g"] for u in used])),
            "median_ratio": (float(np.median(ratios))
                             if target in LOG_TARGETS else None),
            "n_months": len(used)}


# ---------------------------------------------------------------------------
# Test 7 -- the multivariate shift, with a circular-shift null
# ---------------------------------------------------------------------------
def multivariate_shift(rows, draw=N_PERM, seed=SEED):
    """Is the post residual vector shifted as a whole?

    Residuals of M3- (temperature and season, no area, no post) are fitted on the
    **pre** data only, standardised by their pre sd, and stacked into one matrix
    (complete cases only -- a missing gamma2 drops the row, it is never read as
    zero). The null is a circular shift of that matrix against the event
    boundary: the serial dependence and the vector's internal correlation stay
    exactly as observed, only the alignment with the event is randomised.
    """
    spec = {"temp": True, "season": "harmonic"}
    columns, keep, names = [], None, []
    for target in TARGETS:
        prepared = prepare(rows, target, spec)
        if prepared is None:
            continue
        pre_mask = np.array([row["post"] == 0 for row in prepared["rows"]])
        if pre_mask.sum() < 10:
            continue
        fit = stats_mod.ols_fit(prepared["y"][pre_mask], prepared["X"][pre_mask],
                           prepared["names"])
        beta = np.array(fit["beta"], float)
        resid = prepared["y"] - prepared["X"] @ beta
        sd = float(resid[pre_mask].std(ddof=1)) or 1.0
        series = np.full(len(rows), np.nan)
        index = {id(row): i for i, row in enumerate(rows)}
        for i, row in enumerate(prepared["rows"]):
            series[index[id(row)]] = resid[i] / sd
        columns.append(series)
        names.append(target)
    if not columns:
        return {"available": False}
    matrix = np.column_stack(columns)
    complete = ~np.isnan(matrix).any(axis=1)
    matrix = matrix[complete]
    post_mask = np.array([row["post"] == 1 for row, ok in zip(rows, complete) if ok])
    n_pre = int((~post_mask).sum())
    if matrix.shape[0] < 20 or n_pre < 10:
        return {"available": False}

    def statistic(values, mask):
        pre = values[~mask]
        post = values[mask]
        return float(((post.mean(axis=0) - pre.mean(axis=0)) ** 2).sum())

    observed = statistic(matrix, post_mask)
    rng = np.random.default_rng(seed)
    n = matrix.shape[0]
    hits = 0
    for _ in range(draw):
        shift = int(rng.integers(1, n))
        hits += statistic(np.roll(matrix, shift, axis=0), post_mask) >= observed
    return {"available": True, "n": int(n), "n_pre": n_pre,
            "n_post": int(post_mask.sum()), "dims": names,
            "statistic": observed, "p": float((1 + hits) / (1 + draw)),
            "draws": int(draw), "seed": int(seed)}


# ---------------------------------------------------------------------------
# Test 8 -- placebos, month matching and the common-support subset
# ---------------------------------------------------------------------------
def signature_strength(rows, spec):
    """Sum of squared standardised post effects across the five targets.

    One number per split, so a real event can be ranked against pseudo events
    that cannot have the same group sizes. The standardisation is by the sd of
    the transformed outcome, so the number is comparable across splits.
    """
    total, parts = 0.0, {}
    for target in TARGETS:
        prepared = prepare(rows, target, spec)
        if prepared is None or prepared["y"].size < 12:
            continue
        fit = stats_mod.ols_fit(prepared["y"], prepared["X"], prepared["names"])
        column = prepared["names"].index("post")
        sd = float(prepared["y"].std(ddof=1)) or 1.0
        z = float(fit["beta"][column]) / sd
        parts[target] = z
        total += z * z
    return {"statistic": float(total), "z": parts, "dims": len(parts)}


def _standardised_effect(rows, target, spec):
    """|beta_post| in units of the outcome's sd -- comparable across splits."""
    prepared = prepare(rows, target, spec)
    if prepared is None or prepared["y"].size < 10:
        return None
    beta, *_ = np.linalg.lstsq(prepared["X"], prepared["y"], rcond=None)
    scale = float(prepared["y"].std(ddof=1)) or 1.0
    return abs(float(beta[prepared["names"].index("post")])) / scale


def placebo_target(rows, target, spec, n_events=N_PLACEBO, draws=20, seed=SEED):
    """Pseudo events inside the pre period, compared at the **same group size**.

    A pseudo event carved out of the pre period rarely has the real post-group
    size, and a larger group is more power, not more event. Both sides are
    therefore cut to ``m = min(pseudo_post, real_post)`` in every draw, and the
    real effect at that size is what gets ranked against the pseudo effects.
    """
    pre = sorted([row for row in rows if row["post"] == 0], key=lambda r: r["t"])
    post = [row for row in rows if row["post"] == 1]
    if len(pre) < 24 or len(post) < 3:
        return {"available": False}
    real_full = _standardised_effect(rows, target, spec)
    rng = np.random.default_rng(seed)
    events, real_matched = [], []
    for share in np.linspace(0.40, 0.80, n_events):
        index = max(2, min(len(pre) - 6, int(round(share * len(pre)))))
        cut = pre[index]["t"]
        pseudo_base = [{**row, "post": 0} for row in pre if row["t"] < cut]
        pseudo_post = [{**row, "post": 1} for row in pre if row["t"] >= cut]
        size = min(len(pseudo_post), len(post))
        if size < 3 or len(pseudo_base) < 8:
            continue
        pseudo_values, real_values = [], []
        for _ in range(draws):
            taken_pseudo = rng.choice(len(pseudo_post), size=size, replace=False)
            taken_real = rng.choice(len(post), size=size, replace=False)
            pseudo = pseudo_base + [pseudo_post[i] for i in taken_pseudo]
            real = pre + [post[i] for i in taken_real]
            pseudo_values.append(_standardised_effect(pseudo, target, spec))
            real_values.append(_standardised_effect(real, target, spec))
        pseudo_values = [v for v in pseudo_values if v is not None]
        real_values = [v for v in real_values if v is not None]
        if not pseudo_values or not real_values:
            continue
        events.append({"share": float(share), "cut": cut.isoformat(),
                       "size": int(size),
                       "z": float(np.median(pseudo_values))})
        real_matched.append(float(np.median(real_values)))
    if not events:
        return {"available": False}
    real_z = float(np.median(real_matched))
    rank = 1 + sum(1 for event in events if event["z"] >= real_z)
    return {"available": True, "n_post": len(post), "n_pre": len(pre),
            "real_z_full": real_full, "real_z_matched": real_z,
            "events": events, "rank_of_real": rank, "of": len(events),
            "draws_per_event": draws, "seed": seed, "matched": True}


def placebo_panel(rows, specs, n_events=N_PLACEBO, seed=SEED):
    """Signature-level pseudo events, matched to the real post-group size."""
    pre = sorted([row for row in rows if row["post"] == 0], key=lambda r: r["t"])
    post = [row for row in rows if row["post"] == 1]
    if len(pre) < 24 or len(post) < 3:
        return {"available": False}
    rng = np.random.default_rng(seed)
    out = {"available": True, "n_pre": len(pre), "n_post": len(post),
           "by_spec": {}}
    for label, spec in specs.items():
        real_full = signature_strength(rows, spec)
        events, real_matched = [], []
        for share in np.linspace(0.40, 0.80, n_events):
            index = max(2, min(len(pre) - 6, int(round(share * len(pre)))))
            cut = pre[index]["t"]
            pseudo_base = [{**row, "post": 0} for row in pre if row["t"] < cut]
            pseudo_post = [{**row, "post": 1} for row in pre if row["t"] >= cut]
            size = min(len(pseudo_post), len(post))
            if size < 3:
                continue
            taken_pseudo = rng.choice(len(pseudo_post), size=size, replace=False)
            taken_real = rng.choice(len(post), size=size, replace=False)
            pseudo = pseudo_base + [pseudo_post[i] for i in taken_pseudo]
            real = pre + [post[i] for i in taken_real]
            events.append({"share": float(share), "cut": cut.isoformat(),
                           "size": int(size),
                           "statistic": signature_strength(pseudo, spec)["statistic"]})
            real_matched.append(signature_strength(real, spec)["statistic"])
        if not events:
            continue
        real_value = float(np.median(real_matched))
        stronger = sum(1 for event in events
                       if event["statistic"] >= real_value)
        out["by_spec"][label] = {
            "real": real_full, "real_matched": real_value, "events": events,
            "matched_post_size": int(min(len(post), max(e["size"] for e in events))),
            "placebo_at_least_as_strong": int(stronger),
            "rank_of_real": int(stronger) + 1,
            "rank_of_real_matched": int(stronger) + 1,
            "of": int(len(events)),
        }
    return out


def common_support(rows, share=5.0):
    """The pre p5-p95 box of log A and T; returns the rows inside it."""
    pre = [row for row in rows if row["post"] == 0]
    area = np.array([math.log(row["values"]["A"]) for row in pre
                     if row["values"]["A"]], float)
    temp = np.array([row["T"] for row in pre if row["T"] is not None], float)
    if area.size < 5 or temp.size < 5:
        return None, {"available": False}
    lo_a, hi_a = np.percentile(area, [share, 100 - share])
    lo_t, hi_t = np.percentile(temp, [share, 100 - share])
    kept, dropped = [], 0
    for row in rows:
        value = row["values"]["A"]
        if value is None or value <= 0 or row["T"] is None:
            dropped += 1
            continue
        inside = (lo_a <= math.log(value) <= hi_a and lo_t <= row["T"] <= hi_t)
        if inside:
            kept.append(row)
        else:
            dropped += 1
    post_total = sum(1 for row in rows if row["post"] == 1)
    post_kept = sum(1 for row in kept if row["post"] == 1)
    meta = {"available": True, "box": {"logA": [float(lo_a), float(hi_a)],
                                        "T": [float(lo_t), float(hi_t)]},
            "n": {"total": len(rows), "kept": len(kept), "dropped": dropped,
                  "post_total": post_total, "post_kept": post_kept},
            "post_share_kept": float(post_kept / post_total) if post_total else None}
    return kept, meta


# ---------------------------------------------------------------------------
# The verdict, from the pair of ladders
# ---------------------------------------------------------------------------
def final_model(steps):
    """The richest model of a ladder that was actually fitted, and its label."""
    for label in ("M3", "M2", "M1", "M0"):
        fit = steps.get(label)
        if fit and fit.get("available"):
            return fit, label
    return None, None


def excludes_null(target, fit):
    """Does the block-bootstrap CI of the post effect exclude no effect?"""
    if not fit or not fit.get("available"):
        return None
    low, high = fit["ci"]
    return bool(low > 0.0 or high < 0.0)


def verdict(target, adjusted, free):
    """Text for one cell, from whether each ladder excludes no effect."""
    a = excludes_null(target, adjusted)
    f = excludes_null(target, free)
    if a is None or f is None:
        return "not estimable"
    if a and f:
        if (adjusted["beta"] > 0) != (free["beta"] > 0):
            return ("sign flips between the ladders -- co-movement with the "
                    "observed area, not a stable direction")
        return ("survives area, temperature and season "
                "(both ladders exclude no effect, same sign)")
    if (not a) and f:
        return ("co-moves with the observed area "
                "(A-free survives, A-adjusted does not)")
    if a and (not f):
        return "unstable (A-adjusted survives, A-free does not)"
    return "not distinguishable from no effect after the covariates"


# ---------------------------------------------------------------------------
# The whole report
# ---------------------------------------------------------------------------
def analyse(table=None, boot=BOOT, placebo=True, permutation=True):
    table = table or TABLE
    rows = load_rows(table)
    specs = specs_for_ladders()
    report = {
        "schema": SCHEMA, "site": SITE, "label": LABEL, "tier": TREE,
        "table": {"path": os.path.relpath(table, REPO),
                  "sha256": sha256_of(table)},
        "event": EVENT.isoformat(),
        "targets": list(TARGETS), "log_targets": list(LOG_TARGETS),
        "strata": {},
        "rules": {
            "area_is": "explanatory variable; never a target in the same model",
            "word_for_attenuation": "co-moves with the observed area",
            "temperature": "T0 primary, T_smooth30 sensitivity; T24/T72 do not "
                           "exist in this record (one temperature per overpass)",
            "block": BLOCK, "blocks": list(BLOCKS), "boot": boot,
            "permutations": N_PERM if permutation else 0,
            "placebos": N_PLACEBO if placebo else 0,
            "seed": SEED,
            "strata_mode": STRATUM_MODE,
            "estimability": f"n_post >= {FITS_PER_COLUMN} x k, else refused",
            "families": {"A": "the post effects of one stratum, Holm",
                         "B": "post x geometry",
                         "C": "exploratory, never mixed into A or B"},
        },
    }
    main_spec = {"area": "linear", "temp": True, "season": "harmonic"}
    if GEOMETRY_COVARIATE:
        main_spec = {**main_spec, "angle": True}
    geometry_spec = {**main_spec, "angle": True, "post_angle": True}
    collected = {}
    display = stratum_display(rows)
    report["stratum_map"] = display
    for stratum in ordered_strata(rows):
        sub = rows_of(rows, stratum)
        pre = [row for row in sub if row["post"] == 0]
        post = [row for row in sub if row["post"] == 1]
        block = {
            "counts": counts(sub),
            "missing": {target: {"pre": sum(1 for r in pre if r["values"][target] is None),
                                 "post": sum(1 for r in post if r["values"][target] is None)}
                        for target in TARGETS},
            "support": {
                "logA": support_block([r["values"]["A"] for r in pre],
                                      [r["values"]["A"] for r in post], use_log=True),
                "T": support_block([r["T"] for r in pre], [r["T"] for r in post]),
                "T_smooth30": support_block([r["T_smooth30"] for r in pre],
                                            [r["T_smooth30"] for r in post]),
            },
            "raw": {}, "ladders": {}, "attenuation": {}, "nested": {},
            "residual": {}, "block_sensitivity": {}, "diagnostics": {},
            "month_matched": {}, "common_support": {}, "geometry": {},
            "verdict": {}, "final_models": {},
        }
        for target in TARGETS:
            block["raw"][target] = raw_test(sub, target, boot=boot)
            block["ladders"][target] = {
                "A_adjusted": ladder(sub, target, specs["A_adjusted"], boot=boot),
                "A_free": ladder(sub, target, specs["A_free"], boot=boot),
                "sensitivity": ladder(sub, target, specs["sensitivity"], boot=boot),
            }
            adjusted = block["ladders"][target]["A_adjusted"]
            free = block["ladders"][target]["A_free"]
            final_adj, label_adj = final_model(adjusted)
            final_free, label_free = final_model(free)
            block["final_models"][target] = {"A_adjusted": label_adj,
                                             "A_free": label_free}
            block["attenuation"][target] = {
                "by_area": attenuation(adjusted.get("M0"), adjusted.get("M1")),
                "by_area_temp_season": attenuation(adjusted.get("M0"), final_adj),
                "free_by_temp_season": attenuation(free.get("M0"), final_free),
            }
            block["nested"][target] = {
                "temperature": nested_test(sub, target, {"area": "linear"},
                                           {"area": "linear", "temp": True}),
                "season": nested_test(sub, target, {"area": "linear", "temp": True},
                                      main_spec),
                "season_without_area": nested_test(sub, target, {"temp": True},
                                                   {"temp": True, "season": "harmonic"}),
            }
            block["residual"][target] = residual_test(sub, target, boot=boot)
            block["block_sensitivity"][target] = block_sensitivity(sub, target,
                                                                   main_spec, boot=boot)
            block["month_matched"][target] = month_matched(sub, target)
            prepared = prepare(sub, target, main_spec)
            diag = residual_diagnostics(prepared) if prepared else {"available": False}
            diag.pop("resid", None)
            block["diagnostics"][target] = diag
            block["verdict"][target] = verdict(target, final_adj, final_free)
        pvals = []
        for target in TARGETS:
            fit, _label = final_model(block["ladders"][target]["A_adjusted"])
            pvals.append(fit.get("p_boot") if fit else None)
        adjusted_p, family = holm(pvals)
        for target, p in zip(TARGETS, adjusted_p):
            fit, _label = final_model(block["ladders"][target]["A_adjusted"])
            if fit is not None:
                fit["holm_p"] = p
            collected[(display.get(stratum, stratum), target)] = (
                fit.get("p_boot") if fit else None)
        block["holm"] = {"family_A": "the five targets in this stratum",
                         "family_size": family,
                         "alpha": 0.05}
        kept, meta = common_support(sub)
        block["common_support"] = meta
        if meta.get("available") and kept:
            block["common_support"]["models"] = {
                target: (bootstrap_post(prepare(kept, target, main_spec), boot=boot)
                         if prepare(kept, target, main_spec) else {"available": False})
                for target in TARGETS}
        if permutation:
            block["multivariate"] = multivariate_shift(sub)
        if placebo:
            block["placebo"] = placebo_panel(
                sub, {"M0": {}, "M3_free": {"temp": True, "season": "harmonic"}})
        report["strata"][display.get(stratum, stratum)] = block
    if STRATUM_MODE == "angle":
        report["geometry_joint_asc"] = {
            target: joint_asc(rows, target, geometry_spec,
                              terms=("post", "angle", "post_x_angle"), boot=boot)
            for target in TARGETS}
    keys = sorted(collected)
    all_p, size = holm([collected[key] for key in keys])
    report["holm_all_strata"] = {"family": "all target x stratum post effects",
                                 "size": size,
                                 "values": {f"{s}|{t}": p for (s, t), p in zip(keys, all_p)}}
    return report


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
def _f(value, digits=3, signed=False):
    if value is None:
        return "—"
    return f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"


def _ratio(value, digits=3):
    if value is None:
        return "—"
    return f"{math.exp(value):.{digits}f}"


def _ci(pair, digits=3, ratio=False):
    if not pair:
        return "—"
    fmt = (lambda v: f"{math.exp(v):.{digits}f}") if ratio else (lambda v: f"{v:+.{digits}f}")
    return f"[{fmt(pair[0])}, {fmt(pair[1])}]"


def _p(value):
    if value is None:
        return "—"
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def _effect(target, fit, digits=3, ratio=False):
    """One post effect, as a ratio for the log targets and as pixels for S."""
    if not fit:
        return "—"
    beta = fit.get("beta")
    if beta is None:
        return "—"
    if ratio and target in LOG_TARGETS:
        return _ratio(beta, digits)
    return _f(beta, digits, signed=True)


def to_markdown(report):
    lines = []
    add = lines.append
    add(f"# {report['label']} — decomposition of the post-collapse signature")
    add("")
    add(f"* table `{report['table']['path']}` (sha256 `{report['table']['sha256'][:16]}…`)")
    add(f"* event **{report['event'][:10]}**")
    add(f"* strata: {', '.join(report['strata'])}")
    add(f"* targets: {', '.join(report['targets'])} — `A` is the explanatory variable, never a target")
    add(f"* moving-block bootstrap: block {report['rules']['block']} "
        f"(ladder {', '.join(str(b) for b in report['rules']['blocks'])}), "
        f"{report['rules']['boot']} draws, seed {report['rules']['seed']}")
    add("")
    add("## How to read the two ladders")
    add("")
    add("`+A` and `A-free` are two parallel ladders over the same data:")
    add("")
    add("```")
    add("A-adjusted   M0 (state) -> M1 (+log A) -> M2 (+T) -> M3 (+temperature and season)")
    add("A-free       M0 (state) -> M2- (+T)   -> M3- (+temperature and season)")
    add("```")
    add("")
    add("A post effect that **attenuates under A adjustment while staying put in the "
        "A-free ladder** *co-moves with the observed area*. That is a statement about "
        "co-movement, not about cause: the area is itself moved by the collapse "
        "(DESC area falls by 77 %), so conditioning on it can remove the very "
        "change being measured — and, for a mediator, can also *create* "
        "association where the raw data show none. Both directions occur in this "
        "record and are flagged in the verdict column.")
    add("")
    add("## The minimal table")
    add("")
    add("| stratum | target | raw | +A | A-free | +T/season | final 95% CI | Holm |")
    add("|---|---|---|---|---|---|---|---|")
    for stratum, block in report["strata"].items():
        for target in report["targets"]:
            raw = block["raw"].get(target) or {}
            adjusted = block["ladders"][target]["A_adjusted"]
            free = block["ladders"][target]["A_free"]
            labels = (block.get("final_models") or {}).get(target) or {}
            m1 = adjusted.get("M1")
            m3a = adjusted.get(labels.get("A_adjusted") or "") or {}
            m3f = free.get(labels.get("A_free") or "") or {}
            cell_raw = (_ratio(math.log(raw["ratio"]), 3) if raw.get("ratio")
                        else _f(raw.get("median_shift_g"), 2, signed=True))
            add(f"| {stratum} | `{target}` | {cell_raw} | {_effect(target, m1, 3, True)} | "
                f"{_effect(target, m3f, 3, True)} | {_effect(target, m3a, 3, True)} | "
                f"{_ci(m3a.get('ci') if m3a else None, 3, target in LOG_TARGETS)} | "
                f"{_p(m3a.get('holm_p') if m3a else None)} |")
    add("")
    add("Ratio columns are `exp(beta_post)` for `gamma2, P, D, F` (adjusted post "
        "level as a fraction of pre); `S` is in pixels and shows the signed "
        "shift. `+T/season` is the richest A-adjusted model this site could fit "
        "(see the tier and the estimability rule above); the Holm column is "
        "family A: this site's targets **within** a stratum.")
    add("")
    for stratum, block in report["strata"].items():
        add(f"## {stratum}")
        add("")
        add(f"*n* = {block['counts']['pre']} pre / {block['counts']['post']} post. "
            "Missing values (never imputed): "
            + ", ".join(f"`{t}` {v['pre']}/{v['post']}"
                        for t, v in block["missing"].items() if v["pre"] or v["post"])
            + (" — none" if not any(v["pre"] or v["post"]
                                    for v in block["missing"].values()) else ""))
        add("")
        add("### Common support (pre 5th-95th percentile)")
        add("")
        add("| variable | pre p5 | pre p95 | post range | post outside | overlap |")
        add("|---|---|---|---|---|---|")
        for name, support in block["support"].items():
            if not support.get("available"):
                continue
            add(f"| `{name}` | {_f(support['pre_p5_p95'][0], 2)} | "
                f"{_f(support['pre_p5_p95'][1], 2)} | "
                f"[{_f(support['post_range'][0], 2)}, {_f(support['post_range'][1], 2)}] | "
                f"**{100 * support['post_outside_share']:.1f} %** | "
                f"{_f(support['overlap'], 2)} |")
        add("")
        kept = block["common_support"].get("n", {})
        add(f"Common-support subset (the pre p5-p95 box of `log A` **and** `T`): "
            f"{kept.get('kept')} of {kept.get('total')} rows kept, "
            f"{kept.get('post_kept')} of {kept.get('post_total')} post rows "
            f"({100 * (block['common_support'].get('post_share_kept') or 0):.0f} %).")
        add("")
        add("### Raw, unadjusted (reference only)")
        add("")
        add("| target | n pre/post | median shift (g) | ratio | Cliff's delta [CI] | exp(beta) M0 [CI] |")
        add("|---|---|---|---|---|---|")
        for target in report["targets"]:
            raw = block["raw"].get(target) or {}
            if not raw.get("available"):
                add(f"| `{target}` | — | — | — | — | — |")
                continue
            delta = raw.get("delta") or {}
            model = raw.get("model") or {}
            add(f"| `{target}` | {raw['n']['pre']}/{raw['n']['post']} | "
                f"{_f(raw['median_shift_g'], 3, signed=True)} | "
                f"{(lambda r: '—' if r is None else f'{r:.3f}')(raw.get('ratio'))} | "
                f"{_f(delta.get('delta'), 2, signed=True)} {_ci(delta.get('ci'), 2)} | "
                f"{_effect(target, model, 3, True)} "
                f"{_ci(model.get('ci'), 3, target in LOG_TARGETS)} |")
        add("")
        _stratum_details(add, stratum, block, report)
    _global_sections(add, report)
    return "\n".join(lines)


def _stratum_details(add, stratum, block, report):
    targets = report["targets"]
    for ladder_name in ("A_adjusted", "A_free"):
        add(f"### {ladder_name.replace('_', ' ')} ladder")
        add("")
        add("| target | model | post effect | 95% CI | p (bootstrap) | n pre/post | n |")
        add("|---|---|---|---|---|---|---|")
        for target in targets:
            steps = block["ladders"][target][ladder_name]
            for label, fit in steps.items():
                if not fit.get("available"):
                    add(f"| `{target}` | {label} | not estimable | — | — | — | — |")
                    continue
                add(f"| `{target}` | {label} | {_effect(target, fit, 3, True)} | "
                    f"{_ci(fit['ci'], 3, target in LOG_TARGETS)} | {_p(fit.get('p_boot'))} | "
                    f"{fit.get('n_pre', '—')}/{fit.get('n_post', '—')} | {fit['n']} |")
        add("")
    add("### Attenuation (descriptive, no threshold)")
    add("")
    add("| target | raw -> +A | raw -> +A+T+season | raw -> A-free final | verdict |")
    add("|---|---|---|---|---|")
    for target in targets:
        att = block["attenuation"][target]
        add(f"| `{target}` | {_f(att.get('by_area'), 2)} | "
            f"{_f(att.get('by_area_temp_season'), 2)} | "
            f"{_f(att.get('free_by_temp_season'), 2)} | "
            f"{block['verdict'][target]} |")
    add("")
    add("### Nested tests (partial F)")
    add("")
    add("| target | +temperature | +season (A-adjusted) | +season (A-free) |")
    add("|---|---|---|---|")
    for target in targets:
        nested = block["nested"][target]
        cells = []
        for key in ("temperature", "season", "season_without_area"):
            test = nested.get(key) or {}
            if not test.get("available"):
                cells.append("—")
            else:
                cells.append(f"F {_f(test.get('F'), 2)} (df {test.get('dfn')}), "
                             f"p {_p(test.get('p'))}")
        add(f"| `{target}` | " + " | ".join(cells) + " |")
    add("")
    add("### Block-length stability of the final post effect")
    add("")
    add("| target | block 3 | block 5 | block 7 | same verdict |")
    add("|---|---|---|---|---|")
    for target in targets:
        sense = block["block_sensitivity"][target]
        if not sense.get("available"):
            add(f"| `{target}` | — | — | — | — |")
            continue
        by_block = sense["by_block"]
        cells = []
        for length in (3, 5, 7):
            entry = by_block[str(length)]
            point = (_ratio(entry["beta"]) if target in LOG_TARGETS
                     else _f(entry["beta"], 3, signed=True))
            cells.append(f"{point} "
                         f"{_ci(entry['ci'], 3, target in LOG_TARGETS)}")
        add(f"| `{target}` | " + " | ".join(cells) + " | "
            + ("yes" if sense.get("same_verdict") else "**no — unstable**") + " |")
    add("")
    add("### Residual check (model fitted on pre only)")
    add("")
    add("| target | pre median resid | post median resid | delta [CI] |")
    add("|---|---|---|---|")
    for target in targets:
        resid = block["residual"][target]
        if not resid.get("available"):
            add(f"| `{target}` | — | — | — |")
            continue
        delta = resid.get("delta") or {}
        add(f"| `{target}` | {_f(resid['pre_median_residual'], 3, True)} | "
            f"{_f(resid['post_median_residual'], 3, True)} | "
            f"{_f(delta.get('delta'), 2, signed=True)} {_ci(delta.get('ci'), 2)} |")
    add("")
    add("### Residual diagnostics (final model)")
    add("")
    add("| target | lag-1 | runs p | acf 1..4 |")
    add("|---|---|---|---|")
    for target in targets:
        diag = block["diagnostics"][target]
        if not diag.get("available"):
            add(f"| `{target}` | — | — | — |")
            continue
        acf = diag.get("acf") or {}
        add(f"| `{target}` | {_f(diag.get('lag1'), 2, True)} | {_p(diag.get('runs_p'))} | "
            + ", ".join(_f(acf.get(str(lag)), 2, True) for lag in (1, 2, 3, 4)) + " |")
    add("")
    add("### Same calendar month, post versus pre")
    add("")
    add("| target | months used | median ratio (log targets) / shift (S) |")
    add("|---|---|---|")
    for target in targets:
        month = block["month_matched"][target]
        if not month.get("available"):
            add(f"| `{target}` | — | — |")
            continue
        value = (f"{month['median_ratio']:.3f}" if month.get("median_ratio") is not None
                 else _f(month.get("median_shift_g"), 2, signed=True))
        add(f"| `{target}` | {month['n_months']} | {value} |")
    add("")


def _category(text):
    for prefix, label in (("survives", "survives"), ("co-moves", "co-moves with A"),
                          ("sign flips", "sign flip"),
                          ("not distinguishable", "no effect"),
                          ("unstable", "unstable")):
        if (text or "").startswith(prefix):
            return label
    return "n/a"


def _global_sections(add, report):
    if report.get("geometry_joint_asc"):
        _geometry_section(add, report)
    add("## The whole vector (secondary)")
    add("")
    add("Standardised residuals of `M3-` (temperature and season, no area, no "
        "post) fitted on the pre data, with a circular-shift null that keeps the "
        "serial dependence and the vector's internal correlation.")
    add("")
    add("| stratum | n | n pre/post | statistic | p (shift) | dimensions |")
    add("|---|---|---|---|---|---|")
    for stratum, block in report["strata"].items():
        multi = block.get("multivariate") or {}
        if not multi.get("available"):
            add(f"| {stratum} | — | — | — | — | — |")
            continue
        add(f"| {stratum} | {multi['n']} | {multi['n_pre']}/{multi['n_post']} | "
            f"{_f(multi['statistic'], 3)} | {_p(multi['p'])} | "
            f"{', '.join(multi['dims'])} |")
    add("")
    _age_rest(add, report)


def _geometry_section(add, report):
    add("## Geometry — the joint ASC model (family B)")
    add("")
    add("Both ASC tracks in one model: `post`, `angle` and `post x angle`, with "
        "`log A`, temperature and season. With two angle values the interaction "
        "is the stratum contrast, so this adds no information beyond the two "
        "stratified fits — it is reported as its own family, not as a headline.")
    add("")
    add("| target | post | angle | post x angle | F (interaction) | p |")
    add("|---|---|---|---|---|---|")
    for target, geometry in (report.get("geometry_joint_asc") or {}).items():
        if not geometry.get("available"):
            add(f"| `{target}` | — | — | — | — | — |")
            continue
        terms = geometry.get("terms") or {}
        def term(name):
            fit = terms.get(name)
            if not fit:
                return "—"
            return f"{_f(fit['beta'], 3, True)} {_ci(fit['ci'], 3)}"
        test = geometry.get("f_interaction") or {}
        add(f"| `{target}` | {term('post')} | {term('angle')} | {term('post_x_angle')} | "
            f"{_f(test.get('F'), 2)} | {_p(test.get('p'))} |")
    add("")


def _age_rest(add, report):
    add("## Placebo events (exploratory, family C)")
    add("")
    add("Pseudo events inside the pre period, run through the same estimator. The "
        "pseudo-post group is smaller than the real one, so the size-matched rank "
        "is the fair comparison; both are shown.")
    add("")
    add("| stratum | spec | real | rank (raw) | matched size | rank (matched) | of |")
    add("|---|---|---|---|---|---|---|")
    for stratum, block in report["strata"].items():
        panel = block.get("placebo") or {}
        if not panel.get("available"):
            add(f"| {stratum} | — | — | — | — | — | — |")
            continue
        for label, entry in panel["by_spec"].items():
            add(f"| {stratum} | {label} | {_f(entry['real']['statistic'], 3)} | "
                f"{entry['rank_of_real']} | {entry['matched_post_size']} | "
                f"{entry.get('rank_of_real_matched') or '—'} | {entry['of']} |")
    add("")
    add("A rank of 1 means the real event produced a stronger signature than every "
        "pseudo event; a rank of n+1 means all pseudo events were stronger.")
    add("")
    add("## Holm across all target x stratum post effects (sensitivity)")
    add("")
    holm_all = report.get("holm_all_strata") or {}
    add(f"Family size {holm_all.get('size')}; values in `report['holm_all_strata']`.")
    add("")
    add("| stratum | target | Holm p (all 15) |")
    add("|---|---|---|")
    for key, p in (holm_all.get("values") or {}).items():
        stratum, target = key.split("|")
        add(f"| {stratum} | `{target}` | {_p(p)} |")
    add("")
    add("## Decision logic, applied")
    add("")
    add("| stratum | survives | co-moves with A | sign flip | no effect | unstable |")
    add("|---|---|---|---|---|---|")
    for stratum, block in report["strata"].items():
        buckets = {}
        for target, text in block["verdict"].items():
            buckets.setdefault(_category(text), []).append(target)
        add(f"| {stratum} | " + " | ".join(
            f"{len(buckets.get(label, []))}" + (f" ({', '.join(buckets[label])})"
                                                if buckets.get(label) else "")
            for label in ("survives", "co-moves with A", "sign flip", "no effect",
                          "unstable")) + " |")
    add("")
    add("The verdict is about the **covariates**, not about multiplicity: a cell can "
        "survive the area/temperature/season control and still miss its Holm "
        "threshold (see the last column of the minimal table). Read the two "
        "together.")
    add("")
    add("## What this does not say")
    add("")
    add("* `A` is not a neutral confounder here. The area is part of what the "
        "collapse changes, so the A-adjusted column measures **co-movement with "
        "the observed area**, not a correction of an acquisition artefact. Where "
        "the two ladders disagree in sign, the A-adjusted estimate is the one to "
        "distrust: conditioning on a variable that the event moves can create "
        "association as easily as it removes it.")
    add("* `T24` and `T72` do not exist in this record (one temperature per "
        "overpass). `T0` is the primary temperature and `T_smooth30` the "
        "sensitivity; neither is a measurement of the deck temperature.")
    add("* `Post` is a step at a single date. These numbers describe this one "
        "collapse; they are not a causal effect of collapsing and not a forecast.")
    add("* The placebos live inside the pre period and cannot reproduce the real "
        "group sizes; the size-matched rank is reported beside the raw one.")
    add("* Descending geometry has 56 % of its post area below the pre 5th "
        "percentile: its A-adjusted estimate is an **extrapolation** along the "
        "area relation, and the common-support subset is the number to read.")
    add("")


# ---------------------------------------------------------------------------
# Output, verification, self-test
# ---------------------------------------------------------------------------
def write_outputs(report, json_path=None, md_path=None, quiet=False):
    json_path = json_path or OUT_JSON
    md_path = md_path or OUT_MD
    for path in (json_path, md_path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1)
        handle.write("\n")
    with open(md_path, "w") as handle:
        handle.write(to_markdown(report))
        handle.write("\n")
    if not quiet:
        print(f"  wrote {os.path.relpath(json_path, REPO)}")
        print(f"  wrote {os.path.relpath(md_path, REPO)}")
    return json_path, md_path


def _close(left, right, tol=1e-9):
    if left is None or right is None:
        return left is None and right is None
    try:
        return abs(float(left) - float(right)) <= tol
    except (TypeError, ValueError):
        return left == right


def verify(report=None, table=None):
    """Recompute the deterministic numbers (point estimates, n, support)."""
    table = table or TABLE
    if report is None:
        with open(OUT_JSON) as handle:
            report = json.load(handle)
    if report["table"]["sha256"] != sha256_of(table):
        print("FAILED — the table on disk is not the one this report came from")
        return 1
    rows = load_rows(table)
    failures, checked = [], 0
    inverse = {label: key for key, label
               in (report.get("stratum_map") or {}).items()}
    for label in report["strata"]:
        stratum = inverse.get(label, label)
        sub = rows_of(rows, stratum)
        block = report["strata"][label]
        if counts(sub) != block["counts"]:
            failures.append(f"{label}: counts differ")
        checked += 1
        for target in TARGETS:
            raw = block["raw"][target]
            values = [transform(target, row["values"][target]) for row in sub]
            pre = [v for v, row in zip(values, sub) if v is not None and not row["post"]]
            post = [v for v, row in zip(values, sub) if v is not None and row["post"]]
            if raw.get("available"):
                checked += 1
                shift = float(np.median(post) - np.median(pre))
                if not _close(shift, raw["median_shift_g"]):
                    failures.append(f"{stratum}/{target}: median shift")
            for ladder_name in ("A_adjusted", "A_free"):
                for label, fit in block["ladders"][target][ladder_name].items():
                    if not fit.get("available"):
                        continue
                    prepared = prepare(sub, target, fit["spec"])
                    if prepared is None:
                        failures.append(f"{stratum}/{target}/{label}: not estimable")
                        continue
                    fresh = bootstrap_post(prepared, boot=4, seed=SEED)
                    checked += 2
                    if not _close(fresh["beta"], fit["beta"]):
                        failures.append(f"{stratum}/{target}/{label}: beta "
                                        f"{fresh['beta']} != {fit['beta']}")
                    if fresh["n"] != fit["n"]:
                        failures.append(f"{stratum}/{target}/{label}: n")
    print(f"  {checked} field(s) recomputed from {os.path.relpath(table, REPO)}")
    if failures:
        print(f"FAILED — {len(failures)} mismatch(es): {failures[:5]}")
        return 1
    print("  all recomputed values match the JSON")
    return 0


def _synthetic(kind, n_pre=60, n_post=80, seed=1):
    """Synthetic overpasses where the answer is known before the fit.

    ``kind`` is 'area' (the post drop is produced by the area alone), 'season'
    (by the seasonal cycle), 'genuine' (a post step that no covariate carries) or
    'null' (nothing).
    """
    rng = np.random.default_rng(seed)
    rows = []
    start = dt.datetime(2022, 9, 1)
    total = n_pre + n_post
    for index in range(total):
        post = 1 if index >= n_pre else 0
        if kind == "season":
            when = (dt.datetime(2024, 9, 20) + dt.timedelta(days=1 * (index - n_pre))
                    if post else start + dt.timedelta(days=10 * index))
        else:
            when = start + dt.timedelta(days=10 * index)
        doy = float(when.timetuple().tm_yday)
        log_a = rng.normal(5.2, 0.35)
        if kind == "area" and post:
            log_a -= 0.9
        temperature = 10.0 + 12.0 * math.sin(2 * math.pi * (doy - 100) / 365.25)
        seasonal = math.sin(2 * math.pi * doy / 365.25)
        if kind == "area":
            log_y = 1.0 * log_a + 0.0 * post + rng.normal(0, 0.06)
        elif kind == "season":
            log_y = 0.6 * seasonal + rng.normal(0, 0.06)
        elif kind == "genuine":
            log_y = -0.4 * post + rng.normal(0, 0.06)
        else:
            log_y = rng.normal(0, 0.1)
        value = math.exp(log_y)
        rows.append({
            "t": when, "doy": doy, "month": when.month, "year": when.year,
            "stratum": "SYN", "post": post, "angle": 36.87 if index % 2 == 0 else 44.83,
            "T": temperature, "T_smooth30": temperature, "dT30": 0.0,
            "values": {"gamma2": value, "P": value, "D": value, "F": value,
                       "S": value, "A": math.exp(log_a)},
        })
    return rows


def selftest():
    """Four planted worlds, the bootstrap, the guards and the wording rule."""
    configure("carola")
    failures = []

    def check(label, condition, detail=""):
        print(f"  [{'ok  ' if condition else 'FAIL'}] {label}"
              + (f" — {detail}" if detail else ""))
        if not condition:
            failures.append(label)

    specs = specs_for_ladders()

    print("\n1. The effect is produced by the area alone")
    rows = _synthetic("area")
    ladder_rows = ladder(rows, "D", specs["A_adjusted"], boot=200)
    free = ladder(rows, "D", {"M0": {}, "M3": specs["A_free"]["M3"]}, boot=200)
    attenuation_area = attenuation(ladder_rows["M0"], ladder_rows["M1"])
    check("adding log A attenuates the post effect",
          attenuation_area is not None and attenuation_area > 0.7,
          f"attenuation {attenuation_area:.2f}")
    check("the A-adjusted CI covers no effect",
          ladder_rows["M1"]["ci"][0] <= 0 <= ladder_rows["M1"]["ci"][1],
          str([round(v, 3) for v in ladder_rows["M1"]["ci"]]))
    check("the A-free ladder still shows a large effect",
          abs(free["M3"]["beta"]) > 0.3,
          f"beta {free['M3']['beta']:.2f}")
    check("and the verdict says co-movement",
          verdict("D", ladder_rows["M3"], free["M3"]).startswith("co-moves"),
          verdict("D", ladder_rows["M3"], free["M3"]))

    print("\n2. The effect is produced by the season")
    rows = _synthetic("season")
    adjusted = ladder(rows, "D", specs["A_adjusted"], boot=200)
    check("the raw model sees a post effect", abs(adjusted["M0"]["beta"]) > 0.3,
          f"M0 beta {adjusted['M0']['beta']:.2f}")
    check("the seasonal model removes it", abs(adjusted["M3"]["beta"]) < 0.2,
          f"M3 beta {adjusted['M3']['beta']:.2f}")

    print("\n3. A genuine post effect that no covariate carries")
    rows = _synthetic("genuine")
    adjusted = ladder(rows, "D", specs["A_adjusted"], boot=200)
    free = ladder(rows, "D", {"M0": {}, "M3": specs["A_free"]["M3"]}, boot=200)
    check("the effect survives at about -0.4",
          abs(adjusted["M3"]["beta"] + 0.4) < 0.12,
          f"beta {adjusted['M3']['beta']:.3f}")
    check("its CI excludes no effect",
          adjusted["M3"]["ci"][1] < 0,
          str([round(v, 3) for v in adjusted["M3"]["ci"]]))
    check("the verdict says it survives",
          verdict("D", adjusted["M3"], free["M3"]).startswith("survives"),
          verdict("D", adjusted["M3"], free["M3"]))

    print("\n4. Nothing happened")
    rows = _synthetic("null")
    adjusted = ladder(rows, "D", specs["A_adjusted"], boot=200)
    check("no model excludes no effect",
          all(fit["ci"][0] <= 0 <= fit["ci"][1] for fit in adjusted.values()),
          str({k: [round(v, 2) for v in fit["ci"]] for k, fit in adjusted.items()}))

    print("\n5. The bootstrap")
    rows = _synthetic("genuine")
    prepared = prepare(rows, "D", specs["A_adjusted"]["M3"])
    first = bootstrap_post(prepared, boot=200, seed=7)
    again = bootstrap_post(prepared, boot=200, seed=7)
    check("same seed, same interval", first["ci"] == again["ci"], str(first["ci"]))
    index = block_index(10, 5, np.random.default_rng(1))
    check("a block resample stays inside the group",
          index.min() >= 0 and index.max() <= 9 and index.size == 10)
    groups = [np.arange(0, 10), np.arange(10, 30)]
    rng = np.random.default_rng(2)
    drawn = [np.concatenate([g[block_index(g.size, 5, rng)] for g in groups])
             for _ in range(20)]
    check("no resample mixes the two states",
          all(set(d[:10]) <= set(range(10)) and set(d[10:]) <= set(range(10, 30))
              for d in drawn))

    print("\n6. The guards")
    check("T24 cannot be requested, so it cannot be invented",
          prepare(rows, "D", {"temp": True, "temp_column": "T24"}) is None)
    check("T72 likewise",
          prepare(rows, "D", {"temp": True, "temp_column": "T72"}) is None)
    report = analyse(boot=20, placebo=True, permutation=True)
    text = to_markdown(report)
    hits = [phrase for phrase in FORBIDDEN if phrase in text.lower()]
    check("the report never uses a causal phrase", not hits, str(hits))
    check("and it does say co-moves", "co-moves with" in text)
    check("no state label is written", "state label" not in text)

    print("\n7. The estimability rule and the block lengths")
    check("n_post = 42 supports k = 6", fits(42, 6))
    check("n_post = 11 does not support k = 4", not fits(11, 4))
    check("n_post = 11 supports k = 3", fits(11, 3))
    check("n_post = 9 supports k = 2 and not k = 4",
          fits(9, 2) and not fits(9, 4))
    check("block length stops at five", block_for(81) == 5 and block_for(42) == 5,
          f"{block_for(81)}, {block_for(42)}")
    check("and shrinks for a thin post group",
          block_for(11) == 3 and block_for(9) == 3,
          f"{block_for(11)}, {block_for(9)}")
    thin = _synthetic("genuine", n_pre=60, n_post=9, seed=3)
    refused = ladder(thin, "D", {"M0": {}, "M3": specs["A_adjusted"]["M3"]},
                     boot=50)
    check("a model the post group cannot carry is refused, not fitted",
          refused["M3"].get("available") is False
          and "refused" in (refused["M3"].get("reason") or ""),
          str(refused["M3"].get("reason")))
    check("the small model still fits", refused["M0"].get("available") is True)
    check("the three trees differ in what they allow",
          TREES["full"]["A_adjusted"] != TREES["minimal"]["A_adjusted"])

    print("\n8. The robustness guards")
    thin = _synthetic("genuine", n_pre=60, n_post=9, seed=5)
    fitted = bootstrap_post(prepare(thin, "D", {"area": "linear"}), block=None,
                           boot=50, seed=7)
    check("the adaptive block length reaches the fit",
          fitted["block"] in ([5, 3], 3), str(fitted["block"]))
    padded = _synthetic("genuine", n_pre=80, n_post=20, seed=6)
    placebo = placebo_target(padded, "D", {"area": "linear"}, draws=3)
    check("the placebo is matched to the smaller group on both sides",
          bool(placebo.get("available")) and placebo["n_post"] == 20
          and placebo["draws_per_event"] == 3
          and all(event["size"] <= 20 for event in placebo["events"]),
          str({key: placebo.get(key)
               for key in ("n_post", "of", "rank_of_real", "draws_per_event")}))
    check("a real effect ranks high against its own placebos",
          placebo["rank_of_real"] <= 3,
          f"rank {placebo['rank_of_real']}/{placebo['of']}")

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0




def write_summary(reports, path=None):
    """The cross-site table and the tally -- the generalisation question."""
    path = path or os.path.join(REPO, "data", "decompose_summary.md")
    lines = []
    add = lines.append
    add("# Post-collapse decomposition across the sites")
    add("")
    add("The frozen per-candidate verdicts, after the block-length and "
        "size-matched-placebo gate, are in `data/robustness_candidates.md` and "
        "`data/analysis_freeze.md`.")
    add("")
    add("One row per site x stratum x target. `A-free` is the richest model the "
        "site could fit without the area, `A-adjusted` the same model with "
        "`log A`. Both come from the site's tier, and the tier is stated: a "
        "`thin` label marks a post group below twenty overpasses.")
    add("")
    add("| site | tier | stratum | target | n post | raw | A-free | A-adjusted "
        "| 95% CI | Holm | category |")
    add("|" + "---|" * 11)
    tally = {}
    for site, report in reports.items():
        for stratum, block in report["strata"].items():
            for target in report["targets"]:
                raw = block["raw"].get(target) or {}
                labels = (block.get("final_models") or {}).get(target) or {}
                adjusted_steps = block["ladders"][target]["A_adjusted"]
                free_steps = block["ladders"][target]["A_free"]
                adjusted = adjusted_steps.get(labels.get("A_adjusted") or "") or {}
                free = free_steps.get(labels.get("A_free") or "") or {}
                n_post = (adjusted.get("n_post")
                          or (raw.get("n") or {}).get("post")
                          or block["counts"]["post"])
                thin = " (thin)" if (n_post or 0) < 20 else ""
                cell_raw = (_ratio(math.log(raw["ratio"]), 2) if raw.get("ratio")
                            else _f(raw.get("median_shift_g"), 2, signed=True))
                add(f"| {site} | {report.get('tier')} | {stratum} | `{target}` | "
                    f"{n_post}{thin} | {cell_raw} | "
                    f"{_effect(target, free, 2, True)} | "
                    f"{_effect(target, adjusted, 2, True)} | "
                    f"{_ci(adjusted.get('ci'), 2, target in LOG_TARGETS)} | "
                    f"{_p(adjusted.get('holm_p'))} | "
                    f"{_category(block['verdict'][target])} |")
                bucket = tally.setdefault(target, {})
                bucket.setdefault(_category(block["verdict"][target]), []).append(
                    f"{site}/{stratum} (n={n_post})")
    add("")
    add("## The tally: does the Carola pattern generalise?")
    add("")
    categories = ("survives", "co-moves with A", "sign flip", "no effect", "n/a")
    add("| target | " + " | ".join(categories) + " |")
    add("|" + "---|" * (len(categories) + 1))
    for target in ("gamma2", "P", "D", "F", "S"):
        if target not in tally:
            continue
        row = [", ".join(tally[target].get(category, [])) or "—"
               for category in categories]
        add(f"| `{target}` | " + " | ".join(row) + " |")
    add("")
    add("`survives` = the A-adjusted interval excludes no effect and keeps the "
        "sign of the raw move; `co-moves with A` = the raw move is largely gone "
        "once `log A` is in the model; `no effect` = neither ladder separates "
        "from no effect. The categories are about the covariates, not about "
        "multiplicity -- the Holm column is the multiplicity answer -- and a "
        "cell marked `thin` rests on fewer than twenty post overpasses.")
    add("")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    parser = argparse.ArgumentParser(
        description=("Decompose a site's post-collapse signature into area, "
                     "geometry, temperature/season and a residual post effect"))
    parser.add_argument("--site", default="carola",
                        choices=("all",) + tuple(SITES),
                        help="which site to analyse (default: carola)")
    parser.add_argument("--table", default=None)
    parser.add_argument("--json", default=None)
    parser.add_argument("--md", default=None)
    parser.add_argument("--boot", type=int, default=BOOT,
                        help=f"moving-block bootstrap draws (default {BOOT})")
    parser.add_argument("--no-placebo", action="store_true")
    parser.add_argument("--no-permutation", action="store_true")
    parser.add_argument("--robustness", action="store_true",
                        help="the pre-registered candidate check (needs the "
                             "site JSONs from a previous --site all run)")
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()
    wanted = list(SITES) if args.site == "all" else [args.site]
    if args.robustness:
        result = robustness(boot=args.boot)
        json_path = os.path.join(REPO, "data", "robustness_candidates.json")
        with open(json_path, "w") as handle:
            json.dump(result, handle, indent=1)
            handle.write(chr(10))
        md_path = write_robustness(result)
        for cell in result["candidates"]:
            if cell.get("available"):
                print(f"  {cell['site']:8s} {cell['stratum']:9s} "
                      f"{cell['target']:7s} {cell['verdict']}")
        print(f"  wrote {os.path.relpath(json_path, REPO)}")
        print(f"  wrote {os.path.relpath(md_path, REPO)}")
        return 0
    if args.verify:
        configure(wanted[0], table=args.table)
        return verify()
    reports = {}
    for name in wanted:
        configure(name, table=args.table if len(wanted) == 1 else None)
        print(f"{LABEL} — decomposition of the post-collapse signature "
              f"(tier {TREE}, strata {STRATUM_MODE})")
        report = analyse(boot=args.boot,
                         placebo=not args.no_placebo,
                         permutation=not args.no_permutation)
        write_outputs(report,
                      args.json if len(wanted) == 1 else None,
                      args.md if len(wanted) == 1 else None,
                      quiet=args.quiet)
        reports[name] = report
        if not args.quiet:
            for stratum, block in report["strata"].items():
                print(f"  {stratum:9s} {block['counts']['pre']:3d} pre / "
                      f"{block['counts']['post']:3d} post")
                for target in TARGETS:
                    labels = block["final_models"][target]
                    adjusted = block["ladders"][target]["A_adjusted"].get(
                        labels["A_adjusted"] or "") or {}
                    free = block["ladders"][target]["A_free"].get(
                        labels["A_free"] or "") or {}
                    if not adjusted.get("available"):
                        print(f"    {target:7s} refused")
                        continue
                    print(f"    {target:7s} A-adj {adjusted['beta']:+7.3f} "
                          f"[{adjusted['ci'][0]:+.3f},{adjusted['ci'][1]:+.3f}] "
                          f"({labels['A_adjusted']})  A-free "
                          f"{free.get('beta', float('nan')):+7.3f} "
                          f"({labels['A_free']})  "
                          f"{block['verdict'][target][:34]}")
    if len(reports) > 1:
        path = write_summary(reports)
        print(f"  wrote {os.path.relpath(path, REPO)}")
    return 0



# ---------------------------------------------------------------------------
# The robustness pass over the named candidates -- and then the analysis stops
# ---------------------------------------------------------------------------
BLOCKS_ROBUST = (2, 3, 5)

#: Pre-registered candidates: (site, stratum, target, model label, why it is here).
#: ``prereg`` is a label that was fixed **before** the run and outranks the gates.
CANDIDATES = (
    ("morandi", "ALL", "gamma2", "M1",
     "correction: gamma2 was missing only through a column mapping error", None),
    ("morandi", "ALL", "D", "M1", "strongest independent candidate", None),
    ("morandi", "ALL", "F", "M1", "ladder already disagreed; not rescuable",
     "not_rescuable"),
    ("haltern", "ALL", "gamma2", "M1", "thin, no geometry adjustment", None),
    ("carola", "ASC 36.9", "gamma2", "M3",
     "robust, but the placebos say not event-specific", None),
    ("carola", "DESC", "P", "M3", "direction-consistent, never Holm-firm", None),
    ("cts", "ASC", "P", "M3", "boundary case", None),
)

GATE_ORDER = ("1 ladder (A-free and A-adjusted both exclude no effect)",
              "2 blocks (2, 3 and 5 all exclude no effect)",
              "3 placebo (size-matched rank 1 of 8; rank <= 2 counts as "
              "supporting)",
              "4 Holm (p < 0.05 within the site's family)")


def robustness(candidates=CANDIDATES, blocks=BLOCKS_ROBUST, boot=BOOT, draws=20):
    """Block sweep and size-matched placebo for each pre-registered candidate.

    Nothing here searches: the spec of each candidate is the one the site report
    already used, the block lengths are fixed, and a cell that fails a gate is
    labelled -- not re-specified.
    """
    cache = {}
    out = {"schema": SCHEMA + "-robustness", "blocks": list(blocks),
           "boot": int(boot), "draws_per_event": int(draws),
           "gate_order": list(GATE_ORDER),
           "rule": "generalisation-grade = all four gates; supporting = ladder "
                   "and block gates plus placebo rank <= 2; everything else "
                   "keeps a failure label and is not promoted",
           "candidates": []}
    for site, stratum, target, label, note, prereg in candidates:
        if site not in cache:
            configure(site)
            cache[site] = {"rows": load_rows(), "report": json.load(open(OUT_JSON))}
        else:
            configure(site)
        entry = cache[site]
        inverse = {name: key for key, name
                   in (entry["report"].get("stratum_map") or {}).items()}
        sub = rows_of(entry["rows"], inverse.get(stratum, stratum))
        specs = specs_for_ladders()
        spec = specs["A_adjusted"][label]
        prepared = prepare(sub, target, spec)
        if prepared is None:
            out["candidates"].append({"site": site, "stratum": stratum,
                                      "target": target, "model": label,
                                      "available": False, "note": note,
                                      "prereg": prereg})
            continue
        n_post = int(sum(1 for row in prepared["rows"] if row["post"]))
        by_block = {}
        for block in blocks:
            fit = bootstrap_post(prepared, block=block, boot=boot, seed=SEED)
            by_block[str(block)] = {"beta": fit["beta"], "ci": fit["ci"],
                                    "excludes_null": bool(fit["ci"][0] > 0
                                                          or fit["ci"][1] < 0)}
        placebo = placebo_target(sub, target, spec, draws=draws)
        frozen = entry["report"]["strata"][stratum]
        adjusted = frozen["ladders"][target]["A_adjusted"].get(label) or {}
        free_label = frozen["final_models"][target]["A_free"]
        free = frozen["ladders"][target]["A_free"].get(free_label) or {}
        gates = {
            "ladder": bool(excludes_null(target, adjusted)
                           and excludes_null(target, free)),
            "blocks": bool(all(value["excludes_null"]
                               for value in by_block.values())),
            "placebo_rank_1": bool(placebo.get("available")
                                   and placebo["rank_of_real"] == 1),
            "placebo_rank_le_2": bool(placebo.get("available")
                                      and placebo["rank_of_real"] <= 2),
            "holm": bool((adjusted.get("holm_p") is not None)
                         and adjusted["holm_p"] < 0.05),
        }
        # The first failing gate decides, in the pre-registered order; a label
        # fixed before the run outranks all of them.
        if prereg == "not_rescuable":
            verdict = "not rescued (pre-registered)"
        elif not gates["ladder"]:
            verdict = "not rescued (ladder gate)"
        elif not gates["blocks"]:
            verdict = "not robust (block gate)"
        elif not gates["placebo_rank_1"]:
            verdict = ("supporting (placebo rank 2)" if gates["placebo_rank_le_2"]
                       else "not event-specific (placebo gate)")
        elif not gates["holm"]:
            verdict = "candidate, not Holm-firm"
        else:
            verdict = "generalisation-grade"
        out["candidates"].append({
            "site": site, "stratum": stratum, "target": target,
            "model": label, "spec": spec, "available": True,
            "n_post": n_post, "n_pre": int(prepared["y"].size - n_post),
            "by_block": by_block, "placebo": placebo, "gates": gates,
            "holm_p": adjusted.get("holm_p"), "beta": adjusted.get("beta"),
            "ci": adjusted.get("ci"), "verdict": verdict, "note": note,
            "prereg": prereg,
        })
    return out


def write_robustness(robust, path=None):
    """The candidate table and the freeze statement."""
    path = path or os.path.join(REPO, "data", "robustness_candidates.md")
    lines = []
    add = lines.append
    add("# Robustness of the named candidates -- block lengths and size-matched "
        "placebos")
    add("")
    add("Pre-registered: the model of each candidate is the one its site report "
        "already used, the block lengths are 2, 3 and 5, and the placebo "
        "subsamples each pseudo-post to the **real** post-group size. Nothing "
        "here searches for a specification that works.")
    add("")
    add("**Gate order** (a cell must pass in this order; a failure is final):")
    add("")
    for gate in robust["gate_order"]:
        add(f"* {gate}")
    add("")
    add(f"*generalisation-grade* = all four gates. *supporting* = ladder and "
        f"block gates plus a placebo rank of at most 2. Blocks: "
        f"{', '.join(str(b) for b in robust['blocks'])}; {robust['boot']} "
        f"bootstrap draws; {robust['draws_per_event']} subsamples per pseudo "
        f"event.")
    add("")
    add("| site | stratum | target | model | n post | block 2 | block 3 | block 5 "
        "| placebo rank | Holm | gate | verdict |")
    add("|" + "---|" * 12)
    for cell in robust["candidates"]:
        if not cell.get("available"):
            add(f"| {cell['site']} | {cell['stratum']} | `{cell['target']}` | "
                f"{cell['model']} | — | — | — | — | — | — | — | not estimable |")
            continue
        target = cell["target"]
        blocks = " | ".join(
            f"{_f(value['beta'], 2, True)} {_ci(value['ci'], 2, target in LOG_TARGETS)}"
            for value in (cell["by_block"][str(b)] for b in robust["blocks"]))
        placebo = cell["placebo"]
        rank = (f"{placebo['rank_of_real']}/{placebo['of']}"
                if placebo.get("available") else "—")
        failed = [name for name, ok in cell["gates"].items() if not ok]
        add(f"| {cell['site']} | {cell['stratum']} | `{target}` | {cell['model']} | "
            f"{cell['n_post']} | {blocks} | {rank} | {_p(cell['holm_p'])} | "
            f"{', '.join(failed) or 'all pass'} | **{cell['verdict']}** |")
    add("")
    add("## What each candidate may now be called")
    add("")
    for cell in robust["candidates"]:
        if not cell.get("available"):
            continue
        add(f"* **{cell['site']} / {cell['stratum']} / `{cell['target']}`** — "
            f"{cell['verdict']} ({cell['note']}).")
    add("")
    add("A placebo rank of 1 of 8 is an ordering, not a p-value: with eight "
        "pseudo events the floor is about 0.11. It says no pseudo event produced "
        "a stronger signature, nothing more.")
    add("")
    add("## Freeze")
    add("")
    add("After this table the analysis is frozen. **Not** searched, and not to be "
        "searched later: further covariates, change points, alternative "
        "windows, alternative block lengths, further sites, or a different "
        "specification for a candidate that failed a gate. A failed cell stays "
        "in this table and in the site report with its label; it is not "
        "removed.")
    add("")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return path

if __name__ == "__main__":
    sys.exit(main())
