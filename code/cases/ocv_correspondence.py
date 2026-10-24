#!/usr/bin/env python3
"""OCV correspondence — does physical proximity imply behaviour proximity?

The LUMO and Bautzen figures compare a physical plane with an OCV plane.  This
module turns that comparison into a *test*, with one rule frozen in code:

    the campaign / visit is the independent unit -- never the single overpass.

Several SAR acquisitions inside one structural campaign improve the estimate of
(gamma2, P); they do **not** create further independent measurements of (f,
zeta).  Every inference below therefore resamples and permutes **campaigns**;
date-level resampling appears only inside a campaign and only as conditional
measurement noise (`ci_measurement`, reported separately from `ci_inference`).

The four questions, in the order they have to be answered
---------------------------------------------------------
A  Is the OCV response **repeatable**?  Two visits to the same configuration,
   distance / cosine / angle between their OCV vectors, with a paired date
   bootstrap.  Decision rule fixed in advance: ``R = d_visit / median distance
   between states``.  ``R >= 1`` means the visits differ as much as the states,
   i.e. a physical mapping is *not identifiable* -- and C and D are then only
   descriptive.
B  Does the **acquisition geometry** change the mapping?  ``dy = y_ASC - y_DESC``
   per campaign; the statistic is the share of ``dy`` variance explained by the
   state label, tested by permuting the six labels.
C  Do **physical and OCV distances correspond**?  Spearman rho between
   ``d_phys`` and ``d_OCV`` over the campaign pairs, per orbit, tested by a
   Mantel-style permutation of the campaign labels.  ``d_phys`` uses the raw
   signed relative changes (that is the physical meaning); a z-scored variant is
   reported as a **sensitivity only**, never as a way to align the planes.
D  Do configurations that are physically **similar** still differ in the OCV
   residuals?  Fit ``y = a + B x`` per orbit, then the pre-declared contrast
   ``r(DAM 3) - r(DAM 4)`` -- the two states whose frequency change is identical
   on the mode grid.  No p-value is possible at this size; the contrast is an
   effect size with a bootstrap interval.

Nested models, and why M2 is not fitted
---------------------------------------
``M0: y ~ orbit`` (2 parameters), ``M1: y ~ orbit + x`` (4 parameters, two per
orbit), ``M2: M1 + orbit x x`` (6 parameters).  With six campaigns M2 is
**saturated** -- zero residual degrees of freedom -- so it is reported as *not
identifiable* instead of being forced with a ridge penalty, which would produce
numbers without producing information.  M0 vs M1 is compared by
leave-one-campaign-out prediction over the six campaigns plus a permutation of
the x assignment; in-sample R2 is not used as a decision criterion.

Bautzen can only carry part of this, and says so
------------------------------------------------
Bautzen has three disjoint state blocks (RS, DS1, DS2), **one campaign each**,
and independent f/zeta only for RS and DS1 (the raw OpenBridgeLAB data for DS2
is not in this repository).  So: A is not testable (no repeat), C is undefined
(one physical pair), D and M1 are not identifiable (two physical points), B is
descriptive over three contrasts.  DS2 stays visible in the OCV plane with its
physical vector marked *not available* -- no trajectory is invented from state
labels.  The OCV vector is the **paired** bridge-minus-control difference per
date, z-scored against the RS distribution of the same paired difference.

Usage
-----
  python3 code/cases/ocv_correspondence.py --asset lumo
  python3 code/cases/ocv_correspondence.py --asset bautzen
  python3 code/cases/ocv_correspondence.py --asset all
  python3 code/cases/ocv_correspondence.py --selftest
  python3 code/cases/ocv_correspondence.py --build-bautzen-reference

Outputs
-------
  data/cases/ocv_correspondence_lumo.json     + figures/cases/fig_ocv_correspondence_lumo.md|.jpg
  data/cases/ocv_correspondence_bautzen.json  + figures/cases/fig_ocv_correspondence_bautzen.md
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(REPO, "data")
FIGDIR = os.path.join(REPO, "figures", "cases")

SCHEMA = "ocv-correspondence-v1"
SEED = 11
N_PERM_MAX = 20000      # exact enumeration while n! <= N_PERM_MAX
N_BOOT = 2000
ALPHA = 0.05

ORBITS = ("ASC", "DESC")
ORBIT_OF = {"ASC": "ASCENDING", "DESC": "DESCENDING"}

LUMO_CSV = os.path.join(DATA, "lumo", "support", "lumo_ocv_channels.csv")
LUMO_DYN = os.path.join(DATA, "lumo", "reference", "lumo_state_dynamics.json")
BAUTZEN_CSV = os.path.join(DATA, "bautzen", "bautzen_ocv_channels.csv")
BAUTZEN_DYN = os.path.join(DATA, "bautzen", "reference",
                           "bautzen_state_dynamics.json")
BAUTZEN_FREQ_SRC = ("/home/projects/bautzen_openlab_analysis/"
                    "bautzen_frequencies.json")
BAUTZEN_DAMP_SRC = ("/home/projects/bautzen_openlab_analysis/"
                    "bautzen_damping.json")

# LUMO keeps the two physical states apart by the H1 mode; Bautzen's primary mode
# is the one with the strongest median relative power in RS (chosen without
# looking at the damage response), everything else is a sensitivity.
BAUTZEN_PRIMARY_RULE = "largest median rel_power in RS"


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def median(values):
    return float(np.median(values))


def mad(values, centre=None):
    centre = median(values) if centre is None else centre
    return median([abs(v - centre) for v in values])


def dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def cosine(a, b):
    na, nb = math.hypot(*a), math.hypot(*b)
    if na == 0.0 or nb == 0.0:
        return None
    return (a[0] * b[0] + a[1] * b[1]) / (na * nb)


def angle_deg(a, b):
    c = cosine(a, b)
    return None if c is None else math.degrees(math.acos(max(-1.0, min(1.0, c))))


def ranks(values):
    """Average ranks, ties shared."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = rank
        i = j + 1
    return out


def spearman(x, y):
    """Spearman rho with average ranks; None when either side is constant."""
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0.0 or syy == 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def zscore(vectors):
    """Standardise both axes across the campaigns (sensitivity only)."""
    a = np.asarray(vectors, dtype=float)
    mean = a.mean(axis=0)
    sd = a.std(axis=0, ddof=1)
    sd[sd == 0.0] = 1.0
    return (a - mean) / sd


def percentile_ci(values, lo=2.5, hi=97.5):
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return [None, None]
    i = int(lo / 100.0 * (len(vals) - 1))
    j = int(hi / 100.0 * (len(vals) - 1))
    return [vals[i], vals[j]]


def assignment_space(n, seed=SEED, limit=N_PERM_MAX):
    """Every relabelling of n campaigns, or a seeded sample when n! is large."""
    total = math.factorial(n)
    if total <= limit:
        return [list(p) for p in itertools.permutations(range(n))], "exact"
    rng = random.Random(seed)
    seen = set()
    out = []
    while len(out) < limit:
        perm = list(range(n))
        rng.shuffle(perm)
        key = tuple(perm)
        if key in seen:
            continue
        seen.add(key)
        out.append(perm)
    return out, "sampled"


def holm(pvalues):
    """Holm adjustment inside one pre-declared family."""
    items = [(i, p) for i, p in enumerate(pvalues) if p is not None]
    items.sort(key=lambda kv: kv[1])
    m = len(items)
    out = [None] * len(pvalues)
    running = 0.0
    for rank, (i, p) in enumerate(items):
        value = min(1.0, (m - rank) * p)
        running = max(running, value)
        out[i] = running
    return out


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def num(value):
    text = (value or "").strip()
    return float(text) if text else None


# ---------------------------------------------------------------------------
# The unit: transform, reference, and the y vector of one campaign-orbit
# ---------------------------------------------------------------------------
TRANSFORMS = {
    "log": lambda v: math.log(v) if v > 0 else None,
    "identity": lambda v: v,
}


def y_of(values, reference):
    """z of the median channel of one campaign, against its own reference.

    ``values`` are the per-date pairs in the *transformed* scale, ``reference``
    the centre/scale of the reference campaign for the same orbit and channel.
    The median over dates is the campaign estimate -- the dates are not treated
    as replicates of (f, zeta).
    """
    out = []
    for j, channel in enumerate(("gamma2", "P")):
        col = [v[j] for v in values if v[j] is not None]
        ref = reference[channel]
        if not col or not ref["scale"]:
            out.append(None)
            continue
        out.append((median(col) - ref["centre"]) / ref["scale"])
    if any(v is None for v in out):
        return None
    return out


def reference_of(rows, transform):
    """Centre and robust scale of the reference campaign, per orbit and channel."""
    out = {}
    for orbit in ORBITS:
        block = {}
        for j, channel in enumerate(("gamma2", "P")):
            vals = [r[orbit][j] for r in rows
                    if r[orbit] is not None and r[orbit][j] is not None]
            centre = median(vals)
            block[channel] = {"n": len(vals), "centre": centre,
                              "scale": mad(vals, centre),
                              "transform": transform}
        out[orbit] = block
    return out


# ---------------------------------------------------------------------------
# Adapter 1: LUMO -- six campaigns, H1 physical vector
# ---------------------------------------------------------------------------
def _lumo_rows(path):
    rows = []
    with open(path, newline="") as handle:
        for r in csv.DictReader(handle):
            if num(r.get("coherence_gamma2")) is None:
                continue
            rows.append({
                "burst_id": r["burst_id"],
                "date": r["acquisition_date"],
                "orbit": r["orbit_direction"],
                "state": r["damage_label"],
                "gamma2": num(r.get("coherence_gamma2")),
                "P": num(r.get("P")),
            })
    rows.sort(key=lambda r: (r["date"], r["burst_id"]))
    return rows


def _runs(rows):
    runs = []
    for r in rows:
        if not runs or r["state"] != runs[-1]["state"]:
            runs.append({"state": r["state"], "rows": []})
        runs[-1]["rows"].append(r)
    for i, run in enumerate(runs, 1):
        run["run"] = i
        run["first"] = run["rows"][0]["date"]
        run["last"] = run["rows"][-1]["date"]
        run["n"] = len(run["rows"])
    return runs


def load_lumo(csv_path=LUMO_CSV, dyn_path=LUMO_DYN):
    rows = _lumo_rows(csv_path)
    runs = _runs(rows)
    with open(dyn_path, encoding="utf-8") as handle:
        dyn = json.load(handle)
    x_by_state = {}
    for block in dyn["states"]:
        h1 = {m["name"]: m for m in block["modes"]}["first_higher"]
        x_by_state[block["state"]] = [h1["df_rel"], h1["zeta_rel"]]

    transform = "log"
    healthy = []
    for r in rows:
        if r["state"] != "healthy":
            continue
        entry = {"ASC": None, "DESC": None}
        for orbit, native in ORBIT_OF.items():
            if r["orbit"] == native:
                entry[orbit] = [TRANSFORMS[transform](r["gamma2"]),
                                TRANSFORMS[transform](r["P"])
                                if r["P"] is not None else None]
        healthy.append(entry)
    reference = reference_of(healthy, transform)

    seen, subjects = {}, []
    for run in runs:
        state = run["state"]
        if state == "healthy":
            continue
        seen[state] = seen.get(state, 0) + 1
        subject = {
            "id": f"{state.replace(' ', '')} v{seen[state]}",
            "state": state,
            "state_short": state.replace(" ", ""),
            "visit": seen[state],
            "first": run["first"], "last": run["last"],
            "x": x_by_state.get(state),
            "date_values": {}, "dates": {}, "n_dates": {}, "y": {},
            "n_channel": {},
        }
        for orbit, native in ORBIT_OF.items():
            pairs, dates, n_ch = [], [], {"gamma2": 0, "P": 0}
            for r in run["rows"]:
                if r["orbit"] != native:
                    continue
                pair = [TRANSFORMS[transform](r["gamma2"]),
                        TRANSFORMS[transform](r["P"])
                        if r["P"] is not None else None]
                for j, channel in enumerate(("gamma2", "P")):
                    if pair[j] is not None:
                        n_ch[channel] += 1
                if any(v is not None for v in pair):
                    pairs.append(pair)
                    dates.append(r["date"])
            subject["date_values"][orbit] = pairs
            subject["dates"][orbit] = dates
            subject["n_dates"][orbit] = len(dates)
            subject["n_channel"][orbit] = n_ch
            subject["y"][orbit] = y_of(pairs, reference[orbit]) if pairs else None
        subjects.append(subject)

    return {
        "site": "lumo",
        "label": "LUMO lattice tower",
        "unit": "campaign (visit)",
        "unit_rule": ("one campaign = one state run; the two visits of a state "
                      "are the only replicate, the overpasses inside a campaign "
                      "are not"),
        "reference_label": "healthy",
        "transform": transform,
        "x_definition": ("H1 mode: signed df_H1/f_H1,0 and dzeta_H1/zeta_H1,0 "
                         "against the healthy block of the same campaign"),
        "y_definition": ("z of the campaign median log(gamma2) and log(P) "
                         "against the healthy median/MAD of the same orbit"),
        "reference": reference,
        "subjects": subjects,
        "n_rows": len(rows),
        "n_runs": len(runs),
        "csv": csv_path, "csv_sha256": sha256(csv_path),
        "dynamics": dyn_path, "dynamics_sha256": sha256(dyn_path),
    }


# ---------------------------------------------------------------------------
# Adapter 2: Bautzen -- one campaign per state, paired bridge-minus-control
# ---------------------------------------------------------------------------
BAUTZEN_STATES = ("RS", "DS1", "DS2")


def build_bautzen_reference(out_path=BAUTZEN_DYN, freq_path=BAUTZEN_FREQ_SRC,
                            damp_path=BAUTZEN_DAMP_SRC):
    """Committed physical reference for Bautzen (RS / DS1, DS2 documented gap)."""
    with open(freq_path, encoding="utf-8") as handle:
        freq = json.load(handle)
    with open(damp_path, encoding="utf-8") as handle:
        damp = json.load(handle)

    pairs = damp["mode_pair_maps"]["RS->DS1"]
    modes = []
    for rs_key, ds1_key in pairs.items():
        rs_mode = next(m for m in freq["states"]["RS"]["modes"]
                       if f"{m['f_hz']:.3f}" == rs_key)
        ds1_mode = next(m for m in freq["states"]["DS1"]["modes"]
                        if f"{m['f_hz']:.3f}" == ds1_key)
        rs_damp = damp["states"]["RS"]["modes"].get(rs_key, {})
        ds1_damp = damp["states"]["DS1"]["modes"].get(ds1_key, {})
        zeta_rs = rs_damp.get("zeta_median")
        zeta_ds1 = ds1_damp.get("zeta_median")
        resolvable = bool(rs_damp.get("n_resolvable")) and \
            bool(ds1_damp.get("n_resolvable"))
        modes.append({
            "pair": f"{rs_key}->{ds1_key}",
            "rs_hz": rs_mode["f_hz"],
            "state_hz": ds1_mode["f_hz"],
            "df_rel": ds1_mode["f_hz"] / rs_mode["f_hz"] - 1.0,
            "df_pct_reported": ds1_mode["f_hz"] / rs_mode["f_hz"] * 100.0 - 100.0,
            "rel_power_rs": rs_mode.get("rel_power"),
            "zeta_rs": zeta_rs,
            "zeta_state": zeta_ds1,
            "zeta_rel": (None if not zeta_rs or zeta_ds1 is None
                         else (zeta_ds1 - zeta_rs) / zeta_rs),
            "n_zeta_rs": rs_damp.get("n_resolvable"),
            "n_zeta_state": ds1_damp.get("n_resolvable"),
            "zeta_rs_values": rs_damp.get("zeta_values"),
            "zeta_state_values": ds1_damp.get("zeta_values"),
            "resolvable_both": resolvable,
        })
    usable = [m for m in modes
              if m["resolvable_both"] and m["rel_power_rs"] is not None
              and m["zeta_rel"] is not None]
    primary = max(usable, key=lambda m: m["rel_power_rs"])
    for m in modes:
        m["role"] = ("primary" if m is primary else
                     "sensitivity" if m["resolvable_both"] else "context")

    doc = {
        "schema": "bautzen-state-dynamics-v1",
        "method": ("Frequency: Welch PSD per channel, peaks in the analysis band, "
                   "mode consensus across 56 recordings per state "
                   "(`bautzen_frequencies.json`). Damping: half-power bandwidth, "
                   "nperseg 60 s, SNR gate 1.5 (`bautzen_damping.json`). One "
                   "campaign per state: RS 2025-04-29/30, DS1 2025-05-12/13, "
                   "DS2 2025-09-29/30."),
        "params": {
            "n_recordings_per_state": 56,
            "primary_rule": BAUTZEN_PRIMARY_RULE,
            "reference": "RS",
            "n_modes_matched": len(modes),
        },
        "provenance": {
            "sources": {
                "frequencies": {"path": freq_path, "sha256": sha256(freq_path)},
                "damping": {"path": damp_path, "sha256": sha256(damp_path)},
            },
            "note": ("The OpenBridgeLAB analysis lives outside this repository; "
                     "the per-recording damping values are copied verbatim so "
                     "that this JSON is self-contained."),
        },
        "prepared_by": "code/cases/ocv_correspondence.py --build-bautzen-reference",
        "states": [
            {"state": "RS", "role": "reference", "campaign": "2025-04-29/30",
             "n_recordings": 56, "physical": [0.0, 0.0], "modes": modes},
            {"state": "DS1", "role": "damaged", "campaign": "2025-05-12/13",
             "n_recordings": 56,
             "physical": [primary["df_rel"], primary["zeta_rel"]],
             "primary_pair": primary["pair"], "modes": modes},
            {"state": "DS2", "role": "no-independent-physics",
             "campaign": "2025-09-29/30", "n_recordings": None,
             "physical": None, "modes": [],
             "note": ("DS2 appears on the OCV side only. The raw OpenBridgeLAB "
                      "data for DS2 is not in this repository, so no independent "
                      "f/zeta exists and no physical vector is invented.")},
        ],
    }
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(doc, handle, indent=1)
        handle.write("\n")
    return doc


def _bautzen_rows(path):
    rows = []
    with open(path, newline="") as handle:
        for r in csv.DictReader(handle):
            native = (r["orbit_direction"] or "").strip().upper()
            rows.append({
                "state": r["state"],
                "series": r["series"],
                "bridge": "Br\u00fccke" in r["series"],
                "orbit": "ASC" if native.startswith("ASC") else "DESC",
                "date": r["acquisition_date"],
                "gamma2": num(r.get("gamma2_band_floor_corrected")),
                "P": num(r.get("P")),
            })
    return rows


def load_bautzen(csv_path=BAUTZEN_CSV, ref_path=BAUTZEN_DYN):
    rows = _bautzen_rows(csv_path)
    with open(ref_path, encoding="utf-8") as handle:
        ref = json.load(handle)
    transform = "identity"

    paired, pairing = {}, {}
    for state in BAUTZEN_STATES:
        for orbit, native in ORBIT_OF.items():
            by_date = {}
            for r in rows:
                if r["state"] != state or r["orbit"] != orbit:
                    continue
                slot = by_date.setdefault(r["date"], {})
                slot["bridge" if r["bridge"] else "control"] = r
            vals, skipped = [], 0
            for date, slot in sorted(by_date.items()):
                bridge, control = slot.get("bridge"), slot.get("control")
                if bridge is None or control is None:
                    skipped += 1
                    continue
                pair = []
                for channel in ("gamma2", "P"):
                    if bridge[channel] is None or control[channel] is None:
                        pair = None
                        break
                    pair.append(bridge[channel] - control[channel])
                if pair is None:
                    skipped += 1
                else:
                    vals.append((date, pair))
            paired[(state, orbit)] = vals
            pairing[f"{state}.{orbit}"] = {
                "n_paired_dates": len(vals),
                "n_dates_without_complete_pair": skipped,
                "n_dates_total": len(by_date),
            }

    reference = {}
    for orbit in ORBITS:
        block = {}
        for j, channel in enumerate(("gamma2", "P")):
            vals = [v[j] for _d, v in paired[("RS", orbit)]]
            centre = median(vals)
            block[channel] = {"n": len(vals), "centre": centre,
                              "scale": mad(vals, centre), "transform": transform}
        reference[orbit] = block

    physical = {s["state"]: s.get("physical") for s in ref["states"]}
    subjects = []
    for state in BAUTZEN_STATES:
        subject = {
            "id": state, "state": state, "state_short": state, "visit": 1,
            "is_reference": state == "RS",
            "x": physical.get(state),
            "date_values": {}, "dates": {}, "n_dates": {}, "y": {},
            "n_channel": {},
        }
        for orbit in ORBITS:
            vals = paired[(state, orbit)]
            subject["date_values"][orbit] = [v for _d, v in vals]
            subject["dates"][orbit] = [d for d, _v in vals]
            subject["n_dates"][orbit] = len(vals)
            subject["n_channel"][orbit] = {
                "gamma2": sum(1 for _d, v in vals if v[0] is not None),
                "P": sum(1 for _d, v in vals if v[1] is not None),
            }
            subject["y"][orbit] = (y_of([v for _d, v in vals], reference[orbit])
                                   if vals else None)
        subjects.append(subject)

    return {
        "site": "bautzen",
        "label": "Bautzen OpenLabs bridge",
        "unit": "campaign (visit)",
        "unit_rule": ("one campaign per state; the dates inside a campaign are "
                      "paired bridge-minus-control acquisitions, not replicates"),
        "reference_label": "RS",
        "transform": transform,
        "x_definition": ("primary matched mode pair: signed df/f0 and "
                         "dzeta/zeta0 against the RS night; DS2 has no "
                         "independent physics"),
        "y_definition": ("z of the median paired difference "
                         "(bridge - control) of gamma2 and P against the RS "
                         "distribution of the same paired difference, same orbit; "
                         "identity scale because the floor-corrected gamma2 is "
                         "not strictly positive"),
        "reference": reference,
        "subjects": subjects,
        "pairing": pairing,
        "n_rows": len(rows),
        "csv": csv_path, "csv_sha256": sha256(csv_path),
        "dynamics": ref_path, "dynamics_sha256": sha256(ref_path),
        "gamma2_nonpositive": sum(1 for r in rows
                                  if r["gamma2"] is not None and r["gamma2"] <= 0),
    }


ADAPTERS = {
    "lumo": lambda: load_lumo(),
    "bautzen": lambda: load_bautzen(),
}


# ---------------------------------------------------------------------------
# Resampling: dates inside a campaign (measurement), campaigns (inference)
# ---------------------------------------------------------------------------
def resample_y(subject, orbit, reference, rng):
    """Rebuild one campaign-orbit y from resampled dates (measurement noise)."""
    values = subject["date_values"][orbit]
    if not values:
        return None
    draw = [values[rng.randrange(len(values))] for _ in range(len(values))]
    return y_of(draw, reference[orbit])


def subjects_with_physics(asset, require_orbit=None):
    """Campaigns that carry a physical vector (the reference is not a subject)."""
    out = []
    for s in asset["subjects"]:
        if s.get("is_reference") or s.get("x") is None:
            continue
        if require_orbit is not None and s["y"].get(require_orbit) is None:
            continue
        out.append(s)
    return out


def by_state(subjects):
    groups = {}
    for s in subjects:
        groups.setdefault(s["state"], []).append(s)
    for group in groups.values():
        group.sort(key=lambda s: s.get("visit") or 0)
    return groups


# ---------------------------------------------------------------------------
# A - repeatability of the OCV response
# ---------------------------------------------------------------------------
def test_a(asset, n_boot=N_BOOT, seed=SEED):
    ref = asset["reference"]
    subjects = [s for s in asset["subjects"] if not s.get("is_reference")]
    groups = by_state(subjects)
    repeats = {state: g for state, g in groups.items() if len(g) >= 2}
    out = {"defined": bool(repeats), "orbits": {}}
    if not repeats:
        out["reason"] = ("no state was visited twice: every configuration has "
                         "exactly one campaign, so a repeatability estimate "
                         "would have to come from dates inside one campaign, "
                         "which are not independent visits of the tower")
        return out
    rng = random.Random(seed)
    for orbit in ORBITS:
        state_mean, rows = {}, []
        for state, group in repeats.items():
            ys = [s["y"][orbit] for s in group if s["y"][orbit]]
            if len(ys) >= 2:
                state_mean[state] = [sum(c) / len(c) for c in zip(*ys)]
        between = [dist(state_mean[a], state_mean[b])
                   for i, a in enumerate(sorted(state_mean))
                   for b in sorted(state_mean)[i + 1:]]
        for state, group in sorted(repeats.items()):
            s1, s2 = group[0], group[1]
            y1, y2 = s1["y"][orbit], s2["y"][orbit]
            if y1 is None or y2 is None:
                continue
            boots = []
            for _ in range(n_boot):
                b1 = resample_y(s1, orbit, ref, rng)
                b2 = resample_y(s2, orbit, ref, rng)
                if b1 and b2:
                    boots.append(dist(b1, b2))
            d = dist(y1, y2)
            entry = {
                "state": state,
                "visits": [s1["id"], s2["id"]],
                "n_dates": [s1["n_dates"][orbit], s2["n_dates"][orbit]],
                "y": [y1, y2],
                "distance": d,
                "cosine": cosine(y1, y2),
                "angle_deg": angle_deg(y1, y2),
                "ci_measurement": percentile_ci(boots),
                "ratio_to_between_state": (d / median(between)) if between else None,
            }
            rows.append(entry)
        out["orbits"][orbit] = {
            "pairs": rows,
            "state_mean": state_mean,
            "between_state_distances": sorted(between),
            "median_between_state": median(between) if between else None,
        }
    ratios = [r["ratio_to_between_state"] for o in ORBITS
              for r in out["orbits"][o]["pairs"]
              if r["ratio_to_between_state"] is not None]
    out["median_ratio"] = median(ratios) if ratios else None
    out["n_cells"] = len(ratios)
    out["n_cells_below_one"] = sum(1 for r in ratios if r < 1.0)
    out["identifiable_by_rule"] = (None if not ratios
                                   else all(r < 1.0 for r in ratios))
    out["verdict"] = (None if not ratios else
                      "identifiable in every cell" if out["identifiable_by_rule"]
                      else f"identifiable in {out['n_cells_below_one']} of "
                           f"{out['n_cells']} state-orbit cells")
    return out


# ---------------------------------------------------------------------------
# B - does the acquisition geometry change the mapping?
# ---------------------------------------------------------------------------
def test_b(asset, n_boot=N_BOOT, seed=SEED):
    ref = asset["reference"]
    subs = [s for s in asset["subjects"] if not s.get("is_reference")
            and all(s["y"].get(o) is not None for o in ORBITS)]
    rows = []
    rng = random.Random(seed)
    for s in subs:
        dy = [s["y"]["ASC"][j] - s["y"]["DESC"][j] for j in (0, 1)]
        boots = []
        for _ in range(n_boot):
            a = resample_y(s, "ASC", ref, rng)
            d = resample_y(s, "DESC", ref, rng)
            if a and d:
                boots.append(math.hypot(a[0] - d[0], a[1] - d[1]))
        rows.append({
            "id": s["id"], "state": s["state"], "v": dy,
            "norm": math.hypot(*dy),
            "angle_deg": math.degrees(math.atan2(dy[1], dy[0])),
            "n_dates": [s["n_dates"]["ASC"], s["n_dates"]["DESC"]],
            "ci_measurement_norm": percentile_ci(boots),
        })
    n = len(rows)
    out = {"defined": n >= 3, "contrasts": rows, "n_campaigns": n}
    if n < 3:
        out["reason"] = (f"{n} non-reference campaigns carry both geometries; "
                         "a between-state share needs at least three")
        return out
    states = sorted({r["state"] for r in rows})
    counts = {st: sum(1 for r in rows if r["state"] == st) for st in states}
    out["states"] = states
    out["counts_by_state"] = counts
    # A label explains everything when every campaign is its own state, and
    # nothing when the label is constant: both are uninformative, not a result.
    out["informative"] = len(states) >= 2 and min(counts.values()) >= 2

    grand = [sum(r["v"][j] for r in rows) / n for j in (0, 1)]

    def share(labels):
        sst = sum(sum((r["v"][j] - grand[j]) ** 2 for j in (0, 1)) for r in rows)
        if sst == 0.0:
            return None
        ssb = 0.0
        for label in sorted(set(labels)):
            idx = [i for i, lab in enumerate(labels) if lab == label]
            mean = [sum(rows[i]["v"][j] for i in idx) / len(idx) for j in (0, 1)]
            ssb += len(idx) * sum((mean[j] - grand[j]) ** 2 for j in (0, 1))
        return ssb / sst

    base = [r["state"] for r in rows]
    observed = share(base)
    perms, kind = assignment_space(n)
    null = [share([base[p] for p in perm]) for perm in perms]
    null = [v for v in null if v is not None]
    p = (sum(1 for v in null if v >= observed) / len(null)) if null else None
    out.update({
        "statistic": "between-state share of the dy variance",
        "observed": observed,
        "p_perm": p,
        "null_kind": kind,
        "null_n": len(null),
        "null_median": median(null) if null else None,
    })
    return out


# ---------------------------------------------------------------------------
# C - do physical and OCV distances correspond?
# ---------------------------------------------------------------------------
def test_c(asset, standardise=False, n_boot=N_BOOT, seed=SEED):
    subs = subjects_with_physics(asset)
    out = {"defined": False, "standardised": standardise, "orbits": {}}
    for orbit in ORBITS:
        units = [s for s in subs if s["y"].get(orbit) is not None]
        n = len(units)
        block = {"n_campaigns": n}
        if n < 3:
            block["reason"] = ("fewer than three campaigns carry both a physical "
                               "and an OCV vector, so there are no distance "
                               "matrices to compare")
            out["orbits"][orbit] = block
            continue
        xs = [list(s["x"]) for s in units]
        if standardise:
            xs = [[float(v) for v in row] for row in zscore(xs)]
        ys = [s["y"][orbit] for s in units]
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        d_phys = [dist(xs[i], xs[j]) for i, j in pairs]
        d_ocv = [dist(ys[i], ys[j]) for i, j in pairs]
        rho = spearman(d_phys, d_ocv)

        perms, kind = assignment_space(n)
        null = []
        for perm in perms:
            xp = [xs[perm[k]] for k in range(n)]
            dp = [dist(xp[i], xp[j]) for i, j in pairs]
            value = spearman(dp, d_ocv)
            if value is not None:
                null.append(value)
        p = (sum(1 for v in null if v >= rho) / len(null)) if null else None

        rng = random.Random(seed)
        boot = []
        for _ in range(n_boot):
            idx = [rng.randrange(n) for _ in range(n)]
            bx = [xs[k] for k in idx]
            by = [ys[k] for k in idx]
            dp = [dist(bx[i], bx[j]) for i, j in pairs]
            do = [dist(by[i], by[j]) for i, j in pairs]
            value = spearman(dp, do)
            if value is not None:
                boot.append(value)

        block.update({
            "rho": rho,
            "p_mantel": p,
            "null_kind": kind,
            "null_n": len(null),
            "null_median": median(null) if null else None,
            "ci_inference": percentile_ci(boot),
            "pairs": [{"ids": [units[i]["id"], units[j]["id"]],
                       "d_phys": d_phys[k], "d_ocv": d_ocv[k]}
                      for k, (i, j) in enumerate(pairs)],
            "d_phys_range": [min(d_phys), max(d_phys)],
            "d_ocv_range": [min(d_ocv), max(d_ocv)],
        })
        out["orbits"][orbit] = block
        if rho is not None:
            out["defined"] = True
    return out


# ---------------------------------------------------------------------------
# D - residuals, and the pre-declared contrast between the two closest states
# ---------------------------------------------------------------------------
def design(xs):
    return np.column_stack([np.ones(len(xs)), np.asarray(xs, dtype=float)])


def fit(xs, ys):
    """Least squares y = a + B x, returned as a 3 x 2 matrix."""
    matrix, *_ = np.linalg.lstsq(design(xs), np.asarray(ys, dtype=float),
                                 rcond=None)
    return matrix


def residuals(xs, ys):
    matrix = fit(xs, ys)
    predicted = design(xs) @ matrix
    return np.asarray(ys, dtype=float) - predicted


def closest_pair(subjects):
    """The two states whose physical vectors are closest (data-driven)."""
    states = sorted({s["state"] for s in subjects})
    best = None
    for i, a in enumerate(states):
        for b in states[i + 1:]:
            xa = next(s["x"] for s in subjects if s["state"] == a)
            xb = next(s["x"] for s in subjects if s["state"] == b)
            d = dist(xa, xb)
            if best is None or d < best[0]:
                best = (d, a, b)
    return None if best is None else {"distance": best[0], "states": [best[1],
                                                                     best[2]]}


def _state_mean_x(subjects, state):
    xs = [s["x"] for s in subjects if s["state"] == state]
    return [sum(c) / len(c) for c in zip(*xs)]


def test_d(asset, n_boot=N_BOOT, seed=SEED):
    """The closest physical pair -- and why the residual contrast cannot work.

    ``x`` is a property of the **state**, so all campaigns of a state share it.
    With three damage states there are three distinct ``x`` vectors, and
    ``y = a + B x`` has three parameters per channel: the fit therefore
    reproduces the three state means exactly, the state-level residuals are zero
    **by construction** and ``r(a) - r(b)`` is identically zero for any data.
    Reporting it as a result would be an artefact of saturation, not a finding.
    What can be measured is the *raw* OCV contrast at nearly equal ``x``,
    judged against the visit-to-visit noise of the same states.
    """
    subs = [s for s in asset["subjects"] if s["x"] is not None]
    out = {"defined": False, "identifiable": False, "orbits": {}}
    pair = closest_pair(subs)
    out["contrast_pair"] = pair
    if pair is None or len({s["state"] for s in subs}) < 2:
        out["reason"] = ("fewer than two states with a physical vector: the "
                         "contrast between physically similar states is not "
                         "defined")
        return out
    a, b = pair["states"]
    distinct_x = len({tuple(round(v, 12) for v in s["x"]) for s in subs})
    out.update({
        "n_campaigns": len(subs),
        "distinct_x_vectors": distinct_x,
        "parameters_per_channel": 3,
        "reason": (f"{distinct_x} distinct physical vectors for 3 parameters per "
                   "channel: y = a + B x is saturated at the state level, so the "
                   "residuals of a state sum to zero and the residual contrast "
                   "would be zero for any data. The raw contrast is reported "
                   "instead."),
        "physical_distance": dist(_state_mean_x(subs, a), _state_mean_x(subs, b)),
    })
    for orbit in ORBITS:
        units = [s for s in subs if s["y"].get(orbit) is not None]
        block = {"n_campaigns": len(units)}
        left = [s for s in units if s["state"] == a]
        right = [s for s in units if s["state"] == b]
        if not left or not right:
            block["reason"] = "one of the two contrast states is absent"
            out["orbits"][orbit] = block
            continue
        mean_a = list(np.asarray([s["y"][orbit] for s in left], float).mean(axis=0))
        mean_b = list(np.asarray([s["y"][orbit] for s in right], float).mean(axis=0))
        contrast = [float(mean_a[j] - mean_b[j]) for j in (0, 1)]
        rng = random.Random(seed)
        norms, vectors = [], []
        for _ in range(n_boot):
            draw_a = [resample_y(s, orbit, asset["reference"], rng) for s in left]
            draw_b = [resample_y(s, orbit, asset["reference"], rng) for s in right]
            if any(y is None for y in draw_a + draw_b):
                continue
            va = [float(np.mean([y[j] for y in draw_a])) for j in (0, 1)]
            vb = [float(np.mean([y[j] for y in draw_b])) for j in (0, 1)]
            value = [va[j] - vb[j] for j in (0, 1)]
            vectors.append(value)
            norms.append(math.hypot(*value))
        # reference for "is this contrast larger than visit noise?"
        noise = []
        for state in sorted({s["state"] for s in subs}):
            visits = [s for s in units if s["state"] == state]
            if len(visits) >= 2:
                noise.append(dist(visits[0]["y"][orbit], visits[1]["y"][orbit]))
        block.update({
            "states": [a, b],
            "y_mean": {a: mean_a, b: mean_b},
            "raw_contrast": contrast,
            "raw_contrast_norm": math.hypot(*contrast),
            "ci_measurement_norm": percentile_ci(norms),
            "ci_measurement_components": [percentile_ci([v[j] for v in vectors])
                                          for j in (0, 1)],
            "visit_noise_same_orbit": noise,
            "median_visit_noise": median(noise) if noise else None,
            "contrast_over_visit_noise": (math.hypot(*contrast) / median(noise)
                                          if noise else None),
            "residual_contrast": [0.0, 0.0],
            "residual_contrast_note": ("zero by construction: the state-level "
                                       "residuals vanish because the design is "
                                       "saturated at the state level"),
        })
        out["orbits"][orbit] = block
        out["defined"] = True
    return out


# ---------------------------------------------------------------------------
# Nested models: M0 (orbit) vs M1 (orbit + x); M2 is not identifiable
# ---------------------------------------------------------------------------
M2_PARAMS_PER_CHANNEL = 6      # intercept, orbit, 2 x, 2 interactions


def models(asset, n_boot=N_BOOT, seed=SEED):
    subs = subjects_with_physics(asset)
    out = {
        "defined": False,
        "m0": {"formula": "y ~ orbit", "params_per_channel": 2},
        "m1": {"formula": "y ~ orbit + x", "params_per_channel": 4},
        "m2": {
            "formula": "y ~ M1 + orbit x x",
            "params_per_channel": M2_PARAMS_PER_CHANNEL,
            "identifiable": False,
            "reason": ("six independent campaigns and six scalar parameters per "
                       "channel: the two orbits of a campaign share x, so the "
                       "twelve campaign-orbit cells are not twelve independent "
                       "units and M2 has no residual degrees of freedom. It is "
                       "reported, not fitted -- a ridge fit would produce "
                       "numbers without producing information."),
        },
        "orbits": {},
    }
    for orbit in ORBITS:
        units = [s for s in subs if s["y"].get(orbit) is not None]
        n = len(units)
        block = {"n_campaigns": n}
        if n < 4:
            block["reason"] = (f"{n} campaigns: comparison of M0 and M1 needs at "
                               "least four (three parameters plus one)")
            out["orbits"][orbit] = block
            continue
        xs = [s["x"] for s in units]
        ys = [s["y"][orbit] for s in units]

        def loco(x_assign, y_vals):
            err = {"M0": [], "M1": []}
            for k in range(n):
                train = [i for i in range(n) if i != k]
                mean = np.mean([y_vals[i] for i in train], axis=0)
                matrix = fit([x_assign[i] for i in train],
                             [y_vals[i] for i in train])
                err["M0"].append(np.asarray(y_vals[k], float) - mean)
                err["M1"].append(np.asarray(y_vals[k], float)
                                 - design([x_assign[k]])[0] @ matrix)

            def rmse(rows):
                stacked = np.concatenate([np.asarray(r, float) ** 2
                                          for r in rows])
                return float(np.sqrt(np.mean(stacked)))

            return rmse(err["M0"]), rmse(err["M1"])

        rmse_m0, rmse_m1 = loco(xs, ys)
        perms, kind = assignment_space(n)
        null_m1, null_delta = [], []
        for perm in perms:
            xp = [xs[perm[k]] for k in range(n)]
            m0p, m1p = loco(xp, ys)
            null_m1.append(m1p)
            null_delta.append(m0p - m1p)
        p = (sum(1 for v in null_m1 if v <= rmse_m1) / len(null_m1)) \
            if null_m1 else None

        rng = random.Random(seed)
        coefs = []
        for _ in range(n_boot):
            drawn = [resample_y(s, orbit, asset["reference"], rng) for s in units]
            if any(y is None for y in drawn):
                continue
            coefs.append(fit(xs, drawn).ravel().tolist())
        ci = ([percentile_ci([c[j] for c in coefs]) for j in range(6)]
              if coefs else [None] * 6)

        block.update({
            "loco": {"M0": rmse_m0, "M1": rmse_m1,
                     "delta": rmse_m0 - rmse_m1, "folds": n},
            "p_perm_m1": p,
            "null_kind": kind,
            "null_n": len(null_m1),
            "null_median_m1": median(null_m1) if null_m1 else None,
            "null_median_delta": median(null_delta) if null_delta else None,
            "coefficients": fit(xs, ys).ravel().tolist(),
            "coefficient_ci_measurement": ci,
            "df_residual_m1": n - 3,
            "in_sample_note": ("the table above is leave-one-campaign-out only: "
                               "in-sample R2 would be meaningless because x "
                               "takes one value per state, so M1 reproduces the "
                               "three state means exactly and only the "
                               "within-state variation is left as residual."),
        })
        out["orbits"][orbit] = block
        out["defined"] = True
    return out


# ---------------------------------------------------------------------------
# Pins
# ---------------------------------------------------------------------------
def check(checks, name, got, want):
    entry = {"check": name, "ok": bool(got == want)}
    if not entry["ok"]:
        entry["recomputed"] = repr(got)[:300]
        entry["committed"] = repr(want)[:300]
    checks.append(entry)


def pins_lumo(asset, results):
    checks = []
    subjects = asset["subjects"]
    check(checks, "lumo.n_rows", asset["n_rows"], 178)
    check(checks, "lumo.n_runs", asset["n_runs"], 13)
    check(checks, "lumo.n_campaigns", len(subjects), 6)
    check(checks, "lumo.n_states", len({s["state"] for s in subjects}), 3)
    check(checks, "lumo.n_visits_per_state",
          {st: len(g) for st, g in by_state(subjects).items()},
          {"DAM 3": 2, "DAM 4": 2, "DAM 6": 2})
    check(checks, "lumo.n_dates_by_campaign_orbit",
          {s["id"]: [s["n_dates"]["ASC"], s["n_dates"]["DESC"]]
           for s in subjects},
          {"DAM3 v1": [5, 12], "DAM3 v2": [3, 5], "DAM4 v1": [3, 5],
           "DAM4 v2": [2, 6], "DAM6 v1": [2, 6], "DAM6 v2": [2, 6]})
    check(checks, "lumo.reference_n",
          {o: {ch: asset["reference"][o][ch]["n"] for ch in ("gamma2", "P")}
           for o in ORBITS}, {"ASC": {"gamma2": 40, "P": 40},
                              "DESC": {"gamma2": 81, "P": 81}})
    check(checks, "lumo.transform", asset["transform"], "log")
    check(checks, "lumo.contrast_pair_is_dam3_dam4",
          sorted(results["d"]["contrast_pair"]["states"]), ["DAM 3", "DAM 4"])
    check(checks, "lumo.permutation_null_exact_720",
          results["c"]["orbits"]["ASC"].get("null_n"), 720)
    check(checks, "lumo.loco_folds",
          results["models"]["orbits"]["ASC"]["loco"]["folds"], 6)
    check(checks, "lumo.m2_not_identifiable",
          results["models"]["m2"]["identifiable"], False)
    check(checks, "lumo.m2_params_per_channel",
          results["models"]["m2"]["params_per_channel"], 6)
    check(checks, "lumo.a_defined", results["a"]["defined"], True)
    check(checks, "lumo.b_informative", results["b"]["informative"], True)
    check(checks, "lumo.c_defined", results["c"]["defined"], True)
    check(checks, "lumo.d_defined", results["d"]["defined"], True)
    return checks


def pins_bautzen(asset, results):
    checks = []
    subjects = {s["state"]: s for s in asset["subjects"]}
    check(checks, "bautzen.n_rows", asset["n_rows"], 142)
    check(checks, "bautzen.states", sorted(subjects), ["DS1", "DS2", "RS"])
    check(checks, "bautzen.paired_dates",
          {k: v["n_paired_dates"] for k, v in asset["pairing"].items()},
          {"RS.ASC": 6, "RS.DESC": 7, "DS1.ASC": 10, "DS1.DESC": 23,
           "DS2.ASC": 6, "DS2.DESC": 5})
    check(checks, "bautzen.incomplete_pairs",
          {k: v["n_dates_without_complete_pair"]
           for k, v in asset["pairing"].items()},
          {"RS.ASC": 1, "RS.DESC": 0, "DS1.ASC": 13, "DS1.DESC": 0,
           "DS2.ASC": 0, "DS2.DESC": 0})
    check(checks, "bautzen.ds2_has_no_physics", subjects["DS2"]["x"], None)
    check(checks, "bautzen.rs_is_origin", subjects["RS"]["x"], [0.0, 0.0])
    check(checks, "bautzen.transform", asset["transform"], "identity")
    check(checks, "bautzen.gamma2_nonpositive", asset["gamma2_nonpositive"], 72)
    check(checks, "bautzen.a_not_defined", results["a"]["defined"], False)
    check(checks, "bautzen.c_not_defined", results["c"]["defined"], False)
    check(checks, "bautzen.d_defined", results["d"]["defined"], True)
    check(checks, "bautzen.d_contrast_pair",
          sorted(results["d"]["contrast_pair"]["states"]), ["DS1", "RS"])
    check(checks, "bautzen.d_not_identifiable",
          results["d"]["identifiable"], False)
    check(checks, "bautzen.d_distinct_x", results["d"]["distinct_x_vectors"], 2)
    check(checks, "bautzen.models_not_defined",
          results["models"]["defined"], False)
    check(checks, "bautzen.b_not_defined", results["b"]["defined"], False)
    return checks


PINNERS = {"lumo": pins_lumo, "bautzen": pins_bautzen}


def descriptive(asset, n_boot=N_BOOT, seed=SEED):
    """Every OCV vector with its measurement interval (no test attached)."""
    rng = random.Random(seed)
    out = []
    for s in asset["subjects"]:
        for orbit in ORBITS:
            if s["y"].get(orbit) is None:
                continue
            draws = []
            for _ in range(n_boot):
                value = resample_y(s, orbit, asset["reference"], rng)
                if value is not None:
                    draws.append(value)
            out.append({
                "id": s["id"], "state": s["state"], "visit": s.get("visit"),
                "orbit": orbit,
                "x": s["x"],
                "y": s["y"][orbit],
                "ci_measurement": [percentile_ci([d[j] for d in draws])
                                   for j in (0, 1)],
                "n_dates": s["n_dates"][orbit],
                "is_reference": bool(s.get("is_reference")),
            })
    return out


def compute(asset_name):
    asset = ADAPTERS[asset_name]()
    results = {
        "a": test_a(asset),
        "b": test_b(asset),
        "c": test_c(asset),
        "c_standardised": test_c(asset, standardise=True),
        "d": test_d(asset),
        "models": models(asset),
        "descriptive": descriptive(asset),
    }
    family = [("B.geometry", results["b"].get("p_perm"))]
    for orbit in ORBITS:
        family.append((f"C.rho.{orbit}",
                       results["c"]["orbits"][orbit].get("p_mantel")))
    for orbit in ORBITS:
        family.append((f"M1.loco.{orbit}",
                       results["models"]["orbits"][orbit].get("p_perm_m1")))
    adjusted = holm([p for _name, p in family])
    results["multiplicity"] = {
        "family": [{"test": name, "p": p, "p_holm": adjusted[i]}
                   for i, (name, p) in enumerate(family)],
        "note": ("Holm inside this five-test family; the standardised distance "
                 "variant is a sensitivity and stays outside it."),
    }
    checks = PINNERS[asset_name](asset, results)
    n_ok = sum(1 for c in checks if c["ok"])
    return {
        "schema": SCHEMA,
        "asset": asset_name,
        "label": asset["label"],
        "unit": asset["unit"],
        "unit_rule": asset["unit_rule"],
        "reference_label": asset["reference_label"],
        "transform": asset["transform"],
        "x_definition": asset["x_definition"],
        "y_definition": asset["y_definition"],
        "reference": asset["reference"],
        "subjects": [{k: v for k, v in s.items() if k != "date_values"}
                     for s in asset["subjects"]],
        "pairing": asset.get("pairing"),
        "sources": {"csv": asset["csv"], "csv_sha256": asset["csv_sha256"],
                    "dynamics": asset["dynamics"],
                    "dynamics_sha256": asset["dynamics_sha256"]},
        "tests": results,
        "cross_checks": checks,
        "n_cross_checks": len(checks),
        "n_cross_checks_failed": len(checks) - n_ok,
    }


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------
def _f(value, digits=2, plus=False):
    if value is None:
        return "--"
    return f"{value:+.{digits}f}" if plus else f"{value:.{digits}f}"


def _v(vector, digits=2):
    if vector is None:
        return "not available"
    return "(" + ", ".join(_f(v, digits, plus=True) for v in vector) + ")"


def _table(header, rows):
    def cell(value):
        return str(value).replace("|", r"\|")

    out = ["| " + " | ".join(cell(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for row in rows:
        out.append("| " + " | ".join(cell(c) for c in row) + " |")
    return "\n".join(out)


def _head(res):
    L = []
    A = L.append
    A(f"# {res['label']} — correspondence between the physical and the OCV plane")
    A("")
    A("**Unit of inference: the campaign / visit, never the single overpass.** "
      "Several SAR acquisitions inside one structural campaign improve the "
      "estimate of (gamma2, P); they do not create additional independent "
      "measurements of (f, zeta). Every *inference* interval (`ci_inference`, "
      "the permutation p-values, the Holm adjustment) therefore resamples "
      "campaigns; the date-level bootstrap appears only as `ci_measurement` "
      "inside a campaign and never increases n.")
    A("")
    A(_table(["property", "value"], [
        ["asset", res["label"]],
        ["unit", res["unit"]],
        ["reference", res["reference_label"]],
        ["channel transform", res["transform"]],
        ["x", res["x_definition"]],
        ["y", res["y_definition"]],
        ["data", f"`{os.path.relpath(res['sources']['csv'])}` "
                 f"sha256 `{res['sources']['csv_sha256'][:10]}…`"],
        ["physics", f"`{os.path.relpath(res['sources']['dynamics'])}` "
                    f"sha256 `{res['sources']['dynamics_sha256'][:10]}…`"],
    ]))
    A("")
    if res["transform"] == "identity":
        A("`gamma2` uses the **identity** transform, not `log`: the "
          "floor-corrected Bautzen coherence takes non-positive values, so a "
          "logarithm is undefined. The statistic is identical to LUMO; only "
          "this declared per-asset transform differs.")
        A("")
    A("## 1. Campaigns")
    A("")
    rows = []
    for s in res["subjects"]:
        rows.append([s["id"], s["state"], s["visit"],
                     _v(s["x"], 4),
                     _v(s["y"].get("ASC")), _v(s["y"].get("DESC")),
                     f"{s['n_dates'].get('ASC')}/{s['n_dates'].get('DESC')}",
                     "yes" if s.get("is_reference") else ""])
    A(_table(["campaign", "state", "visit", "x (df/f0, dz/z0)", "y ASC",
              "y DESC", "n dates ASC/DESC", "reference"], rows))
    A("")
    A("The reference campaign sits at the origin of both planes by construction "
      "and is not a subject of the tests." if any(
          s.get("is_reference") for s in res["subjects"]) else
      "The healthy reference is not a campaign of its own: it is the healthy "
      "block of each campaign.")
    A("")
    if res.get("pairing"):
        A("**Pairing.** The OCV vector is the per-date difference "
          "bridge − control; a date enters only when **both channels exist in "
          "both series**. Pairing is exact (each bridge acquisition has a "
          "control acquisition on the same date), but the control chips are "
          "missing on some ASC dates, which is why the ASC counts are lower:")
        A("")
        A(_table(["cell", "dates", "complete pairs", "incomplete"],
                 [[k, v["n_dates_total"], v["n_paired_dates"],
                   v["n_dates_without_complete_pair"]]
                  for k, v in sorted(res["pairing"].items())]))
        A("")
    A("### OCV vectors with measurement intervals")
    A("")
    A("`ci_measurement` resamples the dates **inside** a campaign; it is the "
      "uncertainty of that campaign's median, not an increase in the number of "
      "campaigns.")
    A("")
    rows = []
    for item in res["tests"].get("descriptive", []):
        rows.append([item["id"], item["orbit"],
                     _v(item["y"]),
                     f"[{_f(item['ci_measurement'][0][0])}; "
                     f"{_f(item['ci_measurement'][0][1])}]",
                     f"[{_f(item['ci_measurement'][1][0])}; "
                     f"{_f(item['ci_measurement'][1][1])}]",
                     item["n_dates"],
                     _v(item["x"], 4)])
    A(_table(["campaign", "orbit", "y", "ci gamma2", "ci P", "n dates",
              "x"], rows))
    A("")
    return L


def _section_a(res):
    a = res["tests"]["a"]
    L = []
    A = L.append
    A("## 2. A — is the OCV response repeatable?")
    A("")
    if not a["defined"]:
        A(f"**Not defined.** {a.get('reason', '')}")
        A("")
        return L
    rows = []
    for orbit in ORBITS:
        for p in a["orbits"][orbit]["pairs"]:
            rows.append([orbit,
                         p["state"],
                         f"{p['n_dates'][0]}/{p['n_dates'][1]}",
                         _f(p["distance"]), _f(p["cosine"], 3),
                         _f(p["angle_deg"], 1),
                         f"[{_f(p['ci_measurement'][0])}; "
                         f"{_f(p['ci_measurement'][1])}]",
                         _f(p["ratio_to_between_state"])])
    A(_table(["orbit", "state", "n dates", "d_visit", "cos", "angle [deg]",
              "ci_measurement", "d_visit / d_between"], rows))
    A("")
    A(f"Median ratio **R = {_f(a['median_ratio'], 3)}** against the pre-declared "
      f"rule `R < 1`: the mapping is **{a['verdict']}**.")
    A("")
    bad = [(o, p["state"]) for o in ORBITS for p in a["orbits"][o]["pairs"]
           if p["cosine"] is not None and p["cosine"] < 0.9]
    A(f"Visits with a cosine below 0.9: "
      f"{', '.join(f'{s} {o}' for o, s in bad) if bad else 'none'}.")
    A("")
    return L


def _section_b(res):
    b = res["tests"]["b"]
    L = []
    A = L.append
    A("## 3. B — does the acquisition geometry change the mapping?")
    A("")
    if not b["defined"]:
        A(f"**Not defined.** {b.get('reason', '')}")
        A("")
        return L
    rows = [[c["id"], c["state"], _v(c["v"]), _f(c["norm"]),
             _f(c["angle_deg"], 1),
             f"[{_f(c['ci_measurement_norm'][0])}; "
             f"{_f(c['ci_measurement_norm'][1])}]",
             f"{c['n_dates'][0]}/{c['n_dates'][1]}"] for c in b["contrasts"]]
    A(_table(["campaign", "state", "dy = y_ASC - y_DESC", "|dy|",
              "angle [deg]", "ci_measurement |dy|", "n dates"], rows))
    A("")
    if not b["informative"]:
        A(f"**Not informative.** The statistic is {_f(b['observed'], 3)} with "
          f"{b['null_kind']} permutation p = {_f(b['p_perm'], 3)}, but every "
          f"campaign is its own state ({b['counts_by_state']}), so the label "
          f"explains the contrast by construction. With one campaign per state "
          f"this test cannot separate a state-dependent geometry effect from "
          f"the campaign itself.")
    else:
        A(f"Between-state share of the `dy` variance: **{_f(b['observed'], 3)}**, "
          f"{b['null_kind']} permutation over {b['null_n']} relabellings, "
          f"p = **{_f(b['p_perm'], 3)}** (null median "
          f"{_f(b['null_median'], 3)}).")
        A("")
        A("A large share with a small p means the ASC-DESC contrast changes with "
          "the state, i.e. the physical response interacts with the geometry. "
          "With two visits per state it cannot be separated from the visit "
          "scatter itself, which is why the per-contrast intervals are printed "
          "next to it.")
    A("")
    return L


def _section_c(res):
    c = res["tests"]["c"]
    L = []
    A = L.append
    A("## 4. C — do physical and OCV distances correspond?")
    A("")
    if not c["defined"]:
        for orbit in ORBITS:
            A(f"**{orbit}: not defined.** "
              f"{c['orbits'][orbit].get('reason', '')}")
        A("")
        return L
    rows = []
    for orbit in ORBITS:
        block = c["orbits"][orbit]
        std = res["tests"]["c_standardised"]["orbits"][orbit]
        rows.append([orbit, block["n_campaigns"],
                     _f(block["rho"], 3), _f(block["p_mantel"], 3),
                     f"[{_f(block['ci_inference'][0])}; "
                     f"{_f(block['ci_inference'][1])}]",
                     _f(std.get("rho"), 3), _f(std.get("p_mantel"), 3)])
    A(_table(["orbit", "campaigns", "rho (raw)", "p (Mantel)",
              "ci_inference (raw)", "rho (standardised)",
              "p (standardised)"], rows))
    A("")
    A("`d_phys` is the Euclidean distance in (df/f0, dz/z0): the physical "
      "meaning is kept, so the damping axis dominates by its own scale. The "
      "standardised column is the **sensitivity** that asks whether the "
      "conclusion depends on that scale; it is never used to make the planes "
      "look alike. The Mantel test permutes which physical vector belongs to "
      "which campaign, exhaustively where the campaign count allows.")
    A("")
    same = all(c["orbits"][o].get("rho") == res["tests"]["c_standardised"]["orbits"][o].get("rho")
               for o in ORBITS if c["orbits"][o].get("rho") is not None)
    if same:
        A("Here the raw and the standardised rho are **identical**: with three "
          "distinct physical vectors the rank order of the pair distances does "
          "not change when the axes are rescaled, so the conclusion does not "
          "depend on axis scaling.")
        A("")
    pairs = c["orbits"]["ASC"].get("pairs") or []
    if pairs:
        A("All campaign pairs (ASC, raw scale):")
        A("")
        A(_table(["pair", "d_phys", "d_OCV"],
                 [[f"{p['ids'][0]} - {p['ids'][1]}", _f(p["d_phys"], 4),
                   _f(p["d_ocv"])] for p in pairs]))
        A("")
    return L


def _section_d(res):
    d = res["tests"]["d"]
    L = []
    A = L.append
    A("## 5. D — the closest physical pair")
    A("")
    if not d.get("defined"):
        A(f"**Not defined.** {d.get('reason', '')}")
        A("")
        return L
    a, b = d["contrast_pair"]["states"]
    A(f"Contrast pair: **{a} vs {b}** (physical distance "
      f"{_f(d['physical_distance'], 4)}). `x` takes only "
      f"**{d['distinct_x_vectors']} distinct values** while the model has "
      f"**{d['parameters_per_channel']} parameters** per channel, so "
      f"`y = a + B x` is saturated at the state level and the residual "
      f"contrast would be zero for any data. It is reported as **not "
      f"identifiable**; the raw contrast is measured instead.")
    A("")
    rows = []
    for orbit in ORBITS:
        block = d["orbits"][orbit]
        if "raw_contrast" not in block:
            continue
        rows.append([orbit, _v(block["raw_contrast"]),
                     _f(block["raw_contrast_norm"]),
                     f"[{_f(block['ci_measurement_norm'][0])}; "
                     f"{_f(block['ci_measurement_norm'][1])}]",
                     _f(block["median_visit_noise"]),
                     _f(block["contrast_over_visit_noise"])])
    A(_table(["orbit", "raw contrast", "|contrast|", "ci_measurement",
              "median visit noise", "ratio"], rows))
    A("")
    if any(block.get("median_visit_noise") for block in d["orbits"].values()):
        A("A ratio near or below one means the OCV contrast between the two "
          "physically near-identical states is no larger than the scatter "
          "between two visits of the *same* state.")
    else:
        A("This asset has no repeated visit, so there is no within-state yardstick "
          "for the contrast: its own date-level interval is the only uncertainty "
          "it carries, and no repeatability claim can be attached to it.")
    A("")
    return L


def _section_models(res):
    m = res["tests"]["models"]
    L = []
    A = L.append
    A("## 6. Nested models")
    A("")
    A(_table(["model", "formula", "parameters / channel", "status"], [
        ["M0", m["m0"]["formula"], m["m0"]["params_per_channel"], "fitted"],
        ["M1", m["m1"]["formula"], m["m1"]["params_per_channel"],
         "fitted" if m["defined"] else "not identifiable"],
        ["M2", m["m2"]["formula"], m["m2"]["params_per_channel"],
         "**not identifiable**"],
    ]))
    A("")
    A(f"M2: {m['m2']['reason']}")
    A("")
    if not m["defined"]:
        A("No orbit carries enough campaigns for the M0/M1 comparison either.")
        A("")
        return L
    rows = []
    for orbit in ORBITS:
        block = m["orbits"][orbit]
        if "loco" not in block:
            rows.append([orbit, "--", "--", "--", "--",
                         block.get("reason", "")])
            continue
        rows.append([orbit, _f(block["loco"]["M0"], 3),
                     _f(block["loco"]["M1"], 3),
                     _f(block["loco"]["delta"], 3),
                     _f(block["p_perm_m1"], 3),
                     f"{block['loco']['folds']} folds"])
    A(_table(["orbit", "LOCO RMSE M0", "LOCO RMSE M1", "M0 - M1",
              "p (x permuted)", "folds"], rows))
    A("")
    A("`M0 - M1 > 0` means the physical vector improves out-of-campaign "
      "prediction; the p-value permutes the assignment of physical vectors to "
      "campaigns. " + str(m["orbits"]["ASC"].get("in_sample_note", "")))
    A("")
    return L


def _section_multiplicity(res):
    L = []
    A = L.append
    A("## 7. Multiplicity")
    A("")
    family = res["tests"]["multiplicity"]["family"]
    if all(f["p"] is None for f in family):
        A("**No test in this family is defined for this asset**, so there is "
          "nothing to adjust; the report above states for each test why.")
        A("")
        return L
    A(_table(["test", "p", "p (Holm)"],
             [[f["test"], _f(f["p"], 3), _f(f["p_holm"], 3)] for f in family]))
    A("")
    A(res["tests"]["multiplicity"]["note"])
    A("")
    return L


def _section_checks(res):
    L = []
    A = L.append
    A("## 8. Cross-checks")
    A("")
    A(f"**{res['n_cross_checks']} checks, "
      f"{res['n_cross_checks_failed']} failures.**")
    A("")
    if res["n_cross_checks_failed"]:
        A(_table(["check", "recomputed", "committed"],
                 [[c["check"], c.get("recomputed", ""), c.get("committed", "")]
                  for c in res["cross_checks"] if not c["ok"]]))
    else:
        groups = sorted({c["check"].split(".")[0] for c in res["cross_checks"]})
        A(_table(["group", "checks"],
                 [[g, sum(1 for c in res["cross_checks"]
                          if c["check"].split(".")[0] == g)] for g in groups]))
    A("")
    return L


def _section_caveats(res):
    L = []
    A = L.append
    A("## 9. What this does and does not show")
    A("")
    A("- **No p-value in A or D.** Three repeat pairs and one contrast cannot "
      "support a test; both are reported as effect sizes with measurement "
      "intervals.")
    A("- **Few effective physical points.** The campaign pairs are not "
      "independent observations: the visits of a state share `x`, so the "
      "information about the mapping comes from the number of distinct "
      "configurations, not from the number of pairs.")
    A("- **A permutation p answers a narrow question**: is the observed "
      "assignment of physical vectors to campaigns unusual, not whether a "
      "mechanism exists. The permutation is exact for the campaign count used.")
    A("- **A matching direction would not show that (gamma2, P) measures "
      "stiffness or damping.** It would show that two measurement channels "
      "agree about which state is which.")
    A("- **Assets are not pooled.** LUMO and Bautzen are modelled separately; "
      "only the pattern (sign, order of magnitude) is comparable, never "
      "coefficients or vector lengths.")
    A("")
    A("---")
    A("")
    A(f"Generated by `code/cases/ocv_correspondence.py --asset "
      f"{res['asset']}` ({res['n_cross_checks']} cross-checks, "
      f"{res['n_cross_checks_failed']} failures).")
    return L


def report(res):
    lines = (_head(res) + _section_a(res) + _section_b(res) + _section_c(res)
             + _section_d(res) + _section_models(res)
             + _section_multiplicity(res) + _section_checks(res)
             + _section_caveats(res))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Figure (LUMO only: Bautzen has no C/D to draw)
# ---------------------------------------------------------------------------
def make_figure(res, figdir=FIGDIR):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "mathtext.fontset": "stix",
        "font.size": 8.5, "axes.titlesize": 8.2, "axes.labelsize": 8.0,
        "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 7.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.4,
        "figure.dpi": 110, "savefig.dpi": 630, "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    })
    mm = 1.0 / 25.4
    tests = res["tests"]
    fig, axes = plt.subplots(1, 2, figsize=(183 * mm, 92 * mm))

    ax = axes[0]
    colours = {"within": "#4C72B0", "between": "#DD8452"}
    handles = []
    for orbit, marker, filled in (("ASC", "o", True), ("DESC", "s", False)):
        block = tests["c"]["orbits"][orbit]
        for kind in ("within", "between"):
            xs, ys = [], []
            for pair in block.get("pairs", []):
                same = pair["ids"][0].split(" ")[0] == pair["ids"][1].split(" ")[0]
                if (kind == "within") != same:
                    continue
                xs.append(pair["d_phys"])
                ys.append(pair["d_ocv"])
            if xs:
                ax.plot(xs, ys, linestyle="none", marker=marker, ms=5.4,
                        mfc=(colours[kind] if filled else "white"),
                        mec=colours[kind], mew=0.9, alpha=0.95)
    top = max(1.0, max(p["d_ocv"] for o in ORBITS
                       for p in tests["c"]["orbits"][o].get("pairs", [])))
    ax.plot([0, 0], [0, top * 1.05], color="0.6", lw=0.7, ls=":")
    ax.set_xlabel(r"$d_{\mathrm{phys}}$  (raw, relative changes)")
    ax.set_ylabel(r"$d_{\mathrm{OCV}}$  (z units)")
    ax.set_title("(a)  physical vs OCV distance, all campaign pairs", loc="left")
    txt = []
    for orbit in ORBITS:
        block = tests["c"]["orbits"][orbit]
        txt.append(f"{orbit}: rho {_f(block['rho'], 2)}, "
                   f"p {_f(block['p_mantel'], 3)}")
    ax.text(0.03, 0.97, "\n".join(txt), transform=ax.transAxes, va="top",
            fontsize=7.0, color="0.25",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.0))
    ax.text(0.03, 0.03, "the vertical line holds the within-state pairs\n"
                        "(x = 0): the two visits of one state",
            transform=ax.transAxes, va="bottom", fontsize=6.4, color="0.35")
    handles = [Line2D([], [], marker="o", ls="none", color=colours["within"],
                      label="within state (repeat)"),
               Line2D([], [], marker="o", ls="none", color=colours["between"],
                      label="between states"),
               Line2D([], [], marker="o", ls="none", color="0.35", label="ASC"),
               Line2D([], [], marker="s", ls="none", mfc="white",
                      mec="0.35", label="DESC")]

    ax = axes[1]
    a = tests["a"]
    d = tests["d"]
    entries = []
    for orbit in ORBITS:
        block_x = [p for p in a["orbits"][orbit]["pairs"]]
        for pair in block_x:
            entries.append({"label": pair["state"].replace(" ", ""),
                            "value": pair["distance"],
                            "ci": pair["ci_measurement"],
                            "orbit": orbit, "kind": "repeat"})
        if "raw_contrast" in d["orbits"][orbit]:
            entries.append({"label": "D3/D4",
                            "value": d["orbits"][orbit]["raw_contrast_norm"],
                            "ci": d["orbits"][orbit]["ci_measurement_norm"],
                            "orbit": orbit, "kind": "contrast"})
    gap = 0.9
    xpos, cursor = [], 0.0
    for orbit_index, orbit in enumerate(ORBITS):
        if orbit_index:
            cursor += gap
        for entry in [e for e in entries if e["orbit"] == orbit]:
            entry["x"] = cursor
            xpos.append(cursor)
            cursor += 1.0
    for kind, colour, marker in (("repeat", "0.35", "o"),
                                 ("contrast", "#C44E52", "D")):
        block = [e for e in entries if e["kind"] == kind]
        if not block:
            continue
        for orbit, filled in (("ASC", True), ("DESC", False)):
            subset = [e for e in block if e["orbit"] == orbit]
            if not subset:
                continue
            marker_used = marker if kind == "contrast" else (
                "o" if orbit == "ASC" else "s")
            ax.errorbar([e["x"] for e in subset],
                        [e["value"] for e in subset],
                        yerr=[[e["value"] - e["ci"][0] for e in subset],
                              [e["ci"][1] - e["value"] for e in subset]],
                        fmt=marker_used, ms=5.0,
                        mfc=(colour if filled else "white"), mec=colour,
                        ecolor="0.6", elinewidth=0.8, capsize=2.4, ls="none")
    ax.set_xticks(xpos)
    ax.set_xticklabels([e["label"] for e in entries], fontsize=6.3)
    for orbit, label in (("ASC", "ASC"), ("DESC", "DESC")):
        block = [e for e in entries if e["orbit"] == orbit]
        centre = sum(e["x"] for e in block) / len(block)
        ax.text(centre, -0.20, label, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.6, color="0.3")
    between = [a["orbits"][o]["median_between_state"] for o in ORBITS
               if a["orbits"][o]["median_between_state"]]
    for value in between:
        ax.axhline(value, color="#C44E52", lw=0.8, ls="--", alpha=0.6)
    ax.set_ylabel("distance in the OCV plane (z)")
    ax.set_title("(b)  visit repeats and the closest physical pair", loc="left")
    ax.text(0.02, 0.97, "o/s grey: visit repeat\nred D: DAM3-DAM4 contrast\n"
                        "dashed: median between-state distance",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.2,
            color="0.3",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.0))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 0.0), fontsize=6.8)
    fig.suptitle(
        f"{res['label']} — correspondence model, campaign as the unit\n"
        f"{len(res['subjects'])} campaigns, "
        f"{sum(1 for s in res['subjects'] if s['x'] is not None)} with a "
        "physical vector; bars = date-level bootstrap (measurement noise, not "
        "an increase in n)\n"
        "Mantel permutation over the campaign labels; the visit repeats sit at "
        "x = 0 in (a) because both visits share the same physical vector",
        fontsize=7.8, y=0.995)
    fig.tight_layout(rect=(0, 0.075, 1, 0.90))
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, f"fig_ocv_correspondence_{res['asset']}.jpg")
    fig.savefig(out, dpi=630, facecolor="white")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Self-test, including synthetic positive / negative / degenerate controls
# ---------------------------------------------------------------------------
def _synthetic_asset(visits=2, coupling=True, noise=0.05, seed=SEED):
    rng = random.Random(seed)
    states = {"A": [-0.02, -0.10], "B": [-0.06, -0.24], "C": [-0.13, -0.28]}
    reference = {o: {ch: {"n": 20, "centre": 0.0, "scale": 1.0,
                          "transform": "identity"}
                     for ch in ("gamma2", "P")} for o in ORBITS}
    subjects = []
    for state, x in states.items():
        for visit in range(1, visits + 1):
            if coupling:
                base = [2.0 * x[0] + 1.0, -3.0 * x[1] + 0.5]
            else:
                base = [rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)]
            entry = {"id": f"{state} v{visit}", "state": state,
                     "state_short": state, "visit": visit, "x": x,
                     "is_reference": False, "n_dates": {}, "dates": {},
                     "date_values": {}, "y": {}, "n_channel": {}}
            for orbit in ORBITS:
                values = [[base[0] + rng.gauss(0.0, noise),
                           base[1] + rng.gauss(0.0, noise)] for _ in range(4)]
                entry["date_values"][orbit] = values
                entry["dates"][orbit] = [f"d{i}" for i in range(4)]
                entry["n_dates"][orbit] = 4
                entry["n_channel"][orbit] = {"gamma2": 4, "P": 4}
                entry["y"][orbit] = y_of(values, reference[orbit])
            subjects.append(entry)
    return {"site": "synthetic", "label": "synthetic", "unit": "campaign",
            "unit_rule": "campaign", "reference_label": "A",
            "transform": "identity", "x_definition": "x", "y_definition": "y",
            "reference": reference, "subjects": subjects,
            "n_rows": 0, "n_runs": 0, "csv": "", "csv_sha256": "",
            "dynamics": "", "dynamics_sha256": ""}


def selftest():
    fails = []

    def check(label, condition, detail=""):
        if not condition:
            fails.append(f"{label}{': ' + detail if detail else ''}")

    check("spearman.increasing", spearman([1, 2, 3, 4], [1, 2, 3, 4]) == 1.0)
    check("spearman.decreasing", spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0)
    check("spearman.constant", spearman([1, 1, 1], [1, 2, 3]) is None)
    check("ranks.ties", ranks([10, 20, 20, 30]) == [1.0, 2.5, 2.5, 4.0],
          str(ranks([10, 20, 20, 30])))
    check("dist.3-4-5", abs(dist([0, 0], [3, 4]) - 5.0) < 1e-12)
    check("cosine.orthogonal", abs(cosine([1, 0], [0, 1])) < 1e-12)
    check("angle.45", abs(angle_deg([1, 0], [1, 1]) - 45.0) < 1e-9)
    check("holm.known", [round(v, 4) for v in holm([0.01, 0.04, 0.03])]
          == [0.03, 0.06, 0.06], str(holm([0.01, 0.04, 0.03])))
    perms, kind = assignment_space(6)
    check("assignment.exact_6", (len(perms), kind) == (720, "exact"),
          f"{len(perms)} {kind}")
    check("assignment.unique", len({tuple(p) for p in perms}) == 720)

    # Positive control: y is a linear image of x -> the correspondence must show.
    pos = test_c(_synthetic_asset(coupling=True, noise=0.02))
    rho_pos = pos["orbits"]["ASC"]["rho"]
    check("control.positive.rho", rho_pos is not None and rho_pos > 0.7,
          f"rho={rho_pos}")
    check("control.positive.p",
          pos["orbits"]["ASC"]["p_mantel"] < 0.05,
          f"p={pos['orbits']['ASC']['p_mantel']}")
    models_pos = models(_synthetic_asset(coupling=True, noise=0.02))
    check("control.positive.loco",
          models_pos["orbits"]["ASC"]["loco"]["M1"]
          < models_pos["orbits"]["ASC"]["loco"]["M0"],
          str(models_pos["orbits"]["ASC"]["loco"]))

    # Negative control: y independent of x -> no correspondence.
    neg = test_c(_synthetic_asset(coupling=False, noise=0.5))
    rho_neg = neg["orbits"]["ASC"]["rho"]
    check("control.negative.rho", rho_neg is not None and abs(rho_neg) < 0.7,
          f"rho={rho_neg}")

    # Degenerate control: one campaign per state, as in Bautzen.
    degenerate = _synthetic_asset(visits=1)
    check("control.degenerate.a_not_defined", not test_a(degenerate)["defined"])
    check("control.degenerate.d_not_identifiable",
          test_d(degenerate)["identifiable"] is False)
    check("control.degenerate.models_not_defined",
          models(degenerate)["defined"] is False)

    for item in fails:
        print(f"  FAIL {item}")
    print(f"selftest: {'ok' if not fails else str(len(fails)) + ' failures'}")
    return 0 if not fails else 1


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Physical-to-OCV correspondence with the campaign as unit.")
    ap.add_argument("--asset", default="all",
                    choices=("lumo", "bautzen", "all"))
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--build-bautzen-reference", action="store_true",
                    help="rebuild data/bautzen/reference/"
                         "bautzen_state_dynamics.json")
    args = ap.parse_args(argv)

    if args.build_bautzen_reference:
        doc = build_bautzen_reference()
        print(f"wrote {os.path.relpath(BAUTZEN_DYN)} "
              f"({len(doc['states'])} states, "
              f"{doc['params']['n_modes_matched']} matched modes, primary by "
              f"{doc['params']['primary_rule']})")
        return 0

    if args.selftest:
        return selftest()

    assets = ["lumo", "bautzen"] if args.asset == "all" else [args.asset]
    failures = 0
    for name in assets:
        res = compute(name)
        json_out = os.path.join(DATA, "cases",
                                f"ocv_correspondence_{name}.json")
        md_out = os.path.join(FIGDIR, f"fig_ocv_correspondence_{name}.md")
        os.makedirs(os.path.dirname(json_out), exist_ok=True)
        with open(json_out, "w", encoding="utf-8") as handle:
            json.dump(res, handle, indent=1, sort_keys=True)
            handle.write("\n")
        os.makedirs(os.path.dirname(md_out), exist_ok=True)
        with open(md_out, "w", encoding="utf-8") as handle:
            handle.write(report(res))
        print(f"{name}: {res['n_cross_checks']} checks, "
              f"{res['n_cross_checks_failed']} failed -> "
              f"{os.path.relpath(json_out)}, {os.path.relpath(md_out)}")
        for c in res["cross_checks"]:
            if not c["ok"]:
                print(f"    FAIL {c['check']}: {c.get('recomputed')} != "
                      f"{c.get('committed')}")
        if name == "lumo" and not args.no_figures:
            print(f"  {os.path.relpath(make_figure(res, args.figdir))}")
        failures += res["n_cross_checks_failed"]
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
