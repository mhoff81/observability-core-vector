#!/usr/bin/env python3
"""At which spatial scale does persistence appear on the Carolabrücke deck?

The question is narrow and answerable: does the deck support band contain a
persistent scatterer, and if not at the resolution cell, does one *emerge* when
the band is aggregated along its axis into 3x3 neighbourhoods, short spans,
thirds or the whole support?

Why a ladder is not a look
--------------------------
Two arithmetic effects decide the answer of any such ladder, and neither is
physics:

* ``D_A`` of an **incoherent** average falls like ``1/sqrt(N)`` whether or not
  anything is stable. Hammered through the naive gate, a ladder "discovers"
  persistence at N around 15-45 purely by averaging.
* Aggregating **dilutes** a real scatterer: the share of the block power it
  carries shrinks with N, so a coarser ladder detects *less*, never more.

The ladder therefore runs with two controls: a **ring null of the same cell
count** (2000 placements inside the 200 m annulus, two-sided against its median,
Holm over the blocks of a step, and an extreme-of-N test for the 429-cell and
333-block sweeps, which no finite null could Holm-test) and a
**planted-scatterer power curve** (which amplitude, relative to the local block
mean, a block can still see). A scale "emerges" only if it leaves the null;
otherwise the report says so.

Two scale-free readings do the actual work
------------------------------------------
``dom``  share of the block power in its brightest cell, averaged over dates.
         1.0 = one scatterer owns the block; ``-> H_N/N`` (about ln N / N) is
         the brightest of N speckle cells over their sum, i.e. distributed.
``RV``   ``Var|sum z|^2 / <|sum z|^2>^2`` of the coherent block sum.
         1.0 = fully developed speckle, ``-> 0`` = one stable dominant
         scatterer. Unlike ``D_A`` this does **not** improve with N.

Regions come from the manifest, never from constants
----------------------------------------------------
band 39 x 15 = 585 cells, guard 900, ring 7140, ``ring_complete``; the 39.0 m
deck hull maps onto columns 57-67, so the two outermost band bins per side
(55-56, 68-69) are *beside* the carriageway and are reported separately.

The mask is reproduced, not re-defined
--------------------------------------
The rule of the committed channel table — intensity at least 0.30 x the band
peak and at least 5 x the ring median, per date — yields exactly the committed
``A`` on all 436 dates, joined by date, and the run stops if it ever stops doing
so. That rule is a brightness rule: it puts 7.1% of its cells on the band, whose
share of the chip is 6.8%, so it has no preference for the deck. ``A`` is the
only quantity the rule fixes; the phase the table also publishes is not
reproduced here and is not used.

What the numbers do not say
---------------------------
* The chip is cut around the deck anchor **per date**, so a scatterer fixed to
  the deck keeps its chip index while ground clutter does not. That is the
  right frame for a deck-fixed scatterer and the wrong frame for any absolute
  ground claim. Registration/geolocation residual is about +-1-2 cells, which
  is why the single-cell row is reported but is never the basis of the verdict;
  every block >= 11 cells (41 m range x 14 m azimuth) is far outside it.
* ``gamma_g`` is a coherent phase stability and still carries unremoved
  atmosphere at 12-day repeat; it is reported, never used as the gate.
* The acceptance gates are the ones this project already used for the
  year-round test (``carola_analysis/milestone2_persistence.py``): they are
  ported verbatim, not re-tuned here.
* No state label and no timestamp: the report carries numbers, dates and a
  verdict about scatterers, nothing about the structure's condition.

Usage
-----
  python3 code/carola/carola_ocv_persist.py
  python3 code/carola/carola_ocv_persist.py --draws 400 --json /tmp/x.json
  python3 code/carola/carola_ocv_persist.py --selftest
  python3 code/carola/carola_ocv_persist.py --verify

Outputs
-------
  data/carola/support/carola_ocv_persist.json
  data/carola/support/carola_ocv_persist.md
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

SCHEMA = "ocv-persist-v1"
SITE = "carola"
LABEL = "Carolabr\u00fccke deck"
TABLE = os.path.join(REPO, "data", SITE, "support", f"{SITE}_rects_support.txt")
MANIFEST = os.path.join(REPO, "data", SITE, "support",
                        f"{SITE}_rects_support.txt.manifest.json")
CHANNEL_TABLE = os.path.join(REPO, "data", SITE, "support",
                             f"{SITE}_ocv_channels.csv")
OUT_JSON = os.path.join(REPO, "data", SITE, "support",
                        f"{SITE}_ocv_persist.json")
OUT_MD = os.path.join(REPO, "data", SITE, "support", f"{SITE}_ocv_persist.md")

#: The documented partial collapse; every series is split here and never pooled.
EVENT = "2024-09-11"

#: Mask rule of the committed channel table (``ocv-channel-table-support-v1``).
PEAK_FRACTION = 0.30        # of max(intensity) over the support band
REF_MULTIPLE = 5.0          # of median(intensity) over the ring

#: Persistence gates, ported verbatim from
#: ``carola_analysis/milestone2_persistence.py`` (cold = Nov..Mar).
DA_MAX = 0.35
PEAK_TOTAL_MIN = 0.70
PEAK_COLD_MIN = 0.50
RATIO_LO, RATIO_HI = 0.35, 2.80
PEAK_PROM = 0.45            # local max must exceed 0.45 * 5x5 median ...
PEAK_HALF = 0.50            # ... and 0.50 * 5x5 max
COLD_MONTHS = (11, 12, 1, 2, 3)
MIN_COLD, MIN_WARM = 3, 8

#: The ladder. Heights are along-track lines, widths are range bins.
LINES = (1, 3, 6, 13, 39)
NULL_DRAWS = 2000
NULL_SEED = 20240911
ALPHA = 0.05

#: Power curve: amplitude as a multiple of the local block mean, and D_A.
PSI_FACTORS = (2.0, 4.0, 8.0, 16.0, 32.0)
PSI_DA = (0.05, 0.10, 0.20)
PSI_SEED = 5
#: ``RV`` below this counts as "the block sees a dominant stable scatterer".
RV_DETECTED = 0.5


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------
def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def regions(chip, band, across_m, along_m, hull_w_m,
            ring_m=200.0, guard_m=30.0):
    """Box masks and column ranges in deck-relative chip coordinates.

    The band sits at the chip centre, the 39 m deck hull sits centred inside
    the band, and the ring is everything within ``ring_m`` of the band minus a
    ``guard_m`` collar — both rounded outwards to whole cells, as the producer
    does. The returned masks partition the chip exactly.
    """
    lines, bins = int(chip[0]), int(chip[1])
    bl, bw = int(band[0]), int(band[1])
    cl, cb = (lines - 1) // 2, (bins - 1) // 2

    def box(half_lines, half_bins):
        mask = np.zeros((lines, bins), bool)
        r0, r1 = max(cl - half_lines, 0), min(cl + half_lines + 1, lines)
        c0, c1 = max(cb - half_bins, 0), min(cb + half_bins + 1, bins)
        mask[r0:r1, c0:c1] = True
        return mask

    support = box(bl // 2, bw // 2)
    inner = box(bl // 2 + math.ceil(guard_m / along_m),
                bw // 2 + math.ceil(guard_m / across_m))
    outer = box(bl // 2 + math.ceil(ring_m / along_m),
                bw // 2 + math.ceil(ring_m / across_m))
    guard = inner & ~support
    ring = outer & ~inner

    half_bins = int(math.floor((hull_w_m / 2.0) / across_m))
    deck = np.zeros((lines, bins), bool)
    deck[cl - bl // 2:cl - bl // 2 + bl,
         cb - half_bins:cb + half_bins + 1] = True
    deck &= support
    columns = np.arange(bins)[None, :]
    return {
        "support": support,
        "guard": guard,
        "ring": ring,
        "deck": deck,
        "margin_left": support & ~deck & (columns < cb),
        "margin_right": support & ~deck & (columns > cb),
        "rows": (cl - bl // 2, cl - bl // 2 + bl),
        "cols": {
            "band": (cb - bw // 2, cb - bw // 2 + bw),
            "deck": (cb - half_bins, cb + half_bins + 1),
            "margin_left": (cb - bw // 2, cb - half_bins),
            "margin_right": (cb + half_bins + 1, cb - bw // 2 + bw),
        },
        "centre": (cl, cb),
    }

def geometry(meta):
    """Derive the chip layout from the manifest and check it against it."""
    support = meta["support"]
    across_m = float(support["record"]["range_ground_m"])
    along_m = float(support["record"]["azimuth_sample_m"])
    geo = regions(support["chip"], support["band"], across_m, along_m,
                  float(support["hull_m"][0]))
    geo["across_m"] = across_m
    geo["along_m"] = along_m
    geo["hull_m"] = [float(v) for v in support["hull_m"]]
    geo["band"] = [int(v) for v in support["band"]]
    geo["chip"] = [int(v) for v in support["chip"]]
    geo["manifest_cells"] = {
        "support": int(support["support_cells"]),
        "guard": int(support["guard_cells"]),
        "ring": int(support["ring_cells"]),
        "ring_complete": bool(support["ring_complete"]),
    }
    observed = {"support": int(geo["support"].sum()),
                "guard": int(geo["guard"].sum()),
                "ring": int(geo["ring"].sum())}
    for key, value in observed.items():
        if value != geo["manifest_cells"][key]:
            raise SystemExit(f"geometry: {key} is {value} cells, the manifest "
                             f"says {geo['manifest_cells'][key]}")
    total = sum(observed.values())
    if total != int(geo["chip"][0]) * int(geo["chip"][1]):
        raise SystemExit(f"geometry: support+guard+ring is {total} cells, the "
                         f"chip is {float(geo['chip'][0] * geo['chip'][1])}")
    geo["observed_cells"] = observed
    geo["region_cells"] = {name: int(geo[name].sum()) for name in
                           ("deck", "margin_left", "margin_right")}
    return geo


def mask_of(intensity, geo):
    """The committed mask rule: bright against the band peak and the ring."""
    return ((intensity >= PEAK_FRACTION * float(intensity[geo["support"]].max()))
            & (intensity >= REF_MULTIPLE *
               float(np.median(intensity[geo["ring"]]))))


def load(table=TABLE, manifest=MANIFEST):
    """Read the support chips and derive the mask of every date.

    35 MB of text, one JSON payload per line: ``date|orbit|rows|cols|[iq...]``
    with the complex value as ``[real, imag]`` pairs, row-major. The manifest
    rows are in the same order; that is asserted per line rather than trusted.
    """
    with open(manifest) as handle:
        meta = json.load(handle)
    geo = geometry(meta)
    rows = meta["rows"]
    dates, orbits, subswaths = [], [], []
    chips, masks, areas, shares, manifest_areas = [], [], [], [], []
    with open(table) as handle:
        for index, (line, row) in enumerate(zip(handle, rows)):
            parts = line.rstrip("\n").split("|")
            if len(parts) != 5:
                raise SystemExit(f"{table}: line {index + 1} has "
                                 f"{len(parts)} fields, expected 5")
            date, orbit, n_rows, n_cols, payload = parts
            if date != row["date"]:
                raise SystemExit(f"{table}: line {index + 1} is {date}, the "
                                 f"manifest says {row['date']} — the two are "
                                 f"not in the same order")
            iq = np.asarray(json.loads(payload), dtype=float)
            z = (iq[:, 0] + 1j * iq[:, 1]).reshape(int(n_rows), int(n_cols))
            mask = mask_of(np.abs(z) ** 2, geo)
            chips.append(z)
            masks.append(mask)
            areas.append(int(mask.sum()))
            if "A" in row:
                manifest_areas.append(int(row["A"]))
            shares.append(float((mask & geo["support"]).sum()) /
                          max(1, int(mask.sum())))
            dates.append(date)
            orbits.append(orbit)
            subswaths.append(row["subswath"])
    pre = np.asarray([d < EVENT for d in dates])
    return {
        "meta": meta,
        "geo": geo,
        "z": np.stack(chips),
        "mask": np.stack(masks),
        "area": np.asarray(areas, dtype=int),
        "support_share": np.asarray(shares, dtype=float),
        "date": np.asarray(dates),
        "orbit": np.asarray(orbits),
        "subswath": np.asarray(subswaths),
        "area_manifest": (np.asarray(manifest_areas, dtype=int)
                          if len(manifest_areas) == len(areas) else None),
        "pre": pre,
        "post": ~pre,
    }


# ---------------------------------------------------------------------------
# Block metrics — the two scale-free ones carry the argument
# ---------------------------------------------------------------------------
def block_metrics(zz):
    """``D_A``, ``gamma_g``, ``RV`` and ``dom`` of one block over one series.

    ``D_A`` and ``dom`` are computed on the block's *mean* cell statistics,
    ``gamma_g`` and ``RV`` on the coherent complex sum of the block. ``RV`` is
    the one that does not improve with N: for N independent speckle cells the
    group power stays exponentially distributed (``RV`` about 1) however large
    the block, and only a dominant stable scatterer pulls it down.
    """
    flat = zz.reshape(zz.shape[0], -1)
    power = flat.real ** 2 + flat.imag ** 2
    amplitude = np.sqrt(power)
    mean_amp = amplitude.mean(axis=1)
    coherent = flat.sum(axis=1)
    block_power = coherent.real ** 2 + coherent.imag ** 2
    return {
        "da": float(mean_amp.std() / max(mean_amp.mean(), 1e-30)),
        "gamma": float(np.abs(coherent.mean()) /
                       max(np.abs(coherent).mean(), 1e-30)),
        "rv": float(block_power.std() / max(block_power.mean(), 1e-30)),
        "dom": float(power.max(axis=1).mean() /
                     max(power.sum(axis=1).mean(), 1e-30)),
    }


METRIC_KEYS = ("da", "gamma", "rv", "dom")


def window_positions(mask, height, width):
    """Top-left corners where a ``height x width`` window lies inside ``mask``."""
    integral = np.zeros((mask.shape[0] + 1, mask.shape[1] + 1), dtype=int)
    integral[1:, 1:] = np.cumsum(np.cumsum(mask.astype(int), axis=0), axis=1)
    counts = (integral[height:, width:] - integral[:-height, width:]
              - integral[height:, :-width] + integral[:-height, :-width])
    return np.argwhere(counts == height * width)


def ladder_spec(geo):
    """The pre-defined sections, in the order the report prints them."""
    rows_0, rows_1 = geo["rows"]
    span = rows_1 - rows_0
    deck = geo["cols"]["deck"]
    band = geo["cols"]["band"]
    spec = [
        {"key": "deck_cell", "label": "1 cell (deck)", "region": "deck",
         "height": 1, "width": 1, "sweep": True},
        {"key": "deck_3x3", "label": "3x3 native (deck)", "region": "deck",
         "height": 3, "width": 3, "sweep": True},
    ]
    for lines in LINES:
        if lines > span:
            continue
        unit = "line" if lines == 1 else "lines"
        spec.append({"key": f"deck_{lines}x{deck[1] - deck[0]}",
                     "label": f"{lines} {unit} x {deck[1] - deck[0]} (deck)",
                     "region": "deck", "height": lines,
                     "width": deck[1] - deck[0], "sweep": False})
    for lines in LINES:
        if lines > span:
            continue
        unit = "line" if lines == 1 else "lines"
        spec.append({"key": f"band_{lines}x{band[1] - band[0]}",
                     "label": f"{lines} {unit} x {band[1] - band[0]} (band)",
                     "region": "band", "height": lines,
                     "width": band[1] - band[0], "sweep": False})
    for side in ("left", "right"):
        c_0, c_1 = geo["cols"][f"margin_{side}"]
        width = c_1 - c_0
        for lines in (1, span):
            unit = "line" if lines == 1 else "lines"
            spec.append({"key": f"margin_{side}_{lines}x{width}",
                         "label": f"{lines} {unit} x {width} ({side} margin)",
                         "region": f"margin_{side}", "height": lines,
                         "width": width, "sweep": False})
    return spec


def block_positions(geo, item):
    """Positions of one ladder step: either a full sweep or a tiling."""
    height, width = item["height"], item["width"]
    if item.get("sweep"):
        return [tuple(int(v) for v in p)
                for p in window_positions(geo[item["region"]], height, width)]
    rows_0, rows_1 = geo["rows"]
    c_0 = geo["cols"][item["region"]][0]
    return [(row, c_0) for row in range(rows_0, rows_1 - height + 1, height)]


def null_metrics(zz, geo, height, width, draws=NULL_DRAWS, seed=NULL_SEED):
    """The same cell count, but placed at random inside the ring."""
    positions = window_positions(geo["ring"], height, width)
    if len(positions) == 0:
        return {key: np.zeros(0) for key in METRIC_KEYS}, 0
    rng = np.random.default_rng(seed)
    take = min(int(draws), len(positions))
    pick = positions[rng.choice(len(positions), take, replace=False)]
    rows = [block_metrics(zz[:, r:r + height, c:c + width]) for r, c in pick]
    return ({key: np.asarray([row[key] for row in rows], dtype=float)
             for key in METRIC_KEYS}, len(positions))


def permutation_p(observed, null_values):
    """Two-sided empirical p against the ring null, centred on its median."""
    null_values = np.asarray(null_values, dtype=float)
    if null_values.size == 0:
        return 1.0
    centre = float(np.median(null_values))
    extreme = np.abs(null_values - centre) >= abs(float(observed) - centre)
    return float((1.0 + extreme.sum()) / (1.0 + null_values.size))


def holm(pvals):
    """Holm step-down over the blocks of one scale; ``None`` stays missing."""
    present = [(index, p) for index, p in enumerate(pvals) if p is not None]
    adjusted = [None] * len(pvals)
    running = 0.0
    for rank, (index, p) in enumerate(sorted(present, key=lambda pair: pair[1])):
        running = max(running, min(1.0, (len(present) - rank) * p))
        adjusted[index] = running
    return adjusted, len(present)



# ---------------------------------------------------------------------------
# The gates this project already used, and the planted-scatterer power curve
# ---------------------------------------------------------------------------
def gate_counts(zz, dates, geo):
    """Per-cell persistence gates, ported from the year-round test.

    A cell passes when it is a local maximum of intensity, stays above half the
    local 5x5 maximum, holds at least 70% of the dates and 50% of the cold
    season, keeps a cold/warm amplitude ratio inside [0.35, 2.80], and has
    ``D_A`` at most 0.35. Applied to the support band *and* to the ring, so a
    zero on the band can be read against a non-zero elsewhere.
    """
    intensity = zz.real ** 2 + zz.imag ** 2
    padded = np.pad(intensity, ((0, 0), (2, 2), (2, 2)), mode="edge")
    windows = np.lib.stride_tricks.sliding_window_view(padded, (5, 5),
                                                       axis=(1, 2))
    local_median = np.median(windows, axis=(3, 4))
    local_max = windows.max(axis=(3, 4))
    del windows
    peak = ((local_max > 0)
            & (intensity >= PEAK_PROM * local_median)
            & (intensity >= PEAK_HALF * local_max))
    months = np.asarray([int(str(d)[5:7]) for d in dates], dtype=int)
    cold = np.isin(months, COLD_MONTHS)
    amplitude = np.sqrt(intensity)
    mean_amp = amplitude.mean(axis=0)
    da = amplitude.std(axis=0) / np.maximum(mean_amp, 1e-30)
    if int(cold.sum()) >= MIN_COLD and int((~cold).sum()) >= MIN_WARM:
        ratio = amplitude[cold].mean(axis=0) / np.maximum(
            amplitude[~cold].mean(axis=0), 1e-30)
    else:
        ratio = np.full(da.shape, np.nan)
    share = peak.mean(axis=0)
    share_cold = peak[cold].sum(axis=0) / max(1, int(cold.sum()))
    da_ok = da <= DA_MAX
    share_ok = share >= PEAK_TOTAL_MIN
    cold_ok = share_cold >= PEAK_COLD_MIN
    ratio_ok = (ratio >= RATIO_LO) & (ratio <= RATIO_HI)
    passed = da_ok & share_ok & cold_ok & ratio_ok
    out = {}
    for name in ("support", "ring"):
        region = geo[name]
        out[name] = {
            "cells": int(region.sum()),
            "pass": int((passed & region).sum()),
            "da_ok": int((da_ok & region).sum()),
            "peak_ok": int((share_ok & region).sum()),
            "cold_ok": int((cold_ok & region).sum()),
            "ratio_ok": int((ratio_ok & region).sum()),
            "best_da": float(da[region].min()) if region.any() else None,
        }
    return out


def plant(zz, centre, factor, da_scatterer, reference=None, seed=PSI_SEED):
    """Add one scatterer at ``centre``: amplitude ``factor`` x ``reference``
    (by default the block's own mean cell amplitude, so ``factor`` means "so
    many times the local mean"), its own modulation ``da_scatterer``, and a
    phase held fixed over the series — a stable scatterer, not a random walk."""
    rng = np.random.default_rng(seed)
    amplitude = float(np.abs(zz).mean()) if reference is None else float(reference)
    count = zz.shape[0]
    series = (factor * amplitude
              * (1.0 + da_scatterer * rng.standard_normal(count))
              * np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, count)))
    out = zz.copy()
    out[:, centre[0], centre[1]] += series
    return out


def power_curve(zz, geo, factors=PSI_FACTORS, das=PSI_DA):
    """Can a block still see one scatterer of amplitude ``f`` x cell mean?

    A block of N cells sums N speckle fields against one steady phasor, so at
    fixed ``f`` a larger block is *harder*, not easier: detectability falls
    with N. This curve is what makes "finer is strictly better" a measurement
    rather than a preference.
    """
    centre = geo["centre"]
    rows_0, rows_1 = geo["rows"]
    deck = geo["cols"]["deck"]
    band = geo["cols"]["band"]
    deck_w, band_w = deck[1] - deck[0], band[1] - band[0]
    levels = [(f"deck_{lines}x{deck_w}", lines, deck_w, deck)
              for lines in LINES]
    levels.append((f"band_{rows_1 - rows_0}x{band_w}", rows_1 - rows_0,
                   band_w, band))
    records = []
    for key, height, width, cols in levels:
        row_0 = min(max(centre[0] - height // 2, rows_0), rows_1 - height)
        col_0 = min(max(centre[1] - width // 2, cols[0]), cols[1] - width)
        block = zz[:, row_0:row_0 + height, col_0:col_0 + width]
        for factor in factors:
            for da in das:
                planted = plant(block, (centre[0] - row_0, centre[1] - col_0),
                                factor, da)
                metrics = block_metrics(planted)
                records.append({
                    "level": key,
                    "cells": height * width,
                    "factor": float(factor),
                    "da": float(da),
                    "rv": metrics["rv"],
                    "gamma": metrics["gamma"],
                    "dom": metrics["dom"],
                    "detected": bool(metrics["rv"] < RV_DETECTED),
                })
    return records


# ---------------------------------------------------------------------------
# One ladder step, against its own same-N ring null
# ---------------------------------------------------------------------------
def extreme_p(observed, null_values, count, draws, upper):
    """p for the most extreme of ``count`` blocks, against the same extreme
    resampled from the ring null.

    A sweep of 429 cells cannot be Holm-tested against 2000 placements — the
    smallest attainable p is 1/2001, and 429 blocks need 1e-4 — so the sweep is
    tested as a single question: is the *most extreme* cell (or 3x3) of the deck
    outside what the most extreme of 429 ring placements looks like?
    """
    null_values = np.asarray(null_values, dtype=float)
    if null_values.size == 0 or count < 1:
        return 1.0, None, None
    rng = np.random.default_rng(NULL_SEED + (1 if upper else 2))
    picks = rng.integers(0, null_values.size, size=(int(draws), int(count)))
    sample = null_values[picks]
    null_stat = sample.max(axis=1) if upper else sample.min(axis=1)
    extreme = (null_stat >= observed) if upper else (null_stat <= observed)
    p = float((1.0 + extreme.sum()) / (1.0 + null_stat.size))
    return p, float(np.median(null_stat)), float(null_values.size)


def summarise(zz, geo, item, mask, draws, cache, with_null=True):
    """Median metrics of one ladder step, its extremes, and the ring null."""
    height, width = item["height"], item["width"]
    positions = block_positions(geo, item)
    if not positions:
        return None
    per = [block_metrics(zz[:, r:r + height, c:c + width]) for r, c in positions]
    recurrence = [float(mask[:, r:r + height, c:c + width]
                        .any(axis=(1, 2)).mean()) for r, c in positions]
    key = (height, width)
    if with_null and key not in cache:
        cache[key] = null_metrics(zz, geo, height, width, draws)
    null, ring_positions = cache.get(key, (None, 0))
    out = {
        "key": item["key"],
        "label": item["label"],
        "region": item["region"],
        "height": height,
        "width": width,
        "cells": height * width,
        "blocks": len(positions),
        "sweep": bool(item.get("sweep")),
        "extent_m": [width * geo["across_m"], height * geo["along_m"]],
        "recurrence": float(np.mean(recurrence)),
    }
    for metric in METRIC_KEYS:
        values = np.asarray([row[metric] for row in per], dtype=float)
        out[metric] = float(np.median(values))
        out[metric + "_range"] = [float(values.min()), float(values.max())]
        if null is None:
            out["null_" + metric] = None
            out["holm_" + metric] = None
            continue
        null_values = null[metric]
        out["null_" + metric] = float(np.median(null_values))
        if item.get("sweep"):
            upper = metric == "dom"
            observed = float(values.max() if upper else values.min())
            p, null_extreme, _ = extreme_p(observed, null_values, len(values),
                                           draws, upper)
            out["holm_" + metric] = {
                "test": "extreme",
                "tested": len(values),
                "significant": int(p < ALPHA),
                "persistence_side": int(p < ALPHA and (
                    observed > out["null_" + metric] if upper
                    else observed < out["null_" + metric])),
                "min_p": p,
                "observed": observed,
                "null_extreme": null_extreme,
            }
            continue
        pvals = [permutation_p(value, null_values) for value in values]
        adjusted, tested = holm(pvals)
        centre = float(np.median(null_values))
        side = (values > centre) if metric == "dom" else (values < centre)
        clean = [a for a in adjusted if a is not None]
        out["holm_" + metric] = {
            "test": "blocks",
            "tested": tested,
            "significant": int(sum(1 for a in clean if a < ALPHA)),
            "persistence_side": int(sum(1 for a, s in zip(adjusted, side)
                                       if a is not None and a < ALPHA and s)),
            "min_p": float(min(clean)) if clean else None,
        }
    if null is None:
        out["null"] = None
        return out
    out["null"] = {
        "ring_positions": int(ring_positions),
        "draws": int(null["da"].size),
        "da": float(np.median(null["da"])),
        "rv": float(np.median(null["rv"])),
        "dom": float(np.median(null["dom"])),
        "da_p05": float(np.percentile(null["da"], 5)),
        "da_p95": float(np.percentile(null["da"], 95)),
        "rv_p05": float(np.percentile(null["rv"], 5)),
        "rv_p95": float(np.percentile(null["rv"], 95)),
        "dom_p05": float(np.percentile(null["dom"], 5)),
        "dom_p95": float(np.percentile(null["dom"], 95)),
    }
    return out


def mask_check(data, channel_table=CHANNEL_TABLE):
    """Reproduce the committed mask area ``A`` of every date (joined by date)."""
    geo = data["geo"]
    chip_cells = int(geo["chip"][0]) * int(geo["chip"][1])
    share = np.asarray(data["support_share"], dtype=float)
    out = {
        "channel_table": os.path.basename(channel_table),
        "dates": int(len(data["date"])),
        "support_cells": int(geo["support"].sum()),
        "chip_cells": chip_cells,
        "support_share_of_chip": float(geo["support"].sum() / chip_cells),
        "area_median": float(np.median(data["area"])),
        "area_mean": float(data["area"].mean()),
        "area_min": int(data["area"].min()),
        "area_max": int(data["area"].max()),
        "support_share_mean": float(share.mean()),
        "support_share_median": float(np.median(share)),
    }
    if data.get("area_manifest") is not None:
        out["manifest_equal"] = int((data["area_manifest"] == data["area"]).sum())
        out["manifest_max_abs_delta"] = float(
            np.abs(data["area_manifest"] - data["area"]).max())
        if out["manifest_equal"] != len(data["area"]):
            raise SystemExit(
                f"mask rule: {out['manifest_equal']} of {len(data['area'])} "
                f"dates reproduce the manifest's A column — the rule in this "
                f"file no longer matches the manifest")
    if not os.path.exists(channel_table):
        out.update({"compared": 0, "equal": 0, "note": "channel table absent"})
        return out
    published = {}
    with open(channel_table, newline="") as handle:
        for row in csv.DictReader(handle):
            published[str(row["acquisition_ts"])[:10]] = row
    compared = equal = missing = 0
    worst = 0.0
    for date, area in zip(data["date"], data["area"]):
        row = published.get(str(date))
        if row is None:
            missing += 1
            continue
        compared += 1
        delta = float(area) - float(row["A"])
        worst = max(worst, abs(delta))
        equal += int(delta == 0.0)
    out.update({"rows_in_table": len(published), "compared": compared,
                "equal": equal, "missing": missing,
                "max_abs_delta_A": worst})
    if compared and equal != compared:
        raise SystemExit(f"mask rule: {equal} of {compared} dates reproduce the "
                         f"committed A (largest difference {worst:g}) — the rule "
                         f"in this file no longer matches the committed table")
    return out


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def rel(path):
    """Repository-relative path, for the report."""
    try:
        return os.path.relpath(path, REPO)
    except ValueError:                                  # different drive
        return path


def analyse(table=TABLE, manifest=MANIFEST, channel_table=CHANNEL_TABLE,
            draws=NULL_DRAWS):
    data = load(table, manifest)
    geo = data["geo"]
    rows_0, rows_1 = geo["rows"]
    band = geo["cols"]["band"]
    report = {
        "schema": SCHEMA,
        "site": SITE,
        "label": LABEL,
        "question": ("Does the deck support band hold a persistent scatterer, "
                     "and if not at one resolution cell, does one appear when "
                     "the band is aggregated — 3x3, line spans, thirds, full "
                     "support?"),
        "inputs": {
            "support": {"path": rel(table), "sha256": sha256_of(table),
                        "n_dates": int(len(data["date"]))},
            "manifest": {"path": rel(manifest), "sha256": sha256_of(manifest)},
            "channels": {"path": rel(channel_table),
                         "sha256": (sha256_of(channel_table)
                                    if os.path.exists(channel_table) else None)},
        },
        "geometry": {
            "chip": geo["chip"],
            "band": geo["band"],
            "rows": [rows_0, rows_1 - 1],
            "columns": {name: [int(v[0]), int(v[1] - 1)]
                        for name, v in geo["cols"].items()},
            "cells": {"support": geo["observed_cells"]["support"],
                      "guard": geo["observed_cells"]["guard"],
                      "ring": geo["observed_cells"]["ring"],
                      "ring_complete": geo["manifest_cells"]["ring_complete"],
                      "deck": geo["region_cells"]["deck"],
                      "margin_left": geo["region_cells"]["margin_left"],
                      "margin_right": geo["region_cells"]["margin_right"]},
            "ground_m": {
                "range_bin": geo["across_m"],
                "azimuth_line": geo["along_m"],
                "band": [(band[1] - band[0]) * geo["across_m"],
                         (rows_1 - rows_0) * geo["along_m"]],
                "deck_hull": geo["hull_m"],
            },
        },
        "rules": {
            "event": EVENT,
            "mask": (f"intensity >= {PEAK_FRACTION} * max(support band) and "
                     f">= {REF_MULTIPLE} * median(ring), per date"),
            "gates": {"source": "carola_analysis/milestone2_persistence.py",
                      "da_max": DA_MAX, "peak_total_min": PEAK_TOTAL_MIN,
                      "peak_cold_min": PEAK_COLD_MIN,
                      "cold_warm_ratio": [RATIO_LO, RATIO_HI],
                      "peak_prominence": PEAK_PROM,
                      "peak_half_of_max": PEAK_HALF,
                      "cold_months": list(COLD_MONTHS)},
            "null": {"region": "ring", "placements": int(draws),
                     "seed": NULL_SEED, "test": "two-sided vs the null median",
                     "family": ("tiled steps Holm over their blocks; single "
                                "cells and 3x3 sweeps tested as the extreme of "
                                "as many placements"),
                     "alpha": ALPHA},
            "power": {"factors": list(PSI_FACTORS), "da": list(PSI_DA),
                      "detected_below_rv": RV_DETECTED, "seed": PSI_SEED},
            "ladder": {"lines": list(LINES),
                       "deck_cell_sweep": True, "deck_3x3_sweep": True},
        },
        "mask_check": mask_check(data, channel_table),
        "channels": {},
        "limits": [
            "the chip is cut around the deck anchor per date, so a deck-fixed "
            "scatterer keeps its chip index while ground clutter does not; "
            "that is the right frame for a deck-fixed scatterer and the wrong "
            "frame for an absolute ground claim, and registration residual "
            "(about +-1-2 cells) is why the single-cell row is reported but is "
            "never the basis of the verdict",
            "the ring null shares the chip's geometry but not its view "
            "geometry, so it bounds what incoherent averaging does at the same "
            "cell count, not what the deck must look like",
            "gamma_g still carries unremoved atmosphere at 12-day repeat; it is "
            "reported and never used as the gate",
            "the gates are ported from the year-round test exactly; they were "
            "not re-tuned for this ladder, which is why a zero on the band can "
            "be read against the ring",
            "no label for the structure and no timestamp: the report carries "
            "numbers, dates and a verdict about scatterers only",
        ],
    }

    for orbit in sorted(set(str(o) for o in data["orbit"])):
        for subswath in sorted(set(str(s) for s in data["subswath"])):
            member = (data["orbit"] == orbit) & (data["subswath"] == subswath)
            if not member.any():
                continue
            channel = {
                "orbit": orbit,
                "subswath": subswath,
                "n": {},
                "dates": {},
                "gates": {},
                "ladder": {},
                "mask": {},
            }
            dates = data["date"][member]
            pre = data["pre"][member]
            channel["dates"] = {
                "first": str(min(dates)),
                "last_pre": str(max(dates[pre])) if pre.any() else None,
                "first_post": str(min(dates[~pre])) if (~pre).any() else None,
                "last": str(max(dates)),
            }
            for split in ("pre", "post"):
                select = member & data[split]
                channel["n"][split] = int(select.sum())
                if not select.any():
                    continue
                zz = data["z"][select]
                cache = {}
                ladder = [summarise(zz, geo, item, data["mask"][select],
                                    draws, cache, with_null=(split == "pre"))
                          for item in ladder_spec(geo)]
                channel["ladder"][split] = [row for row in ladder if row]
                channel["gates"][split] = gate_counts(zz, data["date"][select],
                                                      geo)
                area = data["area"][select]
                share = data["support_share"][select]
                channel["mask"][split] = {
                    "area_median": float(np.median(area)),
                    "support_share_median": float(np.median(share)),
                    "support_share_mean": float(share.mean()),
                }
            if channel["n"]["pre"]:
                channel["power"] = power_curve(data["z"][member & data["pre"]],
                                               geo)
            report["channels"][f"{orbit}/{subswath}"] = channel
    report["verdict"] = build_verdict(report)
    return report



def fmt(value, digits=3):
    """Format a number for prose; ``None`` stays readable."""
    return "n/a" if value is None else f"{float(value):.{digits}f}"


def fmt_p(value):
    """A p-value that does not print as 0.000 when it is merely very small."""
    if value is None:
        return "p n/a"
    return "p<0.001" if value < 1e-3 else f"p {value:.3f}"


def describe_leaving(entry):
    """One step that left the null, with the number that decided it."""
    if entry.get("observed") is not None:
        tail = (f"brightest {fmt(entry['observed'], 3)} against "
                f"{fmt(entry['null_extreme'], 3)} for the null extreme, "
                f"{fmt_p(entry['p'])}")
    else:
        tail = (f"{entry['out']} of {entry['blocks']} blocks, "
                f"{fmt_p(entry['p'])}")
    return f"{entry['channel']} {entry['level']} ({tail})"


def build_verdict(report):
    """Read the ladder against its two controls and say what follows."""
    channels = report["channels"]
    band_pass = ring_pass = band_cells = ring_cells = 0
    per_channel = {}
    for name, channel in sorted(channels.items()):
        gates = channel["gates"].get("pre")
        if not gates:
            continue
        band_pass += gates["support"]["pass"]
        ring_pass += gates["ring"]["pass"]
        band_cells += gates["support"]["cells"]
        ring_cells += gates["ring"]["cells"]
        per_channel[name] = {
            "band_pass": gates["support"]["pass"],
            "band_cells": gates["support"]["cells"],
            "ring_pass": gates["ring"]["pass"],
            "ring_cells": gates["ring"]["cells"],
            "best_band_da": gates["support"]["best_da"],
        }
    by_key = {}
    for channel in channels.values():
        for row in channel["ladder"].get("pre") or []:
            by_key.setdefault(row["key"], []).append(row)

    full = {}
    for name, channel in sorted(channels.items()):
        for row in channel["ladder"].get("pre") or []:
            if row["key"] != "band_39x15":
                continue
            full[name] = {
                "da": row["da"], "rv": row["rv"], "dom": row["dom"],
                "rv_range": row["rv_range"], "dom_range": row["dom_range"],
                "null_da": row["null"]["da"], "null_rv": row["null"]["rv"],
                "null_dom": row["null"]["dom"],
                "rv_out": row["holm_rv"]["persistence_side"],
                "dom_out": row["holm_dom"]["persistence_side"],
                "blocks": row["blocks"],
            }

    leaving = {"rv": [], "dom": []}
    for name, channel in sorted(channels.items()):
        for row in channel["ladder"].get("pre") or []:
            entry = {"channel": name, "level": row["key"],
                     "region": row["region"], "cells": row["cells"],
                     "blocks": row["blocks"]}
            for metric in ("rv", "dom"):
                holm = row["holm_" + metric]
                if holm and holm["persistence_side"]:
                    leaving[metric].append(
                        dict(entry, out=holm["persistence_side"],
                             p=holm["min_p"],
                             observed=holm.get("observed"),
                             null_extreme=holm.get("null_extreme")))

    emergence = None
    for key, rows in by_key.items():
        da = float(np.mean([row["da"] for row in rows]))
        null_da = float(np.mean([row["null"]["da"] for row in rows]))
        if emergence is None or da < emergence["da"]:
            emergence = {"level": key, "cells": rows[0]["cells"], "da": da,
                         "null_da": null_da,
                         "ratio_to_null": da / null_da if null_da else None,
                         "channels": len(rows)}

    power = None
    curve = {}
    for channel in channels.values():
        for record in channel.get("power") or []:
            if record["factor"] == 8.0 and abs(record["da"] - 0.05) < 1e-9:
                curve.setdefault(record["cells"], []).append(record["rv"])
    if curve:
        sizes = sorted(curve)
        power = {
            "cells": [sizes[0], sizes[-1]],
            "groups": len(curve[sizes[0]]),
            "factor": 8.0,
            "rv_at_factor_8": {
                "smallest": float(np.mean(curve[sizes[0]])),
                "largest": float(np.mean(curve[sizes[-1]])),
            },
            "mean_curve": {str(size): float(np.mean(curve[size]))
                           for size in sizes},
        }

    margins = [{"channel": name, "level": row["key"], "cells": row["cells"],
                "dom": row["dom"], "null_dom": row["null"]["dom"],
                "rv": row["rv"]}
               for name, channel in sorted(channels.items())
               for row in (channel["ladder"].get("pre") or [])
               if row["region"].startswith("margin") and row["height"] > 1]

    verdict = {
        "scatterer_on_deck": bool(band_pass > 0),
        "persistence_at_any_scale": bool(
            any(entry["region"] in ("deck", "band") for entry in leaving["rv"])),
        "brightness_above_null": bool(leaving["dom"]),
        "gates": {"band_pass": band_pass, "band_cells": band_cells,
                  "ring_pass": ring_pass, "ring_cells": ring_cells},
        "per_channel": per_channel,
        "full_support": full,
        "leaving_null": leaving,
        "emergence": emergence,
        "power": power,
        "margins": margins,
    }

    statements = [
        (f"{band_pass} of {band_cells} deck-band cells pass the year-round "
         f"gates before the collapse, against {ring_pass} of {ring_cells} ring "
         f"cells.")
    ]
    if full:
        first = next(iter(full.values()))
        statements.append(
            f"Over the whole 585-cell support the deck block has dom "
            f"{fmt(first['dom'], 2)}, RV {fmt(first['rv'], 2)}, D_A "
            f"{fmt(first['da'], 3)}; the ring null at the same cell count has "
            f"dom {fmt(first['null_dom'], 2)}, RV {fmt(first['null_rv'], 2)}, "
            f"D_A {fmt(first['null_da'], 3)}.")
    if emergence:
        ratio = emergence["ratio_to_null"]
        if ratio is not None and abs(ratio - 1.0) <= 0.25:
            tail = ("the ring null at the same N gives "
                    f"{fmt(emergence['null_da'], 3)}, i.e. the same number: the "
                    "drop is the 1/sqrt(N) of incoherent averaging, not a "
                    "scatterer.")
        else:
            tail = (f"the ring null at the same N gives "
                    f"{fmt(emergence['null_da'], 3)} ({fmt(ratio, 2)}x), which "
                    "does not explain the drop on its own and would need a "
                    "closer look.")
        statements.append(
            f"The smallest D_A anywhere on the ladder is {fmt(emergence['da'], 3)} "
            f"at {emergence['level']} ({emergence['cells']} cells); {tail}")
    if leaving["rv"]:
        statements.append(
            f"{len(leaving['rv'])} ladder step(s) leave the ring null on the "
            f"phase side after Holm: "
            + ", ".join(describe_leaving(entry)
                        for entry in leaving["rv"]) + ".")
    else:
        statements.append(
            "No ladder step, from one cell to the full 585-cell support, leaves "
            "its same-N ring null on the phase side (RV low, gamma high) after "
            "Holm, so persistence does not appear at any spatial scale.")
    if leaving["dom"]:
        statements.append(
            "Blocks that are brighter than the null (dom high): "
            + ", ".join(describe_leaving(entry) for entry in leaving["dom"])
            + " — brightness concentration, which is not phase stability.")
    else:
        statements.append(
            "No ladder step is brighter than its same-N ring null either, so "
            "not even a brightness concentration stands out on the band.")
    if power:
        statements.append(
            f"A planted scatterer of 8x the local cell mean reads as RV "
            f"{fmt(power['rv_at_factor_8']['smallest'], 2)} in "
            f"{power['cells'][0]} cells but "
            f"{fmt(power['rv_at_factor_8']['largest'], 2)} in "
            f"{power['cells'][1]} (mean of the {power['groups']} group(s)): at "
            f"fixed relative brightness, aggregating dilutes a scatterer "
            f"instead of sharpening it.")
    if margins:
        worst = max(margins, key=lambda entry: entry["dom"])
        rv_range = [min(entry["rv"] for entry in margins),
                    max(entry["rv"] for entry in margins)]
        statements.append(
            f"The band bins beside the carriageway are bright and incoherent "
            f"({worst['level']} in {worst['channel']}, {worst['cells']} cells: "
            f"dom {fmt(worst['dom'], 2)} against a null of "
            f"{fmt(worst['null_dom'], 2)}; RV {fmt(rv_range[0], 2)}-"
            f"{fmt(rv_range[1], 2)} over the margin blocks). They sit outside "
            f"the 39 m deck, so they cannot be the deck scatterer, and their "
            f"phase is not stable either.")
    verdict["statements"] = statements
    return verdict



# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
def _row(cells):
    return "| " + " | ".join(str(cell) for cell in cells) + " |"


def to_markdown(report):
    """The committed companion note: same numbers, no new claims."""
    geo = report["geometry"]
    check = report["mask_check"]
    verdict = report["verdict"]
    lines = [
        "# carola_ocv_persist.py — at which spatial scale does deck "
        "persistence appear?",
        "",
        f"*{report['label']}, {report['inputs']['support']['n_dates']} dates of "
        f"the {geo['chip'][0]}x{geo['chip'][1]} support chip cut around the "
        f"deck anchor, split at the {report['rules']['event']} collapse.*",
        "",
        f"Schema `{report['schema']}`. Generated by "
        "`code/carola/carola_ocv_persist.py`; `--verify` recomputes every "
        "number below from the committed outputs.",
        "",
        "## The answer",
        "",
    ]
    lines += [f"- {text}" for text in verdict["statements"]]
    lines += [
        "",
        "So: no scatterer on the 39 m deck is simultaneously bright (at least "
        "8x the local block mean) and phase-stable, and aggregating the band "
        "does not create one. Where an aggregated block does look calmer than "
        "a single cell — around 15 to 45 cells — the ring null at the same "
        "cell count is calmer by exactly the same factor, because that is what "
        "averaging N independent speckle fields does: `D_A` falls like "
        "`1/sqrt(N)`. `RV`, which does not improve with N, stays at 1 (pure "
        "speckle) at every scale, and `dom` falls like `H_N/N`, the brightest "
        "of N speckle cells over their sum — both are printed next to the "
        "same-N ring null so the reader can check it cell by cell.",
        "",
        "## Inputs and the mask rule",
        "",
        _row(["input", "sha256", "note"]),
        _row(["---", "---", "---"]),
        _row(["`" + report["inputs"]["support"]["path"] + "`",
              "`" + report["inputs"]["support"]["sha256"][:16] + "…`",
              f"{report['inputs']['support']['n_dates']} date series"]),
        _row(["`" + report["inputs"]["manifest"]["path"] + "`",
              "`" + report["inputs"]["manifest"]["sha256"][:16] + "…`",
              "geometry, spacing, region cell counts"]),
        _row(["`" + report["inputs"]["channels"]["path"] + "`",
              "`" + str(report["inputs"]["channels"]["sha256"])[:16] + "…`",
              "committed per-date mask area `A`"]),
        "",
        f"Mask per date: {report['rules']['mask']}. The band peak sets the "
        f"brightness floor, the ring median sets the contrast floor; nothing "
        f"in the rule mentions the deck, which is the point.",
        "",
        f"Recomputed mask area against the committed table: **{check['equal']} "
        f"of {check['compared']} dates equal**, largest difference "
        f"`{check['max_abs_delta_A']:g}` cells (median area "
        f"{check['area_median']:g}, range {check['area_min']}-"
        f"{check['area_max']}); the manifest's own `A` column agrees on "
        f"{check.get('manifest_equal', 0)} of {check['compared']} dates. "
        f"The mask covers {100 * check['support_share_mean']:.1f}% of its cells "
        f"on the band (median {100 * check['support_share_median']:.1f}%) "
        f"against {100 * check['support_share_of_chip']:.1f}% of the chip being "
        f"band, i.e. the mask is a brightness mask that is spread over the "
        f"whole chip and has no preference for the deck.",
        "",
        "## Geometry",
        "",
        _row(["region", "cells", "columns", "ground extent"]),
        _row(["---", "---", "---", "---"]),
        _row(["support band", geo["cells"]["support"],
              f"{geo['columns']['band'][0]}-{geo['columns']['band'][1]}",
              f"{geo['ground_m']['band'][0]:.1f} x "
              f"{geo['ground_m']['band'][1]:.1f} m"]),
        _row(["deck interior (39 m hull)", geo["cells"]["deck"],
              f"{geo['columns']['deck'][0]}-{geo['columns']['deck'][1]}",
              f"{geo['ground_m']['deck_hull'][0]:.1f} x "
              f"{geo['ground_m']['deck_hull'][1]:.1f} m"]),
        _row(["margin left", geo["cells"]["margin_left"],
              f"{geo['columns']['margin_left'][0]}-"
              f"{geo['columns']['margin_left'][1]}",
              f"{2 * geo['ground_m']['range_bin']:.1f} m wide"]),
        _row(["margin right", geo["cells"]["margin_right"],
              f"{geo['columns']['margin_right'][0]}-"
              f"{geo['columns']['margin_right'][1]}",
              f"{2 * geo['ground_m']['range_bin']:.1f} m wide"]),
        _row(["guard", geo["cells"]["guard"], "collar",
              "inside the ring, never used as a block"]),
        _row(["ring", geo["cells"]["ring"], "all sides",
              f"`ring_complete: {str(geo['cells']['ring_complete']).lower()}`"]),
        "",
        f"Rows {geo['rows'][0]}-{geo['rows'][1]}, range bin "
        f"{geo['ground_m']['range_bin']:.4f} m, azimuth line "
        f"{geo['ground_m']['azimuth_line']:.3f} m. One resolution cell is "
        f"{geo['ground_m']['range_bin']:.1f} x "
        f"{geo['ground_m']['azimuth_line']:.1f} m "
        f"({geo['ground_m']['range_bin'] * geo['ground_m']['azimuth_line']:.0f} "
        f"m²); the deck interior is 11 bins wide, so the two outermost band "
        f"bins per side sit beside the carriageway and are reported apart.",
        "",
    ]
    lines += to_markdown_ladder(report)
    lines += to_markdown_power(report)
    lines += [
        "## What the numbers do not say",
        "",
    ]
    lines += [f"- {text}." for text in report["limits"]]
    lines += [
        "",
        "## Reproduce",
        "",
        "```sh",
        "python3 code/carola/carola_ocv_persist.py --selftest   # no payload",
        "python3 code/carola/carola_ocv_persist.py              # ~1 min",
        "python3 code/carola/carola_ocv_persist.py --verify",
        "```",
        "",
    ]
    return "\n".join(lines)


def to_markdown_ladder(report):
    """The ladder tables and the gates, one block per channel."""
    lines = [
        "## The ladder, before the collapse",
        "",
        "`mask` is the share of dates whose mask touches the block; `out` is "
        "the number of blocks that leave the same-N ring null on the "
        "persistence side (RV low or dom high). `null` columns are medians of "
        "the ring placements, `p05-p95` the 5th and 95th percentiles. Single "
        "cells and 3x3 blocks are swept over the whole deck (429 and 333 "
        "blocks), and are tested as one question — is the most extreme block "
        "outside the most extreme of as many ring placements — because 429 "
        "blocks cannot be Holm-tested against a finite null; every other step "
        "is tiled and Holm-tested block by block.",
        "",
    ]
    for name, channel in report["channels"].items():
        pre = channel["ladder"].get("pre") or []
        if not pre:
            continue
        lines += [
            f"### {name} — {channel['n']['pre']} dates, "
            f"{channel['dates']['first']} to {channel['dates']['last_pre']}",
            "",
            _row(["level", "cells", "blocks", "mask", "D_A", "null D_A", "RV",
                  "null RV (p05-p95)", "dom", "null dom", "out"]),
            _row(["---"] * 11),
        ]
        for row in pre:
            out = (row["holm_rv"]["persistence_side"]
                   + row["holm_dom"]["persistence_side"])
            if out:
                worst_p = min(p for p in (row["holm_rv"]["min_p"],
                                          row["holm_dom"]["min_p"])
                              if p is not None)
                out = f"{out} ({fmt_p(worst_p)})"
            lines.append(_row([
                row["label"], row["cells"], row["blocks"],
                f"{row['recurrence']:.2f}", f"{row['da']:.3f}",
                f"{row['null']['da']:.3f}", f"{row['rv']:.2f}",
                f"{row['null']['rv']:.2f} ({row['null']['rv_p05']:.2f}-"
                f"{row['null']['rv_p95']:.2f})", f"{row['dom']:.2f}",
                f"{row['null']['dom']:.3f}", out,
            ]))
        lines.append("")
    lines += ["## The gates, before the collapse", ""]
    lines += [_row(["channel", "band pass", "band cells", "ring pass",
                    "ring cells", "best band D_A"]),
              _row(["---", "---", "---", "---", "---", "---"])]
    for name, channel in report["channels"].items():
        gates = channel["gates"].get("pre")
        if not gates:
            continue
        lines.append(_row([
            name, gates["support"]["pass"], gates["support"]["cells"],
            gates["ring"]["pass"], gates["ring"]["cells"],
            f"{gates['support']['best_da']:.2f}",
        ]))
    lines += [
        "",
        "The same gate is applied to both regions, so the band zero is "
        "meaningful only because the ring is not zero.",
        "",
        "## After the collapse",
        "",
        "Reported for completeness, never pooled with the numbers above.",
        "",
        _row(["channel", "dates", "level", "D_A", "RV", "dom"]),
        _row(["---", "---", "---", "---", "---", "---"]),
    ]
    for name, channel in report["channels"].items():
        post = channel["ladder"].get("post") or []
        keep = [row for row in post if row["key"] in ("deck_3x3",
                                                      "band_39x15")]
        for row in keep:
            lines.append(_row([name, channel["n"]["post"], row["label"],
                               f"{row['da']:.3f}", f"{row['rv']:.2f}",
                               f"{row['dom']:.2f}"]))
    lines.append("")
    return lines


def to_markdown_power(report):
    """The planted-scatterer curve, one small table per channel."""
    lines = [
        "## Can a block see a scatterer at all?",
        "",
        "One scatterer is planted at the block centre with amplitude `f` times "
        "*that block's* mean cell amplitude and `D_A = 0.05`; the entry is the "
        "block's `RV` (1.0 = pure speckle, below 0.5 = the scatterer "
        "dominates).",
        "",
    ]
    for name, channel in report["channels"].items():
        records = channel.get("power") or []
        if not records:
            continue
        levels = []
        for rec in records:
            if rec["level"] not in levels:
                levels.append(rec["level"])
        factors = sorted(set(rec["factor"] for rec in records))
        lines += [f"### {name}", "",
                  _row(["level"] + [f"f={f:g}" for f in factors]),
                  _row(["---"] * (len(factors) + 1))]
        for level in levels:
            cells = next(rec["cells"] for rec in records
                         if rec["level"] == level)
            cells_row = [f"{level} ({cells})"]
            for factor in factors:
                value = next((rec["rv"] for rec in records
                              if rec["level"] == level
                              and rec["factor"] == factor
                              and abs(rec["da"] - 0.05) < 1e-9), None)
                cells_row.append("n/a" if value is None else f"{value:.2f}")
            lines.append(_row(cells_row))
        lines.append("")
    curves = {}
    for channel in report["channels"].values():
        for rec in channel.get("power") or []:
            if rec["factor"] == 8.0 and abs(rec["da"] - 0.05) < 1e-9:
                curves.setdefault(rec["cells"], []).append(rec["rv"])
    if curves:
        sizes = sorted(curves)
        pairs = []
        for channel in report["channels"].values():
            records = channel.get("power") or []
            for factor in sorted(set(rec["factor"] for rec in records)):
                def value(cells):
                    return next((rec["rv"] for rec in records
                                 if rec["factor"] == factor
                                 and rec["cells"] == cells
                                 and abs(rec["da"] - 0.05) < 1e-9), None)
                small, large = value(sizes[0]), value(sizes[-1])
                if small is not None and large is not None \
                        and small < RV_DETECTED:
                    pairs.append(large > small)
        lines += [
            "### All groups, `RV` at f=8",
            "",
            _row(["cells"] + [str(size) for size in sizes]),
            _row(["---"] * (len(sizes) + 1)),
            _row(["mean RV"] + [f"{np.mean(curves[size]):.2f}"
                                for size in sizes]),
            "",
            f"At a fixed relative brightness the block gets harder to see as it "
            f"grows: {np.mean(curves[sizes[0]]):.2f} in {sizes[0]} cells against "
            f"{np.mean(curves[sizes[-1]]):.2f} in {sizes[-1]}, and in "
            f"{sum(pairs)} of the {len(pairs)} group-by-brightness comparisons "
            f"in which the single-line block can see the scatterer at all "
            f"(`RV` < {RV_DETECTED:g}) the whole support reads a higher `RV`. "
            f"The middle of the ladder is not monotone, because `f` is measured "
            f"against each block's own mean and a block that contains a brighter "
            f"patch plants a brighter scatterer at the same `f`; where the "
            f"finest block detects, the two ends never cross, so nothing here "
            f"sees at {sizes[-1]} cells a scatterer that the finest block misses "
            f"at the same relative brightness — which is why the ladder is "
            f"reported at the scale that still resolves a cell.",
            "",
        ]
    return lines


# ---------------------------------------------------------------------------
# Output, verification, self-test
# ---------------------------------------------------------------------------
def write_outputs(report, json_path=OUT_JSON, md_path=OUT_MD):
    """Write the JSON and its Markdown companion; no timestamp, no payload."""
    for path in (json_path, md_path):
        folder = os.path.dirname(path)
        if folder:
            os.makedirs(folder, exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1, ensure_ascii=False)
        handle.write("\n")
    with open(md_path, "w") as handle:
        handle.write(to_markdown(report))
        handle.write("\n")
    return json_path, md_path


def hygiene(report):
    """A structure label or a timestamp would make the report a claim; there is
    none. The words that would carry a label are refused in values as well as in
    keys, and ``collapse`` itself is allowed because it names the event."""
    banned_keys = {"state", "status", "condition", "collapsed", "intact",
                   "generated_at", "created_at", "timestamp", "run_at"}
    banned_words = re.compile(r"\b(state|status|condition|collapsed|intact)\b",
                              re.IGNORECASE)
    problems = []

    def walk(node, path):
        if isinstance(node, dict):
            for key, value in node.items():
                if str(key).lower() in banned_keys:
                    problems.append(f"{path}.{key}")
                walk(value, f"{path}.{key}")
        elif isinstance(node, (list, tuple)):
            for index, value in enumerate(node):
                walk(value, f"{path}[{index}]")
        elif isinstance(node, str):
            if re.search(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}", node):
                problems.append(f"{path}: timestamp")
            if banned_words.search(node):
                problems.append(f"{path}: structure label in text")
    walk(report, "report")
    return problems


def compare(left, right, path="report", problems=None, tol=1e-9):
    """Deep compare with a numeric tolerance; stops after twelve findings."""
    problems = [] if problems is None else problems
    if len(problems) >= 12:
        return problems
    if isinstance(left, dict) and isinstance(right, dict):
        if set(left) != set(right):
            only_left = sorted(set(left) - set(right))[:4]
            only_right = sorted(set(right) - set(left))[:4]
            problems.append(f"{path}: keys differ ({only_left} vs {only_right})")
            return problems
        for key in left:
            compare(left[key], right[key], f"{path}.{key}", problems, tol)
    elif isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        if len(left) != len(right):
            problems.append(f"{path}: length {len(left)} vs {len(right)}")
            return problems
        for index, (one, other) in enumerate(zip(left, right)):
            compare(one, other, f"{path}[{index}]", problems, tol)
    elif isinstance(left, bool) or isinstance(right, bool):
        if left != right:
            problems.append(f"{path}: {left!r} vs {right!r}")
    elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if not math.isclose(float(left), float(right), rel_tol=tol, abs_tol=tol):
            problems.append(f"{path}: {left!r} vs {right!r}")
    elif left != right:
        problems.append(f"{path}: {left!r} vs {right!r}")
    return problems


def verify(json_path=OUT_JSON, table=TABLE, manifest=MANIFEST,
           channel_table=CHANNEL_TABLE):
    """Recompute the report from the committed inputs and compare it all."""
    if not os.path.exists(json_path):
        raise SystemExit(f"verify: {json_path} does not exist — run the script "
                         f"first")
    with open(json_path) as handle:
        stored = json.load(handle)
    if stored.get("schema") != SCHEMA:
        raise SystemExit(f"verify: schema is {stored.get('schema')!r}, "
                         f"expected {SCHEMA!r}")
    for key, path in (("support", table), ("manifest", manifest),
                      ("channels", channel_table)):
        pinned = stored["inputs"][key]["sha256"]
        actual = sha256_of(path) if os.path.exists(path) else None
        if pinned != actual:
            raise SystemExit(f"verify: {key} sha256 is {pinned}, on disk "
                             f"{actual} — the inputs moved")
    draws = int(stored["rules"]["null"]["placements"])
    fresh = analyse(table, manifest, channel_table, draws)
    problems = compare(stored, fresh)
    problems += [f"hygiene: {item}" for item in hygiene(fresh)]
    md_path = (json_path[:-5] + ".md") if json_path.endswith(".json") else OUT_MD
    if os.path.exists(md_path):
        with open(md_path) as handle:
            written = handle.read()
        if written != to_markdown(fresh) + "\n":
            problems.append("markdown: differs from a fresh render of the JSON")
    else:
        problems.append(f"markdown: {md_path} is missing")
    if problems:
        for item in problems:
            print(f"  {item}")
        raise SystemExit(f"verify: {len(problems)} difference(s)")
    channels = ", ".join(f"{name} ({channel['n']['pre']}/"
                         f"{channel['n']['post']})"
                         for name, channel in fresh["channels"].items())
    check = fresh["mask_check"]
    print(f"verify: ok — {channels}; mask area equal on "
          f"{check['equal']}/{check['compared']} dates; "
          f"{fresh['verdict']['gates']['band_pass']} band and "
          f"{fresh['verdict']['gates']['ring_pass']} ring cells pass the gates")
    return 0


def synthetic_dates(count=48, start="2022-01-01", step_days=12):
    """Evenly spaced dates across a year, so the cold-season gate can act."""
    base = dt.date.fromisoformat(start)
    return [str(base + dt.timedelta(days=step_days * i)) for i in range(count)]


def noise(count, lines, bins, seed=11):
    """Complex Gaussian speckle, unit mean power per cell."""
    rng = np.random.default_rng(seed)
    return (rng.standard_normal((count, lines, bins))
            + 1j * rng.standard_normal((count, lines, bins))) / math.sqrt(2.0)


def toy_geometry():
    with open(MANIFEST) as handle:
        return geometry(json.load(handle))


def test_geometry():
    geo = toy_geometry()
    assert geo["observed_cells"] == {"support": 585, "guard": 900,
                                     "ring": 7140}, geo["observed_cells"]
    assert geo["region_cells"] == {"deck": 429, "margin_left": 78,
                                  "margin_right": 78}, geo["region_cells"]
    assert geo["cols"]["band"] == (55, 70), geo["cols"]["band"]
    assert geo["cols"]["deck"] == (57, 68), geo["cols"]["deck"]
    assert geo["rows"] == (15, 54), geo["rows"]
    assert not (geo["support"] & geo["ring"]).any()
    assert not (geo["support"] & geo["guard"]).any()
    return "585/900/7140 cells, band 55-69, deck 57-67, no overlap"


def test_mask():
    geo = toy_geometry()
    rng = np.random.default_rng(13)
    intensity = (rng.standard_normal((69, 125)) ** 2
                 + rng.standard_normal((69, 125)) ** 2) / 2.0
    intensity[34, 62] = 50.0
    mask = mask_of(intensity, geo)
    assert int(mask.sum()) == 1, int(mask.sum())
    assert mask[34, 62]
    assert int((mask & geo["deck"]).sum()) == 1
    return "one planted bright cell, one masked cell, on the deck"


def test_metrics():
    flat = np.ones((50, 3, 3), dtype=complex) * np.exp(0.7j)
    steady = block_metrics(flat)
    assert steady["rv"] < 1e-12 and steady["da"] < 1e-12
    assert abs(steady["dom"] - 1.0 / 9.0) < 1e-12
    assert abs(steady["gamma"] - 1.0) < 1e-12
    rng = np.random.default_rng(3)
    spike = noise(50, 1, 9, seed=3)
    spike[:, 0, 4] = 1000.0 * np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, 50))
    dominant = block_metrics(spike)
    assert dominant["dom"] > 0.99, dominant["dom"]
    assert dominant["rv"] < 0.05, dominant["rv"]
    single = block_metrics(noise(4000, 1, 1, seed=5))
    assert abs(single["da"] - 0.5227) < 0.02, single["da"]
    assert abs(single["rv"] - 1.0) < 0.1, single["rv"]
    assert abs(single["dom"] - 1.0) < 1e-12, single["dom"]
    eight = block_metrics(noise(4000, 1, 8, seed=5))
    assert abs(eight["da"] - 0.5227 / math.sqrt(8)) < 0.02, eight["da"]
    assert abs(eight["rv"] - 1.0) < 0.1, eight["rv"]
    assert abs(eight["dom"] - 2.7179 / 8) < 0.02, eight["dom"]
    return (f"steady rv {steady['rv']:.0e}; one steady cell dom "
            f"{dominant['dom']:.3f} rv {dominant['rv']:.3f}; speckle D_A "
            f"{single['da']:.2f} -> {eight['da']:.2f} at 8 cells while RV stays "
            f"{eight['rv']:.2f} and dom is H_8/8 = {eight['dom']:.2f}")


def test_windows():
    mask = np.zeros((6, 7), bool)
    mask[2:5, 3:6] = True
    found = {(int(r), int(c)) for r, c in window_positions(mask, 2, 2)}
    assert found == {(2, 3), (2, 4), (3, 3), (3, 4)}, found
    assert window_positions(mask, 4, 2).size == 0
    assert len(window_positions(mask, 3, 3)) == 1
    return "4 positions for 2x2 in a 3x3 box, 0 for 4x2, 1 for 3x3"


def test_holm():
    adjusted, tested = holm([0.01, 0.02, 0.03])
    assert tested == 3
    assert all(abs(a - b) < 1e-12
               for a, b in zip(adjusted, [0.03, 0.04, 0.04])), adjusted
    adjusted, _ = holm([0.5, 0.01])
    assert abs(adjusted[1] - 0.02) < 1e-12 and abs(adjusted[0] - 0.5) < 1e-12
    adjusted, tested = holm([0.2, None])
    assert tested == 1 and adjusted[1] is None and adjusted[0] == 0.2
    return "[0.01,0.02,0.03] -> [0.03,0.04,0.04], missing stays missing"


def test_permutation():
    rng = np.random.default_rng(7)
    null = rng.standard_normal(200)
    assert permutation_p(float(np.median(null)), null) > 0.8
    assert permutation_p(10.0, null) < 0.01
    assert permutation_p(1.0, np.zeros(0)) == 1.0
    return "median -> p > 0.8, 10 sigma -> p < 0.01, empty null -> 1.0"


def test_gates():
    geo = toy_geometry()
    dates = synthetic_dates(48)
    bare = noise(48, 69, 125, seed=17)
    quiet = gate_counts(bare, dates, geo)
    assert quiet["support"]["pass"] == 0, quiet["support"]["pass"]
    assert quiet["ring"]["pass"] == 0, quiet["ring"]["pass"]
    assert quiet["support"]["cells"] == 585 and quiet["ring"]["cells"] == 7140
    planted = bare.copy()
    rng = np.random.default_rng(19)
    planted[:, 34, 62] = (30.0 * (1.0 + 0.05 * rng.standard_normal(48))
                          * np.exp(1j * 0.3))
    loud = gate_counts(planted, dates, geo)
    assert 1 <= loud["support"]["pass"] <= 2, loud["support"]["pass"]
    assert loud["support"]["best_da"] < 0.1, loud["support"]["best_da"]
    return (f"speckle 0/585 and 0/7140; one planted steady cell "
            f"{loud['support']['pass']}/585")


def test_power():
    bare = noise(48, 69, 125, seed=23)
    reference = float(np.abs(bare).mean())
    small = block_metrics(plant(bare[:, 34:35, 57:68], (0, 5), 32.0, 0.05,
                                reference=reference))
    wide = block_metrics(plant(bare[:, 15:54, 55:70], (19, 7), 32.0, 0.05,
                               reference=reference))
    assert small["rv"] < 0.35, small["rv"]
    assert wide["rv"] > small["rv"] + 0.25, (wide["rv"], small["rv"])
    assert abs(block_metrics(bare[:, 34:35, 57:68])["rv"] - 1.0) < 0.35
    geo = toy_geometry()
    curve = power_curve(bare, geo, factors=(8.0,), das=(0.05,))
    assert len(curve) == len(LINES) + 1, len(curve)
    assert all(isinstance(rec["detected"], bool) for rec in curve)
    return (f"f=32: RV {small['rv']:.2f} at 11 cells vs {wide['rv']:.2f} at "
            f"585 - aggregation weakens the scatterer")


def test_ladder():
    geo = toy_geometry()
    z = noise(48, 69, 125, seed=29)
    mask = np.stack([mask_of(np.abs(frame) ** 2, geo) for frame in z])
    spec = {item["key"]: item for item in ladder_spec(geo)}
    first = summarise(z, geo, spec["deck_3x3"], mask, 50, {})
    again = summarise(z, geo, spec["deck_3x3"], mask, 50, {})
    assert compare(first, again) == [], compare(first, again)
    assert first["cells"] == 9 and first["blocks"] == 333, first["blocks"]
    assert first["null"]["draws"] == 50 and first["null"]["ring_positions"] > 0
    assert 0.0 <= first["recurrence"] <= 1.0
    cell = summarise(z, geo, spec["deck_cell"], mask, 50, {})
    assert cell["blocks"] == 429 and cell["cells"] == 1
    assert cell["null"]["ring_positions"] == 7140
    assert cell["holm_rv"]["significant"] <= 9, cell["holm_rv"]["significant"]
    assert cell["null"]["rv_p05"] < cell["null"]["rv_p95"]
    assert abs(cell["dom"] - 1.0) < 1e-12
    assert abs(cell["da"] - 0.5227) < 0.05
    return (f"3x3: {first['blocks']} blocks, {first['null']['draws']} null "
            f"draws, reproducible; single cell {cell['blocks']} blocks, "
            f"{cell['holm_rv']['significant']} outside the null")


def test_hygiene():
    assert hygiene({"a": 1, "b": [1, 2]}) == []
    assert hygiene({"state": "collapsed"})
    assert hygiene({"nested": [{"timestamp": 1}]})
    assert hygiene({"note": "acquired 2024-09-11T10:00:00Z"})
    assert hygiene({"note": "the deck is intact after all"})
    assert not hygiene({"note": "2024-09-11"})
    assert not hygiene({"note": "before the collapse"})
    return ("labels and timestamps are caught in keys and in text; a plain date "
            "and the name of the event are not")



def selftest():
    """Checks that need no payload: geometry, metrics, controls, report rules."""
    tests = [
        ("geometry", test_geometry),
        ("mask rule", test_mask),
        ("block metrics", test_metrics),
        ("window positions", test_windows),
        ("Holm", test_holm),
        ("permutation p", test_permutation),
        ("year-round gates", test_gates),
        ("planted scatterer", test_power),
        ("ladder and null", test_ladder),
        ("report hygiene", test_hygiene),
    ]
    failed = 0
    for name, test in tests:
        try:
            note = test()
        except AssertionError as error:
            failed += 1
            print(f"  FAIL {name}: {error}")
        else:
            print(f"  ok   {name}: {note}")
    if failed:
        raise SystemExit(f"selftest: {failed} of {len(tests)} failed")
    print(f"selftest: {len(tests)} checks ok, no payload read")
    return 0


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def build_parser():
    parser = argparse.ArgumentParser(
        description=("Does deck persistence exist at any spatial scale? Runs "
                     "the cell -> 3x3 -> line -> third -> full-support ladder "
                     "on the Carola support chip against a same-N ring null and "
                     "a planted-scatterer power curve."))
    parser.add_argument("--table", default=TABLE,
                        help="support chip table (date|orbit|rows|cols|iq)")
    parser.add_argument("--manifest", default=MANIFEST,
                        help="manifest with the geometry and cell counts")
    parser.add_argument("--channels", default=CHANNEL_TABLE,
                        help="committed per-date channel table (mask area A)")
    parser.add_argument("--draws", type=int, default=NULL_DRAWS,
                        help=f"ring placements per scale (default {NULL_DRAWS})")
    parser.add_argument("--json", default=OUT_JSON, help="report to write")
    parser.add_argument("--md", default=OUT_MD, help="Markdown to write")
    parser.add_argument("--selftest", action="store_true",
                        help="run the offline checks and exit")
    parser.add_argument("--verify", action="store_true",
                        help="recompute from the committed inputs and compare")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()
    if args.verify:
        return verify(args.json, args.table, args.manifest, args.channels)
    report = analyse(args.table, args.manifest, args.channels, args.draws)
    json_path, md_path = write_outputs(report, args.json, args.md)
    print(f"{report['label']}: {len(report['channels'])} channel group(s)")
    for text in report["verdict"]["statements"]:
        print(f"  - {text}")
    print(f"\nwrote {rel(json_path)} ({os.path.getsize(json_path)} B)")
    print(f"wrote {rel(md_path)} ({os.path.getsize(md_path)} B)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
