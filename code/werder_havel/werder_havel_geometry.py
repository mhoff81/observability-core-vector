#!/usr/bin/env python3
"""WERDER HAVEL - the derived metre per pixel of the record, and what its 10 % decile buys.

File: ``data/werder_havel/reference/werder_havel_geometry.json`` (derived here, offline,
from the two committed reference files).

The question this file answers
------------------------------
The anchor hands over the structure in metres and *refuses* to turn it into pixels:
``werder_havel_aoi.json`` commits the outline as 704.043 x 35.929 m (25 295.7 m2 at
bearing 59.59 deg) and the record as an **80 x 80 window per acquisition** whose
``coherence_masked_pixels`` is 640 in every one of the 991 rows. That 640 is 10 % of
6400 **by definition** — the brightest decile of the window's mean amplitude — and
``werder_havel_ocv_masks.py`` says so and stops. A decile is only a *share of a window*;
to know what it covers on the ground, one needs the metres per pixel, and the window
carries no pixel-spacing spec. Only one thing in this package can supply it: the
**committed burst footprints** of ``werder_havel_bursts_resolve.py`` — 609 frames in
1218 catalogue records over the 491 dates, each with the catalogue's own ground polygon
and its ``LinesPerBurst`` / ``SamplesPerBurst``. This file is that division.

The grid, and why it is a division and not a lookup
---------------------------------------------------
Sentinel-1 SLC pixel spacings are a *product* property (range ~2.33 m slant, azimuth
~13.9 m), not a ground distance: the ground distance depends on the local incidence
angle, which changes along the swath and with the orbit tube. The family's method
(``haltern_osm_anchor.py``, YWF before it) therefore measures the grid **off the
catalogue's own footprint polygon** rather than quoting a constant: the longest edge of
the footprint is the **range** axis, its collinear pieces are summed, the vertex spreads
along and across that axis give the two ground extents, and

    range_m_per_px   = footprint range extent / SamplesPerBurst
    azimuth_m_per_px = footprint azimuth extent / LinesPerBurst

The result for this site, per pass geometry (medians over that geometry's frames):

    asc_146_iw1   3.77 m/px range x 15.45 m/px azimuth   58.3 m2/px
    asc_44_iw3    3.12 m/px range x 15.30 m/px azimuth   47.7 m2/px
    desc_168_iw1  3.77 m/px range x 15.39 m/px azimuth   58.0 m2/px
    desc_95_iw3   3.11 m/px range x 15.21 m/px azimuth   47.3 m2/px

— i.e. the IW3 geometries see this deck at a higher incidence than the IW1 ones, so a
range pixel is 3.77 m of ground there against 3.12 m here: a 21 % difference that decides
the fill ratio below. The spread is committed beside every median, because it is real and
not noise: the *same* BurstId returns a slightly different footprint on different dates
(the orbit tube), 3.72..3.96 m/px on ascending IW1 and 3.08..3.13 m/px on ascending IW3.

What the grid prices
--------------------
The deck is anisotropic and so is the window, so the extent is projected, not divided:
at bearing 59.59 deg the deck lies 17-42 deg off the range axis of the four geometries,
which is why its range extent is 548..678 m but its azimuth extent only 262..497 m.
Against that:

  * **the deck in pixels** — 434/436 px on the IW1 geometries, 530/535 px on the IW3
    ones. The exported decile is 640 px, so the deck accounts for **0.68 of it on IW1 and
    0.83/0.84 on IW3** — which is what the external calibration fixture claims at 0.86
    ("geometric object fill ratio 0.86 for a 704 m x ~40 m deck"): the claim is an **IW3
    reading**, reproduced to 2.8-3.7 %, and missed by 21 % on IW1. The file commits the
    claim, its quote, its source, both derived readings and the residual;
  * **containment** — it does *not* hold. The window covers 249..302 m of ground along
    range and 1217..1236 m along azimuth; the deck needs 548..678 m along range, so an
    80 x 80 window contains at most about half of the deck's length and 147..219 px per
    side would be needed to contain it. The record's window is therefore a sample of the
    deck's *middle*, not of the deck — a fact the fill ratio cannot show, which is why it
    is committed here as geometry and not inferred later from a mask;
  * **the deck's share of the window** — 6.8 % on IW1, 8.3/8.4 % on IW3. Together with
    the claim this closes the fixture's own arithmetic: a deck of 8.4 % of the window that
    fills 0.86 of the mask makes the mask 9.6-9.7 % of the window, i.e. the "brightest
    decile" the export names. The claim's implied pixel area is 45.96 m2 (3.00 m x 14.73 m
    on one IW3 geometry, 3.02 m x 14.79 m on the other), 2.8-3.7 % finer than the IW3
    footprints give; that residual is committed below rather than smoothed away.

What this file does *not* claim
-------------------------------
``status``: *derived metres per pixel — no pixel layer is claimed*. Three limits are
committed with the numbers rather than left to the reader:

  * the reference and outside windows (``ref_*``, ``refout_*``, 250 + 250 rows) lie
    **outside** the queried ring, and the export carries no coordinate for any row, so
    their own metres per pixel is *not* derivable from this layer. That the grid above
    applies to them is true only as far as "same relative orbit and sub-swath, and within
    the ±4 % the orbit tube moves along the swath";
  * the deck is the **map's** outline, not a survey: the family has measured a -81..+7 m
    shape error (rms 56 m) for a polygon anchor, and 7 m of shape error is 1.9 px at the
    IW1 range spacing - 0.4 % of the deck's 434 px on these grids;
  * the footprint polygons are the catalogue's, rounded to ~1e-5 deg (1.1 m in latitude),
    on an 80 km range extent — 1e-5 of the extent, i.e. the grids are limited by the
    orbit's real variation (the spreads above), not by that rounding.

Modes
-----
    --write   derive the geometry from the committed anchor and burst layer and write it
    --check   re-derive and compare field by field (the default; offline)
    --json    print the derived file as JSON (the numbers the notes above quote)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import werder_havel_aoi as aoi  # noqa: E402
import werder_havel_bursts_resolve as bursts_mod  # noqa: E402
import werder_havel_ocv_masks as masks  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "werder_havel")
AOI_PATH = aoi.OUT
BURSTS_PATH = bursts_mod.OUT
OUT = os.path.join(DATA, "reference", "werder_havel_geometry.json")

SITE = "Havelbruecke Werder (A 10, Brandenburg) - the bridge outline as the anchor"
SRC_AOI = "data/werder_havel/reference/werder_havel_aoi.json"
SRC_BURSTS = "data/werder_havel/reference/werder_havel_bursts.json"
SRC_RECORD = "data/werder_havel/werder_havel_measurements_full.txt"


# The family's equirectangular constants, as the anchor's own ``metrics`` block uses
# them (110540.0 m/deg latitude, 111320.0 m/deg longitude at cos(lat0) = 0.609903).
M_PER_DEG_LAT = 110540.0
M_PER_DEG_LON = 111320.0

WINDOW_SIDE_PX = masks.WINDOW_SIZE                    # 80 - the backend's window here
WINDOW_PX = masks.WINDOW_PX                           # 6400
MASKED_PX = masks.RECORD_EXPECTED["masked_pixels"]    # 640 = 10 % of the window
WINDOW_NOTE = masks.WINDOW_NOTE

# The SLC's own range sample spacing, **slant** - not used to derive anything here, only
# to price what using it would do (see ``area_corrections``).
SLANT_RANGE_M = 2.3296

# The one external claim this file tests. It is *not* in this repository: it is the
# external chain's ASCENDING coherence calibration fixture, quoted here verbatim with its
# path so a reader sees both the claim and the grid it is priced against.
FILL_CLAIM = 0.86
FILL_CLAIM_QUOTE = ("window 80x80 px, mask = brightest decile of the pre-event mean "
                    "amplitude (geometric object fill ratio 0.86 for a 704 m x ~40 m deck)")
FILL_CLAIM_SOURCE = ("/home/projects/havel_werder_analysis/corrected2/fixtures/"
                     "werder_h1_coherence_asc.json (the external chain's ASCENDING "
                     "coherence fixture, 2026-09-21) - outside this repository, mutable "
                     "external input: the quote is committed here, the file is not")

# The load-bearing facts of the anchor this file divides by. They cannot be re-derived
# from the record (the export carries no coordinate), so they are asserted - the family's
# ``*_EXPECTED`` pattern - while the derivation of the numbers stays in
# ``werder_havel_aoi.py``.
ANCHOR_EXPECTED = {
    "osm_way": 413381418,
    "version": 5,
    "changeset": 178106321,
    "n_distinct_nodes": 9,
    "kind": "structure_outline",
    "span_m": 704.043,
    "width_m": 35.929,
    "area_m2": 25295.7,
    "span_bearing_deg": 59.59,
    "wgs84_area_m2": 25075.8,
    "site_point": {"lat": 52.4175279, "lon": 12.9146744},
}

# The facts this file needs from the covering layer, asserted for the same reason.
BURSTS_EXPECTED = {
    "n_dates": 491,
    "n_frames": 609,
    "n_burst_records": 1218,
    "tracks": ["asc_146_iw1", "asc_44_iw3", "desc_168_iw1", "desc_95_iw3"],
    "frames_per_track": {"asc_146_iw1": 126, "asc_44_iw3": 236,
                         "desc_168_iw1": 124, "desc_95_iw3": 123},
    "swath_of_track": {"asc_146_iw1": "IW1", "asc_44_iw3": "IW3",
                       "desc_168_iw1": "IW1", "desc_95_iw3": "IW3"},
}

# The latitude the equirectangular projection is taken at: the anchor's own site point,
# the same latitude its ``metrics`` block used (cos_lat0 = 0.609903).
LAT0 = ANCHOR_EXPECTED["site_point"]["lat"]


# ---------------------------------------------------------------------------
# The committed inputs
# ---------------------------------------------------------------------------
def load_aoi(path=AOI_PATH):
    with open(path) as fh:
        return json.load(fh)


def load_bursts(path=BURSTS_PATH):
    with open(path) as fh:
        return json.load(fh)


def anchor_facts(a):
    """The asserted identity of the committed anchor, read back off its own fields."""
    an = a["anchor"]
    return {
        "osm_way": an["osm_way"],
        "version": an["version"],
        "changeset": an["changeset"],
        "n_distinct_nodes": an["n_distinct_nodes"],
        "kind": an["kind"],
        "span_m": an["span_m"],
        "width_m": an["rectangle_fit"]["width_from_area_perimeter_m"],
        "area_m2": a["metrics"]["flat_area_m2"],
        "span_bearing_deg": an["span_bearing_deg"],
        "wgs84_area_m2": a["metrics"]["wgs84_area_m2"],
        "site_point": a["site_point"],
    }


def deck(a):
    """The structure this file divides by: the anchor's own rectangle, in metres."""
    an = a["anchor"]
    m = a["metrics"]
    return {
        "span_m": an["span_m"],
        "width_m": an["rectangle_fit"]["width_from_area_perimeter_m"],
        "area_m2": m["flat_area_m2"],
        "area_wgs84_m2": m["wgs84_area_m2"],
        "span_bearing_deg": an["span_bearing_deg"],
        "area_source": ("the anchor's shoelace area of the committed outline ring; the "
                        "rectangle with the ring's own area and perimeter is 704.043 x "
                        "35.929 m (rectangle_fit), which is the rectangle used here"),
        "bearing_note": ("the span axis of the ring (its own span_edge_index 6 and the "
                         "two nodes of that edge), as werder_havel_aoi.py derived it: the "
                         "deck enters this file as a bearing and two extents, never as a "
                         "raster"),
        "source": SRC_AOI,
    }


def frames(bursts, expected=BURSTS_EXPECTED):
    """``[(track, date, BurstId, record, n_pol)]`` — one entry per burst frame.

    The catalogue answers the dual-pol pair (VV and VH: two records, one shared
    ``BurstId``) for every frame, so the pair is collapsed here and its two footprints are
    required to be the *same* polygon — if a later catalogue version ever returns two
    geometries under one BurstId, that is raised, not averaged.
    """
    out = []
    for date in sorted(bursts["dates"]):
        e = bursts["dates"][date]
        by_id = {}
        for key, rec in sorted(e["bursts"].items()):
            by_id.setdefault(rec["BurstId"], []).append(rec)
        for bid, recs in sorted(by_id.items()):
            rings = {tuple((round(float(p[0]), 12), round(float(p[1]), 12))
                           for p in bursts_mod.footprint_ring(r)) for r in recs}
            if len(rings) != 1:
                raise SystemExit(f"{date} burst {bid}: the {len(recs)} records of the "
                                 f"dual-pol pair do not share one footprint")
            out.append((e["track"], date, bid, recs[0], len(recs)))
    got = {"n_dates": len({f[1] for f in out}), "n_frames": len(out),
           "n_burst_records": sum(f[4] for f in out),
           "tracks": sorted({f[0] for f in out}),
           "frames_per_track": {t: sum(1 for f in out if f[0] == t)
                                for t in sorted({f[0] for f in out})},
           "swath_of_track": {f[0]: f[3]["SwathIdentifier"] for f in out}}
    problems = [f"bursts.{k}: {got.get(k)!r} != {v!r}"
                for k, v in sorted(expected.items()) if got.get(k) != v]
    return out, problems




# ---------------------------------------------------------------------------
# The grid of the record, from the committed footprints
# ---------------------------------------------------------------------------
def footprint_axes(rec):
    """The footprint's two ground axes, as ``(major, minor, bearing deg)``.

    The **long** axis of an IW burst footprint is *range* (a burst spans one sub-swath
    across range and ~21 km along track), so the axes are separated by extent and not
    guessed: the longest edge of the polygon gives the long direction, collinear pieces of
    that edge are summed (the catalogue may subdivide an edge), and the two axes are the
    vertex spreads along it and across it. That is exact for the parallelogram the
    catalogue outlines, and it avoids the bias a principal-axis fit carries on a ~4:1
    quadrilateral. The bearing returned is that of the **range** axis, i.e. a line.
    """
    xy = [(lon * M_PER_DEG_LON * math.cos(math.radians(LAT0)),
           lat * M_PER_DEG_LAT) for lon, lat in bursts_mod.footprint_ring(rec)]
    n = len(xy)
    edges = [(xy[(i + 1) % n][0] - xy[i][0], xy[(i + 1) % n][1] - xy[i][1])
             for i in range(n)]
    dx, dy = max(edges, key=lambda e: e[0] ** 2 + e[1] ** 2)
    ref = math.hypot(dx, dy)
    if ref == 0.0:
        raise SystemExit(f"degenerate GeoFootprint for {rec['Id']}")
    ux, uy = dx / ref, dy / ref
    ex = ey = 0.0
    for a, b in edges:
        if abs(a * uy - b * ux) <= 1e-3 * ref and a * ux + b * uy > 0.0:
            ex, ey = ex + a, ey + b
    ux, uy = ex / math.hypot(ex, ey), ey / math.hypot(ex, ey)
    proj = [p[0] * ux + p[1] * uy for p in xy]
    perp = [-p[0] * uy + p[1] * ux for p in xy]
    return (max(proj) - min(proj), max(perp) - min(perp),
            math.degrees(math.atan2(ux, uy)) % 180.0)


AXIS_NOTE = ("the two ground axes of a burst are perpendicular: the range axis is "
             "reported as a line (mod 180 deg), the azimuth axis is that line rotated by "
             "90 deg, and which perpendicular sense is the flight direction follows from "
             "the orbit direction of the burst's own record in werder_havel_bursts.json")


def burst_grid(fr):
    """``track -> grid`` of the record: metres per pixel per axis, plus the spread."""
    per = {}
    for track, date, bid, rec, _n in fr:
        major, minor, brg = footprint_axes(rec)
        per.setdefault(track, []).append({
            "range": major / int(rec["SamplesPerBurst"]),
            "azimuth": minor / int(rec["LinesPerBurst"]),
            "bearing": brg, "swath": rec["SwathIdentifier"], "bid": rec["BurstId"],
            "platform": rec["PlatformSerialIdentifier"],
            "lines": int(rec["LinesPerBurst"]),
            "samples": int(rec["SamplesPerBurst"]),
            "major": major, "minor": minor, "date": date,
        })
    out = {}
    for track, gs in sorted(per.items()):
        rng = statistics.median(g["range"] for g in gs)
        azi = statistics.median(g["azimuth"] for g in gs)
        out[track] = {
            "n_frames": len(gs),
            "sub_swath": sorted({g["swath"] for g in gs}),
            "burst_ids": sorted({g["bid"] for g in gs}),
            "platforms": sorted({g["platform"] for g in gs}),
            "range_m_per_px": round(rng, 4),
            "azimuth_m_per_px": round(azi, 4),
            "px_area_m2": round(rng * azi, 3),
            "range_axis_bearing_deg": round(statistics.median(g["bearing"]
                                                              for g in gs), 2),
            "range_m_per_px_spread": [round(min(g["range"] for g in gs), 4),
                                      round(max(g["range"] for g in gs), 4)],
            "azimuth_m_per_px_spread": [round(min(g["azimuth"] for g in gs), 4),
                                        round(max(g["azimuth"] for g in gs), 4)],
            "range_axis_bearing_spread_deg": [round(min(g["bearing"] for g in gs), 2),
                                              round(max(g["bearing"] for g in gs), 2)],
            "footprint_range_extent_m": round(statistics.median(g["major"]
                                                               for g in gs), 1),
            "footprint_azimuth_extent_m": round(statistics.median(g["minor"]
                                                                 for g in gs), 1),
            "lines_per_burst": sorted({g["lines"] for g in gs}),
            "samples_per_burst": sorted({g["samples"] for g in gs}),
        }
    return out


SPREAD_NOTE = ("range_m_per_px_spread / azimuth_m_per_px_spread are over the frames of "
               "the track: the same BurstId comes back with a slightly different footprint "
               "on different dates (the orbit tube), and that variation - not the catalogue "
               "coordinate rounding - is the precision of these metres per pixel")


# ---------------------------------------------------------------------------
# What the grid prices: the deck, in pixels
# ---------------------------------------------------------------------------
def _axis_angle_deg(direction_deg, axis_deg):
    """A direction off an axis, signed and folded into ``[-90, 90)``."""
    return ((direction_deg - axis_deg + 90.0) % 180.0) - 90.0


AREA_NOTE = ("deck_area_px = area_m2 / (range_m_per_px * azimuth_m_per_px) is exact on "
             "these grids: the two axes are perpendicular (AXIS_NOTE) and in the same "
             "plane, so the pixel lattice is orthonormal and the number of pixels a "
             "ground area covers does not depend on its bearing. The *extents* do, which "
             "is why they are projected with the deck's angle off the range axis rather "
             "than divided")


def deck_vs_pixel_grid(grid, d):
    """The deck's extent on each geometry's grid: metres, pixels, and its share."""
    out = {}
    for track, s in sorted(grid.items()):
        ang = _axis_angle_deg(d["span_bearing_deg"], s["range_axis_bearing_deg"])
        rad = math.radians(ang)
        along_range = (d["span_m"] * abs(math.cos(rad)) +
                       d["width_m"] * abs(math.sin(rad)))
        along_azi = (d["span_m"] * abs(math.sin(rad)) +
                     d["width_m"] * abs(math.cos(rad)))
        area_px = d["area_m2"] / s["px_area_m2"]
        out[track] = {
            "deck_vs_range_axis_deg": round(ang, 2),
            "span_along_range_m": round(along_range, 1),
            "span_along_azimuth_m": round(along_azi, 1),
            "span_range_px": round(along_range / s["range_m_per_px"], 1),
            "span_azimuth_px": round(along_azi / s["azimuth_m_per_px"], 1),
            "deck_area_px": round(area_px, 1),
            "deck_frac_of_window": round(area_px / WINDOW_PX, 4),
            "deck_frac_of_masked_decile": round(area_px / MASKED_PX, 4),
            "deck_frac_of_masked_decile_wgs84_area": round(
                d["area_wgs84_m2"] / s["px_area_m2"] / MASKED_PX, 4),
        }
    return out


def deck_vs_window(grid, vs):
    """Whether the record's window can hold the deck, and what would."""
    out = {}
    for track, s in sorted(grid.items()):
        win_rng = WINDOW_SIDE_PX * s["range_m_per_px"]
        win_azi = WINDOW_SIDE_PX * s["azimuth_m_per_px"]
        need_r = vs[track]["span_along_range_m"]
        need_a = vs[track]["span_along_azimuth_m"]
        side_r = 2 * int(math.ceil(need_r / s["range_m_per_px"] / 2.0)) + 1
        side_a = 2 * int(math.ceil(need_a / s["azimuth_m_per_px"] / 2.0)) + 1
        out[track] = {
            "window_range_m": round(win_rng, 1),
            "window_azimuth_m": round(win_azi, 1),
            "window_area_m2": round(WINDOW_PX * s["px_area_m2"], 1),
            "deck_contained": bool(need_r <= win_rng and need_a <= win_azi),
            "deck_span_over_window_range": round(need_r / win_rng, 3),
            "deck_span_over_window_azimuth": round(need_a / win_azi, 3),
            "side_px_to_contain_range": side_r,
            "side_px_to_contain_azimuth": side_a,
            "square_side_px_to_contain": max(side_r, side_a),
        }
    return out




def containment_note(vs):
    """The reading of ``deck_vs_window``, written from its own numbers."""
    r = sorted(v["window_range_m"] for v in vs.values())
    a = sorted(v["window_azimuth_m"] for v in vs.values())
    nr = sorted(v["deck_span_over_window_range"] for v in vs.values())
    na = sorted(v["deck_span_over_window_azimuth"] for v in vs.values())
    side = sorted(v["square_side_px_to_contain"] for v in vs.values())
    return (f"an 80 x 80 window of this record is {r[0]:.0f}-{r[-1]:.0f} m of ground "
            f"across range and {a[0]:.0f}-{a[-1]:.0f} m along azimuth; the deck needs "
            f"{nr[0]:.2f}-{nr[-1]:.2f} of that range span and {na[0]:.2f}-{na[-1]:.2f} of "
            f"that azimuth span, so not one of the four pass geometries contains the "
            f"structure - a square frame of {side[0]}-{side[-1]} px per side at the same "
            f"spacing would. Every mask, amplitude and coherence number of this record is "
            f"therefore measured over the *middle* of the deck; that is what 'containment' "
            f"means in this package, and it is why it is committed as geometry here "
            f"instead of being left to a figure or to a reader's assumption")


def fill_ratio(grid, vs, d):
    """The external 0.86 claim, priced against every pass geometry of the record."""
    implied = d["area_m2"] / (FILL_CLAIM * MASKED_PX)
    per = {}
    for track, s in sorted(grid.items()):
        derived = vs[track]["deck_frac_of_masked_decile"]
        per[track] = {
            "sub_swath": s["sub_swath"],
            "deck_px": vs[track]["deck_area_px"],
            "deck_frac_of_masked_decile": derived,
            "derived_minus_claim": round(derived - FILL_CLAIM, 4),
            "derived_minus_claim_pct": round(100.0 * (derived / FILL_CLAIM - 1.0), 2),
            "implied_px_area_over_derived_pct": round(
                100.0 * (implied / s["px_area_m2"] - 1.0), 2),
            "implied_azimuth_m_per_px": round(implied / s["range_m_per_px"], 4),
            "implied_range_m_per_px": round(implied / s["azimuth_m_per_px"], 4),
        }
    iw3 = [per[t] for t in sorted(per) if grid[t]["sub_swath"] == ["IW3"]]
    iw1 = [per[t] for t in sorted(per) if grid[t]["sub_swath"] == ["IW1"]]
    fine = sorted(v["implied_px_area_over_derived_pct"] for v in iw3)
    low3 = sorted(v["derived_minus_claim_pct"] for v in iw3)
    low1 = sorted(v["derived_minus_claim_pct"] for v in iw1)
    grid3 = sorted(f"{v['implied_azimuth_m_per_px']:.2f} m azimuth / "
                   f"{v['implied_range_m_per_px']:.2f} m range" for v in iw3)
    return {
        "claim": FILL_CLAIM,
        "quote": FILL_CLAIM_QUOTE,
        "source": FILL_CLAIM_SOURCE,
        "reading": ("the deck's ground area in pixels over the exported 640-pixel decile: "
                    "deck_px = area_m2 / (range_m_per_px x azimuth_m_per_px). The decile is "
                    "10 % of the window by definition and the deck is 6.8-8.4 % of the "
                    "window, so the three numbers of the quote close only at one grid"),
        "implied_px_area_m2": round(implied, 3),
        "implied_grid_note": ("what pixel area the claim needs: 550.4 px for a 25 295.7 m2 "
                              "deck, i.e. the derived IW3 range spacing with "
                              f"{grid3[0]} (and {grid3[1]} on the other IW3 geometry) - "
                              "the footprint grid of the record is what the claim is "
                              "priced against, not a quoted product spacing"),
        "per_track": per,
        "verdict": (f"an IW3 reading: the derived ratio is {low3[0]:.2f} and {low3[1]:.2f} "
                    f"% below the claim on the two IW3 geometries and {low1[0]:.2f} and "
                    f"{low1[1]:.2f} % below it on the two IW1 ones. The claim's own grid "
                    f"({implied:.3f} m2/px, {fine[0]:.1f}-{fine[1]:.1f} % finer than the "
                    f"IW3 footprints give) is recorded nowhere - not in this repository "
                    f"and not in the external chain's committed output - so the residual "
                    f"is committed rather than resolved, and both readings stand"),
    }



CORRECTIONS_NOTE = ("three definitions a reader of a fill ratio can get wrong, each priced "
                    "on every geometry of the record: the pixel count of one complex "
                    "sample, the ground-vs-slant range spacing, and the flat-vs-WGS84 "
                    "reading of the deck's own area. None of the three changes the "
                    "committed ratio (the flat, ground, sample reading is the one "
                    "committed); they are committed so that a reader can see how big each "
                    "mistake would be rather than having to guess")


def area_corrections(grid, vs, d):
    """The three definition corrections of this geometry, priced per pass geometry."""
    frac_num, frac_slant = {}, {}
    for track, s in sorted(grid.items()):
        frac_num[track] = round(vs[track]["deck_area_px"] / (2.0 * MASKED_PX), 4)
        frac_slant[track] = round(
            d["area_m2"] / (SLANT_RANGE_M * s["azimuth_m_per_px"]) / MASKED_PX, 4)
    return {
        "complex_sample_not_two_numbers": {
            "note": ("the record's window is 6400 *complex samples*, i.e. 6400 pixels; the "
                     "stored payload writes two numbers per sample (I and Q - 12 800 "
                     "numbers, the ~250 MB of werder_havel_ocv_masks.WINDOW_NOTE), so "
                     "counting the payload's numbers as pixels would double the window and "
                     "its decile. 'px' in this file always means one complex sample."),
            "deck_frac_if_decile_were_1280": frac_num,
        },
        "ground_range_not_slant": {
            "note": ("the derived metres per pixel are *ground* values off the footprint "
                     f"polygons; the SLC's own range sample spacing is {SLANT_RANGE_M} m of "
                     "*slant* range, and slant x azimuth is not a ground area at all. "
                     "Pricing the deck on it would put the ratio above 1 - the mistake "
                     "points the other way and is several times the claim's own residual."),
            "deck_frac_if_slant_range": frac_slant,
        },
        "flat_constants_not_wgs84": {
            "note": ("the anchor commits both readings of the deck's area - "
                     f"{d['area_m2']} m2 on the family's flat constants and "
                     f"{d['area_wgs84_m2']} m2 on the local WGS84 degree lengths, 0.87 % "
                     "apart - and the burst footprints are converted with the flat ones "
                     "(M_PER_DEG_* above), so the flat ratio is the committed one; the "
                     "WGS84 ratio is carried beside it."),
            "deck_frac_if_wgs84_area": {t: v["deck_frac_of_masked_decile_wgs84_area"]
                                        for t, v in sorted(vs.items())},
        },
    }



# ---------------------------------------------------------------------------
# The derived file
# ---------------------------------------------------------------------------
def derive():
    """The whole derived file, from the two committed reference files (offline)."""
    a = load_aoi()
    b = load_bursts()
    fr, bp = frames(b)
    if bp:
        raise SystemExit("the committed burst layer is not the one this file divides "
                         "by:\n  - " + "\n  - ".join(bp))
    an = anchor_facts(a)
    ap = [f"anchor.{k}: {an.get(k)!r} != {v!r}"
          for k, v in sorted(ANCHOR_EXPECTED.items()) if an.get(k) != v]
    if ap:
        raise SystemExit("the committed anchor is not the one this file divides by:\n"
                         "  - " + "\n  - ".join(ap))
    d = deck(a)
    grid = burst_grid(fr)
    vs = deck_vs_pixel_grid(grid, d)
    win = deck_vs_window(grid, vs)
    return {
        "site": SITE,
        "asset_id": masks.ASSET_ID,
        "asset_name": masks.ASSET_NAME,
        "status": ("derived metres per pixel - no pixel layer is claimed; the grid of "
                   "the committed burst footprints, the deck's extent in pixels and the "
                   "containment statement only"),
        "sources": {"anchor": SRC_AOI, "bursts": SRC_BURSTS, "record": SRC_RECORD},
        "window": {
            "side_px": WINDOW_SIDE_PX,
            "px": WINDOW_PX,
            "masked_px": MASKED_PX,
            "masked_frac_of_window": round(MASKED_PX / float(WINDOW_PX), 4),
            "site_point": an["site_point"],
            "note": WINDOW_NOTE,
            "decile_note": ("masked_px is the export's own coherence_masked_pixels column, "
                            "640 in all 991 rows: the brightest decile of the window by "
                            "definition, not a mask over this area and not the echo mask "
                            "of werder_havel_ocv_masks.py (the peak/median rule). "
                            "Everything in this file that compares the deck to 640 "
                            "compares it to this definition"),
        },
        "deck": d,
        "chip_grid": grid,
        "chip_grid_note": SPREAD_NOTE + " " + AXIS_NOTE,
        "deck_vs_pixel_grid": vs,
        "deck_vs_pixel_grid_note": AREA_NOTE,
        "deck_vs_window": win,
        "deck_vs_window_note": containment_note(win),
        "fill_ratio_claim": fill_ratio(grid, vs, d),
        "area_corrections": area_corrections(grid, vs, d),
        "area_corrections_note": CORRECTIONS_NOTE,
        "reference_tracks_note": (
            "the record's four reference chains repeat these 491 dates on the same "
            "relative orbits (ref_asc_146_iw1, ref_desc_168_iw1, refout_asc_146_iw1, "
            "refout_desc_168_iw1 - werder_havel_ocv_masks.RECORD_EXPECTED) but measure "
            "other windows. The export carries no coordinate for any row and the queried "
            "ring does not cover them, so their own metres per pixel is not derivable "
            "here: the grid above is the covering layer's, at the site"),
        "notes": [
            "the deck is divided by the *covering* burst layer, not by the pipeline's own "
            "frames: the export committed no burst id (werder_havel_bursts_resolve.py), "
            "and this file prices that distinction instead of hiding it",
            "the grid is the catalogue's own ground geometry divided by the frame's "
            "committed LinesPerBurst / SamplesPerBurst - no SLC is read and no product "
            "spacing constant is assumed",
            "the 0.86 fill ratio is the external chain's claim, quoted verbatim with its "
            "path; that file is a mutable external input outside this repository, and "
            "only the quote and its pricing are committed here",
            "the anchor's flat constants are the ones the footprints are converted with, "
            "so the flat reading of every area is the committed one and the WGS84 "
            "reading is carried beside it (area_corrections)",
        ],
    }



# ---------------------------------------------------------------------------
# Write / check / main
# ---------------------------------------------------------------------------
def _diff(got, want, path=""):
    """Field-by-field comparison of two JSON trees (floats with a tolerance)."""
    problems = []
    if isinstance(want, dict):
        if not isinstance(got, dict):
            return [f"{path}: not a mapping"]
        for k in sorted(set(got) | set(want)):
            if k not in got:
                problems.append(f"{path}.{k}: missing")
            elif k not in want:
                problems.append(f"{path}.{k}: not in the committed file")
            else:
                problems += _diff(got[k], want[k], f"{path}.{k}")
    elif isinstance(want, list):
        if not isinstance(got, list) or len(got) != len(want):
            problems.append(f"{path}: {got!r} != {want!r}")
        else:
            for j, (a, b) in enumerate(zip(got, want)):
                problems += _diff(a, b, f"{path}[{j}]")
    elif isinstance(want, float):
        if not (isinstance(got, (int, float)) and abs(got - want) <= 1e-9):
            problems.append(f"{path}: {got!r} != {want!r}")
    elif got != want:
        problems.append(f"{path}: {got!r} != {want!r}")
    return problems


def write():
    """Derive the geometry and write the committed reference file."""
    g = derive()
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(g, fh, indent=1, sort_keys=True)
        fh.write("\n")
    fr, w = g["fill_ratio_claim"], g["deck_vs_window"]
    print(f"written: {os.path.relpath(OUT)}\n"
          f"  window {g['window']['side_px']}x{g['window']['side_px']} px = "
          f"{g['window']['px']} complex samples, masked decile "
          f"{g['window']['masked_px']} (definition)\n"
          f"  deck {g['deck']['span_m']} x {g['deck']['width_m']} m = "
          f"{g['deck']['area_m2']} m2 at bearing {g['deck']['span_bearing_deg']} deg\n"
          f"  grid: " + "; ".join(
              f"{t} {s['sub_swath'][0]} {s['range_m_per_px']} m/px range x "
              f"{s['azimuth_m_per_px']} m/px azimuth ({s['px_area_m2']} m2/px, range axis "
              f"at {s['range_axis_bearing_deg']} deg)"
              for t, s in sorted(g["chip_grid"].items())) + "\n"
          f"  deck vs grid: " + "; ".join(
              f"{t} {v['deck_area_px']} px = {v['deck_frac_of_masked_decile']:.4f} of the "
              f"decile" for t, v in sorted(g["deck_vs_pixel_grid"].items())) + "\n"
          f"  claim {fr['claim']}: " + "; ".join(
              f"{t} {v['derived_minus_claim_pct']:+.2f} %"
              for t, v in sorted(fr["per_track"].items())) + "\n"
          f"  contained: " + "; ".join(
              f"{t} {str(v['deck_contained']).lower()} "
              f"({v['square_side_px_to_contain']} px/side would)"
              for t, v in sorted(w.items())))
    return g


def check(committed=None):
    """Re-derive the geometry and compare it with the committed file.

    ``committed`` may be handed in memory - the teeth test tampers with a payload -
    and by default the reference file on disk is read.
    """
    if committed is None:
        with open(OUT) as fh:
            committed = json.load(fh)
    got = derive()
    return _diff(got, committed), got


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="write data/werder_havel/reference/werder_havel_geometry.json")
    ap.add_argument("--json", action="store_true",
                    help="print the derived file as JSON and stop")
    args = ap.parse_args(argv)
    if args.write:
        write()
    if args.json:
        print(json.dumps(derive(), indent=1, sort_keys=True))
        return
    problems, g = check()
    if problems:
        print("deviations:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        raise SystemExit(1)
    frac = sorted(v["deck_frac_of_masked_decile"]
                  for v in g["deck_vs_pixel_grid"].values())
    print(f"ok: {os.path.relpath(OUT)} re-derives field by field - a "
          f"{g['deck']['area_m2']} m2 deck at {g['deck']['span_bearing_deg']} deg over "
          f"{len(g['chip_grid'])} pass geometries, grid " + ", ".join(
              f"{t} {s['range_m_per_px']} x {s['azimuth_m_per_px']} m/px"
              for t, s in sorted(g["chip_grid"].items())) +
          f"; the deck is {frac[0]:.4f}-{frac[-1]:.4f} of the 640-px decile "
          f"(claim {g['fill_ratio_claim']['claim']}); contained: "
          f"{any(v['deck_contained'] for v in g['deck_vs_window'].values())}")


if __name__ == "__main__":
    main()
