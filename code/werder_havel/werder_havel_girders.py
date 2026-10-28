#!/usr/bin/env python3
"""WERDER HAVEL - the eleven spans of the record, as per-track masks.

File: ``data/werder_havel/reference/werder_havel_girders.json`` (derived here, offline, from
the committed anchor, the committed grid and **one** external, quoted fact: the coordinate the
external captures were cropped at).

The question this file answers
------------------------------
``werder_havel_aoi.json`` hands over the structure as 704.043 x 35.929 m (25 295.7 m2) at
bearing 59.59 deg, and ``werder_havel_geometry.json`` prices it in pixels (433.7 / 529.9 /
436.3 / 535.0 px per track, 0.68..0.84 of the exported decile of a window that does not
contain the deck). Both stop before the one thing a **per-girder** layer needs: the deck is
not divided. The record's own unit is the whole structure (one decile per acquisition), so no
channel of the export can ask whether *this opening* changed. The external 400 x 400 captures
can - and the division of the deck into its documented openings is geometry, i.e. this file.

What the mask is
----------------
Eleven **spans** (openings), each a rotated rectangle of the deck's own width in the map
frame, rasterised into each track's rect indices: 2 x 55 m end spans + 9 x 66 m inner spans =
704.0 m, eleven openings (``DIVISION`` carries the source chain). The deck's own committed
edge is 704.043 m, so the 0.043 m residual is carried by the two end spans (55.0215 m each)
while the nine inner spans stay exactly 66.0 m: the division sums to the deck, not to the
map's rounded ``length`` tag.

The masks are **static** - the difference to the family's carola / morandi / deck layers,
where every chip is a fresh window payload: this record's crop rule is one coordinate and one
fixed side, so the same pixel set is the mask on all 491 dates. The *rect* is not quite
constant: five asc_44_iw3 captures are clipped by the burst edge (394..398 of 400 lines) and
are absent from the record (123 captured dates, 118 record rows) - a clipped rect is a
truncated mask, so an extraction must drop it, never stretch it.

The two frames
--------------
*Map*: east/north metres of the family's equirectangular projection about the anchor's site
point (52.4175279, 12.9146744) with the constants the burst footprints are converted with
(``werder_havel_aoi.M_PER_DEG_*``). *Rect*: index (r, c), r = line = azimuth and c = sample =
range, the map point of index (r, c) being

    P(r, c) = CROP_POINT + (c - 200) * range_m_per_px * r_hat
                        + (r - 200) * azimuth_m_per_px * a_hat

with ``r_hat`` the **range direction** (away from the satellite: 78.51 / 80.44 deg on the
ascending geometries, +180 deg on the descending ones) and ``a_hat`` the **flight direction**
(``r_hat`` - 90 deg: northward on the ascending geometries, southward on the descending). The
sense is not assumed here: it is read off the burst layer's own tag
(``werder_havel_bursts.json`` -> ``track_geometries``, ASCENDING/146/IW1 etc.), and the file
prices the 180 deg error it rules out - a mirrored row sense would move the deck's centre by
2 x (200 - r) = 3.4..6.8 px along the rows and flip the sign of its tilt about the range axis,
which is what the ``--place`` gate tests against the captures' own echoes (and what it
currently *fails*: see ``PLACE_STATE_NOTE``).

One external fact, quoted and priced
------------------------------------
The placement needs the coordinate the captures were cropped at. It is **not** in this
repository: it is ``CROP_POINT`` = 52.41734 / 12.91413, quoted from
``/home/projects/havel_werder_analysis/werder_acquisition_census.json`` ("candidate".lat/lon,
whose ``geometry_source`` also names the span division used here) and from the four bridge
runs' ``--lat/--lon`` in ``run_captures.sh`` in the same directory; the crop rule
(``r0 = round(line) - 200``, i.e. the coordinate at index 200, 200) is
``tower/backend/src/input_proxy/copernicus/client.rs::download_burst_rect_ranged``. Both files
are mutable external inputs; what is committed here is the pair, its two sources, the rule, and
everything derived from it. The offset from the crop point to the ring's area centroid is
76.879 m at bearing 58.80 deg, which lands the deck's centre at index

    asc_146_iw1 (201.68, 219.18)   asc_44_iw3  (201.85, 222.90)
    desc_168_iw1 (196.61, 185.00)  desc_95_iw3 (196.70, 181.26)

- rows within 3.5 px of the crop point, columns up to 23 px off it, and the deck (139.5..211.1
px long, 9.0..15.4 px across depending on the track) sits inside all four rects with 85..188
px to the nearest edge. A 10 m error in the external coordinate would move the mask by
1.3..3.0 px (0.03..0.07 of a span): the file therefore quotes the coordinate and its rule,
derives everything from them, and validates the *placed* mask against the captures' own echo
ridges (``--place``) instead of trusting either side.

The conventions, and what would move them
-----------------------------------------
  * **span ids** run 1..11 from the SW end: the ring's edge 6 is node 4147156906 ->
    4147156911, i.e. from (52.4162368, 12.9104842) to (52.4194436, 12.9194436), ascending
    latitude. The map's ring carries no span numbering, so ids 1..11 are this file's convention
    and nothing else's - a per-span channel id is *not* invariant under reversing them, so the
    gate reports a reversal as a *relabelling* (which it is) rather than as a change;
  * **the axis** is the flat-projected span edge (bearing 59.769 deg), not the anchor's
    committed ``span_bearing_deg`` 59.59 deg, which is the *initial great-circle* bearing of
    the same two nodes (``werder_havel_aoi.bearing_deg``). The two differ by 0.180 deg = 2.21 m
    at the far end of the deck = 0.56..0.71 px: small, real, and committed as ``axis_reading``
    rather than folded away. The flat edge is the self-consistent choice because the m/px
    divided by are the flat-projected footprints;
  * **what a wrong division costs** is priced rather than asserted: the partition is also
    rasterised shifted by -66.0 and +66.0 m along the deck axis (the ``variants`` block), which
    is the +-1-span case a channel-level stability test needs - the union is the same deck in
    all three, only the cuts move.

What this file does *not* claim
-------------------------------
``status``: *derived girder masks - geometry only; no echo, no coherence, no OCV is reported
here*. Four limits are committed with the numbers:

  * the **placement is not confirmed**: the crop point is external and quoted, and ``--place``
    measures the placed mask at 0.40..0.90 of its matched nulls' amplitude in the four tracks,
    so everything below is geometry rather than a measurement of an echo (``PLACE_STATE_NOTE``
    keeps the numbers);
  * the division is **documented, not surveyed** - it comes from the source chain above, not
    from a measurement in this repository;
  * the mask is the **map's** rectangle. The family has measured a -81..+7 m (rms 56 m) shape
    error for polygon anchors (``werder_havel_geometry.py``); 5.7 m of it along the axis is
    1.4..1.8 px = 0.08..0.10 of a span, i.e. below the span unit and far below the whole-deck
    one;
  * the **carriageway is not in it**. The record and the map carry one deck, not two: the OSM
    description names two superstructures renewed per carriageway (1995-97 and 2002-04) and the
    map draws their combined 35.929 m. A per-carriageway channel is therefore *not* available
    from this geometry - it cannot be manufactured by halving a rectangle whose internal
    division (touching webs? shared pier line? a gap?) is not in the map. The file says so
    where a reader would look for the per-girder layer's ceiling.

Modes
-----
    --write   derive the girders from the committed anchor and grid and write the file
    --check   re-derive and compare field by field (the default; offline)
    --json    print the derived file as JSON (the numbers the notes above quote)
    --place   the live gate: read one rect line per track from the external captures, test both
              reshape orientations, and require the placed mask to sit on the echo ridge. Not
              confirmed on the rect files read on 2026-10-14 (``PLACE_STATE_NOTE``), so the
              placement stays a quoted external fact rather than a validated one.

Interpreter: ``PYTHON=/usr/bin/python3`` where the PATH python is the numpy-less virtualenv.
"""


from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import werder_havel_aoi as aoi          # noqa: E402
import werder_havel_geometry as geo     # noqa: E402
import werder_havel_ocv_masks as masks  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "werder_havel")
REF = os.path.join(DATA, "reference")
OUT = os.path.join(REF, "werder_havel_girders.json")

BURSTS_PATH = os.path.join(REF, "werder_havel_bursts.json")
AOI_PATH = aoi.OUT
GEOMETRY_PATH = geo.OUT

# ---------------------------------------------------------------------------
# The one external fact: where the captures were cropped
# ---------------------------------------------------------------------------
HARVEST_DIR = "/home/projects/havel_werder_analysis"
# The four *target* tracks' captures, in the harvest directory. They carry the family's rect
# shape ``date|ORBIT|rows|cols|[payload]`` - the same shape the record's own window file has,
# which is why the reader below follows werder_havel_ocv_masks.read_windows.
RECT_FILES = {
    "asc_146_iw1": "werder_rects_asc_146_iw1.txt",
    "asc_44_iw3": "werder_rects_asc_44_iw3.txt",
    "desc_168_iw1": "werder_rects_desc_168_iw1.txt",
    "desc_95_iw3": "werder_rects_desc_95_iw3.txt",
}
RECT_SIDE_PX = 400
CROP_INDEX = RECT_SIDE_PX // 2
CROP_POINT = {"lat": 52.41734, "lon": 12.91413}
CROP_RULE = ("the coordinate sits at index (200, 200): the captures were cut with "
             "r0 = round(line) - 200, c0 = round(pixel) - 200, so the request point is the "
             "rect's centre for any even side")
CROP_RULE_SOURCE = ("tower/backend/src/input_proxy/copernicus/client.rs -> "
                    "download_burst_rect_ranged (a 2 * half + 1 window at half = 200), "
                    "external to this repository")
CROP_POINT_SOURCES = [
    "havel_werder_analysis/werder_acquisition_census.json, \"candidate\".lat 52.41734 / "
    ".lon 12.91413 (its geometry_source also names the 2 x 55 + 9 x 66 m division this file "
    "divides by) - external, mutable",
    "havel_werder_analysis/run_captures.sh, the four bridge runs at --lat 52.41734 "
    "--lon 12.91413 - external, mutable",
]

# The five asc_44_iw3 captures the burst edge clipped (394..398 of 400 lines). They are in the
# rect file and *not* in the record (123 captured dates against 118 record rows), because a
# clipped rect is a truncated mask: excluded by rule, and the rule is committed.
CLIPPED_EXPECTED = {
    "track": "asc_44_iw3",
    "dates": ["2025-05-01", "2025-05-25", "2025-06-30", "2025-07-12", "2025-07-24"],
    "lines": [398, 398, 395, 394, 397],
}

# ---------------------------------------------------------------------------
# The documented division and the variant/null definitions
# ---------------------------------------------------------------------------
SPAN_IDS = [f"span_{i:02d}" for i in range(1, 12)]

# The documented division: 2 x 55 m end spans + 9 x 66 m inner spans = 704.0 m, i.e. eleven
# openings. ``documented_m`` is the sum the reading names; the deck this file divides is
# 704.043 m, so the 0.043 m residual is carried by the two end spans (55.0215 m each) while
# the nine inner spans stay exactly 66.0 m - the division sums to the deck, never to the map's
# rounded ``length`` tag.
DIVISION = {
    "span_ids": SPAN_IDS,
    "n_spans": 11,
    "end_m": 55.0,
    "inner_m": 66.0,
    "documented_m": 704.0,
    "rule": ("the documented span division of the A 10 Havelbruecke Werder: two 55 m end "
             "spans and nine 66 m inner spans, i.e. eleven openings; the residual to the "
             "deck's own 704.043 m outline is split between the two end spans, so the "
             "eleven lengths sum to the deck rather than to the rounded tag"),
    "sources": [
        "havel_werder_analysis/werder_acquisition_census.json, \"candidate\".geometry_source "
        "- external, mutable",
        "data/werder_havel/reference/werder_havel_aoi.json, anchor.length_tag (704 m, the "
        "rounded OSM `length` tag the documented division matches to 0.006 %) - committed",
    ],
    "note": ("the division is documented, not surveyed: the eleven openings come from the "
             "source chain above, not from a measurement in this repository. The "
             "subdivision is this file's own convention - the map's ring carries no span "
             "numbering, so ids 1..11 run from the SW end (edge 6, node 4147156906 -> "
             "4147156911) and nothing else attaches to them"),
}

# The partition-shift variants for the channel-level stability test: the +-1 span case, i.e.
# the whole span lattice moved along the deck axis and clipped to the deck, so that the union
# is the same deck in every variant and only the cuts move.
SHIFT_M = [-66.0, 0.0, 66.0]

# The matched off-deck null: the same mask translated in the rect's index frame, by whole
# pixels, so that its shape and size are identical and only the ground under it moves.
NULL_OFFSETS_PX = [(-6, 0), (6, 0), (0, -6), (0, 6)]
NULL_RULE = ("the same committed mask translated by the offsets above, per rect, clipped to "
             "the rect and re-normalised to its own pixel count: matched in shape and size, "
             "off the deck in ground. The offsets are 3.7..6.0 px of range (22.9..23.0 m) and "
             "6 px of azimuth (91.3..92.7 m) - both well off a 35.929 m deck, and the "
             "extraction must report which is the larger displacement of the two")

# The rasterisation band of the per-track invariant. A rotated rectangle's pixel count is the
# number of pixel *centres* inside it, so it differs from ``area_m2 / px_area_m2`` by the
# boundary, not by the area: the four tracks land at 0.0..1.3 % here (0.0..7.1 px), the most
# tilted grid (asc_44_iw3, -41.9 deg off the range axis) straddling the most partial pixels.
# 2 % is that band with room to spare - a band, not a tolerance on the geometry, which is what
# ``--place`` is for.
RASTER_BAND_FRAC = 0.02


# ---------------------------------------------------------------------------
# The pins
# ---------------------------------------------------------------------------
# The girders divide the *same* anchor and the *same* grid the previous two files committed,
# so their identity is reused rather than restated (the family's ``*_EXPECTED`` pattern), and
# the values the masks were cut against are copied into the committed file's ``inputs`` block:
# a change upstream therefore fails this file twice, once as a moved number and once as a
# moved mask.
ANCHOR_EXPECTED = geo.ANCHOR_EXPECTED
BURSTS_EXPECTED = geo.BURSTS_EXPECTED
RING_IS_A_RECTANGLE = ("the anchor's own rectangle_fit: side chain minus span 0.000 m and end "
                       "edges minus 2 x width 0.000 m - which is what lets a 704.043 m x "
                       "35.929 m rectangle stand in for the ring here")
TRACKS = list(BURSTS_EXPECTED["tracks"])


def load(path):
    with open(path) as fh:
        return json.load(fh)


def orbit_of(track, bursts=None):
    """``ASCENDING`` / ``DESCENDING`` off the burst layer's own ``track_geometries`` tag.

    The sense of the rect's rows and columns is read, not assumed: this is the only committed
    statement of the orbit direction (``werder_havel_bursts.json``), and the 180 deg error it
    rules out is priced in the committed file's ``frame`` block and tested by ``--place``.
    """
    if bursts is None:
        bursts = load(BURSTS_PATH)
    tag = bursts["track_geometries"][track]
    orbit = tag.split("/")[0]
    if orbit not in ("ASCENDING", "DESCENDING"):
        raise SystemExit(f"{track}: track_geometries tag {tag!r} carries no orbit direction")
    return orbit


def local_m(lat, lon, lat0, lon0):
    """East/north metres of the family's equirectangular projection about ``(lat0, lon0)``.

    With the constants of ``werder_havel_aoi`` (``M_PER_DEG_*``), which are the ones the burst
    footprints - and therefore every m/px divided by below - are converted with, so the map
    frame here and the pixel grid are the *same* frame. The 0.66 % / 0.21 % bias against the
    ellipsoid is the anchor's own committed ``METRIC_NOTE``, not re-derived here.
    """
    kx = geo.M_PER_DEG_LON * math.cos(math.radians(lat0))
    return ((lon - lon0) * kx, (lat - lat0) * geo.M_PER_DEG_LAT)


def latlon_of(x, y, lat0, lon0):
    """The inverse of ``local_m`` (the anchor's own round trip)."""
    kx = geo.M_PER_DEG_LON * math.cos(math.radians(lat0))
    return (lat0 + y / geo.M_PER_DEG_LAT, lon0 + x / kx)


def unit_bearing(bearing_deg):
    """Unit vector (east, north) of a compass bearing."""
    rad = math.radians(bearing_deg)
    return (math.sin(rad), math.cos(rad))


def bearing_of(v):
    return math.degrees(math.atan2(v[0], v[1])) % 360.0


def deck_in_map(a, d):
    """The deck in the map frame, re-derived from the committed extract.

    ``d`` is ``werder_havel_geometry.deck``'s own rectangle (the committed ``span_m`` /
    ``width_m`` / ``area_m2`` / ``span_bearing_deg``), reused rather than re-read, so the
    extents this file divides by are *its* committed numbers. Returns the site point, the ring
    in local metres, the area centroid (metres and lat/lon), the axis as a unit vector, and
    **both readings of the axis**: the flat-projected span edge (used below) and the committed
    ``span_bearing_deg`` (a great-circle bearing of the same two nodes), with their difference.
    """
    raw = aoi.load_raw()
    way, pts = aoi.way_geometry(raw)
    coords = {e["id"]: (e["lat"], e["lon"]) for e in raw["elements"]
              if e.get("type") == "node"}
    uniq = list(dict.fromkeys(way["nodes"]))
    site = aoi.node_centroid([coords[n] for n in uniq])
    lat0, lon0 = site
    xy = [local_m(la, lo, lat0, lon0) for la, lo in pts]
    cen = aoi.ring_area_centroid(xy)
    an = a["anchor"]
    isp = an["span_edge_index"]
    axis = aoi.unit(xy[isp], xy[isp + 1])
    flat_deg = bearing_of(axis)
    gc_deg = d["span_bearing_deg"]
    return {
        "site_point": {"lat": lat0, "lon": lon0},
        "ring": [list(p) for p in xy],
        "n_ring_nodes": len(xy),
        "area_centroid_m": list(cen),
        "area_centroid_latlon": list(latlon_of(cen[0], cen[1], lat0, lon0)),
        "axis": list(axis),
        "axis_len_m": math.dist(xy[isp], xy[isp + 1]),
        "flat_bearing_deg": flat_deg,
        "great_circle_bearing_deg": gc_deg,
        "span_bearing_deg": gc_deg,
        "axis_delta_deg": ((flat_deg - gc_deg + 180.0) % 360.0) - 180.0,
        "span_edge_index": isp,
        "span_edge_nodes": list(an["span_edge_nodes"]),
        "span_m": d["span_m"],
        "width_m": d["width_m"],
        "area_m2": d["area_m2"],
        "area_wgs84_m2": d["area_wgs84_m2"],
    }


def rect_frame(track, grid, vs, deck):
    """One track's rect index frame: the two unit vectors, the m/px, and the deck in indices.

    ``range_dir`` points the way the **sample index** grows: the grid's own range axis bearing
    (away from the satellite) on the ascending geometries, that bearing + 180 deg on the
    descending ones. ``azimuth_dir`` points the way the **line index** grows: the flight
    direction, ``range_dir - 90 deg`` (348.51 / 350.44 deg, i.e. northward, ascending;
    191.49 / 189.57 deg, i.e. southward, descending).

    The signed angle between the deck's axis and the range axis is committed here as the flat
    reading (``deck_vs_range_deg``); the committed grid file carries the great-circle one
    (``deck_vs_range_axis_deg``, -18.92 / -20.85 / -41.90 / -39.98), and the two differ by the
    same 0.180 deg ``axis_reading`` prices.
    """
    orbit = orbit_of(track)
    rng_bearing = grid["range_axis_bearing_deg"]
    rng_dir = (rng_bearing + 180.0) % 360.0 if orbit == "DESCENDING" else rng_bearing
    az_dir = (rng_dir - 90.0) % 360.0
    origin = local_m(CROP_POINT["lat"], CROP_POINT["lon"], deck["site_point"]["lat"],
                     deck["site_point"]["lon"])
    delta = (deck["area_centroid_m"][0] - origin[0], deck["area_centroid_m"][1] - origin[1])
    r_hat, a_hat = unit_bearing(rng_dir), unit_bearing(az_dir)
    r_m, a_m = grid["range_m_per_px"], grid["azimuth_m_per_px"]
    centre = (CROP_INDEX + (delta[0] * a_hat[0] + delta[1] * a_hat[1]) / a_m,
              CROP_INDEX + (delta[0] * r_hat[0] + delta[1] * r_hat[1]) / r_m)
    axis = deck["axis"]
    deck_dir = ((axis[0] * a_hat[0] + axis[1] * a_hat[1]) / a_m,
                (axis[0] * r_hat[0] + axis[1] * r_hat[1]) / r_m)
    return {
        "orbit": orbit,
        "range_axis_bearing_deg": rng_bearing,
        "range_dir": list(r_hat),
        "range_dir_bearing_deg": rng_dir,
        "azimuth_dir": list(a_hat),
        "azimuth_dir_bearing_deg": az_dir,
        "range_m_per_px": r_m,
        "azimuth_m_per_px": a_m,
        "px_area_m2": grid["px_area_m2"],
        "origin_m": list(origin),
        "crop_to_centroid_m": list(delta),
        "crop_to_centroid_dist_m": math.hypot(*delta),
        "crop_to_centroid_bearing_deg": bearing_of(delta),
        "centre_index": list(centre),
        "deck_dir_index": list(deck_dir),
        "deck_m_per_pxstep": 1.0 / math.hypot(*deck_dir),
        "deck_vs_range_deg": ((deck["flat_bearing_deg"] - rng_dir + 90.0) % 180.0) - 90.0,
        "deck_area_px_geometry": vs["deck_area_px"],
        "deck_frac_of_masked_decile": vs["deck_frac_of_masked_decile"],
        "deck_vs_range_axis_deg_committed": vs["deck_vs_range_axis_deg"],
        "deck_vs_range_axis_delta_deg": round(
            (((deck["flat_bearing_deg"] - rng_dir + 90.0) % 180.0) - 90.0)
            - vs["deck_vs_range_axis_deg"], 3),
    }


def index_metres(frame, side=RECT_SIDE_PX):
    """The map point (east, north) of every index centre of a ``side x side`` rect.

    Index (r, c) is line r, sample c: the crop's own convention, with the crop point at
    (200, 200) - so the map point of index (r, c) is ``origin + c * range * r_hat +
    r * azimuth * a_hat`` in metres about the origin, with r and c counted from the centre.
    """
    i = np.arange(side, dtype=float) - CROP_INDEX
    rr, cc = np.meshgrid(i, i, indexing="ij")
    r_hat, a_hat = frame["range_dir"], frame["azimuth_dir"]
    east = (frame["origin_m"][0]
            + cc * frame["range_m_per_px"] * r_hat[0]
            + rr * frame["azimuth_m_per_px"] * a_hat[0])
    north = (frame["origin_m"][1]
             + cc * frame["range_m_per_px"] * r_hat[1]
             + rr * frame["azimuth_m_per_px"] * a_hat[1])
    return east, north


def index_of(frame, point):
    """The float ``(row, col)`` index of a map point (the inverse of ``index_metres``)."""
    d = (point[0] - frame["origin_m"][0], point[1] - frame["origin_m"][1])
    row = CROP_INDEX + (d[0] * frame["azimuth_dir"][0] + d[1] * frame["azimuth_dir"][1]) \
        / frame["azimuth_m_per_px"]
    col = CROP_INDEX + (d[0] * frame["range_dir"][0] + d[1] * frame["range_dir"][1]) \
        / frame["range_m_per_px"]
    return (row, col)


# ---------------------------------------------------------------------------
# The deck rectangle, its stations, and the rasterisation
# ---------------------------------------------------------------------------
def division_stations(deck):
    """The documented division on the deck's **own** edge, in metres from its SW end.

    ``2 x 55 + 9 x 66 = 704.0`` m is the documented division and 704.043 m is the deck, so the
    0.043 m residual is split between the two end spans (55.0215 m each) and the nine inner
    spans are exactly 66.0 m. Returns ``(stations, end_m)`` with twelve stations that sum to
    the deck - the map's rounded ``length`` tag is never divided by.
    """
    end = DIVISION["end_m"] + 0.5 * (deck["span_m"] - DIVISION["documented_m"])
    stations = [0.0]
    for i in range(DIVISION["n_spans"]):
        step = end if i in (0, DIVISION["n_spans"] - 1) else DIVISION["inner_m"]
        stations.append(stations[-1] + step)
    return stations, end


def shifted_cuts(stations, span_m, shift):
    """The same lattice moved ``shift`` metres along the axis and clipped to the deck.

    The deck's own ends stay cuts, so the union of the shifted spans is the same deck in every
    variant - only the interior cuts move.
    """
    return ([0.0] + [s + shift for s in stations[1:-1] if 0.0 < s + shift < span_m]
            + [span_m])


def ring_points(deck, frame):
    """The deck's own pixels in one rect: ``(rows, cols, station t, cross-offset u)``.

    The deck is the anchor's rectangle - 704.043 m along its own flat span edge, 35.929 m
    across, centred on the ring's area centroid - placed in the rect's index frame by
    ``index_metres``. ``t`` runs from the SW end of that rectangle to its NE end and ``u`` is
    the signed offset across the axis, so a pixel belongs to exactly one span of any cutting
    that spans ``[0, span_m]``: the partition below is a partition of the pixel set, not an
    approximation of one.
    """
    east, north = index_metres(frame)
    cx, cy = deck["area_centroid_m"]
    dx, dy = deck["axis"]
    nx, ny = -dy, dx
    ex = east - cx
    en = north - cy
    t = ex * dx + en * dy + 0.5 * deck["span_m"]
    u = ex * nx + en * ny
    inside = ((t >= 0.0) & (t <= deck["span_m"])
              & (np.abs(u) <= 0.5 * deck["width_m"]))
    rows, cols = np.nonzero(inside)
    return rows, cols, t[inside], u[inside]


def cut_pixels(points, cuts):
    """Split the deck's pixels at ``cuts`` into row-major pixel lists ``[[row, col], ...]``.

    Half-open on the left and closed on the right, so a pixel on a cut belongs to the *later*
    span and none is lost or double-counted.
    """
    rows, cols, t, _ = points
    out = []
    for j, (a, b) in enumerate(zip(cuts[:-1], cuts[1:])):
        right = (t <= b) if j == len(cuts) - 2 else (t < b)
        sel = (t >= a) & right
        out.append([[int(r), int(c)] for r, c in zip(rows[sel], cols[sel])])
    return out


def span_corners(deck, stations, i, shift=0.0):
    """The four corners of span ``i`` (1-based) in map metres (east, north).

    The rectangle of the span, in the deck's own axes: its corners are what the invariant block
    tests for containment in the 400 x 400 rect, so they are committed rather than left implicit.
    """
    cx, cy = deck["area_centroid_m"]
    dx, dy = deck["axis"]
    nx, ny = -dy, dx
    half_l = 0.5 * deck["span_m"]
    half_w = 0.5 * deck["width_m"]
    cuts = shifted_cuts(stations, deck["span_m"], shift)
    a, b = cuts[i - 1], cuts[i]
    corners = []
    for s, sign in ((a, 1.0), (b, 1.0), (b, -1.0), (a, -1.0)):
        along = s - half_l
        px = cx + along * dx + sign * half_w * nx
        py = cy + along * dy + sign * half_w * ny
        corners.append([px, py])
    return corners


def corners_index(frame, corners):
    """``[[row, col], ...]`` of map-metre corners, rounded to the index they fall in."""
    out = []
    for x, y in corners:
        r, c = index_of(frame, (x, y))
        out.append([int(round(r)), int(round(c))])
    return out







# ---------------------------------------------------------------------------
# The committed file
# ---------------------------------------------------------------------------
SITE = ("Havelbruecke Werder (A 10, Brandenburg) - the deck as its eleven documented "
        "openings, per track, in the captures' own index frame")
SRC_AOI = "data/werder_havel/reference/werder_havel_aoi.json"
SRC_BURSTS = "data/werder_havel/reference/werder_havel_bursts.json"
SRC_GEOMETRY = "data/werder_havel/reference/werder_havel_geometry.json"
SRC_RECORD = "data/werder_havel/werder_havel_measurements_full.txt"

AXIS_READING_NOTE = (
    "the deck's axis enters as the *flat-projected* span edge (bearing 59.769 deg), not as "
    "the anchor's committed span_bearing_deg 59.59 deg, which is the initial great-circle "
    "bearing of the same two nodes (werder_havel_aoi.bearing_deg). The two differ by 0.180 "
    "deg = 2.21 m at the far end of the deck = 0.56..0.71 px = 0.05..0.06 of a span: small, "
    "real, and committed here rather than folded away, because the m/px divided by below are "
    "flat-projected footprints and a per-span channel must not compare two different axes")

STATIC_NOTE = (
    "the masks are static: one crop coordinate and one fixed rect side, so the same pixel set "
    "is the mask on all 491 dates of the record. The *rect* is not quite constant - five "
    "asc_44_iw3 captures are clipped by the burst edge (394..398 of 400 lines) and are absent "
    "from the record (123 captured dates against 118 record rows) - and a clipped rect is a "
    "truncated mask: excluded by rule (inputs.clipped_captures), never stretched")

CEILING_NOTE = (
    "the carriageway is not in it. The record and the map carry one 35.929 m deck, not two: "
    "the OSM description names two superstructures renewed per carriageway (1995-97 and "
    "2002-04) and the map draws their combined width. A per-carriageway channel is therefore "
    "*not* available from this geometry - it cannot be manufactured by halving a rectangle "
    "whose internal division (touching webs? shared pier line? a gap?) is not in the map")

SHAPE_NOTE = (
    "the mask is the map's rectangle. The family has measured a -81..+7 m (rms 56 m) shape "
    "error for polygon anchors (werder_havel_geometry.py): 5.7 m of it along the axis is "
    "1.4..1.8 px = 0.08..0.10 of a span, i.e. below the span unit; across the axis it is what "
    "the matched nulls price rather than a per-span problem")

PLACE_NOTE = (
    "--place is the live gate and the only mode that reads the external captures: it needs "
    "the harvest directory, and it exits non-zero when the placed mask is not the one the "
    "captures' own echoes show. The committed file is derived offline and does not depend on "
    "it")

PLACE_STATE_NOTE = (
    "the placement is priced, not validated. On the rect files read on 2026-10-14 --place "
    "measures an on-mask / matched-null amplitude contrast of 0.40..0.90 in the four tracks - "
    "the deck is *not* the brightest structure in these 400 x 400 rects, the urban land west "
    "of it is - and the mirrored row sense fits two of the four better, so the gate FAILS and "
    "this file's crop point stays a quoted external fact. Everything below is therefore "
    "geometry: the census' coordinate, the documented division, the burst layer's m/px and its "
    "own orbit tag. Nothing in it is a measurement of an echo, and a per-span channel built on "
    "it inherits that limit until --place is green")


def committed_inputs():
    """The identity of the two committed references, asserted, and the deck they hand over.

    The girders reuse ``werder_havel_geometry``'s own reads (``anchor_facts``, ``frames``,
    ``burst_grid``, ``deck_vs_pixel_grid``) so the extents below are *its* committed numbers,
    and they are asserted against the same ``*_EXPECTED`` pins: a moved upstream value fails
    here as well as there. Returns ``(bursts, deck, grid, vs)``.
    """
    a = load(AOI_PATH)
    an = geo.anchor_facts(a)
    problems = [f"anchor.{k}: {an.get(k)!r} != {v!r}"
                for k, v in sorted(ANCHOR_EXPECTED.items()) if an.get(k) != v]
    if a.get("asset_id") != masks.ASSET_ID:
        problems.append(f"anchor.asset_id: {a.get('asset_id')!r} != {masks.ASSET_ID!r}")
    if problems:
        raise SystemExit("the committed anchor is not the one this file divides by:\n  - "
                         + "\n  - ".join(problems))
    b = load(BURSTS_PATH)
    fr, bp = geo.frames(b)
    if bp:
        raise SystemExit("the committed burst layer is not the one this file divides by:\n"
                         "  - " + "\n  - ".join(bp))
    d = geo.deck(a)
    deck = deck_in_map(a, d)
    grid = geo.burst_grid(fr)
    vs = geo.deck_vs_pixel_grid(grid, d)
    return b, deck, grid, vs


def division_block(deck, stations, end):
    """The documented division, its fit to the deck, and the chain the reading came from."""
    lengths = [b - a for a, b in zip(stations[:-1], stations[1:])]
    return {
        "span_ids": SPAN_IDS,
        "n_spans": DIVISION["n_spans"],
        "end_m": DIVISION["end_m"],
        "inner_m": DIVISION["inner_m"],
        "documented_m": DIVISION["documented_m"],
        "deck_m": deck["span_m"],
        "residual_m": round(deck["span_m"] - DIVISION["documented_m"], 4),
        "end_fitted_m": round(end, 4),
        "stations_m": [round(s, 4) for s in stations],
        "lengths_m": [round(x, 4) for x in lengths],
        "sum_m": round(sum(lengths), 4),
        "rule": DIVISION["rule"],
        "sources": DIVISION["sources"],
        "note": DIVISION["note"],
    }


# ---------------------------------------------------------------------------
# The masks: the eleven spans, the shifted variants and the matched nulls
# ---------------------------------------------------------------------------
def span_masks(deck, frame, stations):
    """One track's eleven spans: the pixels, the corners and the extent of each.

    The rectangles are the deck's own axes cut at the twelve stations, so their corners are
    exact (two cuts x the deck's width) and are what the invariant block tests for containment
    in the 400 x 400 rect. They are committed rather than left implicit because a span that
    left the rect would be a silently truncated channel. Returns ``(spans, cut)``, where
    ``cut`` is the same eleven pixel lists in station order.
    """
    cut = cut_pixels(ring_points(deck, frame), stations)
    spans = {}
    for j, sid in enumerate(SPAN_IDS):
        corners = span_corners(deck, stations, j + 1)
        ci = corners_index(frame, corners)
        spans[sid] = {
            "station_m": [round(stations[j], 4), round(stations[j + 1], 4)],
            "length_m": round(stations[j + 1] - stations[j], 4),
            "n_px": len(cut[j]),
            "pixels": cut[j],
            "corners_m": [[round(x, 3), round(y, 3)] for x, y in corners],
            "corners_index": ci,
            "contained_in_rect": bool(all(0 <= r < RECT_SIDE_PX
                                          and 0 <= c < RECT_SIDE_PX for r, c in ci)),
        }
    return spans, cut


def variant_masks(deck, frame, stations):
    """The same span lattice moved ``-66 / 0 / +66`` m along the axis and clipped to the deck.

    This is the +-1-span case a channel-level stability test needs: the deck's own ends stay
    cuts, so the union of a variant's spans is the same deck in all three and only the interior
    cuts move. Only the pixel *counts* are committed - the shifted spans are priced, not used
    as masks.
    """
    out = {}
    for shift in SHIFT_M:
        cuts = shifted_cuts(stations, deck["span_m"], shift)
        px = cut_pixels(ring_points(deck, frame), cuts)
        out[f"shift_{shift:+.0f}m"] = {
            "shift_m": shift,
            "n_spans": len(cuts) - 1,
            "stations_m": [round(s, 4) for s in cuts],
            "n_px": [len(p) for p in px],
            "total_px": sum(len(p) for p in px),
        }
    return out


def matched_nulls(frame, pixels):
    """The committed mask moved by whole pixels: matched in shape, off the deck in ground.

    The offsets are ``NULL_OFFSETS_PX``. What each one means in metres depends on the axes' own
    m/px, so the block reports both displacements and which is larger: the null is only matched
    if the *smaller* of its two still clears the 35.929 m deck.
    """
    out = {}
    for dr, dc in NULL_OFFSETS_PX:
        moved = [[r + dr, c + dc] for r, c in pixels]
        inside = [[r, c] for r, c in moved
                  if 0 <= r < RECT_SIDE_PX and 0 <= c < RECT_SIDE_PX]
        d_rng = abs(dc) * frame["range_m_per_px"]
        d_azi = abs(dr) * frame["azimuth_m_per_px"]
        out[f"dr{dr:+d}_dc{dc:+d}"] = {
            "offset_px": [dr, dc],
            "range_shift_m": round(d_rng, 1),
            "azimuth_shift_m": round(d_azi, 1),
            "n_px": len(moved),
            "n_px_in_rect": len(inside),
            "n_px_clipped": len(moved) - len(inside),
            "larger_displacement": "azimuth" if d_azi > d_rng else "range",
        }
    return out


def derive():
    """The whole committed file, from the two committed reference files (offline)."""
    _b, deck, grid, vs = committed_inputs()
    stations, end = division_stations(deck)
    frames, spans, variants, nulls, invariants = {}, {}, {}, {}, {}
    for track in TRACKS:
        frame = rect_frame(track, grid[track], vs[track], deck)
        frames[track] = frame
        sp, cut = span_masks(deck, frame, stations)
        spans[track] = sp
        variants[track] = variant_masks(deck, frame, stations)
        deck_px = sum(len(p) for p in cut)
        nulls[track] = matched_nulls(frame, [rc for p in cut for rc in p])
        corners = [rc for sid in SPAN_IDS for rc in sp[sid]["corners_index"]]
        invariants[track] = {
            "n_deck_px": deck_px,
            "n_span_px": [sp[sid]["n_px"] for sid in SPAN_IDS],
            "partition_covers_deck": sum(sp[sid]["n_px"] for sid in SPAN_IDS) == deck_px,
            "smallest_span_px": min(sp[sid]["n_px"] for sid in SPAN_IDS),
            "largest_span_px": max(sp[sid]["n_px"] for sid in SPAN_IDS),
            "all_corners_in_rect": all(sp[sid]["contained_in_rect"] for sid in SPAN_IDS),
            "min_corner_margin_px": min(min(r, c, RECT_SIDE_PX - 1 - r,
                                            RECT_SIDE_PX - 1 - c) for r, c in corners),
            "stations_sum_to_deck": abs(stations[-1] - deck["span_m"]) < 1e-9,
            "deck_px_vs_geometry": {
                "mask_px": deck_px,
                "geometry_px": vs[track]["deck_area_px"],
                "delta_px": round(deck_px - vs[track]["deck_area_px"], 1),
            },
        }
    return {
        "site": SITE,
        "asset_id": masks.ASSET_ID,
        "asset_name": masks.ASSET_NAME,
        "status": ("derived girder masks - geometry only; no echo, no coherence, no OCV is "
                   "reported here"),
        "sources": {"anchor": SRC_AOI, "bursts": SRC_BURSTS, "geometry": SRC_GEOMETRY,
                    "record": SRC_RECORD},
        "inputs": {
            "crop": {
                "point": CROP_POINT,
                "index": [CROP_INDEX, CROP_INDEX],
                "side_px": RECT_SIDE_PX,
                "rule": CROP_RULE,
                "rule_source": CROP_RULE_SOURCE,
                "sources": CROP_POINT_SOURCES,
                "external": True,
            },
            "deck": deck,
            "axis_reading": {
                "used": "the flat-projected span edge",
                "flat_bearing_deg": deck["flat_bearing_deg"],
                "great_circle_bearing_deg": deck["great_circle_bearing_deg"],
                "delta_deg": deck["axis_delta_deg"],
                "note": AXIS_READING_NOTE,
            },
            "division": division_block(deck, stations, end),
            "grid": {t: grid[t] for t in TRACKS},
            "vs_pixel_grid": {t: vs[t] for t in TRACKS},
            "pins": {"anchor": ANCHOR_EXPECTED, "bursts": BURSTS_EXPECTED},
            "clipped_captures": CLIPPED_EXPECTED,
        },
        "frames": frames,
        "spans": spans,
        "variants": variants,
        "nulls": nulls,
        "invariants": invariants,
        "notes": [STATIC_NOTE, DIVISION["note"], CEILING_NOTE, SHAPE_NOTE, PLACE_NOTE,
                  PLACE_STATE_NOTE],
    }


def write():
    """Derive the girder masks and write the committed file."""
    g = derive()
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(g, fh, indent=1, sort_keys=True)
        fh.write("\n")
    d = g["inputs"]["division"]
    print(f"written: {os.path.relpath(OUT)}\n"
          f"  {d['n_spans']} spans: "
          + " + ".join(f"{x}" for x in d["lengths_m"])
          + f" m = {d['sum_m']} m (deck {d['deck_m']} m, residual {d['residual_m']} m "
            f"in the two end spans)\n"
          f"  masks: " + "; ".join(
              f"{t} {g['invariants'][t]['n_deck_px']} px in 11 spans "
              f"({g['invariants'][t]['smallest_span_px']}.."
              f"{g['invariants'][t]['largest_span_px']} px per span)"
              for t in TRACKS) + "\n"
          f"  centres: " + "; ".join(
              f"{t} ({g['frames'][t]['centre_index'][0]:.2f}, "
              f"{g['frames'][t]['centre_index'][1]:.2f})" for t in TRACKS) + "\n"
          f"  corners inside the {RECT_SIDE_PX} px rect: "
          + ", ".join(f"{t} {g['invariants'][t]['min_corner_margin_px']} px to the nearest "
                      f"edge" for t in TRACKS))
    return g


# ---------------------------------------------------------------------------
# Check
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


def check(committed=None):
    """Re-derive the masks and compare them with the committed file (offline).

    ``committed`` may be handed in memory - the teeth test tampers with a payload - and by
    default the reference file on disk is read.
    """
    if committed is None:
        if not os.path.exists(OUT):
            return [f"{os.path.relpath(OUT)}: missing - run --write"], None
        with open(OUT) as fh:
            committed = json.load(fh)
    got = derive()
    problems = _diff(got, committed, "girders")

    # The readings below are re-asserted on the *re-derived* file, so they hold even when the
    # committed copy has been tampered with: a partition that no longer covers the deck, a span
    # that left the rect, or a division that stopped summing to the deck is a problem whatever
    # is on disk.
    d = got["inputs"]["division"]
    if abs(d["sum_m"] - d["deck_m"]) > 1e-6:
        problems.append(f"division: the eleven lengths sum to {d['sum_m']} m, the deck is "
                        f"{d['deck_m']} m")
    if abs(d["deck_m"] - d["documented_m"]) > 0.05:
        problems.append(f"division: the deck moved {abs(d['deck_m'] - d['documented_m']):.3f} "
                        f"m off the documented {d['documented_m']} m - re-read the reading")
    for track in TRACKS:
        inv = got["invariants"][track]
        if not inv["partition_covers_deck"]:
            problems.append(f"invariants.{track}: the spans do not cover the deck exactly")
        if not inv["all_corners_in_rect"]:
            problems.append(f"invariants.{track}: a span corner left the {RECT_SIDE_PX} px "
                            f"rect")
        if not inv["stations_sum_to_deck"]:
            problems.append(f"invariants.{track}: the stations do not sum to the deck")
        if abs(inv["deck_px_vs_geometry"]["delta_px"]) > RASTER_BAND_FRAC * \
                inv["deck_px_vs_geometry"]["geometry_px"]:
            problems.append(f"invariants.{track}: the rasterised deck is "
                            f"{inv['deck_px_vs_geometry']['delta_px']} px off the geometry's "
                            f"own area in pixels - beyond the "
                            f"{RASTER_BAND_FRAC:.0%} rasterisation band ("
                            f"{RASTER_BAND_FRAC * inv['deck_px_vs_geometry']['geometry_px']:.1f}"
                            f" px here)")
    return problems, got


# ---------------------------------------------------------------------------
# The live gate: --place
# ---------------------------------------------------------------------------
def rect_rows(track):
    """Yield ``(date, orbit, rows, cols, payload)`` of one track's external rect file.

    The shape is the family's ``date|ORBIT|rows|cols|[payload]`` - the same shape the record's
    own window file has, which is why this reader follows
    ``werder_havel_ocv_masks.read_windows`` - and the reader is also its acceptance test:
    malformed lines are skipped here and counted by the caller.
    """
    path = os.path.join(HARVEST_DIR, RECT_FILES[track])
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("|", 4)
            if len(parts) != 5:
                continue
            yield parts[0], parts[1], int(parts[2]), int(parts[3]), parts[4]


def amplitude(payload, rows, cols):
    """``|z|`` of one capture's payload, reshaped ``(lines, samples)`` by the rect's own shape.

    Row 0 is the first line the crop rule put at index 200: which *ground* that is depends on
    the row sense, which is exactly what ``place`` tests both ways round rather than assuming.
    """
    arr = json.loads(payload)
    z = np.asarray([complex(a[0], a[1]) for a in arr], dtype=complex)
    return np.abs(z.reshape(rows, cols))


def _mean_amplitude(a, pixels, rows, cols):
    """The mean amplitude over a pixel list, with what leaves the rect dropped (never filled).

    Returns ``(mean, n_px_in_rect)``.
    """
    inside = [[r, c] for r, c in pixels if 0 <= r < rows and 0 <= c < cols]
    if not inside:
        return float("nan"), 0
    return float(np.mean([a[r, c] for r, c in inside])), len(inside)


def bright_centroid(a):
    """The amplitude-weighted centre of a rect: a robust reading of where its echo ridge sits.

    The weights are the amplitude above the rect's own median, so the dark water contributes
    nothing and the bright structure - deck, banks, buildings - carries the centroid.
    """
    w = a - np.median(a)
    w = np.where(w > 0.0, w, 0.0)
    if w.sum() <= 0.0:
        return (float("nan"), float("nan"))
    rr = np.arange(a.shape[0], dtype=float)[:, None]
    cc = np.arange(a.shape[1], dtype=float)[None, :]
    return (float((w * rr).sum() / w.sum()), float((w * cc).sum() / w.sum()))


def place():
    """The live gate: the placed mask against the captures' own echoes, both row senses.

    One rect line per track is read from the external harvest directory, and three things are
    required of it: the clipped-capture rule of ``CLIPPED_EXPECTED`` holds; the placed mask is
    brighter than all four matched off-deck nulls of ``NULL_OFFSETS_PX``; and the committed
    row sense beats the mirrored one, whose deck centre sits 2 x (200 - r) px away and whose
    tilt about the range axis has the opposite sign. The mirror test is the 180 deg error the
    docstring prices - and it is what decides the row sense from the data instead of from the
    orbit tag alone.
    """
    if not os.path.isdir(HARVEST_DIR):
        raise SystemExit(f"{HARVEST_DIR}: the external captures are not there - --place "
                         f"reads the rect files, the offline modes do not")
    _b, deck, grid, vs = committed_inputs()
    stations, _end = division_stations(deck)
    problems = []
    print(f"placed mask vs the captures' own echoes ({HARVEST_DIR}):")
    for track in TRACKS:
        rows = list(rect_rows(track))
        clipped = [(d, n) for d, _o, n, _c, _p in rows if n != RECT_SIDE_PX]
        if track == CLIPPED_EXPECTED["track"]:
            got = ([d for d, _n in clipped], [n for _d, n in clipped])
            want = (CLIPPED_EXPECTED["dates"], CLIPPED_EXPECTED["lines"])
            if got != want:
                problems.append(f"{track}: the clipped captures read "
                                f"{list(zip(*got))}, CLIPPED_EXPECTED says {list(zip(*want))}")
        elif clipped:
            problems.append(f"{track}: {len(clipped)} clipped capture(s) "
                            f"{[d for d, _n in clipped]} in a track CLIPPED_EXPECTED calls "
                            f"whole")
        full = next((row for row in rows
                     if row[2] == RECT_SIDE_PX and row[3] == RECT_SIDE_PX), None)
        if full is None:
            problems.append(f"{track}: no {RECT_SIDE_PX} x {RECT_SIDE_PX} capture in "
                            f"{RECT_FILES[track]}")
            continue
        date, _orbit, nr, nc, payload = full
        a = amplitude(payload, nr, nc)
        frame = rect_frame(track, grid[track], vs[track], deck)
        _spans, cut = span_masks(deck, frame, stations)
        mask = [rc for p in cut for rc in p]
        centre = frame["centre_index"]
        cen = bright_centroid(a)
        res = {}
        for name, mk in (("committed", mask),
                         ("mirrored", [[nr - 1 - r, c] for r, c in mask])):
            on_mean, n_in = _mean_amplitude(a, mk, nr, nc)
            offs = [_mean_amplitude(a, [[r + dr, c + dc] for r, c in mk], nr, nc)[0]
                    for dr, dc in NULL_OFFSETS_PX]
            off_max = max(offs)
            res[name] = {
                "on_mean": on_mean,
                "n_px_in_rect": n_in,
                "off_max": off_max,
                "off_mean": sum(offs) / len(offs),
                "contrast": on_mean / off_max if off_max > 0.0 else float("inf"),
                "centre_row": centre[0] if name == "committed" else nr - 1 - centre[0],
            }
        c, m = res["committed"], res["mirrored"]
        d_cen = abs(cen[0] - c["centre_row"])
        d_mir = abs(cen[0] - m["centre_row"])
        if c["contrast"] <= 1.0:
            problems.append(f"{track}: the placed mask is not brighter than its matched "
                            f"nulls (contrast {c['contrast']:.3f}) - the crop coordinate or "
                            f"the frame is wrong")
        if c["contrast"] < m["contrast"]:
            problems.append(f"{track}: the mirrored row sense fits the echoes better "
                            f"({m['contrast']:.3f} against {c['contrast']:.3f}) - "
                            f"azimuth_dir is the wrong way round")
        if d_cen > d_mir:
            problems.append(f"{track}: the bright centroid ({cen[0]:.1f}, {cen[1]:.1f}) is "
                            f"{d_cen:.1f} px from the committed deck centre row and "
                            f"{d_mir:.1f} px from the mirrored one - the row sense is mirrored")
        print(f"  {track} {date}: on-mask {c['on_mean']:.4f} vs matched null "
              f"{c['off_max']:.4f} (contrast x{c['contrast']:.3f}, mirrored "
              f"x{m['contrast']:.3f}); {c['n_px_in_rect']} px; deck centre "
              f"({centre[0]:.2f}, {centre[1]:.2f}) vs bright centroid "
              f"({cen[0]:.1f}, {cen[1]:.1f}) = {d_cen:.1f} px (mirrored {d_mir:.1f} px)")
    if problems:
        print(f"PLACE FAILED - {len(problems)} problem(s)", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        raise SystemExit(1)
    print(f"ok: --place confirms the placement on all {len(TRACKS)} tracks - the committed "
          f"masks sit on the captures' own echo ridge and the mirrored row sense fits none of "
          f"them better")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="derive the girder masks and write the committed file")
    ap.add_argument("--check", action="store_true",
                    help="re-derive and compare field by field (the default; offline)")
    ap.add_argument("--json", action="store_true",
                    help="print the derived file as JSON and stop")
    ap.add_argument("--place", action="store_true",
                    help="the live gate: the placed mask against the external captures")
    args = ap.parse_args(argv)
    if args.place:
        return place()
    if args.write:
        write()
    if args.json:
        print(json.dumps(derive(), indent=1, sort_keys=True))
        return 0
    problems, g = check()
    if problems:
        print("deviations:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        raise SystemExit(1)
    d, inv = g["inputs"]["division"], g["invariants"]
    ax = g["inputs"]["axis_reading"]
    print(f"ok: {os.path.relpath(OUT)} re-derives field by field - {d['n_spans']} spans of "
          + "/".join(f"{x}" for x in d["lengths_m"])
          + f" m summing to the deck's {d['deck_m']} m (residual {d['residual_m']} m in the "
            f"end spans); "
          + "; ".join(f"{t} {inv[t]['n_deck_px']} px "
                      f"({inv[t]['smallest_span_px']}..{inv[t]['largest_span_px']} px per "
                      f"span, corners {inv[t]['min_corner_margin_px']} px inside the rect)"
                      for t in TRACKS)
          + f"; axis {ax['flat_bearing_deg']:.3f} deg flat vs "
            f"{ax['great_circle_bearing_deg']} deg committed "
            f"({ax['delta_deg']:+.3f} deg)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
