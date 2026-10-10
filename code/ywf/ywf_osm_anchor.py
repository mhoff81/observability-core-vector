#!/usr/bin/env python3
"""YWF — the OSM anchor of the site: which road the tower fell on, and what that can support.

Files: ``data/ywf/reference/ywf_osm_road.json`` (the committed map extract) and
``data/ywf/reference/ywf_geometry.json`` (the derived anchor, this module).

The question this file answers
------------------------------
Bautzen carries a *deck-edge* layer — a band of brightness across the bridge deck,
anchored to a rectified pixel grid (``range_m_per_px`` 3.37,
``azimuth_m_per_px`` 12.38, an 81-px ROI). It exists because Bautzen has a
committed **rect spec**: the export resolves the deck's pixels, so an edge is a
geometric object *of the record*.

YWF has none of that. Its record is a **7 x 7** complex window per acquisition
around one anchor point (49 px), the export's own ``segment_latitude`` /
``segment_longitude`` columns are **empty for all five mast sections**, and a
window carries no pixel-spacing spec. A band or edge test needs at least a
resolvable ground extent per pixel *along* the band axis. So this module commits
what *is* known and states the limit:

  * the site point (``onboarder.assets.lat/lon``) stands **36.85 m from** the
    carriageway centreline of OSM way ``577440814`` (``highway=tertiary``,
    ``surface=asphalt``) — the public road the collapsed Unit 21 monopole fell
    onto (``distance_from_asset_point_m``, below). The point is the record's
    *measurement* anchor, not a point on the road: the export does not say where
    the tower struck it;
  * the anchor is taken on the road's **carriageway centreline**, never on an
    area/polygon: the Bergsøysund analysis measured a −81..+7 m shape error
    (rms 56 m) when it anchored to a polygon instead;
  * the burst grid of the record is **estimated** from the committed burst
    footprints (metres per pixel, per axis) — no pixel spacing is invented;
  * the road is then compared to that grid as a **number**, not as a layer.

The estimate cross-checks against the family: for IW2/IW3 it returns 3.35-3.45 m
per pixel **across range**, against the 3.37 m of Bautzen's committed rect spec,
and 14.4-15.3 m per pixel **across azimuth** against Bautzen's 12.38 m.

What this module does *not* claim
---------------------------------
``status``: *anchor reference — no edge layer is claimed*, and this file commits
**why** in two independent steps (``road_vs_window``). First, **containment**: the
near carriageway edge lies **8.10 px (ASC/IW3) and 9.33 px (DESC/IW2) across
range** from the window centre, while a 7 x 7 frame reaches **6 px at most** (3 px
from its centre) — the road is *outside* the frame, so not one committed pixel
crosses the carriageway and there is nothing to run a band test on. Second,
**resolution**: even a frame that contained it would still have to resolve the
~7 m carriageway, which is 2.03-2.09 px across range but **0.46-0.49 px** across
azimuth. Containment alone needs ~19-21 px per side (~27-29 px with Bautzen's
+/-4-row ROI) — 361-841 samples, 7.4x to 17.2x the committed payload of 49 px,
because the cost is quadratic in the frame side. The anchor's payoff
therefore stays gated on the **window size**, a backend constant
(``TOWER_WINDOW_SIZE = 7`` in ``insar_monitor.rs``): that is the staged step C of
``code/ywf/README.md``, not this file.

Modes
-----
    --fetch   (re)fetch the way from the public OSM API and write the extract [network]
    --write   derive the geometry from the committed extract and write it
    --check   re-derive and compare field by field (the default)
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

import ywf_ocv_masks as masks  # noqa: E402
import ywf_bursts_resolve as bursts_mod  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "ywf")
RAW = os.path.join(DATA, "reference", "ywf_osm_road.json")
OUT = os.path.join(DATA, "reference", "ywf_geometry.json")

OSM_API = "https://api.openstreetmap.org/api/0.6/way/{wid}/full.json"
LICENSE = ("\u00a9 OpenStreetMap contributors, ODbL 1.0 "
           "(https://www.openstreetmap.org/copyright)")

WAY_ID = 577440814
# The structural identity of the committed extract (it cannot be re-derived from
# the record, so it is asserted — Carola's ``OSM_EXPECTED`` pattern).
OSM_EXPECTED = {
    "way_id": WAY_ID,
    "version": 8,
    "tags": {"highway": "tertiary", "surface": "asphalt"},
    "n_nodes": 346,
}

SITE = "Yeongdeok Wind Farm (Samgye-ri)"
ASSET_ID = "af683c7e-416a-4395-86cc-094ee83d9497"
ASSET_NAME = "Yeongdeok Wind Turbine (Samgye-ri)"
ASSET_LAT, ASSET_LON = 36.42356, 129.42031

# The carriageway width is an **assumption**, flagged as such in the file: the
# way carries no ``width`` tag. 7.0 m is the documented two-lane carriageway
# figure used elsewhere in this family (Carola's B 170 deck is 18 m with 4 lanes
# + tram; a Korean tertiary road between fields is the narrow case).
CARRIAGEWAY_WIDTH_M = 7.0
CARRIAGEWAY_WIDTH_NOTE = ("assumption, not a survey: the way carries no width "
                          "tag; a two-lane carriageway is taken as 7.0 m")

M_PER_DEG_LAT = 110540.0
M_PER_DEG_LON = 111320.0


# ---------------------------------------------------------------------------
# The committed map extract
# ---------------------------------------------------------------------------
def fetch():
    """Fetch the way from the public OSM API and write the raw extract."""
    import requests

    url = OSM_API.format(wid=WAY_ID)
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    d = r.json()
    if "error" in d:
        raise SystemExit(f"{url}: {d['error']}")
    d["_source"] = url
    d["_license"] = LICENSE
    os.makedirs(os.path.dirname(os.path.abspath(RAW)), exist_ok=True)
    with open(RAW, "w") as fh:
        json.dump(d, fh, indent=1, sort_keys=True)
        fh.write("\n")
    ways = [e for e in d["elements"] if e.get("type") == "way"]
    print(f"written: {os.path.relpath(RAW)} — way {ways[0]['id']} "
          f"v{ways[0]['version']}, {len(ways[0]['nodes'])} nodes, "
          f"{len(d['elements']) - len(ways)} of them with coordinates")
    return d


def load_raw(path=RAW):
    with open(path) as fh:
        return json.load(fh)


def way_geometry(raw):
    """``(way, [(lat, lon), ...])`` of the committed extract, in node order."""
    ways = [e for e in raw["elements"] if e.get("type") == "way"]
    if len(ways) != 1:
        raise SystemExit(f"{len(ways)} ways in the extract, expected 1")
    way = ways[0]
    coords = {e["id"]: (e["lat"], e["lon"]) for e in raw["elements"]
              if e.get("type") == "node"}
    pts = [coords[n] for n in way["nodes"] if n in coords]
    if len(pts) != len(way["nodes"]):
        raise SystemExit("the extract holds a way node without coordinates: "
                         "fetch again with /full.json")
    return way, pts


# ---------------------------------------------------------------------------
# Planar geometry (metres), local to the site point
# ---------------------------------------------------------------------------
def _local(lat, lon, lat0=ASSET_LAT, lon0=ASSET_LON):
    """Local planar metres east/north of ``(lat0, lon0)`` (equirectangular).

    Adequate at the scale of this question (a < 12 km road and a < 1 km anchor
    neighbourhood): the error of the flat approximation over 12 km at 36.4 N is
    far below the ~7 m carriageway width that is being compared.
    """
    return ((lon - lon0) * M_PER_DEG_LON * math.cos(math.radians(lat0)),
            (lat - lat0) * M_PER_DEG_LAT)


def bearing_deg(a, b):
    """Initial bearing a -> b, degrees clockwise from north (platform degrees)."""
    lat1, lon1 = math.radians(a[0]), math.radians(a[1])
    lat2, lon2 = math.radians(b[0]), math.radians(b[1])
    y = math.sin(lon2 - lon1) * math.cos(lat2)
    x = (math.cos(lat1) * math.sin(lat2) -
         math.sin(lat1) * math.cos(lat2) * math.cos(lon2 - lon1))
    return math.degrees(math.atan2(y, x)) % 360.0


def nearest_point(pts, p=(ASSET_LAT, ASSET_LON)):
    """Closest point of the polyline to ``p``, as metres and an index.

    Returns ``(distance_m, i, (lat, lon), bearing_deg)`` where ``i`` is the index
    of the segment's first node and the bearing is that segment's own, i.e. the
    local direction of the **carriageway centreline** at the anchor.
    """
    best = None
    px, py = _local(*p)
    for i in range(len(pts) - 1):
        ax, ay = _local(*pts[i])
        bx, by = _local(*pts[i + 1])
        dx, dy = bx - ax, by - ay
        seg2 = dx * dx + dy * dy
        t = 0.0 if seg2 == 0 else max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / seg2))
        qx, qy = ax + t * dx, ay + t * dy
        d = math.hypot(px - qx, py - qy)
        if best is None or d < best[0]:
            lat = pts[i][0] + t * (pts[i + 1][0] - pts[i][0])
            lon = pts[i][1] + t * (pts[i + 1][1] - pts[i][1])
            best = (d, i, (lat, lon), bearing_deg(pts[i], pts[i + 1]))
    return best


def polyline_length_m(pts):
    return sum(math.dist(_local(*pts[i]), _local(*pts[i + 1]))
               for i in range(len(pts) - 1))


# ---------------------------------------------------------------------------
# The pixel grid of this record (estimated from the committed burst footprints)
# ---------------------------------------------------------------------------
def _footprint_polygon(rec):
    """``[(lon, lat), ...]`` of the burst's ``GeoFootprint``, closure removed.

    The catalogue returns the footprint as a GeoJSON polygon
    (``{"type": "Polygon", "coordinates": [[[lon, lat], ...]]}``); anything else
    is refused rather than guessed.
    """
    g = rec["GeoFootprint"]
    if not (isinstance(g, dict) and g.get("type") == "Polygon"
            and g.get("coordinates")):
        raise SystemExit(f"unexpected GeoFootprint for {rec['Id']}: "
                         f"{type(g).__name__} — the resolver writes the "
                         f"catalogue's own GeoJSON polygon")
    ring = [(float(p[0]), float(p[1])) for p in g["coordinates"][0]]
    return ring[:-1] if len(ring) > 1 and ring[0] == ring[-1] else ring


def footprint_axes(rec):
    """The footprint's two ground axes, as ``(major, minor, bearing°)``.

    The **long** axis of an IW burst footprint is *range* (a burst covers one
    sub-swath across range, ~21 km along track), so the axes are separated by
    extent and not guessed: the longest edge of the polygon gives the long
    direction, collinear pieces of that edge are summed (the catalogue may
    subdivide an edge), and the two axes are the vertex spreads along it and
    across it. That is exact for the rectangle the catalogue outlines, and it
    avoids the bias a principal-axis fit carries on a ~4:1 parallelogram. The
    bearing returned is that of the **range** axis, i.e. a line (see
    ``AXIS_NOTE``).
    """
    xy = [(lon * M_PER_DEG_LON * math.cos(math.radians(ASSET_LAT)),
           lat * M_PER_DEG_LAT) for lon, lat in _footprint_polygon(rec)]
    n = len(xy)
    edges = [(xy[(i + 1) % n][0] - xy[i][0], xy[(i + 1) % n][1] - xy[i][1])
             for i in range(n)]
    dx, dy = max(edges, key=lambda e: e[0] ** 2 + e[1] ** 2)
    ref = math.hypot(dx, dy)
    if ref == 0.0:
        raise SystemExit(f"degenerate GeoFootprint for {rec['Id']}")
    ux, uy = dx / ref, dy / ref
    # sum the pieces of the long edge (collinear to within ~1e-3 of its length)
    ex = ey = 0.0
    for a, b in edges:
        if abs(a * uy - b * ux) <= 1e-3 * ref and a * ux + b * uy > 0.0:
            ex, ey = ex + a, ey + b
    ux, uy = ex / math.hypot(ex, ey), ey / math.hypot(ex, ey)
    proj = [p[0] * ux + p[1] * uy for p in xy]
    perp = [-p[0] * uy + p[1] * ux for p in xy]
    return (max(proj) - min(proj), max(perp) - min(perp),
            math.degrees(math.atan2(ux, uy)) % 180.0)


def burst_grid(bursts):
    """``orbit -> metres per pixel per axis``, from the 33 committed footprints."""
    per_orbit = {}
    for rec in bursts["bursts"].values():
        major, minor, bearing = footprint_axes(rec)
        per_orbit.setdefault(rec["OrbitDirection"], []).append({
            "range": major / rec["SamplesPerBurst"],
            "azimuth": minor / rec["LinesPerBurst"],
            "bearing": bearing,
            "swath": rec["SwathIdentifier"],
            "lines": rec["LinesPerBurst"],
            "samples": rec["SamplesPerBurst"],
        })
    out = {}
    for orbit, gs in sorted(per_orbit.items()):
        out[orbit] = {
            "n_bursts": len(gs),
            "sub_swath": sorted({g["swath"] for g in gs}),
            "range_m_per_px": round(statistics.median(g["range"] for g in gs), 3),
            "azimuth_m_per_px": round(statistics.median(g["azimuth"]
                                                        for g in gs), 3),
            "range_axis_bearing_deg": round(statistics.median(g["bearing"]
                                                              for g in gs), 2),
            "lines_per_burst": sorted({g["lines"] for g in gs}),
            "samples_per_burst": sorted({g["samples"] for g in gs}),
        }
    return out


AXIS_NOTE = ("the two ground axes of a burst are perpendicular: the range axis is "
             "reported as a line (mod 180°), the azimuth axis is that line rotated "
             "by 90°, and which of the two perpendicular senses is the flight "
             "direction follows from the orbit direction reported per chip in "
             "ywf_bursts.json")

# The window of this record, and what the export does *not* carry.
WINDOW_PX = 49
WINDOW_NOTE = ("one 7x7 complex window per acquisition at the asset point; the "
               "export committed no per-section coordinates (its "
               "segment_latitude / segment_longitude columns are empty for all "
               "five mast sections), so this package has no second anchor and no "
               "rect spec")

# The frame of that window as a pixel *index* box, and the two facts the
# containment question below turns on: a 7x7 frame reaches (7 - 1) px from one of
# its pixels to the opposite one, i.e. 3 px from its centre when the anchor sits
# at the centre. ``ROI_HALF_PX`` is Bautzen's own ROI half-width (its rect spec
# drives a band of +/-4 rows around its anchor row); it is carried here only to
# price the comparison, not because this record has such an ROI.
WINDOW_SIDE_PX = masks.WINDOW_SIZE           # 7 — the backend's TOWER_WINDOW_SIZE
WINDOW_REACH_PX = float(WINDOW_SIDE_PX - 1)  # 6 px at most, 3 px from the centre
ROI_HALF_PX = 4


def road_offset_px(offset_m, offset_bearing_deg, g):
    """The two pixel components of a ground offset, on one orbit's grid.

    ``g`` is one ``chip_grid`` entry (metres per pixel per axis and the bearing
    of the range axis). The burst's two ground axes are perpendicular, so a
    ground vector of ``offset_m`` metres at ``offset_bearing_deg`` decomposes
    into them by projection: its range component over ``range_m_per_px``, its
    azimuth component (the other axis) over ``azimuth_m_per_px``. Only the
    magnitudes are returned: which sense of an axis is positive depends on the
    pass geometry, and the question here — does the frame reach the carriageway?
    — is a question about extent.
    """
    da = math.radians(offset_bearing_deg - g["range_axis_bearing_deg"])
    return (abs(offset_m * math.cos(da)) / g["range_m_per_px"],
            abs(offset_m * math.sin(da)) / g["azimuth_m_per_px"])


def containment_note(vs_window):
    """The reading of ``road_vs_window``, written from its own numbers."""
    near = sorted(v["near_edge_range_px"] for v in vs_window.values())
    xsr = sorted(v["cross_section_range_px"] for v in vs_window.values())
    xsa = sorted(v["cross_section_azimuth_px"] for v in vs_window.values())
    side = sorted(v["side_px_to_contain"] for v in vs_window.values())
    side_roi = sorted(v["side_px_with_roi"] for v in vs_window.values())
    return (f"the road lies outside the frame: the near carriageway edge is "
            f"{near[0]:.2f}-{near[-1]:.2f} px across range from the window "
            f"centre, while a {WINDOW_SIDE_PX}x{WINDOW_SIDE_PX} frame reaches "
            f"{WINDOW_REACH_PX:.0f} px at most ({WINDOW_REACH_PX / 2.0:.0f} px "
            f"from its centre), so no committed pixel crosses the carriageway - "
            f"and a frame that did contain it would still have to resolve the "
            f"{CARRIAGEWAY_WIDTH_M:.1f} m cross-section ({xsr[0]:.2f}-{xsr[-1]:.2f} "
            f"px across range, {xsa[0]:.2f}-{xsa[-1]:.2f} px across azimuth). "
            f"Containment needs {side[0]}-{side[-1]} px per side "
            f"({side_roi[0]}-{side_roi[-1]} px with Bautzen's +/-4-row ROI): "
            f"a frame of {side[0]} px per side is {side[0] ** 2} samples "
            f"({side[0] ** 2 / WINDOW_PX:.1f}x the committed {WINDOW_PX} px) and "
            f"{side_roi[-1]} px per side {side_roi[-1] ** 2} "
            f"({side_roi[-1] ** 2 / WINDOW_PX:.1f}x) — the cost of containment "
            f"is quadratic in the side, and this is the containment half of "
            f"blocker 4 of code/ywf/README.md")


# ---------------------------------------------------------------------------
# The committed anchor
# ---------------------------------------------------------------------------
def derive(raw=None, bursts=None):
    """The anchor file ``data/ywf/reference/ywf_geometry.json``, from committed data."""
    raw = load_raw() if raw is None else raw
    if bursts is None:
        with open(bursts_mod.OUT) as fh:
            bursts = json.load(fh)
    way, pts = way_geometry(raw)
    dist, i, near, brg = nearest_point(pts)
    grid = burst_grid(bursts)
    road = {}
    for orbit, g in sorted(grid.items()):
        road[orbit] = {
            "carriageway_width_range_px": round(
                CARRIAGEWAY_WIDTH_M / g["range_m_per_px"], 2),
            "carriageway_width_azimuth_px": round(
                CARRIAGEWAY_WIDTH_M / g["azimuth_m_per_px"], 2),
        }
    # The containment half of ``status``: the same ground facts measured inside
    # the window's own frame. The near carriageway edge is the anchor's own
    # offset reduced by half the carriageway, taken at the anchor-to-road
    # bearing (the road's normal and that direction differ by well under a
    # degree, i.e. far below the 0.01 px this file rounds to); the cross-section
    # is the carriageway width along that normal, the centreline bearing + 90°.
    near_edge_m = dist - CARRIAGEWAY_WIDTH_M / 2.0
    offset_brg = bearing_deg((ASSET_LAT, ASSET_LON), near)
    vs_window = {}
    for orbit, g in sorted(grid.items()):
        r_off, a_off = road_offset_px(dist, offset_brg, g)
        r_near, a_near = road_offset_px(near_edge_m, offset_brg, g)
        r_xs, a_xs = road_offset_px(CARRIAGEWAY_WIDTH_M, brg + 90.0, g)
        side = 2 * int(math.ceil(r_near)) + 1
        vs_window[orbit] = {
            "road_offset_range_px": round(r_off, 2),
            "road_offset_azimuth_px": round(a_off, 2),
            "near_edge_range_px": round(r_near, 2),
            "near_edge_azimuth_px": round(a_near, 2),
            "contained": r_near <= WINDOW_REACH_PX,
            "outside_by_px": round(r_near - WINDOW_REACH_PX, 2),
            "side_px_to_contain": side,
            "side_px_with_roi": side + 2 * ROI_HALF_PX,
            "cross_section_range_px": round(r_xs, 2),
            "cross_section_azimuth_px": round(a_xs, 2),
        }
    return {
        "site": SITE,
        "asset_id": ASSET_ID,
        "asset_name": ASSET_NAME,
        "asset_point": {"lat": ASSET_LAT, "lon": ASSET_LON},
        "road": {
            "osm_way": int(way["id"]),
            "version": way.get("version"),
            "timestamp": way.get("timestamp"),
            "tags": dict(sorted(way.get("tags", {}).items())),
            "n_nodes": len(way["nodes"]),
            "centreline_length_km": round(polyline_length_m(pts) / 1000.0, 3),
            "nearest_point": {"lat": round(near[0], 7), "lon": round(near[1], 7)},
            "nearest_segment_node_index": i,
            "distance_from_asset_point_m": round(dist, 2),
            "centreline_bearing_deg": round(brg, 2),
            "carriageway_width_m": CARRIAGEWAY_WIDTH_M,
            "carriageway_width_source": CARRIAGEWAY_WIDTH_NOTE,
            "anchor": ("carriageway centreline of the OSM way — never an "
                       "area/polygon (the Bergsøysund analysis measured a "
                       "−81..+7 m shape error, rms 56 m, for a polygon anchor)"),
        },
        "chip_grid": grid,
        "chip_grid_axes_note": AXIS_NOTE,
        "window": {"px": WINDOW_PX, "note": WINDOW_NOTE},
        "road_vs_pixel_grid": road,
        "road_vs_window": {
            "anchor_offset_m": round(dist, 2),
            "anchor_offset_bearing_deg": round(offset_brg, 2),
            "near_edge_offset_m": round(near_edge_m, 2),
            "window_side_px": WINDOW_SIDE_PX,
            "window_px": WINDOW_PX,
            "window_half_extent_px": (WINDOW_SIDE_PX - 1) / 2.0,
            "window_reach_px": WINDOW_REACH_PX,
            "roi_half_px": ROI_HALF_PX,
            "per_orbit": vs_window,
            "note": containment_note(vs_window),
        },
        "sources": {
            "osm_way": OSM_API.format(wid=WAY_ID),
            "osm_extract": "data/ywf/reference/ywf_osm_road.json",
            "bursts": "data/ywf/reference/ywf_bursts.json",
            "license": LICENSE,
        },
        "status": ("anchor reference — no edge layer is claimed: a 7x7 window per "
                   "acquisition, no rect/pixel-spacing spec and no per-section "
                   "coordinates (see code/ywf/README.md)"),
        "note": ("this anchor supports a statement about *where the tower stood and "
                 "which road it fell on*, not a band/edge layer: such a layer needs "
                 "a resolvable ground extent per pixel along the band axis, and "
                 "``road_vs_pixel_grid`` is the number that says how far this "
                 "record is from that"),
    }


def write():
    """Derive the anchor from the committed extract and write it."""
    g = derive()
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(g, fh, indent=1, sort_keys=True)
        fh.write("\n")
    r = g["road"]
    print(f"written: {os.path.relpath(OUT)}\n"
          f"  way {r['osm_way']} v{r['version']}  {r['tags']}  "
          f"{r['n_nodes']} nodes  {r['centreline_length_km']} km\n"
          f"  asset point is {r['distance_from_asset_point_m']} m from the "
          f"centreline, which runs at {r['centreline_bearing_deg']}° there\n"
          f"  grid: " + "; ".join(
              f"{o} {s['sub_swath']} {s['range_m_per_px']} m/px range, "
              f"{s['azimuth_m_per_px']} m/px azimuth (range axis at "
              f"{s['range_axis_bearing_deg']}°)"
              for o, s in sorted(g["chip_grid"].items())) + "\n"
          f"  road vs grid: " + "; ".join(
              f"{o} {v['carriageway_width_range_px']} px across range, "
              f"{v['carriageway_width_azimuth_px']} px across azimuth"
              for o, v in sorted(g["road_vs_pixel_grid"].items()))
          + f"\n  road vs window: " + "; ".join(
              f"{o} near edge {v['near_edge_range_px']} px across range "
              f"(contained: {str(v['contained']).lower()}; "
              f"{v['side_px_to_contain']} px/side would contain it)"
              for o, v in sorted(g["road_vs_window"]["per_orbit"].items())))
    return g


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


def check():
    """Re-derive the anchor and compare it with the committed file."""
    with open(OUT) as fh:
        committed = json.load(fh)
    raw = load_raw()
    way, _ = way_geometry(raw)
    facts = {"way_id": int(way["id"]), "version": way.get("version"),
             "tags": dict(way.get("tags", {})), "n_nodes": len(way["nodes"])}
    problems = [f"osm.{k}: {facts.get(k)!r} != {v!r}"
                for k, v in sorted(OSM_EXPECTED.items()) if facts.get(k) != v]
    problems += _diff(derive(raw), committed)
    return problems, facts


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch", action="store_true",
                    help="(re)fetch the way from the public OSM API")
    ap.add_argument("--write", action="store_true",
                    help="write data/ywf/reference/ywf_geometry.json")
    args = ap.parse_args(argv)
    if args.fetch:
        fetch()
    if args.write:
        write()
    problems, facts = check()
    if problems:
        print("deviations:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        raise SystemExit(1)
    print(f"ok: {os.path.relpath(RAW)} is way {facts['way_id']} "
          f"v{facts['version']} ({facts['tags']}, {facts['n_nodes']} nodes) and "
          f"{os.path.relpath(OUT)} re-derives field by field")


if __name__ == "__main__":
    main()

