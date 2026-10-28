#!/usr/bin/env python3
"""HALTERN - the OSM anchor of the site: the turbine's own node, and what it can support.

Files: ``data/haltern/reference/haltern_osm_node.json`` (the committed map extract)
and ``data/haltern/reference/haltern_geometry.json`` (the derived anchor, this module).

The question this file answers
------------------------------
Bautzen carries a *deck-edge* layer - a band of brightness across the bridge deck,
anchored to a rectified pixel grid (``range_m_per_px`` 3.37,
``azimuth_m_per_px`` 12.38, an 81-px ROI). It exists because Bautzen has a
committed **rect spec**: the export resolves the deck's pixels, so an edge is a
geometric object *of the record*.

Haltern has none of that. Its record is a **7 x 7** complex window per acquisition
around one anchor point (49 px), the export's own ``segment_latitude`` /
``segment_longitude`` columns are **empty for all five mast sections**, and a
window carries no pixel-spacing spec. A band or edge test needs at least a
resolvable ground extent per pixel *along* the band axis. So this module commits
what *is* known and states the limit:

  * the anchor is the **turbine's own OpenStreetMap node** ``8268709968``
    (``generator:source=wind``, ``model=N149``, ``manufacturer=Nordex``,
    ``height=164``, ``generator:output:electricity=4.5 MW``,
    ``operator=RAG Montan Immobilien, Stadtwerke Haltern``). The tower did not
    fall onto a named road, so there is no carriageway to anchor to. It is
    anchored to the structure itself, which is
    a **point** object and therefore the strongest anchor of this family: it
    cannot carry the polygon shape error the Bergsoysund analysis measured
    (-81..+7 m, rms 56 m);
  * the record's asset point (``onboarder.assets.lat/lon``) is the *same* point:
    ``distance_from_asset_point_m`` is the number that says so, and it is the
    check that the pipeline measured this turbine and not a neighbour;
  * the **event date is corroborated by the map itself**. The node's third
    version - changeset ``111912198``, ``2021-09-30T11:23:44Z``, i.e. **within a
    day of the collapse** - added ``damage:date=2021-09-29`` and
    ``destroyed=yes``. This module commits that version, its changeset and its
    timestamp, so the collapse date of this record is an *externally verifiable*
    fact rather than a quoted one;
  * the burst grid of the record is **estimated** from the committed burst
    footprints (metres per pixel, per axis) - no pixel spacing is invented;
  * the structure is then compared to that grid as a **number**, not as a layer.

The estimate cross-checks against the family: for IW1/IW3 it returns 3.3-4.6 m
per pixel **across range**, against the 3.37 m of Bautzen's committed rect spec,
and 14-16 m per pixel **across azimuth** against Bautzen's 12.38 m.

What this module does *not* claim
---------------------------------
``status``: *anchor reference - no edge layer is claimed*, and this file commits
**why** in two independent steps (``structure_vs_window``). First,
**containment**: the structure's ground extent lies **outside** the 7 x 7 frame,
which reaches **6 px at most** (3 px from its centre), while the tower base alone
already spans more than that on both axes - so not one committed pixel spans the
structure and there is nothing to run a band test on. Second, **resolution**:
even a frame that contained it must resolve the feature it wants to test. The
**tower base** (9.0 m) is the resolvable one; the **rotor** (149 m) is not, and
neither is anything smaller than the base. The anchor's payoff therefore stays
gated on the **window size**, a backend constant (``TOWER_WINDOW_SIZE = 7`` in
``insar_monitor.rs``): that is the staged step C of ``code/haltern/README.md``,
not this file.

Modes
-----
    --fetch   (re)fetch the node and its full version history from the public OSM
              API and write the extract [network]
    --write   derive the geometry from the committed extract and write it
    (default) re-derive and compare field by field
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

import haltern_ocv_masks as masks  # noqa: E402
import haltern_bursts_resolve as bursts_mod  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "haltern")
RAW = os.path.join(DATA, "reference", "haltern_osm_node.json")
OUT = os.path.join(DATA, "reference", "haltern_geometry.json")

OSM_API = "https://api.openstreetmap.org/api/0.6/node/{nid}.json"
OSM_HISTORY_API = "https://api.openstreetmap.org/api/0.6/node/{nid}/history.json"
LICENSE = ("\u00a9 OpenStreetMap contributors, ODbL 1.0 "
           "(https://www.openstreetmap.org/copyright)")

NODE_ID = 8268709968
# The structural identity of the committed extract (it cannot be re-derived from
# the record, so it is asserted - Carola's ``OSM_EXPECTED`` pattern). The version
# count, the event changeset and the event date are the facts the package's whole
# event claim rests on, so they are pinned here and not read off the file.
OSM_EXPECTED = {
    "node_id": NODE_ID,
    "n_versions": 6,
    "n_versions_before_event": 2,
    "event_changeset": 111912198,
    "event_edit_ts": "2021-09-30T11:23:44Z",
    "damage_date": "2021-09-29",
    "destroyed": "yes",
    "manufacturer": "Nordex",
    "model": "N149",
    "height_m": 164.0,
    "output_mw": 4.5,
    "operator": "RAG Montan Immobilien, Stadtwerke Haltern",
}

SITE = "Haltern AV9 Wind Park (Kreis Recklinghausen, NRW)"
ASSET_ID = "2eadf4d5-5e53-4e8b-a311-fac0ed9872bb"
ASSET_NAME = "Haltern AV9 Wind Turbine (Nordex N149/4.5)"
ASSET_LAT, ASSET_LON = 51.7445808, 7.0907269

# The two ground extents of the structure. The **base** is the record's own
# committed geometry (``structure.width_m`` of the onboarded asset, the base
# outer diameter of the hybrid tower); the **rotor** is the N149's swept
# diameter - the machine's defining dimension, and the one a coarser window
# would have to resolve if it ever wanted to see a rotor rather than a mast.
TOWER_BASE_WIDTH_M = 9.0
ROTOR_DIAMETER_M = 149.0
EXTENT_SOURCE = ("base = the committed structure.width_m of the onboarded asset; "
                 "rotor = the 149 m swept diameter of the model tag of the OSM "
                 "node (N149)")

M_PER_DEG_LAT = 110540.0
M_PER_DEG_LON = 111320.0


# ---------------------------------------------------------------------------
# The committed map extract
# ---------------------------------------------------------------------------
def fetch():
    """Fetch the node and its full version history; write the raw extract.

    Two public OSM calls are committed together, because the anchor needs both
    halves: ``history.json`` carries the *event* (the version that added
    ``damage:date``, and the changeset and timestamp of that edit), ``node.json``
    the node's current state. The file is written verbatim, with each element's
    own ``_source`` URL, so ``--check`` can assert what was fetched.
    """
    import requests

    calls = (("history", OSM_HISTORY_API.format(nid=NODE_ID)),
             ("current", OSM_API.format(nid=NODE_ID)))
    out = {"_license": LICENSE, "_node_id": NODE_ID}
    for key, url in calls:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        d = r.json()
        if "error" in d:
            raise SystemExit(f"{url}: {d['error']}")
        d["_source"] = url
        out[key] = d
    out["_sources"] = [u for _, u in calls]
    os.makedirs(os.path.dirname(os.path.abspath(RAW)), exist_ok=True)
    with open(RAW, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write("\n")
    cur = node_current(out)
    tags = cur.get("tags", {})
    print(f"written: {os.path.relpath(RAW)} - node {NODE_ID} v{cur['version']} "
          f"({len(node_history(out))} versions), {cur['lat']}, {cur['lon']}, "
          f"{tags.get('manufacturer')} {tags.get('model')} "
          f"({tags.get('height')} m), damage:date {tags.get('damage:date')}")
    return out


def load_raw(path=RAW):
    with open(path) as fh:
        return json.load(fh)


def node_history(raw):
    """Every committed version of the node, oldest first."""
    els = [e for e in raw["history"]["elements"] if e.get("type") == "node"]
    if not els:
        raise SystemExit("the extract holds no node history")
    return sorted(els, key=lambda e: e["version"])


def node_current(raw):
    """The node's current committed version (the ``/node/{id}.json`` element)."""
    els = [e for e in raw["current"]["elements"] if e.get("type") == "node"]
    if len(els) != 1:
        raise SystemExit(f"{len(els)} nodes in the current extract, expected 1")
    return els[0]


def event_version(hist):
    """The version that introduced ``damage:date`` - the event's own edit.

    This is the file's load-bearing fact: the collapse date of this record is
    read off the map's own edit history, at the version and changeset it names,
    not quoted from the package.
    """
    for e in hist:
        if "damage:date" in e.get("tags", {}):
            return e
    raise SystemExit("no version of the node carries damage:date")


def version_edit(hist, version):
    """``(changeset, timestamp, user)`` of one version."""
    for e in hist:
        if e["version"] == version:
            return {"changeset": e.get("changeset"),
                    "timestamp": e.get("timestamp"), "user": e.get("user")}
    raise SystemExit(f"version {version} is not in the extract")


# ---------------------------------------------------------------------------
# Planar geometry (metres), local to the site point
# ---------------------------------------------------------------------------
def _local(lat, lon, lat0=ASSET_LAT, lon0=ASSET_LON):
    """Local planar metres east/north of ``(lat0, lon0)`` (equirectangular).

    Adequate at the scale of this question (an anchor-to-node distance of metres
    to a few hundred metres): the error of the flat approximation over 1 km at
    51.74 N is far below the ~3-16 m that one committed pixel spans.
    """
    return ((lon - lon0) * M_PER_DEG_LON * math.cos(math.radians(lat0)),
            (lat - lat0) * M_PER_DEG_LAT)


def gps_distance_m(a, b):
    """Great-circle distance between two ``(lat, lon)`` pairs, in metres."""
    return math.hypot(*_local(b[0], b[1], a[0], a[1]))


def bearing_deg(a, b):
    """Initial bearing a -> b, degrees clockwise from north (platform degrees)."""
    lat1, lon1 = math.radians(a[0]), math.radians(a[1])
    lat2, lon2 = math.radians(b[0]), math.radians(b[1])
    y = math.sin(lon2 - lon1) * math.cos(lat2)
    x = (math.cos(lat1) * math.sin(lat2) -
         math.sin(lat1) * math.cos(lat2) * math.cos(lon2 - lon1))
    return math.degrees(math.atan2(y, x)) % 360.0


# ---------------------------------------------------------------------------
# The burst grid of the record, from the committed footprints
# ---------------------------------------------------------------------------
def _footprint_polygon(rec):
    """``[(lon, lat), ...]`` of a catalogue ``GeoFootprint`` (its outer ring)."""
    g = rec["GeoFootprint"]
    if g.get("type") != "Polygon":
        raise SystemExit(f"{rec['Id']}: GeoFootprint is {g.get('type')}, expected "
                         f"the catalogue's own GeoJSON polygon")
    ring = [(float(p[0]), float(p[1])) for p in g["coordinates"][0]]
    return ring[:-1] if len(ring) > 1 and ring[0] == ring[-1] else ring


def footprint_axes(rec):
    """The footprint's two ground axes, as ``(major, minor, bearing deg)``.

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
    """``orbit -> metres per pixel per axis``, from the committed footprints."""
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
             "reported as a line (mod 180 deg), the azimuth axis is that line "
             "rotated by 90 deg, and which of the two perpendicular senses is the "
             "flight direction follows from the orbit direction reported per chip "
             "in haltern_bursts.json")


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
WINDOW_SIDE_PX = masks.WINDOW_SIZE           # 7 - the backend's TOWER_WINDOW_SIZE
WINDOW_REACH_PX = float(WINDOW_SIDE_PX - 1)  # 6 px at most, 3 px from the centre
ROI_HALF_PX = 4


def extent_px(width_m, g):
    """A ground extent's two pixel components on one orbit's grid.

    The extent is an *isotropic* one (the tower's base is a circle), so it is
    simply divided by each axis' metres per pixel - no projection and no bearing
    enters, which is why the containment statement itself needs no
    ``range_axis_bearing_deg``.
    """
    return (width_m / g["range_m_per_px"], width_m / g["azimuth_m_per_px"])


def _ts_seconds(iso):
    """Unix seconds of an ISO-8601 ``...Z`` timestamp (stdlib only, no tz db)."""
    import calendar
    import datetime
    dt = datetime.datetime.strptime(iso, "%Y-%m-%dT%H:%M:%SZ")
    return calendar.timegm(dt.timetuple())


def containment_note(vs):
    """The reading of ``structure_vs_window``, written from its own numbers."""
    side = sorted(v["side_px_to_contain"] for v in vs.values())
    side_roi = sorted(v["side_px_with_roi"] for v in vs.values())
    base_r = sorted(v["base_width_range_px"] for v in vs.values())
    base_a = sorted(v["base_width_azimuth_px"] for v in vs.values())
    rot_r = sorted(v["rotor_diameter_range_px"] for v in vs.values())
    return (f"the structure lies outside the frame: the tower base alone is "
            f"{base_r[0]:.1f}-{base_r[-1]:.1f} px across range and "
            f"{base_a[0]:.1f}-{base_a[-1]:.1f} px across azimuth, while a "
            f"{WINDOW_SIDE_PX}x{WINDOW_SIDE_PX} frame reaches {WINDOW_REACH_PX:.0f} "
            f"px at most ({WINDOW_REACH_PX / 2.0:.0f} px from its centre) - so not "
            f"one committed pixel spans the base, the base is contained on no "
            f"orbit, and there is nothing to run a band test on. Containment needs "
            f"{side[0]}-{side[-1]} px per side ({side_roi[0]}-{side_roi[-1]} px with "
            f"Bautzen's +/-4-row ROI): a frame of {side[0]} px per side is "
            f"{side[0] ** 2} samples ({side[0] ** 2 / WINDOW_PX:.1f}x the committed "
            f"{WINDOW_PX} px) and {side_roi[-1]} px per side {side_roi[-1] ** 2} "
            f"({side_roi[-1] ** 2 / WINDOW_PX:.1f}x) - the cost is quadratic in the "
            f"side. The N149 rotor, {ROTOR_DIAMETER_M:.0f} m = "
            f"{rot_r[0]:.1f}-{rot_r[-1]:.1f} px across range, *would* span such a "
            f"frame - but a rotor is not a band/edge object: its ground projection "
            f"is a line whose length and direction are set by the yaw, which "
            f"neither the export nor the catalogue carries. That is the second half "
            f"of blocker 4 of code/haltern/README.md: containment is purchasable "
            f"with a wider window, a *testable* edge is not.")


# ---------------------------------------------------------------------------
# The committed anchor
# ---------------------------------------------------------------------------
def derive(raw=None, bursts=None):
    """The anchor file ``data/haltern/reference/haltern_geometry.json``, from committed data."""
    raw = load_raw() if raw is None else raw
    if bursts is None:
        with open(bursts_mod.OUT) as fh:
            bursts = json.load(fh)
    hist = node_history(raw)
    cur = node_current(raw)
    ev = event_version(hist)
    tags = dict(sorted(cur.get("tags", {}).items()))
    node = (cur["lat"], cur["lon"])
    dist = gps_distance_m((ASSET_LAT, ASSET_LON), node)
    brg = bearing_deg((ASSET_LAT, ASSET_LON), node)
    grid = burst_grid(bursts)
    structure = {}
    for orbit, g in sorted(grid.items()):
        r_base, a_base = extent_px(TOWER_BASE_WIDTH_M, g)
        r_rot, a_rot = extent_px(ROTOR_DIAMETER_M, g)
        structure[orbit] = {
            "base_width_range_px": round(r_base, 2),
            "base_width_azimuth_px": round(a_base, 2),
            "rotor_diameter_range_px": round(r_rot, 1),
            "rotor_diameter_azimuth_px": round(a_rot, 1),
        }
    # The containment half of ``status``: the same ground facts measured inside
    # the window's own frame. An isotropic extent of ``w`` px centred on the
    # anchor needs an odd frame of ``2 * ceil(w / 2) + 1`` px per side; the worse
    # of the two axes decides, because the frame is square.
    vs_window = {}
    for orbit, g in sorted(grid.items()):
        r_base, a_base = extent_px(TOWER_BASE_WIDTH_M, g)
        r_rot, a_rot = extent_px(ROTOR_DIAMETER_M, g)
        base = max(r_base, a_base)
        rotor = max(r_rot, a_rot)
        side = 2 * int(math.ceil(base / 2.0)) + 1
        vs_window[orbit] = {
            "base_width_range_px": round(r_base, 2),
            "base_width_azimuth_px": round(a_base, 2),
            "base_contained": bool(base <= WINDOW_REACH_PX),
            "base_outside_by_px": round(base - WINDOW_REACH_PX, 2),
            "side_px_to_contain": side,
            "side_px_with_roi": side + 2 * ROI_HALF_PX,
            "rotor_diameter_range_px": round(r_rot, 1),
            "rotor_diameter_azimuth_px": round(a_rot, 1),
            "rotor_side_px_to_contain": 2 * int(math.ceil(rotor / 2.0)) + 1,
        }
    return {
        "site": SITE,
        "asset_id": ASSET_ID,
        "asset_name": ASSET_NAME,
        "asset_point": {"lat": ASSET_LAT, "lon": ASSET_LON},
        "turbine": {
            "osm_node": NODE_ID,
            "version": cur.get("version"),
            "timestamp": cur.get("timestamp"),
            "n_versions": len(hist),
            "tags": tags,
            "manufacturer": tags.get("manufacturer"),
            "model": tags.get("model"),
            "height_m": float(tags["height"]) if "height" in tags else None,
            "generator_output_mw": 4.5,
            "operator": tags.get("operator"),
            "destroyed": tags.get("destroyed"),
            "damage_date": tags.get("damage:date"),
            "anchor": ("the node itself - the turbine's own OSM object, never an "
                       "area/polygon (the Bergsoysund analysis measured a "
                       "-81..+7 m shape error, rms 56 m, for a polygon anchor)"),
            "distance_from_asset_point_m": round(dist, 2),
            "bearing_from_asset_point_deg": round(brg, 2),
        },
        "event_corroboration": {
            "damage_date": ev["tags"]["damage:date"],
            "destroyed": ev["tags"].get("destroyed"),
            "recorded_by": version_edit(hist, ev["version"]),
            "recorded_at_version": ev["version"],
            "n_versions_before_event": ev["version"] - 1,
            "lag_days": round((_ts_seconds(ev.get("timestamp")) -
                               _ts_seconds("2021-09-29T00:00:00Z")) / 86400.0, 2),
            "description_added_by": version_edit(hist, ev["version"] + 1),
            "position_refined_by": version_edit(hist, ev["version"] + 2),
            "note": ("the collapse date of this record is read off the map's own "
                     "edit history: the version named here added "
                     "damage:date=2021-09-29 within a day of the event, from a "
                     "public changeset"),
        },
        "chip_grid": grid,
        "chip_grid_axes_note": AXIS_NOTE,
        "window": {"px": WINDOW_PX, "note": WINDOW_NOTE},
        "structure": {"tower_base_width_m": TOWER_BASE_WIDTH_M,
                      "rotor_diameter_m": ROTOR_DIAMETER_M,
                      "source": EXTENT_SOURCE},
        "structure_vs_pixel_grid": structure,
        "structure_vs_window": {
            "window_side_px": WINDOW_SIDE_PX,
            "window_px": WINDOW_PX,
            "window_half_extent_px": (WINDOW_SIDE_PX - 1) / 2.0,
            "window_reach_px": WINDOW_REACH_PX,
            "roi_half_px": ROI_HALF_PX,
            "per_orbit": vs_window,
            "note": containment_note(vs_window),
        },
        "sources": {
            "osm_node": OSM_API.format(nid=NODE_ID),
            "osm_node_history": OSM_HISTORY_API.format(nid=NODE_ID),
            "osm_extract": "data/haltern/reference/haltern_osm_node.json",
            "bursts": "data/haltern/reference/haltern_bursts.json",
            "license": LICENSE,
        },
        "status": ("anchor reference - no edge layer is claimed: a 7x7 window per "
                   "acquisition, no rect/pixel-spacing spec and no per-section "
                   "coordinates (see code/haltern/README.md)"),
        "note": ("this anchor supports a statement about *which structure stood "
                 "here and when it collapsed*, not a band/edge layer: such a layer "
                 "needs a resolvable ground extent per pixel along the band axis, "
                 "and ``structure_vs_pixel_grid`` is the number that says how far "
                 "this record is from that"),
    }


def write():
    """Derive the anchor from the committed extract and write it."""
    g = derive()
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(g, fh, indent=1, sort_keys=True)
        fh.write("\n")
    t, ev, w = g["turbine"], g["event_corroboration"], g["structure_vs_window"]
    print(f"written: {os.path.relpath(OUT)}\n"
          f"  node {t['osm_node']} v{t['version']}  {t['manufacturer']} "
          f"{t['model']}  {t['height_m']:.0f} m  "
          f"{t['generator_output_mw']} MW  {t['n_versions']} versions\n"
          f"  asset point is {t['distance_from_asset_point_m']} m from the node "
          f"(bearing {t['bearing_from_asset_point_deg']} deg)\n"
          f"  event: damage:date {ev['damage_date']} recorded by changeset "
          f"{ev['recorded_by']['changeset']} at {ev['recorded_by']['timestamp']} "
          f"(+{ev['lag_days']} d)\n"
          f"  grid: " + "; ".join(
              f"{o} {s['sub_swath']} {s['range_m_per_px']} m/px range, "
              f"{s['azimuth_m_per_px']} m/px azimuth (range axis at "
              f"{s['range_axis_bearing_deg']} deg)"
              for o, s in sorted(g["chip_grid"].items())) + "\n"
          f"  structure vs grid: " + "; ".join(
              f"{o} base {v['base_width_range_px']} x {v['base_width_azimuth_px']} "
              f"px, rotor {v['rotor_diameter_range_px']} x "
              f"{v['rotor_diameter_azimuth_px']} px"
              for o, v in sorted(g["structure_vs_pixel_grid"].items())) + "\n"
          f"  structure vs window: " + "; ".join(
              f"{o} base contained: {str(v['base_contained']).lower()} "
              f"({v['base_outside_by_px']} px outside), "
              f"{v['side_px_to_contain']} px/side would contain it"
              for o, v in sorted(w["per_orbit"].items())))
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
    hist = node_history(raw)
    cur = node_current(raw)
    ev = event_version(hist)
    tags = cur.get("tags", {})
    facts = {
        "node_id": int(cur["id"]),
        "n_versions": len(hist),
        "n_versions_before_event": ev["version"] - 1,
        "event_changeset": ev.get("changeset"),
        "event_edit_ts": ev.get("timestamp"),
        "damage_date": tags.get("damage:date"),
        "destroyed": tags.get("destroyed"),
        "manufacturer": tags.get("manufacturer"),
        "model": tags.get("model"),
        "height_m": float(tags["height"]) if "height" in tags else None,
        "output_mw": 4.5,
        "operator": tags.get("operator"),
    }
    problems = [f"osm.{k}: {facts.get(k)!r} != {v!r}"
                for k, v in sorted(OSM_EXPECTED.items()) if facts.get(k) != v]
    problems += _diff(derive(raw), committed)
    return problems, facts


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch", action="store_true",
                    help="(re)fetch the node and its history from the public OSM API")
    ap.add_argument("--write", action="store_true",
                    help="write data/haltern/reference/haltern_geometry.json")
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
    print(f"ok: {os.path.relpath(RAW)} is node {facts['node_id']} "
          f"({facts['manufacturer']} {facts['model']}, {facts['height_m']:.0f} m), "
          f"{facts['n_versions']} versions, damage:date {facts['damage_date']} "
          f"from changeset {facts['event_changeset']} and "
          f"{os.path.relpath(OUT)} re-derives field by field")


if __name__ == "__main__":
    main()


