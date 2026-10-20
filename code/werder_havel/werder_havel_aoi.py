#!/usr/bin/env python3
"""WERDER — the OSM anchor of the site: the bridge's own outline, and what it can support.

Files: ``data/werder_havel/reference/werder_havel_osm_bridge.json`` and
``..._osm_bridge_history.json`` (the committed map extracts) and
``data/werder_havel/reference/werder_havel_aoi.json`` (the derived anchor, this module).

The question this file answers
------------------------------
Haltern anchors to the turbine's own OSM **node** — a point, the strongest anchor of
this family — and YWF to the carriageway **centreline** of the road the tower fell onto
— a line. Werder can do neither, because its record carries **no coordinate at all**:
the export's 28 columns hold the asset identity, the acquisition and the 80 x 80 complex
window, but no latitude, no longitude, no segment index and no burst id. (Haltern's and
YWF's ``segment_latitude``/``segment_longitude`` columns at least existed, empty.) So
the anchor of this record cannot be read out of the record. It is taken from the
**structure itself**, and the structure the record names is its own ``asset_name``:
*Havelbrücke Werder (A 10)*.

The anchor is OSM way ``413381418`` — ``name=Havelbrücke``, ``man_made=bridge``,
``bridge=yes``, ``bridge:structure=beam``, ``material=steel``, ``layer=1``,
``wikidata=Q19963809``. Committed here are the ``/full.json`` response **and** the way's
own ``/history.json``, because both are read below. Three facts of that extract are used:

  * it is a **closed outline**, 10 node refs / 9 distinct nodes, which is the rectangle
    **704.043 m x 35.929 m** (area 25 295.7 m2) — the run of the whole bridge and its
    width over *both* carriageways. That is a mapper's rectangle, not a survey, and it is
    checked against the map's own tag, the only external check this record allows: the way
    carries ``length=704``, and the family's flat metric puts the derived span 0.043 m
    (0.006 %) from it, while the WGS84 span 706.33 m sits 2.33 m (0.33 %) above it — both
    committed, see ``metrics``;
  * its ``description`` names its two superstructures ("Neubau Überbau der
    Richtungsfahrbahn Nauen 1997; Neubau Überbau der Richtungsfahrbahn Werder 2004") and
    ``start_date=1995;2004`` — **construction** years, kept apart from the event below;
  * its five versions (2016-04-25, 2018-12-11, 2021-07-30, 2026-01-30, 2026-02-04) never
    name the 2024/2025 renewal of this record: the two 2026 edits add the structural tags
    above (v4 also adds the vertex that splits the outline's long side) and both are
    *after* the works — and v5 exists only to restore an "N" that v4 lost from its own
    ``description`` ("1997;eubau" -> "1997;Neubau"). The same hand, five days before,
    typed ``length=704`` and the 1997/2004 years: the transcript is a mapper's, and this
    file reads it as such. ``event.from_map`` is therefore **false**.

What this module does *not* claim
---------------------------------
``status``: *anchor reference — no pixel layer is claimed*. The anchor is a **polygon**,
and this family has already measured what that costs: the Bergsøysund analysis found a
-81..+7 m shape error (rms 56 m) when a polygon was used as an anchor. So the derived file
claims **containment and axis** — which orbit/sub-swath footprints cover the structure and
along which bearing it runs — and *not* the position of one pixel on a carriageway. The
cross-axis share is the number that matters: 7 m of shape error is 19 % of the 35.929 m
width, 81 m of it 11 % of the 704.043 m run.

The event of this record is 2024-07-01 .. 2025-09-05 (Autobahn GmbH: expansion joints,
>1800 m of guardrails, mastic asphalt, steel corrosion protection), carried by the
external chain and corroborated only by the record's own strata (``W3_sanierung`` /
``W4_post``, cut on exactly those two dates in ``werder_havel_ocv_masks.py``). The map is
**silent** about it, and the coincidence that the way's own description also speaks of
renewal — of the *other* kind, in 1997 and 2004 — is exactly the reading this file exists
to prevent. ``event`` below commits that in the file, not in prose.

The pixel grid is *not* in this file. The window is 80 x 80 samples (6400) and the
record's ``coherence_masked_pixels`` is 640 = 10 % of it in every one of the 991 rows — a
definition, not a mask. Turning 25 295.7 m2 and 704.043 m into pixels needs the metres per
pixel of the covering bursts, i.e. the committed footprints of
``werder_havel_bursts_resolve.py``; that comparison is ``werder_havel_geometry.json``, the
next step of this package, not this file.

Modes
-----
    --fetch   (re)fetch the way and its history from the public OSM API [network]
    --write   derive the anchor from the committed extracts and write it
    --check   re-derive and compare field by field (the default)
"""


from __future__ import annotations

import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import werder_havel_ocv_masks as masks  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "werder_havel")
REF = os.path.join(DATA, "reference")
RAW = os.path.join(REF, "werder_havel_osm_bridge.json")
RAW_HIST = os.path.join(REF, "werder_havel_osm_bridge_history.json")
OUT = os.path.join(REF, "werder_havel_aoi.json")

OSM_API = "https://api.openstreetmap.org/api/0.6/way/{wid}/full.json"
OSM_HIST_API = "https://api.openstreetmap.org/api/0.6/way/{wid}/history.json"
LICENSE = ("\u00a9 OpenStreetMap contributors, ODbL 1.0 "
           "(https://www.openstreetmap.org/copyright)")

WAY_ID = 413381418

# The record's event, as the family carries it (never a precursor): the documented
# change the strata are cut on. It is NOT in the map — see the module docstring.
EVENT_KIND = "documented change (Autobahn GmbH, A 10 renewal)"

# ---------------------------------------------------------------------------
# The structural identity of the committed extract
# ---------------------------------------------------------------------------
# None of this can be re-derived from the record (which carries no coordinate and no
# name), so it is asserted (the family's ``OSM_EXPECTED`` pattern), and the geometry is
# asserted against the map's own ``length`` tag — the one external check this anchor
# allows. If the way is edited upstream, ``--check`` fails here first and loudly rather
# than silently re-basing the anchor of a committed record.
OSM_EXPECTED = {
    "way_id": WAY_ID,
    "version": 5,
    "timestamp": "2026-02-04T18:11:39Z",
    "changeset": 178106321,
    "user": "Waldamsel",
    "n_node_refs": 10,
    "n_distinct_nodes": 9,
    "closed": True,
    "n_versions": 5,
    "version_timestamps": [
        "2016-04-25T16:34:52Z",
        "2018-12-11T20:03:10Z",
        "2021-07-30T16:59:13Z",
        "2026-01-30T15:56:19Z",
        "2026-02-04T18:11:39Z",
    ],
    "version_node_refs": [9, 9, 9, 10, 10],
    "tagged_length_m": 704.0,
    "tags": {
        "bridge": "yes",
        "bridge:structure": "beam",
        "description": ("Neubau \u00dcberbau der Richtungsfahrbahn Nauen 1997;"
                        "Neubau \u00dcberbau der Richtungsfahrbahn Werder 2004"),
        "layer": "1",
        "length": "704",
        "man_made": "bridge",
        "material": "steel",
        "name": "Havelbr\u00fccke",
        "start_date": "1995;2004",
        "wikidata": "Q19963809",
        "wikimedia_commons": "Category:Havelbr\u00fccke Werder",
        "wikipedia": "de:Havelbr\u00fccke Werder",
    },
}

# ---------------------------------------------------------------------------
# Two metrics, both committed (a finding of this file, not a footnote)
# ---------------------------------------------------------------------------
# ``_local`` is the family's projection (``ywf_osm_anchor.py``): equirectangular about
# the anchor, M_PER_DEG_LAT 110540 / M_PER_DEG_LON 111320. It is used for every geometry
# number below *because the burst footprints of the next step are converted with the same
# constants*, so the metres-per-pixel comparison of the pixel grid stays self-consistent.
# But it is NOT the ellipsoid: at 52.4175 N the WGS84 degree lengths are 111275.2 m
# (meridional) and 68037.2 m (parallel), i.e. the family's constants are 0.66 % and
# 0.21 % low. On a 704 m span that is 2.29 m (0.325 %), on 25 295.7 m2 0.87 % — so the
# WGS84 figures are committed next to the flat ones rather than left implicit, and the
# map's own ``length=704`` tag can be read against both. Both are true; they answer
# different questions (grid consistency vs the ellipsoid).
M_PER_DEG_LAT = 110540.0
M_PER_DEG_LON = 111320.0

METRIC_NOTE = (
    "flat = the family's equirectangular constants (ywf_osm_anchor.py), used for every "
    "'*_m' / '*_m2' field because the burst footprints that the pixel-grid audit divides "
    "by are converted with them; wgs84 = the same quantity against the local WGS84 "
    "degree lengths (111275.2 m/deg meridional, 68037.2 m/deg parallel at the anchor). "
    "The family's constants are 0.66 % / 0.21 % low, i.e. 2.29 m (0.325 %) on this span."
)

# ---------------------------------------------------------------------------
# Geometry (all of it defined here, none of it inherited from a previous session)
# ---------------------------------------------------------------------------


def meridional_m_per_deg(lat):
    """WGS84 meridional degree length (standard series, mm-accurate for |lat| < 80)."""
    p = math.radians(lat)
    return (111132.92 - 559.82 * math.cos(2 * p) + 1.175 * math.cos(4 * p)
            - 0.0023 * math.cos(6 * p))


def parallel_m_per_deg(lat):
    """WGS84 parallel degree length at ``lat`` (prime-vertical radius times cos)."""
    a, f = 6378137.0, 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    p = math.radians(lat)
    n = a / math.sqrt(1.0 - e2 * math.sin(p) ** 2)
    return n * math.cos(p) * math.pi / 180.0


def _local(lat, lon, lat0, lon0):
    """Local planar metres east/north of ``(lat0, lon0)`` (equirectangular).

    The family's projection, verbatim (``ywf_osm_anchor.py``): adequate at the scale of
    this question because the outline is 0.7 km, three orders below the scale where the
    flat approximation matters — and used for every metric field here so that the
    metres-per-pixel comparison of the pixel-grid audit, which converts the burst
    footprints with the same constants, stays self-consistent. Its known 0.32 % bias on
    this span is committed separately (``*_wgs84_*`` fields, ``METRIC_NOTE``).
    """
    return ((lon - lon0) * M_PER_DEG_LON * math.cos(math.radians(lat0)),
            (lat - lat0) * M_PER_DEG_LAT)


def geodesic_m(a, b):
    """The same edge against the local WGS84 degree lengths, at the edge's mid-lat."""
    lat_mid = 0.5 * (a[0] + b[0])
    dn = (b[0] - a[0]) * meridional_m_per_deg(lat_mid)
    de = (b[1] - a[1]) * parallel_m_per_deg(lat_mid)
    return math.hypot(dn, de)


def bearing_deg(a, b):
    """Initial great-circle bearing a -> b in [0, 360), the family's formula."""
    la1, lo1 = math.radians(a[0]), math.radians(a[1])
    la2, lo2 = math.radians(b[0]), math.radians(b[1])
    y = math.sin(lo2 - lo1) * math.cos(la2)
    x = (math.cos(la1) * math.sin(la2)
         - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1))
    return math.degrees(math.atan2(y, x)) % 360.0


def node_centroid(pts):
    """Arithmetic mean of the distinct outline nodes: the record's own site point."""
    n = len(pts)
    return (sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n)


def ring_area_m2(xy):
    """Shoelace area of a closed ring of (east, north) metres; the sign is dropped."""
    s = 0.0
    for i in range(len(xy) - 1):
        s += xy[i][0] * xy[i + 1][1] - xy[i + 1][0] * xy[i][1]
    return abs(s) / 2.0


def edge_lengths_m(pts, lat0, lon0):
    """The flat length of every edge of the closed ring, in node order."""
    xy = [_local(la, lo, lat0, lon0) for la, lo in pts]
    return [math.dist(xy[i], xy[i + 1]) for i in range(len(xy) - 1)]


def unit(a, b):
    """Unit vector a -> b, or (0, 0) for a degenerate edge."""
    d = math.hypot(b[0] - a[0], b[1] - a[1])
    return (0.0, 0.0) if d == 0 else ((b[0] - a[0]) / d, (b[1] - a[1]) / d)


def ring_sides_and_ends(pts, lat0, lon0, i_span, width,
                        parallel_cos=math.cos(math.radians(45.0))):
    """Split the ring's edges into the two side chains and the end edges.

    By **direction, not length rank**: of the nine edges of this outline the longest is a
    side (704.043 m), but the side opposite it is drawn as *two collinear edges*
    (662.458 + 41.585 m) because the mapper split it with an extra vertex — so ranking by
    length would file a 41.585 m piece of the long side as an end, and any such split
    would still leave the two chains ungrouped, because both of them lie on the same side
    of the *line through the span edge*. An edge is therefore a side when it runs within
    45 deg of the span edge, an end otherwise, and the sides are separated by how far
    their midpoint sits from that line: the span chain sits on it (0 m) and the chain
    opposite it a full width away, so ``width / 2`` is the natural cut between them.
    Returns ``(sides, ends, offsets)`` with the perpendicular offsets for every edge.
    """
    xy = [_local(la, lo, lat0, lon0) for la, lo in pts]
    ux, uy = unit(xy[i_span], xy[i_span + 1])
    ox, oy = xy[i_span]
    sides, ends, offsets = [[], []], [], []
    for j in range(len(xy) - 1):
        dx = xy[j + 1][0] - xy[j][0]
        dy = xy[j + 1][1] - xy[j][1]
        d = math.hypot(dx, dy)
        par = abs(dx * ux + dy * uy) / d if d else 0.0
        mx = 0.5 * (xy[j][0] + xy[j + 1][0]) - ox
        my = 0.5 * (xy[j][1] + xy[j + 1][1]) - oy
        off = abs(mx * uy - my * ux)
        offsets.append(off)
        if par < parallel_cos:
            ends.append(j)
        else:
            sides[0 if off <= 0.5 * width else 1].append(j)
    return sides, ends, offsets


def ring_area_centroid(xy):
    """(east, north) of the ring's area centroid, in the same local metres as ``xy``."""
    a2 = cx = cy = 0.0
    for i in range(len(xy) - 1):
        cr = xy[i][0] * xy[i + 1][1] - xy[i + 1][0] * xy[i][1]
        a2 += cr
        cx += (xy[i][0] + xy[i + 1][0]) * cr
        cy += (xy[i][1] + xy[i + 1][1]) * cr
    if a2 == 0.0:
        return (0.0, 0.0)
    return (cx / (3.0 * a2), cy / (3.0 * a2))


def rectangle_from_area_perimeter(area_m2, perimeter_m):
    """(long, short) side of the rectangle with that area and perimeter, or None.

    Independent of the node decomposition — it uses only the shoelace area and the ring
    perimeter — so agreement with the longest-edge reading is a real cross-check, and for
    this outline the fit is exact: the two side chains are each 704.0426 m long
    (662.4581 + 41.5845 on one, 704.0426 on the other) and the six end edges sum to
    71.8585 m, i.e. 2 x the 35.9292 m width to 0.04 mm — the mapper drew a rectangle whose
    sides are subdivided by extra nodes.
    """
    half = perimeter_m / 2.0
    disc = half * half - 4.0 * area_m2
    if disc < 0:
        return None
    root = math.sqrt(disc)
    return ((half + root) / 2.0, (half - root) / 2.0)


# ---------------------------------------------------------------------------
# The committed map extracts
# ---------------------------------------------------------------------------
def fetch():
    """Fetch the way and its version history from the public OSM API and write both."""
    import requests

    written = []
    for path, url in ((RAW, OSM_API.format(wid=WAY_ID)),
                      (RAW_HIST, OSM_HIST_API.format(wid=WAY_ID))):
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        d = r.json()
        if "error" in d:
            raise SystemExit(f"{url}: {d['error']}")
        d["_source"] = url
        d["_license"] = LICENSE
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w") as fh:
            json.dump(d, fh, indent=1, sort_keys=True)
            fh.write("\n")
        written.append(d)
        print(f"written: {os.path.relpath(path)} — {len(d['elements'])} elements")
    raw, hist = written
    way = [e for e in raw["elements"] if e.get("type") == "way"][0]
    print(f"         way {way['id']} v{way['version']}: {len(way['nodes'])} node refs, "
          f"{len(hist['elements'])} versions")
    return raw, hist


def load_raw(path=RAW):
    with open(path) as fh:
        return json.load(fh)


def load_history(path=RAW_HIST):
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


def way_versions(hist):
    """The way's versions, oldest first, as the small dicts the anchor commits.

    The ``tags`` of every version are kept because they are the evidence that the map
    never names the record's event: the 2026 versions add structural tags and repair a
    typo, they do not record the 2024/2025 renewal.
    """
    out = []
    for e in hist["elements"]:
        if e.get("type") != "way":
            continue
        out.append({"version": e["version"],
                    "timestamp": e["timestamp"],
                    "changeset": e["changeset"],
                    "user": e.get("user", ""),
                    "n_node_refs": len(e["nodes"]),
                    "tags": dict(sorted(e.get("tags", {}).items()))})
    out.sort(key=lambda v: v["version"])
    return out


# ---------------------------------------------------------------------------
# The anchor
# ---------------------------------------------------------------------------
def derive(raw=None, hist=None):
    """The committed anchor, derived from the two map extracts and the committed record.

    Pure: every field is a function of ``RAW``, ``RAW_HIST`` and the committed
    measurements file, so ``--check`` rebuilds it and compares it field by field.
    """
    raw = load_raw() if raw is None else raw
    hist = load_history() if hist is None else hist
    c = masks.census()
    way, pts = way_geometry(raw)
    refs = way["nodes"]
    coords = {e["id"]: (e["lat"], e["lon"]) for e in raw["elements"]
              if e.get("type") == "node"}
    uniq = list(dict.fromkeys(refs))
    site = node_centroid([coords[n] for n in uniq])
    lat0, lon0 = site
    xy = [_local(la, lo, lat0, lon0) for la, lo in pts]
    cen = ring_area_centroid(xy)
    cen_ll = (lat0 + cen[1] / M_PER_DEG_LAT,
              lon0 + cen[0] / (M_PER_DEG_LON * math.cos(math.radians(lat0))))
    d_cen = math.dist(_local(site[0], site[1], lat0, lon0), cen)
    edges = [math.dist(xy[i], xy[i + 1]) for i in range(len(xy) - 1)]
    perimeter = sum(edges)
    area = ring_area_m2(xy)
    order = sorted(range(len(edges)), key=lambda j: -edges[j])
    i_span = order[0]
    span = edges[i_span]
    width = area / span
    sides, ends, offsets = ring_sides_and_ends(pts, lat0, lon0, i_span, width)
    side_sums = [sum(edges[j] for j in ch) for ch in sides]
    side_offsets = [[round(min(offsets[j] for j in ch), 3),
                     round(max(offsets[j] for j in ch), 3)] for ch in sides if ch]
    end_sum = sum(edges[j] for j in ends)
    fit = rectangle_from_area_perimeter(area, perimeter)
    span_wgs = geodesic_m(pts[i_span], pts[i_span + 1])
    area_wgs = area * ((M_PER_DEG_LAT / meridional_m_per_deg(lat0))
                       * (M_PER_DEG_LON * math.cos(math.radians(lat0))
                          / parallel_m_per_deg(lat0)))
    tagged = OSM_EXPECTED["tagged_length_m"]
    lats = [p[0] for p in pts]
    lons = [p[1] for p in pts]
    return {
        "asset": masks.ASSET,
        "asset_id": masks.ASSET_ID,
        "asset_name": masks.ASSET_NAME,
        "asset_type": masks.ASSET_TYPE,
        "status": ("anchor reference — no pixel layer is claimed; containment and "
                   "axis only (the record carries no coordinate of its own)"),
        "site_point": {"lat": round(site[0], 7), "lon": round(site[1], 7)},
        "site_point_kind": "arithmetic mean of the outline's 9 distinct nodes",
        "site_point_area_centroid": {"lat": round(cen_ll[0], 7),
                                     "lon": round(cen_ll[1], 7)},
        "site_point_vs_area_centroid_m": round(d_cen, 3),
        "site_point_note": (
            "from the map extract, not from the record: the export carries no latitude "
            "and no longitude, so there is no exported point to measure against "
            "(Haltern's ``distance_from_asset_point_m`` has no counterpart here). The "
            "*node* mean is the anchor because the mapped vertices are what the extract "
            "actually commits and the mean is hand-checkable; the ring's area centroid "
            "is committed beside it with the distance between them, because the nodes "
            "are not evenly spaced along the sides and a reader should see that the "
            "choice is worth a few metres rather than take it on trust. Neither point "
            "is a surveyed point of the record."),
        "anchor": {
            "kind": "structure_outline",
            "why": (
                "the record's own asset_name names the bridge and the record carries no "
                "coordinate, so the structure itself is the anchor. The family measured "
                "-81..+7 m (rms 56 m) of shape error for a polygon anchor "
                "(Bergsøysund), which is why the derived file claims containment and "
                "axis only: 7 m of shape error is 19 % of the 35.929 m width, 81 m of it "
                "11 % of the 704.043 m run"),
            "osm_way": way["id"],
            "version": way["version"],
            "timestamp": way["timestamp"],
            "changeset": way["changeset"],
            "user": way.get("user", ""),
            "tags": dict(sorted(way.get("tags", {}).items())),
            "n_node_refs": len(refs),
            "n_distinct_nodes": len(uniq),
            "closed": refs[0] == refs[-1],
            "node_ids": refs,
            "edge_lengths_m": [round(e, 3) for e in edges],
            "side_edge_indexes": sides,
            "side_chain_lengths_m": [round(s, 3) for s in side_sums],
            "side_chain_offsets_m": side_offsets,
            "end_edge_indexes": ends,
            "end_edges_m": [round(edges[j], 3) for j in ends],
            "end_edges_sum_m": round(end_sum, 3),
            "ring_perimeter_m": round(perimeter, 3),
            "area_m2": round(area, 1),
            "area_ha": round(area / 1e4, 3),
            "span_m": round(span, 3),
            "span_edge_index": i_span,
            "span_edge_nodes": [refs[i_span], refs[i_span + 1]],
            "span_bearing_deg": round(bearing_deg(pts[i_span], pts[i_span + 1]), 2),
            "width_m": round(width, 3),
            "width_is": ("area / span — the mean width of the whole outline, i.e. both "
                         "carriageways and the median between them, not one deck"),
            "rectangle_fit": {
                "span_from_area_perimeter_m": round(fit[0], 3),
                "width_from_area_perimeter_m": round(fit[1], 3),
                "span_minus_area_over_span_m": round(fit[0] - span, 3),
                "width_minus_area_over_span_m": round(fit[1] - width, 3),
                "side_chain_minus_span_m": round(max(abs(s - span) for s in side_sums), 3),
                "end_edges_sum_minus_2width_m": round(end_sum - 2.0 * width, 3),
                "note": (
                    "from the shoelace area and the ring perimeter alone — no node "
                    "decomposition — the rectangle with those two invariants is "
                    "704.04 x 35.93 m, i.e. the longest-edge reading to under 0.01 m; and "
                    "the two side chains, split by direction and separated by their "
                    "offset from the span edge's line, are each the span while the end "
                    "edges sum to twice the width. That is what makes the outline a "
                    "rectangle for the purposes of this anchor, and both readings are "
                    "re-derived here rather than asserted."),
            },
            "length_tag": {
                "tagged_m": tagged,
                "added_in_version": 4,
                "span_flat_minus_tagged_m": round(span - tagged, 3),
                "span_wgs84_minus_tagged_m": round(span_wgs - tagged, 2),
                "note": (
                    "the way's own integer ``length`` tag is the only external check "
                    "this anchor allows, and it was typed by hand in the same 2026-01-30 "
                    "changeset as a ``description`` that had to be corrected five days "
                    "later (v4 '1997;eubau …' -> v5 '1997;Neubau …'): it carries the "
                    "*scale* of the outline, it is not a survey, and the residual "
                    "against each metric is committed rather than smoothed away."),
            },
        },
        "metrics": {
            "flat": {
                "m_per_deg_lat": M_PER_DEG_LAT,
                "m_per_deg_lon": M_PER_DEG_LON,
                "cos_lat0": round(math.cos(math.radians(lat0)), 6),
            },
            "flat_span_m": round(span, 3),
            "flat_area_m2": round(area, 1),
            "wgs84_span_m": round(span_wgs, 2),
            "wgs84_area_m2": round(area_wgs, 1),
            "wgs84_deg_lengths_m": {"meridional": round(meridional_m_per_deg(lat0), 1),
                                    "parallel": round(parallel_m_per_deg(lat0), 1)},
            "wgs84_minus_flat_span_m": round(span_wgs - span, 2),
            "wgs84_minus_flat_span_pct": round(100.0 * (span_wgs - span) / span, 3),
            "note": METRIC_NOTE,
        },


        "aoi": {
            "crs": "EPSG:4326",
            "bbox": {"west": round(min(lons), 7), "south": round(min(lats), 7),
                     "east": round(max(lons), 7), "north": round(max(lats), 7)},
            "bbox_extent_m": {"ns": round((max(lats) - min(lats)) * M_PER_DEG_LAT, 1),
                              "ew": round((max(lons) - min(lons)) * M_PER_DEG_LON
                                          * math.cos(math.radians(lat0)), 1)},
            "polygon_lonlat": [[round(p[1], 7), round(p[0], 7)] for p in pts],
            "polygon_role": (
                "the /Intersects() polygon of the CDSE OData burst query in "
                "werder_havel_bursts_resolve.py — the committed extract's own node ring, "
                "in node order and closed, lon/lat (GeoJSON axis order)"),
            "bbox_note": (
                "the bbox is the ring's own extent, not a buffered box: the query AOI is "
                "the structure itself, and the resolver reports which bursts cover it "
                "rather than guessing a margin"),
        },
        "event": {
            "start": masks.EVENT_START,
            "end": masks.EVENT_END,
            "kind": EVENT_KIND,
            "source": masks.EVENT_SOURCE,
            "from_map": False,
            "map_note": (
                "none of the way's five versions names the 2024/2025 renewal: v4/v5 "
                "(2026) add the structural tags, and the ``description`` names the "
                "*construction* of the two superstructures (1997 / 2004 — the "
                "'Richtungsfahrbahn Nauen' and 'Richtungsfahrbahn Werder' decks). The "
                "event is the external chain's, corroborated by the record's own strata "
                "only (W3_sanierung / W4_post are cut on exactly these two dates)"),
        },
        "versions": way_versions(hist),
        "record": {
            "measurements": c["measurements"],
            "measurements_sha256": c["measurements_sha256"],
            "rows": c["rows"],
            "columns": c["columns"],
            "n_dates": c["n_dates"],
            "first_date": c["first_date"],
            "last_date": c["last_date"],
            "n_tracks": c["n_tracks"],
            "tracks": c["tracks"],
            "roles": c["roles"],
            "state_windows": c["state_windows"],
            "state_window_ranges": c["state_window_ranges"],
            "window_size": c["window_size"],
            "window_px": c["window_px"],
            "masked_pixels": c["masked_pixels"],
            "masked_pixels_frac_of_window": c["masked_pixels_frac_of_window"],
        },
        "notes": [
            "the anchor is a polygon and a mapper's: containment and axis only, per the "
            "docstring — the record itself carries no coordinate to check it against",
            "coherence_masked_pixels = 640 in all 991 rows is 10 % of the 6400-sample "
            "window by definition, not a mask over this area "
            "(werder_havel_ocv_masks.py, WINDOW_NOTE)",
            "the pixel grid is not in this file: an 80 x 80 window carries no pixel "
            "spacing spec, so the metres-per-pixel audit is "
            "data/werder_havel/reference/werder_havel_geometry.json after "
            "werder_havel_bursts_resolve.py commits the covering burst footprints",
            "the extract is mutable external input (ODbL): if way 413381418 is edited "
            "upstream, --check fails on OSM_EXPECTED and on the derived file rather than "
            "re-basing the anchor of a committed record",
        ],
    }


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


def check(raw=None, hist=None):
    """Re-verify the committed anchor offline. Returns ``(problems, facts)``."""
    raw = load_raw() if raw is None else raw
    hist = load_history() if hist is None else hist
    way, pts = way_geometry(raw)
    refs = way["nodes"]
    versions = way_versions(hist)
    problems = []

    got = {
        "way_id": way["id"],
        "version": way["version"],
        "timestamp": way["timestamp"],
        "changeset": way["changeset"],
        "user": way.get("user", ""),
        "n_node_refs": len(refs),
        "n_distinct_nodes": len(set(refs)),
        "closed": refs[0] == refs[-1],
        "n_versions": len(versions),
        "version_timestamps": [v["timestamp"] for v in versions],
        "version_node_refs": [v["n_node_refs"] for v in versions],
        "tagged_length_m": float(way["tags"]["length"]),
        "tags": dict(sorted(way.get("tags", {}).items())),
    }
    for k, want in OSM_EXPECTED.items():
        if got.get(k) != want:
            problems.append(f"osm.{k}: {got.get(k)!r} != {want!r}")

    # ``event.from_map: False`` is only true while no version of the way names the
    # event: if a mapper ever writes the renewal into it, this fails and the claim has
    # to be re-read rather than left standing as a silent contradiction.
    years = sorted({masks.EVENT_START[:4], masks.EVENT_END[:4]})
    for v in versions:
        blob = " ".join(f"{k} {val}" for k, val in v["tags"].items())
        hit = [y for y in years if y in blob]
        if hit:
            problems.append(f"version {v['version']} of the way names {hit}: "
                            f"event.from_map must be re-read, not assumed false")

    facts = derive(raw, hist)
    if not os.path.exists(OUT):
        problems.append(f"{os.path.relpath(OUT)}: missing — run --write")
        return problems, facts
    with open(OUT) as fh:
        committed = json.load(fh)
    problems += _diff(committed, facts, "aoi")

    # the geometry is still the rectangle the anchor claims, in both readings
    a, rf = facts["anchor"], facts["anchor"]["rectangle_fit"]
    if abs(rf["end_edges_sum_minus_2width_m"]) > 0.05:
        problems.append("the end edges no longer sum to twice the width: "
                        f"{rf['end_edges_sum_minus_2width_m']}")
    if rf["side_chain_minus_span_m"] > 0.1:
        problems.append("the two side chains are no longer each the span: "
                        f"{rf['side_chain_minus_span_m']}")
    if abs(rf["span_minus_area_over_span_m"]) > 0.1:
        problems.append("the outline is no longer the rectangle the anchor claims: "
                        f"{rf['span_minus_area_over_span_m']}")
    if abs(committed["metrics"]["flat_span_m"] - a["span_m"]) > 1e-9:
        problems.append("the committed metrics block disagrees with the anchor span")
    return problems, facts


def write():
    """Derive the anchor from the committed extracts and write it."""
    facts = derive(load_raw(), load_history())
    os.makedirs(REF, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(facts, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(f"written: {os.path.relpath(OUT)}")
    return facts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fetch", action="store_true",
                    help="(re)fetch the way and its history from the public OSM API")
    ap.add_argument("--write", action="store_true",
                    help="derive the anchor from the committed extracts and write it")
    ap.add_argument("--check", action="store_true",
                    help="re-derive and compare field by field (the default)")
    ap.add_argument("--json", action="store_true", help="print the derived anchor")
    args = ap.parse_args(argv)

    if args.fetch:
        fetch()
    if args.fetch or args.write:
        facts = write()
        if args.json:
            print(json.dumps(facts, indent=1, sort_keys=True))
        return 0

    problems, facts = check()
    a = facts["anchor"]
    if problems:
        print(f"CHECK FAILED — {len(problems)} problem(s)")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"OK — way {a['osm_way']} v{a['version']}: {a['n_node_refs']} node refs "
          f"({a['n_distinct_nodes']} distinct, closed {a['closed']}); "
          f"span {a['span_m']} m on bearing {a['span_bearing_deg']} deg, "
          f"width {a['width_m']} m, area {a['area_m2']} m2")
    print(f"     site point {facts['site_point']['lat']}, {facts['site_point']['lon']} "
          f"({facts['site_point_vs_area_centroid_m']} m from the ring's area centroid); "
          f"tagged length {a['length_tag']['tagged_m']:.0f} m: flat residual "
          f"{a['length_tag']['span_flat_minus_tagged_m']:+.3f} m, WGS84 residual "
          f"{a['length_tag']['span_wgs84_minus_tagged_m']:+.2f} m")
    print(f"     event {facts['event']['start']} .. {facts['event']['end']}, "
          f"from_map {facts['event']['from_map']}; {len(facts['versions'])} map "
          f"versions; {facts['record']['rows']} record rows over "
          f"{facts['record']['n_dates']} dates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
