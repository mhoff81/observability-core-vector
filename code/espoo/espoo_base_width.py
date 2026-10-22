#!/usr/bin/env python3
"""Measure the Espoo mast base width in a Street View pano, and adopt it.

`data/espoo/asset_evidence.json` records `base_width_m` as
``unconfirmed_default`` (12.0 m, the registry default for a lattice tower) with
the reason "nobody read a width". This script is that reading, in a form that
cannot silently pass:

* the width is computed from **pixels in a named pano at a known distance and a
  declared horizontal field of view** -- every input is recorded;
* **two independent readings** must agree within ``TOLERANCE_M``, otherwise
  nothing is adopted;
* ``--adopt`` rewrites `asset_evidence.json` only after that check passes, and
  only with the readings and the formula beside the value.

Why the tolerance is 1.5 m: the base width only sizes the **support rectangle**
inside the chip (the chip itself is set by the 200 m reference ring), and the
30 m guard means a 1.5 m error cannot clip the mast. See the measurement task in
``code/espoo/README.md``.

Usage
-----
  # what is needed
  python3 code/espoo/espoo_base_width.py --panos

  # one reading: pano id, measured pixels, horizontal FOV of the crop
  python3 code/espoo/espoo_base_width.py --pano AEYr7sSKMbOVLnwPVDRY3g --px 116 --fov 60

  # two readings, adopt if they agree
  python3 code/espoo/espoo_base_width.py \
      --pano AEYr7sSKMbOVLnwPVDRY3g --px 116 --fov 60 \
      --pano SbCRupSUGyADKLwb7veNGQ --px 148 --fov 60 --adopt

  python3 code/espoo/espoo_base_width.py --selftest
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
EVIDENCE = os.path.join(REPO, "data", "espoo", "asset_evidence.json")

#: The panos of ``espoo_site_geometry.md``, by priority (distance, bearing, date).
PANOS = {
    "AEYr7sSKMbOVLnwPVDRY3g": {"distance_m": 63.1, "bearing_deg": 75.5,
                               "date": "2022-10", "note": "best view"},
    "SbCRupSUGyADKLwb7veNGQ": {"distance_m": 49.0, "bearing_deg": 35.4,
                               "date": "2022-10"},
    "GveCM2wS-gN-zWrNFmII_w": {"distance_m": 49.9, "bearing_deg": 46.5,
                               "date": "2022-10"},
    "lSdEAOCz4_FlRTNG352ziQ": {"distance_m": 53.0, "bearing_deg": 56.9,
                               "date": "2022-10"},
    "W9HFZy8cxmXKeMXtRohb3A": {"distance_m": 54.1, "bearing_deg": 13.5,
                               "date": "2022-10"},
    "vvhBpOg7YlLxmWTeFPMLbQ": {"distance_m": 68.2, "bearing_deg": 82.9,
                               "date": "2022-10"},
    "eggBJwZBZH1iEGBuWAXE3w": {"distance_m": 50.2, "bearing_deg": 23.6,
                               "date": "2022-10"},
    "f2MnW9O3tbDOVvzhfNaBLg": {"distance_m": 60.1, "bearing_deg": 4.7,
                               "date": "2022-10"},
    "9c6Glb7byX0fkTMNHBHhqA": {"distance_m": 74.0, "bearing_deg": 88.9,
                               "date": "2022-10"},
}

#: The crop the readings are defined on: 640 px wide, 60 deg horizontal.
CROP_PX = 640
DEFAULT_FOV = 60.0
TOLERANCE_M = 1.5
REGISTRY_DEFAULT_M = 12.0


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------
def width_from_pixels(distance_m, px, fov_h_deg, crop_px=CROP_PX):
    """Width of a target that spans ``px`` in a pano at ``distance_m``.

    ``width = 2 * d * tan(half_angle)`` with ``half_angle = (px / 2) / px_per_deg``
    and ``px_per_deg = crop_px / fov_h_deg`` -- the same construction
    ``asset_evidence.json`` uses for its "12 m = 110 px at 64.5 m" note.
    """
    if distance_m <= 0 or px <= 0 or fov_h_deg <= 0 or crop_px <= 0:
        return None
    px_per_deg = crop_px / float(fov_h_deg)
    half_angle = math.radians((px / 2.0) / px_per_deg)
    return 2.0 * distance_m * math.tan(half_angle)


def pixels_for_width(distance_m, width_m, fov_h_deg, crop_px=CROP_PX):
    """Inverse of :func:`width_from_pixels` -- used by the self-test."""
    half_angle = math.atan((width_m / 2.0) / distance_m)
    px_per_deg = crop_px / float(fov_h_deg)
    return 2.0 * math.degrees(half_angle) * px_per_deg


def reading(pano, px, fov_h_deg):
    """One recorded reading: pano id, pixels, FOV, distance, resulting width."""
    meta = PANOS.get(pano)
    if meta is None:
        raise SystemExit(f"unknown pano {pano!r}; see --panos")
    width = width_from_pixels(meta["distance_m"], px, fov_h_deg)
    return {"pano": pano, "date": meta["date"],
            "distance_m": meta["distance_m"], "bearing_deg": meta["bearing_deg"],
            "px_measured": float(px), "fov_h_deg": float(fov_h_deg),
            "crop_px": CROP_PX, "width_m": None if width is None else round(width, 3)}


def consistency(readings, tolerance=TOLERANCE_M):
    """Do the readings agree? Returns ``(ok, adopted, spread, why)``."""
    values = [r["width_m"] for r in readings if r["width_m"] is not None]
    if len(values) < 2:
        return (False, None, None,
                "at least two independent readings are required")
    spread = max(values) - min(values)
    mean = sum(values) / len(values)
    if spread > tolerance:
        return (False, None, spread,
                f"the readings differ by {spread:.2f} m, more than the "
                f"{tolerance} m tolerance")
    return (True, round(mean, 2), spread, "ok")


def adopt(readings, adopted, spread, path=None, quiet=False):
    """Write the confirmed width into ``asset_evidence.json`` (keeps a backup)."""
    path = path or EVIDENCE
    with open(path) as handle:
        evidence = json.load(handle)
    block = evidence["dimensions"]["base_width_m"]
    if block.get("status") == "confirmed_by_reading":
        raise SystemExit("base_width_m is already confirmed; refusing to overwrite")
    backup = path + ".bak-width"
    with open(backup, "w") as handle:
        json.dump(evidence, handle, indent=1)
        handle.write(chr(10))
    block.update({
        "value_m": adopted,
        "uncertainty_m": round(spread / 2.0, 2),
        "status": "confirmed_by_reading",
        "method": ("pano pixel measurement: width = 2 d tan(px / (2 px_per_deg)), "
                   f"{CROP_PX} px = 1 FOV"),
        "source_type": "Street View pano, human pixel reading",
        "reason": f"two readings within {TOLERANCE_M} m of each other",
        "readings": readings,
        "backup": os.path.basename(backup),
    })
    with open(path, "w") as handle:
        json.dump(evidence, handle, indent=1)
        handle.write(chr(10))
    if not quiet:
        print(f"  adopted base_width_m = {adopted} m (spread {spread:.2f} m)")
        print(f"  wrote {os.path.relpath(path, REPO)} (backup {os.path.basename(backup)})")
    return path


# ---------------------------------------------------------------------------
# Self-test and CLI
# ---------------------------------------------------------------------------
def selftest():
    failures = []

    def check(label, condition, detail=""):
        print(f"  [{'ok  ' if condition else 'FAIL'}] {label}"
              + (f" — {detail}" if detail else ""))
        if not condition:
            failures.append(label)

    print("\n1. The geometry round-trips")
    for width in (11.0, 12.0, 13.0, 15.0):
        px = pixels_for_width(63.1, width, 60.0)
        back = width_from_pixels(63.1, px, 60.0)
        check(f"{width} m at 63.1 m round-trips", abs(back - width) < 1e-9,
              f"{px:.1f} px -> {back:.3f} m")
    px = pixels_for_width(63.1, 12.0, 60.0)
    check("12 m at 63.1 m is ~116 px (as the evidence note says)",
          abs(px - 116) < 1.0, f"{px:.1f} px")

    print("\n2. The tolerance rule")
    good = [reading("AEYr7sSKMbOVLnwPVDRY3g", 116, 60),
            reading("SbCRupSUGyADKLwb7veNGQ", 148, 60)]
    ok, adopted, spread, why = consistency(good)
    check("two agreeing readings are adopted", ok and adopted is not None,
          f"adopted {adopted} m, spread {spread:.2f} m, {why}")
    bad = [reading("AEYr7sSKMbOVLnwPVDRY3g", 116, 60),
           reading("SbCRupSUGyADKLwb7veNGQ", 180, 60)]
    ok, adopted, spread, why = consistency(bad)
    check("two disagreeing readings are refused", not ok and adopted is None,
          f"spread {spread:.2f} m, {why}")
    ok, _adopted, _spread, why = consistency(good[:1])
    check("one reading is not enough", not ok, why)

    print("\n3. The adoption guard")
    block = json.load(open(EVIDENCE))["dimensions"]["base_width_m"]
    check("the evidence still says unconfirmed_default",
          block.get("status") == "unconfirmed_default", block.get("status"))
    check("value_m is still null", block.get("value_m") is None, str(block.get("value_m")))

    print()
    if failures:
        print(f"FAILED — {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


def print_panos():
    print("panos of espoo_site_geometry.md, by priority:")
    print(f"{'pano':26s} {'distance':>9s} {'bearing':>8s} {'12 m':>7s} {'11-13 m':>11s}")
    for pano, meta in PANOS.items():
        px = pixels_for_width(meta["distance_m"], 12.0, DEFAULT_FOV)
        lo = pixels_for_width(meta["distance_m"], 11.0, DEFAULT_FOV)
        hi = pixels_for_width(meta["distance_m"], 13.0, DEFAULT_FOV)
        print(f"{pano:26s} {meta['distance_m']:8.1f} m {meta['bearing_deg']:7.1f} "
              f"{px:6.0f} px {lo:5.0f} … {hi:4.0f} px")
    print()
    print(f"measure: outer edge to outer edge of the two outer legs at ground "
          f"level, at {DEFAULT_FOV:.0f} deg horizontal over {CROP_PX} px.")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Measure and adopt the Espoo mast base width")
    parser.add_argument("--pano", action="append", default=[],
                        help="pano id (repeat once per reading)")
    parser.add_argument("--px", action="append", type=float, default=[],
                        help="measured pixels (one per --pano)")
    parser.add_argument("--fov", action="append", type=float, default=[],
                        help=f"horizontal FOV of the crop (default {DEFAULT_FOV})")
    parser.add_argument("--panos", action="store_true", help="list the panos")
    parser.add_argument("--adopt", action="store_true",
                        help="write the value into asset_evidence.json if the "
                             "readings agree")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args(argv)
    if args.selftest:
        return selftest()
    if args.panos or not args.pano:
        print_panos()
        return 0
    readings = []
    for index, pano in enumerate(args.pano):
        if index >= len(args.px):
            raise SystemExit(f"missing --px for pano {pano}")
        fov = args.fov[index] if index < len(args.fov) else DEFAULT_FOV
        readings.append(reading(pano, args.px[index], fov))
    for entry in readings:
        print(f"  {entry['pano']:26s} {entry['distance_m']:5.1f} m  "
              f"{entry['px_measured']:5.0f} px  {entry['fov_h_deg']:4.0f} deg  "
              f"-> {entry['width_m']:.2f} m")
    ok, adopted, spread, why = consistency(readings)
    print(f"  consistency: {why}" + (f" (spread {spread:.2f} m)" if spread is not None else ""))
    if not ok:
        print("  nothing adopted")
        return 1
    print(f"  adopted value would be {adopted} m")
    if args.adopt:
        adopt(readings, adopted, spread)
    return 0


if __name__ == "__main__":
    sys.exit(main())
