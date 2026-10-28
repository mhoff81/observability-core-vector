#!/usr/bin/env python3
"""The units of the OCV geometry channels, checked on every support chip.

The channel ``A`` is a **count of cells on a raster**, and the raster has two
sizes that are easy to confuse:

    support_cells   the support band -- a geometry constant of the chip
                    (CTS 225, Carola 585, Haltern 545, ...)
    chip cells      rows x cols of the whole chip (CTS 39 x 145 = 5 655)

``A`` counts the mask over the **whole chip**, so it may -- and routinely does --
exceed ``support_cells``: CTS 74 rows, Carola 27, Espoo 45, Haltern 19. Reading
``A`` as "cells of the support band" makes the record look inconsistent (225 vs
247) when nothing is inconsistent; the two numbers are different quantities of
the same raster. This module is the formal check of that reading, and it derives
the unit from the data instead of assuming it:

  (i)   ``support + guard + ring == rows * cols`` where the meta quantifies the
        chip (the three bands partition it);
  (ii)  ``A`` is an integer, ``0 <= A <= rows * cols`` (zero is an empty mask,
        which Bautzen's control series really has -- ``has_mask`` false);
  (iii) ``A / D == bbox_area`` -- exact up to 1e-9 relative, which fixes ``D``
        as ``A`` divided by the mask's bounding box and fixes ``A`` as the
        numerator of that quotient (skipped where the table exports no bbox);
  (iv)  ``A`` against the two columns it is mistaken for: the export's own
        masked-pixel count and ``A_committed``. Agreement with the export can be
        legitimate -- LUMO re-ran the export's rule on the same strip and gets
        the same count on 178 of 178 rows -- but agreement with *both* the
        export and the committed column on every row means the column was
        copied, and that is a defect;
  (v)   the mask is measured on the chip the record declares -- the table's
        ``mask_rows``/``mask_cols`` are the chip's, for every series. (Where no
        chip is declared the shapes are reported, not policed: LUMO cuts one
        rect per width and Bautzen's ASC and DESC chips differ.)

What it does not claim: ``A`` is not comparable between sites (each site has its
own chip) and not comparable with a committed table computed on an 80 x 80
export (``analysis_freeze.md``). The rule it enforces for reports is the one this
repository had to learn: quote ``A`` as "cells of the chip (<= N)", never beside
``support_cells`` without saying which raster each number counts.

Usage:
  python3 code/audit_ocv_units.py                      # audit, print, no write
  python3 code/audit_ocv_units.py --write-reference    # (re)write data/ocv_channel_units.json
  python3 code/audit_ocv_units.py --verify             # audit vs. committed reference
  python3 code/audit_ocv_units.py --selftest           # the checker on synthetic rows
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REFERENCE = os.path.join(DATA, "ocv_channel_units.json")
TOL = 1e-9


def read_table(path):
    """Rows of a channel table, comma- or pipe-delimited (Werder writes pipes)."""
    with open(path, encoding="utf-8") as fh:
        head = fh.readline()
        fh.seek(0)
        if "|" in head and head.count("|") > head.count(","):
            return list(csv.DictReader(fh, delimiter="|"))
        return list(csv.DictReader(fh))


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def audit_rows(rows, chip):
    """The row rules of one table: returns the observed numbers and violations."""
    chip_cells = (chip["rows"] * chip["cols"]
                  if chip["rows"] and chip["cols"] else None)
    obs = {"rows": len(rows), "rows_with_A": 0, "rows_with_export_column": 0,
           "rows_with_committed_column": 0, "rows_with_A_zero": 0,
           "min_A": None, "max_A": None, "min_D": None, "max_D": None,
           "A_gt_support_cells": 0, "A_eq_export_masked_pixels": 0,
           "A_eq_committed_column": 0, "mask_dims": None, "violations": {}}

    def bad(rule, where, detail):
        example = "%s: %s" % (where, detail)
        obs["violations"].setdefault(rule, [])
        if len(obs["violations"][rule]) < 3:
            obs["violations"][rule].append(example)

    dims = {}
    group_column = next((c for c in ("series", "record_class", "track", "role")
                         if rows and c in rows[0]), None)
    for r in rows:
        where = (r.get("id") or r.get("acquisition_date")
                 or r.get("record_key") or "?")
        a, d = _num(r.get("A")), _num(r.get("D"))
        if _num(r.get("mask_rows")) and _num(r.get("mask_cols")):
            group = r.get(group_column) if group_column else "site"
            dims.setdefault(group, set()).add(
                (int(r["mask_rows"]), int(r["mask_cols"])))
        if a is None:
            continue          # a row without a chip (the migrated LUMO table has them)
        obs["rows_with_A"] += 1
        obs["min_A"] = a if obs["min_A"] is None else min(obs["min_A"], a)
        obs["max_A"] = a if obs["max_A"] is None else max(obs["max_A"], a)
        if a == 0:
            obs["rows_with_A_zero"] += 1

        # (ii) A is a count of cells of this chip: an integer in [0, rows x cols].
        # Zero is a legitimate answer -- an empty mask (Bautzen's control series has
        # such dates, ``has_mask`` false) -- and is not a defect of the unit.
        if a != int(a) or a < 0:
            bad("ii_A_is_a_cell_count", where, "A = %r, not an integer >= 0" % a)
        if chip_cells is not None and a > chip_cells:
            bad("ii_A_is_a_cell_count", where,
                "A = %g > %d, the chip's %d cells"
                % (a, chip_cells, chip_cells))
        if chip["support"] and a > chip["support"]:
            obs["A_gt_support_cells"] += 1

        # (iii) D = A / bbox_area wherever the table carries the bounding box
        # (Werder exports D without it, so the rule is recorded as not applicable).
        bbox = _num(r.get("bbox_area"))
        if a == 0:
            if bbox not in (None, 0.0) or d not in (None, 0.0):
                bad("iii_A_over_D_is_bbox_area", where,
                    "A = 0 but bbox_area = %r and D = %r" % (bbox, d))
        elif bbox is None:
            if "bbox_area" in r:
                bad("iii_A_over_D_is_bbox_area", where,
                    "empty bbox_area with A = %g" % a)
        elif bbox != int(bbox) or bbox <= 0:
            bad("iii_A_over_D_is_bbox_area", where, "bbox_area = %r" % bbox)
        elif d is None or d <= 0 or d > 1:
            bad("iii_D_is_a_fraction", where,
                "D = %r with A = %g and bbox_area = %g" % (d, a, bbox))
        else:
            obs["min_D"] = d if obs["min_D"] is None else min(obs["min_D"], d)
            obs["max_D"] = d if obs["max_D"] is None else max(obs["max_D"], d)
            if abs(a / d - bbox) > TOL * max(1.0, abs(bbox)):
                bad("iii_A_over_D_is_bbox_area", where,
                    "A/D = %.12g but bbox_area = %g" % (a / d, bbox))

        # (iv) A against the two columns it could be confused with: the export's
        # own masked-pixel count and the committed value of the same dimension.
        export = _num(r.get("coherence_masked_pixels"))
        if export is not None:
            obs["rows_with_export_column"] += 1
            if export == a:
                obs["A_eq_export_masked_pixels"] += 1
        committed = _num(r.get("A_committed"))
        if committed is not None:
            obs["rows_with_committed_column"] += 1
            if committed == a:
                obs["A_eq_committed_column"] += 1

    shapes = {group: [list(v) for v in values]
              for group, values in ((g, sorted(s)) for g, s in dims.items())}
    obs["mask_dims"] = shapes
    if chip["rows"] and chip["cols"]:
        declared = (chip["rows"], chip["cols"])
        for group, values in sorted(shapes.items()):
            for value in values:
                if tuple(value) != declared:
                    bad("v_mask_is_measured_on_the_chip", str(group),
                        "mask geometry %d x %d, the record declares %d x %d"
                        % (value[0], value[1], declared[0], declared[1]))
    # A copy of the committed column is the defect (iv) is there to catch: the
    # support migration must re-measure, so agreement with *both* the export and
    # the committed column on every row cannot be a recomputation.
    if (obs["rows_with_committed_column"] > 1
            and obs["A_eq_committed_column"] == obs["rows_with_committed_column"]):
        bad("iv_A_is_not_the_committed_column", "table",
            "A equals A_committed on all %d rows -- the column was copied, "
            "not re-measured" % obs["rows_with_committed_column"])
    return obs



def audit_site(site, table_path):
    """Audit one support table against its meta and return its reference block."""
    meta_path = table_path.replace("_ocv_channels.csv", "_ocv_channels_meta.json")
    meta = {}
    if os.path.exists(meta_path):
        with open(meta_path, encoding="utf-8") as fh:
            meta = json.load(fh)
    rows = read_table(table_path)
    chip = chip_of(meta, rows)
    obs = audit_rows(rows, chip)

    # (i) support + guard + ring partition the chip
    cells = (chip["rows"] * chip["cols"]
             if chip["rows"] and chip["cols"] else None)
    message = bands_partition(chip)
    if message:
        obs["violations"].setdefault("i_bands_partition_the_chip", []).append(message)

    block = {
        "table": os.path.relpath(table_path, ROOT),
        "rows": obs["rows"],
        "rows_with_A": obs["rows_with_A"],
        "rows_with_A_zero": obs["rows_with_A_zero"],
        "chip": [chip["rows"], chip["cols"]] if cells else None,
        "chip_cells": cells,
        "chip_source": chip["source"],
        "support_cells": chip["support"],
        "guard_cells": chip["guard"],
        "ring_cells": chip["ring"],
        "mask_dims": obs["mask_dims"],
        "bbox_area_column": bool(rows and "bbox_area" in rows[0]
                                 and any(r.get("bbox_area") for r in rows)),
        "A_range": [obs["min_A"], obs["max_A"]],
        "D_range": [_round(obs["min_D"]), _round(obs["max_D"])],
        "A_gt_support_cells": obs["A_gt_support_cells"],
        "rows_with_export_masked_pixels": obs["rows_with_export_column"],
        "A_eq_export_masked_pixels": obs["A_eq_export_masked_pixels"],
        "rows_with_A_committed": obs["rows_with_committed_column"],
        "A_eq_A_committed": obs["A_eq_committed_column"],
        "checks": {},
    }
    block["checks"] = {rule: len(items) for rule, items in obs["violations"].items()}
    block["examples"] = {rule: items for rule, items in obs["violations"].items()}
    return block


def chip_of(meta, rows):
    """``{rows, cols, support, guard, ring, source}`` -- where the chip comes from.

    Three declarations exist in the record and all three are used, in this order:
    the ``chip`` block of the support metas, the ``mask_rule.window`` of the sites
    that carry no band table (Werder: ``80 x 80``), and the ``mask_rows``/
    ``mask_cols`` columns of the table itself when they are constant.
    """
    chip = meta.get("chip") or {}
    if chip.get("rows") and chip.get("cols"):
        return {"rows": int(chip["rows"]), "cols": int(chip["cols"]),
                "support": chip.get("support_cells"),
                "guard": chip.get("guard_cells"),
                "ring": chip.get("ring_cells"),
                "ring_complete": chip.get("ring_complete"),
                "source": "meta.chip"}
    window = (meta.get("mask_rule") or {}).get("window") or ""
    match = re.search(r"(\d+)\s*x\s*(\d+)", window)
    if match:
        return {"rows": int(match.group(1)), "cols": int(match.group(2)),
                "support": None, "guard": None, "ring": None,
                "ring_complete": None, "source": "meta.mask_rule.window"}
    dims = {(int(r["mask_rows"]), int(r["mask_cols"])) for r in rows
            if _num(r.get("mask_rows")) and _num(r.get("mask_cols"))}
    if len(dims) == 1:
        r, c = dims.pop()
        return {"rows": r, "cols": c, "support": None, "guard": None,
                "ring": None, "ring_complete": None, "source": "table.mask_rows"}
    return {"rows": None, "cols": None, "support": None, "guard": None,
            "ring": None, "ring_complete": None, "source": "undeclared"}


def bands_partition(chip):
    """Rule (i): the band message when support + guard + ring != rows * cols."""
    if None in (chip["support"], chip["guard"], chip["ring"]):
        return None
    if not (chip["rows"] and chip["cols"]):
        return None
    total = int(chip["support"]) + int(chip["guard"]) + int(chip["ring"])
    cells = chip["rows"] * chip["cols"]
    if total == cells:
        return None
    return ("table: %d support + %d guard + %d ring = %d cells, but %d x %d = %d"
            % (chip["support"], chip["guard"], chip["ring"], total,
               chip["rows"], chip["cols"], cells))


def _round(v):
    return None if v is None else round(v, 12)


CONTRACT = {
    "A": "count of cells the site's mask marks on its chip -- a pixel count of a "
         "raster, not a band size and not comparable between sites",
    "D": "A / bbox_area of the mask -- a fraction in (0, 1]",
    "F": "number of 8-connected components of the mask",
    "S": "distance from the mask centroid to the site's peak reference (px)",
    "P": "mean pixel frequency inside the own mask, one map per chip shape",
    "bbox_area": "area of the mask's bounding box, in cells (the denominator of D)",
    "support_cells": "size of the support band, in cells -- a geometry constant of "
                     "the chip; A may exceed it and usually does",
    "guard_cells": "size of the guard band, in cells",
    "ring_cells": "size of the reference ring, in cells",
    "coherence_masked_pixels": "the export's own masked-pixel column -- another "
                               "measurement, never A",
}

RULE = (
    "A is a count of cells of the chip; support_cells is the size of the support "
    "band. The two are different quantities of the same raster, so A > "
    "support_cells is normal (CTS 74 of 178 rows, Carola 27 of 436, Espoo 45 of "
    "75, Haltern 19 of 107) and never an inconsistency. The identity that fixes "
    "the unit is A / D == bbox_area, with 0 <= A <= rows x cols (0 = empty mask) "
    "and D in (0, 1]. Report A as 'cells of the chip (<= N)', never beside "
    "support_cells without saying which raster each number counts, and never "
    "across sites -- each site has its own chip (analysis_freeze.md).")


def build():
    """The whole reference: the contract, every support table, the totals."""
    sites = {}
    pattern = os.path.join(DATA, "*", "support", "*_ocv_channels.csv")
    for table in sorted(glob.glob(pattern)):
        site = os.path.basename(os.path.dirname(os.path.dirname(table)))
        sites[site] = audit_site(site, table)
    totals = {"sites": len(sites),
              "rows": sum(b["rows"] for b in sites.values()),
              "rows_with_A": sum(b["rows_with_A"] for b in sites.values()),
              "violations": sum(sum(b["checks"].values()) for b in sites.values()),
              "sites_with_A_over_support_cells":
                  sum(1 for b in sites.values() if b["A_gt_support_cells"]),
              "A_over_support_cells":
                  sum(b["A_gt_support_cells"] for b in sites.values())}
    return {"schema": "ocv-channel-units-v1", "tolerance": TOL,
            "rule": RULE, "contract": CONTRACT,
            "sites": sites, "totals": totals}


def report_lines(ref):
    """One line per site: the numbers that fixed the unit, then the violations."""
    lines = []
    for site in sorted(ref["sites"]):
        b = ref["sites"][site]
        a_range = ("%g..%g" % tuple(b["A_range"])
                   if b["A_range"][1] is not None else "-")
        chip = ("%d x %d = %d" % (b["chip"][0], b["chip"][1], b["chip_cells"])
                if b["chip"] else "undeclared")
        lines.append("%-12s %4d rows  A %-9s chip %-14s support %-5s "
                     "A>support %3d  A=export %d/%d  %s"
                     % (site, b["rows_with_A"], a_range, chip,
                        b["support_cells"], b["A_gt_support_cells"],
                        b["A_eq_export_masked_pixels"],
                        b["rows_with_export_masked_pixels"],
                        "ok" if not b["checks"] else "FAILED"))
        for rule, items in sorted(b["examples"].items()):
            for item in items:
                lines.append("    %-32s %s" % (rule, item))
    return lines


def diff(left, right, path=""):
    """Differences between two reference blocks (numbers within TOL)."""
    problems = []
    if isinstance(left, tuple):
        left = list(left)
    if isinstance(right, tuple):
        right = list(right)
    if isinstance(left, dict) and isinstance(right, dict):
        for key in sorted(set(left) | set(right)):
            if key not in left or key not in right:
                problems.append("%s%s: present on one side only" % (path, key))
            else:
                problems.extend(diff(left[key], right[key], "%s%s." % (path, key)))
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            problems.append("%s: %d vs %d entries"
                            % (path.rstrip("."), len(left), len(right)))
        else:
            for i, (a, b) in enumerate(zip(left, right)):
                problems.extend(diff(a, b, "%s%d." % (path, i)))
    elif isinstance(left, bool) or isinstance(right, bool):
        if left != right:
            problems.append("%s%r != %r" % (path, left, right))
    elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if abs(float(left) - float(right)) > TOL * max(1.0, abs(float(right))):
            problems.append("%s%r != %r" % (path, left, right))
    elif left != right:
        problems.append("%s%r != %r" % (path, left, right))
    return problems


def selftest():
    """The rules on synthetic rows: each must fire, and a consistent row pass."""
    checks = []

    def check(name, ok, detail):
        checks.append((bool(ok), name, detail))

    cts_chip = {"rows": 39, "cols": 145, "support": 225, "guard": 420,
                "ring": 5010, "ring_complete": True, "source": "test"}
    good = {"id": "t|2021-06-24|target", "A": "67", "D": repr(67 / 3458),
            "bbox_area": "3458", "mask_rows": "39", "mask_cols": "145",
            "coherence_masked_pixels": "640"}
    obs = audit_rows([good], cts_chip)
    check("a consistent row passes every rule, at A < support_cells",
          not obs["violations"] and obs["A_gt_support_cells"] == 0,
          "A = 67 of 5 655 chip cells, 225 support cells, D = 67/3 458")

    bigger = dict(good, A="247", D=repr(247 / 3458))
    obs = audit_rows([bigger], cts_chip)
    check("A above support_cells is normal, and counted",
          not obs["violations"] and obs["A_gt_support_cells"] == 1,
          "A = 247 > 225 support cells, still inside the chip and consistent: "
          "no violation, one row counted")

    div = dict(good, D="0.5")
    obs = audit_rows([div], cts_chip)
    check("a D that is not A / bbox_area is a violation",
          "iii_A_over_D_is_bbox_area" in obs["violations"],
          "D = 0.5 with A = 67, bbox_area = 3 458")
    over = dict(good, A="5656", D=repr(5656 / 5656), bbox_area="5656")
    obs = audit_rows([over], cts_chip)
    check("an A larger than the chip is a violation",
          "ii_A_is_a_cell_count" in obs["violations"],
          "A = 5 656 > 39 x 145 = 5 655")
    frac = dict(good, bbox_area="27", D="2.5")
    obs = audit_rows([frac], cts_chip)
    check("a D outside (0, 1] is a violation",
          "iii_D_is_a_fraction" in obs["violations"], "D = 2.5")
    empty = dict(good, A="0", D="", bbox_area="")
    obs = audit_rows([empty], cts_chip)
    check("an empty mask (A = 0, no bbox, no D) passes",
          not obs["violations"] and obs["rows_with_A_zero"] == 1,
          "A = 0 with empty bbox_area and D -- Bautzen's control series")
    copied = [dict(good, A="640", D=repr(640 / 6400), bbox_area="6400",
                   A_committed="640")] * 2
    obs = audit_rows(copied, cts_chip)
    check("an A copied from A_committed on every row is a violation",
          "iv_A_is_not_the_committed_column" in obs["violations"],
          "A = A_committed = 640 twice: copied, not re-measured")
    recomputed = [dict(good, A="5", D=repr(5 / 1152), bbox_area="1152",
                       coherence_masked_pixels="5", A_committed="7")] * 2
    obs = audit_rows(recomputed, cts_chip)
    check("A equal to the export column, but not to A_committed, passes",
          not obs["violations"] and obs["A_eq_export_masked_pixels"] == 2
          and obs["A_eq_committed_column"] == 0,
          "LUMO's case: the same rule on the same strip gives the same count "
          "(5 = 5), while the committed value differs (7)")
    obs = audit_rows([good, dict(good, mask_rows="80", mask_cols="80")],
                     cts_chip)
    check("a mask measured off the declared chip is a violation",
          "v_mask_is_measured_on_the_chip" in obs["violations"],
          "80 x 80 in the table, 39 x 145 declared")

    check("the CTS bands partition the chip",
          bands_partition(cts_chip) is None,
          "225 + 420 + 5 010 = 5 655 = 39 x 145")
    broken = dict(cts_chip, ring=5009)
    check("a band sum that misses a cell is a violation",
          bands_partition(broken) is not None, "225 + 420 + 5 009 != 5 655")

    failed = sum(1 for ok, _, _ in checks if not ok)
    for ok, name, detail in checks:
        print("%s  %-60s %s" % ("PASS" if ok else "FAIL", name, detail))
    print("\n%d checks, %d failed" % (len(checks), failed))
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="audit_ocv_units.py",
        description="Check the unit of A (cells of the chip) against every support "
                    "table: the bands partition the chip, A is a count inside it, "
                    "A / D is the bounding-box area, the mask sits on the declared "
                    "chip and A was not copied from the export or committed column.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write-reference", action="store_true",
                       help="write the audited reference to data/ocv_channel_units.json")
    group.add_argument("--verify", action="store_true",
                       help="compare the audit with the committed reference")
    group.add_argument("--selftest", action="store_true",
                       help="check the rules on synthetic rows, read nothing")
    args = parser.parse_args(argv)
    if args.selftest:
        return selftest()

    ref = build()
    for line in report_lines(ref):
        print(line)
    totals = ref["totals"]
    print("\n%d sites, %d rows (%d with a chip), %d violations; A > support_cells "
          "on %d rows in %d sites"
          % (totals["sites"], totals["rows"], totals["rows_with_A"],
             totals["violations"], totals["A_over_support_cells"],
             totals["sites_with_A_over_support_cells"]))
    code = 1 if totals["violations"] else 0

    if args.write_reference:
        os.makedirs(os.path.dirname(REFERENCE), exist_ok=True)
        with open(REFERENCE, "w", encoding="utf-8") as fh:
            json.dump(ref, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print("wrote %s" % os.path.relpath(REFERENCE, ROOT))
    elif args.verify:
        with open(REFERENCE, encoding="utf-8") as fh:
            committed = json.load(fh)
        problems = diff(ref, committed)
        for problem in problems[:20]:
            print("  %s" % problem)
        print("%s  the audit reproduces %s"
              % ("OK  " if not problems else "FAIL",
                 os.path.relpath(REFERENCE, ROOT)))
        code = code or (1 if problems else 0)
    return code


if __name__ == "__main__":
    sys.exit(main())

