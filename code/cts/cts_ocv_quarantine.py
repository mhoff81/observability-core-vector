#!/usr/bin/env python3
"""CTS -- the quarantined asset metadata, checked against the record itself.

``data/ocv_metadata_quarantine.json`` says that ``asset_name``/``asset_type`` of
the CTS export are wrong (they carry the Werder bridge's labels while
``asset_id`` is CTS's) and that identity comes from ``track``/``role``/``day``.
This module makes that statement checkable instead of asserted:

  A. the register's pins   -- sha256 and row count of every file it names
  B. the register's claim  -- the quarantined columns really hold those values,
                              on every row of both CTS tables
  C. the identity rule     -- ``identity_key`` is unique over the whole table,
                              and rewriting the quarantined columns leaves every
                              key unchanged (the poison test)
  D. the cited evidence    -- ``cts_segments.txt`` labels the same ``asset_id``
                              "Champlain Towers South", ``cts_asset.json`` names
                              the site with ``asset_id: null``, and Werder pairs
                              the same name with a *different* id
  E. the code guard        -- no literal read of a quarantined column anywhere
                              in ``code/cts/`` (AST scan)

Scope of E, stated honestly: the scan sees literal column names
(``row["asset_name"]``, ``row.get("asset_name")``). A generic loop over
``SOURCE_TEXT`` is invisible to it -- and that generic copy is exactly how the
column reaches the CSV. E therefore guards against *using* the column, not
against copying it; C is the guard that matters, because it proves that no
identity derived here depends on it.

Usage:
  python3 code/cts/cts_ocv_quarantine.py            # check the record
  python3 code/cts/cts_ocv_quarantine.py --selftest # and the checker itself
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cts_ocv_core as core    # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
REGISTER = core.QUARANTINE_PATH
SUPPORT_CSV = os.path.join(ROOT, "data", "cts", "support", "cts_ocv_channels.csv")
SEGMENTS = os.path.join(ROOT, "data", "cts", "cts_segments.txt")
WERDER_META = os.path.join(ROOT, "data", "werder_havel", "support",
                           "werder_havel_ocv_channels_meta.json")

#: The register's column names, taken from the register / from ``core`` instead
#: of being typed here: this module must not read a quarantined column by
#: literal name, or its own AST scan would -- rightly -- flag it.
ID_COLUMN = "asset_id"                       # correct value, still not identity
NAME_COLUMN = core.QUARANTINED_COLUMNS[0]    # quarantined (wrong value)
TYPE_COLUMN = core.QUARANTINED_COLUMNS[1]    # quarantined (wrong value)
REQUEST_COLUMN = core.MISNAMED_COLUMNS[0]    # right value, wrong name


def sha256(path):
    """Hex digest of a file, read in chunks (the repository's pin)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_pipe(path):
    """A pipe-delimited table with a header row (``#`` comments skipped).

    Rows are padded to the header width. ``cts_segments.txt`` omits a trailing
    empty field on every row (``...|label|``), and the strict reader used by the
    site package therefore returns **no** rows from it -- the table is carried as
    provenance only, which is exactly why nobody noticed. This reader pads, so
    the evidence the register cites can be read at all.
    """
    rows, header = [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("|")
            if header is None:
                header = parts
                continue
            if len(parts) <= len(header):
                parts += [""] * (len(header) - len(parts))
                rows.append(dict(zip(header, parts)))
    return rows


def read_csv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def literal_reads(source):
    """``[(lineno, column)]`` of quarantined names read as literals in ``source``.

    Two forms are counted, because they are the two ways a column is read by
    name: a subscript with a string key (``row["asset_name"]``) and a ``.get``
    call whose first argument is a string (``row.get("asset_name")``).
    """
    wanted = set(core.QUARANTINED_COLUMNS) | set(core.MISNAMED_COLUMNS)
    found = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Subscript):
            key = node.slice
            if isinstance(key, ast.Constant) and key.value in wanted:
                found.add((node.lineno, key.value))
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == "get" and node.args:
            arg = node.args[0]
            if isinstance(arg, ast.Constant) and arg.value in wanted:
                found.add((node.lineno, arg.value))
    return sorted(found)


def scan_package(folder=None):
    """``[(relpath, lineno, column)]`` of literal reads over ``code/cts/*.py``."""
    folder = folder or HERE
    hits = []
    for name in sorted(os.listdir(folder)):
        if not name.endswith(".py"):
            continue
        path = os.path.join(folder, name)
        with open(path, encoding="utf-8") as fh:
            for line, column in literal_reads(fh.read()):
                hits.append((os.path.relpath(path, ROOT), line, column))
    return hits


def run(selftest=False):
    """Run every check, print PASS/FAIL lines, return the exit code."""
    checks = []

    def check(name, ok, detail):
        checks.append((bool(ok), name, detail))
        return bool(ok)

    with open(REGISTER, encoding="utf-8") as fh:
        reg = json.load(fh)
    entry = reg["entries"][0]
    observed = entry["observed"]
    cts_id, cts_name = observed[ID_COLUMN], observed[NAME_COLUMN]

    # A. the register's pins -------------------------------------------------
    pinned = [entry["origin"]] + entry["files"]
    bad = []
    for item in pinned:
        path = os.path.join(ROOT, item["path"])
        if not os.path.exists(path):
            bad.append("%s: missing" % item["path"])
            continue
        if "sha256" in item and sha256(path) != item["sha256"]:
            bad.append("%s: sha256 differs from the register" % item["path"])
        if "rows" in item:
            n = (len(read_pipe(path)) if path.endswith(".txt")
                 else len(read_csv(path)))
            if n != item["rows"]:
                bad.append("%s: %d rows, register says %d"
                           % (item["path"], n, item["rows"]))
    check("the files the register names are byte-identical to it",
          not bad,
          "; ".join(bad) or "%d files, sha256 and row count match" % len(pinned))

    # B. the register's claim -----------------------------------------------
    txt = read_pipe(os.path.join(ROOT, entry["origin"]["path"]))
    tables = [("cts_measurements_full.txt", txt),
              ("cts_ocv_channels.csv", read_csv(core.CSV_PATH)),
              ("support/cts_ocv_channels.csv", read_csv(SUPPORT_CSV))]
    for label, rows in tables:
        triple = (cts_id, cts_name, observed[TYPE_COLUMN])
        held = [r for r in rows
                if (r.get(ID_COLUMN), r.get(NAME_COLUMN), r.get(TYPE_COLUMN))
                == triple]
        check("%s: every row carries the quarantined triple" % label,
              rows and len(held) == len(rows),
              "%d/%d rows: %s + %r + %r"
              % (len(held), len(rows), cts_id, cts_name, triple[2]))

    # C. the identity rule ---------------------------------------------------
    rows = tables[1][1]
    sup = tables[2][1]
    keys = [core.identity_key(r) for r in rows]
    check("identity_key is unique over the table",
          len(set(keys)) == len(keys),
          "%d rows, %d distinct track|day|role keys"
          % (len(keys), len(set(keys))))
    check("identity_key is unique over the support table",
          len(set(core.identity_key(r) for r in sup)) == len(sup),
          "%d rows, %d distinct keys"
          % (len(sup), len(set(core.identity_key(r) for r in sup))))
    poisoned = []
    for r in rows:
        q = dict(r)
        for column in ("asset_id",) + core.QUARANTINED_COLUMNS:
            q[column] = "poisoned"
        poisoned.append(core.identity_key(q))
    check("rewriting the asset columns leaves every identity unchanged",
          poisoned == keys,
          "poison test: %d keys before and after %s = 'poisoned'"
          % (len(keys), "/".join(("asset_id",) + core.QUARANTINED_COLUMNS)))
    clash = [k for k in ("track", "day", "role")
             if k in set(core.QUARANTINED_COLUMNS) | set(core.MISNAMED_COLUMNS)]
    check("no identity component is itself quarantined", not clash,
          "identity = track|day|role, quarantined = %s, misnamed = %s"
          % (", ".join(core.QUARANTINED_COLUMNS), ", ".join(core.MISNAMED_COLUMNS)))

    # D. the cited evidence ---------------------------------------------------
    seg = read_pipe(SEGMENTS)
    seg_label = seg[0].get(NAME_COLUMN) if seg else None
    check("cts_segments.txt labels that asset_id as Champlain Towers South",
          bool(seg) and seg[0].get(ID_COLUMN) == cts_id
          and seg_label == "Champlain Towers South",
          "asset_id %s -> %r" % (cts_id, seg_label))
    with open(core.REF_ASSET, encoding="utf-8") as fh:
        asset = json.load(fh)
    check("cts_asset.json names the site and asserts no id from the export",
          asset.get("label") == "Champlain Towers South"
          and asset.get(ID_COLUMN) in (None, ""),
          "label %r, %s %r" % (asset.get("label"), ID_COLUMN,
                               asset.get(ID_COLUMN)))
    with open(WERDER_META, encoding="utf-8") as fh:
        werder = json.load(fh)
    check("the same name is Werder's, under a different id",
          werder.get(NAME_COLUMN) == cts_name
          and werder.get(ID_COLUMN) not in (None, "", cts_id),
          "Werder: %r + %s; CTS rows carry %s"
          % (werder.get(NAME_COLUMN), werder.get(ID_COLUMN), cts_id))

    # E. the code guard -------------------------------------------------------
    hits = scan_package()
    modules = [n for n in os.listdir(HERE) if n.endswith(".py")]
    check("no literal read of a quarantined column in code/cts/", not hits,
          "; ".join("%s:%d %s" % h for h in hits) or
          "AST scan of %d modules: literal subscripts and .get() calls"
          % len(modules))
    drift = [r for r in rows if r.get(REQUEST_COLUMN) != r.get("track")]
    check("request_id is the track under a misleading name", not drift,
          "%d rows, %d differ from `track`" % (len(rows), len(drift)))

    # F. the checker's own controls (selftest only) ---------------------------
    if selftest:
        found = dict(literal_reads('row["asset_name"]\n'
                                   'x = row.get("asset_type")\n'
                                   'y = row.get("track")\n'
                                   'z = row["day"]\n'))
        check("the scan flags a literal read of a quarantined column",
              found.get(1) == "asset_name" and found.get(2) == "asset_type",
              "synthetic row[\"asset_name\"] -> %r, row.get(\"asset_type\") -> %r"
              % (found.get(1), found.get(2)))
        check("the scan leaves the identity columns alone",
              found.get(3) is None and found.get(4) is None,
              "synthetic row.get(\"track\") and row[\"day\"] not flagged")
        try:
            core.identity_key({"track": "", "day": "2021-06-24", "role": "target"})
            raised = False
        except KeyError:
            raised = True
        check("identity_key refuses a row without a footprint", raised,
              "no track -> KeyError, never a silent fallback to the asset columns")
        try:
            core.assert_identity_unique([{"track": "t", "day": "d", "role": "r"}] * 2)
            raised = False
        except ValueError:
            raised = True
        check("assert_identity_unique refuses a duplicated row", raised,
              "two identical keys -> ValueError")

    failed = sum(1 for ok, _, _ in checks if not ok)
    for ok, name, detail in checks:
        print("%s  %-60s %s" % ("PASS" if ok else "FAIL", name, detail))
    print("\n%d checks, %d failed" % (len(checks), failed))
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="cts_ocv_quarantine.py",
        description="Check the asset-metadata quarantine against the record: the "
                    "register's pins, the wrong values it names, the identity rule "
                    "that replaces them, and the code that must not read them.")
    parser.add_argument("--check", action="store_true",
                        help="run the checks and print them (the default)")
    parser.add_argument("--selftest", action="store_true",
                        help="also exercise the checker's own controls")
    args = parser.parse_args(argv)
    return run(selftest=args.selftest)


if __name__ == "__main__":
    sys.exit(main())
