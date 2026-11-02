#!/usr/bin/env python3
"""The typography and the print width of every figure, checked from the source
and from the written file.

A figure can be wrong in two independent ways, so two things are audited:

  (i)  **what a generator declares.**  Every ``code/**/*.py`` that writes a
       figure is parsed (AST, not grep): the ``rcParams`` block, the
       ``rcParams[...] = ...`` assignments and every ``fontsize=`` /
       ``labelsize=`` / ``titlesize=`` keyword are collected together with the
       call they belong to (``set_title`` -> title, ``set_xlabel`` -> label,
       ``tick_params`` -> tick, ``legend`` -> legend, ``text`` -> annotation)
       and with the line number.  Values are resolved where they are
       resolvable: a literal, a module-level constant of the same file, or
       ``FS["tick"]`` from :mod:`figure_style`.  Each role then has to sit in
       its band (title 8.5-9.0, label 7.0-8.0, tick and legend 6.5-7.0,
       annotations >= 6.0 pt, ``font.size`` 8.0).
  (ii) **what the file really is.**  Point sizes are authored, not printed: a
       canvas that ``bbox_inches="tight"`` widens after the layout shrinks
       every glyph with it.  So the saved width of every artefact is measured
       (PNG ``pHYs`` / JPEG ``JFIF`` / PDF ``/MediaBox`` / EPS
       ``%%HiResBoundingBox``, stdlib only, see :mod:`figure_style`) and has to
       be 170 +/- 0.5 mm.  The effective size at 170 mm is printed next to each
       declared value, so a 9.0 pt title on a 375 mm canvas is visible as the
       4.1 pt a reader would actually see.

A plate is one row, not one per container: the raster, the ``.pdf`` and the
``.eps`` that :func:`figure_style.save` writes under one stem are listed
together, and **every** container is judged on its own width -- a plate can be
170.0 mm as a JPEG and cropped as an EPS.  A vector container is judged on its
box, never on a density it cannot have.

Scope.  The site directories in ``SCOPE_DIRS`` are the figures of this project
and are failed on.  The wall plates (``WALL_PLATES``, empty now that the
oversized panels are retired) are measured and listed but not failed while
their final width is an open decision, and ``figures/concept/`` and
``figures/sds/`` are documentation plates that are measured and reported only.
``--strict`` fails on those two groups as well.

One exception runs the other way.  The plates of ``PAGE_PLATES`` are authored on
their own page -- 183 mm, three panels of one isotropic unit -- and are
deliberately *not* the 170 mm print plate.  They are not exempt from the width
check, only from that one number: every container of such a plate has to agree
with the declared page width, so a cropped EPS or a padded JPEG still fails, and
the row says which page it was measured against.

Exit status: 0 when every check passes, 1 on any violation (the convention of
``audit_ocv_units.py``).

Usage:
  python3 code/audit_figure_style.py                 # audit, print, no write
  python3 code/audit_figure_style.py --only kdlo     # one site (repeatable)
  python3 code/audit_figure_style.py --strict        # wall plates fail too
  python3 code/audit_figure_style.py --json /tmp/style.json   # machine-readable
  python3 code/audit_figure_style.py --selftest      # the rules on synthetic input
"""
from __future__ import annotations

import argparse
import ast
import json
import os
import sys

import figure_style as style

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGURES = os.path.join(ROOT, "figures")

#: Directories whose figures are failed on (the site figures of this repo).
SCOPE_DIRS = ("bautzen", "carola", "cases", "cts", "espoo", "haltern", "kdlo",
              "lumo", "morandi", "werder_havel")

#: Documentation plates: measured and reported, never failed (unless --strict).
DOC_DIRS = ("concept", "sds")

#: Plates whose final width is an open decision (the oversized wall panels).
#: Empty: the wall panels were retired with the twelve-plate scope, so every
#: remaining plate is judged on its width.  The hook stays -- ``--strict`` and
#: the per-plate row layout are defined in terms of it.
WALL_PLATES = ()

#: Plates authored on their own page width, not on the print width, keyed by
#: figure stem and given that page width in millimetres.
#:
#: The two concept plates are 183 mm by design: three equal panels whose one
#: length unit is 55.632 mm, so the panels stay isotropic in both directions.
#: Cutting them to 170 mm would either squash the geometry or shrink the 7 pt
#: floor the plate is built on, so the page -- not the print width -- is the
#: promise here.  It is still a promise: the group is failed on like a site
#: plate, only against its own width instead of 170 mm.
PAGE_PLATES = {
    "fig_ocv_schematic": 183.0,
    "fig_ocv_measurement": 183.0,
}

IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg")

#: Everything ``figure_style.save`` writes: the raster plus its vector sidecars.
ARTEFACT_SUFFIXES = IMAGE_SUFFIXES + (".pdf", ".eps")

#: ``fontsize=`` on one of these calls carries that role.
CALLEE_ROLES = {
    "set_title": "title", "title": "title", "suptitle": "suptitle",
    "set_suptitle": "suptitle",
    "set_xlabel": "label", "set_ylabel": "label", "xlabel": "label",
    "ylabel": "label", "supxlabel": "label", "supylabel": "label",
    "legend": "legend", "figlegend": "legend",
    "text": "annotation", "annotate": "annotation", "table": "annotation",
    "set_xticklabels": "tick", "set_yticklabels": "tick",
    "set_ticklabels": "tick", "tick_params": "tick",
}

#: Calls where ``labelsize=`` means the tick labels (otherwise: the axis label).
TICK_CALLEES = ("tick_params", "set_xticklabels", "set_yticklabels",
                "set_ticklabels")


# ---------------------------------------------------------------------------
# (i) What a generator declares
# ---------------------------------------------------------------------------
def _callee_name(func):
    """The bare name of a call's target (``ax.set_title`` -> ``set_title``)."""
    if isinstance(func, ast.Attribute):
        return func.attr
    if isinstance(func, ast.Name):
        return func.id
    return ""


def _number(node):
    """A literal, a simple product/quotient of literals, or ``None``."""
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool):
            return None
        if isinstance(node.value, (int, float)):
            return float(node.value)
        return None
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        inner = _number(node.operand)
        if inner is None:
            return None
        return -inner if isinstance(node.op, ast.USub) else inner
    if isinstance(node, ast.BinOp):
        left, right = _number(node.left), _number(node.right)
        if left is None or right is None:
            return None
        if isinstance(node.op, ast.Div):
            return left / right if right else None
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
    return None


def _module_constants(tree):
    """Numeric module-level constants of one file (``TICK_PT = 6.8``)."""
    consts = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        value = _number(node.value)
        if value is not None:
            consts[target.id] = value
    return consts


def _resolve(node, consts):
    """The point size behind an expression, or ``None`` if it is symbolic.

    Literals and the file's own numeric constants are resolved, and so is
    ``FS["tick"]`` (or ``figure_style.FS["tick"]``): a generator that reads the
    shared ladder is checked against the ladder's own value.
    """
    value = _number(node)
    if value is not None:
        return value
    if isinstance(node, ast.Name) and node.id in consts:
        return consts[node.id]
    if isinstance(node, ast.Subscript):
        base = node.value
        name = base.id if isinstance(base, ast.Name) else _callee_name(base)
        key = node.slice
        if name == "FS" and isinstance(key, ast.Constant) \
                and isinstance(key.value, str) and key.value in style.FS:
            return style.FS[key.value]
    return None


def _role_of_kwarg(kwarg, callee):
    """The role a ``fontsize``/``labelsize``/``titlesize`` keyword carries."""
    if kwarg == "titlesize":
        return "title"
    if kwarg == "labelsize":
        return "tick" if callee in TICK_CALLEES else "label"
    return CALLEE_ROLES.get(callee, "other")


def scan_source(path):
    """Every declared font size of one generator, by role and line.

    Returns ``{"roles": {role: [entry, ...]}, "unresolved": [...],
    "writes_figure": bool}``; an entry is
    ``{"pt": float|None, "line": int, "where": str, "text": str}``.
    """
    with open(path, encoding="utf-8") as fh:
        source = fh.read()
    tree = ast.parse(source, filename=path)
    consts = _module_constants(tree)
    roles, unresolved, writes_figure = {}, [], False

    def add(role, pt, node, where, text):
        entry = {"role": role, "pt": pt, "line": getattr(node, "lineno", 0),
                 "where": where, "text": text}
        if pt is None:
            unresolved.append(entry)
        else:
            roles.setdefault(role, []).append(entry)

    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            # rcParams["axes.labelsize"] = 7.6
            for target in node.targets:
                if (isinstance(target, ast.Subscript)
                        and _callee_name(target.value) == "rcParams"
                        and isinstance(target.slice, ast.Constant)
                        and target.slice.value in style.RC_ROLES):
                    key = target.slice.value
                    add(style.RC_ROLES[key], _resolve(node.value, consts), node,
                        "rcParams[%r]" % key, "rcParams[%r]" % key)
            continue
        if not isinstance(node, ast.Call):
            continue
        callee = _callee_name(node.func)
        if callee == "savefig":
            writes_figure = True
        elif (callee in ("save", "save_eps") and isinstance(node.func, ast.Attribute)
              and _callee_name(node.func.value) in ("style", "figure_style")):
            # A migrated generator hands the write to the shared module
            # (``style.save(fig, out)``), so it has no ``savefig`` of its own --
            # without this the audit would skip exactly the plates that follow
            # the ladder, and check only the ones that do not.
            writes_figure = True
        if (callee == "update" and isinstance(node.func, ast.Attribute)
                and _callee_name(node.func.value) == "rcParams"
                and node.args and isinstance(node.args[0], ast.Dict)):
            for key, value in zip(node.args[0].keys, node.args[0].values):
                if isinstance(key, ast.Constant) and key.value in style.RC_ROLES:
                    add(style.RC_ROLES[key.value], _resolve(value, consts), node,
                        "rcParams[%r]" % key.value, "rcParams[%r]" % key.value)
        if callee in ("apply", "rcparams") and isinstance(node.func, ast.Attribute):
            if _callee_name(node.func.value) in ("style", "figure_style"):
                for role in style.CHECKED_ROLES + ("base",):
                    add(role, style.FS[role], node,
                        "figure_style.%s()" % callee, "the shared ladder")
        for kw in node.keywords:
            if kw.arg in ("fontsize", "labelsize", "titlesize"):
                role = _role_of_kwarg(kw.arg, callee)
                add(role, _resolve(kw.value, consts), node,
                    "%s(%s=)" % (callee or "?", kw.arg),
                    "%s(%s=...)" % (callee or "?", kw.arg))
    return {"roles": roles, "unresolved": unresolved,
            "writes_figure": writes_figure}


def generators():
    """Relative paths of every source under ``code/`` that writes a figure."""
    found = []
    for dirpath, dirnames, filenames in os.walk(HERE):
        dirnames[:] = sorted(d for d in dirnames if d != "__pycache__")
        for name in sorted(filenames):
            if not name.endswith(".py"):
                continue
            path = os.path.join(dirpath, name)
            try:
                writes = scan_source(path)["writes_figure"]
            except (SyntaxError, OSError):
                writes = False
            if writes:
                found.append(os.path.relpath(path, ROOT))
    return sorted(found)


def role_verdict(role, pt):
    """Why *pt* is outside *role*'s band, or ``None`` when it is inside."""
    if role == "base":
        if abs(pt - style.FS["base"]) < 1e-9:
            return None
        return "base is the anchor %.1f pt, not a free value" % style.FS["base"]
    if role not in style.BANDS:
        return None  # e.g. a fontsize on a call with no print role
    low, high = style.BANDS[role]
    if low is not None and pt < low - 1e-9:
        return "below %.1f pt" % low
    if high is not None and pt > high + 1e-9:
        return "above %.1f pt" % high
    return None


# ---------------------------------------------------------------------------
# (ii) What the written files really are
# ---------------------------------------------------------------------------
def artefacts(only=()):
    """Every written plate, with its group and the width of every container.

    One row per *stem*, not per file: ``jpg``, ``pdf`` and ``eps`` are three
    containers of one plate, and the report has to say so -- otherwise the plate
    list triples the moment the vector sidecars are written.
    """
    rows = []
    for dirpath, dirnames, filenames in os.walk(FIGURES):
        dirnames.sort()
        site = os.path.relpath(dirpath, FIGURES).split(os.sep)[0]
        if only and site not in only:
            continue
        plates = {}
        for name in sorted(filenames):
            stem, suffix = os.path.splitext(name)
            if suffix.lower() not in ARTEFACT_SUFFIXES:
                continue
            plates.setdefault(stem, []).append(
                (suffix.lower(), os.path.join(dirpath, name)))
        for stem in sorted(plates):
            if stem in WALL_PLATES:
                group = "wall"
            elif stem in PAGE_PLATES:
                group = "page"
            elif site in DOC_DIRS:
                group = "doc"
            else:
                group = "scope"
            row = {"site": site, "stem": stem, "group": group,
                   "expect_mm": PAGE_PLATES.get(stem, style.FINAL_WIDTH_MM),
                   "path": None, "formats": {}}
            for suffix, path in plates[stem]:
                entry = {"path": os.path.relpath(path, ROOT)}
                try:
                    info = style.measure(path)
                except (ValueError, OSError) as exc:
                    entry.update({"mm": None, "dpi": None, "px": None,
                                  "vector": suffix in (".pdf", ".eps"),
                                  "error": str(exc)})
                else:
                    entry.update({"mm": info["mm"], "dpi": info["dpi"],
                                  "px": info["px"], "vector": info["vector"]})
                row["formats"][suffix] = entry
                raster = not entry["vector"]
                if row["path"] is None or raster:
                    row["path"] = entry["path"]
            rasters = [s for s in sorted(row["formats"])
                       if not row["formats"][s]["vector"]]
            shown = row["formats"][rasters[0] if rasters
                                   else sorted(row["formats"])[0]]
            row["mm"] = shown.get("mm")
            row["error"] = shown.get("error")
            rows.append(row)
    return rows


def _width_verdict(mm, expect_mm=style.FINAL_WIDTH_MM):
    """Why *mm* is off *expect_mm* (the print width by default), or ``None``."""
    if mm is None:
        return "no print width in the file"
    off = mm - expect_mm
    if abs(off) > style.WIDTH_TOL_MM:
        return "%+.1f mm off %.0f mm" % (off, expect_mm)
    return None


def _plate_verdict(row):
    """Why a plate is off its width, naming the container, or ``None``.

    Every container is judged on its own: a JPEG at 170.0 mm says nothing about
    the EPS beside it, which can be cropped, padded, or -- in the worst case --
    carry no box at all, so its width is not unknown but wrong.  A page plate
    (``PAGE_PLATES``) is judged the same way, against the page it declares.
    """
    expect = row.get("expect_mm", style.FINAL_WIDTH_MM)
    for suffix, info in sorted(row["formats"].items()):
        name = suffix.lstrip(".")
        if info.get("error"):
            return "%s: %s" % (name, info["error"])
        problem = _width_verdict(info.get("mm"), expect)
        if problem:
            return "%s %s" % (name, problem)
    return None


def _formats_summary(formats):
    """``jpg+pdf+eps``, or ``jpg 170.0 + eps 165.9`` when they disagree."""
    def rank(suffix):
        order = (".png", ".jpg", ".jpeg", ".pdf", ".eps")
        return order.index(suffix) if suffix in order else len(order)

    groups = {}
    for suffix in sorted(formats, key=rank):
        info = formats[suffix]
        width = None if info.get("mm") is None else round(info["mm"], 1)
        groups.setdefault(width, []).append(suffix.lstrip("."))
    parts = []
    for width in sorted(groups, key=lambda w: (w is not None, w)):
        names = "+".join(groups[width])
        parts.append(names if width is None or len(groups) == 1
                     else "%s %.1f" % (names, width))
    return " + ".join(parts)


def _role_order(role):
    """The print order of the roles (the ladder's own order)."""
    order = ("base", "title", "suptitle", "label", "tick", "legend",
             "annotation", "other")
    return order.index(role) if role in order else len(order)


def audit(only=(), strict=False):
    """The whole audit: every generator's type, every artefact's width."""
    sources = generators()
    by_stem = {}
    for rel in sources:
        by_stem.setdefault(os.path.splitext(os.path.basename(rel))[0], rel)
    rows, totals = [], {"generators": 0, "type_violations": 0,
                        "artefacts": 0, "containers": 0,
                        "width_violations": 0, "unresolved": 0}
    for rel in sources:
        site = os.path.dirname(rel).split(os.sep)[-1]
        if only and site not in only:
            continue
        base = os.path.basename(rel)
        if base.startswith("audit_") or rel == "code/figure_style.py":
            continue  # the audit's own helpers, not a figure generator
        scan = scan_source(os.path.join(ROOT, rel))
        verdicts = []
        for role in sorted(scan["roles"], key=_role_order):
            for entry in scan["roles"][role]:
                verdicts.append({"role": role, "pt": entry["pt"],
                                 "line": entry["line"],
                                 "where": entry["where"],
                                 "problem": role_verdict(role, entry["pt"])})
        totals["generators"] += 1
        totals["unresolved"] += len(scan["unresolved"])
        if any(v["problem"] for v in verdicts):
            totals["type_violations"] += 1
        rows.append({"path": rel, "group": "source", "verdicts": verdicts,
                     "unresolved": scan["unresolved"]})
    arts = artefacts(only=only)
    for row in arts:
        # A page plate is held to its own page, so it is judged like a site
        # plate -- under --strict too, where it would otherwise fall through
        # the "wall or doc" branch and be reported without a verdict.
        fails = row["group"] in ("scope", "page") or (
            strict and row["group"] in ("wall", "doc"))
        row["problem"] = _plate_verdict(row) if fails else None
        if row["problem"]:
            totals["width_violations"] += 1
        row["generator"] = by_stem.get(row["stem"])
    totals["artefacts"] = len(arts)
    totals["containers"] = sum(len(row["formats"]) for row in arts)
    return {"sources": rows, "artefacts": arts, "totals": totals}


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def _compact(role, entries):
    """``tick 6.8``, or ``tick 6.5-7.0`` when a role carries several sizes."""
    values = sorted({round(e["pt"], 2) for e in entries})
    if len(values) == 1:
        return "%s %.1f" % (role, values[0])
    return "%s %.1f-%.1f" % (role, values[0], values[-1])


def report_lines(result, scans=None):
    """The audit as printable lines (PASS/FAIL per generator and artefact)."""
    scans = scans or {}
    out = ["== the declared type (parsed from the sources) =="]
    for row in result["sources"]:
        verdicts = row["verdicts"]
        bad = [v for v in verdicts if v["problem"]]
        status = "FAIL" if bad else ("NOTE" if not verdicts else "PASS")
        parts = []
        for role in ("base", "title", "suptitle", "label", "tick", "legend",
                     "annotation"):
            entries = [v for v in verdicts if v["role"] == role]
            if entries:
                parts.append(_compact(role, entries))
        out.append("%s  %-50s %s" % (status, row["path"],
                                     "  ".join(parts) if parts
                                     else "no size declared"))
        for v in bad:
            out.append("      L%-5d %-10s %.1f pt  %-24s %s!"
                       % (v["line"], v["role"], v["pt"], v["where"],
                          v["problem"]))
        for u in row["unresolved"]:
            out.append("      L%-5d %-10s symbolic  %s"
                       % (u["line"], u["role"], u["where"]))
    out.append("")
    out.append("== the saved width, and the type a reader sees at %.0f mm =="
               % style.FINAL_WIDTH_MM)
    for row in result["artefacts"]:
        mm = row.get("mm")
        label = os.path.relpath(row["path"], "figures")
        if mm is None:
            out.append("FAIL  %-50s %s"
                       % (label, row.get("error") or "no print width"))
            continue
        status = "FAIL" if row.get("problem") else (
            "PASS" if row["group"] == "scope" else "note")
        page = ""
        if row.get("expect_mm", style.FINAL_WIDTH_MM) != style.FINAL_WIDTH_MM:
            page = "  page %.0f mm" % row["expect_mm"]
        eff = ""
        scan = scans.get(row.get("generator"))
        if scan:
            pieces = []
            for role in ("title", "label", "tick", "legend"):
                entries = scan["roles"].get(role)
                if entries:
                    pt = style.effective_pt(entries[0]["pt"], mm)
                    pieces.append("%s %.1f" % (role, pt))
            eff = "  -> " + "  ".join(pieces) if pieces else ""
        out.append("%s  %-50s %6.1f mm  %-16s%s%s"
                   % (status, label, mm, _formats_summary(row["formats"]),
                      page, eff))
        if row.get("problem") and row["group"] != "scope":
            out.append("      %s -- %s plate, %s"
                       % (row["problem"], row["group"],
                          "held to its own page" if row["group"] == "page"
                          else "reported only"))
    totals = result["totals"]
    out.append("")
    out.append("%s, %d plate(s) in %d file(s): %d type violation(s), "
               "%d width violation(s), %d symbolic size(s)"
               % (_plural(totals["generators"], "generator"),
                  totals["artefacts"], totals["containers"],
                  totals["type_violations"], totals["width_violations"],
                  totals["unresolved"]))
    return out


def _plural(count, word):
    """``1 generator`` / ``24 generators``."""
    return "%d %s%s" % (count, word, "" if count == 1 else "s")


# ---------------------------------------------------------------------------
# The self-test: the rules on synthetic sources and synthetic files
# ---------------------------------------------------------------------------
def _synth_png(width_px, dpi):
    """A 1-pixel-high PNG that declares *dpi* in its ``pHYs`` chunk."""
    import struct
    import zlib

    def chunk(kind, payload):
        crc = zlib.crc32(kind + payload) & 0xFFFFFFFF
        return (struct.pack(">I", len(payload)) + kind + payload
                + struct.pack(">I", crc))

    ppm = int(round(dpi / 0.0254))
    ihdr = struct.pack(">IIBBBBB", width_px, 1, 8, 0, 0, 0, 0)
    phys = struct.pack(">IIB", ppm, ppm, 1)
    raw = zlib.compress(b"\x00" + b"\x00" * width_px, 9)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"pHYs", phys)
            + chunk(b"IDAT", raw) + chunk(b"IEND", b""))


def _synth_jpeg(width_px, height_px, dpi, unit=1):
    """A header-only JPEG (APP0 + SOF0 + SOS): enough for the reader."""
    import struct
    app0 = (b"JFIF\x00" + b"\x01\x01" + bytes([unit])
            + struct.pack(">HH", dpi, dpi) + b"\x00\x00")
    sof = (b"\x08" + struct.pack(">HH", height_px, width_px) + b"\x03"
           + b"\x01\x11\x00\x02\x11\x00\x03\x11\x00")
    sos = b"\x03\x01\x00\x02\x11\x03\x11\x00\x3f\x00"
    return (b"\xff\xd8"
            + b"\xff\xe0" + struct.pack(">H", len(app0) + 2) + app0
            + b"\xff\xc0" + struct.pack(">H", len(sof) + 2) + sof
            + b"\xff\xda" + struct.pack(">H", len(sos) + 2) + sos
            + b"\xff\xd9")


def _synth_eps(width_pt, height_pt=None, hires=True):
    """A header-only EPS: exactly the DSC lines the reader looks for."""
    height_pt = width_pt if height_pt is None else height_pt
    head = ["%!PS-Adobe-3.0 EPSF-3.0"]
    if hires:
        head.append("%%%%HiResBoundingBox: 0 0 %r %r" % (width_pt, height_pt))
    head.append("%%%%BoundingBox: 0 0 %d %d"
                % (round(width_pt), round(height_pt)))
    return ("\n".join(head) + "\n%%EndComments\nshowpage\n").encode("latin-1")


def _synth_pdf(width_pt, height_pt=None, key="/MediaBox"):
    """A PDF whose first page object carries *key* with the given box."""
    height_pt = width_pt if height_pt is None else height_pt
    return ("%%PDF-1.4\n1 0 obj\n<< /Type /Page %s [0 0 %r %r] >>\nendobj\n"
            "trailer\n<< /Root 1 0 R >>\n%%%%EOF\n"
            % (key, width_pt, height_pt)).encode("latin-1")


def selftest():
    """Check the reader and the classifier on input whose answer is known."""
    import tempfile

    checks = []

    def check(name, ok, detail=""):
        checks.append((bool(ok), name, detail))

    tmpdir = tempfile.mkdtemp(prefix="figure_style_selftest_")
    png = os.path.join(tmpdir, "a.png")
    with open(png, "wb") as fh:
        fh.write(_synth_png(4016, 600))
    info = style.measure(png)
    check("a PNG's pHYs gives the print width", abs(info["mm"] - 170.0) < 0.05,
          "4016 px at 600 dpi = %.2f mm" % info["mm"])
    try:
        style.assert_width(png)
        ok, detail = True, "170 +/- 0.5 mm"
    except RuntimeError as exc:
        ok, detail = False, str(exc)[:58]
    check("the width assertion passes on a 170 mm file", ok, detail)

    jpg = os.path.join(tmpdir, "b.jpg")
    with open(jpg, "wb") as fh:
        fh.write(_synth_jpeg(4200, 1000, 630))
    info = style.measure(jpg)
    check("a JPEG's JFIF density gives the print width",
          abs(info["mm"] - 169.3) < 0.06 and info["px"] == (4200, 1000),
          "4200 px at 630 dpi = %.2f x %.2f mm"
          % (info["mm"], info["height_mm"]))
    wide = os.path.join(tmpdir, "c.jpg")
    with open(wide, "wb") as fh:
        fh.write(_synth_jpeg(5000, 1000, 630))
    try:
        style.assert_width(wide)
        ok, detail = False, "no exception"
    except RuntimeError as exc:
        ok, detail = "overflows the print width" in str(exc), "200 mm canvas"
    check("a 200 mm canvas is rejected by the width assertion", ok, detail)

    src = os.path.join(tmpdir, "gen.py")
    with open(src, "w", encoding="utf-8") as fh:
        fh.write(
            "TICK_PT = 5.8\n"
            "def draw(ax, fig):\n"
            "    plt.rcParams.update({'font.size': 8.0,\n"
            "                         'axes.titlesize': 8.8,\n"
            "                         'axes.labelsize': 7.6,\n"
            "                         'xtick.labelsize': 6.8,\n"
            "                         'legend.fontsize': 6.5})\n"
            "    ax.set_title('a panel')\n"
            "    ax.set_xlabel('day', fontsize=7.6)\n"
            "    ax.tick_params(labelsize=TICK_PT)\n"
            "    ax.legend(fontsize=FS['legend'])\n"
            "    ax.text(0.1, 0.9, 'x', fontsize=6.2)\n"
            "    fig.savefig('out.png')\n")
    scan = scan_source(src)
    check("the scanner sees the rcParams block",
          [e["pt"] for e in scan["roles"]["title"]] == [8.8]
          and [e["pt"] for e in scan["roles"]["base"]] == [8.0],
          "title 8.8, base 8.0")
    check("the scanner resolves a module constant",
          5.8 in [e["pt"] for e in scan["roles"]["tick"]], "TICK_PT = 5.8")
    check("the scanner resolves FS['legend'] from the shared ladder",
          sorted(e["pt"] for e in scan["roles"]["legend"])
          == [6.5, style.FS["legend"]],
          "6.5 from rcParams and %.1f from FS" % style.FS["legend"])
    check("tick_params(labelsize=) is a tick size, not a label",
          5.8 not in [e["pt"] for e in scan["roles"]["label"]]
          and [e["pt"] for e in scan["roles"]["label"]] == [7.6, 7.6],
          "both label sizes are the 7.6 of set_xlabel and rcParams")
    check("text(fontsize=) is an annotation",
          [e["pt"] for e in scan["roles"]["annotation"]] == [6.2], "6.2 pt")
    check("a call that writes a figure is recognised",
          scan["writes_figure"] is True, "fig.savefig")
    check("a 5.8 pt tick is a violation",
          role_verdict("tick", 5.8) == "below 6.5 pt", "below 6.5 pt")
    check("an 8.8 pt title is inside the band",
          role_verdict("title", 8.8) is None, "8.5-9.0 pt")
    check("a 13.4 pt title is a violation",
          role_verdict("title", 13.4) == "above 9.0 pt", "above 9.0 pt")
    check("a 5.5 pt annotation is a violation",
          role_verdict("annotation", 5.5) == "below 6.0 pt", "floor 6.0 pt")
    check("a 7.1 pt legend is a violation",
          role_verdict("legend", 7.1) == "above 7.0 pt", "6.5-7.0 pt")
    check("169.8 mm is inside the width tolerance",
          _width_verdict(169.8) is None, "170 +/- 0.5 mm")
    check("a canvas without a print width is a violation",
          _width_verdict(None) is not None, "no print width in the file")

    eps = os.path.join(tmpdir, "d.eps")
    with open(eps, "wb") as fh:
        fh.write(_synth_eps(481.889764))
    info = style.measure(eps)
    check("an EPS's HiResBoundingBox gives the print width",
          info["vector"] and abs(info["mm"] - 170.0) < 0.01,
          "481.889764 pt = %.2f mm" % info["mm"])
    check("a vector container is judged on its box, not on a dpi",
          info["dpi"] is None and _width_verdict(info["mm"]) is None,
          "no dpi, %.2f mm" % info["mm"])
    coarse = os.path.join(tmpdir, "e.eps")
    with open(coarse, "wb") as fh:
        fh.write(_synth_eps(482, hires=False))
    info = style.measure(coarse)
    check("an EPS without a HiRes box falls back to the whole-point box",
          abs(info["mm"] - 170.0) < 0.05, "482 pt = %.2f mm" % info["mm"])
    narrow = os.path.join(tmpdir, "f.eps")
    with open(narrow, "wb") as fh:
        fh.write(_synth_eps(400))
    try:
        style.assert_width(narrow)
        ok, detail = False, "no exception"
    except RuntimeError as exc:
        ok, detail = "overflows the print width" in str(exc), "141 mm canvas"
    check("a cropped EPS is rejected by the width assertion", ok, detail)
    boxless = os.path.join(tmpdir, "g.eps")
    with open(boxless, "wb") as fh:
        fh.write(b"%!PS-Adobe-3.0 EPSF-3.0\n%%EndComments\nshowpage\n")
    try:
        style.measure(boxless)
        ok, detail = False, "no exception"
    except ValueError as exc:
        ok, detail = "no bounding box" in str(exc), "no DSC box at all"
    check("an EPS with no bounding box has no width at all", ok, detail)

    pdf = os.path.join(tmpdir, "h.pdf")
    with open(pdf, "wb") as fh:
        fh.write(_synth_pdf(481.889764))
    info = style.measure(pdf)
    check("a PDF's /MediaBox gives the print width",
          info["vector"] and abs(info["mm"] - 170.0) < 0.01,
          "481.889764 pt = %.2f mm" % info["mm"])
    crop = os.path.join(tmpdir, "i.pdf")
    with open(crop, "wb") as fh:
        fh.write(_synth_pdf(595.276, 841.89, key="/CropBox"))
    info = style.measure(crop)
    check("a PDF without /MediaBox falls back to /CropBox",
          abs(info["mm"] - 210.0) < 0.01
          and abs(info["height_mm"] - 297.0) < 0.01,
          "A4 crop box = %.1f x %.1f mm" % (info["mm"], info["height_mm"]))
    check("a plate verdict names the container that is off width",
          _plate_verdict({"formats": {".jpg": {"mm": 169.9},
                                     ".eps": {"mm": 160.0}}})
          == "eps -10.0 mm off 170 mm", "the jpg is fine, the eps is not")
    check("a plate lists its containers raster first",
          _formats_summary({".eps": {"mm": 170.0}, ".jpg": {"mm": 170.0},
                            ".pdf": {"mm": 170.0}}) == "jpg+pdf+eps",
          "jpg+pdf+eps")
    check("a page plate is measured against its own page",
          _width_verdict(183.0, 183.0) is None
          and _width_verdict(183.0) == "+13.0 mm off 170 mm",
          "183 mm passes on a 183 mm page and fails on the print width")
    check("a page plate's containers are judged on that page",
          _plate_verdict({"expect_mm": 183.0,
                          "formats": {".jpg": {"mm": 183.0},
                                      ".eps": {"mm": 170.0}}})
          == "eps -13.0 mm off 183 mm"
          and _plate_verdict({"expect_mm": 183.0,
                              "formats": {".jpg": {"mm": 183.0},
                                          ".pdf": {"mm": 183.0}}}) is None,
          "a 170 mm EPS is still wrong on a 183 mm page")
    check("the page plates are the two 183 mm concept plates",
          PAGE_PLATES == {"fig_ocv_schematic": 183.0,
                          "fig_ocv_measurement": 183.0}
          and all(mm > style.FINAL_WIDTH_MM + style.WIDTH_TOL_MM
                  for mm in PAGE_PLATES.values()),
          "183 mm is a declared page, not a width violation")

    shared = os.path.join(tmpdir, "shared.py")
    with open(shared, "w", encoding="utf-8") as fh:
        fh.write("import figure_style as style\n"
                 "def draw(fig):\n"
                 "    style.apply()\n"
                 "    fig.savefig('x.png')\n")
    scan = scan_source(shared)
    check("style.apply() is audited as the shared ladder",
          [e["pt"] for e in scan["roles"]["tick"]] == [style.FS["tick"]]
          and all(role_verdict(role, style.FS[role]) is None
                  for role in ("title", "label", "tick", "legend")),
          "the ladder passes its own bands")
    handed = os.path.join(tmpdir, "shared_save.py")
    with open(handed, "w", encoding="utf-8") as fh:
        fh.write("import figure_style as style\n"
                 "def draw(fig):\n"
                 "    style.apply()\n"
                 "    style.save(fig, 'x.jpg')\n")
    check("style.save() is audited as a figure writer",
          scan_source(handed)["writes_figure"] is True,
          "no savefig of its own, only the shared write")

    failed = sum(1 for ok, _, _ in checks if not ok)
    for ok, name, detail in checks:
        print("%s  %-56s %s" % ("PASS" if ok else "FAIL", name, detail))
    print("\n%d checks, %d failed" % (len(checks), failed))
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="audit_figure_style.py",
        description="Audit the typography of every figure generator and the "
                    "print width of every written figure: titles 8.5-9.0 pt, "
                    "labels 7.0-8.0 pt, ticks and legends 6.5-7.0 pt, "
                    "annotations >= 6.0 pt, canvas %.0f +/- %.1f mm."
                    % (style.FINAL_WIDTH_MM, style.WIDTH_TOL_MM))
    parser.add_argument("--only", action="append", default=[], metavar="SITE",
                        help="restrict the audit to one site directory "
                             "(repeatable)")
    parser.add_argument("--strict", action="store_true",
                        help="fail on the wall plates and the documentation "
                             "plates as well (a page plate in PAGE_PLATES is "
                             "always held to its own declared width)")
    parser.add_argument("--json", metavar="PATH",
                        help="write the audit as JSON to PATH")
    parser.add_argument("--selftest", action="store_true",
                        help="check the reader and the classifier on synthetic "
                             "input, read no figure")
    args = parser.parse_args(argv)
    if args.selftest:
        return selftest()

    result = audit(only=tuple(args.only), strict=args.strict)
    scans = {}
    for row in result["sources"]:
        scans[row["path"]] = scan_source(os.path.join(ROOT, row["path"]))
    for line in report_lines(result, scans):
        print(line)
    totals = result["totals"]
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print("\nwrote %s" % args.json)
    return 1 if (totals["type_violations"] or totals["width_violations"]) else 0


if __name__ == "__main__":
    sys.exit(main())

