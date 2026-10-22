#!/usr/bin/env python3
"""Carola entry point for the shared decomposition module.

The method lives in ``code/orbits/ocv_decompose.py``, where the site registry
keeps one description per site (table, event, geometry rule, targets, model
tree). This wrapper only pins the site, so the commands and the output paths
documented in ``code/carola/README.md`` keep working:

  python3 code/carola/carola_ocv_decompose.py            # writes JSON + MD
  python3 code/carola/carola_ocv_decompose.py --selftest
  python3 code/carola/carola_ocv_decompose.py --verify

Everything else -- the two ladders, the moving-block bootstrap, the estimability
rule ``n_post >= 3 k``, the Holm families -- is the shared module's.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "code", "orbits"))

import ocv_decompose as core  # noqa: E402


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if "--site" not in arguments:
        arguments = ["--site", "carola"] + arguments
    return core.main(arguments)


if __name__ == "__main__":
    sys.exit(main())
