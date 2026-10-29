#!/usr/bin/env bash
# Runs code/concept/fig_ocv_schematic.py from the repository root, so its defaults
# (no input file at all; data/concept/fig_ocv_schematic.json,
# figures/concept/fig_ocv_schematic.md, .jpg and .pdf out) resolve correctly.
# Any argument is passed through, e.g. --figdir /tmp/x or --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/concept/fig_ocv_schematic.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/concept/fig_ocv_schematic.py "$@"
