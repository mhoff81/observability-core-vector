#!/usr/bin/env bash
# Runs code/cases/ocv_trajectory.py --site haltern from the repository root, so its
# defaults (data/haltern/support/haltern_ocv_channels.csv in;
# data/haltern/fig_haltern_ocv_trajectory.json, figures/haltern/
# fig_haltern_ocv_trajectory.md and the plate -- .jpg, .pdf and .eps -- out)
# resolve correctly. The plate is one shared frame, 170 x 156 mm, from
# code/figure_style.py's TRAJECTORY table, and its .pdf/.eps are the same
# drawing as vector. Any argument is passed through, e.g. --verify or
# --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/haltern/fig_haltern_ocv_trajectory.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/cases/ocv_trajectory.py --site haltern "$@"
