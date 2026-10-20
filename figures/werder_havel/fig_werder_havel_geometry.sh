#!/usr/bin/env bash
# Runs code/werder_havel/fig_werder_havel_geometry.py from the repository root, so its
# defaults (data/werder_havel/reference/werder_havel_geometry.json and ..._aoi.json in;
# data/werder_havel/fig_werder_havel_geometry.json,
# figures/werder_havel/fig_werder_havel_geometry.md and
# figures/werder_havel/fig_werder_havel_geometry.png out) resolve correctly. Any
# argument is passed through, e.g. --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed where
# the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/werder_havel/fig_werder_havel_geometry.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/werder_havel/fig_werder_havel_geometry.py "$@"
