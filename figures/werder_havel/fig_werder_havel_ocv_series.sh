#!/usr/bin/env bash
# Runs code/werder_havel/fig_werder_havel_ocv_series.py from the repository root, so its
# defaults (data/werder_havel/werder_havel_measurements_full.txt,
# data/werder_havel/werder_havel_segments.txt and
# data/werder_havel/reference/werder_havel_aoi.json in;
# data/werder_havel/fig_werder_havel_ocv_series.json,
# figures/werder_havel/fig_werder_havel_ocv_series.md and
# figures/werder_havel/fig_werder_havel_ocv_series.png out) resolve correctly. Any
# argument is passed through, e.g. --check (self-checks and pins, nothing written) or
# --no-figures (the JSON and the .md only).
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed where
# the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/werder_havel/fig_werder_havel_ocv_series.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/werder_havel/fig_werder_havel_ocv_series.py "$@"
