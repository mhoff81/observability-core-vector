#!/usr/bin/env bash
# Runs code/cases/ocv_trajectory.py --site espoo from the repository root, so its
# defaults (data/espoo/support/espoo_ocv_channels.csv in;
# data/espoo/fig_espoo_ocv_trajectory.json, figures/espoo/
# fig_espoo_ocv_trajectory.md and .jpg out) resolve correctly. Any argument is
# passed through, e.g. --verify.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/espoo/fig_espoo_ocv_trajectory.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/cases/ocv_trajectory.py --site espoo "$@"
