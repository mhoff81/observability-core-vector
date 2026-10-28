#!/usr/bin/env bash
# Runs code/cases/ocv_trajectory.py --site carola from the repository root, so its
# defaults (data/carola/support/carola_ocv_channels.csv in;
# data/carola/fig_carola_ocv_trajectory.json, figures/carola/
# fig_carola_ocv_trajectory.md and .jpg out) resolve correctly. Any argument is
# passed through, e.g. --verify or --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/carola/fig_carola_ocv_trajectory.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/cases/ocv_trajectory.py --site carola "$@"
