#!/usr/bin/env bash
# Runs code/espoo/fig_espoo_ocv_paper.py from the repository root, so its defaults
# (data/espoo/espoo_ocv_channels.csv in; data/espoo/fig_espoo_ocv_paper.json,
# figures/espoo/fig_espoo_ocv_paper.md and figures/espoo/fig_espoo_ocv_paper_[A-E].png
# out) resolve correctly. Any argument is passed through, e.g. --quick.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/espoo/fig_espoo_ocv_paper.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/espoo/fig_espoo_ocv_paper.py "$@"
