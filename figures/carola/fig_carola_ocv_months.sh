#!/usr/bin/env bash
# Runs code/carola/fig_carola_ocv_months.py from the repository root, so its
# defaults (data/carola/carola_ocv_channels.csv in; data/carola/fig_carola_ocv_months.json,
# figures/carola/fig_carola_ocv_months.md and figures/carola/fig_carola_ocv_months.png
# out) resolve correctly. Any argument is passed through, e.g. --figdir /tmp/x.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/carola/fig_carola_ocv_months.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/carola/fig_carola_ocv_months.py "$@"
