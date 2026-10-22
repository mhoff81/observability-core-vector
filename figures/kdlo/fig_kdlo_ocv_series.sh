#!/usr/bin/env bash
# Runs code/kdlo/fig_kdlo_ocv_series.py from the repository root, so its defaults
# (data/kdlo/kdlo_ocv_channels.csv, data/kdlo/fig_kdlo_0_raw.json and
# data/kdlo/reference/kdlo_s1_availability.json in; data/kdlo/fig_kdlo_ocv_series.json,
# figures/kdlo/fig_kdlo_ocv_series.md and figures/kdlo/fig_kdlo_ocv_series.jpg
# out) resolve correctly. Any argument is passed through, e.g. --figdir /tmp/x.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/kdlo/fig_kdlo_ocv_series.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/kdlo/fig_kdlo_ocv_series.py "$@"
