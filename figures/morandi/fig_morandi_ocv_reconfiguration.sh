#!/usr/bin/env bash
# Runs code/morandi/fig_morandi_ocv_reconfiguration.py from the repository root,
# so its defaults (data/morandi/support/morandi_ocv_channels.csv in;
# data/morandi/ocv_reconfiguration_morandi.json,
# data/morandi/ocv_reconfiguration_morandi_table.csv,
# figures/morandi/fig_morandi_ocv_reconfiguration.jpg and .md out) resolve
# correctly. Any argument is passed through, e.g. --figdir /tmp/x or --selftest.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/morandi/fig_morandi_ocv_reconfiguration.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/morandi/fig_morandi_ocv_reconfiguration.py "$@"
