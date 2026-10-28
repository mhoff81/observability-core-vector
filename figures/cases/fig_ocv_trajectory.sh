#!/usr/bin/env bash
# The six OCV trajectory plates in one run — one figure per site, nothing pooled:
#   figures/<site>/fig_<site>_ocv_trajectory.jpg   the 3 x 4 plate
#   figures/<site>/fig_<site>_ocv_trajectory.md    its report
#   data/<site>/fig_<site>_ocv_trajectory.json     every plotted number
# Any argument is passed through, e.g. --verify, --no-figures or --quiet.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/cases/fig_ocv_trajectory.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/cases/ocv_trajectory.py --all "$@"
