#!/usr/bin/env bash
# Runs code/cases/ocv_prospective.py from the repository root, so its defaults
# (data/cases/eas_policy_pin.json in; data/morandi/fig_morandi_ocv_trajectory_prospective.json,
# figures/morandi/fig_morandi_ocv_trajectory_prospective.md and .jpg, data/cases/fig_ocv_prospective_comparison.json and
# figures/cases/fig_ocv_prospective_comparison.md out) resolve correctly. Any argument is passed
# through, e.g. --verify or --selftest.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one -- needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/morandi/fig_morandi_ocv_trajectory_prospective.md.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/cases/ocv_prospective.py "$@"
