#!/usr/bin/env bash
# Runs code/lumo/fig_lumo_ocv_physics.py from the repository root, so its
# defaults (data/lumo/support/lumo_ocv_channels.csv and
# data/lumo/reference/lumo_state_dynamics.json in;
# data/lumo/fig_lumo_ocv_physics.json,
# figures/lumo/fig_lumo_ocv_physics.md and fig_lumo_ocv_physics.jpg out) resolve
# correctly. Any argument is passed through, e.g. --chip committed, --selftest
# or --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_physics.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/lumo/fig_lumo_ocv_physics.py "$@"
