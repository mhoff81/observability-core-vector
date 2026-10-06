#!/usr/bin/env bash
# Runs code/bautzen/fig_bautzen_ocv_amp_phase.py from the repository root, so its
# defaults (data/bautzen/bautzen_ocv_channels.csv and the four complex rect files
# in; data/bautzen/fig_bautzen_ocv_amp_phase.json,
# figures/bautzen/fig_bautzen_ocv_amp_phase.md and
# figures/bautzen/fig_bautzen_ocv_amp_phase.png out) resolve correctly. Any
# argument is passed through, e.g. --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/bautzen/fig_bautzen_ocv_amp_phase.py "$@"
