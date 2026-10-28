#!/usr/bin/env bash
# Runs code/werder_havel/fig_werder_havel_ocv_trajectory.py from the repository root, so
# its defaults (data/werder_havel/support/werder_havel_ocv_channels.csv,
# data/werder_havel/support/werder_havel_ocv_channels_meta.json and
# data/werder_havel/reference/werder_havel_aoi.json in;
# data/werder_havel/fig_werder_havel_ocv_trajectory.json,
# figures/werder_havel/fig_werder_havel_ocv_trajectory.md and
# figures/werder_havel/fig_werder_havel_ocv_trajectory.jpg out) resolve correctly.
# The documented works interval is read from the committed channel metadata; the
# external operator record (EVENTS_SOURCE) may only agree with it or be reported as a
# note -- it never supplies a date. Any argument is passed through, e.g. --selftest
# (the registered rules and pins, nothing written), --verify (re-derive and compare with
# the committed JSON, report and plate) or --no-figures (the JSON and the .md only).
#
# Interpreter: `python3` from PATH, unless PYTHON names another one -- needed where
# the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/werder_havel/fig_werder_havel_ocv_trajectory.py "$@"
