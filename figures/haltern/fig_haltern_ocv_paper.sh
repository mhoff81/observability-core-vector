#!/usr/bin/env bash
# Runs code/haltern/fig_haltern_ocv_paper.py from the repository root, so its
# defaults (data/haltern/haltern_ocv_channels.csv and
# data/haltern/reference/haltern_ocv_findings.json in; data/haltern/fig_haltern_ocv_paper.json,
# figures/haltern/fig_haltern_ocv_paper.md and figures/haltern/fig_haltern_ocv_paper_[A-E].png out)
# resolve correctly. Any argument is passed through, e.g. --write-reference.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/haltern/fig_haltern_ocv_paper.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/haltern/fig_haltern_ocv_paper.py "$@"
