#!/usr/bin/env bash
# Runs code/cts/fig_cts_ocv_paper.py from the repository root, so its defaults
# (data/cts/cts_ocv_channels.csv in; data/cts/fig_cts_ocv_paper.json,
# figures/cts/fig_cts_ocv_paper.md and figures/cts/fig_cts_ocv_paper_[A-F].png
# out) resolve correctly. Any argument is passed through, e.g. --quick.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy / matplotlib / sklearn,
# e.g.
#   PYTHON=/usr/bin/python3 bash figures/cts/fig_cts_ocv_paper.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/cts/fig_cts_ocv_paper.py "$@"
