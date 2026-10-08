#!/usr/bin/env bash
# Runs code/ywf/fig_ywf_ocv_paper.py from the repository root, so its defaults
# (data/ywf/ywf_ocv_channels.csv and data/ywf/reference/ywf_ocv_findings.json in;
# data/ywf/fig_ywf_ocv_paper.json, figures/ywf/fig_ywf_ocv_paper.md and
# figures/ywf/fig_ywf_ocv_paper_[A-E].png out) resolve correctly. Any argument is
# passed through, e.g. --write-reference.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/ywf/fig_ywf_ocv_paper.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/ywf/fig_ywf_ocv_paper.py "$@"
