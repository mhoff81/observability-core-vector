#!/usr/bin/env bash
# Runs code/haltern/fig_haltern_ocv_series.py from the repository root, so its
# defaults (data/haltern/haltern_ocv_channels.csv, data/haltern/haltern_windows_full.txt,
# data/haltern/haltern_windows_mask_cache.txt and its manifest,
# data/haltern/haltern_ocv_channels_meta.json, data/haltern/fig_haltern_ocv_paper.json and
# data/haltern/reference/haltern_ocv_findings.json in; data/haltern/fig_haltern_ocv_series.json,
# figures/haltern/fig_haltern_ocv_series.md and figures/haltern/fig_haltern_ocv_series.png out)
# resolve correctly. Any argument is passed through, e.g. --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/haltern/fig_haltern_ocv_series.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/haltern/fig_haltern_ocv_series.py "$@"
