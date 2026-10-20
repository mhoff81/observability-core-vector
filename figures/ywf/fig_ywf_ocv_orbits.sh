#!/usr/bin/env bash
# Runs code/ywf/fig_ywf_ocv_orbits.py from the repository root, so its defaults
# (data/ywf/ywf_ocv_channels.csv, data/ywf/ywf_windows_full.txt,
# data/ywf/ywf_windows_mask_cache.txt and its manifest,
# data/ywf/ywf_measurements_full.txt, data/ywf/ywf_segments.txt,
# data/ywf/ywf_ocv_channels_meta.json, data/ywf/fig_ywf_ocv_paper.json and
# data/ywf/reference/ywf_ocv_findings.json + ywf_bursts.json + ywf_geometry.json
# + ywf_osm_road.json in; data/ywf/fig_ywf_ocv_orbits.json,
# figures/ywf/fig_ywf_ocv_orbits.md and figures/ywf/fig_ywf_ocv_orbits.png out)
# resolve correctly. Any argument is passed through, e.g. --no-figures.
#
# Interpreter: `python3` from PATH, unless PYTHON names another one — needed
# where the default python3 is a virtualenv without numpy/matplotlib, e.g.
#   PYTHON=/usr/bin/python3 bash figures/ywf/fig_ywf_ocv_orbits.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${PYTHON:-python3}" code/ywf/fig_ywf_ocv_orbits.py "$@"
