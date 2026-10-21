#!/usr/bin/env bash
# OCV time series per dimension and orbit, for every support-sized table.
# Run from the repository root: bash figures/orbits/fig_ocv_orbits.sh
#
# The script needs numpy + matplotlib. Whoever runs it may not be the one who
# installed them, so the interpreter is chosen by trying candidates in order
# (override with PYTHON=...). The chosen one is printed, so a report can say
# which interpreter drew it.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=""
for candidate in "${PYTHON:-}" python3 /home/projects/erlking_mvp/.venv/bin/python3; do
  if [ -n "$candidate" ] && "$candidate" -c 'import matplotlib, numpy' >/dev/null 2>&1; then
    PY="$candidate"
    break
  fi
done
if [ -z "$PY" ]; then
  echo "no interpreter with matplotlib + numpy found; set PYTHON=" >&2
  exit 2
fi
echo "interpreter: $PY"
"$PY" code/orbits/fig_ocv_orbits.py --verify "$@"
