# observability-core-vector

Two standalone, offline figure/report packages for two Sentinel-1 observability
sites, side by side in one repository. Each answers the same question about its
own **observability core vector**

```
OCV = [gamma2, P, D] = [interferometric coherence, mask persistence, mask density]
```

for its own structural states — and neither claims an answer: both **recompute**
the committed numbers of the site's analysis project from a small input table,
render their figures and report, and **pin every number against the committed
reference JSONs**: a single deviation makes the run exit with a non-zero status.
The claim of the repository is therefore not "we assert", but "we reproduce and
read off".

| site | folders | states | mask layer | results |
| --- | --- | --- | --- | --- |
| LUMO lattice tower (Hannover) | `code/lumo/`, `data/lumo/`, `figures/lumo/` | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`) | [`figures/lumo/fig_lumo_ocv_paper.md`](figures/lumo/fig_lumo_ocv_paper.md) |
| Bautzen OpenLabs bridge | `code/bautzen/`, `data/bautzen/`, `figures/bautzen/` | RS, DS1, DS2 | deck-edge mask (geometrically anchored deck band) | [`figures/bautzen/`](figures/bautzen/) (figure script pending) |

The Bautzen package is a **strict mirror** of the LUMO one — same module names,
same file layout, same JSON/Markdown shapes, same determinism — and only the
mask layer (and with it the site constants) differ, so the two packages can be
read as a diff of one another.

## What the LUMO figures show

| figure | question | content |
| --- | --- | --- |
| A | what do the channels look like? | raw distributions of gamma2, A, D, F, S, P over the 178 coherence overpasses, per state |
| B | which component carries a difference? | in-sample effect sizes: Cliff's delta, bootstrap CI, SDS = 10&#124;delta&#124;, for all six components |
| C | is the state recoverable at all? | 4-class LDA, leave-one-out, shrinkage grid, confusion matrix, label-permutation null |
| D | is the separation state-specific? | out-of-fold scores of the six binary state pairs for four feature sets (+ paired bootstrap vs. gamma2 alone) |
| E | is it the state driving it? | severity trend (Jonckheere-Terpstra, Spearman), season/orbit strata, mask brightness and wind covariates |

`P` belongs to the core vector; `A`, `F` and `S` are controls for the
*morphology* of the same echo mask (they are all derived from the identical mask
as `P` and `D`).

## What the Bautzen figure shows

The same five panels for the bridge: the six channels per series (ASC/DESC
bridge and control 500 m east), effect sizes, state recovery, pairwise
separation and the controls — but on the **deck-edge mask** instead of the echo
mask, and with the site's own two documented caveats:

* only the **ASC bridge** series passes the deck-edge detection gate (a global
  echo rule of the LUMO kind is not usable at this site: the scene clutter pulls
  the reference median up and the mask collapses to ~0 pixels);
* **state and season are perfectly confounded** at this site (RS is spring, DS1
  summer, DS2 autumn), so the environment covariates are a counter-check, not a
  control.

Read [`figures/bautzen/fig_bautzen_ocv_paper.md`](figures/bautzen/fig_bautzen_ocv_paper.md)
for the verdict; it is written by the script and contains no claim that the data
does not support.

## Folder structure

```
observability-core-vector/
├── code/
│   ├── lumo/      LUMO analysis + plotting scripts — see code/lumo/README.md
│   ├── bautzen/   the Bautzen mirror (deck-edge mask) — see code/bautzen/README.md
│   └── README.md
├── data/
│   ├── lumo/      the copied LUMO overpass table, the extended table generated
│   │              from it, and the five committed reference JSONs
│   ├── bautzen/   the four rect files, the generated channel table, the two
│   │              committed reference JSONs and the hourly weather series
│   └── README.md
├── figures/
│   ├── lumo/      the .sh wrapper, the .md report and the .png figures
│   ├── bautzen/   the same for the bridge site
│   └── README.md
├── .gitignore
└── README.md
```

Per site: `code/<site>/` reads `data/<site>/…`, writes `data/<site>/*.json` and
`figures/<site>/{*.md,*.png}`. See the three READMEs for the file-by-file lists.

## Quick start

```bash
cd research/observability-core-vector

# LUMO: figures + report + JSON, verified against the five committed references
bash figures/lumo/fig_lumo_ocv_paper.sh            # ~5 min (two permutation tests)
bash figures/lumo/fig_lumo_ocv_paper.sh --quick    # ~2 min (permutation tests copied)

# Bautzen: the same for the bridge — **script pending**, see code/bautzen/README.md
# bash figures/bautzen/fig_bautzen_ocv_paper.sh

# if `python3` of the PATH is a virtualenv without numpy/matplotlib:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_paper.sh

# rebuild the channel tables (LUMO needs the local burst cache, Bautzen only
# the committed rect files)
python3 code/lumo/lumo_ocv_channels_csv.py
python3 code/bautzen/bautzen_ocv_channels_csv.py
```

Everything is offline and deterministic: fixed seeds (RNG seed 11 and 2000
bootstrap draws for LUMO, 7 / 42 with 10000 bootstrap draws for Bautzen, 20000
permutation draws in the port), no timestamps in any output, so repeated runs
produce identical JSON, Markdown and PNG/JPEG files. Requirements: Python 3 with
numpy, matplotlib and scipy (scipy only for the p-values of the control tests).
