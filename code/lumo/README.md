# `code/lumo/` — LUMO site package: analysis and plotting scripts

Six Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values of the control
tests; without it the scripts still run, just with fewer p-values). This is the
LUMO half of the repository; the Bautzen bridge site lives in `../bautzen/` and
mirrors it module for module — only the mask layer is swapped.

All paths below are relative to the **repository root** and the scripts are
always started from there (the wrappers `../../figures/lumo/fig_lumo_ocv_paper.sh`
and `../../figures/lumo/fig_lumo_ocv_orbits.sh` do that): inputs and outputs in
`data/lumo/`, figures in `figures/lumo/`.

### What the state labels mean

The dataset's own description of its reversible damage mechanisms is positional,
not ordinal: **DAM 1 … DAM 6 are six attachment positions on the tower**, and each
campaign removes the mechanisms at one or more of them —

* `DAM 3` = *all damage mechanisms removed at one strut* (the dataset's wording),
* `DAM 4` = the same pattern at a **different strut**,
* `DAM 6` = *all damage mechanisms removed* (the most extensive state).

The measurement levels are `ML 1 … ML 10` (ten accelerometer heights, ~1 m apart,
which is why the tower is modelled as 3 segments × 3 m).  Consequences that matter
for every number in this package:

* `DAM 3` and `DAM 4` are **different locations with the same mechanism pattern**,
  and they produce the *same* first-higher-mode drop (−1.00 Hz, −7.41 %) — two
  different damages that global modal frequency cannot tell apart.  Grouping them
  as one "dose level" (as `lumo_state_frequencies.json` does, `L1`) is therefore a
  grouping by *measured frequency*, never by damage type or location.
* The dataset publishes **six** damaged campaigns — the three positions above in
  the pattern `111` *and* in a second pattern `010` — of which this repository uses
  only the `111` sets (`02_DAM6_111`, `04_DAM4_111`, `06_DAM3_111`).  The `010`
  sets are the same three positions with a different set of mechanisms removed;
  they are not part of the committed tables, and their campaign dates have not been
  checked against the SAR record.

### `lumo_ocv_stats.py` — statistics layer (verbatim port)

`_stats`, `welch_mw_test`, `cliffs_delta`, `delta_sds`, `bootstrap_delta_ci`,
`monotonicity` (Jonckheere-Terpstra + Spearman), `strata`, `loo_cv`,
`permutation_test`, `oof_scores`, `permutation_p`, `paired_bootstrap_diff`,
`feature_dataset`, `mask_definition`, `by_state`, `vals`. Ported line by line
from the LUMO project's analysis scripts so that the numbers stay comparable with
the committed results — no reimplementation, no "improved" variant.

### `lumo_ocv_masks.py` — mask layer (verbatim port)

The echo-mask rule of `analyze_lumo_echo_mask_vector.py`: threshold at
`0.30 x peak` and `5 x median`, at least 2 pixels, then
`A = n_masked`, `D = A / bbox_area`, `F = 8-neighbourhood components`,
`S = ||mask centroid - peak||`, plus the persistence `P` (mean pixel frequency
across the 177 majority-shape masks).

**Audit note (2026-10).** `A`, `D`, `F`, `S`, `centroid_r/c` and `bbox_area`
reproduce **exactly** (deviation 0, 178 of 178 overpasses) when the rule is
replayed against `support/lumo_rects_ring.txt`, which carries every chip cell —
see `physical_model-observability/code/mask_replay.py`.  The two columns
`peak_r` / `peak_c` do **not**: they agree neither with the strip argmax used by
`S`, nor with the brightest cell of the mask, nor with the argmax of the 11 bin
columns, nor with the chip-wide argmax.  `S` itself reproduces, so the peak used
in the computation is the strip argmax; the two audit columns were evidently
written from something else.  Nothing in the LUMO pipeline reads them, but they
should not be quoted as "the peak pixel" until the producer is identified.

### `lumo_ocv_core.py` — data layer and constants

Paths, the channel/feature order (`FEATURES = [gamma2, P, D, A, F, S]`), the
three feature sets of Figure D (`MODELS`), the shrinkage/seed/permutation
constants, the label maps (`SHORT`, `SHORT_ALC`, `LABEL_EN`) and the loaders:
`load_csv()` (the extended overpass table), `load_csv_meta()`,
`load_reference()` (the five committed JSONs), `feature_dataset()`,
`n_by_state()`, `sha256()`.

### `lumo_ocv_channels_csv.py` — generator of the extended overpass table

Reads `data/lumo/lumo_channels.csv` (the copy of the committed overpass table)
and the machine-local burst cache four levels above this folder
(`<repo>/../../lumo_dam6_analysis/monthly_bursts/`, a pipeline artefact, ~100 MB,
not committed, overridable with `--burst-dir`), computes the mask geometry of the
178 coherence overpasses with `lumo_ocv_masks.py` and appends the columns
`A, D, F, S, P` (plus the auditable ingredients `mask_rows, mask_cols,
bbox_area, centroid_r/c, peak_r/c`) to `data/lumo/lumo_ocv_channels.csv`. It
verifies while writing: every overpass of the committed table must be found in
the cache, and `A` must equal the committed `n_masked`.

```bash
python3 code/lumo/lumo_ocv_channels_csv.py                # build + verify
python3 code/lumo/lumo_ocv_channels_csv.py --verify-only  # verify only
```

### `fig_lumo_ocv_paper.py` — figures A–E, JSON result, report, pin block

The main script. `compute()` recalculates every number of the five figures from
`data/lumo/lumo_ocv_channels.csv`; `pin_all()` then compares each number with the
committed reference JSONs in `data/lumo/reference/`, key by key and row by row
(2982 checks in a full run, 2892 in `--quick`; the German documentation strings of
the references are translated through `core.LABEL_EN`). A single deviation ends
the run with exit code 1. `main()` writes the five figures to `figures/lumo/`,
the result JSON to `data/lumo/fig_lumo_ocv_paper.json` and the English report to
`figures/lumo/fig_lumo_ocv_paper.md`.

```bash
python3 code/lumo/fig_lumo_ocv_paper.py            # full run, ~5 min
python3 code/lumo/fig_lumo_ocv_paper.py --quick    # ~2 min: the two permutation
                                                   # tests are copied from the refs
python3 code/lumo/fig_lumo_ocv_paper.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

Wall-clock times are only printed to stdout, never written into an artifact: with
the fixed seeds (`RNG_SEED = 11`, 2000 bootstrap draws, 1000 / 500 permutation
draws) two runs produce identical JSON, Markdown and PNG files.

### `fig_lumo_ocv_orbits.py` — 6 x 2 time-series companion

The paper figure summarises the six channels per state; this companion asks the
same channels as a **time series** and separates the two pass geometries. It
draws six rows in the formal OCV order (`gamma2`, `P`, `D`, `F`, `A`, `S`) and
two columns (ASC / DESC): twelve panels on `figures/lumo/fig_lumo_ocv_orbits.jpg`
(183 mm wide, 630 dpi JPEG). The x axis is the **actual local overpass instant**:
the coherence rows carry no timestamp, so `compute()` joins the
`record_class="acquisition"` rows by `burst_id`, converts UTC to
`Europe/Berlin` with `zoneinfo`, sorts by acquisition instant and cuts the
**13 consecutive `damage_label` runs** (36, 8, 6, 8, 54, 17, 7, 8, 4, 8, 5, 8,
9 points). Bands follow those runs, not fixed windows; twelve dashed lines mark
the transitions. `P`, `D`, `F` and `A` use log axes; the one missing `P`
(2021-07-14 ASC) is drawn as a gap and never interpolated.

The run writes `data/lumo/fig_lumo_ocv_orbits.json` and
`figures/lumo/fig_lumo_ocv_orbits.md` as well. It carries **144 cross-checks**
against `lumo_ocv_channels_meta.json`, `fig_lumo_ocv_paper.json` and
`reference/lumo_tower_coherence_states.json`: row and date counts, the 482-row
acquisition join, the 13-run sequence and lengths, sub-swath counts,
per-state/per-orbit counts, row-by-row `burst_id` / date / orbit / label /
`gamma2` / `n_masked` parity, `A == n_masked`, channel completeness, strict
positivity on the log axes, local-time ranges / UTC offsets and every committed
sha256. Any failed check exits non-zero. No output carries a timestamp, so
repeated runs are byte-identical.

```bash
python3 code/lumo/fig_lumo_ocv_orbits.py              # figure + JSON + report
python3 code/lumo/fig_lumo_ocv_orbits.py --no-figures # checks only
```

### `fig_lumo_ocv_physics.py` — physical dynamics against OCV behaviour

Two panels, one state colour scheme, asking whether states that sit close
together in an **independently measured physical** plane also sit close together
in the OCV behaviour plane:

| panel | x | y |
| --- | --- | --- |
| A physical | `df_H1 / f_H1,0` (signed, frequency loss to the left) | `dzeta_H1 / zeta_H1,0` |
| B OCV | `z(gamma2)` | `z(P)` |

Both origins are the healthy state. Panel A has **one point per state** — each
state was measured in one controlled campaign of five healthy plus five damaged
recordings, so there is no physical replicate; it carries a 95 % bootstrap
interval over those ten recordings (seeded, 2000 draws) and a ±0.125 Hz
quantisation bar on the frequency axis. Panel B has **two visits per state** (the
two SAR state runs); colour = state, marker shape = orbit (circle ASC, square
DESC), fill = visit (filled first, open repeat), marker area ∝ n, a grey ring
marks n < 3. The orbits are never pooled; the dashed arrow is the state-level
(visit-pooled) vector, which is the reproducible object.

The report (`figures/lumo/fig_lumo_ocv_physics.md`) gives the direction and
magnitude of every vector, the visit-to-visit repeatability (cosine and angle,
with the per-visit n), the state ordering by `|v|` in both planes and the
cross-plane cosine **as a descriptive statistic only** — a matching direction
would be evidence of consistency between two measurement channels, not evidence
that `(gamma2, P)` measures stiffness or damping.

Data: `data/lumo/reference/lumo_state_dynamics.json` (H1 frequency change from
the committed Welch mode consensus, H1 damping from the half-power bandwidth of
30 recordings; built and verified by `lumo_state_dynamics_json.py`) and the
**support** chip `data/lumo/support/lumo_ocv_channels.csv`. `--chip committed`
repeats the whole panel on `data/lumo/lumo_ocv_channels.csv`: the two chips
define `gamma2` differently (the support re-cut raises it by an order of
magnitude), and the report shows how far that moves the OCV points.

```bash
python3 code/lumo/fig_lumo_ocv_physics.py                 # support chip
python3 code/lumo/fig_lumo_ocv_physics.py --chip committed
python3 code/lumo/fig_lumo_ocv_physics.py --selftest       # geometry helpers
python3 code/lumo/fig_lumo_ocv_physics.py --no-figures     # checks only
```

### `lumo_state_dynamics_json.py` — builder of the damping reference

The damping half of `v_phys` lives outside this repository (LUMO SHM project).
This script reads the three sources, pins their sha256, copies the thirty
per-recording damping values verbatim into
`data/lumo/reference/lumo_state_dynamics.json` and cross-checks the medians
against `lumo_dam6_analysis/lumo_damping_gamma2_modulation.json`. Without
`--build` it only verifies the committed file (107 checks: medians, spreads,
relative changes, campaign names, agreement with the committed frequency
reference).

```bash
python3 code/lumo/lumo_state_dynamics_json.py --build   # rebuild (needs the SHM data)
python3 code/lumo/lumo_state_dynamics_json.py           # verify only
```

Note that the module and site names keep their `lumo_` prefix although the folder
already says `lumo/`: the import statements (`import lumo_ocv_core as core`),
the reference file names and every command line in the reports stay unchanged
between the two sites, which is what makes the mirror auditable as a diff.
