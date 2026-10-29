# `code/` — the eight site packages

One subfolder per site, all standalone, offline and built the same way (stats
layer, data layer, generator, figure script). They differ in what the observability
core vector contains: LUMO masks the **echo** of the tower, Bautzen masks the
**deck edge** of the bridge, KDLO consumes the mask size `A` that the
committed record delivers (it has no mask layer of its own), and Carola, Morandi
and CTS each mask both the **echo of every chip** and the **per-girder /
per-track deck** of a collapsed structure. CTS is the one
**difference-in-differences** site: four footprints of the *same* orbit (target
plus two control towers) share strata, dates and counts, so its package carries a
first-class DiD layer (`cts_ocv_did.py`) next to the mask layer. Espoo is the
one site with **no event at
all**: 150 acquisitions with an empty `damage_label`, `structural_state` and
`condition_label` on every row and **no mask raster** either (56 of its 78 columns
are empty), so its OCV is `[gamma2, A]` with `A` = the delivered mask *size*,
**copied and proved** rather than recomputed, and the only contrast left is the
orbit geometry (ASC / DESC) — the family's event-free control.

| folder | site | states | mask layer | package README |
| --- | --- | --- | --- | --- |
| [`lumo/`](lumo/README.md) | LUMO lattice tower, Hannover | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`, ≥ 2 px) | `lumo/README.md` |
| [`bautzen/`](bautzen/README.md) | Bautzen OpenLabs bridge | RS, DS1, DS2 | deck-edge mask (geodetically anchored deck band) | `bautzen/README.md` |
| [`kdlo/`](kdlo/README.md) | KDLO media tower, Garden City SD | pre-collapse, during-rebuild (epochs) | echo mask recomputed from a committed 64-strip full-dwell cache (LUMO rule, upper median) | `kdlo/README.md` |
| [`carola/`](carola/README.md) | Carolabrücke, Dresden (collapsed 2024-09-11) | pre-collapse (healthy), post-collapse | two layers: the echo mask of every chip **and** the per-girder deck mask, both from a committed window-payload cache (`0.30 x peak`, `5 x np.median`, ≥ 2 px) | `carola/README.md` |
| [`morandi/`](morandi/README.md) | Ponte Morandi / Polcevera, Genova (collapsed 2018-08-14) | pre-collapse (healthy), post-collapse | two layers: the echo mask of every chip **and** the per-track deck mask, both from a committed window-payload cache (`0.30 x peak`, `5 x np.median`, ≥ 2 px) | `morandi/README.md` |
| [`cts/`](cts/README.md) | Champlain Towers South, Surfside (collapsed 2021-06-24) | pre-collapse (healthy), post-collapse | two layers over four footprints of one orbit (target + two controls + a beach reference), the echo mask **and** the per-track deck mask from a committed window cache (`0.30 x peak`, `5 x np.median`, ≥ 2 px), **plus a difference-in-differences layer** | `cts/README.md` |
| [`haltern/`](haltern/) | Haltern AV9 wind turbine, Kreis Recklinghausen (collapsed 2021-09-29) | pre-collapse (healthy), post-collapse | **one layer**: the echo mask of the committed 7 x 7 window payload (`0.30 x peak`, `5 x np.median`, ≥ 2 px) — the record is not segment-resolved, so a per-mast-section layer would be dishonest, and no edge/band layer is claimed because the tower base spans more than the frame reaches; two further references are *resolved* from public catalogues (the CDSE burst pass of every chip, the turbine's **own OSM node**) and committed beside it | `haltern_ocv_core.py` docstring (the package carries no separate README) |
| [`espoo/`](espoo/README.md) | Espoo Kurttila mast, Espoo, Finland (**event-free** 2024-09-05 … 2026-09-07) | none — the record has no damage label at all; the contrast is the orbit geometry (ASCENDING n = 70 / DESCENDING n = 80) | **no mask layer**: the record delivers only the mask *size* `A` (`coherence_masked_pixels`) and no raster (56 of 78 columns empty), so the OCV is `[gamma2, A]` and `A` is **copied and proved**, never recomputed | `espoo/README.md` |

### `orbits/` — one OCV time-series figure per support-sized table

`[orbits/fig_ocv_orbits.py](orbits/fig_ocv_orbits.py)` draws the twelve-panel
companion (six where only one orbit has data — KDLO is DESC-only, so its plate
has one column and the report names the missing channel) for the six sites
migrated to the support-sized chip (CTS, KDLO, Carolabrücke, Haltern, Ponte
Morandi, and Espoo): one panel per **dimension
x observation channel**, in the order and with the axes the LUMO and Bautzen
figures use (`gamma2, P, D, F, A, S`; `P, D, F, A` logarithmic).

Two states for the five labelled sites: **pre-collapse and post-collapse**,
except KDLO, whose record carries **pre-rebuild and post-rebuild**. A state with a
single overpass (KDLO's `single`) is dropped and named in the markdown. A channel
with no data is not drawn empty — the panel carries the reason (KDLO's support
cut holds 60 overpasses, all DESCENDING, so its ASC column has nothing to draw).

**One state is allowed too**, for a site whose record labels no damage at all:
Espoo registers `states = {"all": "event-free"}` and `label_column = None`, so
every overpass takes the single state, the band is the one run, the median is the
level of the whole series, and the title reads *across* instead of *before and
after*. The two-state sites keep every number: their JSON gains exactly one key,
`out_name`, and no value, run or median changes. Everything state-dependent is now driven by the
registry: the per-state loops, the counts column, the legend and the markdown
tables all iterate `spec["states"]`, and `STATE_COLORS/FILL/MARKERS` carry an
`all` slot (a neutral fill, because one state has no contrast to colour).

Espoo also carries `out_name` (`fig_espoo_ocv_orbits_chip`), because the site
already owns `fig_espoo_ocv_orbits.jpg` — the *pipeline* figure (A vs S, no chip,
see [`espoo/README.md`](espoo/README.md)). Same base name, different figure: the
`.png` is the pipeline record, the `.jpg`/`.json`/`.md` of the `_chip` stem are
the shared six-dimension family on the support chip.

```sh
bash figures/orbits/fig_ocv_orbits.sh            # all six, with --verify
bash figures/orbits/fig_ocv_orbits.sh --site cts
```

Outputs: `figures/<site>/fig_<site>_ocv_orbits.{jpg,md}` and
`data/<site>/support/fig_<site>_ocv_orbits.json`, or the registry's `out_name`
where a site already owns that stem (Espoo's `fig_espoo_ocv_orbits_chip.*`).
`--verify` re-derives every plotted point, run and median from the table and
fails on a mismatch — the figure computes nothing the JSON does not carry.

### `cases/` — the six-case overview page

`[cases/fig_cases_overview.py](cases/fig_cases_overview.py)` reads the six cases on
the **support-sized chip** (`data/sds/sds_support.json` for the state contrasts,
Espoo's chip figure for the event-free control) and draws one page: the case map on the event spectrum, the
dimension-response matrix, and six schematic time lines. A matrix cell carries an
arrow only when that case's own interval excludes zero **and** both sides carry at
least six units (`○` = tested, nothing resolved; `?` = no evidence; fill = the SDS
band inside the cell), so Haltern's single resolved dimension is `F ↑` and
Morandi's row is four circles on purpose. `--selftest` checks the gate rules and
the three adapters; `--verify` re-derives every cell and every source hash.

```sh
python3 code/cases/fig_cases_overview.py            # JSON + report + 630 dpi JPEG
python3 code/cases/fig_cases_overview.py --selftest # 11 checks
python3 code/cases/fig_cases_overview.py --verify   # every cell, every sha256
```

### `cases/` — the standard OCV trajectory plate

`[cases/ocv_trajectory.py](cases/ocv_trajectory.py)` draws the **same 3 x 4 plate**
for every case site, on the support-sized chip, with the post-event tail **drawn**
instead of cut off — the timeline question the six cases are read for, asked once
and answered six times in one design. The reference timeline is defined **case by
case** (R2), each of the two blocks of `ocv_reconfigure.BLOCKS` is standardised on
**its own** reference set (R3), the z is robust and clipped (R4), and the months are
bins with holes drawn, never spanned, with one seeded moving-block interval for the
contrast (R5).

| site | reference window | chips | tail | contrast |
| --- | --- | --- | --- | --- |
| carola | 2022-09-03 → 2024-09-09 | 184 | 252 post-collapse | resolves `gamma2, P, S` |
| cts | 2015-09-21 → 2021-06-21 | 136 | 42 post-collapse | resolves `P, A, D` |
| morandi | 2016-08-24 → 2018-08-08 | 198 | 9 post-collapse | resolves `D, R_G` |
| kdlo | 2022-06-02 → 2022-12-11 | 16 | 43 rebuild | resolves `P, A, S, R_B, R_G` |
| haltern | 2021-01-01 → 2021-09-17 | 96 | 11 post-collapse | resolves nothing |
| espoo | 2024-09-05 → 2026-09-07 | 75 | none | no contrast (event-free) |

Two states are the rule; KDLO carries three and its `single` acquisition (the
collapse day itself) is excluded before any statistic, and Espoo carries one and
standardises **per orbit**, because the pass geometry is the only contrast its
record has. The full rule set, the rejected windows and the panel inventory are in
[`ocv_trajectory_spec.md`](ocv_trajectory_spec.md).

```sh
python3 code/cases/ocv_trajectory.py --site carola   # one plate (+ JSON + report)
python3 code/cases/ocv_trajectory.py --all           # the six
python3 code/cases/ocv_trajectory.py --selftest      # 201 checks, writes nothing
python3 code/cases/ocv_trajectory.py --site kdlo --verify
```

### `cases/` — the correspondence model (campaign as the unit)

`[cases/ocv_correspondence.py](cases/ocv_correspondence.py)` tests the question the
two-panel figure only *shows*: do states that sit close in the physical plane also
sit close in the OCV plane? One rule is frozen in the code and in every report:
**the campaign / visit is the independent unit, never the single overpass** —
several SAR acquisitions inside one structural campaign improve the estimate of
`(gamma2, P)` but do not create additional independent measurements of `(f,
zeta)`. Inference intervals resample **campaigns**; the date-level bootstrap
appears only as `ci_measurement` inside a campaign and never increases n.

Four tests in the order the question forces: **A** repeatability of the OCV
response across visits (with the pre-declared rule `R = d_visit / d_between < 1`),
**B** a state-dependent ASC-DESC contrast (permutation of the six labels), **C**
the Spearman correlation between `d_phys` and `d_OCV` over the campaign pairs
(Mantel permutation of the campaign labels, exact at 720 relabellings; the
z-scored variant is a sensitivity, not a way to align the planes) and **D** the
pre-declared contrast between the two physically closest states. Two findings are
built into the code instead of being glossed over: `x` takes one value per state,
so **M2 (`orbit x x`) is saturated and is reported, not fitted**, and the residual
contrast of D is **zero by construction** for exactly that reason, which is why D
reports the *raw* contrast against the visit noise.

`--asset lumo` runs the six LUMO campaigns (17 pins: 178 overpasses, 13 runs, the
per-campaign/per-orbit date counts, the healthy reference counts, the exact 720
relabellings, six LOCO folds, the DAM3-DAM4 contrast pair); `--asset bautzen`
runs the degenerate arm (16 pins) where the same code reports A and C as *not
defined* and B as *not informative*, and where D degenerates to the single
RS→DS1 contrast that the asset can carry — Bautzen has one campaign per state and
no independent f/zeta for DS2. `--build-bautzen-reference` rebuilds
`data/bautzen/reference/bautzen_state_dynamics.json` from the OpenBridgeLAB
analysis; `--selftest` runs the synthetic positive, negative and degenerate
controls against the real statistics.

```sh
python3 code/cases/ocv_correspondence.py --asset all   # both JSONs + reports + JPEG
python3 code/cases/ocv_correspondence.py --selftest    # 15 checks incl. controls
python3 code/cases/ocv_correspondence.py --build-bautzen-reference
```

### `sds/` — the cross-site SDS report (the one script that reads two sites)

`[sds/sds_support.py](sds/sds_support.py)` asks the question the new chip has to
answer: **do the OCV dimensions still separate the labelled states, and for which
state pair?** It is a cross-site analysis, not a site package — it reads
`data/lumo/` and `data/bautzen/` (committed table and support-sized table each)
and writes `data/sds/`.

The metric is the project's own: `cliffs_delta` and `bootstrap_delta_ci` are
imported from [`lumo/lumo_ocv_stats.py`](lumo/lumo_ocv_stats.py), so
`SDS = 10 * |delta|` and the percentile bootstrap (10 000 resamples, seed 7) are
the definitions the committed LUMO reference was built with. Before it reports
anything it **pins** `data/lumo/reference/lumo_echo_mask_vector.json`: the five
`sds_pooled_all_damaged` values must come out identically or the run exits 1.

| level | role |
| --- | --- |
| `per_channel_per_state` | **primary analysis** — one SDS per observation channel (ASC / DESC), dimension and state pair. ASC and DESC are never pooled into a primary number |
| `per_channel_max_pairwise` | the compact summary; every entry carries its `winning_pair`, that pair's signed delta, **that pair's** interval and `n1 / n2` |
| `per_channel_pooled_damage` | healthy versus all damage states, one channel at a time |
| `pooled_all_damage` | both channels pooled — the **compatibility pin and context**, never the primary evidence |

Two rules are enforced in code rather than trusted to a reader: the writer
refuses to emit the phrase "CI of max" (an interval belongs to one comparison,
not to a selection over comparisons), and a maximum without its winning pair,
delta, interval or counts raises instead of leaving a blank cell.

```sh
python3 code/sds/sds_support.py            # both sites, 10 000 resamples
python3 code/sds/sds_support.py --selftest # 25 offline checks
```

### `sds/ocv_reconfigure.py` — do the two experiments reconfigure the same way?

The second cross-site script, and the one that asks the research question of
Paper 2: **does the magnitude of OCV reconfiguration grow with an independently
defined state dose, and is the direction the same in both experiments?** Every
number is a distance from that experiment's *own* reference state — LUMO against
`healthy`, Bautzen against `RS` on the paired bridge-minus-control difference, the two
bridges and the two towers against their own pre-collapse state — so the two assets
never share a scale and their magnitudes are never compared. Three registered pairs run
through it: the dose pair (`lumo,bautzen`, the default), the bridge pair and the tower
pair.

`code/sds/ocv_reconfigure_spec.md` fixes the design and was written **before** the
first run: the dose (LUMO's measured first-higher-mode drop, 0 / 7.41 / 12.96 %,
which makes `DAM 3` and `DAM 4` one level because they sit one FFT bin apart), the
transforms (`log` for positive dimensions, identity for `S` and Bautzen's signed
γ²), the zero rule (drop the unit for that dimension, count it), the blocks
(`behaviour = [γ², P]` (cross-acquisition echo behaviour) and `geometry = [A, D, F, S]` (spatial mask geometry, with `D` kept off the plane because it partly inherits `A` and `S` annotated because a shift is position, not morphology);
both blocks carry their dimensions into the magnitude), the tests and
the three pre-registered predictions.

| step | what it does |
| --- | --- |
| **H1, per asset** | Jonckheere–Terpstra on the **per-unit** magnitudes (the reference's median magnitude is 0 by construction, so a test on medians would be degenerate), **Holm within each asset's own family — two cells per channel the record carries**, moving-block bootstrap with the reference re-estimated inside every draw |
| **verdict, per asset** | each asset gets its own verdict **and a bound** — the largest top-step change its interval still covers — because neither alone is a null result, only a bound |
| **comparison** | cosine of the standardised response vectors at the largest state, per block and channel, with a random-direction null and a sign-agreement count. **Descriptive**: it asks whether the same metric applied twice answers the same way, not whether a response transfers. γ² is also reported without γ², because it is a different quantity in the two experiments |

The 2026-10-22 amendment in the spec separates the two assets. Two consequences a
reader should know: the reference level carries **no** interval (`R = 0` by
construction, its noise belongs to `ΔR`), and where the reference state has fewer
than ten units — Bautzen has 6 (ASC) and 7 (DESC) complete pairs — the **pooled
MAD** is the standardising scale, because a MAD from six numbers is itself too
unstable. `scale_used` travels with every cell, in the JSON and in the table.

```sh
python3 code/sds/ocv_reconfigure.py                        # the dose pair, ~40 s
python3 code/sds/ocv_reconfigure.py --sites carola,morandi  # the bridge pair
python3 code/sds/ocv_reconfigure.py --sites haltern,kdlo    # the tower pair
python3 code/sds/ocv_reconfigure.py --sites kdlo            # one asset, one channel, two states
python3 code/sds/ocv_reconfigure.py --selftest              # 27 checks, the last ten offline
python3 code/sds/ocv_reconfigure.py --verify                # re-derives every table cell
```

`--sites` chooses the assets. **Carola** (`--sites carola`) is the two-state case:
pre-collapse against post-collapse on the support-sized chip, where only `Girder A`
survives the cut (436 of the committed 1 717 chips). With two states the ordered
test is its two-group case (a one-sided rank test), there is **no ordering claim**,
and the H2 comparison is skipped because a single asset has no partner — the report
says all three things rather than leaving blanks. Its post signature is a time
contrast without a control track, which the closing section states.

**Carolabrücke and Ponte Morandi** (`--sites carola,morandi`) are the second
registered pair (spec amendment 4): two bridges read across a collapse instead of a
dose, on support-sized chips produced by the same window rule (207 chips of
Morandi's committed 232, 436 of Carola's 1 717). There is no dose, so the report
opens with `## There is no dose` and a per-state `n` table; H2 keeps `largest` only
(no `middle`, no 2 × 2 matrix, both declared under `not_computed`); and the report
closes with **Two verdicts, read together** — Carola resolves a large coherent DESC
reconfiguration (both blocks at Holm < 0.001), Morandi resolves nothing, and exactly
one of the four cosine cells shares a direction (DESC geometry, 4/4 signs, null
`p = 0.083` — a guard, not a test).

The **thin-cell rule** applies to every plate and report: a cell whose `n` does not
exceed the bootstrap block (3) is marked `✱` wherever it is quoted and named in
`bootstrap_caveats` in the JSON — its interval then carries the reference's
uncertainty only. Morandi's ASC post-collapse cell (n = 3) is the one such cell in
the bridge pair; the dose pair's list is empty.

**Haltern and KDLO** (`--sites haltern,kdlo`) are the third registered pair (spec
amendment 5): two towers read across an event that removed the structure — Haltern's
AV9 turbine came down on 2021-09-29, KDLO's mast collapsed between 2022-12-11 and
2022-12-23 **and was rebuilt**. That difference is what the pair is for, and it is
also what stops the run being a copy of the bridge pair: KDLO's second state is a new
structure, not a damaged one. Three consequences are handled by data rather than by
prose. First, **channels come from the record**: KDLO was acquired DESCENDING only, so
its cells are DESC cells, its Holm family has **two** members instead of four, its
`ASC` channel is declared under `channels_absent` / `not_computed`, and the plate draws
two empty ASC panels that say so — an empty cell would have read as a null. Second,
**the level names come from the registry**: KDLO's ladder is `pre → rebuild`, so the
verdict says "pre-to-rebuild" and no sentence hardcodes `post`. Third, **one
acquisition is excluded and counted**: KDLO's `single` state is the acquisition at the
collapse itself (2022-12-23), one date that cannot carry a median; it is removed before
any statistic and reported in the audit (`{"single": 1}`), never folded into a
neighbour.

Outputs `data/sds/ocv_reconfigure.{json,md}`,
`data/sds/ocv_reconfigure_table.csv` and three 630 dpi plates with one panel guide:
`fig_sds_ocv_reconfigure.jpg` (the juxtaposition), `..._lumo.jpg` and
`..._bautzen.jpg` (each asset standing alone). Each registered pair writes its own
stem — `data/sds/ocv_carola_morandi.{json,md}`, `..._table.csv` and the plates
`fig_sds_ocv_reconfigure_carola_morandi.jpg`, `..._carola.jpg`, `..._morandi.jpg`;
`data/sds/ocv_haltern_kdlo.{json,md}`, `..._table.csv` and
`fig_sds_ocv_reconfigure_haltern_kdlo.jpg`, `..._haltern.jpg`, `..._kdlo.jpg` — and
never touches the dose files.


The scripts are started from the **repository root** (each wrapper in
`../figures/<site>/` does that), so that their defaults resolve to
`data/<site>/…` in and `data/<site>/…`, `figures/<site>/…` out. No package:
the script's own folder is on `sys.path` (a script run puts its directory first),
which is why the modules import each other by their `lumo_…` / `bautzen_…` /
`kdlo_…` / `carola_…` / `morandi_…` / `cts_…` / `espoo_…` names.

Adding another site means copying one package and swapping the mask layer (if the
record has one) plus the `DATA`/`FIGDIR` lines of `…_ocv_core.py` — nothing
outside the new `code/<site>/`, `data/<site>/` and `figures/<site>/` changes. The
KDLO package is the example of a site *without* a mask layer: the two OCV
components are `gamma2` and the delivered mask size. The Espoo package is the
example of a site *without an event*: the same two components (`gamma2` and the
delivered `coherence_masked_pixels`), no damage label on any row, no mask raster
behind the size, and the orbit geometry (ASC / DESC) as the only contrast the
record supports — so its `espoo_ocv_channels_csv.py` **proves** the copied `A`
against the committed column and the committed reference statistics instead of
recomputing a mask, and its reports read every effect size against the record's own
no-event controls (stationarity, season, weather, coverage). The Carola package is the
example of a site *with two* mask layers over one chip: the echo mask of the
whole 80 x 80 window and the same rule on the girder's own window. The Morandi
package is the same shape one more time, with the two Sentinel-1 tracks
(`A_asc` / `A_des`) in place of the girders; because only 11 of its 232 chips
survive the collapse, its figures are an honest null rather than a contrast. The
CTS package keeps the same shape but adds what no single-structure site can: the
four footprints of its one orbit (the collapsed tower as `target`, two standing
towers as `control_ctn` / `control_cte`, the beach as `reference_beach`) let it
difference the target against a control that shares the atmosphere, so
`code/cts/cts_ocv_did.py` pins the three committed
`data/cts/reference/cts_reference_did_*.json` files (series and markdown, 1e-9).
Its figure script `code/cts/fig_cts_ocv_paper.py` (run by
`figures/cts/fig_cts_ocv_paper.sh`) mirrors the other sites' panels A–E and adds
panel **F**, the DiD forest — the interaction t per channel for both controls and
the placebo, the differencing assumption and the pre-trend parallelism — so the
whole contrast, and its honest null, is readable in one place.

### `audit_ocv_units.py` — what `A` counts, checked on every chip

The one script at this level, and the only one that reads every site at once for a
single number: [`audit_ocv_units.py`](audit_ocv_units.py) asks what `A` **is**.

`A` is a count of cells of a raster, and the raster has two sizes: the band
(`support_cells`, e.g. CTS 225) and the chip (`rows x cols`, CTS 39 x 145 = 5 655).
`A` counts the mask over the whole chip, so `A > support_cells` is normal (CTS 74
of 178 rows, Carola 27 of 436, Espoo 45 of 75, Haltern 19 of 107) and reading it as
"cells of the support band" makes a consistent record look inconsistent. The
script derives the unit from the data instead of assuming it:

| check | what it fixes |
| --- | --- |
| `support + guard + ring == rows x cols` | the three bands partition the chip |
| `A` integer, `0 <= A <= rows x cols` | `A` is a cell count of *this* raster; 0 is an empty mask, which Bautzen's control series really has (`has_mask` false) |
| `A / D == bbox_area` (1e-9) | `D` is `A` over the mask's bounding box, so `A` is that quotient's numerator |
| `A` vs `coherence_masked_pixels` and `A_committed` | the export column is another measurement (Werder: constant 640, unequal on all 991 rows); agreement with *both* on every row means the column was copied, the one defect this catches — LUMO agrees with the export on 178 of 178 rows but with `A_committed` on 5, i.e. the same rule on the same strip, recomputed |
| `mask_rows`/`mask_cols` vs the declared chip | the mask is measured on the chip the record declares |

Output `data/ocv_channel_units.json`: the contract, and per site the chip, the band
sizes, `A_range`, `D_range`, the `A > support_cells` count and the mask geometry.
The rule it enforces for every report: quote `A` as "cells of the chip (`<= N`)",
never beside `support_cells` without saying which raster each number counts, and
never across sites — each site has its own chip.

```sh
python3 code/audit_ocv_units.py                      # audit, print, write nothing
python3 code/audit_ocv_units.py --write-reference    # (re)write data/ocv_channel_units.json
python3 code/audit_ocv_units.py --verify             # audit vs. the committed reference
python3 code/audit_ocv_units.py --selftest           # the rules on synthetic rows
```


### `concept/fig_ocv_schematic.py` — the definition, drawn

The only script in `code/` that reads no data.
[`concept/fig_ocv_schematic.py`](concept/fig_ocv_schematic.py) draws the three-panel
schematic of the observability core vector — asset → SAR observation → measurement
support → echo mask → the six components — at a fixed 183 × 61 mm, and writes the
630 dpi JPEG, the vector PDF, the JSON build record and the Markdown report behind
it. It imports nothing from the site packages and opens nothing under `data/`: the
plate is the definition itself.

Because there is no record to validate it against, the script validates the
**drawing** instead:

| family | what it holds fixed |
| --- | --- |
| canvas | 183 × 61 mm at 630 dpi, isotropic units, every panel's contents inside its own box |
| text inventory | 26 text objects, each matched against a written whitelist |
| legibility | every string ≥ 7 pt and no two drawn strings overlapping |
| vocabulary | no `damage`, `health`, … and no `R_B`/`R_G` anywhere on the plate |
| set algebra | `W_{a,c} ⊂ guard ⊂ ring ⊂ C_{a,c}` with a guard ≥ 0.015 units and a tinted band ≥ 0.03, and `W` identical in (b) and (c) |
| the annulus | the composed canvas is sampled: `W`'s interior is the blue fill over white and the outer band the ring tint, each within 2/255 |
| the two frames | `X, Y, Z` and `az, r` share no symbol, the mapping statement carries `Π_c` and both frames, and `C_{a,c}` sits on the chip's own edge |
| the channel | every channel-dependent symbol carries the `c` subscript (`LOS_c`, `θ_c`, `Π_c`, `W_{a,c}`, `C_{a,c}`) and no frame symbol does |
| raster | four components of 8/6/4/3 cells, one closed hole, centroid and peak ≥ 1.5 cells apart (measured, never drawn), density 5–20 % |
| vector row | six components in order, each brace spanning its own group's measured extent, inter-group gap solved from those widths |
| determinism | two builds render identical pixels, and both files are byte-identical across runs |

A failure exits 2 *after* the artefacts are written, so a broken plate is reported
rather than silently delivered.

```sh
python3 code/concept/fig_ocv_schematic.py                 # both images + JSON + report
python3 code/concept/fig_ocv_schematic.py --no-figures    # build and check, write no image
python3 code/concept/fig_ocv_schematic.py --figdir DIR --json J --md M
python3 code/concept/fig_ocv_schematic.py --help          # the flags
```

