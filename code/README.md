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
first-class DiD layer (`cts_ocv_did.py`) next to the mask layer. YWF is the one
site whose record is too small and too duplicated for a contrast at all: one
mask layer, echo **coverage** instead of effect sizes, and a reference file that
the package itself defines and then pins. Espoo is the one site with **no event at
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
| [`ywf/`](ywf/README.md) | Yeongdeok Wind Farm, Unit 21, Samgye-ri (collapsed 2026-02-02) | pre-collapse (healthy), post-collapse | **one layer**: the echo mask of the committed 7 x 7 window payload (`0.30 x peak`, `5 x np.median`, ≥ 2 px) — the record is not segment-resolved, so a per-mast-section layer would be dishonest; its second cross-cut (the pass geometry, ASCENDING / DESCENDING) is reported as a pinned, explicitly descriptive echo-coverage split; two further references are *resolved* from public catalogues (the CDSE burst pass of every chip, the OSM road anchor) and committed beside it | `ywf/README.md` |
| [`espoo/`](espoo/README.md) | Espoo Kurttila mast, Espoo, Finland (**event-free** 2024-09-05 … 2026-09-07) | none — the record has no damage label at all; the contrast is the orbit geometry (ASCENDING n = 70 / DESCENDING n = 80) | **no mask layer**: the record delivers only the mask *size* `A` (`coherence_masked_pixels`) and no raster (56 of 78 columns empty), so the OCV is `[gamma2, A]` and `A` is **copied and proved**, never recomputed | `espoo/README.md` |

### `orbits/` — one OCV time-series figure per support-sized table

`[orbits/fig_ocv_orbits.py](orbits/fig_ocv_orbits.py)` draws the twelve-panel
companion (six where only one orbit has data — KDLO is DESC-only, so its plate
has one column and the report names the missing channel) for the seven sites
migrated to the support-sized chip (CTS, KDLO,
Carolabrücke, Haltern, YWF, Ponte Morandi, and Espoo): one panel per **dimension
x observation channel**, in the order and with the axes the LUMO and Bautzen
figures use (`gamma2, P, D, F, A, S`; `P, D, F, A` logarithmic).

Two states for the six labelled sites: **pre-collapse and post-collapse**, except
KDLO, whose record carries **pre-rebuild and post-rebuild**. A state with a
single overpass (KDLO's `single`) is dropped and named in the markdown. A channel
with no data is not drawn empty — the panel carries the reason (YWF's ASC track
has one healthy and no damaged overpass).

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
Morandi's row is four circles on purpose. YWF is **not** a row: nineteen of its 33
acquisitions cannot be cut on the support chip (four of the six post-collapse
dates among them), which is fewer than an interval — the report states it instead
of leaving a silent hole. `--selftest` checks the gate rules and the three
adapters; `--verify` re-derives every cell and every source hash.

```sh
python3 code/cases/fig_cases_overview.py            # JSON + report + 630 dpi JPEG
python3 code/cases/fig_cases_overview.py --selftest # 11 checks
python3 code/cases/fig_cases_overview.py --verify   # every cell, every sha256
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
`healthy`, Bautzen against `RS` on the paired bridge-minus-control difference —
so the two assets never share a scale and their magnitudes are never compared.

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
| **H1, per asset** | Jonckheere–Terpstra on the **per-unit** magnitudes (the reference's median magnitude is 0 by construction, so a test on medians would be degenerate), **Holm within each asset's own four-test family**, moving-block bootstrap with the reference re-estimated inside every draw |
| **verdict, per asset** | each asset gets its own verdict **and a bound** — the largest top-step change its interval still covers — because neither alone is a null result, only a bound |
| **comparison** | cosine of the standardised response vectors at the largest state, per block and channel, with a random-direction null and a sign-agreement count. **Descriptive**: it asks whether the same metric applied twice answers the same way, not whether a response transfers. γ² is also reported without γ², because it is a different quantity in the two experiments |

The 2026-10-22 amendment in the spec separates the two assets. Two consequences a
reader should know: the reference level carries **no** interval (`R = 0` by
construction, its noise belongs to `ΔR`), and where the reference state has fewer
than ten units — Bautzen has 6 (ASC) and 7 (DESC) complete pairs — the **pooled
MAD** is the standardising scale, because a MAD from six numbers is itself too
unstable. `scale_used` travels with every cell, in the JSON and in the table.

```sh
python3 code/sds/ocv_reconfigure.py                 # the pair, ~40 s
python3 code/sds/ocv_reconfigure.py --sites carola  # one asset, two states
python3 code/sds/ocv_reconfigure.py --selftest      # 11 checks on synthetic units
python3 code/sds/ocv_reconfigure.py --verify        # re-derives every table cell
```

`--sites` chooses the assets. **Carola** (`--sites carola`) is the two-state case:
pre-collapse against post-collapse on the support-sized chip, where only `Girder A`
survives the cut (436 of the committed 1 717 chips). With two states the ordered
test is its two-group case (a one-sided rank test), there is **no ordering claim**,
and the H2 comparison is skipped because a single asset has no partner — the report
says all three things rather than leaving blanks. Its post signature is a time
contrast without a control track, which the closing section states.

Outputs `data/sds/ocv_reconfigure.{json,md}`,
`data/sds/ocv_reconfigure_table.csv` and three 630 dpi plates with one panel guide:
`fig_sds_ocv_reconfigure.jpg` (the juxtaposition), `..._lumo.jpg` and
`..._bautzen.jpg` (each asset standing alone).


The scripts are started from the **repository root** (each wrapper in
`../figures/<site>/` does that), so that their defaults resolve to
`data/<site>/…` in and `data/<site>/…`, `figures/<site>/…` out. No package:
the script's own folder is on `sys.path` (a script run puts its directory first),
which is why the modules import each other by their `lumo_…` / `bautzen_…` /
`kdlo_…` / `carola_…` / `morandi_…` / `cts_…` / `ywf_…` / `espoo_…` names.

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
whole contrast, and its honest null, is readable in one place. The YWF package is
the family's smallest and the only one with a single mask layer: its 7 x 7 window
record holds 33 unique chips (24 windows are byte-identical duplicates across the
requested mast sections, `segment_resolved = false`), only 17 of them carry an
echo at all, and 1 of the 6 post-event chips does — so `fig_ywf_ocv_paper.py` fits
no classifier, leads with the echo-coverage table and its Fisher test, adds §A2
with the same coverage split by pass geometry (10 of 20 ascending vs. 7 of 13
descending chips, p = 1.0000 — pinned, but read as a description of where the
echoes sit), and
*defines* `data/ywf/reference/ywf_ocv_findings.json` with `--write-reference`
before pinning it on every later run. Two companions read the same record on a
time axis without touching that pin contract: `code/ywf/fig_ywf_ocv_weeks.py`
opens the pre-collapse window, and `code/ywf/fig_ywf_ocv_series.py` (run by
`figures/ywf/fig_ywf_ocv_series.sh`) draws the whole record — all 33 chips on 32
acquisition days, the 9 steps longer than the 6-day cadence as bands and the
collapse inside the 12-day step that holds it.

Two further references are not *defined* but *resolved* against public catalogues
and committed beside it — `reference/ywf_bursts.json` (the CDSE pass of all 33
chips: sub-swath, polarisation, relative orbit, azimuth time, with ASCENDING =
IW3 / rel-orbit 54 at 09:2x UT and DESCENDING = IW2 / rel-orbit 61 at 21:2x UT, VV
and VH interleaved in both) and `reference/ywf_osm_road.json` +
`reference/ywf_geometry.json` (the OSM way the tower fell on, its centreline and
the estimated burst grid) — each with its own offline `--check`
(`ywf_bursts_resolve.py`, `ywf_osm_anchor.py`) and its sha256 digest in the
figures' pin block. The package also states what the record **cannot** do: a
7 x 7 px window per acquisition (`TOWER_WINDOW_SIZE`), no rect spec of its own, no
per-section coordinates and no deck edge — the road lies *outside* the frame (its
near edge is 8.10-9.33 px across range from the window centre against a 6 px
reach, so containment needs 19-21 px per side) and a ~7 m carriageway would still
be 2.0 px across range but 0.5 px across azimuth — which is why the echo mask
stays the only layer it claims, and why the mask rule itself is only ever
*studied*: `ywf_mask_sensitivity.py` sweeps it offline on the committed payload
(the rule is family-wide — identical in the backend's `sarsolve` and in Carola's
`gamma2`).

