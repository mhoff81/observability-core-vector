# `code/` — the nine site packages

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
orbit geometry (ASC / DESC) — the family's event-free control. Werder Havel is the
family's **positive control in the other direction**: a record with documented
works instead of a collapse, whose three strata (`W2_bestand`, `W3_sanierung`,
`W4_post`) are cut on the two dates of the 2024/2025 renewal.

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
| [`werder_havel/`](werder_havel/) | Havelbrücke Werder, A 10 (documented works 2024-07-01 → 2025-09-05) | `W2_bestand`, `W3_sanierung`, `W4_post` | echo mask of the 80 x 80 complex window over the whole record, on the support-sized chip | `werder_havel_ocv_channels_csv.py` docstring (the package carries no separate README) |

### `cases/` — the standard OCV trajectory plate

`[cases/ocv_trajectory.py](cases/ocv_trajectory.py)` draws the **same 3 x 4 plate**
for every case site in its registry, on the support-sized chip, with the post-event
tail **drawn** instead of cut off — the timeline question the five cases are read
for, asked once and answered five times in one design. The sixth trajectory plate
of `figures/` (Werder Havel) is
`[werder_havel/fig_werder_havel_ocv_trajectory.py](werder_havel/fig_werder_havel_ocv_trajectory.py)`,
which implements the same template on an interval-specific registry. The reference
timeline is defined **case by
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

Two states are the rule; KDLO carries three and its `single` acquisition (the
collapse day itself) is excluded before any statistic. Standardisation is pooled
over the orbits at every site. The full rule set, the rejected windows and the panel
inventory are in [`ocv_trajectory_spec.md`](ocv_trajectory_spec.md).

```sh
python3 code/cases/ocv_trajectory.py --site carola   # one plate (+ JSON + report)
python3 code/cases/ocv_trajectory.py --all           # the five
python3 code/cases/ocv_trajectory.py --selftest      # 221 checks, writes nothing
python3 code/cases/ocv_trajectory.py --site kdlo --verify
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
(`behavior = [γ², P]` (cross-acquisition echo behavior) and `geometry = [A, D, F, S]` (spatial mask geometry, with `D` kept off the plane because it partly inherits `A` and `S` annotated because a shift is position, not morphology);
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
python3 code/sds/ocv_reconfigure.py --selftest              # 40 checks, the last ten offline
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
`data/sds/ocv_reconfigure_table.csv` and one plate with its panel guide:
`fig_sds_ocv_reconfigure.jpg` (the juxtaposition of both assets). The per-asset
companions of earlier revisions are no longer written: a pair emits the
juxtaposition only, and an asset named alone on `--sites` still gets its own plate.
A plate is authored at 170 mm — the
repository's text width — as four rows (ASC above DESC, per experiment) by two
columns (behavior, geometry), written from that canvas as `.jpg` (630 dpi), `.pdf`
and `.eps`, with a `.md` guide whose JSON blocks audit the files that were written;
spec amendment 6 records why the pre-registered 2 × 4 grid had to give way (four
columns leave 30 mm per plane, less than a 6.2 pt state label). The printed block
name is American English (`behavior`), while the key stays `behaviour` in `BLOCKS`,
the JSON and the audit, because the other generators index it by name — amendment 7
and the `BLOCK_DISPLAY` table are where that is fixed. Each registered pair
writes its own
stem — `data/sds/ocv_carola_morandi.{json,md}`, `..._table.csv` and the plate
`fig_sds_ocv_reconfigure_carola_morandi.jpg`; `data/sds/ocv_haltern_kdlo.{json,md}`,
`..._table.csv` and `fig_sds_ocv_reconfigure_haltern_kdlo.jpg` — and never touches
the dose files. Only the dose pair's stem is committed; the other two pairs are
runs of the same code and are not shipped.


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
Its shipped plate is the trajectory one (`figures/cts/fig_cts_ocv_trajectory.md`,
drawn by `code/cases/ocv_trajectory.py --site cts`), and the DiD numbers keep their
pinned home in `code/cts/cts_ocv_did.py`: the interaction t per channel for both
controls and the placebo, the differencing assumption and the pre-trend
parallelism, all recomputed from the committed references on every run.

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

### `figure_style.py` — one type ladder and one print width

Every figure in `figures/` is printed at the same width (**170 mm**) and has to
carry the same point sizes at that width, but until this module existed each
generator carried its own `COL_WIDTH = 183 * MM` and its own `rcParams` block.
Because a point size is *authored* and not *printed*, a canvas that
`bbox_inches="tight"` widens after the layout -- a long panel title, a wide
legend, a colour bar -- shrinks every glyph with it: the same class of panel came
out between 3.2 pt and 13.4 pt depending on which generator drew it and how much
text it had to wrap. [`figure_style.py`](figure_style.py) fixes both ends:

| role | what it is | value | band |
| --- | --- | --- | --- |
| `base` | `font.size`, the anchor | 8.0 pt | 8.0 (not free) |
| `title` | `axes.titlesize`, the panel title | 8.8 pt | 8.5 – 9.0 |
| `suptitle` | `figure.titlesize`, the one line above the grid | 9.0 pt | 8.5 – 9.0 |
| `label` | `axes.labelsize`, the axis label | 7.6 pt | 7.0 – 8.0 |
| `tick` | `xtick`/`ytick.labelsize` | 6.8 pt | 6.5 – 7.0 |
| `legend` | `legend.fontsize` | 6.8 pt | 6.5 – 7.0 |
| `annotation` | `text`/`annotate` inside a panel | 6.2 pt | >= 6.0 (a floor) |

and the canvas is `FINAL_WIDTH_MM = 170` with a tolerance of ± 0.5 mm, checked
on the **written file**, not on the intent:

```python
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
for path in (os.path.join(REPO, "code"),):      # the same preamble as sds_support.py
    if path not in sys.path:
        sys.path.insert(0, path)

from figure_style import FS, COL_WIDTH, apply, assert_fits, grid_rect, save  # noqa: E402
apply()                                          # the rcParams block, once
fig = plt.figure(figsize=(COL_WIDTH, 0.9 * COL_WIDTH))
...
leg = fig.legend(...)                            # a figure artist, outside the grid
fig.tight_layout(rect=grid_rect(fig, bottom=leg))
assert_fits(fig)                                 # canvas is 170 mm, tight box fits
save(fig, path)                                  # write at that canvas, measure it
```

A generator that already inserts its own directory (`sys.path.insert(0, HERE)`)
adds one line instead: `sys.path.insert(1, os.path.dirname(HERE))`.

`save()` writes the canvas **uncropped** and then measures the file it wrote: a
cropping save (`bbox_inches="tight"`) makes the file as wide as the layout
happens to be, which is the defect itself -- the reader would get the ladder at
a size nobody authored. `assert_fits()` runs first, so a layout that needs more
than the canvas (too long a title, a legend that will not fit, a grid that has to
be re-cut) is refused instead of cropped, and the width assertion then checks
the written file at 170 ± 0.5 mm.

### `figure_style.py` — the same plate in three containers

Beside the raster, under the same stem, `save()` writes the plate as `.pdf` and
as `.eps` (`VECTOR_FORMATS`) and measures **all three**: a sidecar that is
cropped or padded is as much a defect as a cropped raster. The PDF is the copy
for LaTeX -- full vector, transparency intact, no `/CreationDate` or `/ModDate`
so two runs are byte-identical. PostScript has no alpha at all, so the EPS has to
*resolve* it, and which way it does is the plate's decision (`eps=`):

| `eps=` | what the EPS carries | when |
| --- | --- | --- |
| `"blend"` (default) | the tint the PDF shows **on white**: `alpha*C + (1-alpha)*background` per channel, then alpha 1. Text stays vector | the plate's translucent artists are tints over the page -- fills, faint guides, annotation arrows |
| `"raster"` | the PDF itself through poppler's `pdftops`: pixel-identical to the PDF, but the page becomes one 300 dpi bitmap | alpha **is** the reading -- overlapping translucent marks, where a density has to accumulate (the orbit lines at 0.55: `code/lumo/fig_lumo_ocv_orbits.py`) |
| `"opaque"` | matplotlib's own flattening, i.e. the raw colour at full strength | never in this repository; kept as the documented default of the PS backend |

`"blend"` is the default because matplotlib's PS backend does the third row
whatever you ask: it writes the *unmixed* colour opaque, so an arrow at alpha 0.75
comes out at full strength next to a PDF that shows the tint. `blend_alpha(fig)`
does the mixing instead, and `save()` then re-reads the figure and **refuses to
write** if anything is still translucent -- a blend that silently missed an
artist would be the same defect one layer down:

```
  eps 170.00 mm  blend, 74 artist(s) blended
```

Two artists have to be reached by name, because they are in no artist tree:
`Text.get_children()` and `Annotation.get_children()` return `[]`, so a label's
`bbox_patch` and an annotation's `arrow_patch` are invisible to `fig.findobj()`
although the backend draws them like any other patch. Missing the arrow is not
cosmetic: its alpha sits *in the arrow*, so the EPS would carry the raw colour
while the PDF shows a tint (that is the 18-colour leftover the guard exists for).

A blend cannot accumulate: two translucent marks on top of each other become two
opaque marks at their single-layer composite. That is why a density plate takes
`"raster"`, and why the choice is made per plate at the call site rather than in
this module. Nor does it look *through* a mark: a translucent white label box
over a line becomes the white box, where the PDF still shows the line at 22 %.
Both containers are gitignored (`figures/**/*.pdf|.eps`): they are regenerated
binaries, like the rasters.

### `figure_style.py` — the layout a fixed canvas needs

A fixed canvas changes one more thing that has nothing to do with type: the
`bbox_inches="tight"` save was also what made a hand-tuned
`tight_layout(rect=...)` look right. With a fixed canvas the same rect leaves a
band of white as tall as the caption, because `tight_layout` **already reserves
the suptitle** — a rect that ends below the suptitle counts it twice:

```python
fig.tight_layout(rect=(0, 0, 1, 0.900))                 # 37 mm of white under the caption
fig.tight_layout(rect=style.grid_rect(fig, top=leg))    # the grid fills the canvas
```

`grid_rect(fig, bottom=..., top=...)` returns that rect and answers the second
half of the problem: `tight_layout` knows only about *axes* decorations, so a
`fig.legend()` strip under the grid or a `fig.text()` note above it is invisible
to it and gets drawn across the axis labels unless it is declared. Pass the
artist — its drawn bbox is measured, in millimetres — or an explicit height;
`None` means "the canvas edge".

Two rules follow for a migrated generator:

* place every figure-level artist **inside** the canvas. `y=1.040` was legal only
  while a tight bbox grew the canvas upward; on a fixed canvas it is clipped, so
  the header block (`suptitle`, legend strip, note) moves inside and
  `grid_rect(..., top=leg)` hands the rest of the height to the grid;
* a caption that does not fit is a caption to **re-wrap**, not a canvas to
  widen. The three-line Haltern note at 6.2 pt is 130 mm wide; the one-line form
  it replaced was 206 mm — that 206 mm *was* the 235 mm file.

### Migrating a generator onto `figure_style`

Five edits per generator, then a run, then `audit_figure_style.py`:

1. `sys.path.insert(1, os.path.dirname(HERE))` + `import figure_style as style`
   next to the site imports;
2. the local `plt.rcParams.update({...})` block becomes `style.apply(dpi=...)`,
   keeping only the drawing convention (spines, grid, mathtext, on-screen dpi) as
   overrides — `savefig.bbox` is *gone*, not overridden;
3. `COL_WIDTH = 183 * MM` → `COL_WIDTH = style.COL_WIDTH` (and a
   `P.cmp("figure.width_mm", ...)` pin moves to `style.FINAL_WIDTH_MM`);
4. `fontsize=` / `labelsize=` / `titlesize=` literals → `style.FS["title" |
   "label" | "tick" | "legend" | "annotation"]`; where a size was picked to fit a
   narrow panel, re-cut the *text* (a shorter axis label), not the ladder;
5. `fig.savefig(...)` → `style.save(fig, path)` and a hand-tuned
   `tight_layout(rect=...)` → `style.grid_rect(...)`; then re-cut whatever
   overflows, which `assert_fits` reports in millimetres.

A *family* of plates — several sites on one template — takes two more edits, both
in `figure_style.py` rather than in the generators:

* the frame is a **table**, not a copy: `TRAJECTORY = {width_mm, height_mm, grid,
  legend}` holds the canvas, the grid margins and the band the key may occupy, so
  two generators drawing the same template cannot deliver two different plates.
  The generators read it (`fig.add_gridspec(3, 4, **style.TRAJECTORY["grid"])`);
* the **pen** goes with the type: `style.LW` is the line ladder (series 1.4 pt,
  plane/rule/spine 1.0, grid 0.6, marker edge 0.9) and `style.MS` the three marker
  sizes. They are overrides, not `rcParams`, so a plate already written keeps its
  pens when it adopts the ladder.

And one thing a plate should *lose* on the way: its titles. A plate that carries a
header and a descriptive title per panel is a plate whose text has to be re-flowed
every time the canvas changes. The six trajectory plates now carry their bare tags
(`a1`…`c4`), one legend, and a `## Caption` section in the report — the caption is
set by the journal, at the journal's width, and the plate is left to the data.

Done so far (both artefacts measured from the written file):

| generator | plate | width | containers |
| --- | --- | --- | --- |
| `lumo/fig_lumo_ocv_orbits.py` | `figures/lumo/fig_lumo_ocv_orbits.jpg` | 170.0 mm | jpg + pdf + eps (`raster`) |
| `espoo/fig_espoo_ocv_chips.py` | `figures/espoo/fig_espoo_ocv_chips.jpg` | 170.0 mm | jpg + pdf + eps |
| `sds/ocv_reconfigure.py` | `figures/sds/fig_sds_ocv_reconfigure.jpg` | 170.0 mm | jpg + pdf + eps (`blend`, 60) |
| `cases/ocv_trajectory.py` (five sites) | `figures/{kdlo,carola,morandi,cts,haltern}/fig_*_ocv_trajectory.jpg` | 170.0 mm | jpg + pdf + eps (`blend`, ~150-175) |
| `werder_havel/fig_werder_havel_ocv_trajectory.py` | `figures/werder_havel/fig_werder_havel_ocv_trajectory.jpg` | 170.0 mm | jpg + pdf + eps (`blend`, 167) |
| `bautzen/fig_bautzen_ocv_orbits.py` | `figures/bautzen/fig_bautzen_ocv_orbits.jpg` | 170.0 mm | jpg + pdf + eps (`blend`, 217) |

The six trajectory plates share `style.TRAJECTORY` — 170.0 × 156.0 mm, a 3 × 4 grid
with a reserved legend band — because they are one template drawn twice, so a
per-generator `figsize` was exactly the thing that let them drift to 190.2 … 193.3 mm.
They also carry no header and no descriptive panel titles any more: a bare tag
(`a1`…`c4`) per panel, one legend inside the band, and the sentence a reader needs
in the report's `## Caption`. The census a header used to print is a table in the
report, which is the only place it can be re-flowed without touching the canvas.

`bautzen/fig_bautzen_ocv_orbits.py` came last, and it was the one plate off the
ladder: it drew a four-column grid at its own column width — 180.9 mm — and set
type below the ladder's floor (base 8.5 pt against the 8.0 pt anchor, titles
8.0 pt). It now runs on the same frame as the rest (`style.apply()`,
`style.COL_WIDTH`, `grid_rect()`, `style.save()`), so all ten print plates are
170.0 mm wide — the six trajectory plates on the shared template and Bautzen on
its own 6 x 2 grid. The concept pair is the declared exception: it draws its own
183 x 61 mm page (three panels of one isotropic unit) with its own type tables,
so the audit lists its declared sizes as *symbolic* and
`audit_figure_style.PAGE_PLATES` holds it to that page instead of to 170 mm.

`figure_style.py` imports matplotlib only inside `apply()`, `assert_fits()` and
`save()`, so the constants and the image reader (PNG `pHYs`, JPEG `JFIF`/`SOF`,
stdlib only) are usable on an interpreter without matplotlib:

```sh
python3 code/figure_style.py figures/lumo/fig_lumo_ocv_orbits.jpg
```

### `audit_figure_style.py` — the ladder and the width, checked on every figure

[`audit_figure_style.py`](audit_figure_style.py) is the gate on the section
above, and it checks the two things independently, because a figure can be wrong
in two independent ways:

| what | how | fails when |
| --- | --- | --- |
| the declared type | every `code/**/*.py` that writes a figure is parsed (AST, not grep) -- a `savefig` of its own or the shared `figure_style.save`/`save_eps`: the `rcParams` block, `rcParams[...] = ...` and every `fontsize=`/`labelsize=`/`titlesize=` keyword, attributed to its call (`set_title` -> title, `set_xlabel` -> label, `tick_params` -> tick, `legend` -> legend, `text` -> annotation) and to its line. Values are resolved through the file's own numeric constants and through `FS["tick"]` | a value sits outside its band |
| the print width | the saved width of every artefact is measured from the file itself (PNG `pHYs`, JPEG `JFIF`, EPS `%%HiResBoundingBox`, PDF `/MediaBox` with `/CropBox` as fallback, stdlib only) | any container of the plate is not its expected width (170 ± 0.5 mm, or the page it declares in `PAGE_PLATES`), or a raster declares no density at all |

A plate is one row, not one row per file: the raster, the PDF and the EPS are the
same drawing in three containers, so `artefacts()` groups them by stem and prints
the containers side by side -- `jpg+pdf+eps` when they agree, and
`png 179.8 + jpg 183.0` when they do not. The verdict is **per container**: one
off-width container fails the plate and is named, because a PDF at 170 mm next to
a 235 mm PNG is exactly the defect the fixed canvas is there to prevent. A vector
container is judged on its page box in points, never on a dpi: it has no pixels.

The report prints both at once with the **effective** size at 170 mm, so a 9.0 pt
title on a 375 mm canvas is visible as the 4.1 pt a reader would actually see:

```
PASS  code/bautzen/fig_bautzen_ocv_orbits.py             base 8.0  title 8.8  suptitle 9.0  label 7.6  tick 6.8  legend 6.8  annotation 6.2
PASS  bautzen/fig_bautzen_ocv_orbits.jpg                  170.0 mm  jpg+pdf+eps       -> title 8.8  label 7.6  tick 6.8  legend 6.8
PASS  carola/fig_carola_ocv_trajectory.jpg                170.0 mm  jpg+pdf+eps
note  concept/fig_ocv_schematic.jpg                       183.0 mm  jpg+pdf+eps       page 183 mm
```

The `->` clause is the **effective** size at 170 mm. Bautzen is the plate that
made the argument: saved on its old 180.9 mm canvas with an 8.0 pt title, the same
file at the print width showed a reader 7.5 pt. It now declares 8.8 pt on a
170.0 mm canvas, so declared and effective agree — and a plate that drifts back to
a wider canvas is reported at the size it would really print, not at the size it
asked for. That is the whole case for the fixed canvas in one line.

Scope: the site directories in `SCOPE_DIRS` are failed on, and `WALL_PLATES` is
empty now that the oversized panels are retired — the hook stays, but no plate is
exempt from the width check any more. `figures/concept/` and `figures/sds/` are
documentation plates that are measured and reported only — `--strict` fails on
those two groups as well.

The one exception runs the other way: `PAGE_PLATES` registers the two concept
plates as **183 mm pages** (three panels of one 55.632 mm unit), the deliberate
departure from 170 mm. They are not exempt from the width check, only from that
number: every container has to agree with the declared page, the row prints
`page 183 mm` next to it, and a 170 mm EPS of a 183 mm page still fails. So
`--strict` exits 0 with 12 plates in 36 files — 10 at 170.0 mm, 2 at their own
183.0 mm.

```sh
python3 code/audit_figure_style.py                    # audit, print, write nothing
python3 code/audit_figure_style.py --only kdlo        # one site (repeatable)
python3 code/audit_figure_style.py --strict           # wall and documentation plates fail too
python3 code/audit_figure_style.py --json /tmp/style.json   # machine-readable
python3 code/audit_figure_style.py --selftest         # the rules on synthetic input
```

Exit status 0 on a clean audit, 1 on any violation (the convention of
`audit_ocv_units.py`).


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

