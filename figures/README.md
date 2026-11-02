# `figures/` — the twelve plates

One folder per site (plus `concept/`, `sds/`, `cases/`), and inside a folder one
**stem** per plate:

```
figures/<folder>/<stem>.sh    the wrapper: runs that plate's generator from the repo root
figures/<folder>/<stem>.md    the generated English report — pins, panels, caveats
figures/<folder>/<stem>.jpg   the raster at 630 dpi
figures/<folder>/<stem>.pdf   the same drawing as a vector page
figures/<folder>/<stem>.eps   the same page in PostScript (alpha pre-blended over white)
data/<folder>/<stem>.json     every number the plate draws
```

`cases/` holds no plate of its own: it holds the driver that runs five of the
trajectory plates (Figs 5, 6, 10–12) in one go.

The list is **closed**. Twelve plates are shipped here; a plate that is not in
this table is not in this repository.

| # | plate | folder | generator | renders |
| --- | --- | --- | --- | --- |
| 1 | `fig_ocv_schematic` | `concept/` | `code/concept/fig_ocv_schematic.py` | `.jpg` `.pdf` `.eps` |
| 2 | `fig_ocv_measurement` | `concept/` | `code/concept/fig_ocv_measurement.py` | `.jpg` `.pdf` `.eps` |
| 3 | `fig_espoo_ocv_chips` | `espoo/` | `code/espoo/fig_espoo_ocv_chips.py` | `.jpg` `.pdf` `.eps` |
| 4 | `fig_sds_ocv_reconfigure` | `sds/` | `code/sds/ocv_reconfigure.py` | `.jpg` `.pdf` `.eps` |
| 5 | `fig_kdlo_ocv_trajectory` | `kdlo/` | `code/cases/ocv_trajectory.py --site kdlo` | `.jpg` `.pdf` `.eps` |
| 6 | `fig_carola_ocv_trajectory` | `carola/` | `code/cases/ocv_trajectory.py --site carola` | `.jpg` `.pdf` `.eps` |
| 7 | `fig_werder_havel_ocv_trajectory` | `werder_havel/` | `code/werder_havel/fig_werder_havel_ocv_trajectory.py` | `.jpg` `.pdf` `.eps` |
| 8 | `fig_lumo_ocv_orbits` | `lumo/` | `code/lumo/fig_lumo_ocv_orbits.py` | `.jpg` `.pdf` `.eps` |
| 9 | `fig_bautzen_ocv_orbits` | `bautzen/` | `code/bautzen/fig_bautzen_ocv_orbits.py` | `.jpg` `.pdf` `.eps` |
| 10 | `fig_morandi_ocv_trajectory` | `morandi/` | `code/cases/ocv_trajectory.py --site morandi` | `.jpg` `.pdf` `.eps` |
| 11 | `fig_cts_ocv_trajectory` | `cts/` | `code/cases/ocv_trajectory.py --site cts` | `.jpg` `.pdf` `.eps` |
| 12 | `fig_haltern_ocv_trajectory` | `haltern/` | `code/cases/ocv_trajectory.py --site haltern` | `.jpg` `.pdf` `.eps` |

Two plates (`3`, `4`) are run through their generator directly — they have no
`.sh` because they were never shell-driven; every other plate has its wrapper.

## Three kinds of plate

* **A definition plate** (Figs 1–2) is computed from the *definition* the records
  are measured with, not from a record. It has no table input and therefore no
  pin to a table: it carries structural self-checks instead.
* **A contrast plate** (Figs 3, 4, 8, 9) compares *states* — two pass geometries at
  Espoo, the controlled states of LUMO and Bautzen, the documented states of the
  same structure — and standardises inside the asset so that two panels can be
  compared in shape.
* **A timeline plate** (Figs 5–7, 10–12) answers *where in the record* the vector
  moves: one reference state on the left, the post state drawn, never cut off, and
  a monthly median per dimension with its `n`.

## The rules every plate keeps

1. **One cut, named and hashed.** Every plate that reads a record reads exactly
   one committed table, named in its report together with its sha256, and never
   mixes two chips: the support-sized chip and the committed export are two
   measurements.
2. **Every number is pinned.** Each generator recomputes every number it draws
   from its table and the frozen inputs, and refuses to write — or, where it has
   one, fails its `--verify` — on a single deviation. Two runs of one generator
   produce byte-identical files: the bootstrap is seeded and no wall-clock value
   enters any output.
3. **One page, one measured width.** Every raster is a 630 dpi JPEG, and beside it
   stand the same drawing as `.pdf` and as `.eps` (`figure_style.save`), all three
   measured after they are written. The ten print plates are authored at the
   170 mm text width of `code/figure_style.py` — the six trajectory plates on the
   one shared frame of **170.0 × 156.0 mm**, `figure_style.TRAJECTORY` — and are
   never cropped to their content: a layout that does not fit fails the run rather
   than growing the canvas. The concept pair (Figs 1–2) is the declared exception:
   a **183 mm page** of three equal panels, registered in
   `audit_figure_style.PAGE_PLATES` and held to that page in every container
   instead of to 170 mm. A plate carries no title and no panel titles of its own:
   it carries its bare tags (`a1`…`c4`), one legend in the band its grid reserves,
   and a `## Caption` in its report.
4. **No verdict.** A plate reports a departure, a direction and an interval. It
   does not report damage, failure, a precursor or a cause, and the word it does
   print (`resolves …`, `sustained … departure`) is a function of frozen rules
   stated in the report, never of the shape of the curve.

## 1–2. The concept pair — `concept/`

**`fig_ocv_schematic`** — the observability core vector, drawn: three panels that
carry the reader from the physical asset, through the SAR observation geometry and
the asset-specific measurement support, to the radar echo representation and the
six components. No input file at all; 86 structural checks stand in for the pin.

**`fig_ocv_measurement`** — one chip, three stages, one number set: the native read
chip `C`, the measurement support `W`, the echo mask `M` and finally `OCV(t)`, each
stage admitting less than the one before. 148 checks, 24 text objects, 20 echo
cells in 4 regions, and a stated pixel size (55.632 mm) for the drawn chip.

Both are read as the *definition* of the other ten plates — no site's numbers
enter either of them.

## 3. `espoo/fig_espoo_ocv_chips` — one chip site, six dimensions

The Espoo Kurttila mast on the **support-sized chip**: 75 chips
(2024-09-05 … 2026-09-07, ASCENDING 31 / DESCENDING 44, three shared dates), one
170.00 × 191.77 mm plate. It is the only plate that reads a site's chip table with
no state column and no event: what it has instead is the **pass geometry**.

* **Two tracks.** ASC is relative orbit 160, DESC 51; the other two tracks place
  the mast on the far range edge, where the 200 m ring runs past the swath.
* **No pairing.** The tracks share three dates, so the orbit contrast is a
  two-sample comparison, and the report says so.
* **A tautology, stated.** `D = A / bbox_area` exactly, and `gamma2_mask` carries
  `n = A` in its denominator; `A` differs by a factor 7 between the orbits, so both
  dimensions are dominated by mask size.
* **A second instrument.** These are chip values, not the pipeline columns: the
  two agree only loosely (r(log A, log `coherence_masked_pixels`) = −0.24,
  r(`gamma2_mask`, `coherence_gamma2`) = 0.13). Neither is wrong; they measure
  different things. The pipeline figure that produced the second reading is no
  longer shipped; its frozen lock stays at `data/espoo/fig_espoo_ocv_orbits.json`.

## 4. `sds/fig_sds_ocv_reconfigure` — two assets, one metric

The cross-site plate: the **controlled-dose** experiments of LUMO (healthy, DAM 3,
DAM 4, DAM 6 on a 0 / 7.4 / 13 % frequency drop) and of Bautzen (`RS`, `DS1`,
`DS2`) side by side, drawn as one juxtaposition — rows = experiment × channel (ASC
above DESC), columns = block (behaviour `[gamma2, P]` and geometry
`[A, D, F, S]`), every panel on its own plane.

In every panel the origin is that asset's **reference state** and the arrows end at
its standardised responses (`z = (median − median_ref) / scale`), so the geometry is
within-asset and no scale is shared between assets. The inset repeats the block
magnitude `R` per level with its moving-block 95 % interval — the quantity the
hypotheses test, each asset inside its own Holm family.

Its reading is deliberately a **null about ordering**: no dose ordering survives in
either block at either asset (`P1` does not hold: LUMO's geometry magnitude goes
0 < 0.581 < 0.34 ASC, 0 < 0.419 < 0.225 DESC), the same-dose states `DAM 3` /
`DAM 4` differ more than the cross-dose levels do in the geometry block (`P2`
holds for geometry, partly for behaviour), and the behaviour block need not follow
the dose at all (`P3` consistent: 0 of 4 behaviour tests reach Holm < 0.05). The
cross-asset comparison is descriptive — the same metric applied twice, not a
transferability test.

## 8–9. The orbit plates — one panel per dimension × orbit

**`lumo/fig_lumo_ocv_orbits`** — six dimensions × two orbit directions, twelve
panels over the tower's 178 coherence rows joined to its 482 acquisition rows. The
x coordinate is the **actual local overpass instant** (`Europe/Berlin`, joined by
`burst_id`), not midnight, and the background bands follow the **13 consecutive
`damage_label` runs** of the time-sorted cohort, with the twelve run transitions
dashed — the states are read where the record puts them.

**`bautzen/fig_bautzen_ocv_orbits`** — six rows in the formal OCV order
(`gamma2`, `P`, `D`, `F`, `A`, `S`) × two columns (ASC IW2, DESC IW3) over the
bridge's `RS` / `DS1` / `DS2` states, again on the actual overpass instant, because
the state boundaries fall mid-day. The control rectangle 500 m east is **not**
drawn: its ASC `gamma2` range (−0.94 … +0.93) is an order of magnitude wider than
the bridge's (−0.14 … +0.24) and would compress the signal the plate exists to
show. It stays in `data/bautzen/bautzen_ocv_channels.csv`.

Both plates are the *time-series* reading of a table the site's own package builds;
neither computes a pooled contrast.

## 5–7, 10–12. The trajectory family — six plates, two generators, one design

Six plates over six sites, one design: five of them are emitted by
`code/cases/ocv_trajectory.py` from a registry of sites, the sixth (Fig 7, the
Werder crossing) by `code/werder_havel/fig_werder_havel_ocv_trajectory.py`, which
implements the same template with an interval-specific registry. The question is
the timeline one — *where in the record does the vector move, and does it move
where the event is?* Every other plate of this repository answers a contrast
question; these answer "when".

### The rules the report states (R1–R6)

* **R1 — one cut.** The plate reads the **support-sized** chip
  (`data/<site>/support/<site>_ocv_channels.csv`), the same cut the SDS plate
  standardises with. The committed export is never mixed in: two chips are two
  measurements.
* **R2 — the case's own reference, and the window is chosen per asset.** The
  reference timeline is the case's own reference state (`pre` /
  `pre-collapse (healthy)` / the crossing's `W2_bestand`), restricted to the window
  the registry declares for that site. A healthy chip outside that window is a
  third class: drawn, in the reference frame, and standardising nothing.
* **R3 — per-block standardisation.** `behaviour = [gamma2, P]` against the
  reference chips that carry an echo; `geometry = [A, D, F, S]` against every
  reference chip. Pooled over the orbits at every site.
* **R4 — a robust z.** Transform first (log for `gamma2, P, A, D, F`, identity for
  `S`), then `z = (x − median_ref) / (1.4826 · MAD_ref)`, falling back to the IQR
  and then to the SD, clipped at ±6; a block magnitude is the RMS of its
  standardised dimensions.
* **R5 — no interpolation.** A month with no chip is a hole: the line breaks there
  and never spans it, and neither does the event. The monthly median is drawn where
  a month holds a chip, with its `n`.
* **R6 — the outside class is tested too.** The healthy chips the window leaves out
  are drawn and read against the reference with the same interval, so the plate
  never shows a chip it does not also measure.

### The twelve panels

`a1` z(gamma2) · `a2` z(P) · `a3` behaviour plane · `a4` R_B · `b1` z(A) · `b2` z(F)
· `b3` geometry plane · `b4` R_G · `c1` weather · `c2` observation cadence · `c3`
echo coverage · `c4` the trajectory matrix (six dimensions, one row each).

### The six plates

| # | site | kept chips | reference (chosen) | tail | what the frozen rules print |
| --- | --- | --- | --- | --- | --- |
| 5 | KDLO media tower | 59 | `pre`, 16 (2022-06-02 → 2022-12-11) | `rebuild`, 43 | resolves `P`, `A`, `S`, `R_B`, `R_G` |
| 6 | Carolabrücke, Dresden | 436 | `pre-collapse (healthy)`, 121 of 184 (window 2022-09 … 2023-12) | `post-collapse`, 252 | resolves `gamma2`, `P`, `S` · held-out class: isolated unusual departure |
| 7 | Havelbrücke Werder, A 10 | 491 of 991 rows | `W2_bestand`, 119 | `W3_sanierung` 189 / `W4_post` 61 | sustained unusual geometry departure, confined to the intervention interval |
| 10 | Ponte Morandi, Genova | 207 | `pre-collapse (healthy)`, 198 | `post-collapse`, 9 | resolves `D`, `R_G` |
| 11 | Champlain Towers South | 178 | `pre-collapse (healthy)`, 136 | `post-collapse`, 42 | resolves `P`, `A`, `D` |
| 12 | Haltern AV9 turbine | 107 | `pre-collapse (healthy)`, 96 | `post-collapse`, 11 | resolves nothing |

What each plate adds to its row:

* **The frame is one table, not six choices.** All six are drawn on
  `figure_style.TRAJECTORY` — 170.0 × 156.0 mm, three rows of four, one legend band
  under the grid — and written by `figure_style.save` as `.jpg` (630 dpi) + `.pdf` +
  `.eps`, each measured against that frame before the run is allowed to finish. The
  type is the repository ladder (panel tag 8.8 pt, axis label 7.6 pt, tick and legend
  6.8 pt) and the pen is `figure_style.LW` (series 1.4 pt, rule 1.0, grid 0.6). The
  plate carries **no header and no descriptive panel titles**: each panel is tagged
  with its bare key (`a1`…`c4`), the single legend names every mark, and the census
  the header used to carry is a table in the report, under `## Caption`.

* **Fig 5, KDLO.** The tail is the **rebuild**, not a collapse aftermath, and the
  plate's caption says so rather than a footnote: the contrast answers a different
  question. The reference is 16 chips only, and the acquisition at the
  collapse itself (2022-12-23) is excluded by the registry because one acquisition
  is not a state.
* **Fig 6, Carola.** The window is a *choice* (2022-09 … 2023-12, 121 of the 184
  healthy chips), so 63 held-out chips are drawn and measured in their own block.
  The verdict is an **isolated unusual departure**, and the plate also prints the
  one directional statement it can carry: `z(P)` below the reference envelope in 5
  of 9 held-out months, with the joint `z(gamma2) > 0, z(P) < 0` pattern holding for
  four consecutive months — direction, not size, and not a verdict input.
* **Fig 7, Werder Havel.** Not a collapse but a **documented works interval**
  (2024-07-01 → 2025-09-05, 431 d) over three strata of the export
  (`W2_bestand` / `W3_sanierung` / `W4_post`), with the held-out pre-intervention
  year as the outside class. It answers nine registered questions in order and
  keeps three statements apart: directional departure, departure magnitude and
  sustained departure. `R_G` sits above the monthly reference q95 in 10 of the 15
  interval months and for three consecutive months (2025-02 → 2025-04); one
  acquisition of the operator record is cross-checked against
  `/home/projects/havel_werder_analysis/havel_werder_events.json` **when that file
  exists** — the cross-check is skipped, not failed, when it does not.
* **Fig 10, Morandi.** 9 post chips against 198 reference chips: the tail is drawn
  and its contrast reported with the interval it can carry, which is wide by
  construction.
* **Fig 11, CTS.** One orbit only (ASCENDING, relative orbit 48): the orbit columns
  are empty on the other direction, so pooling is the whole record here, and the
  reference is the site's six-year pre-trend, deliberately not cut.
* **Fig 12, Haltern.** 13 of the record's 120 dates are refused at the raster edge,
  all ascending: the cut's own manifest carries them, and the plate does not
  pretend the missing 11 % are months without data.

## Running them

```bash
bash figures/concept/fig_ocv_schematic.sh
bash figures/concept/fig_ocv_measurement.sh
python3 code/espoo/fig_espoo_ocv_chips.py            # no .sh wrapper
python3 code/sds/ocv_reconfigure.py                  # no .sh wrapper
bash figures/kdlo/fig_kdlo_ocv_trajectory.sh
bash figures/carola/fig_carola_ocv_trajectory.sh
bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh
bash figures/lumo/fig_lumo_ocv_orbits.sh
bash figures/bautzen/fig_bautzen_ocv_orbits.sh
bash figures/morandi/fig_morandi_ocv_trajectory.sh
bash figures/cts/fig_cts_ocv_trajectory.sh
bash figures/haltern/fig_haltern_ocv_trajectory.sh
```

The five trajectory plates also run in one go — the same generator, one process per
site, nothing pooled across sites:

```bash
bash figures/cases/fig_ocv_trajectory.sh              # all five
bash figures/cases/fig_ocv_trajectory.sh --verify     # re-derive, write nothing
```

`--verify` is the check to reach for after touching a table or a generator: it
re-derives every number of the JSON and the report from the table and fails on one
deviation. `--no-figures` computes and writes the report only; `--quiet` prints
nothing but the failures. Every wrapper passes its arguments through, and every one
of them honours `PYTHON=/usr/bin/python3` where the `python3` of `PATH` is a
virtualenv without numpy/matplotlib.

## The audits that watch these plates

* `python3 code/audit_figure_style.py` — the type ladder and the width of every
  generator and of every written file (`--strict` also fails the documentation
  plates of `concept/` and `sds/`). The ten print plates read
  `PASS … 170.0 mm jpg+pdf+eps`: the width line fails a plate that came out 190 mm
  wide, the type line fails a 7.2 pt panel title, and the `jpg+pdf+eps` column fails
  a plate that wrote only the raster. The concept pair reads
  `note … 183.0 mm jpg+pdf+eps   page 183 mm`: it is registered in
  `audit_figure_style.PAGE_PLATES` and held to its own page in every container
  rather than to the 170 mm print width.
* `python3 code/audit_ocv_units.py` — what `A` counts, checked on every chip of
  every committed table.

## What is not here any more

The per-site "paper" plates (panels A–E / A–F), the month companions, the orbit,
geometry, series and weeks companions, the case-overview page and the
prospective/correspondence analyses were **retired** with this scope. No generator,
wrapper, report or image for them ships in `figures/`; the plates above are the
whole set. Their committed numbers stay under `data/` as frozen provenance — a
retired plate's JSON remains the record of what it drew, and two of those locks are
still read: `data/lumo/fig_lumo_ocv_paper.json` by the LUMO orbit plate and
`data/haltern/reference/haltern_ocv_findings.json` by the Haltern mask-sensitivity
script. The nine site packages still build, pin and report the numbers behind their
remaining plates (see `code/README.md`).
