# observability-core-vector

Nine offline site packages, one cross-site plate and a drawn definition of the
measurement itself: **twelve figures in total**, each with its own report, its own
committed inputs and its own pin block. They all measure the same thing on every
site — the **observability core vector**

```
OCV(t) = [gamma2, P, D, F, A, S]
```

| symbol | dimension | paper definition | reads as | block |
| --- | --- | --- | --- | --- |
| `gamma2` (drawn `γ²`) | coherence | consistency of the radar response across acquisitions | echo consistency | behavior `[gamma2, P]` |
| `P` | persistence | fraction of observations with a valid echo | observation continuity | behavior |
| `D` | density | fraction of the observation region occupied by echo-mask pixels (`D = A / bbox_area`, so it partly inherits `A` algebraically) | spatial concentration | geometry `[D, F, A, S]` |
| `F` | fragmentation | number of connected components in the echo mask (the deck-mask sites carry `1 − largest component / A` instead) | spatial organisation | geometry |
| `A` | area | spatial extent of the echo mask | observable extent | geometry |
| `S` | shift | displacement of the echo mask relative to a reference position | spatial displacement | geometry |

The six are what every plate draws and standardises, in the order the two
definition plates pin (Figs 1–2) — and none of them claims an answer.
`[gamma2, P, D]`, the first three, is the *analysis slice* the site packages
call the OCV proper — all six are tabulated and pinned in their channel tables.
Each package **recomputes** the committed
numbers of its site from a small input table, renders its plate and report, and
**pins every number against the committed JSONs**: a single deviation makes the run
exit non-zero. The claim of the repository is not "we assert", but "we reproduce
and read off".

## The twelve plates

| # | plate | what it answers | report |
| --- | --- | --- | --- |
| 1 | OCV schematic | what the vector *is*, drawn from the definition | [`figures/concept/fig_ocv_schematic.md`](figures/concept/fig_ocv_schematic.md) |
| 2 | one chip, three stages | what exactly is measured, and where | [`figures/concept/fig_ocv_measurement.md`](figures/concept/fig_ocv_measurement.md) |
| 3 | Espoo support chips | one site, two pass geometries, six dimensions | [`figures/espoo/fig_espoo_ocv_chips.md`](figures/espoo/fig_espoo_ocv_chips.md) |
| 4 | SDS reconfiguration | do two controlled experiments reconfigure the same way? | [`figures/sds/fig_sds_ocv_reconfigure.md`](figures/sds/fig_sds_ocv_reconfigure.md) |
| 5 | KDLO trajectory | where the vector moves across a rebuild | [`figures/kdlo/fig_kdlo_ocv_trajectory.md`](figures/kdlo/fig_kdlo_ocv_trajectory.md) |
| 6 | Carola trajectory | where it moves across a bridge collapse | [`figures/carola/fig_carola_ocv_trajectory.md`](figures/carola/fig_carola_ocv_trajectory.md) |
| 7 | Werder Havel trajectory | where it moves across documented works | [`figures/werder_havel/fig_werder_havel_ocv_trajectory.md`](figures/werder_havel/fig_werder_havel_ocv_trajectory.md) |
| 8 | LUMO orbit series | the six dimensions through four damage labels | [`figures/lumo/fig_lumo_ocv_orbits.md`](figures/lumo/fig_lumo_ocv_orbits.md) |
| 9 | Bautzen orbit series | the six dimensions through three deck states | [`figures/bautzen/fig_bautzen_ocv_orbits.md`](figures/bautzen/fig_bautzen_ocv_orbits.md) |
| 10 | Morandi trajectory | where the vector moves across a collapse | [`figures/morandi/fig_morandi_ocv_trajectory.md`](figures/morandi/fig_morandi_ocv_trajectory.md) |
| 11 | CTS trajectory | the same, on the two-layer DiD site | [`figures/cts/fig_cts_ocv_trajectory.md`](figures/cts/fig_cts_ocv_trajectory.md) |
| 12 | Haltern trajectory | the same, on the one-layer turbine record | [`figures/haltern/fig_haltern_ocv_trajectory.md`](figures/haltern/fig_haltern_ocv_trajectory.md) |

`figures/README.md` is the file-by-file index: stems, wrappers, containers, run
commands and the four rules every plate keeps (one hashed cut, every number pinned,
one measured page, no verdict).

## The nine site packages

| site | package | states | mask layer | chip the plate reads |
| --- | --- | --- | --- | --- |
| LUMO lattice tower (Hannover) | `code/lumo/` | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`) | 11 x 400 dwell strip |
| Bautzen OpenLabs bridge | `code/bautzen/` | RS, DS1, DS2 | deck-edge mask (geometrically anchored deck band) | 80 x 80 window (the band is inside it) |
| KDLO media tower (Garden City SD) | `code/kdlo/` | pre-collapse, during-rebuild (**epochs**) | echo mask recomputed from a committed 64-strip full-dwell cache | 35 x 475 support chip |
| Carolabrücke (Dresden) | `code/carola/` | pre-collapse (healthy), post-collapse | **two layers** over the committed 80 x 80 window: the echo mask of the window and the per-girder deck mask | 69 x 125 support chip |
| Ponte Morandi / Polcevera (Genova) | `code/morandi/` | pre-collapse (healthy), post-collapse | two layers as Carola, the two tracks in place of the girders | 51 x 137 support chip |
| Champlain Towers South (Surfside) | `code/cts/` | pre-collapse (healthy), post-collapse | two layers as Morandi, plus a **difference-in-differences** layer against two control towers | 39 x 145 support chip |
| Haltern AV9 wind turbine | `code/haltern/` | pre-collapse (healthy), post-collapse | **one layer**: the echo mask of the committed 7 x 7 window payload | 35 x 229 support chip |
| Espoo Kurttila mast | `code/espoo/` | **no event** — the control | the delivered mask size, copied and proved: `A` is never recomputed **in the committed table** — the chip plate reads the cut's own `A` | 35 x 129 support chip (75 chips, the two in-swath tracks) |
| Havelbrücke Werder (A 10) | `code/werder_havel/` | W2_bestand, W3_sanierung, W4_post | echo mask of the export's 80 x 80 complex window, one window per acquisition — the window **is** the chip (6 400 cells, no band cut; `support/` holds the *role* cut) | 80 x 80 window |

The last column is what each site's plate actually measures — the grid on which
`A`, `D`, `F` and `S` are counted, and so the grid the six dimensions describe.
Five of the six trajectory plates read a **support-sized chip**,
`data/<site>/support/<site>_ocv_channels.csv`: the band the asset's own geometry
supports, plus its 30 m guard and its 200 m reference ring (KDLO 35 x 475,
Carola 69 x 125, Morandi 51 x 137, CTS 39 x 145, Haltern 35 x 229 cells). Werder
is the exception, and the one place where `support/` does **not** mean a band:
the export cut one **80 x 80** complex window per acquisition at the site point,
so the window *is* the chip and `A` counts over all 6 400 cells, while
`support/` names the *role* cut of that same window (491 `target` of 991 rows).
Where a site also holds a **committed 80 x 80 layer** beside its chip (Carola,
Morandi, CTS), that layer is a second instrument: the two are reported side by
side and never mixed, and the trajectory plate measures the chip.

Each package README lists its scripts, its inputs, the numbers it pins and what its
numbers do *not* say. Six of them are the collapse, rebuild and works cases of the
trajectory plates; Espoo is the site without an event; LUMO and Bautzen are the two
records whose states are *controlled* rather than observed.

## What the two orbit plates show (Figs 8–9)

LUMO and Bautzen are the records that carry a **controlled** state column — four
damage labels at the tower, three deck states at the bridge — so their plates are
the two places where the vector is read *without* a collapse in the way.

**LUMO** draws one panel per dimension × orbit over 178 coherence rows joined to
482 acquisition rows, on the **actual local overpass instant** (`Europe/Berlin`,
joined by `burst_id`) and with the background bands taken from the 13 consecutive
`damage_label` runs rather than from fixed windows: the states are read where the
record puts them, and the two pass geometries are drawn side by side so that a
movement of one panel is only read against the other.

**Bautzen** draws the same twelve panels (ASC IW2, DESC IW3) over `RS` / `DS1` /
`DS2`. It deliberately leaves the control rectangle 500 m east out of the plate:
the control's ASC `gamma2` range (−0.94 … +0.93) is an order of magnitude wider
than the bridge's (−0.14 … +0.24) and would flatten the very signal the plate
exists to show. The control stays in the channel table, where the package's own
numbers read it.

Neither plate computes a pooled contrast: they are the *time-series* reading of the
table, and the pooled readings live in the site packages.

## What the cross-site plate shows (Fig 4)

The controlled-dose plate puts the two controlled experiments side by side — LUMO
at a 0 / 7.4 / 13 % frequency drop, Bautzen at `RS` → `DS1` → `DS2` — one row set
per experiment × channel, one column set per block, every panel on its own plane
and standardised against its own asset's reference. Its reading is a **null about
ordering**, and it is the honest version of the question "does the vector track the
dose?":

* `P1` — the tower's geometry magnitude rises with the dose — **did not hold**
  (ASC 0 < 0.581 < 0.34; DESC 0 < 0.419 < 0.225);
* `P2` — `DAM 3` and `DAM 4` share one dose, so their magnitudes should agree —
  **holds for geometry, partly for behaviour**: the same-dose gap is 0.036 (ASC) and
  0.108 (DESC) against cross-dose steps of 0.241 and 0.194 in geometry, but 0.204 /
  0.195 against 0.051 / 0.495 in behaviour;
* `P3` — the behaviour block need not follow the dose — **consistent**: 0 of 4
  behaviour tests reach Holm < 0.05.

The cross-asset comparison is descriptive: the same metric applied twice, not a
transferability test.

## What the chip plate shows (Fig 3)

Espoo is the site with **no event**, and its plate reads the support-sized chip
instead of a state contrast: 75 chips (ASCENDING 31 / DESCENDING 44, three shared
dates) over two tracks, standardised **per orbit** because the orbit geometry is
the only contrast the record carries. Two things are stated as facts rather than
caveats: `D = A / bbox_area` exactly and `gamma2_mask` carries `n = A` in its
denominator, so the mask *size* dominates both dimensions (`A` differs by a factor
7 between the orbits); and the plate is a **second instrument** beside the site's
pipeline export, agreeing with it only loosely (r(log A, log
`coherence_masked_pixels`) = −0.24, r(`gamma2_mask`, `coherence_gamma2`) = 0.13).

## What the six trajectory plates show (Figs 5–7, 10–12)

Every plate in this section reads its site's own table
(`data/<site>/support/<site>_ocv_channels.csv`): the support chip for the five
case sites, and the export's 80 x 80 window for Werder, whose site has no band.
None of them reads the committed 80 x 80 layer that Carola, Morandi and CTS also
carry — the nine-site table at the top says on which grid each plate counts, and
that grid is the chip.

Six sites, one template: the same twelve panels (`a1`–`a4` behaviour, `b1`–`b4`
geometry, `c1`–`c4` context), the same standardisation (`z = (x − median_ref) /
1.4826 · MAD_ref`, log for `gamma2, P, D, F, A`, identity for `S`, clipped at ±6,
block magnitudes as the RMS of their standardised dimensions), the same rule that a
month without a chip is a **hole** and never a line across it, and the same
possibility of a third **held-out** class — the healthy chips the chosen reference
window leaves out, drawn and measured, standardising nothing. The word printed at
the end is a function of two frozen persistence rules, never of the shape of the
curve.

| # | site | kept chips | reference | tail | what the rules print |
| --- | --- | --- | --- | --- | --- |
| 5 | KDLO media tower | 59 | `pre`, 16 | `rebuild`, 43 | resolves `P`, `A`, `S`, `R_B`, `R_G` |
| 6 | Carolabrücke, Dresden | 436 | `pre-collapse (healthy)`, 121 of 184 | `post-collapse`, 252 | resolves `gamma2`, `P`, `S`; held-out class: isolated unusual departure |
| 7 | Havelbrücke Werder, A 10 | 491 of 991 rows | `W2_bestand`, 119 | `W3_sanierung` 189 / `W4_post` 61 | sustained unusual geometry departure, confined to the intervention interval |
| 10 | Ponte Morandi, Genova | 207 | `pre-collapse (healthy)`, 198 | `post-collapse`, 9 | resolves `D`, `R_G` |
| 11 | Champlain Towers South | 178 | `pre-collapse (healthy)`, 136 | `post-collapse`, 42 | resolves `P`, `A`, `D` |
| 12 | Haltern AV9 turbine | 107 | `pre-collapse (healthy)`, 96 | `post-collapse`, 11 | resolves nothing |

Read as five statements about what a *timeline* can carry:

* **The record with the smallest reference is the one that resolves most.** KDLO
  has 16 reference chips against a 43-chip `rebuild` tail and resolves five
  quantities — but its tail is a **rebuild**, not a collapse aftermath, so the
  contrast answers a different question, and the plate's caption says so rather than
  leaving it to a footnote.
* **A chosen window buys a held-out class.** Carola restricts its reference to
  2022-09 … 2023-12 (121 of 184 healthy chips), so the other 63 are drawn and
  measured; the verdict is an *isolated unusual departure*, and the one directional
  statement the plate carries is that `z(P)` lies below the reference envelope in 5
  of 9 held-out months, with the joint `z(gamma2) > 0, z(P) < 0` pattern holding
  for four consecutive months — direction, not size, and not a verdict input.
* **Not every case is a collapse.** Werder Havel is a **documented works interval**
  (2024-07-01 → 2025-09-05) over the crossing's three strata, with the
  pre-intervention year held out, and it answers nine registered questions in
  order. `R_G` sits above the monthly reference q95 in 10 of the 15 interval months
  and for three consecutive months (2025-02 → 2025-04) — a *sustained* geometry
  departure, confined to the interval, with no behaviour departure to match.
* **A short tail is drawn, not dropped.** Morandi publishes 9 post chips against
  198 reference chips and reports the contrast with the interval it can carry —
  wide by construction, and said so.
* **An empty column is said out loud.** CTS has one orbit direction only, so its
  pooling is the whole record; Haltern has 13 of its 120 dates refused at the
  raster edge, all ascending, carried in the cut's own manifest rather than
  silently dropped. Haltern's verdict is the honest one for a one-layer record:
  **nothing resolves**.

None of the six claims damage, a failure precursor or a cause, and none compares
its own magnitudes with another's: two plates are comparable in *shape* and in
*where the movement sits*, not in size.

## The definition, drawn (Figs 1–2)

The first two plates have no site and no input table. `fig_ocv_schematic` carries
the reader from the physical asset, through the SAR observation geometry and the
asset-specific measurement support, to the radar echo representation and the six
components; `fig_ocv_measurement` takes **one** chip and draws it three times — the
native read chip, the measurement support, the echo mask — each stage admitting
less, ending at `OCV(t)`. Both are checked structurally instead of pinned to a
table (86 and 148 checks), and they are the definition the other ten plates are
read against.

## Folder structure

```
observability-core-vector/
├── code/          one package per site, plus the cross-site scripts
│   ├── <site>/    the nine records above — see code/<site>/README.md
│   ├── cases/     the trajectory plate (five sites, one registry) + its spec
│   ├── sds/       the cross-site reconfiguration report (the one script that
│   │              reads two sites) + its spec
│   ├── concept/   the definition, drawn — two plates, no input
│   ├── audit_ocv_units.py   the unit audit over every channel table
│   ├── figure_style.py      the one type ladder and the one print width
│   ├── audit_figure_style.py  the audit of both, over every figure
│   └── README.md
├── data/
│   ├── <site>/    the committed inputs of that site (payload or mask cache,
│   │              reference JSONs), the generated channel table and the JSON of
│   │              every plate
│   ├── ocv_channel_units.json, ocv_metadata_quarantine.json, …   the cross-site
│   │              audit files
│   ├── concept/   the build record of the two concept plates (no site)
│   ├── sds/       the cross-site numbers
│   └── README.md
├── figures/
│   ├── <site>/    that site's wrapper, report and plate
│   ├── sds/       the cross-site plate + its panel guide
│   ├── concept/   the definition plates (.jpg + vector .pdf + report)
│   ├── cases/     the driver that runs the five trajectory plates
│   └── README.md
├── .gitignore
└── README.md
```

Per site: `code/<site>/` reads `data/<site>/…` and writes `data/<site>/*.json` plus
`figures/<site>/{*.md,*.jpg}`. A plate migrated onto `code/figure_style.save()`
writes its drawing **three times under one stem** — the raster, a vector `.pdf` for
LaTeX and an `.eps`, which resolves the alpha PostScript cannot express
(`eps="blend"` composites each tint; `eps="raster"` takes the PDF through
`pdftops` where the alpha has to accumulate) — and the audit judges the width of
each container separately. All three are regenerated binaries and git-ignored.
`code/README.md`, `data/README.md`, `figures/README.md` and the site READMEs carry
the file-by-file lists.

## Quick start

```bash
cd research/observability-core-vector

# 1-2. the definition, drawn (no input at all)
bash figures/concept/fig_ocv_schematic.sh        # jpg + vector pdf + JSON + report
bash figures/concept/fig_ocv_measurement.sh      # jpg + vector pdf + JSON + report

# 3. the chip site — run through its generator, no .sh
python3 code/espoo/fig_espoo_ocv_chips.py        # plate + report + JSON
python3 code/espoo/fig_espoo_ocv_chips.py --verify

# 4. the cross-site reconfiguration plate — no .sh either
python3 code/sds/ocv_reconfigure.py              # plate + report + JSON + table
python3 code/sds/ocv_reconfigure.py --verify

# 5-7, 10-12. the six trajectory plates, one wrapper each
bash figures/kdlo/fig_kdlo_ocv_trajectory.sh
bash figures/carola/fig_carola_ocv_trajectory.sh
bash figures/werder_havel/fig_werder_havel_ocv_trajectory.sh
bash figures/morandi/fig_morandi_ocv_trajectory.sh
bash figures/cts/fig_cts_ocv_trajectory.sh
bash figures/haltern/fig_haltern_ocv_trajectory.sh

# 8-9. the two orbit plates
bash figures/lumo/fig_lumo_ocv_orbits.sh         # seconds, pins included
bash figures/bautzen/fig_bautzen_ocv_orbits.sh   # seconds, pins included

# the five registry sites in one run (nothing pooled across them)
bash figures/cases/fig_ocv_trajectory.sh         # a minute or two
bash figures/cases/fig_ocv_trajectory.sh --verify    # re-derive, write nothing
python3 code/cases/ocv_trajectory.py --selftest      # 221 checks, writes nothing

# the audits
python3 code/audit_figure_style.py               # the type ladder + the print width
python3 code/audit_ocv_units.py                  # what `A` counts, every chip

# if `python3` of the PATH is a virtualenv without numpy/matplotlib:
PYTHON=/usr/bin/python3 bash figures/lumo/fig_lumo_ocv_orbits.sh

# rebuild a channel table (each site README lists what its table reads)
python3 code/lumo/lumo_ocv_channels_csv.py
python3 code/carola/carola_ocv_channels_csv.py
```

## Determinism and requirements

Every generator is offline and deterministic: seeded bootstraps (`2000`
moving-block draws, block 3, seed 7 for the trajectory family), permutation nulls
with their own fixed seeds where a package needs them, no timestamps in any output
— so repeated runs produce identical JSON, Markdown and JPEG files, and `--verify`
can compare bytes-and-numbers against what is committed.

Python 3 with numpy, matplotlib and scipy (scipy only for the p-values of some
control tests); the trajectory generator itself runs on the standard library alone.
Only the KDLO package can talk to a database, and only when it is explicitly asked
to re-extract columns (`--fetch-db-extra`); no figure script ever does. Neither the
Carola package nor its plate touches the 182 MB payload file — only the one-time
`carola_ocv_masks.py --extract` does — and the Morandi package reads only its
committed window cache. Espoo has neither a cache nor a payload: its package reads
only `data/espoo/`, so it runs from a fresh clone.
