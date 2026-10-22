# `code/kdlo/` — KDLO site package: analysis and plotting scripts

Eight Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values of the control
tests; without it the scripts still run, just with fewer p-values). This is the
KDLO half of the repository; the LUMO tower site lives in `../lumo/` and Bautzen
in `../bautzen/`.

The mask layer follows the LUMO rule (`0.30 x peak`, `5 x median`) and carries the
same 6D echo-mask vector `x_KDLO = [gamma2, P, D, A, F, S]`, but the KDLO backend
stored only the mask *size* (`A = coherence_masked_pixels`). The masks are
therefore recomputed from a **committed** cache of the 64 full-dwell strips the
pipeline decoded (`data/kdlo/strips/`, 11 x 400 complex samples each, ~70 KB per
file), produced by `tower/backend/src/cli/analysis/kdlo_strip_capture.rs` and read
here by `kdlo_ocv_masks.py`. The OCV proper is `[gamma2, P, D]`; `A` and the
masked coherence are pinned to the pipeline's own `coherence_masked_pixels` and
`coherence_gamma2` on every acquisition.

All paths below are relative to the **repository root** and the scripts are always
started from there (the wrapper `../../figures/kdlo/fig_kdlo_ocv_paper.sh` does
that): inputs and outputs in `data/kdlo/`, figures in `figures/kdlo/`.

### Two chips, two figure sets (additive, never replacing)

KDLO exists twice in this package, and the difference is the **search window**, not
the record:

| chip | table | cut by |
| --- | --- | --- |
| `committed` | `data/kdlo/kdlo_ocv_channels.csv` | `kdlo_ocv_channels_csv.py`, the fixed window of the original pipeline |
| `support` | `data/kdlo/support/kdlo_ocv_channels.csv` | `erlking_mvp/tools/ocv_migrate.py`, the support-sized chip (guard 30 m, reference ring 200 m, chip 35 x 475 cells, ring complete) |

The support table carries the same 60 acquisitions and the committed values beside
the recomputed ones (`A_committed` … `gamma2_committed`), so the two can be compared
row by row — 300 of 300 committed values are carried exactly.

The three figure packages take `--table` and default to the committed chip; the
support run appends `_support` to every output stem, so nothing is overwritten:

```sh
python3 code/kdlo/fig_kdlo_ocv_paper.py                 # committed: 3035 pins
python3 code/kdlo/fig_kdlo_ocv_paper.py --table data/kdlo/support/kdlo_ocv_channels.csv
python3 code/kdlo/fig_kdlo_ocv_months.py --table data/kdlo/support/kdlo_ocv_channels.csv
python3 code/kdlo/fig_kdlo_ocv_series.py --table data/kdlo/support/kdlo_ocv_channels.csv
```

On the support chip the pin block is a **different, smaller set** (`pin_support`):
the chip itself (`generated_by`, rows/cols, `ring_complete`, guard, ring), the same
60 ids, the carry-over audit, and a findings lock
(`data/kdlo/reference/kdlo_support_findings.json`, written once with
`--write-reference`). The 123 committed-value comparisons stay with the committed
run, where they belong. The series' record-side self-checks are skipped on the
support chip for the same reason, and the run says so.

**Why this matters.** Same acquisitions, same code, same estimators — only the
window differs, and the effect sizes move by more than a rounding:

| channel | Δ committed | Δ support | sign flip | SDS |
| --- | ---: | ---: | --- | --- |
| `gamma2` | +0.529 | +0.529 | no | 5.29 → 5.29 |
| `P` | +0.235 | −0.628 | **yes** | 2.35 → **6.28** |
| `S` | −0.145 | −0.640 | no | 1.45 → **6.40** |
| `A` | −0.281 | +0.499 | **yes** | 2.81 → **4.99** |
| `F` | −0.267 | +0.121 | **yes** | 2.67 → 1.21 |
| `D` | −0.070 | −0.201 | no | 0.70 → 2.01 |

`gamma2` is the pipeline coherence and does not depend on the chip, which is why it
is the one row that does not move. On the committed chip it is the largest effect;
on the support chip **`S` and `P` overtake it** and three of the five mask dimensions
flip sign. The support report therefore names its headline from the table instead of
naming `gamma2` in advance.

### `kdlo_ocv_stats.py` — statistics layer (verbatim port)

`_stats`, `welch_mw_test`, `cliffs_delta` (+ `cliffs_delta_bruteforce`),
`bootstrap_delta_ci`, `delta_block`, `sign_test`, `sds_from`, `by_state`, `vals`,
`season_of`, `strata_pair`, `fit_lda`, `predict_lda`, `loo_cv`,
`permutation_test`, `fit_linear_score`, `score_row`, `oof_scores`, `delta_sds`,
`permutation_p`, `paired_bootstrap_diff`. Copied line by line from
`../lumo/lumo_ocv_stats.py`, which is itself a verbatim port of the LUMO project
analysis scripts, so a difference in the numbers cannot come from a difference in
the statistics. The only KDLO edits are naming: the state field is `state`
(LUMO: `damage_label`) and the two states are **epochs** (`pre`, `rebuild`)
instead of damage labels.

### `kdlo_ocv_core.py` — data layer, constants and the committed aggregates

Paths, the analysis channels (`OCV = [gamma2, P, D]`), the 6D feature vector
(`FEATURES = [gamma2, P, D, A, F, S]`, `DIMS`, `DIM_LABEL`) and the six controls
(`CONTROLS`, `WEATHER`), the epoch windows (`STATE_WINDOWS`, `SINGLETON_DAY`), the
pipeline thresholds (`COMPACT_MASK_PX = 10`, `MASK_CLUTTER_PX = 25`), the
model/seed constants (`MODELS`, `HEADLINE_PAIR`, `LAMBDA_GRID`,
`HEADLINE_LAMBDA = 0.5`, `RNG_SEED = 7`, `N_BOOT = 10000`) and the loaders
(`load_csv`, `load_csv_meta`, `load_raw_parts`, `load_db_extra`, `load_reference`).
`MASK_COLUMNS` lists the echo-mask columns the generator appends, so that
`CSV_COLUMNS` (and the CSV header) carry the whole 6D vector.

It also carries the two aggregate blocks the committed channel JSONs were written
with — `summ`, `row_brief`, `echo_modes_of`, `phase_observable_of`,
`part_stats` — and the two pin-only mirrors of the files that cannot be
recomputed from the record: `availability_checks` (188 structural checks of
`kdlo_s1_availability.json` — the record's two windows plus the 2015-2021
ascending IW2 era as `figure_0` and the geometry break that separates it) and
`fem_checks` (67 checks of `kdlo_fem_frequencies.json`).

### `kdlo_ocv_masks.py` — the echo-mask layer (port of `../lumo/lumo_ocv_masks.py`)

`read_complex_bin`, `connected_components_8`, `vector_from_strip` (returns the mask
plus `A`, `D = A / bbox_area`, `F`, `S`, the bbox/centroid/peak geometry and the
masked `gamma2`), `add_persistence` (`P` = mean pixel frequency of the own mask
over all dates of the majority strip shape), `load_vectors` and `check`.

The one deliberate difference from the LUMO port is the background median of the
mask rule: it is the **upper** median `sorted(v)[len // 2]`, the value the Rust
pipeline uses (`sarsolve::coherence::coherence_mast_echo`), because `A` has to
reproduce `coherence_masked_pixels` bit for bit. `np.median` is kept as the
diagnostic `A_np_median` / `gamma2_np_median` pair and never reported.

```bash
python3 code/kdlo/kdlo_ocv_masks.py           # acceptance test on all strips
python3 code/kdlo/kdlo_ocv_masks.py --json    # same, plus the derived vectors
```

The check recomputes the mask of every committed strip and requires
`A == coherence_masked_pixels`, `abs(gamma2_mask - coherence_gamma2) <= 1e-12` and
a 400 x 11 strip shape. It is a self-contained acceptance test of the mask layer
that needs neither database nor network — the per-strip sha256 values of
`strips/manifest.json` are the pins behind it.

### `kdlo_ocv_channels_csv.py` — generator of the acquisition table

Merges the two committed raw API dumps (`data/kdlo/fig_kdlo_1_raw.json`,
`fig_kdlo_2_raw.json`) with the **one-time, read-only** Postgres extraction of
nine columns that the API rows do not contain (`kdlo_db_extra.json`: phase SNR,
phase observability, snow, precipitation, gusts, intensity, incidence angle),
derives `part`, `state`, `month`, `year`, `season` and `echo_mode`, attaches the
echo-mask columns of the 6D vector (`A`, `D`, `F`, `S`, `P`, the mask geometry and
the recomputed `gamma2_mask`) from the committed strip cache via
`kdlo_ocv_masks.attach_masks`, and writes `data/kdlo/kdlo_ocv_channels.csv`
(60 rows, one per acquisition, 58 columns) plus `kdlo_ocv_channels_meta.json`
(provenance: file sha256s, the windows, the definitions, the counts, and the
`strip_cache` block with the mask rule and the manifest sha256).

The build **fails** if a row has no echo mask or if any `A`/`gamma2_mask`
deviates from the record, so the CSV cannot be written with a lookalike mask. The
verification repeats that comparison and additionally re-hashes every committed
strip (`--no-strip-verify` skips only the hashing).

The **2015-2021 ascending IW2 era** (`part0`: 162 acquisition days, rel. orbit 136,
~00:2x UTC) is *not* a row source here. Its echo-mask geometry would have to be
recomputed from full-dwell strips and `data/kdlo/strips/` holds only the 60
acquisitions of parts 1/2, so an era row would carry no `P`/`D`/`F`/`S`. Its
decoded rows are exported instead (`kdlo_export_raw.py --part 0` in the KDLO
project) to `data/kdlo/fig_kdlo_0_raw.json` — committed *evidence*, whose sha256
the availability reference pins. The **generator never reads it** (an era row
cannot enter the table: it would carry no `P`/`D`/`F`/`S`), and the era-to-collapse
series companion `fig_kdlo_ocv_series.py` is its only reader: it re-verifies that
sha256 before using a single value and draws only the two era columns that exist
(`gamma2`, `A`). The era is therefore present as catalogue facts (`figure_0`, the
162 days, the 69-day blind gap), as backend accounting and as the series figure's
era stretch — never as pooled rows.

```bash
python3 code/kdlo/kdlo_ocv_channels_csv.py                 # build + verify
python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only   # verify only
python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only --no-strip-verify
# one-time re-extraction of the nine DB columns (needs the `tower-postgres`
# container; the result is vendored, so this is not part of a normal run):
python3 code/kdlo/kdlo_ocv_channels_csv.py --fetch-db-extra
```

### `kdlo_capture_requests.py` — the request list of the strip cache

Reads the 64 measurement rows of `onboarder.insar_measurements` for the KDLO asset
(read-only, `docker exec psql`) and writes `data/kdlo/strips/capture_requests.json`:
site point, strip spec and, per acquisition, the CDSE burst Id, orbit direction,
relative orbit, sub-swath, polarisation, product name and content start. The Rust
binary turns that list into the strips; this script is what makes the capture
reproducible from the record.

### `fig_kdlo_ocv_paper.py` — figures A–E, JSON result, report, pin block

The main script. `compute()` recalculates every number of the five figures from
`data/kdlo/kdlo_ocv_channels.csv`; `pin_all()` then compares each number with the
four committed reference JSONs in `data/kdlo/reference/`, row by row (3045 checks
in a full run, 3035 in `--quick`; 188 + 67 of them are the structural mirrors of
the two pin-only files), verifies the echo-mask layer (the `strip_cache` block of
the metadata, every committed strip's sha256 through the manifest sha256, and per
row `A == coherence_masked_pixels`, `abs(gamma2_mask - coherence_gamma2) <= 1e-12`,
`D == A / bbox_area`, `P` in [0, 1]). A single deviation ends the run with exit
code 1. `main()` writes the five figures to `figures/kdlo/`, the result JSON to
`data/kdlo/fig_kdlo_ocv_paper.json` and the English report to
`figures/kdlo/fig_kdlo_ocv_paper.md`.

```bash
python3 code/kdlo/fig_kdlo_ocv_paper.py            # full run, ~1 min
python3 code/kdlo/fig_kdlo_ocv_paper.py --quick    # without the permutation null
python3 code/kdlo/fig_kdlo_ocv_paper.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

Wall-clock times are only printed to stdout, never written into an artifact: with
the fixed seeds (RNG seed 7, 10000 bootstrap draws, 2000 for the strata, 1000
label permutations) two runs produce identical JSON, Markdown and PNG files.

### `fig_kdlo_ocv_months.py` — the echo-mask dimensions month by month over the record

The **standalone companion** of the paper figure, run by
`figures/kdlo/fig_kdlo_ocv_months.sh` from the repository root. Where the A–E
figures pool the two **epochs** (`pre-collapse` vs. `during-rebuild`), this script
keeps the *time axis*, so the reader sees **when** the mask moved, not only **that**
it moved. It reads only `data/kdlo/kdlo_ocv_channels.csv`, recomputes each
calendar month's block with `kdlo_ocv_stats`, draws one figure and writes a small
JSON + English report. It adds **no** number to the A–E pin contract and touches
none of it (it never reads `data/kdlo/reference/`).

Why the axis is **calendar months over the whole record** and not the 24 rolling
windows of the other sites: the KDLO record is 60 acquisitions in two **disjoint**
epochs (`pre` 2022-06-02..2022-12-11 16 acquisitions, `rebuild`
2023-05-04..2024-09-19 43 acquisitions), so there is no dense 24-month pre-stretch
to roll; the axis is the whole record by **calendar month** (28 buckets, 4 of them
empty) and the acquisition gap is drawn, not filled. The CSV carries no `day`
column — the month is taken from `acquisition_ts` through `core.day_of`.

```
pre     2022-06-02..2022-12-11   standing mast      16 acquisitions
single  2022-12-23               in neither epoch    1 acquisition
rebuild 2023-05-04..2024-09-19   during the rebuild 43 acquisitions
```

Two properties of the record are kept visible rather than smoothed away. The
collapse month **2022-12 is mixed** (it holds the last `pre` acquisitions, to
2022-12-11, and the `single` 2022-12-23 that belongs to neither epoch), so the
collapse cut is drawn *inside* that bucket (a red dot and label at a fractional
x) instead of on a boundary; and the **4-month acquisition gap** between the
epochs (`2023-01 … 2023-04`, no acquisitions) stays on the axis as four empty
buckets with a dotted gap line. Panel (f) is therefore read against the first
record month (`2022-06`) as the 1.0 line, and the epoch tints (`pre` blue,
`rebuild` orange) mark which side of the cut a month belongs to.

The figure is 2 x 3: (a)–(d) one dot per acquisition for `gamma2`, `P`, `D`, `A`
with the month median and IQR (`D` on a symlog y axis, `linthresh=0.01`, so the
near-zero cluster stays visible); (e) the `(P, D)` plane, one dot per acquisition
coloured by month; (f) the 6D fingerprint — each month's median divided by the
`2022-06` baseline median on a log axis. `main()` writes the PNG to
`figures/kdlo/fig_kdlo_ocv_months.jpg`, the JSON to
`data/kdlo/fig_kdlo_ocv_months.json` and the report to
`figures/kdlo/fig_kdlo_ocv_months.md`:

```bash
python3 code/kdlo/fig_kdlo_ocv_months.py            # ~11 s
bash figures/kdlo/fig_kdlo_ocv_months.sh            # the same, via the wrapper
python3 code/kdlo/fig_kdlo_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

It is deterministic (fixed seeds, no timestamps), so two runs are byte-identical,
and it is **descriptive**: the 6 pure `pre` months carry 16 acquisitions (two to
three per month), the collapse month 2, so every interval is very wide and the
month axis is read as the record's own co-variation over two disjoint windows —
the reconstruction, the season, the acquisition geometry and the pre/post mix move
with the collapse date as well — **never** as damage.

### `fig_kdlo_ocv_series.py` — the acquisition series, from the first scene to the collapse

The **second standalone companion**, and the only script of the package that reads
more than the record's own stretch. It puts the record back into the seven years
that precede it and draws the echo-mask vector on **one dated axis**,
`2015-03-20 … 2022-12-11` (the collapse cut). Its four committed inputs are the
era export `data/kdlo/fig_kdlo_0_raw.json` (committed evidence; its sha256 is
re-verified against the availability reference before a single value is used), the
era's **own** strip cache `data/kdlo/strips_era/` (manifest + capture requests,
162 full-dwell strips, hashed and re-decoded here), the record's
`kdlo_ocv_channels.csv` and
`reference/kdlo_s1_availability.json`, from which the four stretches of the axis
are **read, not typed**: the ascending IW2 era (162 days), the blind gap, the
catalogue-only descending stretch and the record's 16 pre-collapse days.

All **six** dimensions exist on both sides of the break. `gamma2` and the mask
size `A` are the same columns by name (`coherence_gamma2`,
`coherence_masked_pixels`, which the era's own export carries; the mask size comes
out of the era strips exactly and the coherence to `1e-12`). `P`, `D`, `F` and `S`
are recomputed from committed full-dwell strips with the pipeline's own mask rule:
from the strip set the channel table uses for the record, and from a **second**
cache `data/kdlo/strips_era/` for the era, one strip per era acquisition day. The
era cache is never taken on its own say-so — the mask rule recomputed on those 162
strips has to reproduce the export's two columns on every era day, the same
acceptance `kdlo_ocv_masks.py` runs on the record's cache. `P` is a mask
*frequency*, so it is only defined against the set it is counted over: the era's
over its own 162 strips (a mixed-polarisation set — that track carries both VV and
VH, one burst per day), the record's over the channel table's set. The two
frequency maps are never added, averaged or used to rescale one another. Each
stretch is drawn in units of **its own** median (1.0 drawn, log where a dimension
spans decades), because the two stretches sit on different tracks, swaths and
times of day: what each stretch *does* is the comparable quantity, not the level
it sits at. No statistic is computed across the break.

The figure is 3 x 2 (183 mm wide): (a) `gamma2` and (b) `A`, era dots (blue) vs.
record dots + line (orange) on a log y axis, with the era's per-year medians as
bars in all six panels; (c)-(f) `P`, `D`, `F`, `S`, drawn the same way on both
stretches. Every panel carries the blind gap,
the catalogue-only stretch, the collapse cut and the first collapsed acquisition
(`2022-12-23`) drawn hollow and labelled as outside the axis. `main()` writes the
PNG to `figures/kdlo/fig_kdlo_ocv_series.jpg`, the JSON to
`data/kdlo/fig_kdlo_ocv_series.json` and the report to
`figures/kdlo/fig_kdlo_ocv_series.md`:

```bash
python3 code/kdlo/fig_kdlo_ocv_series.py            # <1 s
bash figures/kdlo/fig_kdlo_ocv_series.sh            # the same, via the wrapper
python3 code/kdlo/fig_kdlo_ocv_series.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

The era's own clock is **re-derived rather than asserted**: the reference's
per-year table labels the hour block `00`, the export times its rows on the
pipeline's `05:00Z` `morning` slot, and the one metadata-anomaly day (`2020-08-13`)
sits on the `afternoon` slot at `16:00Z` because its catalogue burst name carries
the wrong hour. All of that is self-checked and none of it enters a drawn value.
The figure carries **125 self-checks** that recompute every drawn number from the
four committed inputs — the era strips among them, re-decoded and held to the
export — and abort before any file is written on drift. Like the
months companion it adds **no** number to the A–E pin contract and touches none of
it, and it is byte-deterministic (no timestamps, fixed inputs).


The module and site names keep their `kdlo_` prefix although the folder already
says `kdlo/`: the import statements (`import kdlo_ocv_core as core`), the
reference file names and every command line in the reports stay parallel to the
other two sites, which is what makes the three packages readable as one diff.
