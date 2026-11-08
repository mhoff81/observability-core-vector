# LUMO OCV channel table — support-sized chip, disjoint ring reference

This directory holds a second channel table for the LUMO lattice tower. It keeps
the committed pipeline and time base, and changes two coupled things: the chip is
sized from the **support** (the 11 x 400 dwell strip plus guard plus ring, in
metres, odd and anchor-centred) instead of the committed search window, and the
echo mask is no longer thresholded against the median *inside* that window but
against a ring that is disjoint from the support and from the support's own
guard band.

The directory is named for the **sizing method** (`support`), the rect file for
the **reference rule** (`lumo_rects_ring.txt`): both belong to the same cut, and
`tools/lumo_rects_s3.py` still writes that name.

The committed table (`../lumo_ocv_channels.csv`, generator
`code/lumo_ocv_channels_csv.py`) is untouched and still reproduces its own
figure exactly. Nothing here supersedes it; the two are meant to be read
side by side.

## The rule

```
mask = intensity >= 0.30 * peak(support)  AND  intensity >= 5.0 * median(ring)
       with at least 2 cells
```

* **support** — the 11 x 400 dwell strip (4 411 cells). The published chip is a
  *search* window, not a mast-centred patch: the eleven-row dwell is what
  smears the vibration response, so a compact support around the mast would
  measure the wrong thing.
* **guard** — 30 m, so the ring never touches the support.
* **ring** — 30-200 m from the support, ground distance, 35 128 / 39 486 /
  40 348 cells depending on the chip.

## The chips

| sub-swath | chip (rows x cols) | dates | support | guard | ring |
| --- | --- | --- | --- | --- | --- |
| IW1 | 431 x 107 = 46 117 | 60 | 4 411 | 6 578 | 35 128 |
| IW2 | 431 x 119 = 51 289 | 33 | 4 411 | 7 392 | 39 486 |
| IW2 | 431 x 121 = 52 151 | 85 | 4 411 | 7 392 | 40 348 |

`support + guard + ring == chip` holds element-exactly on all 178 dates, which
is what `tools/lumo_table.py` asserts per row.

## Why the reference moved

The committed threshold sits at `5 x median` of the window the measurement
itself fills, so the bar tracks the signal: the more the tower vibrates, the
higher the threshold. Two measurements decided the alternative:

* **Window sensitivity (L0).** Re-cutting the same date with a chip 8x and 16x
  larger moves `D` by 8x and `S` by 16x, while the published configuration is
  flat across states in `D` and `F`. A quantity that is flat under a 16x change
  of the window and jumps 8x under a state change is measuring the window, not
  the tower.
* **Peak != mast.** Template matching puts the mast at 52.380037 / 9.719982
  (ASC) and 52.379975 / 9.719994 (DESC), NCC 1.0000, element-exact on 3 of 3
  committed chips; the geolocation grid predicts the window centre to
  +-0.3 cells. But the matched peak wanders ~370 m across dates and is *not*
  the mast, so a mast-centred support (option A) was rejected. The ring keeps
  the support where the committed measurement put it and only moves the
  reference.

## What changed

Median per state, committed -> ring. `A/D/F/S` are **not comparable** between
the two sides: `A` is now the size of the coherent core inside a 431-line chip,
and `S` is measured from a support whose peak is no longer the window peak.

| orbit | dimension | healthy | DAM 3 | DAM 4 | DAM 6 |
| --- | --- | --- | --- | --- | --- |
| ASC | A | 12.5 -> 10 | 15.5 -> 12.5 | 10 -> 11 | 13 -> 5.5 |
| ASC | D | 0.0055 -> 0.0042 | 0.0059 -> 0.0045 | 0.0039 -> 0.0050 | 0.0064 -> 0.502 |
| ASC | F | 5.5 -> 4 | 9 -> 3 | 5 -> 4 | 7 -> 1.5 |
| ASC | S | 50.1 -> 110.9 | 67.6 -> 138.7 | 77.9 -> 108.2 | 57.8 -> 73.1 |
| ASC | P | 0.0155 -> 0.0935 | 0.0136 -> 0.0622 | 0.0130 -> 0.0791 | 0.0130 -> 0.1001 |
| ASC | gamma2 | 0.0673 -> 0.1241 | 0.0621 -> 0.1289 | 0.0673 -> 0.0683 | 0.0467 -> 0.3738 |
| DESC | A | 16 -> 4 | 10 -> 4 | 16 -> 4 | 12 -> 4 |
| DESC | D | 0.0275 -> 0.0102 | 0.0149 -> 0.0107 | 0.0108 -> 1.0 | 0.0297 -> 0.0071 |
| DESC | F | 7 -> 2 | 5 -> 2 | 7 -> 1 | 6 -> 2 |
| DESC | S | 33.8 -> 14.2 | 28.2 -> 13.3 | 41.4 -> 0.7 | 42.9 -> 28.5 |
| DESC | P | 0.0148 -> 0.1917 | 0.0140 -> 0.1333 | 0.0150 -> 0.2400 | 0.0170 -> 0.1869 |
| DESC | gamma2 | 0.0340 -> 0.6509 | 0.0439 -> 0.6086 | 0.0279 -> 0.7066 | 0.0267 -> 0.4556 |

Three consequences worth stating plainly:

1. **gamma2 rises by 10-20x** (0.03-0.07 -> 0.46-0.71). The committed value
   was dominated by the median of the strip the smearing sits on; with a
   disjoint reference the coherent core is what is left.
2. **The state ordering disappears.** Committed DESC runs 0.034 -> 0.044 ->
   0.028 -> 0.027; under the ring it runs 0.651 -> 0.609 -> 0.707 -> 0.456, no
   monotone DAM trend. Any argument that leans on the low committed values is
   leaning on the self-referential threshold.
3. **The mask collapses to the minimum.** DESC `A = 4` in every state, and
   `D = 1.0` at DAM 4 (ASC `D = 0.50`): a 2x2 cluster fills its own bounding box
   exactly. That is the `min_n_masked = 2` floor, i.e. the mask now measures a
   core, not an extent.

## The pooled SDS of the whole tower

The number the paper quotes for this site is a **pooled** one: `healthy` against
all three damaged campaigns at once, with both orbits in the same sample. It is
also the aggregation the committed EMV reference used (`sds_pooled_all_damaged`),
which is why it is reproduced here and not replaced. `gamma2` is the whole-tower
coherence of this record — one scatterer, so the record splits into no segments —
and pooling ASC with DESC therefore pools the whole tower, not a part of it.

```
SDS      = 10 * |Cliff's delta|
delta    = (#{y > x} - #{y < x}) / (n1 * n2)      ties count for neither side
95 % CI  = percentile bootstrap, both groups resampled independently,
           10 000 resamples, seed 7 — one interval per comparison
```

Committed window (left) and this ring chip (right), healthy versus all damaged:

| dimension | committed SDS | delta | 95 % CI | ring SDS | delta | 95 % CI |
| --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 0.10 | +0.010 | [-0.165, +0.185] | 0.26 | +0.026 | [-0.157, +0.212] |
| `P` | 1.05 | -0.105 | [-0.292, +0.089] | 0.54 | -0.054 | [-0.237, +0.128] |
| `D` | 0.86 | -0.086 | [-0.272, +0.098] | 0.20 | +0.020 | [-0.159, +0.198] |
| `F` | 0.67 | -0.067 | [-0.244, +0.107] | 0.51 | -0.051 | [-0.219, +0.119] |
| `A` | 0.81 | -0.081 | [-0.255, +0.090] | 0.40 | -0.040 | [-0.215, +0.138] |
| `S` | 0.10 | -0.010 | [-0.196, +0.172] | 0.22 | -0.022 | [-0.205, +0.162] |

n = 121 healthy / 57 damaged on both sides, 178 overpasses. `P` is missing for
the healthy 2021-07-14 ASC overpass in the committed table, so that cell rests on
120 / 57. **Every interval includes zero**, on both chips.

Three readings:

1. **The pooled score separates nothing, and the ring does not change that.**
   The whole-tower coherence moves from SDS 0.10 to 0.26 and its interval covers
   zero either way. Pooling two observation geometries was always the
   compatibility pin and the context, never the primary evidence; that is now
   true of the ring chip as well.
2. **The mask geometry flattens towards its floor.** `P` 1.05 -> 0.54, `D` 0.86 ->
   0.20, `A` 0.81 -> 0.40, and the medians say why: `A` 15 -> 12 cells commits to
   `A` 4 -> 4 cells, `D` 0.0128 -> 0.0088 to 0.0071 -> 0.0077. With
   `min_n_masked = 2` the pooled mask is a core in every state, so pooling the
   states leaves nothing to rank — the opposite of a threshold that tracks the
   signal. Pooled is doing work in that sentence: the medians are weighted
   towards DESC, and on ASC alone the mask keeps 10-11 cells (see *By orbit*).
3. **Only `gamma2` carries its meaning across the change.** `A/D/F/S` are not
   comparable between the two sides and `P` cannot cross chips at all (it is a
   per-chip-shape pixel-frequency map), so the column pair is a comparison of
   *rules*, not of two measurements of the same quantity.

The committed column is pinned, which is what keeps this section a measurement
instead of a retyping: `gamma2` = 0.1029433 is the pooled score the paper's
reference transcribes (Table 4 and section 4.1: "the pooled
Healthy-versus-all-damaged comparison yielded SDS = 0.10"), and `A D F S P`
reproduce `sds_pooled_all_damaged` of
`data/lumo/reference/lumo_echo_mask_vector.json` to 1e-9. The ring column is the
same number `data/sds/sds_support.md` section 3 reports for `lumo / support`,
counts included.

```sh
python3 code/sds/sds_support.py --site lumo --out /tmp/sds_lumo
```

recomputes all four levels from the two tables alone. Section 3 of
`/tmp/sds_lumo/sds_support.md` (*Compatibility pin and context*) carries the same
twelve numbers, with `n1 / n2` and the CI flags spelled out per row.

### By orbit

The same contrast one orbit at a time — the level the cross-site report calls
`per_channel_pooled_damage`, and exactly the aggregation its reading rules keep
out of the primary number. Ring chip:

| dimension | ASC SDS | delta | 95 % CI | DESC SDS | delta | 95 % CI |
| --- | --- | --- | --- | --- | --- | --- |
| `gamma2` | 0.62 | +0.062 | [-0.276, +0.400] | 0.25 | -0.025 | [-0.252, +0.199] |
| `P` | 1.29 | -0.129 | [-0.453, +0.188] | 0.61 | -0.061 | [-0.283, +0.169] |
| `D` | 1.93 | +0.193 | [-0.146, +0.512] | 0.32 | -0.032 | [-0.243, +0.180] |
| `F` | 1.71 | -0.171 | [-0.472, +0.140] | 0.17 | +0.017 | [-0.178, +0.215] |
| `A` | 0.38 | -0.038 | [-0.363, +0.300] | 0.45 | +0.045 | [-0.163, +0.259] |
| `S` | 0.59 | +0.059 | [-0.279, +0.403] | 0.10 | +0.010 | [-0.212, +0.236] |

n = 40 healthy / 17 damaged ascending, 81 / 40 descending — the two orbits
partition the 121 / 57 above. **All twelve intervals include zero.** Splitting
the orbit halves the sample (40 / 17 on ASC) and buys width, not a separation.

Three things the split shows that the pooled number hides:

1. **The point estimates differ by geometry, without supported directional
   effects.** `gamma2` (+0.062 vs -0.025), `D` (+0.193 vs -0.032), `F` (-0.171
   vs +0.017) and `A` (-0.038 vs +0.045) have opposite delta signs on the two
   orbits. Every one of those intervals is wide and centred near zero on both
   sides, so the difference is a property of the geometry each orbit samples:
   neither orbit supports "damaged differs", in either direction.
2. **The collapse to the `min_n_masked` floor is predominantly DESC-driven.**
   Ring `A` keeps a median of 10 cells healthy / 11 damaged ascending against
   4 / 4 descending, and `D` 0.0042 / 0.0047 against 0.0102 / 0.0106. The
   support-ring echo population is materially larger on ASC: mean `A` 10.3
   against 4.6, and 19 % of ASC rows sit at 5 cells or fewer against 81 % of
   DESC rows. DESC supplies 121 of the 178 rows, so the pooled median (4, mean
   6.4) falls on the DESC value — the collapse quoted under *What changed* is
   mostly this population, not a property of the chip as such.
3. **The one non-straddling orbit result of the whole set does not survive this
   chip.** On the committed window, ascending `P` separates: SDS 4.98, delta
   -0.498, CI [-0.768, -0.192], n 39 / 17 (39 because of the missing 2021-07-14
   ASC overpass). On the ring the same contrast gives 1.29 and an interval that
   includes zero. It is the only orbit-specific interval in either table that
   excludes zero, out of 24, and the ring drops it — which is context, not a
   separation: one cell, on the chip that is being replaced.

The committed counterpart of the two columns, for scale only:

| dimension | committed ASC (40 / 17) | committed DESC (81 / 40) |
| --- | --- | --- |
| `gamma2` | 1.35, -0.135, [-0.435, +0.176] | 0.85, +0.085, [-0.123, +0.294] |
| `P` | 4.98, -0.498, [-0.768, -0.192] (39 / 17) | 0.79, +0.079, [-0.148, +0.302] |
| `D` | 0.18, -0.018, [-0.332, +0.300] | 1.11, -0.111, [-0.338, +0.121] |
| `F` | 2.15, +0.215, [-0.101, +0.512] | 1.56, -0.156, [-0.367, +0.060] |
| `A` | 0.85, +0.085, [-0.232, +0.396] | 1.40, -0.140, [-0.348, +0.070] |
| `S` | 0.21, -0.021, [-0.332, +0.303] | 0.12, -0.012, [-0.241, +0.214] |

Both columns are the `per_channel_pooled_damage` level of the committed report,
`data/sds/sds_support.json` →
`sites.lumo.levels.{committed,support}.per_channel_pooled_damage`, produced by
`code/sds/sds_support.py` (`level_per_channel_pooled_damage`, line 315). The
level is carried in the JSON only — `data/sds/sds_support.md` prints per state
(section 1), the maximum (section 2) and the pooled pair (section 3).

## Verification

`fig_lumo_ocv_orbits.md` / `.json` / `.jpg` are produced by the committed
figure script, pointed at this table:

```sh
python3 code/lumo/fig_lumo_ocv_orbits.py \
    --csv data/lumo/support/lumo_ocv_channels.csv \
    --meta data/lumo/support/lumo_ocv_channels_meta.json
```

**123 cross-checks, 0 failures.** Row, date and burst counts, the 482-row
acquisition join, the 13 runs, the sub-swath counts, `A == n_masked`,
`A <= chip`, strictly positive log channels and the local-time offsets all hold
exactly as for the committed table.

**24 pins are listed as `superseded`, not failed.** They describe the committed
rule — the median inside the window, the 400 x 11 chip, channel completeness,
and the `n_masked` / `gamma2` values the paper reported — so a ring table
cannot satisfy them. Both sides are kept in `superseded` in the JSON and in
section 5 of the markdown instead of being deleted, so the size of the change
stays visible.

## Reproduce from scratch

```sh
cd /home/projects/erlking_mvp
# anchor: exact template match against two committed chips, plus the grid check
.venv/bin/python3 tools/lumo_anchor.py --match <chip> --second <chip> \
    --grid <grid> --json /tmp/lumo_anchor.json
# 178 chips: dwell support, 30 m guard, 30..200 m ring, one chip per sub-swath
.venv/bin/python3 tools/lumo_rects_s3.py --cut --support dwell \
    --guard-m 30 --ring-m 200 \
    --out <repo>/data/lumo/support --json /tmp/lumo_cut.json
# the table: A/D/F/S against the ring, P per chip shape, meta with the pins
.venv/bin/python3 tools/lumo_table.py \
    --rects <repo>/data/lumo/support/lumo_rects_ring.txt \
    --report /tmp/lumo_cut.json \
    --out <repo>/data/lumo/support/lumo_ocv_channels.csv
cd <repo> && python3 code/lumo/fig_lumo_ocv_orbits.py \
    --csv data/lumo/support/lumo_ocv_channels.csv \
    --meta data/lumo/support/lumo_ocv_channels_meta.json
```

`lumo_rects_ring.txt` (82 MB) is the full chip series; the manifest beside it
records the point, the chips, the support/guard/ring sizes and that 0 of 178
dates were refused. `lumo_table.py` reads the non-mask columns from the
committed `../lumo_ocv_channels.csv` and rewrites only what the ring changes.
