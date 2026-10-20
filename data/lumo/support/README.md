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
