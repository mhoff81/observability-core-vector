# The support-sized migration — six sites, one chip model

Every site's chip is now sized from its **asset** and its **channel** instead of
the committed search window:

```
asset polygon (WGS84)  ->  projected base ring + top ring (height shift)
                       ->  uncertainty budget -> margin
                       ->  support band (odd, in cells)
                       ->  chip = support + guard (30 m) + ring (30..200 m)
                       ->  mask = peak(support) and 5 x median(ring)
```

The asset point is the **polygon's centroid**, never a typed coordinate: the
registry used to carry hand-entered points and two of them were wrong by 330 m
(CTS) and ~400 m (Carola), i.e. by more than half the chip's azimuth extent.

The committed tables are untouched. Each migrated site gets
`data/<site>/support/` with three artefacts:

| file | content |
| --- | --- |
| `<site>_rects_support.txt` | the chips, one line per `(date, orbit)` |
| `…txt.manifest.json` | the cut: chip, band, support/guard/ring cells, margins, height shift, assumptions, anchors with residuals, every refusal with its reason, and the burst id per row |
| `<site>_ocv_channels.csv` + `_meta.json` | the channel table (all non-mask columns copied from the committed table, `A D F S P gamma2` recomputed, `*_committed` beside each) |

## Status

| site | band (lines x bins) | chip | support / guard / ring cells | cut | table rows | refusals |
| --- | --- | --- | --- | --- | --- | --- |
| CTS | 9 x 25 | 39 x 145 | 225 / 420 / 5 010 | 178 / 178 | 178 | 0 |
| KDLO | 5 x 379 | 35 x 475 | 1 895 / 2 450 / 12 280 | 60 / 60 | 60 | 0 |
| Carolabrücke | 39 x 15 | 69 x 125 | 585 / 900 / 7 140 | 436 / 436 | 436 | 0 |
| Haltern | 5 x 109 | 35 x 229 | 545 / 852 / 6 618 | 107 / 120 | 107 | 13 |
| Ponte Morandi | 21 x 41 | 51 x 137 | 861 / 678 / 5 448 | 207 / 232 | 207 | 25 (alle ASC) |

`support + guard + ring == chip` holds exactly in every bundled manifest, and the
ring is complete in all five chips.

## What the chip did to the numbers (median per site, new / committed)

| site | `A` | `D` | `F` | `S` | `gamma2` | states |
| --- | --- | --- | --- | --- | --- | --- |
| CTS | 197,5 / 29 | 0,0429 / 0,0094 | 58 / 14 | 19,1 / 20,2 | 0,0083 / 0,0511 | pre 136, post 42 |
| KDLO | 75 / 84,5 | 0,0561 / 0,0306 | 10,5 / 49,5 | 15,7 / 60,7 | 0,0380 / 0,0095 | pre 16, single 1, rebuild 43 |
| Carolabrücke | 205 / 10 | 0,0243 / 0,0078 | 65,5 / 3 | 16,4 / 14,3 | 0,0029 / 0,2086 | pre 184, post 252 |
| Haltern | 376 / 3 | 0,0475 / 0,5000 | 189 / 2 | 19,0 / 0,95 | 0,0038 / 0,6006 | pre 96, post 11 |
| Ponte Morandi | 26 / 11 | 0,0243 / 0,0071 | 14 / 4 | 6,68 / 8,56 | 0,0204 / 0,0794 | pre 198, post 9 |

Two effects are entangled and must not be read as one: the chip is much larger
(more cells pass the threshold, so `A` and `F` rise), and the reference moved out
of the measurement (so `gamma2` drops by one to two orders of magnitude). Neither
number is "better" than the committed one — they are different measurements, and
that is why every mask column is written twice.

## The refusal finding: an ASC track can be too short for the chip

Haltern refuses 13 of 120 dates and Morandi 25 of 232 — **all of them
ASCENDING**, all with the same reason:

```
ValueError: the window ends at row 15137, past the raster's 15120 rows
```

The anchor sits 17 rows before the end of the burst raster and the chip needs 17
rows of half-height, so the window cannot be cut from that burst. The tools now
try the **neighbouring burst of the same acquisition** (`candidates_by_date`, a
date is not a burst), but for these dates the catalogue returns exactly one
candidate — there is no neighbour to fall back to. The refusals are recorded per
date with orbit and reason in the manifest.

This is a consequence of the larger cell, not a defect: the committed 80 x 80
(Haltern) window fits at that edge, a 35-line support with a height shift does
not. Haltern therefore keeps 107 of its 120 dates.

## Ponte Morandi: the deck comes from the project's own record

The OSM query at the site returns the **replacement** bridge
(`(Ponte san Giorgio) A10 ...`, way 732436706, 611 m), not the pre-collapse
viaduct the record measured, so OSM was not used. `morandi_edge_online.json`
carries the geometry itself:

| quantity | value | provenance |
| --- | --- | --- |
| deck centre | 44.425100, 8.889550 | midpoint of the two geodetic centres A and B |
| deck length | 244.1 m | `\|A->B\|` |
| deck width | 18.0 m | `thresholds.deck_width_m` |
| bearing | 316.85 deg | `atan2(A->B)` |

The rectangle is built in UTM from those four numbers and nothing else. The
reconstructed A and B land on the recorded centres to 1e-6 deg, and the deck's
axis in the radar frame comes out at **152.8 deg** against the project's own
`axis_px_bearing_deg` of **156.6 deg** — a 3.8 deg difference between a segment
through the two centres and a per-chip ridge fit, which is the kind of agreement
that makes the rectangle trustworthy rather than merely plausible.

The same edge effect as Haltern hits Morandi's ASC track: 25 of 232
dates are refused because the 51-line chip runs past the raster. DESC is
complete; the refusals are listed per date in the manifest.

## Espoo: the one site whose asset geometry needs an image

Espoo is the event-free control and the only site with **no mask layer** (its `A`
is the record's mask size, copied and proved), so it cannot be migrated to a
support chip at all. What is missing there is `base_m`: OSM carries the mast as a
point, the registry default for a lattice tower is 12 m, and the documented height
is 36 m.

A Google Maps Platform key works: `streetview/metadata` answers `OK` at the site
(pano `wk6r5QGt3DnY8yBZhL_KDw`, captured 2026-07, 49.2 m from the repository
coordinate). The reviewer opened the frame at pano `eyKXqe8bXRKVcBnJdsPogA`
(2026-07, 64.5 m away, yaw 97.46 deg) and reports the mast **on the right**, with
a **plausible height** against the documented 36 m.

`data/espoo/asset_evidence.json` records exactly that and nothing more. The
asset point is **OSM node 2233589804** (the reviewer's decision): it carries the
height the record documents (36 m) and the lattice/communication tags, and it is
the coordinate the site's SAR reference already uses. A point reported 50.7 m north
is recorded as *not adopted* (no corroboration — Overpass was unavailable) and a
fallback 3 869.8 m away as *rejected* (that is not the site the record measured).

* `mast_height_m = 36.0` — source: the committed reference, **confirmed
  qualitatively** by the reviewer (a 36 m mast subtends 29.2 deg at 64.5 m, so it
  is unmistakable).
* `base_width_m = null` — nobody read a width. At 64.5 m a 12 m base subtends
  10.5 deg and a 1.5 m car 1.3 deg, so a reading is possible; it was not made. The
  registry default therefore stands, **unconfirmed**.
* `position_check.status = explained_as_another_structure` — at yaw 97.46 deg a
  mast at the repository coordinate sits **26.4 deg left** of centre, not right.
  The repository point equals **OSM node 2233589804** to 0.0 m, and that node
  carries `height=36`, `man_made=mast`, `tower:construction=lattice`,
  `tower:type=communication`. So the OSM mast is the left one; what sits on the
  right is another antenna on the same telecom site. The asset is the OSM node,
  and nothing is anchored to a structure identified by eye in a panorama.

The imagery itself stays outside the repository (`/tmp/espoo_streetview/`): the
Maps Platform terms forbid storing it, and the model that runs these tools cannot
see images, so only a human reading can be recorded — with its provenance and
without a pretence of measurement.

## Reproduce one site

```sh
cd /home/projects/erlking_mvp
.venv/bin/python3 tools/ocv_migrate.py --site <site> --probe 3 --env-file <env>
.venv/bin/python3 tools/ocv_migrate.py --site <site> --cut   --env-file <env>
.venv/bin/python3 tools/ocv_table.py    --site <site>
```

`--geometry` runs the support arithmetic offline on one annotation XML, and
`tests/test_ocv_migrate.py` derives every expected number by hand.

## What this does not claim

- The two tables are **two measurements**, not a better and a worse instrument.
- A taller asset costs range: KDLO's 519.7 m mast is a 158.85-bin (≈ 527 m) range
  profile, and its chip is 475 bins wide because the support is.
- The refusals are geometry, not quality: a date that cannot hold the ring is not
  a date with "no signal".
- Nothing here is comparable with the committed per-state numbers until the state
  labels and the counts are read side by side — that is the SDS report's job.
