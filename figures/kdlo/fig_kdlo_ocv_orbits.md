# KDLO media tower — OCV trajectories before and after the rebuild

Trajectories of the six OCV dimensions before and after the rebuild, one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), ASC / DESC across the columns. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

## States

Two states, as the record defines them: **pre-collapse** and **after rebuild** (table labels `pre` / `rebuild`).

Excluded: `single` (1 row(s)) — a state that carries a single overpass has no contrast and is not drawn as one.

## Rows and refusals

| channel | pre-collapse | after rebuild | runs |
| --- | --- | --- | --- |
| ASC | 0 | 0 | 0 |
| DESC | 16 | 43 | 2 |

Table `data/kdlo/support/kdlo_ocv_channels.csv` (59 overpasses, sha256 `9c6bb9197197`), cut manifest `data/kdlo/support/kdlo_rects_support.txt.manifest.json`: chip [35, 475], band [5, 379], support 1895 cells, ring 12280 cells, guard 30.0 m, ring 200.0 m, **0 of 60 dates refused**.

## Per-state medians

**ASC** — no overpasses of this channel in the table, so the panel is empty on purpose.

**DESC**

| dimension | pre-collapse | after rebuild |
| --- | --- | --- |
| `gamma2` | 0.0431 (n=16) | 0.03748 (n=43) |
| `P` | 0.5753 (n=16) | 0.4895 (n=43) |
| `D` | 0.07017 (n=16) | 0.05185 (n=43) |
| `F` | 10 (n=16) | 11 (n=43) |
| `A` | 54 (n=16) | 79 (n=43) |
| `S` | 22.93 (n=16) | 7.163 (n=43) |

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
