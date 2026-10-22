# Espoo mast (event-free control) — OCV trajectories across the event-free record (ASC R160 / DESC R51)

Trajectories of the six OCV dimensions across the event-free record (ASC R160 / DESC R51), one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), ASC / DESC across the columns. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

## States

**One state: event-free.** The record labels no damage state, so there is no pre/post contrast: the bands are runs of the single label and the medians are the level of the whole series.

## Rows and refusals

| channel | event-free | runs |
| --- | --- | --- |
| ASC | 31 | 1 |
| DESC | 44 | 1 |

Table `data/espoo/support/espoo_ocv_channels.csv` (75 overpasses, sha256 `82a8cc04208e`), cut manifest `data/espoo/support/espoo_rects_support.txt.manifest.json`: chip [35, 129], band [5, 33], support 165 cells, ring 3976 cells, guard 30.0 m, ring 200.0 m, **1 of 76 dates refused**.

## Per-state medians

**ASC**

| dimension | event-free |
| --- | --- |
| `gamma2` | 0.01639 (n=31) |
| `P` | 0.234 (n=31) |
| `D` | 0.021 (n=31) |
| `F` | 19 (n=31) |
| `A` | 76 (n=31) |
| `S` | 25.72 (n=31) |

**DESC**

| dimension | event-free |
| --- | --- |
| `gamma2` | 0.0007045 (n=44) |
| `P` | 0.1472 (n=44) |
| `D` | 0.1168 (n=44) |
| `F` | 133 (n=44) |
| `A` | 527.5 (n=44) |
| `S` | 14.39 (n=44) |

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
