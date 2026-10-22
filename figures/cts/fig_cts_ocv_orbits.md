# Champlain Towers South (target footprint) — OCV trajectories before and after the collapse

Trajectories of the six OCV dimensions before and after the collapse, one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), `ASC` across the column. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

## States

Two states, as the record defines them: **pre-collapse** and **post-collapse** (table labels `pre-collapse (healthy)` / `post-collapse`).

## Rows and refusals

| channel | pre-collapse | post-collapse | runs |
| --- | --- | --- | --- |
| ASC | 136 | 42 | 2 |
| DESC | 0 | 0 | 0 |

Table `data/cts/support/cts_ocv_channels.csv` (178 overpasses, sha256 `3965a15b2053`), cut manifest `data/cts/support/cts_rects_support.txt.manifest.json`: chip [39, 145], band [9, 25], support 225 cells, ring 5010 cells, guard 30.0 m, ring 200.0 m, **0 of 178 dates refused**.

## Per-state medians

**Channel not drawn:** `DESC` carries no overpass at this site, so the plate has one column for `ASC`. The empty channel is kept in the tables below, where its n = 0 is the record.

**ASC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.009122 (n=136) | 0.007826 (n=42) |
| `P` | 0.2772 (n=136) | 0.2438 (n=42) |
| `D` | 0.04191 (n=136) | 0.04967 (n=42) |
| `F` | 56 (n=136) | 68.5 (n=42) |
| `A` | 193 (n=136) | 238 (n=42) |
| `S` | 19.95 (n=136) | 16.7 (n=42) |

**DESC** — no overpasses of this channel in the table, so it is not drawn as a panel.

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
