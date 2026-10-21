# Champlain Towers South (target footprint) — OCV trajectories before and after the collapse

Trajectories of the six OCV dimensions across the collapse, one panel per **dimension x observation channel**, drawn from the support-sized chip: `gamma2, P, D, F, A, S` down the rows (the order the LUMO and Bautzen figures use), ASC / DESC across the columns. `P, D, F` and `A` are logarithmic, `gamma2` and `S` linear.

The thin dark line in each panel is a **rolling median over 5 consecutive overpasses** (centred, at least 3 valid values). It summarises local evolution — a dip in it is a summary, not a state, and whether it is a drift or a level shift is a question this figure does not answer.

## States

Two states, as the record defines them: **pre-collapse** and **post-collapse** (table labels `pre-collapse (healthy)` / `post-collapse`).

## Rows and refusals

| channel | pre | post | runs |
| --- | --- | --- | --- |
| ASC | 136 | 42 | 2 |
| DESC | 0 | 0 | 0 |

Table `data/cts/support/cts_ocv_channels.csv` (178 overpasses, sha256 `3965a15b2053`), cut manifest `data/cts/support/cts_rects_support.txt.manifest.json`: chip [39, 145], band [9, 25], support 225 cells, ring 5010 cells, guard 30.0 m, ring 200.0 m, **0 of 178 dates refused**.

## Per-state medians

**ASC**

| dimension | pre-collapse | post-collapse |
| --- | --- | --- |
| `gamma2` | 0.009122 (n=136) | 0.007826 (n=42) |
| `P` | 0.2772 (n=136) | 0.2438 (n=42) |
| `D` | 0.04191 (n=136) | 0.04967 (n=42) |
| `F` | 56 (n=136) | 68.5 (n=42) |
| `A` | 193 (n=136) | 238 (n=42) |
| `S` | 19.95 (n=136) | 16.7 (n=42) |

**DESC** — no overpasses of this channel in the table, so the panel is empty on purpose.

## What this does not claim

- The panels show the **support-sized chip**, not the committed table: the chip is larger and the reference is a disjoint ring, so these values are not the committed ones.
- A band with few overpasses is a handful of points, not a distribution; the median carries its n.
- The rolling median is a summary over five overpasses, not a state and not a trend test; it is drawn so that the local evolution is visible beside the global pre/post difference.
- The state boundary is drawn between runs, not at an event date: the record labels the runs, and the band follows the labels.
