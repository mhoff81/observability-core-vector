# `code/cases/ocv_trajectory_spec.md` — the standard OCV trajectory plate

**Status:** implemented by `code/cases/ocv_trajectory.py`, drawn by
`figures/<site>/fig_<site>_ocv_trajectory.sh` (six sites) and
`figures/cases/fig_ocv_trajectory.sh` (all six in one run). `--selftest` asserts
every number below against the six support tables, so this document cannot drift
away from the code without a failing check.

## The question the plate answers

Every other figure of this repository answers a **contrast** question: does the
state differ, and by how much? The trajectory plate answers the **timeline**
question the six cases are actually read for:

1. *Where in the record does the vector move?*
2. *Does it move where the event is — and does it stay moved?*

and it answers it with **one design for all six sites**, because the comparison a
reader makes is between plates. A per-site design would turn that comparison into
a comparison of designs.

## R1 — one cut, the support-sized chip

The plate reads `data/<site>/support/<site>_ocv_channels.csv`, the cut the shared
orbit plate (`code/orbits/fig_ocv_orbits.py`), the SDS reconfiguration
(`code/sds/ocv_reconfigure.py`) and `code/cases/fig_cases_overview.py` read. The
committed chip is **never** mixed in: two chips are two measurements, and the
support cut is the one every site can produce.

## R2 — the case's own reference timeline, and the tail drawn

The reference is the case's own reference state, from its **first** to its **last**
chip; the tail is **every** chip of the post state. Nothing is trimmed to make the
windows look alike — the windows differ because the records differ. The
post-event tail is **drawn**, never cut off: a trajectory that stops at the event
cannot show whether the vector returned, stayed or kept moving.

| site | reference state | reference window | chips | ASC | DESC | event | tail | chips |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| carola | `pre-collapse (healthy)` | 2022-09-03 → 2024-09-09 | **184** | 124 | 60 | 2024-09-11 | `post-collapse` 2024-09-13 → 2026-08-23 | **252** |
| cts | `pre-collapse (healthy)` | 2015-09-21 → 2021-06-21 | **136** | 136 | 0 | 2021-06-24 | `post-collapse` 2021-07-03 → 2022-12-25 | **42** |
| morandi | `pre-collapse (healthy)` | 2016-08-24 → 2018-08-08 | **198** | 91 | 107 | 2018-08-14 | `post-collapse` 2018-08-14 → 2018-09-13 | **9** |
| kdlo | `pre` | 2022-06-02 → 2022-12-11 | **16** | 0 | 16 | 2022-12-11 | `rebuild` 2023-05-04 → 2024-09-19 | **43** |
| haltern | `pre-collapse (healthy)` | 2021-01-01 → 2021-09-17 | **96** | 44 | 52 | 2021-09-29 | `post-collapse` 2021-10-01 → 2021-10-17 | **11** |
| espoo | *all* (no state column) | 2024-09-05 → 2026-09-07 | **75** | 31 | 44 | none | none | 0 |

Two states are the rule, three occur once, one occurs once:

* **KDLO** carries `pre` / `single` / `rebuild`; `single` is **excluded before any
  statistic** — one acquisition at the collapse itself (2022-12-23) is not a state —
  and the exclusion is counted in the JSON (`cut.excluded`) and printed in the
  report. Its tail is the **rebuild**, not a collapse aftermath: the contrast there
  answers a different question, and the plate's own subtitle says so.
* **Espoo** has **no event and no state column**: the whole record is reference
  *and* trajectory, and there is no contrast row at all.

### Windows that were considered and rejected

| site | rejected alternative | why not |
| --- | --- | --- |
| carola | the last 12 months before the collapse (90 chips) | it would cut the record's own start, thin the baseline by half and remove no confound — the season is inside every window of this record either way |
| cts | the last 24 months (59 chips) | the six-year pre-trend **is** the site's reference; cutting it would hide the drift the CTS package measures, and the DiD layer is defined on the full pre period |
| morandi | a 12-month reference (94 chips) | 9 post chips against 94 or 198 leaves the same conclusion, but the longer window is the record's own and needs no justification |
| kdlo | folding `single` into `pre` | it is the acquisition of the collapse day itself; a state of one is not a state, and the orbit plate already drops it |
| haltern | starting at the first *echo-bearing* chip | every chip of this cut carries an echo (see R3), so the filter would change nothing while suggesting it had |
| espoo | splitting the record in halves and calling the first a reference | the record labels no event; a split would invent the very boundary the control exists to deny |

## R3 — per-block standardisation, each with its own reference set

The two response blocks of `ocv_reconfigure.BLOCKS` are used unchanged:

* **behaviour** `(gamma2, P)` — cross-acquisition echo behaviour;
* **geometry** `(A, D, F, S)` — spatial echo-mask geometry.

Each block is standardised on **its own** reference set:

* the behaviour block against the reference chips that **carry an echo**
  (`A >= 2`), because the two channels only exist where the mask rule holds;
* the geometry block against **every** reference chip.

The standardisation is **pooled over the orbits**, except at **Espoo**, whose only
contrast *is* the pass geometry: there the reference is per orbit (ASCENDING
n = 31, DESCENDING n = 44) and never pooled.

On the six support cuts the echo filter is a **no-op**: every chip of every one of
the six carries an echo (`A >= 2`), which `--selftest` asserts per site. The rule
is still applied — it is the rule the committed chips need, where 54 of Haltern's
120 chips carry none — and the report prints `n` per dimension, so a reader can see
the two sets coincide here.

## R4 — a robust z on the family's own scale

    t(x) = log(x) for gamma2, P, A, D, F      identity for S
    z    = (t(x) - median_ref) / (1.4826 * MAD_ref),  clipped at +-6

The transform is `ocv_reconfigure.TRANSFORM` (asserted by `--selftest` on the four
sites both registries know); the scale falls back to `IQR / 1.349`, then to the
sample standard deviation, then to `1.0` with the source named `degenerate`, and
the report prints which one was used. A block magnitude is the **RMS** of its
standardised dimensions — `ocv_reconfigure.magnitude` verbatim, `sqrt(mean(z^2))`
over the dimensions actually available.

A block magnitude is computed **per chip** and the plate and the report show the
**median of those chip values** — the typical joint displacement of a single
overpass. It is therefore not the RMS of the two panel medians beside it, and it
does not vanish when those medians are near zero: a month whose chips scatter in
opposite directions has a zero median per dimension and a large per-chip
magnitude. The contrast's `R_B` / `R_G` rows use the same per-chip values.

The co-variates of panel c1 are standardised the same way on the identity scale.
Where a column carries only a handful of distinct values it is a climatology, not a
per-acquisition measurement (Espoo's `temperature_c` carries 5 values on 31
ascending chips), and the report says so instead of letting the panel imply a
resolution it does not have.

## R5 — monthly bins, no interpolation, a seeded interval

* A month is a **bin**, not a date: its ordinate is the **median** of its chips'
  z, with its `n` printed in the report and in panel c2.
* A month with no chip is a **hole**: the line breaks, and the event month is
  drawn with its own marker and breaks it too. Nothing is ever interpolated over.
* The intervals in the plate are both built this way: the **contrast** — post minus
  reference on the z scale — and R6's **outside class**, from a **circular
  moving-block bootstrap**: `BOOT = 2000` draws, `BLOCK = 3` chips, `SEED = 7`,
  one resample per side, percentile interval at 2.5 / 97.5. The same construction
  the orbit plate uses, so the interval means the same thing in both.
* No wall-clock value enters any output: two runs are byte-identical.

## R6 — the outside class is tested, not only drawn

R2 declares a **third class**: healthy chips of the reference state that the chosen
reference window leaves out. The plate draws them — grey, hollow, in the reference's
own frame — and standardises nothing on them. On Carola that is **63** chips
(2024-01-01 → 2024-09-09, 9 months), which a reader would otherwise see in every
panel of the behaviour block and have no number for.

So the report tests them: the class's median against the reference's, per dimension
and per block magnitude, with **R5's own interval** — `BOOT` draws, `BLOCK` chips,
`SEED`, one resample per side — and **its own seeded stream**, so the contrast's
committed draws cannot move when the same construction is reused. The result is the
JSON's `outside` block and one table in the report, `n_ref` and `n_out` beside every
row.

Where a site has nothing outside — R2's window is the state's own span there, which
is five of the six — the block says `applies: false` with the reason, and the report
prints that sentence instead of the table. The section exists on every site, so the
six reports keep one structure.

Two consequences are worth stating, because the plate is easy to misread without
them:

* **A block magnitude carries no sign.** A class that pushes `z(γ²)` up while
  pulling `z(P)` down moves *inside* the reference frame and can leave `R_B` at the
  reference's own level. Carola's outside class does: `R_B` is **+0.01**
  [−0.15, +0.25] while `z(γ²)` is **+0.20** [−0.09, +0.60]. The magnitude is the
  quietest row of that table, and it is not evidence of stability.
* **No row separating is itself the finding.** Every row of Carola's table spans
  zero: read in the reference's own scales, the 63 excluded chips *are* the
  reference, and the January → September spread inside the class is what makes its
  intervals as wide as they are. A subset chosen **after** looking at the curve (the
  Jun–Aug summer chips, say) can come out tighter — but a post-hoc window is not a
  pin, so it stays out of the report, and the class R2 declared is what is printed.

## The twelve panels

| panel | block | what it draws |
| --- | --- | --- |
| a1 | behaviour | z(γ²) — the monthly median over the reference and the tail |
| a2 | behaviour | z(P) — the monthly median |
| a3 | behaviour | the behaviour plane (z(γ²), z(P)), month by month |
| a4 | behaviour | R_B — the behaviour block magnitude |
| b1 | geometry | z(A) — monthly median |
| b2 | geometry | z(F) — monthly median |
| b3 | geometry | the geometry plane (z(A), z(F)); D and S annotated, drawn in c4 |
| b4 | geometry | R_G — the geometry block magnitude |
| c1 | context | the record's own co-variates (weather), same z |
| c2 | context | the observation cadence per orbit, per month |
| c3 | context | echo coverage: solid = an echo, hollow = a mask without one |
| c4 | context | the trajectory matrix — six dimensions × every month |

Reference months are filled circles with a solid line, post months hollow squares
with a dashed line, a month that holds both a diamond, the reference window a grey
band, and the event a dashed vertical line. The legend carries only the entries a
given site can produce: an event-free site has no post month and no event line.

## Outputs, verification, determinism

```
figures/<site>/fig_<site>_ocv_trajectory.jpg    the 3 x 4 plate, 630 dpi
figures/<site>/fig_<site>_ocv_trajectory.md     the report (every number)
data/<site>/fig_<site>_ocv_trajectory.json      the plotted values + the pins
```

* `--selftest` — 201 checks, no file touched: the registry, the documented
  windows of R2 (counts **and** dates), the echo filter, log positivity, the four
  shared `TRANSFORM` entries, `BLOCKS`/`GEOMETRY_PLANE`, "the reference median
  standardises to zero", "no |z| exceeds the clip", "every monthly median is the
  median of its own bin", the RMS definition, the bootstrap's determinism, and
  R6's five rules: the outside block applies exactly where R2 chooses a window,
  it holds exactly the chips the window leaves out (per site, against the census
  above), its delta is the difference of the two medians, its bootstrap is seeded
  and reproducible, and its interval brackets its own delta.
* `--verify` — re-derives the JSON and the report from the table and fails on any
  disagreement; the table's sha256 is carried in both.
* Two consecutive runs are byte-identical (no wall clock, seeded bootstrap).

```bash
bash figures/cases/fig_ocv_trajectory.sh            # all six plates
bash figures/carola/fig_carola_ocv_trajectory.sh    # one site
python3 code/cases/ocv_trajectory.py --selftest     # the template's rules
bash figures/cases/fig_ocv_trajectory.sh --verify   # re-derive and compare
```

## What this plate does not claim

* **No damage verdict.** Every ordinate is a deviation from the case's *own*
  reference, in that reference's robust scales — two plates are not comparable in
  magnitude, so the comparison they support is of *shape* and of *where the
  movement sits*, never of effect size.
* **No significance on the curve.** Both intervals — the contrast's and the outside
  class's — are differences of medians, not bands on the curve; the monthly medians
  are drawn with their `n` and without one.
* **No cross-site pooling.** The plate is per site by construction, and the driver
  runs the six one after the other; nothing is averaged across sites.
* **No orbit contrast.** R3 pools the orbits (except at Espoo, where it
  standardises per orbit to make the two comparable). The orbit question belongs to
  `code/orbits/fig_ocv_orbits.py` and `code/haltern/fig_haltern_ocv_orbits.py`,
  which answer it with their own construction.
* **No test of the mask rule.** c3 documents coverage on this cut; because every
  support chip carries an echo, it documents rather than tests. The rule itself is
  studied offline in each site's `*_mask_sensitivity.py`.
* **No claim for KDLO's tail as an aftermath.** Its post state is the rebuild; the
  contrast is a tower before and after a reconstruction, which is a different
  question and is labelled as such.
* **No event line where there is no event**, and no contrast row on the control.


