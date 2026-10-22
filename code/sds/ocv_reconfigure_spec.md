# Pre-registration — OCV reconfiguration against the controlled state dose

*Written before the first run of `code/sds/ocv_reconfigure.py`. Every choice below is
fixed here, not in the code; if a choice changes, it changes in this file first and
the change is dated.*

## The question

Do the two controlled experiments increase the **magnitude** of OCV reconfiguration
as their independently defined state dose increases — and is the **direction** of
that reconfiguration the same in both?

* **H1 (within experiment)** — does the reconfiguration magnitude grow with the
  dose, inside each experiment, on each observation channel?
* **H2 (across experiments)** — are the standardised response *directions* similar,
  block by block, when each experiment is compared with its own reference?

The quantity is the **distance from each experiment's own reference state**. Raw
OCV values are never compared between the two assets.

## What the dose is

**LUMO** — measured, not asserted. `data/lumo/reference/lumo_state_frequencies.json`
carries the natural frequencies of the three damage states from the tower's own
controlled excitation campaigns (Welch PSD, mode consensus, n = 5 recordings per
state). The first-higher mode is the one that moves:

| state | first higher mode | drop vs. its own healthy | dose |
|---|---|---|---|
| healthy | 13.501 Hz | 0.00 % | 0 |
| DAM 3 | 12.501 Hz | −7.41 % | 7.41 |
| DAM 4 | 12.5009 Hz | −7.41 % | 7.41 |
| DAM 6 | 11.751 Hz | −12.96 % | 12.96 |

The fundamental (2.750 Hz) is unchanged in every state, and the second-higher mode
confirms the DAM 6 step (−9.63 % against +1.17 % / +0.62 %). **DAM 3 and DAM 4 are
one FFT bin apart (0.25 Hz resolution): they are one dose level, not two.** The
ordered test therefore has **three** levels, `L0 = healthy`, `L1 = {DAM 3, DAM 4}`,
`L2 = DAM 6`. A four-group Jonckheere–Terpstra test (healthy < DAM 3 < DAM 4 <
DAM 6) is **not** licensed by the physics and is reported only as a sensitivity
analysis, labelled as such.

**Bautzen** — ordinal only. `RS < DS1 < DS2` is the documented ladder of the
OpenLabs bridge experiment; there is no measured frequency per state, so **no
dose-regression test is run on Bautzen** and its groups stay ordinal.

## Reference, response, scale

One unit of observation:

* **LUMO** — one `record_class='coherence'` row (an overpass), both orbits as
  separate channels.
* **Bautzen** — one complete `(date, orbit)` **bridge − control** pair. The states
  are also the calendar, so the paired difference is the only contrast that can be
  read as structural.

For a dimension `j`, channel `c` and state `s`, on the transformed scale `g` below,
the **response** is the difference of medians and the **standardised response** is
that difference in units of the reference state's own variability:

    r[j,c,s] = median(g(X[j,c,s])) − median(g(X[j,c,ref]))
    z[j,c,s] = r[j,c,s] / scale[j,c,ref]

For `g = log` this is exactly the log-ratio of the medians; for `g = identity` it is
the difference of medians. Both forms of the plan are covered by one definition, and
`r` is reported in the transformed unit it was computed in.

`scale` is the **MAD of the reference state**, with a fixed fallback chain when it is
zero: MAD → IQR / 1.349 → standard deviation → the dimension is dropped for that cell
and the drop is reported. Which step of the chain was used is a field of the JSON,
never a footnote.

A second scale — the **pooled MAD over all states of the same experiment** — is
reported beside every `z` as a sensitivity analysis. The reference MAD is
conservative by construction (LUMO's healthy variability is wide), and the two scales
must never be mixed inside one comparison.

## Amendment — 2026-10-22: the two assets are reported separately (option C)

*Decided after the first run and before the second. Two of the four changes alter
the published numbers; both follow from a property of the record that the first run
exposed, not from the direction of its result.*

**1. Holm is computed per asset, not over the two experiments together.** The joint
eight-test family expressed a joint claim; from now on each asset carries its own
four-test family (2 channels × 2 blocks) and its own verdict. There is no combined
`reading` label any more.

**2. When the reference state has fewer than ten units, the pooled MAD is the
primary scale.** Bautzen's reference is 6 complete pairs (ASC) and 7 (DESC); a MAD
from six numbers is itself an estimate with a wide spread and came out two to four
times tighter than the pooled spread of the same dimension (`A`: 0.22 against 0.59
on the log scale), which inflates every `z` of that asset by construction. The rule,
fixed before the second run: `n_ref < 10` → `scale = pooled MAD`, `scale_used =
"pooled_mad"`, with the reference MAD kept beside it as a sensitivity. LUMO keeps
its reference MAD (n_ref = 40 / 81, `scale_used = "mad"`). The consequence is that
Bautzen's `R` values **shrink** and their intervals narrow; that is the intended
correction of a scale that was estimated from too few points, and it is recorded
here rather than applied silently.

**3. The reference level carries no interval.** `R(reference) = 0` by construction,
because the same median appears in numerator and denominator. The bootstrap of that
row (previously `[0, 1.02]` for LUMO and `[0, 135.5]` for Bautzen) is a draw of the
reference against a second draw of itself; the noise it shows belongs to `ΔR`, where
both levels share one draw and it cancels. The row is reported as `0` with
`by_construction: true` and no interval.

**4. H2 becomes a descriptive comparison.** The cosine stays — with its interval, its
random-direction null and the sign agreement — but it is no longer a hypothesis: it
moves behind the two per-asset sections and answers "is this the same metric applied
to two assets?", not "is the response transferable?". Each asset section ends with
its own verdict **and a bound**: what that asset's interval excludes.

The claim the pair can carry is therefore: *one metric, defined per asset and applied
twice; the mask geometry answers in opposite directions; neither asset's own interval
licenses a dose gradient.* Each half stands or falls on its own numbers.

---

## Transforms


| dimension | LUMO | Bautzen |
|---|---|---|
| `A` (masked cells) | log | log |
| `D` (= A / bbox_area) | log | log |
| `F` (connected components) | log | log |
| `P` (persistence) | log | log |
| `S` (centroid↔peak distance, px) | identity | identity |
| `gamma2` | log | **identity** (floor-corrected, can be negative) |

**Zeros.** A value `X ≤ 0` makes `log` undefined; such a unit is **dropped for that
dimension only**, and the count is reported per cell. `S` and Bautzen's `gamma2` are
on the identity scale, so their zeros are values, not gaps. Known at
pre-registration time: LUMO `gamma2` 1 zero of 178; Bautzen `A` 8 and `F` 8 of 142
(and one side of a pair at zero kills that pair for that dimension).

## Blocks

`A`, `D`, `F` are not three independent votes: `D = A / bbox_area` holds in
**178 of 178** LUMO rows and in every Bautzen row (`A == n_band` too). The primary
mask block is therefore

    mask   = [A, F]        (two dimensions)
    other  = [gamma2, P, S]

and the literal three-dimension block `[A, D, F]` is reported as a sensitivity
analysis. The magnitude of a block is the root mean square of its `z` values:

    R_block = sqrt( mean_j( z_j^2 ) )

`R` is a descriptive magnitude. It is **not** an evidence count: nothing here treats
the three components as independent tests.

**γ² is not the same quantity in the two experiments** (LUMO: whole-tower coherence
in (0, 1]; Bautzen: floor-corrected band coherence, 72 of 128 values ≤ 0). It stays
in the `other` block, but the cross-experiment cosine of that block is reported
**with and without γ²**, and the γ² row is the weakest link in H2.

## Tests

**H1 — ordered, per experiment, per channel, per block.** Jonckheere–Terpstra on the
**per-unit magnitudes** `R_i` (each overpass or pair gets its own `z_i` against the
reference median and scale; the reference's *median* magnitude is 0 by construction,
so a test on state medians would be degenerate). 20 000 permutations of the group
labels, seed 7, one-sided `p_increasing`. Groups: three dose levels for LUMO, three
ordinal states for Bautzen. ASC and DESC are never pooled.

Uncertainty: moving-block bootstrap, 2 000 draws, seed 7, blocks of consecutive units
within a state that never cross a state boundary; the reference median **and** scale
are re-estimated inside every draw, so the interval carries the reference's own
uncertainty. Reported: CI on `R` per state, and CI on `ΔR = R(L2) − R(L1)`.

**H2 — direction, not magnitude.** Cosine between the standardised response vectors
of the two experiments at their largest state (LUMO `DAM 6`, Bautzen `DS2`), per block
and per channel; the same at the intermediate state (LUMO `{DAM 3, DAM 4}` pooled at
one dose, Bautzen `DS1`); plus the full 2 × 3 matrix of state pairs as a descriptive
map, and the sign-agreement count (0–3) because three dimensions admit only a handful
of distinct directions. CI by the same paired bootstrap; null by 20 000 random unit
vectors in three dimensions.

**LUMO repeatability.** Each damage state is visited twice, separated by a return to
healthy, so the response of visit 1 is compared with visit 2 of the same state (sign
agreement and `ΔR` with its CI). This is not a test of H1; it is what makes LUMO's H1
interpretable, and it is the only place where the record itself can reproduce a state
effect.

## Multiplicity

Eight primary tests: 2 experiments × 2 channels × 2 blocks for H1. Holm within that
family, and the adjusted value travels with every quoted `p`. H2 is descriptive:
cosines are reported with their interval and their null percentile, never as
"significant".

## What would refute the claims

* **Common magnitude continuum** — `R` increases with dose in both experiments on the
  same channel, and `R(LUMO, DAM 6) > R(LUMO, L1) > 0`.
* **Common response profile** — the block cosines are high and the sign patterns agree
  after each experiment's own standardisation.
* **No continuum** — neither magnitude nor direction shows a progression.
* **Partial continuum** — only the mask block scales with dose, while `γ², P, S` stay
  experiment-specific.

Pre-registered predictions, written before the run:

* **P1** LUMO's mask magnitude rises from L0 to L1 to L2 (`R(L2) > R(L1) > 0`).
* **P2** `R(DAM 3) ≈ R(DAM 4)` — they are one dose level. A large difference between
  them is a season/visit artefact, because DAM 3 exists only in summer while DAM 4 is
  season-balanced, and it bounds that artefact for every other number.
* **P3** The `other` block need not follow the dose; if it does not, the honest reading
  is "the mask geometry carries the dose, the coherence dimensions do not".

## Refusals

* No pooling of the two experiments before standardisation, and no comparison of the
  *magnitudes* `R` between them — only directions.
* No severity derived from an OCV response, in either direction.
* No use of the failure cases in the construction of the metric; they are projected
  with frozen definitions afterwards, if at all.
* No "thin" cell is quoted without its `n`: LUMO ASC `DAM 4` (n = 5) and `DAM 6`
  (n = 4), Bautzen `DS2` (11 pairs).
* Bautzen's state is the calendar; nothing here claims a seasonal control for it
  beyond what the paired difference removes.

## Outputs

* `data/sds/ocv_reconfigure.json` — every number, the drop counts, the scale chain.
* `data/sds/ocv_reconfigure.md` — the table and its reading.
* `data/sds/ocv_reconfigure_table.csv` — `Experiment, State, Channel, n, R_mask,
  R_other, Direction, Support`.
* `figures/sds/fig_sds_ocv_reconfigure.jpg` — 630 dpi: rows = experiment, columns =
  ASC-mask, ASC-other, DESC-mask, DESC-other, arrows from the reference state to each
  state in standardised space, bootstrap intervals as thin bars.
