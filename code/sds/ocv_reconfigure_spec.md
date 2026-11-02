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

## Amendment 2 — 2026-10-22: the blocks are behavior and geometry

*The first amendment separated the two assets. This one splits the response the
other way — by what the dimension *is*, not by how it is observed.*

**The two blocks, and why.**

| block | dimensions | meaning |
|---|---|---|
| `behavior` | `γ²`, `P` | **cross-acquisition echo behavior** — `γ²` is the consistency of the echo across acquisitions, `P` its persistence/continuity across acquisitions. `P` belongs here even though its computation depends on valid echo detections: it describes how the echo behaves *between* acquisitions, not how the mask looks. |
| `geometry` | `A`, `D`, `F`, `S` | **spatial echo-mask geometry** — `A` extent (how much is observable), `D` occupancy, `F` fragmentation (how spatially organised the mask is), `S` the spatial shift. |

`D` was previously paired with `A` and `F` in a "mask" block that omitted `S`;
`S` sat with `γ²` and `P` in an "other" block. That cut did not follow the physics:
`S` is a spatial property, and `A`–`D`–`F`–`S` are one family. The literal
`[A, D, F]` sensitivity block is therefore dropped — `D` is no longer an extra vote
inside a two-dimension block, it is one of four geometry dimensions.

**Consequences, fixed before the run.**

* the magnitudes are `R_behavior = rms(z_γ², z_P)` and
  `R_geometry = rms(z_A, z_D, z_F, z_S)`; the two are **not** comparable with each
  other (two dimensions against four) and are never compared across assets;
* `R_geometry` keeps all four dimensions even though the figure's plane draws only
  `A` × `F`: `D` is kept out of the axes because it partly inherits `A`
  algebraically, and `S` is annotated because a shift describes position, not
  morphology;
* the Holm family stays **four tests per asset** (2 channels × 2 blocks), so the
  multiplicity does not change;
* the pre-registered predictions are re-stated per block: **P1** and **P2** are asked
  of the `geometry` block (the dose hypothesis is a mask-geometry hypothesis), and
  **P3** becomes "the `behavior` block need not follow the dose";
* the H2 comparison keeps one cosine per block and channel; the former
  "other block without γ²" variant is replaced by the two single-dimension `z`
  readings of the behavior block, because a one-dimension "vector" has no cosine.

The order of the blocks in the table, the figure and the report is
`behavior`, then `geometry`. The keys themselves keep the original spelling
(`behaviour`) — see amendment 7.



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

## Amendment 3 — 2026-10-22: a single asset with two states (Carola)

*The module was written for the LUMO/Bautzen pair. This amendment lets it run one
asset alone, and adds Carola to the registry.*

**Why Carola, and what changes.** Carola has exactly **two** states — `pre-collapse
(healthy)` and `post-collapse`, separated by the 2024-09-11 collapse — so:

* the ordered test collapses to its two-group case: a one-sided rank test between
  the two states, still through the same Jonckheere–Terpstra statistic, so `ΔR` is
  the post-minus-pre step;
* there is **no ordering claim**: with two states no ladder is tested, and the report
  says "the pre-to-post contrast" instead of "a dose ordering";
* **H2 is not computed** in a single-asset run (there is no partner), and the report
  says so instead of printing an empty table;
* the predictions P1–P3 belong to LUMO and are absent from a Carola run;
* the Holm family stays four tests (2 channels × 2 blocks), as in every other run.

**The table.** `data/carola/support/carola_ocv_channels.csv` — the support-sized
chip, **436 chips, `Girder A` only** (the committed table carries 1 717 chips over
Girders A/B/C, so nothing here speaks for B and C), 184 pre / 252 post, with
`gamma2_mask` as the coherence channel (log scale, like LUMO) and the `*_committed`
carry-over columns.

**What Carola cannot carry.** The two states are a **time split without a control
track**: Bautzen's bridge-minus-control difference removes what both tracks share,
Carola has nothing comparable, so its post signature is confounded with the
calendar. That is written into the report's closing section, together with the
girder limitation and the no-ordering statement.

**Outputs.** A single-asset run writes its own stem — `data/sds/ocv_carola.{json,md}`,
`data/sds/ocv_carola_table.csv`,
`figures/sds/fig_sds_ocv_reconfigure_carola.jpg|.md` — and never touches the pair's
files. `--sites` selects the assets; the default is still `lumo,bautzen`.


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

## Amendment 4 — 2026-10-26: the bridge pair (Carola, Morandi)

*Amendment 3 let the metric run on one two-state asset. This one runs it on two at
once, so the two collapses can be read side by side — the second pair, added
without a new hypothesis class.*

**Why a second pair.** Carolabrücke and Ponte Morandi are the two experiments in the
corpus that are (a) a bridge, (b) read across a collapse, and (c) available as a
support-sized chip produced by the same window rule. Comparing them asks the
question the dose pair asks — does the same metric, applied twice, answer the same
way? — in the one setting where the answer is not confounded with a dose ladder,
because neither asset has one.

**The registry entry.** `morandi` is registered exactly like `carola`: the same
eight dimensions (`A, D, F, S`, `P`, `gamma2_mask`), `log` for `A, D, F, P, γ²` and
`identity` for `S`, states `pre-collapse (healthy) → post-collapse`, reference
`pre-collapse (healthy)`, largest state `post-collapse`, the two-level ladder
`pre → post`, and `data/morandi/support/morandi_ocv_channels.csv` — **207 chips of
the committed 232**, one platform track, both orbit directions, zero drops for any
dimension. Nothing else in the registry changes; the dose pair's entries are
untouched.

**What H2 becomes.** The comparison keeps the cosine, its moving-block interval, the
sign agreement and the random-direction null, and loses what a two-state pair cannot
support:

* `largest` only — both assets have exactly one response state;
* no `middle` (no intermediate state exists) and no 2 × 2 matrix (only `post|post`
  carries a response); both omissions are declared in the JSON under `not_computed`
  rather than left as absent keys;
* the null stays coarse by construction: 2 dimensions give 2 random sign patterns,
  4 dimensions give 16 — so a 4-D `p` is a guard, never a significance claim. For
  this pair the one shared direction (DESC geometry, `A− D− F− S+`, cosine
  0.83, 4/4 signs) carries null `p = 0.083` and an interval that still covers zero.

**What the pair may claim.** That Carola resolves a large, coherent DESC
reconfiguration (both blocks, Holm < 0.001) and Morandi resolves nothing — its
post-collapse cell holds 3 ASC and 6 DESC chips, the ASC one at the bootstrap block,
its largest point estimate (DESC behavior, `R` = 0.98) sits at Holm 0.797 and every
other cell at 1.000 — and that exactly one of the four cells shares a direction. It
may **not** claim that the two bridges reorganised alike, and no magnitude crosses
between them — each `R` is on its own reference's scale.

**The thin-cell rule, now written down.** Amendment 3 quoted cells with `n = 3`
without marking them. Any cell whose `n` does not exceed the bootstrap block
(`BLOCK_BOOT = 3`) is now flagged everywhere it is quoted: `✱` beside its `n` in the
plate labels and a legend line in the suptitle, a `✱` in the report's state table,
and a `bootstrap_caveats` list in the JSON naming site, channel, block, level and
`n`. Such a cell is resampled as a single block, so its interval carries the
*reference's* uncertainty only; the interval is still reported, with that caveat.
The flag is computed per run, so the dose pair's reports carry an empty list.

Two thresholds coexist and are not the same thing: `SMALL_N = 6` dates from the first
pre-registration and writes `thin` in the table's `Support` column, while `✱` marks
the narrower case in which the resample block itself exceeds the cell
(`n ≤ BLOCK_BOOT = 3`) — the only case in which the interval is the reference's
alone. Bautzen's `DS2` cell (n = 5) is therefore `thin` in the table and unmarked on
the plate.

**One key added to the dose block.** `dose.applies` is `true` only when LUMO is in
the run. A bridge run still records the frequency ladder for provenance — it is the
same JSON schema — but flags it `applies: false`, prints `## There is no dose`
instead of the ladder, and puts the per-state `n` in its place.

**Outputs.** A registered pair keeps its own stem: `data/sds/ocv_carola_morandi.{json,md}`,
`data/sds/ocv_carola_morandi_table.csv`, and the plates
`figures/sds/fig_sds_ocv_reconfigure_carola_morandi.jpg` (2 × 4 juxtaposition) plus
the two single-asset plates. `--sites carola,morandi` selects it; the default is
still `lumo,bautzen`, and a single asset is still `--sites carola`.

**Self-test and verify.** Two checks are added to `--selftest` (17 total): a
registered pair is recognised in either order while a single asset is not, and a
cell at or below the bootstrap block is flagged while a larger cell is not.
`--verify` re-derives the expected caveat set from the units and fails if the
committed `bootstrap_caveats` names a different set of cells.

## Amendment 5 — 2026-10-26: the tower pair (Haltern, KDLO)

*Amendment 4 let two bridges be read across a collapse. This one runs the same metric
on two towers, and one of them was rebuilt — which is where a second pair stops being
a copy of the first.*

**Why this pair.** Haltern (an AV9 wind turbine that came down on 2021-09-29) and KDLO
(a media tower whose mast collapsed between 2022-12-11 and 2022-12-23 and was rebuilt)
are the two experiments in the corpus that are (a) a tower, (b) read across an event
that removed the structure, and (c) available as a support-sized chip produced by the
same window rule as the two bridges (`erlking_mvp/tools/ocv_table.py`, `30 %` of the
peak and `5 ×` the ring median, 200 m ring, 30 m guard). They differ from the bridges
in exactly the way that matters for the reading: Haltern's second state is a ruin and
KDLO's is a *replacement*.

**The registry entries.** Both are registered like `carola`: the same six dimensions
(`A, D, F, S`, `P`, `gamma2_mask`), `log` for `A, D, F, P, γ²` and `identity` for `S`,
zero drops for any dimension on either chip. `haltern`: states `pre-collapse (healthy)
→ post-collapse`, reference `pre-collapse (healthy)`, largest `post-collapse`, ladder
`pre → post`, `data/haltern/support/haltern_ocv_channels.csv` (**107 of the 120
committed overpasses**, both orbit directions: 44 ASC and 52 DESC pre-collapse, 5 ASC
and 6 DESC post). `kdlo`: states `pre → rebuild`, reference `pre`, largest `rebuild`,
ladder `pre → rebuild`, `data/kdlo/support/kdlo_ocv_channels.csv` (**60 of 60** — the
support is the whole record, one orbit direction: 16 pre and 43 rebuild, all DESC).

**One acquisition excluded, and counted.** KDLO's `single` state is the one acquisition
at the collapse itself (2022-12-23). It is removed by the registry before any statistic
and counted in the JSON audit (`audit.excluded = {"single": 1}`) and in the report, not
folded into `pre` or into `rebuild`: one date cannot carry a median of its own, and
registering it would invent an `n = 1` cell. The rule is general — an excluded state is
never a registered state, and a state that is neither stops the run.

**What changes in the code, and why it is data-driven.** KDLO has one orbit direction,
so a naive run would have emitted four empty ASC cells, a four-test Holm family with two
`None`s, and an `✱` on a cell with `n = 0`. Instead:

* **channels come from the record.** `channels_present(units)` reads the orbit directions
  the table actually carries; every cell, test, caveat, table row and panel iterates
  that tuple. A channel the record lacks is *declared* (`channels_absent`), carries no
  cell, contributes no test to the Holm family (`2 × channels`, so 2 for KDLO and 4 for
  Haltern) and is drawn as an empty panel that says so — never as a null;
* **H2 runs on the shared channels.** A cosine needs a vector on each side, so KDLO's
  run carries one channel, two cosine cells instead of four, and an `ASC` key under
  `not_computed`;
* **the state names come from the registry.** KDLO's level is `rebuild`, so the report
  says "pre-to-rebuild", the contrast is named from `LEVELS`, and `post` is never
  hardcoded;
* **the pair's prose comes from the registry.** `PAIR_NOTE` holds the headline, the
  noun, the unit, the window-rule sentence and the H2 rule per pair; the closing
  paragraph's lead clause is chosen from the data (both, one or neither asset resolves),
  and the shared-direction count is now a claim only when a cell has *both* a positive
  cosine and all signs equal.

**What the pair may claim.** Resolution or its absence, per asset and per block, in that
asset's own Holm family; the count of cells where the two directions agree; and, for
Haltern, the bound its interval licenses. KDLO's rebuilt state is a different structure,
so a resolved contrast there is a statement about the rebuild, not about damage; Haltern
holds 5 ASC and 6 DESC post-collapse overpasses, so a null there is a statement about
what that many acquisitions can resolve. Neither a shared direction nor a bound crosses
between the two assets, and no magnitude is compared.

**Outputs.** `data/sds/ocv_haltern_kdlo.{json,md}`, `data/sds/ocv_haltern_kdlo_table.csv`,
`figures/sds/fig_sds_ocv_reconfigure_haltern_kdlo.jpg` (2 × 4 juxtaposition, two empty
ASC panels in KDLO's row), plus the plates for `--sites haltern` and `--sites kdlo`.

**Self-test and verify.** `--selftest` grows from 17 to 27 checks: the registry covers
every asset in every table, the reference is the first level and the largest state the
last, an excluded state is never registered, the tower pair keeps its own stem, KDLO
loads 59 units with `{"single": 1}` excluded, Haltern loads 107 with both channels and
no drops, channels are read off the units, an absent channel contributes no caveat cell
and an empty caveat list no sentence, and a pair's cosine runs only on the channels both
records carry. `--verify` re-derives `channels_present`, the exclusion count, the level
sizes (which must sum to the units in those states) and the family size
(`2 × channels`) as well as every cell.



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

## Amendment 6 — 2026-10-30: the plate, authored at 170 mm

*Every amendment so far changed what is computed. This one changes only the sheet the
numbers are drawn on — and it does so because the pre-registered grid could not be
built at the repository's type floor.*

**What the pre-registration asked for.** A 2 × 4 plate: rows = experiment, columns =
channel × block (ASC-behavior, ASC-geometry, DESC-behavior, DESC-geometry), one
arrow per state from the reference into standardised space, and the bootstrap
interval of `R` as a thin bar.

**Why that grid cannot carry the type.** `figure_style.FINAL_WIDTH_MM` is 170 mm and
the ladder's floor is the annotation rung, 6.2 pt. Four columns inside 170 mm leave
`(170 − 14 − 3 − 3 × 11) / 4 = 30.0 mm` per square z-plane. At that rung the longest
state label measured on the real type is `DAM 6 (+2.4; n=17)` at 24.3 mm and
`rebuild (0.56; n=43)` at 24.8 mm; the median label is 21.6 mm. Add the marker and
its 0.8 mm clearance and a 30 mm plane has no side on which a label fits without
sitting on the data or leaving its slot. A plate that meets the width rule by
shrinking the type to 4 pt is not the same plate, so the width rule wins on width
and the grid gives way on grid.

**What replaces it.** Two columns (the block) and one row per experiment × channel:
ASC above DESC, four rows, eight panels. The reading order of the pre-registered
grid survives — experiment first, then channel, then block — and so does the
juxtaposition, because each experiment still carries both blocks in one band of two
rows. Geometry: plane 44 mm square, slot 71.0 mm, `x_stretch = slot / plane` so the
plane keeps its units and the added width is the margin the state labels sit in.
Canvas 170 × 237 mm for a pair, 170 × 125 mm for a single asset.

**Authored at 170 mm, never rescaled.** The canvas *is* the export box: `.jpg` at
630 dpi, `.pdf` (vector) and `.eps` (alpha blended, text Type 42) written from it,
`bbox_inches` never `tight`. A figure that is later rescaled by LaTeX has a type size
nobody chose; a figure that is cropped to its content has a width nobody measured.
`--verify` re-measures all three containers from their own headers and fails on a
0.5 mm deviation.

**The type ladder on the plate.** Header 9.0 pt, panel title 8.8, axis label 7.6,
tick label 6.8, annotation 6.2 — every one of them read from `style.FS`, so the plate
cannot drift from the rest of the corpus, and nothing on it is set below the
annotation rung.

**Why this plate keeps its headline.** The print plates of the corpus carry *keys,
not titles*: a bare tag per panel (`a1` … `f2`, `figure_style.tag_axes`), a key line
in the legend band and a `## Caption` in the report. This one does not, and the
exception is **named** instead of implied: it is a documentation plate of
`audit_figure_style.DOC_DIRS` (`concept/`, `sds/`), the audit's grammar clause skips
it by that name, and its row prints `grammar: documentation plate, header by name`
(source) and `doc plate (concept/sds)` (artefact). A page mock-up is read as a page,
not as the finished figure; the moment this pair is cut for print it moves under the
grammar — tags, a key line, a caption — rather than carrying the exception with it.

**The `R` strip.** One inset per panel repeats the block magnitude per level with its
moving-block interval, the quantity H1 tests. Its height follows the number of levels
instead of a fixed box: 28 mm wide, one row per level at 3.4 mm (a 6.8 pt label plus
its gap), the level names on the left (10 mm — `rebuild` is 8.6 mm at that rung), the
interval's ticks below (3.6 mm) and the `R` tag above (3.2 mm): 17.0 mm for three
levels, 13.6 mm for two. The strip is placed *before* the state labels, in whichever
of eight windows (four corners, four edge midpoints) is freest of markers, and if the
labels do not fit around the strip's footprint the next window is tried; the labels
themselves step in millimetres off their marker and are measured by the renderer, not
estimated from character counts. Nothing is drawn over anything, and if no window
works the build stops rather than covering a datum.

**What the guards hold.** The header and the ✱ line stay inside the canvas and 0.3 mm
from its edge; the header clears the first row of titles by 1.5 mm and the ✱ line by
0.4 mm; the drawn content stays 0.3 mm inside every edge. The ✱ line names its thin
cells when they fit on one 170 mm row and states the rule alone when they do not (the
guide keeps the list either way). Every plate is then re-read against the report:
`assert_figure_matches_report` compares each drawn `z`, `R`, `n` and interval with the
JSON, so a layout change cannot move a number.

**Self-test and verify.** `--selftest` (40 checks) gains the plate arithmetic: for
1, 2 and 4 row plates the slots, the gap and the gutters sum to exactly
`FINAL_WIDTH_MM`, `slot / x_stretch` is the plane, the height is the sum of its bands,
the grid stays inside the canvas, and the typography floor is the annotation rung.
`--verify` additionally re-measures the width of every container, re-compares the
audit block in each guide with the plate that was built, and re-derives every number
the guide quotes from `data/sds/`.


## Amendment 7 — 2026-10-30: American spelling, and a header that names the states

*Neither change touches a number. The first is the manuscript's variety of English;
the second removes a phrase that read like a claim the design refuses.*

**Every printed name is American English.** The manuscript is set in American
English, so the plate says `behavior` — in the panel titles, the guide, the tables
and the CSV column `R_behavior`. The **key** keeps the spelling this document
registered (`behaviour`): `BLOCKS`, the report's JSON, the H1/H2 cell ids and the
plate audit all key on it, and the other generators index `BLOCKS` by name — one of
them asserts the dict verbatim against its own copy. One table, `BLOCK_DISPLAY` in
`code/sds/ocv_reconfigure.py`, maps key to printed name, and nothing but a label goes
through it — no JSON key, no CSV field of the report itself, no audit field. A reader
therefore never sees two spellings of the same block in one document; a program never
sees two names for one key.

**The pair header.** It read *OCV reconfiguration against the state dose — LUMO and
Bautzen (no shared scale)*. "State dose" reads like one ordered loading axis that
both assets sit on, and then the parenthesis has to take the claim back. There is no
such axis: LUMO's ladder is a measured frequency-drop dose, Bautzen's is an ordinal
rank, and no panel of the plate shares a scale with any other — each arrow starts at
its own reference and is standardised by its own scale. The header now reads

> OCV reconfiguration across controlled LUMO and Bautzen states

and the parenthetical is gone, because the title no longer claims a shared scale for
the guide to deny: the standardisation is stated where it is made, in the panel
geometry and in the guide's first paragraph. The two registered pairs keep their
headers ("across a collapse", "across a collapse and a rebuild"): a collapse *is* a
common event, and its "(no shared scale)" names the one thing a reader might still
assume. The `--sites` help text and the default report's title line follow the same
wording.


## Amendment 8 — 2026-10-30: the headline comes off the plate, and the band goes to the plane

*Amendment 7 settled what the header said. This one removes it. Presentation only: no
number, label, sample size, interval or state order changes — `git diff data/sds/` is
empty for this amendment.*

**What was redundant.** The plate carried one 9.0 pt headline above the grid, *OCV
reconfiguration across controlled LUMO and Bautzen states* (amendment 7). The same
sentence is the report's first line and the caption's first line, and on the plate it
sat 2 mm above a band that already names every panel. A reader who has reached the
figure has read the section; a plate that repeats its own caption spends a band it
needs for the panels on words the page is not missing.

**What the band is now: the sum of what is drawn, not a constant.**
`plate_layout(rows, thin=...)` composes it — 2.0 mm of air, the block identifier's own
line (an 8.8 pt title plus the 2 pt pad it sits above its axes with, 4.7 mm), and only
on a plate that marks a thin cell the ✱ line (3.1 mm) with its 1.5 mm of clearance:
**6.7 mm** on the shipped pair, 11.3 mm on a plate that draws a ✱ line, against the
15.0 mm that was reserved for all of them before. The freed millimetres went into the
plane: **46 mm square on every row instead of 44** (`x_stretch` 1.5435 in the same
71.0 mm slot, 12.5 mm of lateral margin per side where the state labels sit). Canvas
170.0 × 236.7 mm for the pair — 0.3 mm shorter than amendment 6, because the band gave
up 8.3 mm and the planes took 8.0.

**The identifiers, and no bracket.** Each panel is addressed by
`block_identifier(site, channel, block)`:

> LUMO / ASC / behavior · LUMO / DESC / behavior · Bautzen / DESC / geometry

Same 8.8 pt, same title rung, no heading and no pair in brackets. Amendment 6 put
`[gamma2, P]` and `[A, F]` after the block because nothing else on the plate listed
the plane's dimensions. The axes do: `z(gamma2)`, `z(P)`, `z(A)`, `z(F)` sit on the
two edges of every plane, so the bracket was the same two words twice on the same
46 mm square. What the identifier does not say is said by the drawing — and `D` and
`S` appear **only** in the block magnitude `R` (the inset) and in the `D… S…` note at
the geometry panel's top right, never as an axis: `D` partly inherits `A`
algebraically and `S` describes position, not morphology. The panel's plane is still
recorded as `dims` in the audit and in the guide's JSON, so a machine still knows
exactly which two dimensions a plane draws.

**A key's pad, not a title's.** The identifier sits `PANEL_TAG_PAD_PT` (2.0 pt) above
its axes, as a `tag_axes` key does, instead of the 3.0 pt a title took; the band's
4.7 mm is that pad plus one 8.8 pt line, and the guard measures the renderer instead
of trusting the arithmetic. The rung itself is unchanged — nothing on the plate is set
smaller than it was, and nothing is set larger: the freed space is *plane*, not type.

**The guards, rewritten for a band with no second row.** The first row of identifiers
stays `PLATE_MARGIN_MM` (2.0 mm) below the top edge of the canvas and the drawn
content 0.3 mm inside every edge (`get_tightbbox`, so a tight bbox cannot hide a
collision); a plate that draws a ✱ line keeps `PLATE_CLEAR_MM` (1.5 mm) between that
line and the identifiers, and the line still has to fit one 170 mm row and stop
0.3 mm short of the right edge. The export contract is untouched: one canvas, three
containers from it, `ps.fonttype = 42` → `/FontFile2` and `/FontType 42`, no
`bbox="tight"`, no clipping.

**The caption, and the audit clause.** The sentence is not deleted, it is moved: the
report and the guide carry it, and the guide's new `## Caption` is the one place the
plate's sentence is written. The plate's own record says so — the audit block now
holds `"caption"` and `"headline_on_plate": null`, adds `layout.band_mm`, and drops
`suptitle` from `final_font_sizes_pt`, because a plate that draws no headline
declares no rung for one. The audit clause narrows to match:
`audit_figure_style.grammar_verdict(grammar, doc=...)` fails a `suptitle` on **every**
plate, documentation plates included, and exempts `DOC_DIRS` only from the *key* half
— this plate's identifiers are admitted **by name**, and its row prints `grammar:
synthesis plate, block identifiers by name` (a mock-up that keys nothing prints
`grammar: documentation plate, no heading, no key`). The artefact row is unchanged:
`doc plate (concept/sds)`, `170.0 mm`, `h236.7`, `fonts 42`.

**Self-test and verify.** `--selftest` goes 40 → 44: the band is checked to be exactly
the items it draws, on both variants (with and without a ✱ line), to stay below the
8 mm a headline-free band must not exceed and below the 15.0 mm it replaced; the
identifier is checked to be `asset / channel / block` with no bracket; and the plane
is checked to have taken what the band gave up. The audit's `--selftest` goes 42 → 44
with the two halves of the narrowed clause — a documentation plate may name its
panels, and a `suptitle` still fails there. `--verify` follows: the guide's audit
block is compared on `caption`, `layout` and `final_font_sizes_pt` as well, so a stale
guide fails instead of ageing quietly.


## Multiplicity

Eight primary tests: 2 experiments × 2 channels × 2 blocks for H1 — four per asset since
amendment 2026-10-22 separates the two, and **two per asset** when a record carries a
single orbit direction (amendment 5: KDLO has no ASCENDING acquisition, so its family is
2, not 4, and no test is computed on a channel that does not exist). Holm within that
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
  (n = 4), Bautzen `DS2` (11 pairs), Haltern ASC post-collapse (n = 5).
* No cell, test, caveat or panel for a channel a record does not carry (amendment 5):
  KDLO has no ASCENDING acquisition, so its `ASC` cells are declared absent, not
  reported as nulls, and its Holm family has two members.
* No acquisition folded out of sight: a state the registry excludes (KDLO's `single`,
  2022-12-23) is removed before any statistic and *counted* in the audit and in the
  report.
* Bautzen's state is the calendar; nothing here claims a seasonal control for it
  beyond what the paired difference removes. The same holds for every two-state asset:
  its second state is confounded with the calendar, and for KDLO with the rebuild.

## Outputs

* `data/sds/ocv_reconfigure.json` — every number, the drop counts, the scale chain.
* `data/sds/ocv_reconfigure.md` — the table and its reading.
* `data/sds/ocv_reconfigure_table.csv` — `Experiment, State, Channel, n, R_behavior,
  R_geometry, Direction, Scale, Support`.
* `figures/sds/fig_sds_ocv_reconfigure.jpg` — 630 dpi, authored at the 170 mm print
  width, 236.7 mm tall for the pair: rows = experiment × channel (ASC above DESC),
  columns = behavior | geometry (amendment 6), one block identifier per panel
  (amendment 8), arrows from the reference state to each state in standardised
  space, the bootstrap interval of `R` in the panel's inset. Beside the jpg:
  `.pdf` (vector) and `.eps` (alpha blended), and a `.md` guide whose JSON blocks
  carry the audit of the plate that was written.

A registered pair writes its own stem beside those: `ocv_carola_morandi`
(amendment 4) and `ocv_haltern_kdlo` (amendment 5), each with `.json`, `.md`,
`_table.csv`, a `fig_sds_ocv_reconfigure_<pair>.jpg` juxtaposition and one plate per
asset. A single asset writes `ocv_<site>.*` and `fig_sds_ocv_reconfigure_<site>.jpg`.
