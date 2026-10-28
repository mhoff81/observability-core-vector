# Ponte Morandi, Genova — the prospective read, two readers over one frozen reference

**Case** Ponte Morandi, Genova · **event** 2018-08-14 · **frozen window** 2016-09 → 2017-08 (**109** chips, 12 months, selected by a rule written before the read) · **held out** 2017-09-02 → 2018-08-08 (**88** chips, 12 bins) · **control** Carolabrücke, Dresden, read after the case · **policy** erlking-eas 1.0.0, pinned (`data/cases/eas_policy_pin.json`, sha256 `e0557e90cbad29f4`)

```
window rule   the earliest run of twelve consecutive calendar months in which every month holds at least four chips
              reads the acquisition date of every chip; the number of chips in each calendar month
              never reads the event date; any measured quantity of any chip
envelope      monthly medians of the 12 frozen reference months, q05/q50/q95 (empirical)   decision quantile q95
persistence   trajectory rule 2 consecutive months   policy v1.0.0 2 consecutive months   lookback 4 months
policy floors block |z| > 2.5   dimension |z| > 3.5   alert_enabled False
readers       trajectory and live  -- co-primary, never collapsed
```

Two readers look at the same 12 held-out monthly states: the trajectory rule of this repository, whose boundary is the reference months' own q05/q95, and EAS policy v1.0.0, whose boundary is a fixed `|z|` floor. Both use the same persistence rule -- 2 consecutive months -- so the whole difference between their words is a distance between two boundaries, and that distance is measured in Q8.

## Q1 · What window was frozen, and by what rule?

Earliest full cycle, written before the read ran: **the earliest run of twelve consecutive calendar months in which every month holds at least four chips**. The rule reads the acquisition date of every chip and the number of chips in each calendar month and nothing else — in particular it never reads the event date or any measured quantity of any chip, so the size of the held-out spread cannot have steered it.

It was tried from the first occupied month and kept the first start that held:

| start tried | would end | months under 4 chips | verdict |
|---|---|---|---|
| 2016-08 | 2017-07 | 2016-08 (1) | skipped |
| 2016-09 | 2017-08 | none | kept |

So the reference is **2016-09 → 2017-08**, 109 chips over 12 calendar months:

| month | 2016-09 | 2016-10 | 2016-11 | 2016-12 | 2017-01 | 2017-02 | 2017-03 | 2017-04 | 2017-05 | 2017-06 | 2017-07 | 2017-08 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chips | 6 | 10 | 9 | 10 | 9 | 7 | 9 | 10 | 9 | 10 | 9 | 11 |

The guard that keeps this prospective rather than retrospective: the window ends **12 calendar months** and **351 days** before the event month's first day, against a 12-month guard — passed. A rule of 365 days would have refused this window by 14 days while the event month is 12 bins away, which is why the guard is counted in months.

The registry declares no window for this site — its reference state spans 198 chips over 25 months and stops at the event — so this window *narrows* what was published rather than contradicting it. Q11 checks the selection against the data it was not allowed to see.

## Q2 · What was held out, and where did every other chip go?

| bucket | what it is | chips | months | window |
|---|---|---|---|---|
| reference | frozen pre-event reference | 109 | 12 | 2016-09-01 → 2017-08-31 |
| outside | held-out pre-event | 88 | 12 | 2017-09-02 → 2018-08-08 |
| prefix | prefix | 1 | 1 | 2016-08-24 → 2016-08-24 |
| post | post-event | 9 | 2 | 2018-08-14 → 2018-09-13 |
| unassigned | unassigned | 0 | 0 | — |
| **total** | | **207** | | |

Disjoint and complete: one chip, one bucket: chip_class returns exactly one of reference / outside / prefix / post / unassigned. Every kept chip sits in exactly one bucket, chips after the event are kept but never read as evidence, and the one chip before the window is named rather than discarded silently.

## Q3 · Did either reader see anything before the event?

The trajectory rule first, because it is the reader with the tighter boundary. Each of the eight streams is standardised against the frozen reference months and compared with its own q05/q95 envelope, over the 12 held-out bins (2017-09 → 2018-08 (pre), 88 chips, every one of them before the event):

| stream | months outside | longest run | run span | median z | min → max | envelope q05 → q95 |
|---|---|---|---|---|---|---|
| `R_B` | 7 / 12 | 6 | 2017-11 → 2018-04 | 1.11 | 0.89 → 1.42 | 0.59 → 1.04 |
| `R_G` | 4 / 12 | 2 | 2018-01 → 2018-02 | 1.66 | 1.07 → 2.12 | 0.55 → 1.71 |
| `gamma2` | 11 / 12 | 10 | 2017-11 → 2018-08 (pre) | -1.34 | -2.01 → -0.29 | -0.47 → 0.38 |
| `P` | 1 / 12 | 1 | 2018-08 (pre) → 2018-08 (pre) | 0.18 | -0.60 → 0.42 | -0.36 → 0.45 |
| `A` | 10 / 12 | 10 | 2017-11 → 2018-08 (pre) | 1.45 | 0.33 → 3.07 | -0.53 → 0.80 |
| `D` | 0 / 12 | 0 | — | 1.31 | -0.03 → 2.06 | -0.22 → 2.68 |
| `F` | 9 / 12 | 9 | 2017-12 → 2018-08 (pre) | 1.16 | 0.36 → 2.09 | -0.51 → 0.80 |
| `S` | 1 / 12 | 1 | 2018-08 (pre) → 2018-08 (pre) | 0.06 | -0.32 → 0.85 | -0.65 → 0.58 |

Longest runs, in held-out months: 10 on `gamma2`, 10 on `A`, 9 on `F`, 6 on `R_B`. They sit inside one span rather than scattering across the year, and the signs disagree: `gamma2` falls below q05 from 2017-11 while `A` and `F` rise above q95 from 2017-11 and 2017-12 — one departure seen on both sides of the envelope. `D` never leaves its envelope at all, and the quiet stream is named here so that the departure is not read as a property of the whole asset.

## Q4 · Was it sustained, or a single month?

The rule requires **2 consecutive months** outside the envelope, in a block magnitude and in a dimension, before it will say CHANGE; one of the two alone is WATCH:

| stream | block or dimension | longest run | vs 2 | sustained | fired |
|---|---|---|---|---|---|
| `R_B` | block | 6 | ≥ | yes | **yes** |
| `R_G` | block | 2 | ≥ | yes | **yes** |
| `gamma2` | dimension | 10 | ≥ | yes | **yes** |
| `P` | dimension | 1 | < | no | no |
| `A` | dimension | 10 | ≥ | yes | **yes** |
| `D` | dimension | 0 | < | no | no |
| `F` | dimension | 9 | ≥ | yes | **yes** |
| `S` | dimension | 1 | < | no | no |

Every run is a run of *consecutive held-out calendar months*, not a count of flagged months: `R_G` leaves the envelope 4 times, but its longest run without a gap covers 2018-01 → 2018-02; `P` and `S` leave it once each (2018-08 (pre)). Read up to and including the last held-out month, the rule returns **CHANGE** — *sustained multidimensional observability departure*.

## Q5 · How much lead, and under which definition?

Four definitions, four dates. They are printed side by side because 'lead time' without a definition is a slogan:

| definition | first month | bin tick | days before the event |
|---|---|---|---|
| first_flagged | 2017-09 | 2017-09-15 | 333 |
| first_persistence | 2017-12 | 2017-12-15 | 242 |
| first_change | 2017-12 | 2017-12-15 | 242 |
| first_live_fire | — | — | did not happen |

The signed lead is **333 days** (first month outside the envelope) and the sustained lead is **242 days** (first month read as CHANGE); the pinned policy never fires, so it has no lead to report. Dates come from the bin's own tick: the 15th of a whole month; the event is 2018-08-14. The month is the anchor, not the chip — a monthly state is dated by its tick, and dating it by the last chip in the bin would buy or lose days that no measurement supports.

The definitions the rows above are computed from:

* `first_flagged` — the first held-out month in which any of the eight streams lies outside the frozen envelope.
* `first_persistence` — the first month at which some stream had already been outside for 2 consecutive months.
* `first_change` — the first month at which the trajectory rule, read on the months up to it, returns CHANGE.
* `first_live_fire` — the first month at which a pinned policy rule would have fired, both blocks and dimensions considered.

## Q6 · Did the pinned policy fire?

Policy `erlking-eas` v1.0.0, pinned by digest (`ef2805e1c262b980…`) and applied by value: every threshold below is read from `data/cases/eas_policy_pin.json`, not from the policy's own text, so a change to the text cannot move a word here without the digest check in Q9 catching it.

Its six eligibility gates, read on this cut:

| gate | value | threshold | read as | passed |
|---|---|---|---|---|
| min_months | 12 | 6 | held-out bins | yes |
| min_reference_months | 12 | 3 | reference months inside the frozen window | yes |
| min_eligible_share | 1.00 | 0.5 | chips kept over chips offered | yes |
| min_complete_months | 12 | 4 | held-out bins carrying all eight streams | yes |
| max_degenerate_dims | 0 | 2 | dimensions with a missing or zero reference scale | yes |
| max_clip_share | 0.00 | 0.25 | standardised values clipped at the plate's own clip | yes |

Then its rules. `observability_drop` is listed but not applied: the rule reads a level of observability -- the eligible share of a month and its echo share -- and this cut carries chips, not levels: it can neither fire nor clear.

| rule | quantity | boundary | max abs | months over | longest run | vs 2 | fired |
|---|---|---|---|---|---|---|---|
| block_z | `R_B` | abs(R_B) > 2.5 | 1.42 | 0 | 0 | < | no |
| block_z | `R_G` | abs(R_G) > 2.5 | 2.12 | 0 | 0 | < | no |
| dimension_z | `gamma2` | abs(gamma2) > 3.5 | 2.01 | 0 | 0 | < | no |
| dimension_z | `P` | abs(P) > 3.5 | 0.60 | 0 | 0 | < | no |
| dimension_z | `A` | abs(A) > 3.5 | 3.07 | 0 | 0 | < | no |
| dimension_z | `D` | abs(D) > 3.5 | 2.06 | 0 | 0 | < | no |
| dimension_z | `F` | abs(F) > 3.5 | 2.09 | 0 | 0 | < | no |
| dimension_z | `S` | abs(S) > 3.5 | 0.85 | 0 | 0 | < | no |
| observability_drop | — | not applied | — | — | — | — | — |

**NOMINAL** — no rule this cut can read fires; one of the policy's three firing rules is not carried by an observability cut.

Read over every firing rule it names instead, the same policy returns **INSUFFICIENT**: observability_drop cannot be read from this cut. Both words are printed because only one of them is a clearance and neither is a verdict on the asset.

Scope, stated rather than assumed:

* the policy's `rules.ocv` carries five entries: three firing rules (`block_z`, `dimension_z` and `observability_drop`) and two parameters (`persist_months` and `lookback_months`). This cut applies two of the firing rules by value and cannot read `observability_drop`, which wants a per-month level and finds chips, so that rule is carried as missing rather than as clear.
* the lookback: the whole held-out span is the more permissive of the two readings, so a null across it implies a null under the engine's four-month lookback.

## Q7 · What does each reader say?

| reader | word | what the word is |
|---|---|---|
| trajectory rule (this repository) | **CHANGE** | CHANGE needs a sustained block magnitude and a sustained dimension; WATCH either one alone; STABLE neither; NOT_ASSESSABLE a held-out bin without all eight streams |
| erlking-eas v1.0.0, rules this cut carries | **NOMINAL** | no rule this cut can read fires; one of the policy's three firing rules is not carried by an observability cut |
| erlking-eas v1.0.0, every firing rule | **INSUFFICIENT** | observability_drop cannot be read from this cut |

Two vocabularies over one set of numbers, and the record does not collapse them: `concur` is **False** and `collapsed` is **False**. The words are not translations of each other — the trajectory rule prints `STABLE`, `WATCH`, `CHANGE` and `NOT_ASSESSABLE`, each word a reading of the reference months' own envelope, while the policy chooses among its seven words `INSUFFICIENT`, `NOMINAL`, `DEGRADED_OBSERVABILITY`, `WATCH`, `UNRESOLVED`, `CORROBORATED_WATCH` and `ALERT` over a fixed floor. A sentence that said 'the policy agrees' or 'the policy disagrees' would be inventing a mapping neither reader defines.

## Q8 · Why do they differ, and by how much?

Persistence is **not** the difference: the trajectory rule wants 2 consecutive months, the policy wants 2. The difference is entirely the boundary. For each rule, the largest floor at which that rule would still fire over 2 consecutive held-out months is the smallest `|x|` inside its best run — so the margin to the pinned floor is the distance between silence and a word:

| quantity | rule | pinned floor | largest firing floor | run | margin | fires |
|---|---|---|---|---|---|---|
| `R_B` | block_z | 2.5 | 1.325 | 2018-02 → 2018-03 | 1.175 | no |
| `R_G` | block_z | 2.5 | 1.897 | 2018-01 → 2018-02 | 0.603 | no |
| `gamma2` | dimension_z | 3.5 | 1.849 | 2018-02 → 2018-03 | 1.651 | no |
| `P` | dimension_z | 3.5 | 0.248 | 2017-09 → 2017-10 | 3.252 | no |
| `A` | dimension_z | 3.5 | 2.948 | 2018-01 → 2018-02 | 0.552 | no |
| `D` | dimension_z | 3.5 | 1.615 | 2017-12 → 2018-01 | 1.885 | no |
| `F` | dimension_z | 3.5 | 1.994 | 2018-01 → 2018-02 | 1.506 | no |
| `S` | dimension_z | 3.5 | 0.513 | 2018-07 → 2018-08 (pre) | 2.987 | no |

The closest of them is `A`, whose run 2018-01 → 2018-02 would fire at a floor of **2.948** against the pinned **3.5** — a margin of **0.552**. Nothing in the held-out span comes near the policy's floor under either of its block or dimension rules, while the envelope it is compared with is far tighter: `gamma2`'s q95 is 0.38 and `A`'s is 0.80, both against a dimension floor of 3.5. The trajectory rule therefore calls runs of 10 months on `gamma2` and 10 months on `A` a CHANGE while the policy, correctly applying its own floors, has nothing to report.

Read the other way: the two readers are not in conflict, they are answering at two distances. `gamma2`'s sustained 10-month departure is *large* against the reference months' own spread and *small* against a fixed absolute-scale floor — |z| reaches only 2.01 there, and the envelope's own q95 sits at 0.38. Naming that gap is the whole finding; a single word would have hidden it.

## Q9 · Does the finding survive the guards?

Every guard below is computed from the same artefact the answers above are read from, and every one of them passes — `guards.passed` is **True**.

| guard | passed | evidence |
|---|---|---|
| window_event_free | yes | 12 calendar months and 351 days before the event month |
| buckets_disjoint_and_reconciled | yes | 207 chips in, 207 reconciled across 5 buckets |
| no_pooling | yes | one build per site; the control never enters the case's numbers, which is checked by rebuilding the case after the control has been read |
| envelope_from_reference_months_only | yes | 12 reference months, monthly median of each calendar-month bin |
| every_held_out_bin_fielded | yes | 12 of 12 bins carry all eight streams |
| state_recomputes_from_stored_numbers | yes | trajectory CHANGE → CHANGE, live NOMINAL → NOMINAL |
| record_is_max_over_streams | yes | longest run 10, on `gamma2` |
| control_after_case_no_influence | yes | morandi → carola → morandi, identical: True |
| runs_reported_separately_from_magnitude | yes | every stream carries its longest run beside its own median and extremes; no sentence substitutes one for the other |
| magnitudes_at_1e_12 | yes | R_B = rms(gamma2, P); R_G = rms(A, D, F, S) over 414 chips, max difference 0 |
| vocabulary_scan | yes | both readers, case and control, each word in its own reader's list: emitted CHANGE/WATCH (trajectory) and INSUFFICIENT/NOMINAL (live); no stray word |
| no_failure_language | yes | scanned the pair, both states, the leads and the partition |
| hindsight | yes | three perturbations, all identical |
| thresholds_pinned | yes | 3 digests verified against data/cases/eas_policy_pin.json |
| event_is_the_registrys | yes | 2018-08-14 |
| policy_sha_matches_pin | yes | 3 of 3 digests match |

The hindsight guard deletes evidence and checks that nothing moves. It holds out the reference numbers and the held-out numbers separately, because a tail beyond the last held-out bin can only move one of them:

| perturbation | chips dropped | reference moved | held-out moved | bins vanished |
|---|---|---|---|---|
| tail_deleted | 9 | no | no | none |
| tail_absurd | 0 | no | no | none |
| reference_from_window_alone | 109 | no | no | none |

One label moves when the post-event chips are deleted — 2018-08: `2018-08 (pre)` → `2018-08`. That is the split's `(pre)` marker, owned by the split rather than by a number, so it is reported and excluded from the comparison instead of being allowed to pass as a difference.

And the guard that matters most for this file being *prospective*: 88 held-out chips, the latest on **2018-08-08**; the event is **2018-08-14**. Nothing after it is read as evidence.

## Q10 · Does the control behave the same way?

The control is **Carolabrücke, Dresden**, event **2024-09-11**, read *after* the case and never as an input to it. The two windows arrive from different places:

* **control** — the registry declares it, and its descriptive plate is published against it: the registry's window for this asset.
* **case** — the stated rule selected it: the earliest run of twelve consecutive calendar months in which every month holds at least four chips.

That asymmetry is the point: a rule-selected window and a registry-declared one are both legitimate, and the file says which is which instead of treating them as interchangeable.

Two scopes are named in this table, so both are labelled. The reference *state* is the plate's: every chip the table carries for it. The *frozen reference* is the slice this case builds its envelope from — the months inside each asset's window — and it is the smaller of the two.

| | Ponte Morandi, Genova (case) | Carolabrücke, Dresden (control) |
|---|---|---|
| reference-state months | 25 | 25 |
| reference-state chips | 198 | 184 |
| frozen reference months | 12 | 16 |
| frozen reference chips | 109 | 121 |
| held-out months | 12 | 9 |
| held-out chips | 88 | 63 |
| trajectory rule | **CHANGE** — sustained multidimensional observability departure | **WATCH** — isolated unusual departure |
| pinned policy | **NOMINAL** | **NOMINAL** |
| first month outside | 2017-09 (2017-09-15, 333 days) | 2024-01 (2024-01-15, 240 days) |
| first CHANGE | 2017-12 (2017-12-15, 242 days) | did not happen |

Per stream, longest run — the comparison is of *runs*, not of words:

| stream | case run | case months outside | control run | control months outside | differs |
|---|---|---|---|---|---|
| `R_B` | 6 | 7 / 12 | 1 | 2 / 9 | yes |
| `R_G` | 2 | 4 / 12 | 1 | 2 / 9 | yes |
| `gamma2` | 10 | 11 / 12 | 2 | 4 / 9 | yes |
| `P` | 1 | 1 / 12 | 3 | 5 / 9 | yes |
| `A` | 10 | 10 / 12 | 0 | 0 / 9 | yes |
| `D` | 0 | 0 / 12 | 0 | 0 / 9 | no |
| `F` | 9 | 9 / 12 | 1 | 1 / 9 | yes |
| `S` | 1 | 1 / 12 | 2 | 3 / 9 | yes |

The control does not reproduce the case: its longest run in any stream is 3 against the case's 10, its own word is **WATCH** against **CHANGE**, and its first month outside the envelope is 2024-01, which is 240 days before its event — a signed lead with no sustained departure behind it. That is what a control is for: the same machinery, run on a second asset, does not print the same word, so the case's **CHANGE** is a reading of those numbers rather than something the code prints everywhere.

Read in this order: the case, then the control, then the case again. The third reading is compared with the first and must be identical; the control is never an input to the case. The two assets are not pooled and no number here is an average of them.

## Q11 · Was the window selection genuinely event-free?

Three separate claims, each checkable without leaving this file:

1. **The rule reads no measurement.** It reads the acquisition date of every chip and the number of chips in each calendar month. It never reads the event date or any measured quantity of any chip — the resolver's argument list contains no chip measurement and no event date, so there is no path by which a held-out value could have moved the window.
2. **The rule was applied from the first month, not searched.** It starts at the first occupied month (2016-08) and keeps the first start that satisfies it. Only 2 starts were tried, and the first was rejected for a reason that has nothing to do with any measurement: 2016-08 held 1 chip.
3. **The guard is counted in months, not days.** The event month is 2018-08; the window's last month is 2017-08, 12 bins earlier and 351 days before the event date. With a 12-month guard the check passes; a 365-day guard would have refused this window by 14 days even though 12 whole bins separate it from the event.

What this does *not* claim: that the earliest satisfying window is the only defensible one. A later window would satisfy the same rule and would not be selected by it — the rule fixes *one* window in advance, and the record's job is to be honest that it is one, not to pretend it is the only one. The held-out span (88 chips, 12 bins) is not thinner than the reference (109 chips, 12 bins), so the test set is not a convenience sample from a sparse tail.

## Q12 · What would change the answer?

Numbers, not adjectives. Each line names the measurement and the value it would take to move a word:

| to move | the pinned value | it would take | which is |
|---|---|---|---|
| `R_B` fires under block_z | 2.5 | a floor of ≤ 1.325 | 1.175 below today's |
| `R_G` fires under block_z | 2.5 | a floor of ≤ 1.897 | 0.603 below today's |
| `gamma2` fires under dimension_z | 3.5 | a floor of ≤ 1.849 | 1.651 below today's |
| `P` fires under dimension_z | 3.5 | a floor of ≤ 0.248 | 3.252 below today's |
| `A` fires under dimension_z | 3.5 | a floor of ≤ 2.948 | 0.552 below today's |
| `D` fires under dimension_z | 3.5 | a floor of ≤ 1.615 | 1.885 below today's |
| `F` fires under dimension_z | 3.5 | a floor of ≤ 1.994 | 1.506 below today's |
| `S` fires under dimension_z | 3.5 | a floor of ≤ 0.513 | 2.987 below today's |

The tightest of these is `A` at **2.948**: the policy's dimension floor is **3.5**, the strongest single held-out value reaches `|z|` = 3.074, and the weakest month inside the strongest two-month run is 2.948 — so the word would flip only if the floor were lowered by 0.552. Persistence is not the lever: even at one month no held-out value reaches **3.5**, so the policy would still be silent. `alert_enabled` is `False` in the pin, so no rule could have produced ALERT here in any case.

On the other side, the trajectory word moves with the persistence rule, not with the envelope: the runs are `R_B` 6, `R_G` 2, `gamma2` 10, `P` 1, `A` 10, `D` 0, `F` 9, `S` 1, so any persistence up to **6** keeps CHANGE — a block run and a dimension run both clear it — and at **7** the longest block run is exhausted and the word falls to WATCH on dimensions alone. What would *not* move it is a slightly looser envelope: `gamma2` sits outside q05/q95 in 11 of 12 held-out months, so widening the reference spread by a little leaves the same 10-month run.

## Q13 · What does this record not answer?

* **The observability level.** It carries chips, not `(month, level)` records, so `observability_drop` (0.5) cannot be read and one of the policy's three firing rules is absent. The policy's word restricted to the rules this cut carries is `NOMINAL`; over every firing rule it names it is `INSUFFICIENT`. Neither is a verdict on the asset.
* **Four of the policy's seven words.** `DEGRADED_OBSERVABILITY`, `UNRESOLVED`, `CORROBORATED_WATCH` and `ALERT` are unreachable in an observability cut that carries no observability level and no second stream, so a reader should not read their absence as clearance. This cut prints only `NOMINAL`, `WATCH` and `INSUFFICIENT`.
* **ALERT, in particular.** `alert_enabled` is `False` in the pinned policy, so that word could not have been produced here whatever the numbers: the policy requires a human to escalate.
* **Causation.** Neither reader names a cause, and no number here is evidence that anything caused the collapse. The record says a monthly state departed from a frozen reference; it does not say why.
* **Sub-monthly timing.** Every state is a calendar-month median, so a lead counted in months is a count of bins, not of days: a departure that began on the 25th of a month reads as that whole month. The lead figures in Q5 are anchored to the bin's tick for exactly this reason.
* **A population.** One case and one control is not a sample. The quantiles are empirical quantiles of 12 reference months, so each is positioned by one or two observations; the envelope is a description of those months, not an estimated population distribution.
* **The nine post-event chips.** They are kept in the cut and named in the partition, and they are read by nothing except the hindsight guard, which deletes them and shows that the answer does not move.

## The pair, in one paragraph

With the reference frozen at 2016-09 → 2017-08 (109 chips, selected by a rule that cannot see the event) and 12 bins held out before the event, the trajectory rule reads **CHANGE**: `gamma2`, `A` and `F` leave the reference envelope for 10, 10 and 9 held-out months respectively, and the block runs rise with them, from 2017-12. EAS policy v1.0.0 reads **NOMINAL**, because no held-out value reaches its block or dimension floor — the closest run would need a floor of 2.948 against the pinned 3.5. The distance between the two words is that margin, and it is the finding: a departure that is unambiguous against the reference months' own spread is invisible against a fixed absolute-scale floor of the same persistence rule. Both words are printed; neither is retracted and neither is preferred.

## Reproducing this file

```
python3 code/cases/ocv_prospective.py            # write the artefact, the report, the plate
python3 code/cases/ocv_prospective.py --verify   # rebuild and compare with the committed files
python3 code/cases/ocv_prospective.py --selftest # the pinned values, the guards, the vocabulary; writes nothing
```

The artefact beside this page is `data/morandi/fig_morandi_ocv_trajectory_prospective.json`; the comparison is `data/cases/fig_ocv_prospective_comparison.json`. Every number in every answer above is read from that JSON, and `--verify` rebuilds it from the chips and refuses to agree with a stale file.

`schema` **ocv-prospective-v1** · pin `data/cases/eas_policy_pin.json` sha256 `e0557e90cbad29f496f483720fa1ce63f96765cd5c58a3c210be59639d0968d8` (3 digests verified) · policy digests policy_v1.json, rules.py, metric.py · guards all pass
