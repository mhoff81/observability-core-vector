# `code/ywf/` — Yeongdeok Wind Farm, the collapsed Unit 21 monopole

**Site** Yeongdeok Wind Farm / Changpo Wind Power Complex, Samgye-ri
(36.4236 N 129.4203 E) · **asset** `Yeongdeok Wind Turbine (Samgye-ri)`
(`af683c7e-416a-4395-86cc-094ee83d9497`) · **event** the tower **collapsed onto a
public road on 2026-02-02** · **states** `pre-collapse (healthy)` (Aug 2025 -
Jan 2026) and `post-collapse` (Feb - Apr 2026) of the *same* structure.

This is the seventh package of the family (after `lumo/`, `bautzen/`, `kdlo/`,
`carola/`, `morandi/`, `cts/`) and it keeps their shape: a stats layer, a data
layer, a mask layer, a channel-table generator, a figure/report script with a
pin block, and a committed reference under `data/ywf/reference/`. What differs is
the *size* of the record and, with it, what the package is allowed to claim.

```bash
python3 code/ywf/ywf_ocv_masks.py --check     # re-verify the mask rule on the cache
python3 code/ywf/ywf_ocv_masks.py --extract   # rewrite data/ywf/ywf_windows_mask_cache.txt
python3 code/ywf/ywf_ocv_channels_csv.py      # build + verify the channel table
python3 code/ywf/ywf_ocv_core.py              # read-only brief of the committed layer
python3 code/ywf/ywf_ocv_core.py --chips      # one line per unique chip, with its 6D vector
bash    figures/ywf/fig_ywf_ocv_paper.sh      # figures A-E + report + pins (~8 s)
bash    figures/ywf/fig_ywf_ocv_paper.sh --write-reference   # (re)define the reference
```

## The record

| fact | value |
| --- | --- |
| committed windows | **57** on **32** acquisition dates (`data/ywf/ywf_windows_full.txt`, 38 kB, committed **in full**) |
| unique chips | **33** — 24 of the 57 windows are **byte-identical duplicates** that differ only in the requested mast section (index 2 vs. 3, same date, same payload) |
| `segment_resolved` | **false** — the two requests decoded the *same* chip at the asset point |
| mast sections present | `turbine-mast-section-3` only (of the five in `data/ywf/ywf_segments.txt`) |
| chips with an echo | **17** of 33 — **16** of 27 pre-event chips, **1** of 6 post-event chips |
| window size | 7 x 7 px (49 px), one complex window per request |

Because the record is not segment-resolved, this package carries **one** mask
layer — the echo mask of the whole 7 x 7 window. A per-mast-section (per-girder)
layer, the second layer of Carola/Morandi/CTS, would be an artefact of a
duplication and is deliberately absent; figure E d shows the duplication.

## The mask layer

`data/ywf/ywf_windows_full.txt` is the committed payload (this site is small
enough to commit it), and `data/ywf/ywf_windows_mask_cache.txt` is the committed
mask cache the other packages have, written by `ywf_ocv_masks.py --extract`:

```
intensity = |z|^2
mask      = intensity >= 0.30 * max(intensity)  AND  intensity >= 5.0 * np.median(intensity)
A         = count(mask)                    (no echo when A < 2)
gamma2    = |sum(z[mask])|^2 / (A * sum(|z[mask]|^2)),  clamped to [0, 1]
D         = A / bbox_area          F = 8-connected components      S = ||centroid - peak||
P         = mean pixel frequency inside the own mask
```

The cache holds **only the 17 windows that carry an echo** (8 kB); its manifest
records the 16 dropped ones (`n_without_echo`), the 24 de-duplicated payloads and
the section counts. `--check` re-derives every dimension of every cache row from
the committed payload and fails on the first mismatch.

## The channel table

`data/ywf/ywf_ocv_channels.csv` — **one row per unique chip** (33 rows, 81
columns, built by `ywf_ocv_channels_csv.py`, which **refuses to write** unless the
mask layer verifies). `echo_mask_present` marks the 17 echo-bearing rows; on the
other 16 the six mask dimensions are empty, because a mask vector without a mask
does not exist:

| group | columns |
| --- | --- |
| provenance | `id`, `asset_*`, `request_id`, `segment_index`, `acquisition_ts`, `pass_label`, `orbit_direction`, `status`, `…` (34 source columns) |
| derived | `segment`, `section_index`, `state`, `pre_collapse`, `day`, `month`, `year`, `season`, `echo_mode`, `echo_mask_present`, `payload_group_size`, `payload_shared_with_segments` |
| recomputed mask | `A`, `D`, `F`, `S`, `P`, `gamma2_mask`, `mask_bbox_area`, `mask_row_span`, `mask_col_span`, `mask_centroid_row/col`, `mask_peak_row/col`, `mask_n_components`, `mask_largest_component`, `mask_fragmentation` |
| pipeline column | `coherence_masked_pixels`, `coherence_gamma2`, `intensity`, `brightness_ratio`, … (kept as provenance, not as the mask) |

`echo_mode` is cut on this site's own scale (`A <= 2` compact, `3..5`
intermediate, `>= 6` distributed of 49 px) — the project-wide cuts belong to the
80 x 80 bridge chips and would put every YWF echo in one bucket. The cut used is
recorded in the metadata (`echo_mode_cuts`).

### `A` is not `coherence_masked_pixels`

The pipeline's own column was computed on the **pre-fix asset-point chip of the
whole-scene window** (75, 133, … px) while this package recomputes the echo mask
of the committed 7 x 7 window (2..7 px). The two agree on **0 of 17** rows
(`n_coherence_masked_pixels_matched` in the metadata, panel E a), so the CSV keeps
both and reports the coincidence instead of merging them.

## The reference file — defined here, then pinned

The other sites pin against JSONs their upstream analysis project committed.
**This site has no upstream analysis project** (the record *is* the site's own
export), so the reference layer is defined by this package:

* `fig_ywf_ocv_paper.py --write-reference` writes
  `data/ywf/reference/ywf_ocv_findings.json` and proves the round trip (the file is
  re-read and re-compared field by field);
* every later run **without** the flag re-derives the whole file from the
  committed tables and aborts (exit code != 0) on the first deviation.

The reference carries no timestamp and is written with sorted keys, so two runs
are byte-identical; the mask cache, the window payload, the measurements, the
segments, the manifest and the channel CSV are pinned by **sha256** as well. A
`nan` — a test that cannot be computed, here the Welch/Mann-Whitney pair of a
control whose two states are nearly identical (`coherence` = 0.999 throughout) —
is stored as JSON `null`, never as the non-standard `NaN` literal.

## What the figure shows (and what it does not)

| panel | content |
| --- | --- |
| `A` | the six echo-mask channels per state, over the **17 echo-bearing chips** — a box where a state has ≥ 3 of them, a single marker for the one post-event chip |
| `B` | the **headline**: (a) echo coverage per state with its Fisher exact test; (b) the coverage timeline across the collapse (marker area = unique chips of the day) |
| `C` | in-sample effect sizes: (a) Cliff's delta + bootstrap CI95 and SDS of the six mask channels; (b) the same for the co-variate controls |
| `D` | the OCV plane (`gamma2` vs `P`, `gamma2` vs `D`, `P` vs `D`), 16 pre chips and 1 post chip |
| `E` | why the contrast is not a damage signal: (a) the export `A` vs. the recomputed `A`; (b) the echo-mode mix; (c) the season composition; (d) the 24 duplicated windows |

**No classifier is fitted.** Panel D is a scatter, not a decision boundary: one
post-event echo chip is not a class, so a leave-one-out accuracy computed from it
would be a property of the sample size, not of the tower. The report says so
instead of printing a number.

**The site's claim is an echo-coverage statement.** 16 of 27 pre-event chips
carry an echo, 1 of 6 post-event chips does; the Fisher exact test gives
p = 0.085 — not significant, and it cannot be, because that single surviving echo
is in the post side. Had all six post-event chips been echo-free, four chips
would already have sufficed at alpha = 0.05
(`echo_coverage.power.post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha` in
the reference). What the record supports is the *observability* statement: the
echo mask is the quantity that can be observed on this collapsed tower, and its
coverage is the number to extend before any state contrast is attempted — never a
damage verdict.

## The modules

| file | role |
| --- | --- |
| `ywf_ocv_stats.py` | the statistics of the family (Cliff's delta + bootstrap CI, SDS, Welch/Mann-Whitney, sign test, Spearman/Mann-Whitney, shrinkage LDA + LOO + label permutation, block-permutation correlation, OLS) — ported verbatim, unused parts included so the seven packages read as one diff |
| `ywf_ocv_masks.py` | the mask rule, the payload reader, the cache writer (`--extract`), the acceptance test (`--check`), the 6D vector of a chip |
| `ywf_ocv_core.py` | paths, the site/event/state constants, the state and echo-mode definitions, the loaders of the committed tables, `chip_table`, `coverage_block`, `fisher_p`, and a read-only CLI brief |
| `ywf_ocv_channels_csv.py` | the channel-table generator (build + verify) |
| `fig_ywf_ocv_paper.py` | panels A–E, the result JSON (`data/ywf/fig_ywf_ocv_paper.json`), the report (`figures/ywf/fig_ywf_ocv_paper.md`), the reference writer and the pin block |

## Pins

Every run prints its pin summary and exits non-zero on any deviation (1080
checks on the current record):

| pinned against | what |
| --- | --- |
| `data/ywf/reference/ywf_ocv_findings.json` | the echo-coverage table and its Fisher test, the mask dimensions per state, the pipeline column, the co-variate controls, the month/section composition, the 32-row timeline and the payload-digest groups — field by field |
| `data/ywf/ywf_ocv_channels_meta.json` | the table's own definition, guard counts and per-state summaries against this package's recomputation |
| `data/ywf/ywf_windows_mask_cache.txt` | each cached chip's `n_masked` / `n_components` against the channel table's `A` / `F` (17 rows) |
| the committed inputs | their sha256 digests (payload, cache, manifest, measurements, segments) |

Two runs of `fig_ywf_ocv_paper.py` with the fixed seed are byte-identical — the
result JSON, the report and all five PNGs (checked with `sha256sum`, including a
run into a scratch directory); no wall-clock value is written and the runtime
stays on stdout. `--write-reference` is idempotent: re-writing the reference
leaves the committed `data/ywf/reference/ywf_ocv_findings.json` byte-for-byte
unchanged.

## Reading

Before the collapse **16 of 27** unique chips carry an echo, after it **1 of 6**.
With six post-event acquisitions and one surviving echo the site cannot separate
the states, and the single post-event echo chip lies inside the pre-event cloud of
`gamma2` / `P` / `D` (panel D) — so the effect sizes of panel C are reported to
*document* the co-variation of season, weather and acquisition generation with the
state, not as a mechanism. The honest output of this package is the coverage table
of panel B and the power statement that goes with it.

