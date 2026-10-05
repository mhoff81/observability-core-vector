# `code/` — the two site packages

One subfolder per site, both standalone, offline and built the same way (stats
layer, mask layer, data layer, generator, figure script). They differ in exactly
one place: the mask layer. LUMO masks the **echo** of the tower, Bautzen masks
the **deck edge** of the bridge.

| folder | site | states | mask layer | package README |
| --- | --- | --- | --- | --- |
| [`lumo/`](lumo/README.md) | LUMO lattice tower, Hannover | healthy, DAM3, DAM4, DAM6 | echo mask (`0.30 x peak`, `5 x median`, ≥ 2 px) | `lumo/README.md` |
| [`bautzen/`](bautzen/README.md) | Bautzen OpenLabs bridge | RS, DS1, DS2 | deck-edge mask (geodetically anchored deck band) | `bautzen/README.md` |

The scripts are started from the **repository root** (each wrapper in
`../figures/<site>/` does that), so that their defaults resolve to
`data/<site>/…` in and `data/<site>/…`, `figures/<site>/…` out. No package:
the script's own folder is on `sys.path` (a script run puts its directory first),
which is why the modules import each other by their `lumo_…` / `bautzen_…` names.

Adding another site means copying one package and swapping `…_ocv_masks.py` plus
the `DATA`/`FIGDIR` lines of `…_ocv_core.py` — nothing outside the new
`code/<site>/`, `data/<site>/` and `figures/<site>/` changes.
