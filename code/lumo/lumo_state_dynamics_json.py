#!/usr/bin/env python3
"""LUMO state dynamics — builder and verifier of the reference JSON.

``data/lumo/reference/lumo_state_dynamics.json`` is the *physical* counterpart of
``lumo_state_frequencies.json``: for each damage state it stores the relative
change of the two quantities that make up the physical response vector

    v_phys = [df_H1 / f_H1,0 , dzeta_H1 / zeta_H1,0]

i.e. the **relative modal-frequency change** (signed, from the Welch-PSD mode
consensus in ``lumo_state_frequencies.json``) and the **relative damping
change** of the same mode (H1, half-power bandwidth).  The fundamental mode is
carried along as a sensitivity, never as the primary axis.

Where the numbers come from.  The frequency half is already committed; the
damping half lives outside this repository in the LUMO SHM project.  This builder
therefore reads three sources and pins all of them:

  ``data/lumo/reference/lumo_state_frequencies.json``   committed, H1 + fundamental
  ``/home/projects/lumo_shm_data/lumo_damping_ratio.json``
      per-recording half-power bandwidth (30 recordings: 6 campaign blocks x 5)
  ``/home/projects/lumo_dam6_analysis/lumo_damping_gamma2_modulation.json``
      the per-state damping summary, used as an independent cross-check

Each state was measured in **one** controlled campaign: five healthy recordings
followed by five damaged recordings, and the healthy reference is the healthy
block of the *same* campaign (the tower as it stood days before that state).  The
stored ``zeta_*_recordings`` lists are the raw per-recording values, so the
medians, the min/max range and any bootstrap interval can be recomputed from the
committed JSON alone.

Usage:
  python3 code/lumo/lumo_state_dynamics_json.py --build   # rebuild the reference
  python3 code/lumo/lumo_state_dynamics_json.py           # verify only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import lumo_ocv_core as core    # noqa: E402

SHM_ROOT = "/home/projects/lumo_shm_data"
DAMPING_JSON = os.path.join(SHM_ROOT, "lumo_damping_ratio.json")
MODULATION_JSON = os.path.join("/home/projects/lumo_dam6_analysis",
                               "lumo_damping_gamma2_modulation.json")
FREQ_JSON = os.path.join(core.REF, "lumo_state_frequencies.json")
OUT_JSON = os.path.join(core.REF, "lumo_state_dynamics.json")

SCHEMA = "lumo-state-dynamics-v1"

# The damping campaign of a state and the mode name used by the damping
# analysis.  The *frequencies* reference names the modes
# ``first_higher`` / ``fundamental``; the damping source keys them
# ``H1`` / ``fundamental``.  Both spellings are pinned below.
DAMPING_KEY = {"DAM 3": "dam3", "DAM 4": "dam4", "DAM 6": "dam6"}
DAMPING_MODE = {"first_higher": "H1", "fundamental": "fundamental",
                "second_higher": "H2"}
# Which mode carries the analysis: H1 is the primary axis (it supplies df),
# the fundamental mode is the sensitivity, second_higher is context only.
ROLE = {"first_higher": "primary", "fundamental": "sensitivity",
        "second_higher": "context"}
DAMPING_FIELD = {"H1": "zeta_h1", "fundamental": "zeta_fundamental",
                 "H2": "zeta_h2"}

METHOD = (
    "Frequency: Welch PSD per accelerometer channel (4 s windows), averaged "
    "over the 18 acceleration channels, peaks in the 1-50 Hz band, mode "
    "consensus across the five recordings of a campaign (+-0.75 Hz clustering). "
    "Damping: half-power bandwidth on the channel-mean Welch PSD (nperseg 32 s, "
    "parabolic peak interpolation, SNR gate 3x, state-adaptive mode windows), "
    "one value per recording, median over the five recordings of a campaign. "
    "The healthy reference of a state is the healthy block of the same "
    "campaign, so both halves of v_phys refer to the tower as it stood days "
    "before that state."
)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _records(path):
    """Per campaign block: the five H1 / fundamental damping values."""
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    out = {}
    for r in doc["recordings"]:
        key = (r["campaign"], r["state"])
        out.setdefault(key, []).append({
            "file": r["file"],
            "zeta_h1": float(r["H1"]["zeta"]),
            "zeta_fundamental": float(r["fundamental"]["zeta"]),
            "zeta_h2": float(r["H2"]["zeta"]),
        })
    return doc, out


def build():
    with open(FREQ_JSON, encoding="utf-8") as fh:
        freq = json.load(fh)
    with open(MODULATION_JSON, encoding="utf-8") as fh:
        modulation = json.load(fh)
    damping_doc, blocks = _records(DAMPING_JSON)
    mod_by_state = {r["state"]: r for r in modulation["rows"]}

    bin_hz = float(freq["params"]["fft_bin_hz"])
    states = []
    for camp in freq["campaigns"]:
        state = camp["state"]
        key = DAMPING_KEY[state]
        healthy = blocks[(key, "healthy")]
        damaged = blocks[(key, key)]
        row = mod_by_state[key]
        modes = []
        for mode in camp["modes"]:
            name = mode["name"]
            dkey = DAMPING_MODE[name]
            field = DAMPING_FIELD[dkey]
            zh = [r[field] for r in healthy]
            zd = [r[field] for r in damaged]
            med_h = st.median(zh)
            med_d = st.median(zd)
            modes.append({
                "name": name,
                "damping_key": dkey,
                "role": ROLE[name],
                "healthy_hz": mode["healthy_hz"],
                "state_hz": mode["state_hz"],
                "healthy_n": mode["healthy_n"],
                "state_n": mode["state_n"],
                "df_rel": (mode["state_hz"] - mode["healthy_hz"])
                          / mode["healthy_hz"],
                "df_bin_uncertainty_rel": bin_hz / 2.0 / mode["healthy_hz"],
                "zeta_healthy": med_h,
                "zeta_state": med_d,
                "zeta_rel": (med_d - med_h) / med_h,
                "zeta_healthy_recordings": zh,
                "zeta_state_recordings": zd,
                "zeta_healthy_files": [r["file"] for r in healthy],
                "zeta_state_files": [r["file"] for r in damaged],
            })
        states.append({
            "state": state,
            "healthy_campaign": row["healthy_campaign"],
            "damaged_campaign": row["damaged_campaign"],
            "healthy_window": camp["healthy_window"],
            "state_window": camp["state_window"],
            "modes": modes,
        })

    doc = {
        "schema": SCHEMA,
        "method": METHOD,
        "params": {
            "fft_bin_hz": bin_hz,
            "h1_bin_uncertainty_hz": bin_hz / 2.0,
            "n_recordings_per_block": 5,
            "n_recordings_total": len(damping_doc["recordings"]),
            "primary_mode": "first_higher",
            "sensitivity_mode": "fundamental",
            "rel_reference": "healthy block of the same campaign",
        },
        "provenance": {
            "sources": {
                "frequencies": {"path": FREQ_JSON, "sha256": sha256(FREQ_JSON)},
                "damping_recordings": {"path": DAMPING_JSON,
                                       "sha256": sha256(DAMPING_JSON)},
                "damping_summary": {"path": MODULATION_JSON,
                                    "sha256": sha256(MODULATION_JSON)},
            },
            "note": ("The damping summaries live outside this repository; the "
                     "per-recording values are copied verbatim so that this "
                     "JSON is self-contained.  `--verify` recomputes every "
                     "derived number from these lists."),
        },
        "prepared_by": "code/lumo/lumo_state_dynamics_json.py",
        "states": states,
    }
    return doc, modulation, damping_doc


def verify(doc=None, freq_path=FREQ_JSON):
    """Recompute every derived number from the committed file alone."""
    if doc is None:
        with open(OUT_JSON, encoding="utf-8") as fh:
            doc = json.load(fh)
    checks = []

    def check(name, got, want):
        entry = {"check": name, "ok": bool(got == want)}
        if not entry["ok"]:
            entry["recomputed"] = repr(got)[:300]
            entry["committed"] = repr(want)[:300]
        checks.append(entry)

    check("schema", doc["schema"], SCHEMA)
    with open(freq_path, encoding="utf-8") as fh:
        freq = json.load(fh)
    freq_by = {c["state"]: c for c in freq["campaigns"]}
    check("n_states", len(doc["states"]), 3)
    check("state_order", [s["state"] for s in doc["states"]],
          ["DAM 3", "DAM 4", "DAM 6"])
    check("primary_mode", doc["params"]["primary_mode"], "first_higher")
    check("sensitivity_mode", doc["params"]["sensitivity_mode"], "fundamental")

    for s in doc["states"]:
        state = s["state"]
        check(f"{state}.n_modes", len(s["modes"]), 3)
        for mode in s["modes"]:
            tag = f"{state}.{mode['name']}"
            zh = mode["zeta_healthy_recordings"]
            zd = mode["zeta_state_recordings"]
            check(f"{tag}.n_healthy", len(zh), 5)
            check(f"{tag}.n_state", len(zd), 5)
            check(f"{tag}.median_healthy", st.median(zh), mode["zeta_healthy"])
            check(f"{tag}.median_state", st.median(zd), mode["zeta_state"])
            check(f"{tag}.zeta_rel",
                  (st.median(zd) - st.median(zh)) / st.median(zh),
                  mode["zeta_rel"])
            check(f"{tag}.df_rel",
                  (mode["state_hz"] - mode["healthy_hz"]) / mode["healthy_hz"],
                  mode["df_rel"])
            check(f"{tag}.healthy_in_range",
                  min(zh) <= mode["zeta_healthy"] <= max(zh), True)
            check(f"{tag}.state_in_range",
                  min(zd) <= mode["zeta_state"] <= max(zd), True)
            ref = {m["name"]: m for m in freq_by[state]["modes"]}[mode["name"]]
            check(f"{tag}.freq_reference_healthy_hz",
                  mode["healthy_hz"], ref["healthy_hz"])
            check(f"{tag}.freq_reference_state_hz",
                  mode["state_hz"], ref["state_hz"])
            check(f"{tag}.freq_reference_shift_pct",
                  round(mode["df_rel"] * 100.0, 2), ref["shift_pct"])
    n_ok = sum(1 for c in checks if c["ok"])
    return {"checks": checks, "n_checks": len(checks),
            "n_failed": len(checks) - n_ok}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build or verify the LUMO state-dynamics reference.")
    ap.add_argument("--build", action="store_true",
                    help="read the three sources and rewrite the reference JSON")
    ap.add_argument("--out", default=OUT_JSON, help="path of the reference JSON")
    args = ap.parse_args(argv)

    if args.build:
        doc, modulation, _damping = build()
        mod_by_state = {r["state"]: r for r in modulation["rows"]}
        for s in doc["states"]:
            row = mod_by_state[DAMPING_KEY[s["state"]]]
            h1 = [m for m in s["modes"] if m["name"] == "first_higher"][0]
            fund = [m for m in s["modes"] if m["name"] == "fundamental"][0]
            pairs = ((h1["zeta_healthy"], row["zeta_H1_healthy"], "H1 healthy"),
                     (h1["zeta_state"], row["zeta_H1_damaged"], "H1 damaged"),
                     (fund["zeta_healthy"], row["zeta_fund_healthy"],
                      "fundamental healthy"),
                     (fund["zeta_state"], row["zeta_fund_damaged"],
                      "fundamental damaged"))
            for got, want, tag in pairs:
                if abs(got - want) > 1e-12:
                    raise SystemExit(
                        f"median mismatch for {s['state']} {tag}: "
                        f"{got!r} != {want!r}")
            if s["healthy_campaign"] != row["healthy_campaign"]:
                raise SystemExit(f"campaign mismatch for {s['state']}")
            if s["damaged_campaign"] != row["damaged_campaign"]:
                raise SystemExit(f"campaign mismatch for {s['state']}")
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=1)
            fh.write("\n")
        print(f"wrote {os.path.relpath(args.out)}")
        print(f"  {len(doc['states'])} states, "
              f"{doc['params']['n_recordings_total']} recordings")

    res = verify()
    print(f"verify: {res['n_checks']} checks, {res['n_failed']} failed")
    for c in res["checks"]:
        if not c["ok"]:
            print(f"  FAIL {c['check']}: {c.get('recomputed')} != "
                  f"{c.get('committed')}")
    return 0 if res["n_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
