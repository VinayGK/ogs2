# EPFL Task-3.3 — term-ON (p'≥φ_m·Π) vs term-OFF baseline (2026-06-03)

**Setup.** Decision #4 = option 1: the Maxwell-conjugate gate threshold rescaled
from the intrinsic micro Π to the REV-consistent partial stress φ_m·Π = n_S·n_l·Π
(commit `d542d6c08c`; trigger-only, S₁ keeps the full Π). Scaled binary rebuilt
(`maxwell-conjugate-20260602`, 23:56). Same ported PRJs, new folder.

- **term-OFF baseline** (coded gate p'≥Π, never opened): `beacon_wp3_epfl_stresspath_2026-06-02`
- **term-ON** (this folder, gate p'≥φ_m·Π)

## Result: term-ON ≡ term-OFF, bit-for-bit

| run | max |Δ| over ALL fields, ALL steps |
|---|---|
| `path2_P2-1_dsm_mcc` | **0.000e+00** |
| `path2_P2-1_swellingpressure_dsm_mcc` | **0.000e+00** |

`path1_P1-3_dsm_mcc` still diverges at the identical time (918 685 s) — the
gate-independent wetting-front closure limit. column_le / LE / single-structure
MCC carry no DSM term. So the e–σ figure is **identical** to the baseline
(`epfl_both_paths_termON_IDENTICAL.png`).

## Why the rescaled term still does not bite (see `gate_trajectory_path2.png`)

`p'` **hugs the φ_m·Π operating point** along the whole path. The REV gate is open
only where the term cannot act:

| window | p' vs φ_m·Π | why no effect |
|---|---|---|
| start (t=0–10 d) | open (p' 0.18 > φ_m·Π 0.02) | ε_v ≈ 0, and the term ∝ ε_v ⇒ 0 |
| swell + compression (t=16–230 d) | **shut** (p' 5–9 < φ_m·Π 9.5–10.2) | gate closed |
| end (t=235–240 d) | open, marginal (p' 9.7–10.4 vs φ_m·Π 9.4–9.2) | only the last 1–2 steps; the explicit gate (reads the previous converged confining stress) lags below threshold |

This **empirically confirms the note's §6 prediction**: because φ_m·Π = |σ_sw| ≈
the swelling pressure, p' equilibrates *at* the gate rather than crossing it
decisively. The rescale (option 1) is correct and necessary, but the EPFL paths
**do not exercise the term** — they sit at the operating point.

## What it would take to demonstrate the term firing (not yet done)
1. A path loaded **well past** the swelling pressure (φ_m·Π ≈ 8–10 MPa here) with
   live strain — e.g. a higher-load confined oedometer probe (above-gate
   verification, task #3). *BC magnitude is Vinay's call.*
2. And/or fix the **explicit-lag gate** (evaluate on the current iterate, not the
   previous converged stress) so a marginal late crossing is caught.
3. And/or fix `path1_P1-3` convergence (the wetting-front limit), since P1-3
   reaches higher void ratio / larger strains and would load deeper past φ_m·Π.

## Provenance
- scaled binary: `/Users/vinaykumar/git/build/maxwell-conjugate-20260602/bin/ogs` (commit `d542d6c08c`)
- term-off control = the kept baseline folder (coded-gate binary, same source pre-`d542d6c08c`)
- OMP_NUM_THREADS=4
