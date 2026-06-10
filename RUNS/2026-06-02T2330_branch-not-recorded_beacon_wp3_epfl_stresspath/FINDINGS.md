# EPFL Task-3.3 stress-path suite on the Maxwell-conjugate binary — findings (2026-06-02)

**Request:** Vinay, "run ... EPFL test (the stress path dependency), do this first."
Binary chosen: **Maxwell-conjugate only** (`maxwell-conjugate-20260602`, full-Π
closure + the Maxwell-conjugate micro-potential term, residual-only; analytic
Jacobian block not yet wired).

**Port:** inputs copied from `beacon_wp3_epfl_repro_2026-05-31/model` into
`./model`; the only edit was stripping the retired `<accumulate_swelling_contributions>`
tag (switch retired 2026-06-01 — the full-Π closure accumulates the swelling
eigenstress by default). **No parameter changed** — K = 44200 J/kg (Dixon 2023
MX-80, ρ_d = 1500, cited in-PRJ) and all vdW/micro params carry over. Run via
`run_epfl_maxwell.sh`; statuses in `results/EPFL_MAXWELL_SUMMARY.txt`.

## Run outcomes (4/5 complete)

| PRJ | result | reached | vtus |
|---|---|---|---|
| `column_le` | ✅ completed | 200 d | 6 |
| `path1_P1-3_dsm_mcc` | ❌ diverged | 918 685 s (~10.6 d), step 62, min-dt floor | 1 (t=0) |
| `path2_P2-1_dsm_mcc` | ✅ completed | 240 d | 20 |
| `path2_P2-1_swellingpressure_dsm_mcc` | ✅ completed | 240 d | 20 |
| `path2_compression_homogeneous_mcc` | ✅ completed | 40 d | 9 |

## Headline — the gate as written is inert on the real EPFL stress path

On `path2_P2-1` (a wetting+loading path: n_l 1.5e-4 → 0.332, S 0→1, p′ rising to 10.4 MPa):

| quantity | value |
|---|---|
| max p′ = −tr(σ′)/3 | **10.38 MPa** |
| Π = micro_pressure (the coded gate threshold) | 22–28 MPa |
| φ_m·Π = micro_porosity·Π (REV-consistent threshold, #4 option 1) | 9.2–10.2 MPa |
| **coded gate** p′ ≥ Π | **never fires** (10.4 < 27.7) → term inert |
| **REV-consistent gate** p′ ≥ φ_m·Π | **fires** near peak load (10.4 > 9.2) |

This is the **third independent confirmation** (after the confined oedometer probe
and the gate-scale analysis) that the coded gate `p′(REV) ≥ Π(micro)` mixes scales
and is **never crossed on realistic paths**, while the REV-consistent `p′ ≥ φ_m·Π`
activates exactly where physics expects (late compression, micro hydrated). Strong
empirical support for **gate-scale decision #4, option 1**. *(Observed, not yet a
re-run with the rescaled gate — that is the next step once #4 is settled.)*

## path1 divergence — attributed, NOT the Maxwell term

Control: `path1_P1-3` re-run on `dsm_native_hier_wt-release` (same full-Π closure,
**no** Maxwell term) diverges at the **bit-identical** Time 918684.60940596461 s.
⇒ the divergence is intrinsic to the full-Π closure + MFront-MCC on the P1-3 path
(a sharp hydraulic/wetting-front transition; cf. DSM/AGENTS.md "Known Limitation —
Hydraulic Side of OGS RM-DSM"), independent of the Maxwell-conjugate work. Consistent
with the term being inert there (gate firmly closed: Π ≈ 203 MPa at the dry initial state).

## Provenance
- binary: `/Users/vinaykumar/git/build/maxwell-conjugate-20260602/bin/ogs`
- control: `/Users/vinaykumar/git/build/dsm_native_hier_wt-release/bin/ogs`
- inputs: `./model` (ported copy, dead key stripped, params unchanged)
- OMP_NUM_THREADS=4
