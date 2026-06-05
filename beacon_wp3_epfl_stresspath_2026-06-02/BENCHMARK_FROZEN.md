# FROZEN benchmark — current_porosity_split MX-80 calibration + EPFL Task-3.3 verification

**Branch `dsm_cps_mx80`, frozen 2026-06-05.** Disjoining switched to the
energy-consistent micro-only mode `micro_solid_volume_fraction_mode =
current_porosity_split` (nS = 1 − n_l). Reference results below are the freeze.

## Calibration — MX-80 swelling pressure (ANCHORS MS33 Model I, 1-element const. vol.)
Fit K to Dixon (2023) MX-80 σ_swell = 0.003·exp(5.2883·ρ_d) under the new mode:

| ρ_d (kg/m³) | K [J/kg] | model Pₛ | Dixon target | err |
|---|---|---|---|---|
| 1400 | 36795.5 | 4.91 MPa | 4.92 | −0.31% |
| 1600 | 85775.6 | 14.16 | 14.16 | −0.01% |
| 1800 | 227262  | 40.86 | 40.86 | +0.00% |

K barely shifts from the `reference`-mode values (end-state mode-coincidence: at
saturation n_l→φ0 ⇒ nS=1−n_l→1−φ0). Needs `initial_dt`≈0.5 (split mode is stiff at the
S_L→1 transition; the dnS/dn_l local-solve tangent is not yet wired — interim).

## Verification — EPFL Task-3.3 granular MX-80, e–σ_v vs Ferrari et al. (2022)
`figures/epfl_cps_both_paths.png`. K(1500)=56183 (log-interp of the new anchors).
- **P2 (const-volume):** model end e=0.50 @16.5 MPa vs measured C′ 0.57 — close (with a
  current_porosity_split swelling-pressure overshoot in the dry-wetting transient).
- **P1 free-swell:** capped e≈1.0 — the measured B=2.34 (~82% strain) is beyond OGS RM
  small-strain; both E-softening and K-cranking diverge before ~e1.3. NOT reproduced.
- **P1 compression:** MCC run from the IMPOSED measured B=2.34; λ_mcc recalibrated
  0.077→0.31 (fit to D; pc/yield-stress is a weak lever) → endpoint e=0.56 = measured D.
  Intermediate stays above measured C (large-strain MCC shape ≠ log-linear).

## Commits (this freeze)
72b4c78732 recalib · e43e5c0b8a EPFL redo · 184f6dae17 P1 imposed-B · bd945ed471 λ-fit.
Built on `dsm_native_pdisj_maxwell` (f7050d6319 ∂ρ̂/∂ε tangent, e32fde2a4d CTSF rejected).

## Known limitations (carried, not bugs)
1. dnS/dn_l local-solve tangent missing → small dt needed (split mode).
2. Free-swell magnitude needs finite-strain (small-strain caps ~e1.3).
3. λ=0.31 is a fit-to-D in the large-strain MCC, not the textbook index.
4. EPFL K is interpolated from the LE anchors (LE→MCC + block→granular gaps, as before).
