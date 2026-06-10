# BEACON WP5.3 MGR tests (MX-80 pellet+block homogenisation) — provenance extraction

Source deliverables (texts symlinked from ../../beacon_wp3_epfl_repro_2026-05-31/spec/):
- **D5.7** (BEACON Final report, 15/01/2022) §5.3 — UPC/CODE_BRIGHT BExM modelling of MGR22/23/27.
- **D5.6** (WP5 Step-3 specifications, 30/04/2021) §2.2, §3, §4 — test specs + multi-team synthesis.
- Experimental data: **Villar et al. (2021)**; WP4 deliverable **D4.3** (D5.6.txt:432). CIEMAT.

> NOTE: locators are line numbers in D5.7.txt / D5.6.txt (the symlinked copies).
> Figure-only quantities are flagged; do not invent numerics for them.

## The test family
Three constant-volume hydration columns, MX-80, pellet layer + compacted block layer
stacked in an oedometric cell. Homogenisation (do the two densities converge?) is the
target. MGR22/23 = calibration (results known); **MGR27 = blind/predictive**, geometry
INVERTED (block at hydration face). Vinay decision 2026-06-01: do MGR22/23 first.

## 1. Geometry (D5.7.txt:8578-8592)
- Cylinder, diameter 5 cm, total height 10 cm; pellet layer ~5 cm + block layer ~5 cm.
- Axisymmetric, 400 linear quads / 451 nodes; effectively 1-D (no friction in UPC model).
- **MGR22/23:** pellets at BOTTOM (hydration face), block on top (D5.6.txt:509).
- **MGR27:** block at BOTTOM (hydration face), pellets on top — inverted (D5.6.txt:511-513).

## 2. Initial conditions — Table 5.3-2 (D5.7.txt:8777-8819)
| Test | layer | w[%] | h[cm] | rho_d[g/cm3] | Sr[%] | suction[MPa] |
|---|---|---|---|---|---|---|
| MGR22 | block | 13.60 | 4.94 | 1.61 | 55 | 120 |
| MGR22 | pellets | 9.90 | 5.04 | 1.28 | 25 | 115 |
| MGR23 | block | 14.20 | 4.98 | 1.60 | 56 | 120 |
| MGR23 | pellets | 3.50 | 5.00 | 1.30 | 9 | 290 (model used 120*) |
| MGR27 | block | 14.20 | 4.98 | 1.60 | 56 | 120 |
| MGR27 | pellets | 3.50 | 5.00 | 1.30 | 9 | 310 (model used 120*) |
- *Pellet initial suction reset to 120 MPa in the model (D5.7.txt:8775-8776): "effect of
  suction changes at very high values are slight." A modelling choice, flagged.
- Grain density rho_s = 2735 kg/m3 (D5.6.txt:1698-1699).

## 3. Boundary conditions (D5.7.txt:8828-8847)
- **MGR22:** constant water FLOW RATE 0.047 cm3/h at base, after 10 d dry period.
- **MGR23 / MGR27:** constant water PRESSURE 15 kPa at base, ramped over 1 day.
- Mechanical: top vertical displacement = 0 (axial load measured there); constant volume.
- Duration: ~200 d to saturation (MGR22, D5.7.txt:8900); others similar.

## 4. Measured validation targets
### Final state — Table 2-5 MGR22/23 (D5.6.txt:416-426), Table 2-6 MGR27 (D5.6.txt:450-460)
| Test | layer | w_final[%] | rho_d_final[g/cm3] | Sr[%] |
|---|---|---|---|---|
| MGR22 | pellets | 35.3 | 1.35 | 95 |
| MGR22 | block | 30.7 | 1.51 | 106 |
| MGR22 | total | 32.7 | 1.43 | 100 |
| MGR23 | pellets | 35.7 | 1.34 | 95 |
| MGR23 | block | 31.1 | 1.51 | 107 |
| MGR23 | total | 32.7 | 1.43 | 100 |
| MGR27 | pellets | 32.3 | 1.434 | 98.7 |
| MGR27 | block | 30.0 | 1.454 | 94.4 |
| MGR27 | total | 31.01 | 1.44 | 96.4 |
- **Homogenisation (the headline):** initial gap ~0.30-0.33 g/cm3 → final:
  MGR22 0.16 (~50% closure), MGR23 0.17 (~43%), **MGR27 0.02 (~93%, near-complete)**.
  MGR27 homogenises MORE despite being blind/inverted (D5.6.txt:7272-7273).
- Water intake: total ~150 cm3 all three (similar pore space) — D5.6.txt:7150-7151;
  time curves figure-only (D5.7 Figs 5.3-3b/8b/13b).
- Axial/swelling pressure: MGR22 final "correctly reproduced"; MGR23 model +17% over;
  MGR27 grossly over-predicted by model, attributed to neglected LATERAL FRICTION
  (D5.7.txt:9079-9084). Magnitudes figure-only.

## 5. Reference modelling parameters (UPC/CODE_BRIGHT BExM, Table 5.3-1, D5.7.txt:8637-8733)
*** §12.1 FLAG: this is a BExM (CODE_BRIGHT) parameter set, NOT the native Pi-path DSM.
    It may be cited as the Villar/EURAD-MS material-parameter SOURCE per group, but the
    swelling calibration K for a native-DSM build MUST come from Dixon/Villar (§12.1),
    and translating BExM macro/micro params into the Pi-path closure is a MODELLING
    decision for Vinay, not a transcription. ***
Key groups (block / pellets):
- elastic: nu 0.3/0.3; phi 25/25 deg; kappa 7e-4/1.5e-4; kappa_s 0.028/0.03.
- macro BBM: p_c 0.5/0.1 MPa; lambda(sat) 0.25; r 0.65/0.41; lambda 0.25/0.025; k_s 0.01;
  leakage 0.4e-5.
- micro WRC: P0_m 180 Pa; lambda_o 0.072; theta_o 0.70; Slr_m 0.2/0.1.
- **macro WRC (the new ingredient): P0_M 10/1 Pa; theta_o 0.30/0.33; Slr_M 0.0/0.001;
  Pd_M 1500/900 Pa; lambda_d 3.0/2.5.**
- intrinsic permeability k0 1.0e-20/0.6e-18 m2; b 12/6; krel_min 0.18/0.392.
- init: Po* 1.9/0.4 MPa; eM 0.402/0.526; e_m 0.222/0.134; PLm -120/-114 MPa.
- Multi-team ranges (D5.6 Tab 4-2..4-6): phi 15-27 deg, nu 0.25-0.41, kappa 1e-4..0.074,
  lambda 0.088..0.25; pellet k 4 orders of magnitude spread; block k 2 orders.

## 6. Reconciliation with the Task 3.3 phi_M->0 finding (Vinay 2026-06-01)
MGR is NOT the sharp-Pi wall. Two reasons: (a) it runs with a LIVE MACRO WRC (P0_M ~1-10 Pa)
so capillary water populates the macro through the wetting front — the "capillary-first
macro fill" the Task-3.3 macro-dry PRJs lacked; (b) homogenisation here is inter-LAYER
mass redistribution at moderate suction (120 MPa), not interlayer-dominated 200+ MPa.
=> Published UPC DSM captures it "very satisfactorily" (D5.7.txt:9201-9208). The native
   Pi-path DSM has not yet been run with a live macro WRC; that is the build's novelty
   and risk.

OPEN (Vinay, §9): calibration anchoring + how BExM macro/micro WRC maps to Pi-path. See
response 2026-06-01.
