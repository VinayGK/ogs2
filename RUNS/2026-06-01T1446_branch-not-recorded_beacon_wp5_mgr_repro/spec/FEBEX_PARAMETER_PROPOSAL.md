# FEBEX bentonite — proposed native-DSM parameter set (FOR VINAY APPROVAL)

First-of-kind FEBEX parameterisation of the native Pi-path DSM (0 prior FEBEX PRJs).
Every literal cited to a deliverable line or a prior-approved memory value. Vinay
to approve / correct each group before any PRJ is built (§9, §12.1, §1.1).

Locators: D5.2.1.txt / D4.1.txt = the converted deliverables in
beacon_wp3_epfl_repro_2026-05-31/spec/ and beacon_wp5_mgr_repro_2026-06-01/spec/.
Layout-mode line numbers (pdftotext mangles some; key ones re-verified raw).

## DECISION 1 — swelling-pressure law coefficient (the calibration target) ***CONFLICT***
D5.2.1 PRINTS:  Ps = exp(6.0*rho_d - 9.07)  [MPa, g/cm3]  (D5.2.1.txt:4474 layout / 7449 raw)
PRIOR-APPROVED MEMORY (Task 13): ln Ps = 6.77*rho_d - 9.07  (Ps=11.46 @1.7, K=5.8e4)
Checks:
  b=6.0  -> 1.60:1.70 MPa  1.70:3.10  1.45:0.69   (matches NOTHING)
  b=6.77 -> 1.60:5.82 MPa  1.70:11.46 1.45:2.11   (matches prose "~5 MPa" D5.2.1:1465,
            "3-7 MPa" :4276; MGR ~2.0+/-0.5 @1.45 D4.1:3724; Task-13 11.46@1.7)
=> PROPOSED: b = 6.77 (printed "6.0" treated as OCR artifact of 6.77). 
   Source: Villar (2002) / ENRESA (2000), as in prior-approved Task-13 calibration.
   *** Vinay: approve b=6.77, or insist on the literal 6.0? ***
Resulting Ps calibration targets (if 6.77 approved):
  pellet rho_d=1.30 -> 0.76 MPa ; block rho_d=1.60 -> 5.82 MPa ; homog 1.45 -> 2.11 MPa.

## DECISION 2 — block density: 1.60 (D5.7 model) vs 1.70 (D4.1 experiment Table 5-1)?
D4.1 Table 5-3 (experiment): block initial rho_d 1.60-1.62 (MGR21-24). D5.7 model used 1.60.
D4.1 §5.1 prose says blocks "compacted at nominal 1.70" but Table 5-3 finals start ~1.60.
=> PROPOSED: block rho_d = 1.60 (matches Table 5-3 initial + D5.7 model). pellet = 1.30.

## Proposed §12.2 parameter groups

| group | value | source (proposed) | locator | status |
|---|---|---|---|---|
| grain density rho_s | 2700 kg/m3 (Gs 2.70) | ENRESA 2000 | D5.2.1:4263-4264 | verified |
| smectite content | "very high" (92% not in this deliv.) | ENRESA 2000 | D5.2.1:4179 | VERIFIED qual. |
| swelling law | Ps=exp(6.77 rho_d-9.07) MPa | Villar 2002 | D5.2.1:4474 + memory | DECISION 1 |
| K (pellet 1.30) | calibrate to Ps=0.76 MPa | Villar 2002 target | derived | to fit |
| K (block 1.60) | calibrate to Ps=5.82 MPa | Villar 2002 target | derived | to fit |
| phi0 pellet (rho_d 1300, rho_s 2700) | 0.5185 | = 1-1300/2700 | derived | computed |
| phi0 block (rho_d 1600, rho_s 2700) | 0.4074 | = 1-1600/2700 | derived | computed |
| WRC block (vG, 1.60-1.65) | P0=30 MPa, Sr0=0.32, Srmax... | ENRESA 2000 | D5.2.1:5156-5158 | verified |
| WRC modified-vG block | P0=35, lambda 0.30, lam_d 4000 | ENRESA 2000 | D5.2.1:5158 | verified |
| WRC pellet (~1.30, use 1.55-1.59 row) | P0=4.5 MPa (vG5) / 2.0 (mvG6) | ENRESA 2000 | D5.2.1:5160-5162 | verified (low-density row) |
| rel. perm exponent | k_r = Sr^3 (adopted) | ENRESA 2000 | (see RE-VERIFIED note; carry MS33 S_e^3) | OK |
| sat. hydraulic K | log K piecewise (K m/s): pellet 1.3e-12, block 4.9e-14 | ENRESA 2000 | D5.2.1:4714-4716 | VERIFIED (see note 1) |
| A_Hamaker | 2.2e-20 J | Israelachvili-Adams 1978 | literature anchor (fixed) | keep |
| micro EOS (a,b,rho_l0) | as MS33 (1e-16,1,100) | prior-commit | — | keep (FEBEX micro EOS unknown) |
| pellet grain dd | 1.95 g/cm3 | ENRESA 2000 | D5.2.1:3733 | verified (for micro/macro split) |

## Initial / boundary (from D4.1 Table 5-3, experiment-authoritative)
| qty | pellet | block | locator |
|---|---|---|---|
| rho_d initial | 1.30 | 1.60 | D4.1:3823,3825 (MGR23) |
| w initial | 3.5 % | 14.2 % | D4.1:3823,3825 |
| Sr initial | 9 % | 56 % | D4.1:3823,3825 |
| rho_d final (measured) | 1.34 | 1.51 | D4.1:3823,3825 |
| w final | 35.7 % | 31.1 % | D4.1:3823,3825 |
BC: MGR23 = constant water pressure 14 kPa (D4.1:3560) / 15 kPa (D5.7 model); 210 d.
Geometry: cell diameter 10.0 cm (D5.2.1... NB D4.1:3498 says inner diameter 10.0 cm,
  sample length 10 cm); pellet 5 cm bottom + block 5 cm top; hydrate bottom.
  *** NOTE: D4.1:3498 diameter 10.0 cm conflicts with D5.7 "5 cm" — re-verify. ***

## CHECK items — RE-VERIFIED 2026-06-01 (direct PDF read, agent locators corrected)

DECISIONS LOCKED (Vinay 2026-06-01): Ps law b=6.77; FEBEX specific surface 725 m2/g
(rest of micro EOS = MS33); re-verify-all-first (done below).

1. **sat. hydraulic conductivity** — VERIFIED D5.2.1:4714-4716 (eq 1.4, ENRESA 2000):
   log K = -6.00*rho_d - 4.09   (1.30<=rho_d<=1.47, r2=0.97, 8 pts)
   log K = -2.96*rho_d - 8.57   (1.47<=rho_d<=1.84, r2=0.70, 26 pts)   +/-30%.
   K in m/s (hydraulic conductivity). => pellet 1.30: log K=-11.89 -> K=1.3e-12 m/s;
   block 1.60: log K=-13.31 -> K=4.9e-14 m/s. CONVERT to intrinsic k=K*mu/(rho*g)
   = K*1.0e-7 approx (mu=1e-3, rho=1000, g=9.81): pellet k~1.3e-19, block k~5.0e-21 m2.
2. **retention vG table** — VERIFIED D5.2.1:5049-5068 (Table 1-8 wetting, ENRESA 2000):
   1.70-1.75 vG(1): P0=90 MPa, Sr0=0.45, Srmax=1.00.
   1.60-1.65 vG(3): P0=30 MPa, Sr0=0.32, Srmax=1.00 ; mvG(4): P0=35, lambda 0.30, Pd 4000, ld 1.5.
   1.58-1.59 vG(5): P0=4.5 MPa, Sr0=0.17 ; mvG(6): P0=2.0, lambda 0.10, Pd 1000, ld 1.3.
   BLOCK (1.60): use vG(3) P0=30 MPa, Sr0=0.32. PELLET (1.30): BELOW lowest band (1.58)
   -> EXTRAPOLATION; lowest-band vG(5) P0=4.5 MPa is the nearest cited anchor. FLAG: pellet
   macro WRC is an extrapolation of ENRESA Table 1-8 to 1.30; the lowest measured row is 1.58.
3. **smectite %** — only qualitative "very high content of montmorillonite" (D5.2.1:4179)
   in this deliverable; precise 92% NOT here (cosmetic, not a model input). Leave qualitative.
4. **MGR cell geometry** — VERIFIED D4.1:3407-3411 (§5.1.1, verbatim): "The large-scale
   oedometer consists of a cylindrical body ... The body has an inner diameter of 10.0 cm
   and the length of the sample inside was 10 cm." Block compacted first, pellets poured on
   top, "then the cell was overturned" (D4.1:3417) -> pellets at BOTTOM (hydration), block
   at top, for MGR22/23.
   => CELL = 10.0 cm diameter x 10 cm height. Axisym: r=50 mm, H=100 mm; pellet 0-50 mm
   (bottom, hydration face), block 50-100 mm (top). (Self-correction: an earlier draft of
   this line wrote an UNSOURCED "50mm x 50mm" mis-citing the Table 5-2 header line 3463 --
   that was wrong and is retracted. D5.7's "5 cm" was the modeller's reduced representation;
   the EXPERIMENT is 10 x 10 cm. The "35 mm" elsewhere is the unrelated EPFL Task-3.3 cell.)
   Cross-section area = pi*(0.05)^2 = 7.854e-3 m2 -> sets absolute water-intake volume target.

## Corrected phi0 (rho_s=2700 FEBEX)
- pellet rho_d=1300: phi0 = 1-1300/2700 = 0.518519 ; n_s = 0.481481
- block  rho_d=1600: phi0 = 1-1600/2700 = 0.407407 ; n_s = 0.592593

## Ps targets (b=6.77 LOCKED): pellet 1.30 -> 0.76 MPa ; block 1.60 -> 5.82 MPa.
K (vdw prefactor) CALIBRATED to these per material (fit, not asserted; §2-clean).

## CALIBRATION RESULT (Vinay approved 2026-06-01; calibrate_febex_K.py)
- BLOCK  rho_d=1.60: FITTED K = 25268 J/kg -> Ps = 5.784 MPa (target 5.82, -0.6%). S=1.
- PELLET rho_d=1.30: FITTED K =  3979 J/kg -> Ps = 0.752 MPa (target 0.76, -1.0%). S=1.
Sanity: FEBEX block K=25268 < MX-80/Dixon K=71900 at same 1.60 (FEBEX swells less, correct);
pellet/block K ratio 0.16 ~ Ps ratio 0.13 (K-governed, correct).
Constant-volume single element, BishopsSaturationCutoff=1, specific_surface=725 (FEBEX).
=> These two K values are the FEBEX calibration anchors for the MGR column build.

ALL CHECK ITEMS CLEARED + K CALIBRATED. Ready to assemble two-material MGR23 column.

## Open (Vinay):
- DECISION 1 (6.0 vs 6.77) and DECISION 2 (block 1.60 vs 1.70).
- FEBEX micro physics: no FEBEX-specific micro EOS / specific surface in deliverables;
  propose carrying the MS33 micro EOS unchanged (only K + macro WRC + rho_s are FEBEX).
  Is that acceptable, or does FEBEX need its own micro parameterisation?
