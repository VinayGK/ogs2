# MGR material identity — CORRECTION (2026-06-01)

## The MGR tests are FEBEX bentonite, NOT MX-80.

Source: BEACON WP4 D4.1 (Sun, Tanttu, Villar, Wieczorek), §5.1 CIEMAT pellet/block
isochoric tests (the MGR series). File spec/D4.1.txt.

Evidence:
- D4.1 experiment table line 303: MGR row material = "FEBEX (MX-80?)" — the
  parenthetical "?" is the deliverable authors' own annotation.
- D4.1:3513-3534 (§5.1 body, decisive): "The blocks were made of FEBEX bentonite
  compacted at a nominal dry density of 1.7 g/cm3 ... pellets were made of the
  same bentonite ... The bentonite used was FEBEX bentonite, with a smectite
  content of 92 +/- 3 %, ... a specific gravity of 2.70 ... montmorillonite ...".
- D4.1:3941: MGR subsamples compared against "FEBEX samples of the same [density]".

=> CONSEQUENCE: the earlier MGR build staged a Dixon/MX-80 calibration anchor
   (K=71900 @ rho_d=1600, MX-80). That is the WRONG MATERIAL FAMILY. The 14-vs-3.5
   MPa "discrepancy" flagged on the block smoke run is explained: a Dixon-MX-80
   anchor was applied to a FEBEX test. Guardrail §12.1: FEBEX is its own allowed
   source family (Villar 2002; Lloret & Villar 2007). K MUST be calibrated to a
   FEBEX swelling-pressure target.

## Authoritative FEBEX MGR data (D4.1, in text — not figure-only)

Materials / geometry (Table 5-1, D4.1:3486-3534):
- Block: FEBEX, dry density 1.70 g/cm3, height 50 mm, w ~14 %.
- Pellets: FEBEX, dry density 1.30 g/cm3 (poured; pellet-grain dd 1.95 implied by
  packing 0.94), height 50 mm, w ~6 % (MGR21/22) or ~3.5 % (MGR23/24).
- Cell: diameter 50 mm, total height 100 mm.
- System average dry density ~1.45 g/cm3.
- FEBEX: smectite 92+/-3 %, specific gravity 2.70 (=> rho_s = 2700 kg/m3),
  hygroscopic w ~13.5 %.
  NOTE this differs from D5.7's UPC modelling values (block 1.60, rho_s 2735).
  D4.1 is the EXPERIMENT; D5.7 is one team's model setup. Prefer D4.1 for IC/data,
  flag the model-vs-experiment difference.

Swelling-pressure target (D4.1:3540-3542, citing Hahn 2007):
- "With the FEBEX bentonite at its hygroscopic water content, the swelling
  pressure measured at a dry density of 1.45 Mg/m3 was around 2 MPa."
- D4.1:3724: theoretical Ps at MGR final density 1.45 = "2.0 +/- 0.5 MPa".
- D4.1:3612: MGR23 stabilised axial pressure ~3 MPa.
=> FEBEX swelling-pressure calibration target: ~2 MPa at rho_d=1.45 (homogenised),
   block 1.70 and pellet 1.30 endpoints from the FEBEX Ps(rho_d) law (Villar).

Test BCs / durations (Table 5-2, D4.1:3555-3563):
- MGR21: constant pressure 14 kPa, 34 d.
- MGR22: constant flow 0.05 cm3/h, 266 d.
- MGR23: constant pressure 14 kPa, 210 d.   <- primary calibration target
- MGR24: constant pressure 14 kPa, 14 d.
  (D5.7 modelling used 15 kPa / 0.047 cm3/h — model-rounded; prefer D4.1 14 kPa /
   0.05 cm3/h for the experiment.)

## Open decisions (Vinay, §9 / §12.1) before any FEBEX MGR build
1. K calibration: re-anchor to FEBEX Ps (~2 MPa @1.45; or block 1.70 / pellet 1.30
   density-compliant FEBEX values). Source = Villar 2002 / Lloret & Villar 2007 /
   Hahn 2007 as cited in D4.1. NOT Dixon/MX-80.
2. rho_s = 2700 (FEBEX, Gs=2.70) replaces 2780 (MX-80). Resets phi0, n_s.
3. FEBEX retention curve + micro physics: do we have a FEBEX-parameterised native
   DSM anywhere, or is this the first FEBEX build of the native model?
4. Block density 1.70 (D4.1 experiment) vs 1.60 (D5.7 model) — which to target.

The earlier mgr_block_dd1600_le.prj (MX-80 anchor) is SUPERSEDED by this finding
and kept only as a record; not to be used for the FEBEX MGR build.
