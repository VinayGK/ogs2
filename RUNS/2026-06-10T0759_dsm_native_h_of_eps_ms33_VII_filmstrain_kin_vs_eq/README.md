<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-10T0759_dsm_native_h_of_eps_ms33_VII_filmstrain_kin_vs_eq
**Run finished:** 2026-06-10 07:59 (+0200)
**Branch / commit:** `dsm_native_h_of_eps` @ `7ff8861847`
**Binary:** `~/git/build/h_of_eps_20260609/bin/ogs`
**Full raw outputs (untracked, local):** `~/git/ogs/ms33_VII_filmstrain_kin_vs_eq_2026-06-09/out_{off,kinematic,equilibrium}/`
**Report / beamer:** beamer = `~/tex/implementation_and_theory/implementation/maxwell_from_psi.tex` Ch.2 (Steps 13-19, extended 2026-06-10 with this run's numbers). Written report: TODO.

---

# MS33 Model VII — film_strain_coupling: off vs kinematic vs equilibrium

**Date:** 2026-06-09.
**Binary:** `~/git/build/h_of_eps_20260609/bin/ogs`
(branch `dsm_native_h_of_eps` @ `7ff8861847`; design doc
`ProcessLib/RichardsMechanics/DSM/STRAINED_FILM_IMPLEMENTATION.md`).
**Base PRJ:** `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj`
(free swelling under 0.2 MPa for 200 d, then load ramps to 5 MPa and unload to
0.4 MPa, per the PRJ header). Variants differ ONLY by the appended tags
`<film_strain_coupling>off|kinematic|equilibrium</...>` +
`<film_strain_kappa>aggregate</...>`.

All three runs completed cleanly (507 / ~600 / ~530 steps; 1 rejected step
each; off 133 s, kinematic 200 s, equilibrium ~200 s at 6 OMP threads,
concurrent).

## Headline numbers (REV means)

| variant     | e @200 d (end swelling) | e @5 MPa peak | e @end (0.4 MPa) | n_l @end |
|-------------|------------------------|---------------|------------------|----------|
| off         | 1.4995 | 1.4548 | 1.4954 | 0.5993 |
| kinematic   | 1.3534 | 1.3214 | 1.3503 | 0.5743 |
| equilibrium | 1.4995 | 1.4548 | 1.4954 | 0.5993 |

Plots: `e_vs_time.png`, `e_vs_sigma.png`; raw table `filmstrain_compare.csv`;
extraction script `compare_kin_vs_eq.py`.

## Findings

1. **Ceiling-pinned chemistry in ALL THREE runs.** At the swollen state
   `max(n_l − φ) = 0.0` exactly — n_l sits on the total-porosity ceiling
   (`boundedMicroWaterContentCeiling`; this PRJ has no macro_porosity_floor,
   so the cap is φ itself). The micro fills the entire pore space; the
   saturated value of n_l is cap-determined, not μ-equilibrium-determined.
   This is the same over-swelling regime noted for VII (e≈1.50 here vs teams'
   ≈1.09).

2. **kinematic — strong, physically-sensible geometric feedback.** Swelling
   expansion (ε_v ≈ +0.4) thickens the films: w_eff = n_l(1+κε_v) > n_l ⇒
   Π(w_eff) < Π(n_l) ⇒ weaker swelling drive. Result: swelling plateau drops
   e 1.50 → 1.35, σ_sw,zz 15.7 → 13.1 MPa, micro pressure 47.5 → 38.9 MPa.
   The over-swelling moves TOWARD the inter-team value (1.09) — still above,
   but the negative feedback is the right direction. The kinematic difference
   propagates through the eigenstress half (weaker σ_sw → less expansion →
   smaller φ-cap → smaller n_l), with the chemistry cap following.

3. **equilibrium — inert on this path, by scale.** The loaded branch requires
   p_conf > Π(n_l); here Π ≈ 39–48 MPa while p_conf ≤ 5 MPa, so the branch
   point NEVER opens (the same micro-vs-REV gate-scale situation documented in
   MAXWELL_CONJUGATE_IMPLEMENTATION.md §7.2). What remains is the Derjaguin
   load term +b·p_conf/ρ_lR in μ_lR — visible in the output micro pressure
   (43.6 vs 47.9 MPa at peak) but unable to move n_l because n_l is
   ceiling-pinned (finding 1). Hence equilibrium ≈ off in e/φ/σ_sw to ~1e-5.

4. **Load/unload loops (e_vs_sigma.png).** All variants show a modest loop
   (Δe ≈ 0.05 at 5 MPa) with near-complete recovery on unload. CAVEAT: the
   apparent hysteresis includes exchange-kinetics lag (α_M rate limit) — e is
   still relaxing during the 5 d ramps (e.g. e keeps dropping on the first
   unload leg) — so the loop area is NOT a clean dissipation measure on this
   schedule. A slower path (or rate sweep) is needed to separate kinetic lag
   from true ∮ ≠ 0 (the operational form's expected Maxwell defect, design doc
   §9a).

## What this comparison says about the variants

- To make **equilibrium** mode bite on MS33-like stress levels, either the
  branch point fires only in much denser/drier states (Π smaller — NOT the
  case here: Π is tens of MPa), or the gate scale needs the REV-consistent
  reading (p_conf vs φ_m·Π ≈ 9 MPa at peak here — that WOULD have opened at
  5 MPa·... close; φ_m·Π ≈ 0.6·43 ≈ 26 MPa REV — still closed). The
  micro-scale force balance as implemented never engages on this benchmark.
  → Gate-scale remains the open §7.2 decision; this run is direct evidence it
  binds the equilibrium variant too.
- **kinematic** mode is live everywhere strain is nonzero and acts as a
  negative feedback on free swelling — a candidate lever for the VII
  over-swelling (e 1.50 → 1.35 with κ = aggregate, no recalibration). NOTE:
  K was calibrated with frozen geometry, so the kinematic plateau shift also
  signals the predicted K-recalibration coupling (design doc §8.6) — do not
  read e=1.35 as a calibrated result.

## Caveats

- Operational Derjaguin cut (design doc §9a): not Maxwell-exact; loop
  interpretation limited per finding 4.
- Single mesh/schedule; no dt-sensitivity check on the differences.
- n_l ceiling-pinning means the μ-side of the equilibrium variant is untested
  by this benchmark; a confined, lower-Π configuration would expose it.
