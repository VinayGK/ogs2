# Live K(rho_d) — augmentation prefactor at the EVOLVING dry density

Status: first-cut trial implementation (2026-06-11), per Vinay's order
"K(rho_d) try it" (2026-06-10). Branch `dsm_native_h_of_eps`.

## What it is

The existing parse-time K(rho_d) table feature (memory:
`project_dsm_k_of_dry_density`, branch `dsm_native_pdisj_maxwell_kofdd`;
ported to this branch) resolves the augmentation prefactor
K = `potential_augmentation_prefactor` ONCE at parse time from a
piecewise-linear table `potential_augmentation_prefactor_vs_dry_density`
evaluated at the PRJ-supplied initial/target `<dry_density>`. K is then a
per-material constant in time.

This LIVE variant keeps the table un-frozen and re-evaluates

    K = K_table(rho_d),    rho_d = rho_SR * (1 - phi)   [kg/m^3]

at run time, where `rho_SR` = `micro_solid_density_reference` and `phi` is
the CURRENT total porosity at the evaluation site. As the column swells
(phi up, rho_d down) the prefactor tracks the calibrated K(rho_d) curve.

## PRJ interface

```xml
<potential_exchange>
  ...
  <potential_augmentation_prefactor_vs_dry_density>
    <dry_densities>...</dry_densities>
    <prefactors>...</prefactors>
  </potential_augmentation_prefactor_vs_dry_density>
  <dry_density>1400</dry_density>  <!-- optional in live mode: fallback K -->
  <potential_augmentation_prefactor_live_dry_density>true</potential_augmentation_prefactor_live_dry_density>
</potential_exchange>
```

- `potential_augmentation_prefactor_live_dry_density` default `false` →
  parse-time freeze, bit-for-bit the existing behavior (verified: dd1400
  off-mode regression below).
- `true` REQUIRES the table (OGS_FATAL otherwise) and is mutually exclusive
  with the scalar `potential_augmentation_prefactor` (unchanged rule).
- In live mode `<dry_density>` is optional; when given, K(dry_density) is
  stored as the FALLBACK scalar used at evaluation sites that have no
  porosity in scope (see below).

## Implementation

One inline helper, `effectiveAugmentationPrefactor(params, phi)` in
`PotentialExchangeParameters.h`:

- live mode && table && `std::isfinite(phi)` →
  `table->getValue(rho_SR * (1 - phi))`;
- otherwise → the parse-time scalar `potential_augmentation_prefactor`.

`MathLib::PiecewiseLinearInterpolation::getValue` HOLDS the endpoint values
outside `[rho_d_min, rho_d_max]` (verified in
`MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp`), so K is
clamped at the table range ends — no extrapolation.

Threaded through every site in `RichardsMechanicsFEM-impl.h` that feeds
`potential_augmentation_prefactor` into `computeVanDerWaalsMicroPotential` /
`computeStrainedFilmState`:

- micro-potential fold point `applyFilmPressureMicroPotential` and
  `computeActiveMicroPotential`, plus the predictor/coupled
  ScalarReferenceMassStorage local solves — phi from
  `PotentialExchangeLocalSolveContext::phi` (infinity sentinel for
  default-constructed contexts → scalar fallback);
- swelling-eigenstress increment
  `computeReferenceMicroPorositySwellingStressIncrement` /
  `computeSwellingStressIncrement` — new trailing defaulted argument
  `total_porosity` (NaN sentinel → scalar), passed from
  `updateSwellingState` as the stateful total `PorosityData.phi`; one K for
  both the prev and curr Pi evaluations of an increment (mirrors the
  held-fixed p_conf telescoping convention);
- the assembleWithJacobian Maxwell p-u tangent and swelling-Jacobian blocks
  — the local total-porosity `phi` at the integration point.

Sites without any porosity in scope (the GP eigenstress-difference driver,
unit tests with default contexts) fall back to the scalar — in OFF mode this
makes every site bit-for-bit identical to before.

## Provisional choices (first cut — flagged)

- LINEAR interpolation between the calibrated knots (the
  PiecewiseLinearInterpolation shape); the interpolation shape between
  calibration anchors is UNDECIDED (memory: `project_dsm_k_of_dry_density`)
  — linear is the provisional placeholder, Vinay's call pending.
- The analytic dK/drho_d (= dK/dphi → strain chain) tangent contribution is
  OMITTED in this first cut: the residual sees the live K but the Jacobian
  does not carry its derivative. Predicted, not benchmarked: Newton may pay
  extra iterations — accepted for this trial (Vinay's "try it").
- NO double counting: only the K value fed to the UNCHANGED augmentation
  law `mu_aug = sign*K*exp(-h/lambda)` varies with rho_d; no new physical
  term, no second rho_d-dependence is introduced anywhere.

## Verification (2026-06-11)

- Unit tests `RichardsMechanicsLiveKOfRhoD.*` in
  `Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`
  (structural in-test knots, not physical values): off-mode returns the
  scalar; live mode reproduces the table at knots and interior points
  (linear-interpolation identity, derived in-file); endpoint clamping.
  31/31 RichardsMechanics unit tests pass (28 prior + 3 new).
- Off-mode regression: `ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj` (no
  live flag) with this build → final mean sigma_zz = -4.9218 MPa, matching
  the recorded baseline (eurad-anchors run
  `2026-06-10_0841_dsm_native_h_of_eps_successful`, README_ms33.md §2).
- A live-mode production run (e.g. free-swelling Model VII with a Dixon
  K(rho_d) table) has NOT been performed yet — behavior under live K is
  PREDICTED only, pending Vinay's run plan.
