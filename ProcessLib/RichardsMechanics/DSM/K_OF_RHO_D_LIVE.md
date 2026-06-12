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
- [HISTORICAL — superseded 2026-06-12, see "Analytic tangent completion"
  below] The analytic dK/drho_d (= dK/dphi → strain chain) tangent
  contribution was OMITTED in the first cut: the residual saw the live K but
  the Jacobian did not carry its derivative.
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

## Analytic tangent completion (2026-06-12, Vinay approved)

The first-cut omission above is closed. Implemented as specified by Vinay
(2026-06-12): **Jacobian-only — the residual is untouched.**

- Chain: `rho_d = rho_SR*(1-phi)` => `drho_d/dphi = -rho_SR` =>
  `dK/dphi = -rho_SR * (dK/drho_d)` with `dK/drho_d` the EXACT per-segment
  slope of the piecewise-linear table — slope 0 at/outside the clamped
  edges, LEFT-segment slope at interior knots (mirrors `getValue`'s
  lower_bound interval selection). Exposed by `AugmentationPrefactorTable::
  getSegmentSlope` (`PotentialExchangeParameters.h`); NOT MathLib's
  `getDerivative`, which quadratically blends adjacent segment slopes and is
  not the derivative of the clamped value the residual actually uses.
- `effectiveAugmentationPrefactorPhiDerivative(params, phi)` returns
  `dK/dphi` [J/kg per unit phi]; 0 in every case where
  `effectiveAugmentationPrefactor` falls back to the parse-time scalar
  (mode off, no table, sentinel/NaN phi).
- mu-level K-partials: `mu_aug = sign*K*exp(-h/lambda)` is LINEAR in K, so
  `VanDerWaalsMicroPotentialData` carries exact `dmu_lR_dK` and
  `ddmu_lR_dnl_dK` (0 when the disjoining floor clamps).
- Wired into the live p-u augmentation Jacobian block in
  `RichardsMechanicsFEM-impl.h` (`assembleWithJacobian`): an additional
  `dmu_lR_dK_tot * dK/dphi * dphi/deps_v` contribution to d mu_lR/d eps_v,
  with `dphi/deps_v = (alpha - phi)/(1 + w)` derived analytically from
  `PorosityFromMassBalance` (the MS33 porosity law); other porosity laws are
  treated as strain-independent (chain 0, exact for Constant). The
  default-OFF `enable_dsm_swelling_up_jacobian` block carries a NOTE and is
  left without the chain (dead code).
- Off mode / frozen table / clamped edge: `dK/dphi == 0` -> the new block is
  skipped entirely -> bit-for-bit identical Jacobian (off-mode verified
  bitwise, see below).

### Verification (2026-06-12, measured)

- New unit tests `RichardsMechanicsLiveKOfRhoD.AnalyticPhiTangentMatchesFD
  InsideSegment` and `...ClampedEdgesAndKnots` (FD-vs-analytic, derived
  identity; scale-derived tolerances; structural in-test knots): PASS.
  Full RichardsMechanics suite: 41 tests, 39 PASS + 2 designed skips
  (build `Pi_fofnlev_20260611`, incremental).
- Off-mode regression: `ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj` run
  with the reference binary (`h_of_eps_20260609`, pre-tangent) vs this
  build: all 12 output VTUs bitwise-identical field data (points + every
  point/cell array).
- Live-K sanity (2-step truncated `1a_robin_A_Kl` copy,
  `task42_case1_2026-06-12/_diagnostics_1bKl/1a_robin_A_Kl_trunc2_diag.prj`):
  converges; measured Newton iterations 17/2/2 for steps 1-3, EQUAL to the
  pre-tangent full-run counts (no regression, no measured gain on this
  confined-path case).
- CURE TEST (the motivating 1b `*_Kl` step-1 singularity,
  `task42_case1_2026-06-12/_diagnostics_1bKl/README_DIAG.md`): NOT CURED.
  Both `1b_A_Kl` and `1b_B_Kl` still fail in time step #1 ("time stepper
  cannot reduce the time step size further"). Measured CHANGE in the
  Newton trajectory, though: with the tangent, |dx|_uz CONTRACTS for 7
  iterations (14.6 -> 11.4 -> 17.7 -> 14.0 -> 5.9 -> 2.4 -> 1.09) before
  blowing up at iteration 8 (168 -> 4.3e3 -> 5.1e4 -> ...), vs the
  pre-tangent monotonic divergence (14.6 -> 26 -> 44 -> 64 -> inf). The
  missing-tangent hypothesis is therefore PARTIALLY supported (basin
  improved) but the 1b compliant-top failure has an additional, still
  unidentified mechanism — evidence in `out_1b_{A,B}_Kl/run.log`
  (2026-06-12 17:04/17:05 runs); no further patching without Vinay's call.
