# KKT treatment of the micro-water ceiling: implementation note (branch dsm_mass_conservation_v3_kkt_ceiling_2026-09-30)

Status 2026-09-30: implemented as specified by Vinay, **NOT adopted, switchable, default = the shipped clamp (bitwise)**.
Formulation (derivation, weak forms, design, test plan) lives in the record folder
`~/ogs-models/scratch/2026-09-30_kkt_ceiling_impl/` (`DERIVATION.md` = D-n, `WEAK_FORMS.md` = W-n, `DESIGN.md`, `THEORY_FIXES.md`);
this file says only what the code does and where. Comments in the code cite those sections.
Branch parent: `dsm_mass_conservation_v2_strain_term_2026-09-30` (9cd3d00a4d; contains the IC fix 686fcd6ef8 and the B / B' / V2 switches, all default false).

## What it is (DERIVED, D-2.3 to D-2.5, D-4.4)

The ceiling `n_l <= n_max(eps_v) = phi` of the `scalar_micro_macro_mass_storage_mode` local solve is a complementarity condition with a
multiplier `lambda >= 0` (code name `multiplier`, never `lambda`: that symbol is also the vdW film spacing). `rho_lR = rho_lR_EOS(n)` slaves the second local
unknown, so the local problem is the scalar NCP `f(n) + dt*alpha_M*lambda/rho_lR = 0`, `f(n) = R_m(n, rho_lR_EOS(n), 0)` (the residual of the base solve).
Active candidate: `f(n_max) < 0` (no tolerance, no stored state read, ties and NaN inactive). At an active candidate a bracketed scan
(`scanForInteriorCeilingRoot`, `N_dec` nodes per decade from `n_floor = 1e-16` to `n_max`, plus a bisection on `f_n` in every cell with a + to - sign change) decides
whether the wall is the unique solution. Status (`MicroCeilingKktStatus`): 0 Interior (base solve, bitwise), 1 Active, 2 NonMonotone, 3 PremiseViolated
(2, 3: the KKT rule chooses nothing; base solve and base sink, counted, WARN once per process and status).
Active: `n = n_max`, `rhohat = S_s = (c_s*rho_l - rho_l_prev)/dt` (booked storage rate, received by micro AND macro), `r = rhohat_pot - S_s`,
`lambda = rho_lR*r/alpha_M`. The macro sink at an active IP is `-S_s` instead of `-rhohat_pot` (the water the clamp deletes, `r`, stays in the macro).
`s` (strain-term sign) is the existing `micro_mass_strain_term_eulerian` switch; `micro_exchange_source` keeps its meaning `rhohat_pot`
(D-0.2), the new `micro_exchange_received` is what the macro sink uses.

## PRJ tags (inside `<potential_exchange>`, per-medium override inherits the global block)

| tag | values | default | notes |
|---|---|---|---|
| `micro_ceiling_treatment` | `clamp`, `kkt` | `clamp` | `kkt` requires `scalar_micro_macro_mass_storage_mode`; exclusive with `ceiling_micro_storage_exchange` / `ceiling_micro_storage_includes_strain`; requires `film_strain_coupling = off`; `beta_SR != 0` is `OGS_FATAL` at assembly |
| `micro_ceiling_pu_tangent` | `overwritten`, `kkt_active`, `all_exchange` | `overwritten` | Q9. `overwritten`: the line `local_Jac.pu = Kpu/dt` erases every exchange p-u entry (shipped). `kkt_active`: the KKT-active entries are added after it. `all_exchange`: `Kpu/dt` is added to the accumulated exchange entries (Maxwell, film, live-K entries at inactive points revive) and the KKT-active entries are added too. Only with `kkt` |
| `micro_ceiling_fd_check` | `true`, `false` | `false` | route-B debug flag: central-difference check of the assembled element Jacobian, log lines `KKT-FD`; LinearElasticIsotropic solids only. Only with `kkt` |
| `micro_ceiling_scan_nodes_per_decade` | integer >= 2 | 8 | `N_dec`, a PROPOSAL (repo rule 1.2) that needs Vinay's approval. Only with `kkt` |
| `micro_ceiling_trace_elements` | whitespace separated element ids | empty (off) | iteration trace `kkt_trace.csv` in the working directory. Only with `kkt` |

A label line `MASSFIX V3 label: ...` is printed when `kkt` is on.

## Code map

| file | change |
|---|---|
| `PotentialExchangeParameters.h` | enums `MicroCeilingTreatment`, `MicroCeilingPuTangent`, five members at the end, `toString`, `isKktCeiling` |
| `CreateRichardsMechanicsProcess.cpp` | parse, validation, label, aggregate init (appended) |
| `ConstitutiveRelations/MicroCeiling.h` (new) | twelve `StrongType<double>` fields (`micro_ceiling_status`, `_multiplier`, `micro_exchange_received`, `micro_ceiling_rejected_exchange`, `_flips`, `_nonmonotone`, `_premise`, `_eps_seen`, `_inc_last`, `_attempt_t`, `_inc_alt`, `_inc_same`) and `MicroCeilingTraceTag` |
| `ConstitutiveRelations/ConstitutiveData.h` | the twelve types appended to `StatefulData` (reflection makes them IP writers / secondary variables; written only with `kkt`) |
| `RichardsMechanicsFEM-impl.h` | `ReducedMicroLiquidDensityData::drho_lR_drho_LR`; after `solveReferenceMassStorageCoupledState`: `scanForInteriorCeilingRoot`, `solveReferenceMassStorageKktState`, `computeCeilingKktActiveExchangeTangents`, `porosityDerivativeWrtVolumetricStrain`, `updateCeilingIterationDiagnostics`, `writeMicroCeilingTraceLine`; dispatch in `updateMicroscaleHydraulicState` (new arguments `t`, `trace_tag`); assembler blocks in `assemble` (residual only) and `assembleWithJacobian`; Q9 at the `Kpu/dt` line; route-B FD block at the end of `assembleWithJacobian` |
| `RichardsMechanicsFEM.h` | `assembleWithJacobianEvalConstitutiveSetting` gets a trailing `trace_tag` argument |

No existing function or lambda is refactored or shared (the file-scope `#pragma STDC FP_CONTRACT OFF` in the header records that a refactor boundary alone once changed clang's fusions);
the KKT solver carries its own evaluation in the same expression order as `mass_residual`. With `micro_ceiling_treatment = clamp` (or the tag absent) every changed line
is reached with the old values: `kkt_active = false` selects the old expression, `dn_l_dpL` / `dn_l_dK_sw` are evaluated as before, and the `Kpu/dt` line is the old assignment.

## Active tangents (DERIVED, D-5.2)

`d rhohat/d p_L = c_s*phi*(d rho_lR/d rho_LR)*(d rho_LR/d p_L)/dt`, `d rhohat/d eps_v = [-s*rho_l + c_s*phi'*(rho_lR + phi*d rho_lR/dn)]/dt`,
`c_s = 1 - s*Deps`, `rho_l = phi*rho_lR`, `phi' = (alpha - phi)/(1 + Deps)` (0 where the porosity clamp acts). At active IPs `dn_l/dp_L` and `dn_l/dK` are set to 0 (n_l = n_max(eps_v)).

## Deviations from DESIGN.md (all also in the IMPLEMENT log of the record folder)

1. Trace file name: `kkt_trace.csv` in the process working directory instead of `<output_prefix>_kkt_trace.csv` (the output prefix is not reachable from the local assembler).
2. `micro_ceiling_fd_check` accepts `LinearElasticIsotropic` only (DESIGN: refuse MFront); history-dependent solids carry material state that the check does not snapshot.
3. Unit-test batch 1 (DESIGN.md 4.1, file `Tests/ProcessLib/RichardsMechanics/MicroCeilingKkt.cpp`) is NOT part of this change (separate step; the free functions above are written to be testable).
