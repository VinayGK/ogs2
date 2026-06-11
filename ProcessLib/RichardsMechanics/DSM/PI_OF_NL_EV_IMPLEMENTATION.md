# Equipresent disjoining law Pi(n_l, eps_v) + compressible-liquid carrier (D2)
# — IMPLEMENTATION & TEST DESIGN

**Branch:** `dsm_native_Pi_fofnlev` (originally off `dsm_native_h_of_eps` @
`7ff8861847`; REBASED 2026-06-11 onto `23a723cc3c` "live K(rho_d) table mode",
after which `dsm_native_h_of_eps` was retired — fully contained in this
branch, ref deleted per Vinay; its worktree directory remains on disk,
detached, because run snapshots cite its untracked campaign artifacts as
canonical sources).
**Worktree:** `~/git/ogs-worktrees/dsm_native_Pi_fofnlev_wt`.
**Builds on:** `STRAINED_FILM_IMPLEMENTATION.md` (esp. §9a exact-energy route and
decision D2), `MAXWELL_PAIR_RESTORATION.md`, `MAXWELL_CONJUGATE_IMPLEMENTATION.md`,
`MAXWELL_CONJUGATE_REV_REFERENCING.md`.
**Goal owner:** Vinay (2026-06-11): *"Pi does not carry p_conf. All the circus we
did was to consider p_conf in the coupling energy correctly, but the p_conf is
still a bolt-on to the micro, not an energetically compliant bolt-on. Equipresence
says Pi(n_l, <mechanical state>)."* Resolution of the same discussion: the
equipresent argument is the CONFIGURATION variable eps_v (not the force p_conf —
a Helmholtz energy may not take the conjugate force of its own variable);
Pi(n_l, p_conf) is legal only after a Legendre transform to G(p_conf, n_l).

Tags: **[D]** derived/established · **[DECISION]** Vinay's call before coding ·
**[PRED]** predicted, not verified · **TODO(Vinay)** value/approval needed.

**Audience note.** This document is written so that an implementing agent needs
NO further physics derivation: every formula is given with its derivatives and
units; every code touch point is given with file + function + line anchor (line
numbers as of `7ff8861847`); every test is specified in CLAUDE.md §3 format.
The implementing agent makes NO numerical or formulation decisions: everything
marked [DECISION]/TODO(Vinay) blocks on the owner. CLAUDE.md §1.1/§3/§5/§7
apply in full.

---

## 0. Scope and goal

Two deliverables, PRJ-selectable, default-off (bit-for-bit when off):

- **(E) Exact one-Psi film pair** — close STRAINED_FILM_IMPLEMENTATION.md §9a:
  replace the *operational* Derjaguin cut (law evaluated at w_eff plus a
  hand-added `+b*p_conf/rho_lR`) by the EXACT energy route: a single
  Psi_film(n_l, eps_v) whose closed-form strain integral exists for the
  kinematic h-law, with sigma_sw = dPsi/deps_v and
  mu_mech = (1/(nS*rho_lR)) dPsi/dn_l. Maxwell holds identically; the
  reversibility loop (∮dW = 0) closes — the property the operational form is
  expected to fail [PRED, design doc §9a].

- **(L) Compressible-liquid carrier (D2 proper)** — give the squeezed micro
  liquid an energy home: at fixed water content a thinner film densifies the
  liquid; with a finite confined-liquid bulk modulus K_liq this stores
  Psi_liq and produces the over-pressure part of the chemical potential
  (mu = psi + p/rho — the Derjaguin v*P term) WITHOUT any bolt-on. The
  trigger condition of D2 ("revisit if variant A shows promise") is satisfied:
  kinematic mode measurably moved MS33 VII toward the inter-team value
  (e 1.4995 -> 1.3534, run `ms33_VII_filmstrain_kin_vs_eq_2026-06-09`).

Out of scope: dissipative residual / hysteresis flow rule (B2-class, unchanged);
MFront parity (native-only branch, like `dsm_native_h_of_eps`); the equilibrium
variant's inversion machinery (kept as is; see §2.5); any re-calibration of K
(separate step, never combined with verification — CLAUDE.md §2).

## 1. The code as read [D] (anchors at `7ff8861847`)

- `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`
  - `computeVanDerWaalsMicroPotential(n_l, rho_lR, nS, rho_SR, A, Sa, sign,
    K_aug, lambda, dnS_dnl, n_l_floor)` — bare law, no strain argument. Film
    thickness h = n_l/(nS*rho_SR*Sa). Returns mu_lR and d/dn_l, d2/dn_l2,
    d/drho_lR, d/dnS, d/drho_SR.
  - `computeIntegrableMechanicalMicroPotential(...)` (~L540–577) — the shipped
    frozen-h partner mu_mech = -[(Pi + n_l*Pi')*eps_v + 0.5*b*K_d*eps_v^2]/rho_lR
    with analytic d/deps_v, d/dn_l, d/drho_lR. Comment block carries the
    Maxwell identity derivation.
  - `invertDisjoiningPressure(...)` (L609) — Newton+bisection solve of
    Pi(w) = p_target (equilibrium variant B).
  - `computeStrainedFilmState(mode, kappa_mode, n_l, active_nS, eps_v, p_conf,
    ...)` (L683) — returns {w_eff, dw_eff_dnl, dw_eff_deps_v, loaded_branch}.
    Kinematic: w_eff = n_l*(1 + kappa*eps_v), kappa = active_nS | 1, positivity
    guard f >= 1e-6, kappa frozen at GP (B1).
  - `FilmStrainCouplingMode { Off, Kinematic, Equilibrium }`,
    `FilmStrainKappaMode { Aggregate, Unity }` enums + the
    `PotentialExchangeParameters` struct (search `film_strain_coupling` for the
    member block).
- `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`
  - `applyFilmPressureMicroPotential(...)` (L680) — THE equipresent fold point.
    Strained modes: bare law at w_eff (nS chains frozen, dnS_dnl = 0), chain
    rule through dw_eff_dnl, macro-floor cutoff, then the OPERATIONAL load term
    `mu_load = b*p_conf/rho_lR` added UNCUT; replaces (does not add to) the
    shipped integrable partner. Off mode: integrable partner added after
    cutoff (frozen-h path).
  - Micro liquid EOS (L503–513): `rho_lR = rho_l0*exp(-a_rho*omega_l^b_rho)
    + rho_LR`, `omega_l = n_l*rho_lR/(nS*rho_SR)` — packing-driven
    (adsorption structuring), NO pressure and NO strain argument. This is the
    pressure-blindness D2 names.
  - `p_L_m = -rho_lR*mu_lR` (L1273) — micro pressure from the FULL potential,
    after the fold; propagates to every consumer.
  - Eigenstress site (~L1689): `sigma_sw` increment; strained modes evaluate
    Pi at the SAME w_eff (eps_v threaded with NaN-sentinel defaults).
- Tests: `Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`
  (7 tests; includes `ReplacementIsExclusiveAtZeroStrain`).
- PRJ parsing: `CreateRichardsMechanicsProcess.cpp` (search
  `film_strain_coupling`); tag docs under
  `Documentation/ProjectFile/` (grep for `film_strain_coupling` to find the
  exact subtree for `<potential_exchange>` children).

Sign conventions (as in the code; do NOT re-derive): Pi > 0 repulsive
(disjoining); mu_vdW < 0; sigma tension-positive; eps_v > 0 expansion;
p_conf > 0 compression. The augmentation sign is whatever
`potential_sign_factor` produces — all formulas below are written per term
against the h-dependence only and are sign-agnostic (CLAUDE.md §1.1: signs are
protected literals; reuse the code's, never re-pick).

## 2. Physics design [D]

### 2.0 The objection this branch answers

The operational pair (current `applyFilmPressureMicroPotential`, strained
modes):

```
mu_lR    = mu_bare(w_eff) + b*p_conf/rho_lR        (load term BOLTED ON)
sigma_sw = -phi_m_REV*(Pi(w_eff) - b*p_conf)
```

Pi itself never carries the mechanical state beyond the geometric w_eff; the
over-pressure enters as an additive term that is NOT the derivative of any
energy at frozen state — Maxwell defect O(Pi*eps_v) (design doc §9a). The
compliant structure is ONE energy:

```
Psi_micro(n_l, eps_v) = Psi_film(n_l, eps_v)  [disjoining, h live      — (E)]
                      + Psi_liq (n_l, eps_v)  [liquid compression      — (L)]
sigma_sw,m = dPsi_micro/deps_v
mu_lR,mech = (1/((1-phi_M)*rho_lR)) * dPsi_micro/dn_l
```

Maxwell then holds identically; no separate load term anywhere.

### 2.1 (E) Exact film energy — closed forms (kinematic h-law)

Notation: nS = active_nS = (1-phi_M) frozen at the GP (B1); kappa per
`FilmStrainKappaMode`; h0 = n_l/(nS*rho_SR*Sa); xi = h0/lambda_aug;
x = kappa*eps_v; w(e) = n_l*(1 + kappa*e). Split the bare disjoining pressure
at the UNSTRAINED state into its two terms (both returned by the existing law;
Pi_T = -rho_lR * mu_T for each term T):

```
Pi_vdw(n_l)  ~ n_l^-3        (cubic core)
Pi_aug(n_l)  ~ exp(-xi)      (exponential augmentation; sign per code)
```

Strain integral I_T(n_l, eps_v) = ∫_0^{eps_v} Pi_T(w(e)) de, closed form
per term [design doc §9a; re-derived and CONFIRMED here]:

```
I_vdw = Pi_vdw(n_l) * J3(x)/1,   J3(x) = [1 - (1+x)^-2] / (2*kappa)     [Pa·(-)]
I_aug = Pi_aug(n_l) * Jx(x),     Jx(x) = (1 - exp(-xi*x)) / (xi*kappa)  [Pa·(-)]
```

(Both J -> eps_v as kappa -> 0; see §4.6 stable forms.) The film energy and
its two halves, with the transmitted-load work S (see [DECISION Q2]):

```
Psi_film(n_l, eps_v) = -(1-phi_M)*n_l * [ I_vdw + I_aug + S(eps_v) ]   [J/m^3 REV]

S(eps_v) = 0.5*b*K_drained*eps_v^2      (route R3, default; dropped when (L) is ON)
```

**Eigenstress half** (differentiate the integrals back down — exact):

```
sigma_sw,m = dPsi_film/deps_v
           = -(1-phi_M)*n_l * [ Pi_vdw(n_l)*(1+x)^-3 + Pi_aug(n_l)*exp(-xi*x)
                                + b*K_drained*eps_v ]                       [Pa]
           = -(1-phi_M)*n_l * [ Pi(w_eff) + b*K_drained*eps_v ]
```

i.e. EXACTLY the operational eigenstress on the drained line
(p_conf = -K_drained*eps_v) — the eigenstress half is unchanged; only the
mu half changes. [D]

**Micro-potential half.** dPsi_film/dn_l needs d/dn_l of n_l*I_T at fixed
eps_v. Use Pi_vdw ~ n_l^-3 => d(n_l*I_vdw)/dn_l = -2*I_vdw. For the
augmentation, xi depends on n_l (dxi/dn_l = xi/n_l) and
dPi_aug/dn_l = -(xi/n_l)*Pi_aug (B1: rho_lR, nS frozen). Product rule on
n_l * Pi_aug(n_l) * Jx(x; xi(n_l)):

```
d(n_l*I_aug)/dn_l = I_aug * (1 - xi) + n_l*Pi_aug * dJx/dxi * (xi/n_l)
dJx/dxi = [ x*exp(-xi*x) - Jx*kappa ] / (xi*kappa) * kappa
        = ( x*exp(-xi*x) - (1 - exp(-xi*x))/xi ) / (xi*kappa) * ... 
```

— the implementing agent MUST NOT hand-simplify further: implement
dJx/dxi = [x*exp(-xi*x)*xi - (1 - exp(-xi*x))/kappa ] ... STOP — to remove
transcription risk the implementation SHALL compute the n_l-derivative of the
aug term as the derivative of the assembled expression

```
F_aug(n_l, eps_v) := n_l * Pi_aug(n_l) * (1 - exp(-xi(n_l)*kappa*eps_v)) / (xi(n_l)*kappa)
```

by symbolic differentiation ONCE (e.g. sympy, see §5 T-0) and transcribe the
generated C++ with a unit comment, then VERIFY against central finite
differences in the unit test (FD tolerance derived from scale, §5 T-3). The
same procedure applies to the vdW term (whose result must reproduce
-2*I_vdw/(...) analytically — a built-in cross-check). Then

```
mu_mech = (1/((1-phi_M)*rho_lR)) * dPsi_film/dn_l                       [J/kg]
        = -( d(n_l*[I_vdw+I_aug+S])/dn_l ) / rho_lR
```

(the (1-phi_M) cancels exactly as in the shipped partner). Derivative blocks
required for the Jacobian (all analytic, unit-commented per CLAUDE.md §4.2):

```
dmu_mech/deps_v   = -( Pi(w_eff) + n_l*dPi/dw|_{w_eff}*(1+x)... ) — equals
                    (1/((1-phi_M)*rho_lR)) * dsigma_sw/dn_l by Maxwell;
                    IMPLEMENT INDEPENDENTLY, assert equality in test T-3.
dsigma_sw/deps_v  = -(1-phi_M)*n_l * [ dPi/dw|_{w_eff} * n_l*kappa + b*K_drained ]  [Pa/-]
dsigma_sw/dn_l    = exact d/dn_l of the eigenstress half (chain through both
                    n_l prefactor and w_eff = n_l*(1+x))                     [Pa/-]
dmu_mech/dn_l     = from the symbolic derivative of dPsi/dn_l               [J/kg per -]
dmu_mech/drho_lR  = -mu_mech/rho_lR  (density mirroring, as shipped)        [(J/kg)/(kg/m^3)]
```

**Limits [D], all become unit tests (§5):**
- eps_v = 0: I_T = 0, sigma_sw = -(1-phi_M)*n_l*Pi(n_l), mu_mech = 0 + S-part 0.
  Exact reduction to the unstrained state.
- kappa -> 0: J3, Jx -> eps_v and the pair reduces EXACTLY to the shipped
  `computeIntegrableMechanicalMicroPotential` (the frozen-h partner is the
  kappa->0 limit of the exact pair — unlike the OPERATIONAL form, which does
  NOT reduce to it; this resolves the "slopes differ" note under D3 of the
  strained-film doc in favour of the exact route).
- Maxwell: dsigma_sw/dn_l == (1-phi_M)*rho_lR * dmu_mech/deps_v identically.
- Loop: ∮ dW = ∮ (sigma_sw deps_v + (1-phi_M)*rho_lR*mu_mech dn_l) = 0 on any
  closed (eps_v, n_l) path (it is a gradient field) — to FD/quadrature
  tolerance.

### 2.2 (L) Compressible-liquid carrier — the over-pressure made compliant

Kinematics [D]: the strained film at FIXED micro water mass has volume
V_film ∝ n_l*(1 + kappa*eps_v) while the mass is ∝ n_l*rho_ref, so the
geometric density is

```
rho_geom(n_l, eps_v) = rho_EOS(omega_l) / (1 + kappa*eps_v)
eps_liq := ln(rho_geom/rho_EOS) = -ln(1 + kappa*eps_v)   (≈ -kappa*eps_v)
```

Energy (per unit REV volume), with phi_m_REV = (1-phi_M)*n_l and a confined-
liquid bulk modulus K_liq [TODO(Vinay), §8 Q4]:

```
Psi_liq(n_l, eps_v) = 0.5 * K_liq * phi_m_REV * eps_liq^2                [J/m^3]
p_liq  = K_liq * eps_liq                       (liquid over-pressure)    [Pa]
```

Halves (differentiate; deps_liq/deps_v = -kappa/(1+kappa*eps_v)):

```
sigma_liq = dPsi_liq/deps_v
          = -phi_m_REV * p_liq * kappa/(1+kappa*eps_v)                   [Pa]
mu_liq    = (1/((1-phi_M)*rho_lR)) * dPsi_liq/dn_l
          = 0.5 * K_liq * eps_liq^2 / rho_lR                             [J/kg]
            + (p_liq/rho_lR) * (n_l d eps_liq/dn_l term = 0 at fixed kappa)
```

NOTE the structure: mu_liq = psi_liq + p_liq/rho-type split emerges when the
mass-derivative is taken at fixed film volume; with the kinematic pinning
above, d eps_liq/dn_l = 0 (B1: kappa frozen, volume tracks n_l), so the p/rho
term enters mu through the EXCHANGE-relevant route only via the energy share.
**This is exactly the place where the formulation choice lives** — whether the
v*P Derjaguin term is recovered through (a) the fixed-volume mass derivative
(adds p_liq/rho_lR to mu_liq; requires V_film independent of n_l in the
derivative — a different B1 freeze), or (b) the kinematic pinning above (no
p/rho term in mu; expulsion driven by the energy share only). The two differ
at first order. **[DECISION Q3 — Vinay; the doc implements (a) as default
ONLY IF approved, else (b); do not code before the call.]** Magnitude flag
[PRED]: with K_liq at the bulk-water scale and |kappa*eps_v| ~ 0.05,
p_liq ~ O(100 MPa) > Pi ~ O(40 MPa) — the liquid channel would DOMINATE and
relax only through the exchange (alpha_M/Damköhler-limited). Conditioning
note: Jacobian gains K_liq*kappa^2 stiffness; expect smaller time steps until
the exchange relaxes the pressure [PRED — must be observed in T-8 before any
claim].

Exclusivity [D, mirrors D3 mechanism-ownership]: when (L) is ON, BOTH the
operational `mu_load = b*p_conf/rho_lR` AND the S = 0.5*b*K_d*eps_v^2 term in
(E) are REMOVED (the load now reaches mu through the liquid energy; running
both double-counts the over-pressure). When (L) is OFF, (E) keeps S (route
R3) so the exact pair reproduces the shipped behaviour class.

### 2.3 What stays untouched

- The bare law `computeVanDerWaalsMicroPotential` (single-argument Pi(n_l)):
  unchanged; (E) calls it per term at w_eff and at n_l.
- The packing EOS rho_EOS(omega_l) (L503–513): UNCHANGED — (L) layers the
  geometric densification ON TOP (rho_geom), it does not modify the
  adsorption-structuring law. Double-count guard: rho_geom is used ONLY inside
  Psi_liq/p_liq; every existing consumer of rho_lR keeps rho_EOS. Anything
  else risks re-running the micro-vs-bulk density incident
  ([[ogs_rm_dsm_potential_physics]]).
- `invertDisjoiningPressure` and the equilibrium branch logic: unchanged.
  Equilibrium mode is already conservative on its manifold (envelope
  theorem); it interoperates with (E) only through the unloaded branch
  (below the branch point the exact kinematic pair applies if selected).
- Exchange law, storage coupling dphi_M/dn_l = -(1-phi_M), p_L_m fold:
  unchanged in form; p_L_m now sees the new mu (by construction, same fold
  point).

## 3. PRJ interface [DECISION Q5 for names; defaults fixed]

New optional tags inside `<potential_exchange>` (parse in
`CreateRichardsMechanicsProcess.cpp`, same pattern as `film_strain_coupling`):

```xml
<film_energy_route>operational | exact</film_energy_route>
    <!-- default: operational (bit-for-bit current behaviour).
         exact: requires film_strain_coupling == kinematic; create-time
         OGS_FATAL otherwise. -->
<micro_liquid_compression>off | kinematic</micro_liquid_compression>
    <!-- default: off. kinematic: requires film_strain_coupling == kinematic
         AND film_energy_route == exact; create-time OGS_FATAL otherwise. -->
<micro_liquid_bulk_modulus> (Parameter, Pa) </micro_liquid_bulk_modulus>
    <!-- REQUIRED iff micro_liquid_compression != off. NO default value —
         CLAUDE.md §1.1: the literal needs a cited source (§8 Q4). -->
```

Mode matrix (create-time validated; each invalid combination is OGS_FATAL
with a message naming this doc):

| film_strain_coupling | film_energy_route | micro_liquid_compression | behaviour |
|---|---|---|---|
| off         | operational | off       | bit-for-bit `7ff8861847` |
| kinematic   | operational | off       | bit-for-bit `7ff8861847` (operational cut) |
| kinematic   | exact       | off       | (E): exact pair, S = ½bK_d eps² |
| kinematic   | exact       | kinematic | (E)+(L): exact pair, S dropped, liquid carrier ON |
| equilibrium | operational | off       | bit-for-bit `7ff8861847` |
| equilibrium | exact / any | any≠off   | OGS_FATAL (not designed here) |
| off         | exact / any | any≠off   | OGS_FATAL |

## 4. Implementation plan (file by file)

Work top-down; after each numbered step the build MUST compile and ALL
pre-existing tests MUST pass (CLAUDE.md §6.5). Commit message pattern (§7):
"RM/DSM: <step>, per Vinay's 2026-06-11 spec (PI_OF_NL_EV_IMPLEMENTATION.md §<n>)".

1. **Enums + parameters** (`PotentialExchange.h`): add
   `enum class FilmEnergyRoute { Operational, Exact };`
   `enum class MicroLiquidCompressionMode { Off, Kinematic };`
   members in `PotentialExchangeParameters`:
   `FilmEnergyRoute film_energy_route = FilmEnergyRoute::Operational;`
   `MicroLiquidCompressionMode micro_liquid_compression = ...::Off;`
   `ParameterLib::Parameter<double> const* micro_liquid_bulk_modulus = nullptr;`

2. **New helper** (`PotentialExchange.h`, after `computeStrainedFilmState`):

   ```cpp
   struct StrainedFilmEnergyPairData
   {
       double Psi_film = 0.0;          // J/m^3 REV
       double sigma_sw_m = 0.0;        // Pa
       double mu_lR_mech = 0.0;        // J/kg
       double dsigma_sw_dnl = 0.0;     // Pa per unit n_l
       double dsigma_sw_deps_v = 0.0;  // Pa per unit strain
       double dmu_mech_dnl = 0.0;      // J/kg per unit n_l
       double dmu_mech_deps_v = 0.0;   // J/kg per unit strain
       double dmu_mech_drho_lR = 0.0;  // (J/kg)/(kg/m^3)
   };

   inline StrainedFilmEnergyPairData computeStrainedFilmEnergyPair(
       double const n_l, double const eps_v, double const kappa,
       double const biot_b, double const K_drained, bool const include_S,
       double const rho_lR, double const nS, double const rho_SR,
       double const hamaker_constant, double const specific_surface,
       double const potential_sign_factor,
       double const potential_augmentation_prefactor,
       double const potential_augmentation_exponent, double const n_l_floor);
   ```

   Implementation contract: evaluate the bare law TERM-WISE (vdW and aug
   separately — either by two calls with K_aug = 0 / A = 0, or by exposing
   the per-term values from the existing helper; prefer the two-call route, it
   needs no change to the shipped function); assemble §2.1 with the §4.6
   stable forms; every assignment carries a unit comment (§4.2). `include_S`
   selects route R3's S-term (false when (L) is ON).

3. **New helper** (same file):

   ```cpp
   struct MicroLiquidCompressionData
   {
       double Psi_liq = 0.0;       // J/m^3 REV
       double p_liq = 0.0;         // Pa
       double sigma_liq = 0.0;     // Pa
       double mu_liq = 0.0;        // J/kg
       double dsigma_liq_deps_v = 0.0, dsigma_liq_dnl = 0.0;
       double dmu_liq_deps_v = 0.0,  dmu_liq_dnl = 0.0, dmu_liq_drho_lR = 0.0;
   };
   MicroLiquidCompressionData computeMicroLiquidCompression(
       double const n_l, double const eps_v, double const kappa,
       double const K_liq, double const rho_lR, double const nS);
   ```

   Body per §2.2 — BLOCKED on [DECISION Q3] (route (a) vs (b)) — implement the
   approved variant only; eps_liq = -log1p(kappa*eps_v) (use log1p), reuse the
   1e-6 positivity guard convention of `computeStrainedFilmState`.

4. **Fold point** (`RichardsMechanicsFEM-impl.h`,
   `applyFilmPressureMicroPotential`, L680 block): inside the existing
   `film_strain_coupling != Off` branch, add the route switch:
   - `Operational`: existing code path, character-for-character.
   - `Exact` (kinematic only): call `computeStrainedFilmEnergyPair`
     (include_S = micro_liquid_compression == Off); set
     `out.mu_lR = mu_bare(w_eff at zero?)` — NO: in the exact route the FULL
     mu is `mu_bare(n_l) + mu_lR_mech` (the bare adsorption part stays
     evaluated at the TRUE n_l; the strain coupling lives in mu_lR_mech which
     already contains the w_eff physics through the integrals — mirror of the
     shipped Off-path structure `mu = mu_bare + mu_mech`, NOT of the
     operational w_eff-evaluation). Apply the macro-floor cutoff to the bare
     part exactly as the Off path does, then ADD (+=) mu_lR_mech and (if ON)
     mu_liq, with their d/dn_l, d/drho_lR folded the same way the shipped
     partner folds (`out.dmu_lR_dnl += ...` etc.). §4.1: ALL accumulations
     `+=`; never `=` on `out.mu_lR` after the bare part is in.
5. **Eigenstress site** (~L1689 + the local solves at L857/L1033/L1307 call
   chain): when route == Exact, sigma_sw uses `sigma_sw_m` (+ `sigma_liq` when
   (L) ON) from the SAME helper call — one evaluation per GP, both halves from
   one struct (no drift between halves). The existing eps_v threading
   (NaN-sentinel defaults) is reused unchanged.
6. **Jacobian blocks**: wire `dsigma_sw_dnl`, `dsigma_sw_deps_v`,
   `dmu_mech_deps_v` (+ liquid blocks) into the same slots the operational
   mode fills today (grep the consumers of the strained-mode derivatives in
   the local Newton and the global assembly; they are the same slots —
   the struct replaces the source of the numbers, not the wiring).
7. **Create-time parsing + validation** (`CreateRichardsMechanicsProcess.cpp`):
   parse the three tags; enforce the §3 mode matrix with OGS_FATAL.
8. **Tag documentation**: add the three files under the
   `Documentation/ProjectFile/...` subtree where `film_strain_coupling`'s
   t_*.md lives (grep; create `t_film_energy_route.md`,
   `t_micro_liquid_compression.md`, `t_micro_liquid_bulk_modulus.md`,
   one-paragraph each, citing this doc).

### 4.6 Numerically stable forms (MANDATORY — small-kappa/small-strain)

With x = kappa*eps_v, y = xi*kappa*eps_v:

```
J3(x)/eps_v:  for |x| > 1e-5:  [1 - 1/((1+x)*(1+x))] / (2*kappa*eps_v)
              else:            (1 - 1.5*x + 2.0*x*x) ... use series
                               J3 = eps_v*(1 - 1.5*x + 2.0*x*x)   // O(x^3)
Jx:           for |y| > 1e-5:  -std::expm1(-y) / (xi*kappa)
              else:            eps_v*(1 - 0.5*y + y*y/6.0)        // O(y^3)
eps_liq:      -std::log1p(kappa*eps_v)
```

Guard (1 + kappa*eps_v) >= 1e-6 exactly as `computeStrainedFilmState` does;
in the guarded regime freeze the strain derivatives to 0 (same convention).
The 1e-5 switch thresholds are series-accuracy bounds (relative error <
1e-10 at the switch), §1.2-scoped numeric defaults — flag in the PR, no
approval needed.

## 5. Tests (CLAUDE.md §3 — structure binding, values TODO(user))

New file `Tests/ProcessLib/RichardsMechanics/ExactFilmEnergyPair.cpp`
(supplements StrainedFilmPotential.cpp; existing 28 tests must stay green).
Representative state for all tests unless noted: take the SAME parameter set
the existing StrainedFilmPotential.cpp fixtures use (already §12-traceable on
the parent branch) — do NOT introduce new material literals.

```
T-0  Name:           SymbolicDerivativeGeneration (offline, not a gtest)
     What:           sympy script `scripts/dsm/derive_pi_fofnlev.py` (new)
                     generating the C++ for d(n_l*I_aug)/dn_l etc.; output
                     committed as a comment block above the implementation.
     Anchor:         (d) symmetry of the procedure (generated-vs-FD in T-3).

T-1  Name:           OffModeBitForBit
     Physics anchor: (f) regression baseline (dd1400 free swelling, parent
                     branch head 7ff8861847, user-approved baseline).
     Input config:   all new tags at defaults.
     Expected:       bitwise-identical field data, all output VTUs (recipe §6).
     Catches:        accidental default-path drift.

T-2  Name:           ExactPair_ZeroStrainReduction
     Physics anchor: (a) analytical limit.
     Input config:   eps_v = 0, kinematic+exact.
     Expected:       sigma_sw == -(1-phi_M)*n_l*Pi(n_l) and mu_mech == 0,
                     to FD-free analytic equality (tol = 1e-14 relative,
                     scale-derived: exact algebraic identity).
     Catches:        integral-assembly constant errors.

T-3  Name:           ExactPair_MaxwellIdentity
     Physics anchor: (d) symmetry (mixed-partial equality).
     Input config:   grid of (n_l, eps_v) states ON and OFF the guard,
                     kappa in {aggregate, unity}; FD vs analytic for ALL
                     derivative blocks AND the cross identity
                     dsigma_sw/dn_l == (1-phi_M)*rho_lR*dmu_mech/deps_v.
     Expected:       agreement to FD tolerance derived from scale
                     (tol = 1e-6 * |value| + 1e-6 * Pi-scale; in-file derived).
     Catches:        any transcription error in §2.1 derivatives.

T-4  Name:           ExactPair_FrozenHLimit
     Physics anchor: (a) analytical limit (kappa -> 0).
     Input config:   kappa = 1e-9 (sub-threshold series branch), compare to
                     computeIntegrableMechanicalMicroPotential at same state.
     Expected:       relative difference -> O(kappa) (assert <= C*kappa with
                     C documented in-file from the series remainder).
     Catches:        wrong limit structure (the operational form FAILS this
                     test by construction — discriminating).

T-5  Name:           ExactPair_ReversibilityLoop
     Physics anchor: (b) conservation law (∮dW = 0 for a gradient field).
     Input config:   closed rectangular loop in (eps_v, n_l), trapezoid
                     quadrature, N points; BOTH routes evaluated:
                     exact (expect ~0) and operational (expect != 0).
     Expected:       |∮dW|_exact <= quadrature error bound (derived in-file
                     from N and path scale); |∮dW|_operational ABOVE that
                     bound at the same N (the §9a predicted defect — this
                     test VERIFIES the prediction; report the measured value
                     in the PR per §5.1 'measured').
     Catches:        the entire point of the branch.

T-6  Name:           LiquidCarrier_EnergyPressureConsistency   [BLOCKED Q3/Q4]
     Physics anchor: (a) analytical limit (p = K*eps_liq from Psi_liq by FD).
     Input config:   FD of Psi_liq vs analytic p_liq, sigma_liq, mu_liq.
     Expected:       FD agreement (scale-derived tol as T-3). K_liq value:
                     TODO(Vinay) — test reads it from the fixture parameter,
                     no literal in the test body.
     Catches:        sign/chain errors in the carrier.

T-7  Name:           ModeMatrixValidation
     Physics anchor: (f) regression of create-time contract.
     Input config:   each invalid row of §3 matrix in a minimal PRJ.
     Expected:       OGS_FATAL with the documented message.
     Catches:        silent double-counting configurations.

T-8  Name:           ExpulsionProbe_DrainedRamp                [run-level]
     Physics anchor: (a) physical limit (loading a saturated film must
                     decrease n_l past the §2 crossover).
     Input config:   single-element drained oedometer ramp, kinematic+exact
                     (+liquid when Q3/Q4 resolved), past b*K_d > 3*kappa*Pi.
     Expected:       sign only: n_l monotonically decreasing past the
                     crossover. Magnitudes: TODO(Vinay).
     Catches:        reversal not realized at run level.
```

GTEST_SKIP() with TODO for anything blocked on Q3/Q4 — never a trivially
passing assertion (§3.5).

## 6. Verification recipe (commands; §6.8 naming applies to any campaign dirs)

```
# build (out of source, like h_of_eps_20260609)
cmake --preset release -B ~/git/build/Pi_fofnlev_$(date +%Y%m%d) \
      -S ~/git/ogs-worktrees/dsm_native_Pi_fofnlev_wt   # mirror the parent branch's preset
ninja -C ~/git/build/Pi_fofnlev_<date> testrunner ogs

# unit tests (all RichardsMechanics + the new file)
~/git/build/Pi_fofnlev_<date>/bin/testrunner --gtest_filter='*RichardsMechanics*'

# T-1 off-mode regression: two binaries, one input (the parent-branch binary
# ~/git/build/h_of_eps_20260609/bin/ogs is the reference)
#   run dd1400 free-swelling PRJ with both binaries into out_ref/ out_new/
#   compare: python3 -c "vtk field-data byte compare" (reuse the comparison
#   script from the 2026-06-09 verification, see strained-film doc §10)

# T-8 + MS33 VII A/B (only after Q3/Q4): campaign dir
#   ms33_VII_pifofnlev_<date>_<status>/  (status suffix MANDATORY, §6.8)
```

Gate for upgrading any [PRED] in this doc or the beamer to verified/measured:
the corresponding test or run above, cited by name+date (§5).

## 7. Documentation deliverables (same commit as the code, §8)

1. THIS file: §10 status log updated at every step (DONE + date + outcome).
2. `STRAINED_FILM_IMPLEMENTATION.md`: §9a annotated "exact route implemented
   on dsm_native_Pi_fofnlev (<date>), see PI_OF_NL_EV_IMPLEMENTATION.md;
   operational cut retained PRJ-selectable" — annotate, do not delete (§6.3).
   D2/D3 entries cross-referenced likewise.
3. `DSM/AGENTS.md`: worklog entry (created with this doc; update on DONE).
4. Memory: `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/`
   `project_dsm_pi_fofnlev.md` + MEMORY.md pointer (created with this doc;
   update status on verification).
5. Tag docs (§4.8).
6. **Beamer** `~/tex/implementation_and_theory/implementation/maxwell_from_psi.tex`
   — apply ONLY after the corresponding verification (§5 discipline; [PRED]
   stays [PRED] until then). Frame-level plan:
   - NOW (done with this doc): Step 19 footer D2 note + Step 23 open-list item
     point to this branch/doc as "in design".
   - After T-3/T-5 pass: Step 21 ("honesty") — retitle the operational cut as
     the PRJ-selectable fallback; the exact pair gets [DERIVED]+[MEASURED]
     loop-closure (quote |∮dW| from T-5); update Step 23 (§9a item -> DONE,
     date); Step 7b "Missing" row -> "supplied on dsm_native_Pi_fofnlev".
   - After T-8 / VII A/B: new frame "Step 24 — the compliant load path:
     first numbers" (table off/operational/exact(+liquid), [MEASURED], run
     dir cited); references frame gains the K_liq source (Q4).
   - Recompile must stay clean (current baseline: 131 pages, 0 overfull).

## 8. DECISIONS — Vinay's queue (block coding where marked)

- **Q1 — §9a closed forms (E).** The derivation in §2.1 (identical in
  substance to design doc §9a, now with full derivative blocks). APPROVE to
  unblock steps 4.2/4.4–4.6 for the (E) route.  [formulation — §9 ASK USER]
- **Q2 — the S-term carrier when (L) is OFF.** Default here: keep
  S = 0.5*b*K_drained*eps_v^2 (transmitted-load work, skeleton modulus — the
  R3 form, frozen-h-compatible). Alternative: no S at all in the exact route
  (pure geometric film energy; load reaches mu only via (L)). Affects T-4's
  limit statement.  [formulation]
- **Q3 — (L) mass-derivative freeze. PARKED (Vinay, 2026-06-11: "park it,
  this is getting very intricate").** Full decision brief, so the parked
  state is self-contained:
  - Route (a) — FIXED-VOLUME mass derivative: adding water at frozen film
    volume densifies the liquid; product rule on m*psi(m/V) gives
    mu_liq = psi_liq + p_liq/rho_lR — the FULL Derjaguin v*P term recovered
    from a Helmholtz energy. Magnitude at kappa*eps ~ 0.05, K_liq ~ 2 GPa:
    mu-jump ~ 9e4 J/kg — an enormous expulsion drive, exchange-relieved on
    the alpha_M timescale. Price: the freeze "V_film independent of n_l in
    the derivative" sits oddly against the B1 kinematics used everywhere
    else (V ∝ n_l).
  - Route (b) — KINEMATIC PINNING: film volume tracks water content (B1
    convention); d(eps_liq)/dn_l = 0 and mu_liq = psi_liq only =
    0.5*K_liq*eps_liq^2/rho_lR — NO p/rho term, quadratic, sign-symmetric
    in eps_v, ~kappa*eps/2 (~40x at 5% strain) smaller than (a). Gentler;
    arguably deletes the Derjaguin physics the carrier exists for.
  - They differ at FIRST order — a constitutive fork, not a correction.
    Agent's reading (not a decision): (a) is Derjaguin-faithful; (b)'s
    consistency with B1 is cleaner. OWNER'S CALL.
- **Q4 — K_liq value + source. PARKED (Vinay, 2026-06-11, same note).**
  §1.1: no literal without citation. CANDIDATE (labelled, awaiting
  approval): bulk-water isothermal bulk modulus K_T ≈ 2.18 GPa
  (kappa_T ≈ 4.58e-10 1/Pa at 25 °C; Kell 1975, J. Chem. Eng. Data
  20(1):97, via CRC/NIST) — a LOWER anchor; confined interlayer water is
  expected stiffer [PRED, no source yet — a confined value needs a thm-lit
  pass, separate task]. Implication of the candidate: at 5% film strain
  p_liq ≈ 110 MPa, i.e. 2–3x Pi (~40–48 MPa) on MS33 states — the liquid
  channel would DOMINATE the disjoining channel and relax only through the
  exchange (Damköhler-limited); Jacobian gains K_liq*kappa^2 stiffness
  [PRED: smaller time steps until drainage relaxes it].
  - COUPLING Q3<->Q4: route (b) makes K_liq much less consequential
    (quadratic-small mu-drive); route (a) + 2.18 GPa is the maximal-physics
    combination (full v*P at near-bulk stiffness — strongest expulsion,
    stiffest numerics). A cautious first sweep would be (a) with K_liq
    scanned over a decade below the candidate in the T-8 probe — sweep
    values themselves need owner approval before any run.
  - UNPARK BY: answering Q3 (a|b) and Q4 (approve candidate | other value+
    source | order thm-lit pass | drop (L)). Until then (E) stands alone:
    the exact pair is Maxwell-closed with the R3 transmitted-load term
    carrying the load coupling; tests T-6/T-8 stay GTEST_SKIP.
- **Q5 — tag names/defaults** (§3): `film_energy_route`,
  `micro_liquid_compression`, `micro_liquid_bulk_modulus`. Rename freely;
  defaults (off/operational, bit-for-bit) are fixed by §6.5 discipline.

## 9. Guardrail notes for the implementing agent

- §1.1 fires on: K_liq (Q4 — never hard-code, Parameter only), any new
  tolerance not derived in-file, any sign not taken from the existing code.
- §2: never assert on a quantity derived from a parameter fitted in the same
  test; T-5's operational-vs-exact contrast is the discriminating design.
- §5: every runtime claim in comments labelled expected/measured; this doc's
  [PRED] markers upgrade only per §6 gate.
- §6.2/6.3: no deletion of PRJs/meshes/.md; §6.5 every commit compiles + green.
- §7: no authorship claims; no Co-Authored-By trailer.
- §12: any NEW PRJ used for T-8/VII carries the full provenance header
  (inherit the parent PRJ's block, add the new tags with "this doc §3" as the
  source for structural defaults; K_liq cites the Q4-approved source).

## 10. Status log

- 2026-06-11: branch `dsm_native_Pi_fofnlev` + worktree created off
  `7ff8861847`; this design doc written (agent, per Vinay's commission of
  2026-06-11); DSM/AGENTS.md worklog entry added; memory file
  `project_dsm_pi_fofnlev` created; beamer Step 19/23 "in design" markers set.
- 2026-06-11 (later): **(E) IMPLEMENTED** per Vinay's "implement that now"
  (read as resolving Q1, Q2-default keep-S, Q5-default names; (L) remains
  blocked on Q3+Q4). Code: `FilmEnergyRoute` enum + parameter + pure
  validation predicate (PotentialExchangeParameters.h),
  `computeStrainedFilmEnergyPair` (PotentialExchange.h; mu-space closed forms
  with x_over_kappa-stable small-kappa evaluation), exact-route fold branch
  in `applyFilmPressureMicroPotential` (FEM-impl.h; bare at TRUE n_l +
  one-Psi partner, cutoff factor g product-ruled, eigenstress site unchanged
  by design), `film_energy_route` parsing + create-time mode-matrix OGS_FATAL
  (CreateRichardsMechanicsProcess.cpp), tests
  `Tests/ProcessLib/RichardsMechanics/ExactFilmEnergyPair.cpp` (T-2..T-5,
  T-7 helper-level; T-6/T-8 GTEST_SKIP blocked on Q3/Q4), regression tool
  `scripts/dsm/compare_vtu_bitwise.py`.
  VERIFIED 2026-06-11 (build `~/git/build/Pi_fofnlev_20260611`, Release):
  - 36/36 RichardsMechanics unit tests pass (28 pre-existing + 8 new; 2 of
    the new are the designed Q3/Q4 skips).
  - T-3 Maxwell identity: analytic-vs-analytic cross identity to 1e-9
    relative; all derivative blocks FD-confirmed.
  - T-4 frozen-h limit: kappa = 1e-9 reduces to the shipped integrable
    partner within the series bound (10*kappa relative).
  - T-5 loop (measured, N=400 trapezoid): |W_exact|/scale = 8.4e-9 (below
    the 50/N^2 quadrature bound, halving-convergent); |W_operational|/scale
    = 0.93. **The §9a "small ∮ != 0" prediction is CORRECTED by measurement:
    on a Pi-dominant loop with augmentation the operational defect is O(1)
    of the path work, not small** (§5 update applied to the strained-film
    doc as well).
  - T-1 off-mode regression: dd1400
    (`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj`)
    run with the parent-head binary (`h_of_eps_20260609`, 7ff8861847) vs the
    new binary, all defaults: ALL 12 output VTUs bitwise-identical field
    data (points + every point/cell array; two binaries, one input).
  PENDING: (L) implementation (Q3+Q4, PARKED); T-8 expulsion probe + MS33
  VII A/B (run-level); Doxygen tag doc for `film_energy_route` (joint TODO
  with the film_strain tags).
- 2026-06-11 (later still): committed + pushed (origin, vgk2). REBASED onto
  `23a723cc3c` (live K(rho_d)) per Vinay's merge order; re-verified after
  rebase: 39/39 RichardsMechanics unit tests (28 base + 3 live-K + 8 here,
  2 designed skips), dd1400 off-mode bitwise-identical vs the 23a723cc
  reference binary (h_of_eps_20260609, rebuilt-in-place build of
  2026-06-11). `dsm_native_h_of_eps` retired (branch ref deleted; history
  fully contained here).
