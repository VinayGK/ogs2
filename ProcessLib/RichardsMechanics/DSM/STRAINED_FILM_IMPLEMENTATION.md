# Strained-film disjoining law — h(w_m, ε_v) — IMPLEMENTATION DESIGN

**Branch:** `dsm_native_h_of_eps` (off `dsm_native_maxwell_conjugate` @ `d98f5f8324`).
**Builds on:** `MAXWELL_PAIR_RESTORATION.md`, `MAXWELL_CONJUGATE_IMPLEMENTATION.md`,
`MAXWELL_CONJUGATE_REV_REFERENCING.md`.
**Goal owner:** Vinay (2026-06-09): *"if you make the thickness of the layer also
varying with increasing compressive stress, the potential has to be reversed.
There is no way that there is no admissible constitutive model for this. It is a
reversible process — not fully, but definitely reversible."*

Tags: **[D]** derived/established · **[DECISION]** Vinay's call before coding ·
**[PRED]** predicted, not verified.

---

## 0. Scope and goal

Reversible (conservative, one-Ψ) de-swelling under mechanical load. The shipped
integrable partner (`computeIntegrableMechanicalMicroPotential`) imbibes at all
loads because the film geometry is frozen (B1); the only expulsion channel in
the code history is the retired, non-integrable film bolt-on (`∮ ≠ 0`). This
design adds a strain-carrying film thickness h(w_m, ε_v) so that **both** the
swelling eigenstress and a load-reversing micro potential derive from one
Ψ_film. Out of scope: the dissipative residual (fabric-collapse hysteresis,
"not fully reversible") — that stays a future flow rule (B2-type).

Two h-laws are implemented, **PRJ-selectable**, to test against each other:
- **Variant A — kinematic spacing:** h follows the volumetric strain.
- **Variant B — equilibrium spacing:** h tracks the film force balance
  Π(h) = p_conf once the load exceeds the unloaded disjoining pressure.

## 1. The code as read [D]

- `PotentialExchange.h` `computeVanDerWaalsMicroPotential(n_l, rho_lR, nS,
  rho_SR, A, Sa, sign, K, lambda, dnS_dnl, n_l_floor)` — **no strain or stress
  argument**. Film thickness `h = n_l/(nS*rho_SR*Sa)` (comment L91): water
  content only.
- Shipped integrable partner: `mu_lR_mech = -[(Pi + n_l*Pi')*eps_v
  + 0.5*b*K_drained*eps_v^2]/rho_lR`; slope `(2*Pi + b*p_conf)/rho_lR > 0` for
  the vdW cubic core ⇒ **imbibition at every load**.
- Retired expulsion channels: `computeMaxwellConjugateMicroPotential` (sharp
  gate, B1) and the film-pressure delta `+g*b*p_conf/rho_lR` — both expel, both
  non-integrable as wired (`∮ ≠ 0`), both retired/superseded 2026-06-08.
- Micro liquid EOS (`FEM-impl.h` L503–513): `rho_lR = rho_l0*exp(-a_rho *
  omega_l^b_rho) + rho_LR`, `omega_l = n_l*rho_lR/(nS*rho_SR)` — **packing-
  driven, no pressure and no strain argument**. The EOS therefore cannot carry
  a load-pressurization channel in its current form.
- Eigenstress site: `FEM-impl.h` ~L1689 (`sigma_sw.sigma_sw += ...`), with the
  WIP closure `sigma_sw,m = -phi_m*p_film = -(1-phi_M)*n_l*(Pi - b*p_conf)`.

## 2. Structural finding — where the reversal can and cannot live [D]

Write the operational micro potential as μ_lR = −Π(h)/ρ_lR with Π > 0 the
(repulsive) disjoining pressure, and let the squeeze act only through the
geometry, h = h(w_m, ε_v) with ∂h/∂ε_v > 0 (compression thins the film). Then

```
∂μ_lR/∂ε_v = -(1/ρ_lR) · Π'(h) · ∂h/∂ε_v  > 0      for ANY repulsive law
                                                     (Π'(h) < 0 always),
```

i.e. **a pure geometric squeeze of a repulsive disjoining law can never reverse
the potential** — squeezing raises Π, the potential −Π/ρ drops, the film pulls
water IN (it wants to re-thicken). This holds for the vdW cubic core
(contribution +3Π_vdW/ρ per unit strain) and for the exponential augmentation
(+ξ·Π_aug/ρ, ξ = h/λ) alike. [D]

The reversal therefore lives in the **load-transmission (over-pressure) term**:
squeezing confined liquid raises its chemical potential by v·Δp (Derjaguin film
thermodynamics: μ = μ₀ + v·(P − Π(h))). That is exactly the physics the retired
film form carried (+b·p_conf/ρ_lR) — its sin was never the sign, it was
non-integrability at frozen h. **The fix is not "strain the Π law instead of the
load term"; it is "keep both, let h be live, and derive both halves from one
Ψ_film so the pair closes."** With h live, the load term and the Π(h) stiffening
combine; the slope

```
ρ_lR · ∂μ_lR/∂ε_v = -Π'(h)·∂h/∂ε_v·ρ-units + b·∂p_conf/∂ε_v-route
                  ~  +3κΠ(h)  -  b·K_drained          (cubic core, drained line)
```

flips sign when `b·K_drained > 3κ·Π(h)` — an **emergent, smooth gate**: stiff
skeleton / soft film ⇒ expulsion; soft skeleton / hard film ⇒ imbibition. No
Macaulay bracket, no bolted-on switch. [D — symbolic; the precise gate locus
depends on the h-law, §4–5.]

## 3. Common framework — one Ψ_film, differentiate twice [D]

Do NOT postulate the two halves separately (that is what broke integrability
twice). Postulate the **energy**, then read both halves off:

```
Ψ(ε, w_m) = Ψ_skel(ε) + Ψ_film(ε_v, w_m),

σ_sw,m  = ∂Ψ_film/∂ε_v                       (eigenstress half)
μ_lR    = (1/((1-φ_M)·ρ_lR)) ∂Ψ_film/∂w_m    (micro-potential half;
                                              REV referencing per
                                              MAXWELL_CONJUGATE_REV_REFERENCING.md)
```

Maxwell then holds **identically** — no separate check, no retro-fit. The two
variants are two constitutive choices of the film state inside Ψ_film.

Referencing constraint [D]: for ∂Ψ_film/∂ε_v to reproduce the existing
eigenstress scale −φ_m·Π at ε_v = 0, the strain rate of the spacing must carry
the aggregate fraction: ∂h/∂ε_v = κ·h₀ with **κ = (1−φ_M)** (the integrable
completion of the shipped closure). κ = 1 (naive "spacing follows REV strain
one-to-one") rescales the eigenstress by 1/(1−φ_M) — the same REV-vs-aggregate
referencing trap as the corrigendum. → DECISION D1.

## 4. Variant A — kinematic spacing

```
h_A(w_m, ε_v) = h₀(w_m) · (1 + κ·ε_v),    h₀ = w_m/((1-φ_M)·ρ_SR·S_a)
```

Equivalently: evaluate the EXISTING law at the effective content
`w_eff = w_m·(1 + κ·ε_v)` (the law depends on w_m only through h), keeping the
mass bookkeeping at w_m. Implementation = a thin wrapper that computes w_eff and
chain-rules the existing derivatives — minimally invasive, byte-identical at
ε_v = 0 or κ = 0.

Energy and the two halves (cubic core shown; augmentation analogous):

```
Ψ_film^A = -φ_m·[ E_film(h_A) + b·(p_conf-coupling, ½·b·K_d·ε_v² on the
            drained line) ]      → σ_sw,m^A = -φ_m·(Π(h_A) - b·p_conf)
                                 → μ_lR^A   = -(Π(h_A)·(1+κε_v) - b·p_conf)/ρ_lR
                                              + O(κ·ε_v) cross-terms (exact
                                              forms to be finalized in code
                                              review; ALL from ∂Ψ^A, no hand
                                              terms)
```

- Expulsion: through the −b·p_conf channel once `b·K_d > 3κΠ(h_A)` (§2);
  Π(h_A) itself stiffens under compression, sharpening the swelling drive below
  the gate and the snap above it — **smooth crossover**, not a sharp gate.
- At fixed w_m the squeeze implies the liquid densifies (mass conservation in a
  thinner film). The current packing-EOS cannot see this (no strain argument);
  variant A therefore treats ρ_lR as the EOS gives it (per ω_l) and carries the
  strain ONLY in the potential. The (deferred) consistent step would add ε_v to
  ω_l — flagged as A2, out of scope. → DECISION D2.

## 5. Variant B — equilibrium spacing

```
h_B(w_m, p_conf) = min( h₀(w_m),  h_eq(p_conf) ),   Π(h_eq) ≡ p_conf
```

The film sits on its force balance whenever the load can compress it; below
that load it is the unloaded h₀(w_m). For the cubic core the inversion is
closed-form: `w_eff = w_m·(Π(w_m)/p_conf)^(1/3)` on the loaded branch
(p_conf > Π(w_m)); the augmentation requires a scalar Newton (bounded, smooth,
local) or evaluation at the cubic-core inverse — implementation choice flagged
in §7.

- The gate is **emergent**: the min() switches branch exactly at
  p_conf = Π(h₀(w_m)) — the unloaded disjoining pressure of the current water
  content. No gate-scale decision needed (it is the film's own balance; note it
  is the MICRO-scale Π, the §7.2 scale question dissolves because the branch
  point is a property of the film law, not a REV comparison).
- Conservativeness: h_B is the minimizer of the film energy under the load
  constraint; by the envelope theorem the reduced Ψ̃(ε_v, w_m) =
  Ψ_film(h_B(·), ·) is still a single potential — derivatives of Ψ̃ ARE the
  constitutive pair; ∮ = 0 on the C¹ manifold. The branch point is C⁰/C¹ —
  kink, no jump (cf. the sharp-cavitation-handover discussion in the Tuller
  closure: branches sequential, not summed).
- On the loaded branch μ_lR^B = −(Π(h_B) − b·p_conf)/ρ_lR with Π(h_B) = p_conf:
  the disjoining part cancels against the load up to the b-fraction —
  μ_lR^B = −(1−b)·p_conf/ρ_lR; for b = 1 the film is exactly equilibrated
  (μ = bulk) and all further load goes to expulsion through the exchange. [D]

## 6. PRJ selection [DECISION D4]

New optional tag in `<potential_exchange>`:

```xml
<film_strain_coupling>off | kinematic | equilibrium</film_strain_coupling>
```

- `off` (default): bit-for-bit current behavior (shipped integrable partner
  stays active).
- `kinematic` (A) / `equilibrium` (B): the strained-film Ψ replaces BOTH the
  bare-law evaluation point AND the shipped `computeIntegrableMechanicalMicroPotential`
  term — **mutually exclusive, enforced at create-time** (double-counting
  guard: the shipped μ_mech is the frozen-h limit of the same physics; running
  both would count the strain coupling twice). → DECISION D3.

## 7. Implementation plan (code sites)

1. `PotentialExchange.h`: new helper `computeStrainedFilmState(w_m, eps_v,
   p_conf, Pi_law, mode, kappa)` returning `{w_eff, dw_eff_deps_v, dw_eff_dwm,
   dw_eff_dpconf, branch}`; the existing `computeVanDerWaalsMicroPotential` is
   then evaluated at `w_eff` with chain-ruled derivatives. No change to the
   bare law itself.
2. Eigenstress site (`FEM-impl.h` ~L1689): evaluate σ_sw with the SAME strained
   state (consistency — both halves of Ψ_film see one h).
3. Exchange assembly sites: μ_lR from the strained law; remove the
   `computeIntegrableMechanicalMicroPotential` add when mode ≠ off (§6).
4. Jacobian: new ∂μ/∂ε_v and ∂σ_sw/∂w_m blocks from the chain rules — they are
   transposes by construction (one Ψ); assert in the GP test.
5. `CreateRichardsMechanicsProcess.cpp`: parse the tag, validate exclusivity.
6. Docs: this file, AGENTS.md entry, ProcessLib docs page for the tag.

## 8. Verification (anchors per CLAUDE.md §3 — structure only, no values)

1. **off-mode regression** [approved baseline]: bit-for-bit vs branch head on
   dd1400 free swelling.
2. **ε_v = 0 reduction** [analytical limit]: kinematic and equilibrium modes
   reproduce the off-mode state exactly at zero strain / below the B branch
   point.
3. **Maxwell GP test** [derived identity]: FD-vs-analytic
   ∂σ_sw/∂w_m = (1−φ_M)·ρ_lR·∂μ_lR/∂ε_v in BOTH modes, on and off the loaded
   branch.
4. **Reversibility loop** [conservation law]: closed (ε_v, w_m) load-unload
   loop; ∮dW = 0 to tolerance in both modes (the property the retired film
   form fails — discriminating).
5. **Expulsion probe** [physical limit]: drained oedometer ramp past the §2
   crossover (A) / branch point (B): w_m must DECREASE; blind to K. Expected
   magnitudes: TODO(Vinay).
6. **K re-calibration** [PRED]: saturated swelling-pressure equilibrium shifts
   in both modes; re-anchor vs Dixon afterwards — separate step, never combined
   with 1–5.

## 9. DECISIONS — resolved (Vinay, 2026-06-09)

- **D1 — κ (spacing-strain weighting): BOTH, PRJ-selectable.**
  `<film_strain_kappa>aggregate|unity</film_strain_kappa>`, default
  `aggregate` (κ = 1−φ_M, the integrable completion of the eigenstress
  scale). κ frozen at the GP (B1) — no dκ chains.
- **D2 — EOS coupling: deferred (A2).** Variant A carries strain in the
  potential only this iteration; strained ω_l(ε_v) revisited only if A shows
  promise.
- **D3 — replacement: PROVISIONAL (Vinay: "not sure").** Implemented as
  replace-when-ON on mechanism-ownership grounds (both terms claim the strain
  coupling of the film; the global no-double-count rule bars running both).
  NOTE the original "shipped-limit" demonstration plan is RETIRED: the
  implemented operational form (§9a) does NOT reduce to the shipped integrable
  partner at small ε_v (slopes differ: +3κΠ/ρ vs +2Π/ρ on the vdW core at
  zero load) — they are genuinely different constitutive objects. The
  exclusivity is instead asserted structurally in the unit test
  (`ReplacementIsExclusiveAtZeroStrain`). Revisit when §9a is resolved.
- **D4 — PRJ tag:** `<film_strain_coupling>off|kinematic|equilibrium</...>`
  inside `<potential_exchange>`, default `off`.

## 9a. HONESTY NOTE — the implemented cut is the Derjaguin operational form,
NOT yet the exact one-Ψ pair [D]

What is implemented (both halves through one strained state w_eff):

```
mu_lR     = mu_bare(w_eff) + b*p_conf/rho_lR          (fold point)
sigma_sw  = -phi_m*(Pi(w_eff) - b*p_conf)             (eigenstress site)
```

This is variant-faithful (Derjaguin: μ = μ₀ + v(P − Π(h)), with h live) and
both halves see the SAME w_eff — but checking the cross-partials symbolically
shows the pair is NOT exactly the gradient of a single Ψ (the Maxwell defect
is O(Π·ε_v) terms, the same class the shipped partner repaired at frozen h).
The EXACT pair requires the energy route:

```
Ψ_film(ε_v, w_m) = ∫₀^{ε_v} σ_sw,m(w_m, e) de
                 = -(1-φ_M)·w_m·[ ∫₀^{ε_v} Π(w(w_m,e)) de + ½·b·K_d·ε_v² ]
```

whose ∫Π de is CLOSED-FORM for both terms of the law under the kinematic
h-law (vdW core: ∫(1+κe)⁻³de = [1−(1+κε)⁻²]/(2κ); exponential augmentation:
∫e^{−ξ(1+κe)}de = e^{−ξ}(1−e^{−ξκε})/(ξκ)), so the exact μ_mech =
(1/((1−φ_M)ρ_lR))·∂Ψ_film/∂w_m is implementable analytically. That derivation
is the NEXT step and a formulation change — Vinay to review the closed forms
before they are coded. Until then: do NOT cite this branch as "Maxwell-exact";
the reversibility-loop verification (§8.4) is expected to show a small ∮ ≠ 0
of the operational form, which the exact pair will close. [PRED]

> UPDATE 2026-06-11: the exact route IS now implemented on the follow-up
> branch `dsm_native_Pi_fofnlev` (design doc PI_OF_NL_EV_IMPLEMENTATION.md;
> `film_energy_route = exact`, PRJ-selectable; this branch's operational cut
> stays the default). MEASURED there (loop test, N=400): |∮dW|/scale =
> 8.4e-9 for the exact pair vs 0.93 for the operational form — the "small
> ∮ ≠ 0" wording above is corrected by measurement: on a Pi-dominant loop
> with augmentation the operational Maxwell defect is O(1) of the path work,
> NOT small. The [PRED] is thereby resolved: closure verified for the exact
> pair, defect measured (and large) for the operational cut.

## 10. Status

- 2026-06-09: branch + worktree created; design doc written; D1–D4 resolved
  (D3 provisional); IMPLEMENTED: enums/params, computeStrainedFilmState +
  invertDisjoiningPressure (PotentialExchange.h), fold-point rewiring
  (applyFilmPressureMicroPotential), eigenstress threading (eps_v/eps_v_prev
  sentinel args through computeSwellingStressIncrement), PRJ parsing, unit
  tests (Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp).
  VERIFIED 2026-06-09:
  - testrunner builds clean; all 28 RichardsMechanics unit tests PASS
    (21 pre-existing + the 7 new StrainedFilmPotential tests).
  - dd1400 free-swelling off-mode regression PASS: same PRJ run with this
    build vs the unmodified branch head (d98f5f8324, built fresh in
    mc_20260608) — all 12 output VTUs have BITWISE-IDENTICAL field data
    (points + every point/cell array); the only byte difference is the
    embedded OGS_VERSION string (different git-describe between builds).
    Valid isolation: two distinct binaries, one input.
  PENDING: the §9a exact-Ψ closed forms (Vinay review); kinematic/equilibrium
  expulsion probe (§8.5); reversibility loop (§8.4); K re-calibration [PRED].
