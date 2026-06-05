# The (1−φ_M) conjugate referencing — energy-conservation corrigendum to B1

**Branch:** `dsm_native_pdisj_maxwell`.
**Corrects:** `MAXWELL_CONJUGATE_IMPLEMENTATION.md` §2 eq.(★) and `PotentialExchange.h`
`computeMaxwellConjugateMicroPotential` (L257). The B1 conjugate term is wired but
its **magnitude is under-referenced by one factor of (1−φ_M)** — it does not form a
true conjugate pair with the swelling eigenstress, so the `(σ, μ_lR)` response is
**non-conservative** (`∮ ≠ 0`). This file is the spec to close that.

Tags: **[D]** derived · **[DECISION]** Vinay's call · **[PRED]** predicted, not verified.

---

## 0. The catch in one line

The two halves of the Maxwell pair live at **different scales**, and only one was
referenced to the REV:

| channel | quantity | as-coded | REV-referenced? |
|---|---|---|---|
| micro→macro (swelling) | `σ_sw = −φ_m·Π = −(1−φ_M)·n_l·Π` | FEM-impl L1628 | **yes** — carries `(1−φ_M)` |
| macro→micro gate | `p′ ≥ φ_m·Π = (1−φ_M)·n_l·Π` | FEM-impl L2744/3607 | **yes** — opt.1, carries `(1−φ_M)` |
| macro→micro **magnitude** | `μ_lR_mech = S₁·ε_v / ρ_lR` | PotentialExchange.h **L257** | **NO** — missing `1/(1−φ_M)` |

The fix is one reciprocal factor: **swelling carries `×(1−φ_M)` (down), so its
conjugate must carry `÷(1−φ_M)` (up).** Get it wrong and the pair is not the gradient
of any single `Ψ`.

## 1. Why reciprocal — the energy argument [D]

The coupling free energy (volumetric sector, per **REV** volume) is bilinear,

```
Ψ_cpl(ε_v, m_l) = σ_sw,m · ε_v ,
```

with the two physically meaningful conjugate variables: REV strain `ε_v` (work-conjugate
to the macro stress) and the **micro liquid mass per REV volume**

```
m_l = ρ_lR · φ_m = ρ_lR · (1−φ_M) · n_l .            (mass referencing)
```

`μ_lR` is a specific potential **[J/kg]** — its conjugate is `m_l`, not `n_l`. For
`Ψ_cpl` to be a state function (so a closed `(ε_v, m_l)` loop returns zero net work,
`∮ dΨ_cpl = 0`), the two partials must satisfy the Maxwell relation

```
∂σ_sw,m/∂m_l  =  ∂μ_lR_mech/∂ε_v .                   (integrability)
```

Differentiate `Ψ_cpl` both ways, using `∂m_l/∂n_l = ρ_lR·(1−φ_M)` (φ frozen, B1):

```
∂Ψ_cpl/∂ε_v = σ_sw,m(n_l)                            # eigenstress — LEFT, present
∂Ψ_cpl/∂m_l = (∂σ_sw,m/∂n_l) / (∂m_l/∂n_l)
            = S₁ / [ρ_lR·(1−φ_M)]                    # the per-REV-mass conjugate
```

so the correct micro-potential term is

```
   μ_lR_mech = S₁ · ε_v / [ ρ_lR · (1−φ_M) ] ,     S₁ ≡ ∂σ_sw,m/∂n_l = −n_S·(Π + n_l·dΠ/dn_l)   (★′)
```

**This is the existing (★) scaled UP by `1/(1−φ_M)`.** Equivalently: the micro feels
the REV effective stress *concentrated* onto the load-bearing aggregate fraction
`(1−φ_M)`, `σ′_felt = σ′ / (1−φ_M)` — only the aggregate carries the skeleton stress,
the macro pores carry fluid. The swelling stress was *diluted* to the REV by the same
`(1−φ_M)`; the work product `[σ_sw down] · [σ′ up]` is invariant ⇒ energy conserved.

**Why the current code violates it.** `mu_lR_mech = S₁·ε_v/ρ_lR` (L257) references the
conjugate to `ρ_lR·n_l` (the *aggregate*-scale mass), not the REV mass
`ρ_lR·(1−φ_M)·n_l`. Then `∂σ_sw,m/∂m_l = S₁/[ρ_lR(1−φ_M)]` but
`∂μ_lR_mech/∂ε_v = S₁/ρ_lR` — the two cross-partials **differ by `(1−φ_M)`**, the
integrability condition fails, `∮ dΨ_cpl = ρ̂·(...)·(1−φ_M − 1) ≠ 0`. The eigenstress
(LEFT) is REV-scaled; the conjugate (RIGHT) is aggregate-scaled. They are not partners.

> The existing doc's "Maxwell closes by construction: `ρ_lR·∂μ_lR/∂ε = S₁`" is closure
> against the measure `ρ_lR·n_l`. The thermodynamically correct measure is the per-REV
> micro liquid mass `ρ_lR·(1−φ_M)·n_l`, because the exchange residual
> `rhs += N_pᵀ·ρ̂·w` (FEM-impl L2762) is assembled at the **REV** scale. Hence the
> extra `(1−φ_M)`.

## 2. The fix — one factor, three sites [D]

**Helper `computeMaxwellConjugateMicroPotential` (PotentialExchange.h L257–260).**
Divide by the **REV micro-liquid mass density** `ρ_lR·(1−φ_M)` instead of `ρ_lR`.
Cleanest: pass the referencing fraction in, keep the helper pure.

```
// add an argument n_S (= 1 − φ_M); replace rho_lR by rho_lR * n_S in the measure
out.mu_lR_mech       = S1 * eps_v / (rho_lR * n_S);          // was / rho_lR
out.dmu_lR_mech_deps_v = S1 / (rho_lR * n_S);                // was / rho_lR
out.dmu_lR_mech_dnl  = dS1_dnl * eps_v / (rho_lR * n_S);     // was / rho_lR
out.dmu_lR_mech_drho_lR = -out.mu_lR_mech / rho_lR;          // unchanged in form (φ frozen ⇒ ∂n_S/∂ρ_lR = 0)
```

**Both call sites (FEM-impl L2745, L3608)** pass `n_S_mc = 1 − φ_M` (already in scope
as `n_S_mc`, the same fraction used in `S1_mc` and `Pi_gate_mc`). No new state, no new
parameter.

**Do NOT touch** the eigenstress (L1628) or the gate (L2744/3607) — both already carry
`(1−φ_M)`. The asymmetry is entirely in the magnitude denominator.

## 3. Consistent tangent + the symmetry test [D]

The Jacobian block `∂ρ̂/∂ε = α_M·∂μ_lR_mech/∂ε·mᵀ` (MAXWELL_CONJUGATE_IMPLEMENTATION
§3) inherits the same `1/(1−φ_M)`:

```
∂μ_lR_mech/∂ε_v = S₁ / [ρ_lR·(1−φ_M)] .
```

**Maxwell-symmetry GP unit test — state it in the REV metric.** The old test
`ρ_lR·∂μ_lR/∂ε_v == ∂σ_sw,m/∂n_l` must become

```
ρ_lR·(1−φ_M)·∂μ_lR_mech/∂ε_v  ==  ∂σ_sw,m/∂n_l   (= S₁)      [D, calibration-free]
```

With (★′) this is an identity; with the current L257 it fails by `(1−φ_M)`. This is the
falsifiable signature of the bug — run it first, it needs no model run.

**Closed-loop ∮ check [D].** Drive a small cycle in `(ε_v, n_l)` that opens the gate
and returns. Net work `∮ (σ_sw dε_v + μ_lR_mech dm_l)` must be `0` to tolerance with
(★′); with L257 it is `O((1−φ_M)−1)·` (loop area) ≠ 0. This is the direct
energy-conservation test the user's argument names.

## 4. Impact [PRED — verify, do not assert]

- **Term is larger.** `(★′)/(★) = 1/(1−φ_M)`. For the ANCHORS densities
  `φ_M ≈ 0.1–0.4` ⇒ the conjugate (and the expulsion it drives) is **~1.1–1.7× stronger**
  than the wired value. The gate is unchanged, so *when* it fires is unchanged; *how hard*
  it drains is up by `1/(1−φ_M)`.
- **Re-run the confined-expulsion probe** (`ms33_confined_expulsion_dd1400.prj`) after the
  fix — with the gate at `φ_m·Π` (opt.1) and the corrected magnitude, show `n_l` drops
  once `p′ ≥ φ_m·Π`. Previously the gate never opened; with opt.1 it does (~7.8 MPa).
- **K re-calibration** vs Dixon still pending per spec §6 — the stronger conjugate
  re-touches the saturated equilibrium past the gate (not below it, where the term is
  inert and `K` was anchored). Never fit-and-verify in one step.
- **[PRED] dd1600 `#104` SparseLU singular wall** — restoring a *true* transpose pair
  (now symmetric in the REV metric, not just `∝S₁`) is a better regulariser candidate
  than the under-referenced version. Re-run to check; relabel only after.

## 5. Provenance / decisions

- **[D]** The `(1−φ_M)` factor and its direction (swelling ÷ down, conjugate × up) are
  derived from the mass referencing `m_l = ρ_lR(1−φ_M)n_l` and the REV-scale assembly of
  the exchange residual — no new constant, no calibration.
- **[DECISION — Vinay, this session]** The reciprocal `(1−φ_M)` referencing is required
  for energy conservation; the conjugate magnitude must be scaled up by `1/(1−φ_M)` to
  partner the REV-referenced swelling eigenstress.
- Sign, gate (opt.1), φ-freeze (B1), stress measure: unchanged from
  MAXWELL_CONJUGATE_IMPLEMENTATION §6.

## 6. Cross-refs

- `MAXWELL_CONJUGATE_IMPLEMENTATION.md` — the B1 design this corrects (its (★) and L257).
- `MAXWELL_PAIR_RESTORATION.md` — the goal spec (claim B).
- `paper_DSM.tex` §2.4 `sec:admissibility` — the swelling-prefactor note
  (`(1−φ_M)·φ_m`, the same REV referencing on the eigenstress side).
- `truesdell_noll_lecture.tex` pp. 77–82 — the Maxwell-pairs unit (`∮ ≠ 0` = dropped
  conjugate); this file is its quantitative `(1−φ_M)` instance.
