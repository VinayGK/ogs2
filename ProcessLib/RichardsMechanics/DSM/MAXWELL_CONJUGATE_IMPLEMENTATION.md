# Maxwell-conjugate term — IMPLEMENTATION DESIGN (B1)

**Branch:** `dsm_native_maxwell_conjugate` (off `a6dfec842c`).
**Builds on:** `MAXWELL_PAIR_RESTORATION.md` (the goal spec). That file is *what
and why*; this file is the concrete algebra + code plan for the **B1 milestone** —
restore the parameter-free conjugate term against the current sharp-`Π`
eigenstress closure.

Tags: **[D]** derived/established · **[DECISION]** Vinay's call before coding ·
**[PRED]** predicted, not verified.

---

## 0. Scope
Spec established the broken pair: LEFT `σ_sw = −φ_m·Π(n_l)` present, RIGHT
`∂μ_lR/∂ε = 0` absent. This file derives the missing RIGHT term against the real
code (`PotentialExchange.h`; the FEM-impl Π-eigenstress block) and lays out the
Jacobian, the code sites, and the verification.

## 1. The code as read [D]
`PotentialExchange.h`: `μ_lR = sign·(μ_vdW + μ_aug)`,
`μ_vdW = (A·Sa³/6π)·(nS³ρ_SR³)/(n_l³ρ_lR)`; micro pressure `p_L_m = −ρ_lR·μ_lR`;
disjoining `Π = −ρ_lR·ψ_Micro = p_L_m`. **No stress argument** in
`computeVanDerWaalsMicroPotential(...)` → `∂μ_lR/∂ε = 0`. Swelling eigenstress
(FEM-impl): `σ_sw = −φ_m·Π·I` (spherical). Hence
`σ_sw,m(n_l,ε_v) = −φ_m·Π(n_l)`, `∂σ_sw,m/∂n_l ≠ 0` — the broken pair, in code.

## 2. The parameter-free conjugate term [D]
One coupling free energy must generate both halves. With the **additive**
eigenstress entering as `σ = C:ε + σ_sw(n_l)·I_m`, the work-conjugate coupling
energy (volumetric sector) is

```
Ψ_cpl(ε_v, n_l) = σ_sw,m(n_l) · ε_v
```

Differentiating the SAME Ψ_cpl both ways:

```
∂Ψ_cpl/∂ε_v = σ_sw,m(n_l)              # the eigenstress — ALREADY implemented (LEFT)
∂Ψ_cpl/∂n_l = (∂σ_sw,m/∂n_l) · ε_v     # the MISSING conjugate (RIGHT)
```

so the restored micro potential gains, with `S₁ ≡ ∂σ_sw,m/∂n_l`,

```
        μ_lR_mech = (1/ρ_lR) · S₁ · ε_v                                   (★)
```

**Parameter-free:** `S₁` is the n_l-derivative of the eigenstress the
displacement residual already uses — no new constant. **Units:** `S₁` [Pa] × `ε_v`
[–] / `ρ_lR` [kg/m³] = J/kg ✓. Because `ε_v = σ'_m / K_drained` at the skeleton,
(★) **is** the mean-effective-stress dependence the spec asks for — strain and
stress views are one law, linked by the drained stiffness already in the model.

**Maxwell closes by construction:** `ρ_lR·∂μ_lR_mech/∂ε_v = S₁ = ∂σ_sw,m/∂n_l` ✓.

## 3. Consistent-tangent derivatives [D]
```
∂μ_lR_mech/∂ε_v  = S₁ / ρ_lR
∂μ_lR_mech/∂n_l  = (∂S₁/∂n_l)·ε_v / ρ_lR   = (∂²σ_sw,m/∂n_l²)·ε_v / ρ_lR
∂μ_lR_mech/∂ρ_lR = −μ_lR_mech / ρ_lR
```
Exchange `ρ̂ = −α_M·(μ_LR − μ_lR)` → NEW global-Jacobian block (exchange ↔ disp):
```
∂ρ̂/∂ε = +α_M·∂μ_lR/∂ε = (α_M·S₁/ρ_lR)·mᵀ          # m = Voigt identity, ε_v = mᵀε
```

## 4. The structural signature — SYMMETRY [D]
Off-diagonal global-tangent blocks:
```
[ u–u      u–n_l (= ∂σ_sw/∂n_l = S₁·m)   ]
[ n_l–u (= ∂μ_lR/∂ε = (S₁/ρ_lR)·mᵀ, NEW)   n_l–n_l ]
```
The two off-diagonals become a **transpose pair** (both ∝ S₁) → the tangent is
**symmetric**, which is exactly the discrete single-Ψ / integrability condition.
The current code has the upper block and a **zero** lower block → a
non-symmetric, non-conservative tangent.

**[PRED]** A strongly non-symmetric tangent is a candidate cause of the dd1600
`#104` SparseLU `compute()` singular wall (the dense case is the one that crosses
`σ_n → Π`, where the missing block bites). Restoring symmetry may regularise it.
**Verify by re-running dd1600 after (★); do not assert until then.**

## 5. Code sites
1. `PotentialExchange.h` — give the micro-potential helper the arguments `σ'_m`
   and `S₁`; add `μ_lR_mech` (★) and the three §3 derivatives. Gate per §6.
2. `ConstitutiveRelations/{ConstitutiveData,ConstitutiveSetting}` — route the
   current **mean effective stress** `σ'_m` and `S₁` into the exchange evaluation
   (it presently receives only micro state).
3. `RichardsMechanicsFEM-impl.h` — add the `(α_M·S₁/ρ_lR)·mᵀ` exchange↔disp
   Jacobian block; assert it is the transpose of the eigenstress block (§4).

## 6. DECISIONS — Vinay, before the algebra is coded
- **[DECISION] Sign of (★)** — choose so the steady state satisfies `Π(n_l)=σ_n`
  (spec §4.1). Drives whether load *expels* (correct) or imbibes micro water.
- **[DECISION] Effective-stress measure** — mean Bishop effective stress
  `σ'_m = tr(σ')/3` (recommended) vs net/total; confirm compression-positive sign.
- **[DECISION] Gate (spec §5)** — B1-sharp: term active only for `σ'_m ≥ Π` via a
  Macaulay `⟨σ'_m − Π⟩₊` (tangent kink at the gate; inherits the snap-drain) **vs**
  a smoothed ramp of width `δ_gate` (Newton-friendly; `δ_gate` to be **derived
  from problem scale**, value TBD). B2 (smear `g(Π)`) is the follow-on.
- **[DECISION] φ_m(ε_v) coupling** — keep `σ_sw,m`'s `ε_v`-dependence (through
  `φ_m`) inside `S₁`, or freeze `φ_m` at the GP? Needs the exact `φ_m` code
  definition; sets the §3 second derivative.

## 7. Verification — NOT fit-and-verify (anchors per spec §6)
1. **Maxwell-symmetry GP test** [anchor: derived identity]: FD-confirm
   `ρ_lR·∂μ_lR/∂ε_v == ∂σ_sw,m/∂n_l` at a sample state. Calibration-independent.
2. **Below-gate regression** [anchor: approved baseline]: `σ'_m < Π` cases (dd1400)
   reproduce current results bit-for-bit (term inactive there).
3. **Stress-driven exchange** [anchor: physical limit]: hold macro suction, ramp
   `σ'_m` past `Π` → `n_l` must drop (water expelled). Discriminating, blind to K.
4. **Jacobian symmetry** [anchor: integrability]: assembled off-diagonals are
   transposes to tolerance.
5. **[PRED]** dd1600 `#104` resolves — verify by re-run; relabel only after.
6. **K re-calibration** vs Dixon afterward (saturated `σ` shifts; spec estimates
   ~13% drop, **[PRED]**) — a separate step; never fit-and-verify with 1–5.

## 8. Status
**DESIGN — §6 decisions pending Vinay.** No code written; no parameter values
introduced; predicted items labelled. Successor to `MAXWELL_PAIR_RESTORATION.md`.
