# Restore the forgotten Maxwell pair — NEXT IMPLEMENTATION GOAL

**Status (set 2026-06-02).** The DSM is currently admissible only on the
restricted domain `σ_n < Π` (claim **A** — written in `paper_DSM.tex` §2.4
`sec:admissibility`, and taught in `truesdell_noll_lecture.tex` pp. 77–82, the
Maxwell-pairs unit). This file is the spec for claim **B**: restore the forgotten
conjugate so the model is thermodynamically admissible **past** that gate
(high-stress / over-compaction / unload). The leading term is **parameter-free**
(its coefficient is fixed by the existing swelling closure); it re-touches the `K`
calibration. **STANDING RULE: read the current code before coding the algebra —
the swelling closure changed (see below).**

---

## 1. Apply the Maxwell test to this model

From a single free energy `Ψ(ε, n_l)` the Maxwell relation is

```
∂σ/∂n_l  =  ∂²Ψ/∂n_l ∂ε  =  ∂μ_lR/∂ε .
```

The DSM keeps the **left** and zeroes the **right**:

- **LEFT — present.** Swelling enters the macro momentum balance as the
  disjoining-pressure **eigenstress** `σ_sw = −φ_m · Π(n_l)` (commit
  `72f4f3a192`, "Π-eigenstress swelling closure"; full vdW + augmentation
  `p^disj`). So `∂σ/∂n_l ≠ 0` — **isotropic / mean-stress only** (the swelling
  eigenstress is spherical; the deviatoric sector is trivially `0 = 0`).
  *Note:* the legacy eigenstrain form `ε^sw = β_sw (n_l − n_l0)` was **retired**
  in that commit — swelling is now a stress, not a strain.
- **RIGHT — absent.** `μ_lR = μ_lR(n_l, ρ_lR)` only. Code fact
  (`PotentialExchange.h`): the solid-fraction argument is the **aggregate**
  fraction `1 − n_l` (a micro quantity, per the 2026-05-26 physics fix), **not**
  the total porosity `1 − φ` — so `μ_lR` carries **no strain and no stress**:
  `∂μ_lR/∂ε = 0`.

**⇒ broken Maxwell pair.** The closures are not the gradient of any single `Ψ`;
the `(σ, μ_lR)` response field is **non-conservative** — a loop in `(ε, n_l)` that
crosses the gate returns with net work (`∮ ≠ 0`). This is the *unlicensed
equipresence deletion* of the lecture's Maxwell-pairs unit: half a conjugate pair
kept, the other dropped.

## 2. Why it has held so far — the `σ_n < Π` gate

Exact **iff** `σ_n < Π`. Below the disjoining pressure the missing cross-term is
**physically zero**: the confining wall cannot squeeze out water the films already
out-push. Swelling-driven, monotonic loading (the swelling-pressure test, free /
low-stress wetting) lives entirely in this domain — stress is an *output*, no
closed loop excites the broken pair, the hole is invisible. The model is therefore
**exact** there (claim A).

It **breaks past Π** (over-compaction, high-stress constant-volume, unload legs):
no channel for stress to expel micro water. Code symptom: the **force balance
`Π(n_l) = σ_n` is never closed** — the exchange `ρ̂ = −α_M (ψ_M − ψ_m)` enforces
only the **chemical** balance `ψ_M = ψ_m`; under load the micro cannot de-swell
→ `φ_M → 0` → the porosity / saturation crash.

## 3. The goal

Give the micro potential a **mean-effective-stress dependence** — the Maxwell
partner of the swelling eigenstress — so that (i) `(σ, μ_lR)` derive from **one**
`Ψ`, and (ii) the exchange relaxes toward the **full** micro↔macro balance
(chemical **and** mechanical): `Π(n_l) = σ_n` as well as `ψ_M = ψ_m`.

```
μ_lR  =  ∂Ψ^ads/∂n_l  +  [mean-effective-stress term] .
```

The added term **is** the mechanical-expulsion (delamination) channel — the route
by which load past `Π` drives water out of the interlayer.

- **Parameter-free (leading term).** Derive both `σ_sw` and the `μ_lR` partner
  from **one** `Ψ`, so the partner's coefficient is fixed by the existing swelling
  closure — **no new constant**.
- **But re-verify `K`.** The change re-touches the saturated swelling-pressure
  equilibrium, so `potential_augmentation_prefactor` (`K`) must be re-calibrated per
  density vs Dixon afterwards (expect `K` to **drop** — the now-explicit term
  takes over part of `K`'s lumped-proxy job; cf. the full-`p^disj` shift, ~13 %).
- **Exact algebra against the current code.** Derive against the **Π-eigenstress**
  assembly, *not* the retired β_sw form. Read `PotentialExchange.h` (μ_lR) and the
  FEM-impl Π-eigenstress block before coding.

## 4. Implementation plan

1. **`μ_lR` gains the mean-effective-stress argument** (`PotentialExchange.h`).
   Add the term with the sign that makes the steady state satisfy `Π(n_l) = σ_n`.
2. **Analytic Jacobian — the Maxwell partner.** Add `∂μ_lR/∂σ_mean` (and the
   chain `∂ρ̂/∂ε` via `ψ_m`) to the consistent tangent. *The symmetry is the
   point:* the eigenstress half (`∂σ/∂n_l`, already in) and the expulsion half
   (`∂μ_lR/∂ε`, new) must **both** sit in the Jacobian. Reuse the GP
   single-integration-point test harness (the eigenstress+exchange Jacobian test).
3. **Close the force balance in the exchange.** With `σ_n` in `ψ_m`,
   `ρ̂ = −α_M (ψ_M − ψ_m)` now relaxes the mechanical balance too. Verify the
   steady state satisfies **both** balances.
4. **Maxwell-symmetry unit test** (the integrability check the current code
   fails): numerically confirm `∂σ/∂n_l = ∂μ_lR/∂ε` at a Gauss point.

## 5. The sharp-Π tax — two tiers

With the **current sharp `Π`** (a δ-function binding spectrum, `g = δ(Π − Π_0)` on
the stress axis), restoring the term **relocates** the crash: past `Π` the micro
**snap-drains** (all the water releases at once) instead of de-swelling gradually.
So:

- **B1 — restore the term, sharp `Π`.** Closes the Maxwell pair; admissible in
  principle past the gate; inherits the snap-drain. Parameter-free. *(milestone)*
- **B2 — + smear `Π`.** Resolve the binding spectrum `g(Π)` (sub-pools `n_l(Π)`)
  so release is gradual from the first load increment — matches the smooth macro
  collapse the MIP shows at point C, and the smoothness BExM gets free from its
  fitted interaction functions. Heavier; forces wetting/drying hysteresis head-on.
  *(follow-on)*

B1 is the structural fix (the pair); B2 is the smoothness that makes B1 robust.

## 6. Verification (NOT fit-and-verify-in-the-same-test)

- The new term is **not identifiable** from the swelling-pressure test that
  anchored `K` — the calibration is **blind** to it (`σ_n` pinned at `Π` there).
  Do not calibrate-and-assert in one step.
- Identify / verify from: the **swelling-strain-vs-water-content** curve;
  **high-stress (> Π) compression**; the **unload legs**. Empirical anchor:
  Ferrari, Bosch, Baryla & Rosone, *Acta Geotechnica* 2022 (granular MX80
  two-path oedometer; intra-aggregate pore volume grows while inter collapses;
  void ratio path-dependent, compression converges to a unique VCL).
- **Falsification:** unload from point C. Reversible ⇒ a static `g(Π)` suffices
  (B1 enough). Hysteretic ⇒ path-dependent spectrum (micro-cavitation) ⇒ B2
  mandatory.
- After: **re-verify `K`** per density vs Dixon (saturated swelling pressure
  shifts).

## 7. Who kept the term — BExM

BExM's micro effective stress `p̂ = p + s` carries **both** the mechanical `p` and
the suction `s`; its interaction functions `f_c, f_s` **are** the BExM analogue of
our missing cross-term — but **fitted** (two calibrated curves) and smooth from
the start (no snap-drain). The DSM goal is to **derive** (one `Ψ`, parameter-free)
what BExM **fits**, then pay the sharp-`Π` tax (B2) that BExM sidestepped by being
smooth. Cite: Gens & Alonso 1992 (BExM); Coleman & Gurtin 1967 (internal
variable); Collins & Houlsby 1997 (MCC hyperplasticity).

## 8. References

- **Code** (`dsm_native_hierarchical`): `PotentialExchange.h` (μ_lR; exchange ρ̂);
  `RichardsMechanicsFEM-impl.h` (Π-eigenstress swelling-stress block; the
  hierarchical split `φ_M = (φ − n_l)/(1 − n_l)`).
- **Theory:** `truesdell_noll_lecture.tex` pp. 77–82 (Maxwell pairs + "The broken
  pair in the DSM"); `paper_DSM.tex` §2.4 `sec:admissibility` (claim A + the B
  extension); claude memory `dsm-equipresence-audit` (the per-term audit, the
  `g(Π)` spectrum unification, the `σ_n < Π` gate).
