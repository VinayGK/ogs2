# CurrentTotalSolidFraction — REJECTED, do not implement (energy-inconsistent)

**Verdict (2026-06-05, Vinay + verification):** the mode proposed in
`CURRENT_TOTAL_SOLID_FRACTION_IMPLEMENTATION.md` — make the disjoining cubic's
`nS = 1 − φ_total` (live) — is **energy-inconsistent** and must NOT be implemented.
The disjoining `Π` must stay a **micro** quantity; the macro porosity `φ_M` enters
the swelling stress **exactly once**, through the eigenstress contact-area factor.

## The principle (Vinay)
The effective stress acts on the aggregate contact area `1 − φ_M`. For the
`(σ', σ_sw)` pair to derive from one free energy (no `∮dW ≠ 0`), the swelling
stress must act **back through the same `1 − φ_M`** — i.e. scale as `(1−φ_M)¹`.

## Why the doc's mode violates it — three ways, one conclusion
Key identity: `1 − φ_total = (1 − φ_M)(1 − n_l)`, so `nS = 1−φ_total` hides a
factor `(1−φ_M)`, which the cubic cubes.

1. **Scale separation (decisive).** `Π` is the interlayer disjoining — the
   conjugate of the *micro* water; a micro potential cannot depend on a *macro*
   DOF. Exact test (`deliverables/energy_test_solid_fraction.py`, sympy):
   `∂Π/∂φ_M|_{n_l}` = **0** for `Reference` and `CurrentPorositySplit`,
   but `= 3(1−n_l)(1−φ0)²/(n_l³(1+ε_v)²) ≠ 0` for `CurrentTotalSolidFraction`.
   (The two micro-only modes giving exactly 0 is the passing sanity check.)
2. **Contact-area power.** `σ_sw = (1−φ_M)·n_l·Π` scales as `(1−φ_M)¹` for
   `CurrentPorositySplit` (correct) but `(1−φ_M)⁴` for `CurrentTotalSolidFraction`
   — three spurious contact-area powers ⇒ over a `φ_M` cycle the swelling stress
   no longer traces one free energy with `σ'` ⇒ **energy production**.
3. **Don't double-count / output-not-input.** The macro-porosity coupling is
   carried **once** by the eigenstress `(1−φ_M)`; putting it in `Π` too is a second
   copy. And `Pₛ(ρ_d)` is an *emergent OUTPUT* (`Π(n_l)` × contact area `(1−φ_M)` ×
   the micro↔macro equilibrium), not a constitutive input to hard-wire into `Π`.

## The correct path
- Keep `Π` micro-only. **`CurrentPorositySplit`** (`nS = 1 − n_l` = aggregate solid
  fraction = the live micro void ratio) is the physically clean LIVE mode: `Π(n_l)`
  only, `φ_M`-free. `Reference` (frozen) is an energy-safe special case.
- `Pₛ(ρ_d)` must **emerge** from the equilibrium. If it isn't, the lever is the
  uptake law / `K` / the exchange — **not** the swelling-stress fraction.
- Switching a model `Reference → CurrentPorositySplit` changes `nS` from `1−φ0`
  (≈0.63) to `1−n_l` (≈1.0 at hygroscopic IC) ⇒ ~4× in `Π` ⇒ **`K` must be
  re-calibrated** (a separate scientific step; not fit-and-verify).

## Verification provenance
`deliverables/energy_test_solid_fraction.py` (in the Task-13 models tree):
test (1) `∂Π/∂φ_M` [exact, sanity-passed] and (2) `(1−φ_M)`-power
[`CurrentPorositySplit`=1, `CurrentTotalSolidFraction`=4]. A finite-loop `∮` was
also coded but its energy-safe baseline did not isolate — **discarded**; the exact
differential closure (test 1) is the decisive, sanity-checked result.

Supersedes the implementation direction of
`CURRENT_TOTAL_SOLID_FRACTION_IMPLEMENTATION.md`.
