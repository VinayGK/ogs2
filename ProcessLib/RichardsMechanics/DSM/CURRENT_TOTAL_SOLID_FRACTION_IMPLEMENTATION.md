# CurrentTotalSolidFraction — make the swelling pressure track the LIVE dry density (Pₛ(ρ_d))

> **❌ REJECTED — DO NOT IMPLEMENT (2026-06-05).** Energy-inconsistent: it puts the
> macro porosity `φ_M` inside the micro disjoining `Π` (`1−φ_total = (1−φ_M)(1−n_l)`),
> making `σ_sw ∝ (1−φ_M)⁴` instead of the contact-area-conjugate `(1−φ_M)¹` →
> energy production over a `φ_M` cycle (verified: `∂Π/∂φ_M ≠ 0` only for this mode).
> Keep `Π` micro-only (`CurrentPorositySplit`); `Pₛ(ρ_d)` must emerge from the
> equilibrium. Full argument + verification: `CURRENT_TOTAL_SOLID_FRACTION_REJECTION.md`.

**Branch:** `dsm_native_pdisj_maxwell`. **Work in a FRESH worktree off it — do not disturb the active one.**
**Authored:** 2026-06-05 (Vinay decision + agent diagnosis). **Status: NOT yet implemented — this is the spec.**
**Why it matters:** this is the *prerequisite* that makes the whole Maxwell-conjugate apparatus (gate +
magnitude + the (1−φ_M) REV-referencing) actually do anything. Verified inert without it (§1).

---

## 0. Guardrails (announce on firing)
This repo's `CLAUDE.md` (§0–§12) governs. Before acting, re-read it. Specifically:
- **No material parameter / BC magnitude / expected value without a cited source** (paper+locator,
  Vinay's words, or an in-file derivation) → else STOP and ASK.
- **No fit-and-verify in one test** (§2): do NOT re-fit K and then assert on a Pₛ derived from K.
- **Supplement, never replace** (§3): the existing `Reference` and `CurrentPorositySplit` modes stay
  byte-identical; this is an opt-in third mode. Default stays `Reference`.
- **predicted ≠ verified**: label every expected result predicted until a re-run confirms it.
- Announce in-chat the moment any guardrail fires.

---

## 1. Diagnosis — why (read first)
The disjoining/swelling potential is **cubic in the micro solid volume fraction n_S**
(`ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h:180–181`,
`computeVanDerWaalsMicroPotential`):
```
mu_lR = sign * [A·Sₐ³/(6π)] * (nS·nS·nS) * (rho_SR³) / (n_l³ · rho_lR)      and   Π = −rho_lR · mu_lR
```
Since **n_S·ρ_SR = (1−φ)·ρ_s = ρ_d** (ρ_SR≈const, REV-referenced), the formula already contains
**Π ∝ ρ_d³** — i.e. a Pₛ(ρ_d) is structurally present.

**But n_S is not live.** `computeActiveMicroSolidVolumeFraction`
(`RichardsMechanicsFEM-impl.h:396–409`) returns either:
- `Reference` (the **default**, `PotentialExchangeParameters.h:115–116`): a **frozen PRJ constant**
  `micro_solid_volume_fraction_reference` → static in ρ_d; or
- `CurrentPorositySplit`: `1 − n_l` → tracks the **interlayer water** (micro state), **not** the macro
  densification.

Neither equals **1 − φ_total = ρ_d/ρ_s**. So the ρ_d³ in line 180 is disconnected from the current density.

**Consequence (verified 2026-06-05, `/tmp/probe` staged LE restart from B, e=2.36, true K=44200):**
under compression ρ_d rises but n_S (frozen) does not → Π **collapses** (σ_sw: 0.87→~0 MPa) → the
Maxwell term (∝ Π) dies. The probe gate opened decisively (p′→15.6 MPa ≫ φ_m·Π) yet the (1−φ_M)
magnitude made **zero** difference (revref ≡ old-on, max|Δ|=2×10⁻¹⁶). On EPFL the term is also inert
(p′ hugs the gate). **The term is inert in every reachable regime until Π is made density-live. This change
fixes that.**

---

## 2. The fix (design)
Add a **third** `MicroSolidVolumeFractionMode::CurrentTotalSolidFraction` that returns
**n_S = 1 − φ_total (live)**. Then Π ∝ (1−φ_total)³ρ_SR³ = ρ_d³ tracks the live density and the existing
cubic becomes a genuine **Pₛ(ρ_d)** — **parameter-free** (no new constant; reuses the existing K), exactly
like the equipresence restoration.

**No new stored state variable, no new data read-in.** The live total porosity is already in scope:
`PotentialExchangeLocalSolveContext::phi` (`RichardsMechanicsFEM-impl.h:279`), with a kinematic fallback
(`(φ_prev + Δε_v)/(1+Δε_v)`) already coded in `boundedMicroWaterContentCeiling`
(`RichardsMechanicsFEM-impl.h:291–317`). The ONLY new PRJ key is the mode-selector string.

---

## 3. Exact code changes

### 3.1 Enum + toString — `PotentialExchangeParameters.h`
- Add `CurrentTotalSolidFraction` to `enum class MicroSolidVolumeFractionMode` (line ~29).
- Add to `toString(MicroSolidVolumeFractionMode)` (line ~81): `case CurrentTotalSolidFraction: return "current_total_solid_fraction";`

### 3.2 Parser — `CreateRichardsMechanicsProcess.cpp:90–99`
In `parseMicroSolidVolumeFractionMode`, add: `if (s == "current_total_solid_fraction") return MicroSolidVolumeFractionMode::CurrentTotalSolidFraction;`
Keep `reference` and `current_porosity_split` unchanged; keep `Reference` the default.

### 3.3 Helper — factor out the live total porosity (DRY)
Extract the total-porosity logic currently inside `boundedMicroWaterContentCeiling` (lines 291–317) into:
```cpp
inline double liveTotalPorosity(PotentialExchangeLocalSolveContext const& c)
{
    constexpr double cap = 1.0 - 1e-12;
    if (std::isfinite(c.phi)) return std::clamp(std::max(0.0, c.phi), 0.0, cap);
    double const phi_prev = std::clamp(std::max(0.0,c.phi_M_prev)+std::max(0.0,c.phi_m_prev),0.0,cap);
    double const d = 1.0 + (c.volumetric_strain - c.volumetric_strain_prev);
    if (std::isfinite(d) && std::abs(d) > 1e-12)
        return std::clamp((phi_prev + (c.volumetric_strain - c.volumetric_strain_prev))/d, 0.0, cap);
    return phi_prev;
}
```
and have `boundedMicroWaterContentCeiling` call it (no behaviour change — verify byte-identical).

### 3.4 `computeActiveMicroSolidVolumeFraction` — `RichardsMechanicsFEM-impl.h:396–409`
Add, BEFORE the existing returns (leave `Reference`/`CurrentPorositySplit` untouched):
```cpp
if (mode == MicroSolidVolumeFractionMode::CurrentTotalSolidFraction)
    return std::max(1e-16, 1.0 - liveTotalPorosity(local_context));   // 1 − φ_total = ρ_d/ρ_s, live
```

### 3.5 `computePreviousMicroSolidVolumeFraction` — `RichardsMechanicsFEM-impl.h:411+` — **MIRROR (do not skip)**
Add the same case using the **previous** total porosity `phi_M_prev + phi_m_prev`:
```cpp
if (mode == MicroSolidVolumeFractionMode::CurrentTotalSolidFraction)
    return std::max(1e-16, 1.0 - std::clamp(phi_M_prev + phi_m_prev, 0.0, 1.0 - 1e-12));
```
**If the previous-state function is not mirrored, the eigenstress increment
Δσ_sw = n_S(n_l_prev·Π_prev − n_l·Π) mixes a live "current" with a frozen "previous" → spurious
energy per step.**

### 3.6 Consistent tangent (the chain) — required for convergence
The cubic is steep ("steep restoring spring"). Add the missing chain to the analytic Jacobian:
```
∂μ_lR/∂ε_v = (∂μ_lR/∂n_S)·(∂n_S/∂φ)·(∂φ/∂ε_v)
```
`∂μ_lR/∂n_S` already exists (`dmu_lR_dnS`, PotentialExchange.h:185); `∂n_S/∂φ = −1`;
`∂φ/∂ε_v` from `PorosityFromMassBalance` (≈ (1−φ)/(1+ε_v) for small strain — derive in-file, cite).
This overlaps the still-missing exchange↔displacement Jacobian block (see
`MAXWELL_CONJUGATE_IMPLEMENTATION.md` §3/§4). Without it expect Newton stalls on the steep parts;
a finite-difference Jacobian is an acceptable interim **with a logged TODO**, not a silent omission.

---

## 4. Initialization (no explicit step needed — but assert it)
At t=0 the live φ = `phi0` (`initial_porosity`), and the PRJs set
`micro_solid_volume_fraction_reference = 1 − phi0`. So 1−φ(0) = reference **by construction** — the live
n_S starts at the reference value, then φ takes over. Do NOT add an init parameter. **Do add an assert/warn**
that `|(1 − phi0) − micro_solid_volume_fraction_reference| < 1e-6` when the mode is selected (a mismatch
would cause a t=0 Pₛ jump).

---

## 5. Gotchas (do not skip)
1. **Two functions** (active §3.4 AND previous §3.5) — both live, consistently.
2. **n_S is re-read** as the *total* packing (→ macroscopic Pₛ(ρ_d), the goal). Keep it **distinct** from
   the Maxwell-conjugate term's own `n_S = 1 − φ_M` (aggregate fraction, `PotentialExchange.h:247`) — that
   one is already live via the porosity split; do not cross the two.
3. **Default unchanged.** `Reference` stays the default and must be byte-identical (existing calibrated runs).

---

## 6. Verification gates (MANDATORY — §3 working rules)
1. **Build clean**; `./bin/testrunner --gtest_filter='*DSMMicroMacro*'` → 14/14 pass. The `Reference` and
   `CurrentPorositySplit` paths MUST be byte-identical to pre-change (the new mode is opt-in).
2. **Calibration CHECK, not a fit (§2).** Single element, no displacement BC, vary the imposed dry density
   (via `initial_porosity`); with the new mode + the **existing** K (do NOT re-fit), read the resulting Pₛ
   and compare its *shape* to the approved Pₛ(ρ_d) anchors:
   - Dixon (2023) MX-80, EMDD≡ρ_d (WG agreement 2026-05-27): **4.92 / 14.16 / 40.86 MPa** at
     ρ_d = 1.40 / 1.60 / 1.80 g/cm³ (per AGENTS.md 2026-06-01); continuous form Villar FEBEX
     Pₛ = exp(6.77·ρ_d − 9.07) [MPa].
   - **If n_S³ + existing K reproduces the curve within tolerance → the cubic is the right carrier (verified).**
   - **If not → report the gap. Do NOT re-fit K to force the match (that is fit-and-verify, §2).**
3. **Re-run the high-load probe** (`/tmp/probe/beacon_probe_compress.prj`, or rebuild from
   `beacon_t33_path1_P1-3_LE` → restart). Expected with the live mode: Π **rises** under compression
   (not collapses), and **revref vs old-on Δ becomes non-zero** — the (1−φ_M) magnitude finally bites.
   Label predicted→verified once run.
4. **Zero rejected steps** on the canonical Model-I LE reruns.
5. **Docs:** update `DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md` (new step section) and `AGENTS.md`
   known-limitations (mark the static-swelling / dry-density gap → resolved by this mode), and add a unit
   test asserting Pₛ rises monotonically with imposed ρ_d under the new mode.

---

## 7. Acceptance
- Default (`Reference`) and `CurrentPorositySplit` byte-identical; new mode opt-in via one PRJ string.
- With the new mode, the cubic reproduces the Dixon/Villar Pₛ(ρ_d) anchors **without re-fitting K** (or the
  gap is reported, not fitted away).
- The high-load probe shows Π rising under load and the Maxwell term (with the (1−φ_M) magnitude) finally
  producing a non-zero effect.
- All energetic claims **predicted** until the re-runs confirm them.

---

## 8. Provenance / locators (all @ `dsm_native_pdisj_maxwell`)
- Swelling cubic: `PotentialExchange.h:180–181` (Π = −ρ_lR·μ_lR; μ_lR ∝ n_S³ρ_SR³/(n_l³ρ_lR)).
- n_S modes: `computeActiveMicroSolidVolumeFraction` `RichardsMechanicsFEM-impl.h:396–409`;
  previous `:411+`; default `Reference` `PotentialExchangeParameters.h:115–116`.
- Live φ in context: `PotentialExchangeLocalSolveContext::phi` `RichardsMechanicsFEM-impl.h:279`;
  kinematic fallback `boundedMicroWaterContentCeiling` `:291–317`.
- Parser: `CreateRichardsMechanicsProcess.cpp:90–99`.
- Pₛ(ρ_d) anchors: Dixon 2023 MX-80 (EMDD≡ρ_d); Villar FEBEX exp(6.77ρ_d−9.07); AGENTS.md 2026-06-01.
- Inert-term diagnosis: probe 2026-06-05 (`/tmp/probe`, LE staged restart from B, e=2.36); MCC apex-blocked
  at the soft restart (status −1) — use LE for the swell, see `MAXWELL_CONJUGATE_IMPLEMENTATION.md` for the
  Maxwell-term context this unblocks.
