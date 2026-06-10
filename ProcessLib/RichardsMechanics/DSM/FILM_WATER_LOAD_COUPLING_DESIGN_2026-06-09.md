# Design — energetic film-water/load coupling for the swelling closure

**Status:** PROPOSED — not implemented, not run. Formulation proposed by Vinay
(discussion 2026-06-09). Companion discussion deck:
`~/tex/implementation_and_theory/implementation/film_pressure_imbibition_problem.{tex,pdf}`.
Diagnosis it builds on: `~/ogs-models/EBS/Task13/2026_06_05_maxwell_task13/DIAGNOSIS_maxwell_nonconvergence.md`
and this session's `DIAGNOSIS_2026-06-09_full_maxwell` evidence (MCC + LE step-1 blow-up,
α_M-independent, Π≈2e11 Pa at Task-13 IC).

This is a design/derivation note (no numerical literals decided here). All material
values, the exact cross-term, and every test's expected value remain Vinay's call.

---

## 0. Problem (recap)

The always-on integrable Maxwell film-pressure closure routes the confining load as a
*static subtraction*:

    sigma_sw = -phi_m * p_film,   p_film = Pi(n_l) - b*p_conf,   phi_m = n_S*n_l,
    Pi(n_l)  = -rho_lR * mu_lR^vdW(n_l),   dPi/dn_l < 0.

`Pi` is held at the current `n_l` within a step. When `b*p_conf` outpaces `Pi`, `p_film<0`
and the conjugate micro potential drives the exchange **backwards** — `n_l` increases
(imbibition) under load. Ungated, this misbehaves in two observed facets (neither run
literally crosses `sigma'_eff>Pi`; loads are << Pi):

- **Task 13 stg1 — Pi-magnitude facet.** Dry/dense IC, `Pi≈2e11 Pa` → exchange/local
  micro solve runs away in step 1 (`Pi:1e11→1e44`, `p→1e13`); MCC *and* LE diverge,
  dt-independent to 1e-10 s; α_M-independent (1e-13…1e-26 all fail).
- **MS33 III — eps_v facet.** `mu_lR^mech ∝ (Pi + n_l*Pi')*eps_v` grows with swelling
  strain; survives slow wetting, diverges at 50× permeability.
- **MS33 VII — control (NOT a problem).** Loads to 5 MPa but `sigma'_eff ≤ 1.7 MPa`
  vs `Pi≈48 MPa` (`sigma'_eff/Pi≈0.03`); converges. (Its over-swell e≈1.50 is the
  density-decoupled drive — a separate item.)

## 1. Proposed formulation (Vinay, 2026-06-09)

Route the confining **work into the film water content**, not the eigenstress alone:

- `p_conf` does work that **expels** film water (`n_l ↓`).
- Since `Pi = Pi(n_l)` with `dPi/dn_l < 0`, expulsion **raises** `Pi` until it meets the
  load: `Pi(n_l) → b*p_conf`, hence `p_film → 0+`, `sigma_sw → 0` (clean consolidation).
- The `sigma'_eff > Pi` branch is then **never reached** — `n_l` (and `Pi`) adjust before
  `p_film` can flip. The self-limiter is `dPi/dn_l < 0` (the film fights back as it thins).

**Free-energy structure.** Replace the static subtraction by a genuinely cross-coupled
`Psi(eps_v, n_l)` so the micro responds to a **net (effective) stress**, with the Maxwell
pair from one `Psi`:

    sigma = dPsi/d(eps),  mu_lR = dPsi/d(m_l),  d(sigma_sw)/d(n_l) = d(mu_lR)/d(eps_v).

This is the double-structure **micro-effective-stress** idea — micro strain/water driven
by `p_hat = p + s` — established in:
- Gens & Alonso (1992), "A framework for the behaviour of unsaturated expansive clays",
  *Can. Geotech. J.* 29(6).
- Sánchez, Gens, Guimarães & Olivella (2005), "A double structure generalized plasticity
  model for expansive materials", *Int. J. Numer. Anal. Methods Geomech.* 29.

The present `mu_lR^mech ∝ (Pi + n_l*Pi')*eps_v` is integrable but encodes the load→water
response with the **wrong sign** (imbibition).

**OPEN (Vinay's call) — exact cross-term, two candidate routes:**
- (a) **Load in the disjoining argument:** the film's equilibrium `n_l` shifts with a
  micro-effective stress `p_hat` (the film "sees" `p_conf`), so `Pi` re-equilibrates
  through `n_l(p_hat)`.
- (b) **Explicit `Psi_couple(eps_v, n_l)`** term whose derivatives give both the
  load→expulsion and the Maxwell-symmetric stress.

## 2. Implementation plan (code sites)

Files: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`,
`ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`.

1. **Local micro solve** (`scalar_micro_macro_mass_storage_mode`, ~FEM-impl.h L717): let the
   micro equilibrium `n_l` depend on `p_conf` (confining enters the micro water balance) so
   loading expels water. This is the core change — the load→`n_l` channel.
2. **Disjoining `Pi(n_l)`**: unchanged in form; becomes load-responsive *via* `n_l`.
3. **Swelling-stress block** (~FEM-impl.h L1843, `sigma_sw = -phi_m*p_film`): with `n_l`
   load-responsive, `p_film → 0` at equilibrium naturally; remove reliance on the negative
   branch.
4. **Exchange residual + Jacobian**: add the `d(micro eq)/d(p_conf)` coupling and its
   tangent (new `n_l ↔ p_conf` cross-term), kept Maxwell-symmetric with the existing p–u
   swelling-stress tangent (FEM-impl.h ~L4200). FD-Jacobian fallback available
   (`fd_jacobian_for_exchange`) for bring-up.
5. **PRJ flag** (proposed `film_water_load_coupling`, default = Vinay's call): when OFF,
   reproduce the current closure **bit-for-bit** (regression safety). Per §0.1 the default
   is a formulation decision — do not pick silently.
6. Discipline: `+=` on accumulators (§4.1); unit-annotate any new potential/stress term
   (§4.2); update this doc + the affected `AGENTS.md` (§8).

## 3. Tests — STRUCTURE ONLY (expected values are `TODO(Vinay)`, §3)

> Five tests proposed (the §3 batch limit). Each states a physics anchor from the closed
> list; expected values are left for Vinay or proposed with a cited source. New tests
> SUPPLEMENT the existing suite — they do not replace it.

**T1 — single-element consolidation under load past Pi**
- Physics anchor: (a) analytical limit — film disjoining balance `Pi = b*p_conf`.
- Input: saturated micro-macro GP; ramp `p_conf` from `<Pi` to `>Pi`, isothermal.
- Expected: `TODO(Vinay)` — assert `n_l` **decreases monotonically** (expulsion, `Δn_l<0`
  throughout) and `Pi → b*p_conf`, `sigma_sw → 0+`; NO imbibition. (The sign assertion
  needs no literal; the `Pi=b*p_conf` endpoint is the analytical target.)
- Catches: the imbibition-under-load pathology (current code gives `Δn_l>0`, `p_film<0`).
- Overlap: none (new GP test).

**T2 — Maxwell symmetry of the coupled Psi**
- Physics anchor: (d) symmetry.
- Input: GP; finite-difference `d(sigma_sw)/d(n_l)` vs `d(mu_lR)/d(eps_v)`.
- Expected: equality to a tolerance derived from problem scale (§1.2; no material literal).
- Catches: a non-integrable / energy-non-conserving coupling.
- Overlap: none.

**T3 — `sigma'_eff > Pi` unreachability**
- Physics anchor: (a) analytical limit.
- Input: the state where the *current* closure imbibes (`p_conf` driven toward `Pi`).
- Expected: `TODO(Vinay)` — `p_film` stays `≥0` (n_l adjusts first); `sigma_sw` never tensile.
- Catches: residual sign-flip / imbibition branch.
- Overlap: contrasts the current closure (documented behaviour); not a replacement.

**T4 — EBS Task 13 Stage 1a/1b convergence + experiment**
- Physics anchor: (e) published benchmark — EBS Task 13 (Task-13-Description-Stages-1-2-3;
  Villar et al., CIEMAT report; digitized `expt_stg1_05.csv` / `expt_stg1_6.csv`).
- Input: existing Task-13 stg1 PRJs (0.5 / 6 MPa) on the coupled closure.
- Expected: `TODO(Vinay)` — converges (target: to ~466 d, as the non-conjugate binary does);
  swelling pressure / void ratio vs the digitized experiment.
- **§2 FLAG (fit-vs-verify):** `K` is calibrated **upstream** (Dixon / Model I) — assert on
  the Task-13 **transient/consolidation**, NOT on a quantity used to fit `K`. Keep the
  streams separate.
- Catches: the step-1 blow-up.
- Overlap: separate benchmark family (not the MS33 ctest suite).

**T5 — MS33 submission invariance (regression)**
- Physics anchor: (f) regression baseline previously approved by user — the committed
  `mc_20260608` 20× results (III p≈9.63 MPa; IV p≈3.66 MPa, e≈1.09; VII e≈1.50;
  commit `2c3b3bb0d6`).
- Input: `ANCHORS_MS33_Model{I,III,IV,VII}` on the coupled closure (flag ON), 20× perm.
- Expected: `TODO(Vinay)` — match the committed results to tolerance **below the cross-over**
  (the fix must not perturb the safe regime; VII already lives there).
- **§3 OVERLAP FLAG:** SUPPLEMENTS the existing ANCHORS suite (regression guard); does NOT
  replace it. The user's suite is authoritative.
- Catches: an unintended change to the calibrated / submission regime.

## 4. Open decisions for Vinay

- Exact `Psi(eps_v, n_l)` cross-term (route a vs b, §1).
- Does the self-limiting (`dPi/dn_l<0`) remove the need for any gate (stay "always-on"),
  or is a soft cap still wanted as a backstop?
- New Jacobian cross-term + conditioning; whether FD-Jacobian suffices for production.
- PRJ flag name + default state.
- Per-material handling (clay vs pellet) of the coupling — **ask, don't assume**
  (per-material parametrization rule); record the choice once made.

## 5. Predicted impact — NOT verified (§5)

- **Predicted (formulation only, nothing re-run):** Task 13 consolidates under load (no
  step-1 blow-up); MS33 III stable at 50× (the eps_v facet self-limits); MS33 I/IV/VII
  invariant below the cross-over. To be **upgraded to verified only** after the affected
  cases are re-run and the predictions confirmed.
