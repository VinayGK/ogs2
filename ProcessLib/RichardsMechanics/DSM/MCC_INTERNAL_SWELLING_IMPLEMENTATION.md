# MCC internal swelling — porting the DSM eigenstress INSIDE the ModCamClay behaviours — IMPLEMENTATION DESIGN

**Branch:** `dsm_native_h_of_eps` (HEAD `23a723cc3c` at design time; same commit
as the Task-13 diagnostic binary).
**Builds on:** `MAXWELL_CONJUGATE_IMPLEMENTATION.md`,
`STRAINED_FILM_IMPLEMENTATION.md` (template + the eigenstress closure),
`/Users/vinaykumar/git/ogs/task13_absp_2026-06-11/TASK13_MCC_BLOCKAGE.md`
(the motivating record).
**Working name of the new behaviour family:** `ModCamClay_semiExpl_*_intSwell`
(naming → DECISION D4).

Tags: **[D]** derived/established · **[DECISION]** Vinay's call before coding ·
**[PRED]** predicted, not verified.

**STATUS: DESIGN ONLY. Nothing in this file is implemented.** (§11.)

---

## 0. Scope and motivation

Established facts from the Task-13 record (all MEASURED there unless tagged;
`TASK13_MCC_BLOCKAGE.md`, binary `dsm_native_h_of_eps @ 23a723cc3c`):

- Stage-1 dry-density residual: model DD ≈ **1.459** vs experiment **1.4139**.
- The residual is **dry-side plasticity that cannot engage**: at the
  force-balance state p′ = 1.6 MPa, q = 3.2 MPa, q/p′ ≈ 2 > M = 1
  (supercritical side, associated MCC flow ⇒ plastic dilation ⇒ lower DD);
  the state sits only ~25 % inside yield (q_yield ≈ 4.1 MPa at this p′);
  ε_p ≡ 0, p_c frozen at 12 MPa — plasticity DORMANT by construction
  (p_c was *raised* above the 11.46 MPa swelling stress to avoid the crash).
- Both engagement paths hit the **MCC semiExpl integrator wall** when the DSM
  swelling eigenstress is applied EXTERNALLY (as eigenstrain):
  - `constE` + lowered p_c (6 / 4 MPa): MFront status −1 cascades from yield
    onset (t ≈ 10–68 s) without recovery and without ever engaging plasticity;
  - `absP`: status −1 at step #1 at the soft dry IC (σ0 = 0.15 MPa, OCR ≈ 80
    — K ∝ p′/κ ≈ 4.8 MPa, too soft to absorb the swelling-stress build-up),
    independent of p_c.
- It is NOT K (live-K Villar re-fit gives a WORSE 1.491) and NOT (only) the
  irreversible/Route-R channel (Vinay, 2026-06-11).

**The port.** Move the *application site* of the swelling stress from the
process (eigenstrain added to the mechanical strain handed to MFront) into the
behaviour itself (an isotropic stress contribution supplied per increment as
an external state variable). This

(a) removes the external stress shock: Δσ_sw no longer arrives as a strain
    increment C_el⁻¹·Δσ_sw inside `deto`, so the elastic predictor no longer
    jumps far across the yield surface each step;
(b) keeps the stress variable that the absP hypoelastic law works with away
    from the p′→0 singularity (no eigenstrain-induced volumetric strain shock
    through a soft K; optionally, under D5(β), the swelling compression
    actively stiffens K) [PRED until run];
(c) keeps the **entire DSM micro engine unchanged** — Pi-path, Maxwell pair,
    live K(ρ_d), film modes all stay the driver supplier; only where σ_sw is
    *applied* moves.

Out of scope: any change to the DSM eigenstress closure itself; BExM/BBM LC
plasticity (route 3 of the blockage record); the implicit-MCC rewrite
(route 2) — this design is the minimal-surgery variant of route 2's intent.

---

## 1. The code as read [D]

All paths relative to the worktree root
`/Users/vinaykumar/git/ogs-worktrees/dsm_native_h_of_eps_wt/`.

### 1.1 The three ModCamClay behaviours
(`MaterialLib/SolidModels/MFront/ModCamClay_semiExpl{,_absP,_constE}.mfront`)

All three: `@DSL Implicit` (L4 in each), `@Theta 1.0`, `@Epsilon 1e-14`,
`@ModellingHypotheses{".+"}`. Integration variables: `eel` (implicit DSL
default), `lp` (`@StateVariable real lp`), plus `rpc`
(`@IntegrationVariable strain rpc`, reduced p_c) in `semiExpl` (L70) and
`constE` (L75); **`absP` has NO p_c integration variable** — see below.
Auxiliary state: `pc`, `epl_V`, `v` (volume ratio) in all three.

**Algorithm state — SURPRISE, design-relevant:**

| file | `@Algorithm` line | status |
|---|---|---|
| `ModCamClay_semiExpl.mfront` | L36: `NewtonRaphson; //_NumericalJacobian;_LevenbergMarquardt` | NR active, alternatives commented |
| `ModCamClay_semiExpl_absP.mfront` | L35: `LevenbergMarquardt; // _NewtonRaphson;_NumericalJacobian (Step-2 diagnostic 2026-05-28)` | **LM already active** |
| `ModCamClay_semiExpl_constE.mfront` | L34: `LevenbergMarquardt; // 2026-06-06: was NewtonRaphson; LM more robust …` | **LM already active** |

Both LM switches (`372da0aafe` for absP, `2a2409c043` for constE) are
ancestors of the Task-13 binary commit `23a723cc3c` (verified via
`git merge-base --is-ancestor`). **The Task-13 status −1 walls were hit WITH
LevenbergMarquardt already on.** LM alone does not cure the blockage — this
retires the "flip to LM" half of the original D2 question (§7).

**Residual systems (quoted from the files):**

- `semiExpl` (`@Integrator`, L130–188):
  ```
  feel = deel + depl - deto;                                   (L166)
  flp  = f / fchar;            f = q*q + M2*p*(p - pc_new)     (L167, L152)
  frpc = drpc + deplV * the * (rpc_new - rpc_min);             (L168)
  ```
  with `p = -trace(sig)/3` (L147), `the = v0/(la-ka)` (L134), elastic range
  short-circuit `feel -= deto; return true;` (L137–141).
- `constE` (`@Integrator`, L118–181): `@Brick StandardElasticity` (L36)
  supplies the elastic part of `feel`, so the plastic correction is
  `feel += depl;` (L154); `flp = f/fchar` (L155);
  `frpc = drpc + deplV*the*(rpc_new - rpc_min)` (L156);
  `p = -trace(sig)/3 + pamb` (L134; `pamb` default 1e3 Pa, L40–42) — the
  ambient-pressure shift that guarantees an initial elastic range (L110–111).
  `eel` is (re)initialized from `sig` at L100:
  `eel = (1+nu)/young*sig - nu/young*trace(sig)*Stensor::Id();`.
- `absP` (`@Integrator`, L162–234): residuals are **`feel` and `flp` only**
  ```
  feel = deel + depl - deto;                                   (L201)
  flp  = f / fchar;            f = q*q + M2*p*(p - pc_new)     (L202, L186)
  ```
  p_c is updated **explicitly inside the residual evaluation** (no `fpc`):
  `pc_new = (pc - pc_min)*exp(-the*deplV) + pc_min;` (L184) with
  `deplV = trace(depl)` recomputed from the flow direction (L195–196), and
  again, definitively, in `@UpdateAuxiliaryStateVariables`
  (`pc = (pc - pc_min)*exp(-v0/(la-ka)*deplV) + pc_min;`, L243).
  NOTE: the header comment block (L26–28) still *names* an `fpc` residual —
  it documents the family, not this file; the code has none.

**Hypoelastic K — absP vs constE vs semiExpl:**

- `absP`: absolute/integral form. `@InitLocalVariables`:
  `p0 = -trace(sig)/3` (L137, stored L144), `K = v0/ka * p` (L138),
  `young = 3K(1-2nu)` (L140), `pc_min = 0.5e-8*pc_char` (L141).
  `@ComputeStress` (L124–128) calls the file-local helper (L92–119):
  ```
  p = p0 * exp(-v0_ka * deelV);   (L111)
  K = v0_ka * p;                  (L112)   // v0_ka = v0/ka
  G = alpha * K;                  (L113)
  ```
  **The singularity:** K ∝ p; at the Task-13 dry IC, p0 = σ0 = 0.15 MPa makes
  K (≈ 4.8 MPa, measured in the blockage record) far too soft, and the
  external eigenstrain shock C_el⁻¹·Δσ_sw inside `deto`/`deel` drives
  `exp(-v0_ka*deelV)` to extreme values — p collapses toward 0 (where K → 0,
  L112) or explodes. There is **no floor on p anywhere in this file**.
- `constE`: constant `young`, `nu` via `@Brick StandardElasticity` (L36);
  no pressure dependence at all.
- `semiExpl`: incremental form, constant over the step, **with an in-file
  pressure floor**: `K = v0/ka * std::max(p, pc_char);` (L103) — the elastic
  stiffness never drops below the pc_char-scale value. This is the in-repo
  precedent for the absP floor knob (§7, D2).

**`@TangentOperator`** — `semiExpl` L220–236 and `absP` L263–279, identical
pattern (no brick):
```
ELASTIC/SECANT: Dt = dsig_deel;
CONSISTENT:     getPartialJacobianInvert(Je); Dt = dsig_deel * Je;
```
`constE` gets its tangent from the StandardElasticity brick. In every case the
interface exposes **one** tangent block, dσ/dΔε — nothing else (see §1.2).

**Semi-explicit volume-ratio update** (`@UpdateAuxiliaryStateVariables`):
`semiExpl` L196–203 (`v += v0 * detoV;`), `absP` L237–245
(`v += v0 * detoV;` L244, plus the explicit `pc` update L243), `constE`
L184–189 (`v *= exp(trace(deto));` L188). This explicit end-of-step update is
what makes the family "semiExpl".

**`@AdditionalConvergenceChecks`** (all three; semiExpl L205–218, absP
L247–261, constE L191–204): a converged solution with `dlp < 0` is rejected
and re-run elastically.

**Material properties:** `nu, M, ka, la, pc_char` (+ `v0` in semiExpl/absP;
`young` in constE; parameter `pamb` in constE). State variables as listed
above. **No `@ExternalStateVariable` is declared in any of the three** —
the implicit DSL provides only the default temperature ESV.

### 1.2 How OGS threads external state variables into small-strain MFront
(`MaterialLib/SolidModels/MFront/`)

- The generic adapter `MFrontGeneric<DisplacementDim, Gradients, TDynForces,
  ExtStateVars>` (MFrontGeneric.h) is instantiated for plain solid mechanics
  as `MFront<DisplacementDim>` with
  `mp_list<Strain>, mp_list<Stress>, mp_list<Temperature>` (MFront.h L15–25).
- **Hard limitation, both compile- and run-time** (MFrontGeneric.h):
  ```
  static_assert(std::is_same_v<ExtStateVars, mp_list<Temperature>>,
      "Temperature is the only allowed external state variable.");   (L368–370)
  if (_behaviour.esvs[0].name != "Temperature") OGS_FATAL(
      "Only temperature is supported as external state variable.");  (L372–379)
  ```
- ESV values are copied **explicitly, per increment, by position [0]**
  (MFrontGeneric.h L516–532): `s0.external_state_variables[0] =
  variable_array_prev.temperature;` and `s1.… = variable_array.temperature;`
  — flagged in-code as a TODO to unify with gradient handling (L514–515).
- Gradients (Strain) are set from the MPL `VariableArray` by name-mapped
  metadata: `Variable.h` `struct Strain` → `mpl_var =
  &VariableArray::mechanical_strain` (L96–110); `struct Temperature` →
  `&VariableArray::temperature` (L203–211). The `mp_for_each<Gradients>`
  copy happens in `integrateStress` (MFrontGeneric.h L542–547).
- Integration: `mgis::behaviour::integrate` (MFrontGeneric.h L561–562);
  a non-1 status throws `NumLib::AssemblyException("MFront: integration
  failed with status …")` (L563–568) — this is exactly the "status −1"
  surfacing in the Task-13 logs.
- Behaviours are loaded by name from the shared library at create time:
  `mgis::behaviour::load(lib_path, behaviour_name, hypothesis)`
  (CreateMFrontGeneric.cpp L290, L305); the name comes from the PRJ
  `<behaviour>` tag (L378). The loaded `_behaviour.esvs` list is available
  at create time — usable for capability detection (D4).
- New `.mfront` files must be added to the build list
  (`MaterialLib/SolidModels/MFront/CMakeLists.txt`, `_mfront_behaviours`
  L22 ff.; current MCC entries L30–32).
- **Consequence [D]:** in the standard small-strain interface an ESV is a
  *known per-increment input*: it is theta-interpolable inside the behaviour
  (`esv + theta*desv`, MFront Implicit-DSL semantics) but carries **no global
  cross-tangent block** — the interface returns dσ/dΔε only (§1.1
  `@TangentOperator`). Cross-tangents require the generic interface
  (`@DSL DefaultGenericBehaviour` + `@TangentOperatorBlocks`), §1.4.

### 1.3 The external eigenstrain path in RichardsMechanics that must be
disabled in internal mode (`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`)

The chain, with the exact variable handoff:

1. **σ_sw accumulation** — `updateSwellingState` (L2239–2318):
   `sigma_sw.sigma_sw += computeSwellingStressIncrement<…>(…)` (L2299–2308),
   which forwards to `computeReferenceMicroPorositySwellingStressIncrement`
   (L2216–2236 → L2211–2212):
   `delta_sigma_sw += n_S*(n_l_prev*Pi_prev - n_l*Pi_curr)*identity2;`
   — **spherical by construction** (a scalar times `identity2`).
   Sign (settled, L2155–2161): Π > 0 ⇒ σ_sw = −φ_m·Π **compressive** in the
   OGS tension-positive convention.
2. **Eigenstrain conversion (volumetric, for MPL)** — same function,
   L2312–2317:
   ```
   variables.volumetric_mechanical_strain = variables.volumetric_strain
       + identity2.transpose() * C_el_inverse * sigma_sw.sigma_sw;
   ```
   (and the `_prev` analogue). Consumer: the MPL porosity model
   `MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp`.
3. **Eigenstrain conversion (tensorial, for the solid model)** — four sites,
   all of the same shape
   `eps_m = swelling_stress_active ? eps + C_el.inverse()*sigma_sw : eps;`:
   - L2861–2864 (`setInitialConditionsConcrete`, fills `eps_m_prev`),
   - L3181–3183 (`assemble`),
   - L3793–3798 (the `updateConstitutiveRelations` helper used by the
     Jacobian assembly),
   - L5245–5247 (numerical-Jacobian path).
   `swelling_stress_active = solid_phase.hasProperty(swelling_stress_rate)
   || isPotentialExchangeEnabled(…)` (e.g. L3177–3180).
4. **Handoff to MFront** — `variables.mechanical_strain.emplace<…>(eps_m)`
   (L3184–3186 / L3796–3798), then
   `ip_data.updateConstitutiveRelation(variables, …, eps_m, eps_m_prev,
   solid_material, …)` (L3207–3210 / L3817–3820), which in
   `IntegrationPointData.h` (L60–101) builds
   ```
   variable_array_prev.stress            = sigma_eff_prev->sigma_eff;  (L81)
   variable_array_prev.mechanical_strain = eps_m_prev->eps_m;          (L82–84)
   variable_array_prev.temperature       = temperature;                (L85)
   solid_material.integrateStress(variable_array_prev, variable_array, …)
                                                                       (L87–89)
   ```
   and `Strain` is read from `mechanical_strain` inside the adapter
   (§1.2). **This is the exact handoff: the swelling enters MFront ONLY as
   the `C_el⁻¹·σ_sw` addend inside `mechanical_strain`; the returned `sig`
   (= `sigma_eff`) therefore already contains σ_sw** (σ = C:ε_m = C:ε + σ_sw
   in the LE picture).
5. **Held-within-iteration precedent** — the film coupling's confining
   pressure: "p_conf is HELD FIXED across the step …: the p_conf→stress
   coupling then enters only through the outer (Newton) iteration as p_conf
   updates step to step" (L1959–1961; same statement at L2091). This is the
   convergence class the internal-swelling ESV inherits (§5).
6. **PRJ-flag and exclusivity precedent** — `film_strain_coupling`:
   parse with default at `CreateRichardsMechanicsProcess.cpp` L653–656
   (enum `FilmStrainCouplingMode` in `PotentialExchangeParameters.h` L250);
   replace-not-add semantics documented at FEM-impl L694–702 (mode ≠ Off
   REPLACES the frozen-geometry path — the D3 no-double-count precedent of
   the strained-film design, STRAINED_FILM_IMPLEMENTATION.md §6).

### 1.4 The in-repo generic-interface precedent
(`Tests/Data/ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/MFrontBehaviour/BentoniteBehaviour.mfront`)

- `@DSL DefaultGenericBehaviour` (L4) — NOT the standard small-strain
  interface.
- Drivers beyond strain: `@Gradient real p_L; p_L.setEntryName(
  "LiquidPressure");` (L17–18); fluxes beyond stress: `@Flux real Sr;`
  (L22–23); **declared cross-tangents**:
  `@TangentOperatorBlocks{dsig_ddeto, dsig_ddp_L, dSr_ddp_L};` (L25).
  It also carries a true ESV (`@ExternalStateVariable stress pr;`, L45–46).
- Consumed through the general adapter instantiation
  `MFrontGeneric<Dim, mp_list<Strain, LiquidPressure>,
  mp_list<Stress, Saturation>, mp_list<Temperature>>`
  (ProcessLib/ThermoRichardsMechanics/ConstitutiveStressSaturation_StrainPressureTemperature/CreateConstitutiveSetting.cpp
  L24–25; TraitsBase.h L18–19).
- **Why it matters here:** ModCamClay's `@DSL Implicit` + standard
  small-strain interface CANNOT export a `dsig_dds_sw` block; the
  BentoniteBehaviour pattern is the upgrade path if the explicit ESV coupling
  proves too loose (§5). It is deliberately NOT the first step: the RM
  process consumes solid models through `MechanicsBase` /
  `MFront<DisplacementDim>` (§1.2), and switching RM to a general adapter is
  a much larger change than this port.

---

## 2. A — The constitutive statement [D]

### 2.1 Sign conventions (verified from the files)

| quantity | convention | source |
|---|---|---|
| OGS effective stress σ′ (Kelvin) | tension-positive | FEM-impl L2155–2161 comment; momentum residual `B'σ_eff` L3231–3232 |
| DSM swelling eigenstress σ_sw | spherical, **compressive ⇒ s_sw ≤ 0** in tension-positive convention | FEM-impl L2155–2161, L2211–2212 |
| ESV scalar (this design) `s_sw := tr(σ_sw)/3` | tension-positive scalar, ≤ 0 during swelling | definition; σ_sw spherical (L2211–2212) makes the scalar lossless |
| MFront `sig` | tension-positive stensor (same continuum convention as OGS) | adapter copies σ ↔ `sig` with no sign flip (MFrontGeneric.h L552–559, L575 ff.) |
| MCC pressure `p` inside the behaviours | **compression-positive**: `p = -trace(sig)/3` | absP L137/L179; semiExpl L102/L147; constE L134 (`+ pamb`) |
| p_c, pc_char, M, κ, λ | compression-positive p-axis | the files' yield `f = q² + M²·p·(p − pc)` (absP L186 etc.) |

Bridge: a tension-positive mean stress σ_m maps to p = −σ_m. The compressive
swelling scalar s_sw ≤ 0 corresponds to a compression-positive swelling
pressure s̄_sw := −s_sw ≥ 0.

### 2.2 The statement

New behaviour variant (per parent; D3 decides which parents):
**total effective stress returned to OGS**

```
sig = sig_MCC(eel)  +  s_sw_mid · I ,        s_sw_mid = s_sw + theta·ds_sw
```

where `s_sw` is a scalar `@ExternalStateVariable` (stress-typed), supplied by
the RM process per increment from the untouched DSM eigenstress site
(D1 = option (i), §3), theta-interpolated inside the behaviour per the
Implicit-DSL ESV semantics (value at step start `s_sw`, increment `ds_sw`).
`sig_MCC(eel)` is the parent's elastic law, untouched. The **MCC machinery
(elastic predictor, yield, flow, p_c evolution) operates on the MCC effective
stress**

```
σ_MCC = σ − s_sw_mid · I        (tension-positive)
p_MCC = −tr(σ_MCC)/3 = p_tot + s_sw_mid = p_tot − s̄_sw_mid   (compression-positive)
q_MCC = q(σ_MCC)     = q(σ)     [deviator of a spherical tensor is zero]
```

i.e. **the prompt-level statement p′_eff = p_total − s_sw** (tension-positive
means) is, on the compression-positive p-axis of the .mfront files, the
isotropic shift `p ↦ p − s̄_sw`. The yield function and p_c evolution become

```
f      = q² + M² · p_MCC · (p_MCC − pc_new)        (replaces absP L186 / semiExpl L152 / constE L139)
df_dp  = M² · (2·p_MCC − pc_new)                   (replaces absP L187 …)
pc_new : UNCHANGED in form (driven by deplV only;  absP L184, semiExpl/constE frpc)
```

### 2.3 The change is an isotropic shift — return-map structure unchanged [D]

Within one local solve, `s_sw_mid` is a **known constant** (ESVs are inputs;
they do not depend on the integration variables `deel, dlp, drpc`). Hence:

- `∂p_MCC/∂(deel) = ∂p_tot/∂(deel)` — identical to the parent;
- the flow direction `n = df_dsig/‖·‖` and every Jacobian block
  (`dfeel_ddeel, dfeel_ddlp, dflp_ddeel, dflp_ddrpc, …`; absP L227–233,
  semiExpl L176–187, constE L167–180) keep their parent FORM with `p`
  replaced by `p_MCC` — the local Newton/LM matrix is the parent's matrix
  evaluated at a translated stress point;
- the elastic predictor check `f_el < 0` (absP L158–159, semiExpl L122–123,
  constE L114–115) is evaluated at the shifted predictor
  `p_el,MCC = p_el − s̄_sw_mid`;
- **zero-ESV reduction:** `s_sw = ds_sw = 0` reproduces the parent behaviour
  exactly, line for line (verification anchor §8.1).

The two mathematically equivalent implementation styles (they differ in which
pressure the absP K-law sees — that coupling is DECISION D5, §7):

- **Style S1 — internal σ_MCC bookkeeping (minimal diff).** At
  `@InitLocalVariables` subtract the incoming shift from the stored stress
  (`sig -= s_sw·I` — the stress handed in by OGS is the TOTAL σ′, set from
  `sigma_eff_prev` at IntegrationPointData.h L81 via MFrontGeneric.h
  L552–559); run the ENTIRE parent algorithm on σ_MCC unchanged; add
  `(s_sw + ds_sw)·I` back at final-stress time (`@ComputeFinalStress` /
  end of `@ComputeStress` chain). absP's `p0` (L137) then IS p_MCC ⇒ the
  hypoelastic K sees the swelling-free pressure (D5 option α).
- **Style S2 — total-stress bookkeeping.** Keep `sig` total throughout;
  substitute `p → p + s_sw_mid` at every invariant evaluation
  (absP L137/L179, predictor sites, yield/flow lines). absP's `p0` stays
  total ⇒ K sees the swelling-compressed pressure (D5 option β).

### 2.4 Honest semantic note — what changes relative to the EXTERNAL mode [D]

In today's eigenstrain path the behaviour receives
`eps_m = eps + C_el⁻¹·σ_sw` (§1.3 item 3), so its `sig` — and therefore the
yield function — already CONTAINS the swelling compression
(σ = C:ε + σ_sw in the LE picture). The internal-swelling statement of §2.2
deliberately evaluates yield on σ_MCC = σ − s_sw·I, i.e. **the swelling push
is removed from the yield argument**. This is not a pure refactor; it is the
constitutive content of the port: the disjoining eigenstress loads the
momentum balance but is no longer allowed to drive the clay across its own
yield surface, which is precisely what forced `InitialPreConsolidationPressure0`
up to 12 MPa ("Was 1e6 → overrun by the 11.46 MPa swelling at IC, MFront MCC
crash. Now > swelling." — PRJ comment quoted in TASK13_MCC_BLOCKAGE.md) and
made plasticity dormant. With the shift, p_c can return to a physical
(Mont-Terri-candidate-scale) value without the IC overrun. The yield-argument
choice is Vinay's constitutive call and is restated explicitly under D1/D5
(§7) for sign-off — it is the formulation decision of this design, not an
implementation detail.

---

## 3. B — Driver choice [DECISION D1]

What does the RM process hand the behaviour?

**(i) ESV = s_sw directly from the DSM eigenstress site — RECOMMENDED.**
The process computes σ_sw exactly as today (`updateSwellingState`,
L2299–2308, untouched — Pi-path, Maxwell pair, live K(ρ_d), film modes all
stay the single owner of the closure) and publishes the scalar
`s_sw = tr(σ_sw)/3` (lossless: σ_sw is spherical, L2211–2212) into the
`VariableArray` instead of converting it to an eigenstrain.
*Costs:* one new `VariableArray` scalar member (+ enum entry,
`MaterialLib/MPL/VariableType.h` L165–194 — no suitable member exists today),
one `Variable.h` metadata struct, the §1.2 adapter generalization (relax the
static_assert L368–370 and replace the positional copy L516–532 with an
`mp_for_each<ExtStateVars>` name-mapped copy, symmetric to the Gradients
handling L542–547), the new `.mfront` variants. *Wins:* most conservative
change; one closure owner; mtest-able standalone; the eigenstress evolution
(incl. every future DSM refinement) reaches the behaviour with zero further
MFront work.

**(ii) ESV = Π, behaviour forms σ_sw = −φ_m·(Π − b·p_conf) internally.**
Needs φ_m, b, p_conf as ADDITIONAL ESVs (4 scalars), duplicates the
eigenstress closure inside MFront, and every closure change (live K(ρ_d) —
K_OF_RHO_D_LIVE.md; strained film — STRAINED_FILM_IMPLEMENTATION.md;
REV-referencing fixes) would need a parallel MFront edit. Violates
single-mechanism ownership ("name the mechanism and say which model carries
the coupling"). REJECTED unless Vinay wants the behaviour self-contained for
export to non-OGS hosts.

**(iii) ESV = suction p_L (BBM-classical).** The behaviour would carry its
own s_sw(p_L) law ⇒ duplicates retention AND swelling, directly conflicting
with the DSM (the micro engine already owns Π); guaranteed double-count
hazard. REJECTED.

**(iv) (surfaced during code reading) Process-side stress shift, no MFront
change at all.** Keep the parent behaviours; in
`IntegrationPointData.h::updateConstitutiveRelation` subtract `s_sw_prev·I`
from `variable_array_prev.stress` (L81) before the call and add
`s_sw·I` to the returned `sigma_eff` (L97) after it. Mathematically equal to
(i)+S1 with end-of-step (not theta-interpolated) s_sw. Zero MFront and zero
adapter work — but: not exercisable in mtest, forecloses D5(β) and any later
internal use of the ESV, and smears the constitutive statement across two
code layers. Worth keeping in the pocket as a one-day DIAGNOSTIC probe of the
shock-removal hypothesis before the full port is built; NOT the deliverable.

**Recommendation: (i)**, with (iv) optionally run first as a cheap
falsification probe of motivation claim (a) [PRED until run].

---

## 4. C — Tangent treatment [D + honest limitation]

- The standard small-strain interface returns exactly one tangent block,
  dσ/dΔε (§1.1 `@TangentOperator`; §1.2 adapter). An ESV carries **no global
  cross-tangent**: there is no `dsig_dds_sw` for the global Newton.
- Therefore the global Newton sees `s_sw` **explicitly — held within the
  iteration, updated between iterations/steps** — exactly the convergence
  class of the film coupling's `p_conf` today: "the p_conf→stress coupling
  then enters only through the outer (Newton) iteration as p_conf updates
  step to step" (FEM-impl L1959–1961; L2091).
- The block dσ/dΔε itself is UNCHANGED by the port [D]: within the local
  solve `s_sw_mid` is constant, so the consistent tangent of the intSwell
  variant equals the parent's tangent at the shifted state (§2.3). The
  u–p_L global coupling through σ_sw (∂σ_sw/∂n_l ↔ the Maxwell pair) is
  neither better nor worse than today — it lives in the process Jacobian,
  not in the behaviour, and is untouched.
- **[PRED]** Expected consequence: same global convergence class as the
  current film coupling / current eigenstrain mode (the eigenstrain path is
  equally explicit in σ_sw); what the port fixes is the INTERNAL return map
  (the status −1 wall), not the outer coupling. Not verified until the §8
  runs.
- **Upgrade path if the explicit coupling proves too loose:** port the
  intSwell variant to `@DSL DefaultGenericBehaviour` with
  `@TangentOperatorBlocks{dsig_ddeto, dsig_dds_sw}` after the
  BentoniteBehaviour pattern (§1.4), consumed through a general adapter
  instantiation as in TRM CreateConstitutiveSetting.cpp L24–25. That is a
  separate design (RM consumes `MechanicsBase`, §1.2) — do not start there.

---

## 5. D — No-double-count wiring

The cardinal rule: σ_sw must act **exactly once**. Internal mode therefore
switches the application site, never adds one.

1. **PRJ flag [DECISION D4]:** new optional tag in `<potential_exchange>`:
   ```xml
   <swelling_application>eigenstrain | internal</swelling_application>
   ```
   default `eigenstrain` (bit-for-bit today). Precedent for tag + enum +
   default + replace semantics: `film_strain_coupling`
   (CreateRichardsMechanicsProcess.cpp L653–656;
   PotentialExchangeParameters.h L250; FEM-impl L694–702;
   STRAINED_FILM_IMPLEMENTATION.md §6/D3-D4).
2. **When `internal`:**
   - `updateSwellingState` keeps accumulating `sigma_sw` (L2299–2308) — the
     state, output fields and the Maxwell-pair bookkeeping stay;
   - it STOPS adding `C_el⁻¹·σ_sw`: the four `eps_m` sites (L2861–2864,
     L3181–3183, L3793–3798, L5245–5247) use `eps` directly (the
     `swelling_stress_active` ternary gains the mode);
   - it PUBLISHES `s_sw = tr(σ_sw)/3` and `s_sw_prev` into
     `variables` / `variables_prev` (new VariableArray member, §3(i));
   - `IntegrationPointData.h::updateConstitutiveRelation` threads the prev
     value into the `variable_array_prev` it builds locally (L80–85) — this
     needs `sigma_sw_prev` (or the pre-filled member) passed in: small
     signature extension, listed in §9.
   - **Open sub-decision D4b:** the `volumetric_mechanical_strain` published
     for MPL at L2312–2317 (consumer:
     `TransportPorosityFromMassBalance.cpp`). Keep the swelling shift in it
     (porosity evolution unchanged vs external mode) or drop it (purely
     kinematic)? This is a physics call — which strain the porosity balance
     should see is not decided by where the stress is applied. ASK Vinay;
     default proposal: KEEP (unchanged porosity evolution isolates the port's
     effect to the return map).
3. **Create-time exclusivity check** (precedent: film_strain_coupling D3
   exclusivity, STRAINED_FILM_IMPLEMENTATION.md §6):
   - `internal` mode + a behaviour that does NOT declare the s_sw ESV →
     OGS_FATAL (the swelling would be applied NOWHERE);
   - `eigenstrain` mode + a behaviour that DOES declare it → OGS_FATAL (risk
     of double application if the PRJ also wires the ESV, and the behaviour
     would otherwise integrate with a NaN/zero ESV silently);
   - `internal` + `solid_phase.hasProperty(swelling_stress_rate)` (the MPL
     swelling law in the same ternary, L3177–3180) → OGS_FATAL (two swelling
     owners).
   The behaviour's declared ESV list is available at create time
   (`_behaviour.esvs`, loaded at CreateMFrontGeneric.cpp L290) — detection is
   by **declared-ESV capability, by name**, which is robust against behaviour
   renames; the PRJ flag stays the explicit intent switch and the two are
   cross-checked (D4).

---

## 6. E — Robustness companions [DECISION D2]

Both as separate, PRJ/parameter-visible knobs, **default off**, so the port's
effect is attributable (one mechanism per switch):

1. **`@Algorithm` — the LM half of this question is MOOT (finding §1.1):**
   LevenbergMarquardt is ALREADY the active algorithm in absP (L35, since
   2026-05-28) and constE (L34, since 2026-06-06), and both predate the
   Task-13 binary (`23a723cc3c`). The −1 walls were hit *with* LM. The new
   intSwell variants inherit LM from their parents; no flip to decide.
   Remaining sub-question: keep the commented `NumericalJacobian` (semiExpl
   L36 comment) as a documented DIAGNOSTIC-only switch (catches hand-Jacobian
   errors in the ported residuals during §8.3 FD checks) — recommend YES as a
   comment, never as shipped default.
2. **p-floor / minimum bulk modulus for absP:** absP's K has no floor
   (helper L111–112 `p = p0*exp(-v0_ka*deelV); K = v0_ka*p;`, init L138
   `K = v0/ka * p`). The sibling `semiExpl` carries the in-file precedent
   `K = v0/ka * std::max(p, pc_char);` (L103) — floor scale = `pc_char`, an
   EXISTING material property (no new literal; §1.1-style provenance
   preserved). Proposal: optional `@MaterialProperty stress p_floor`
   (default 0 = off ⇒ bit-compatible) clamping `p` at exactly the two K
   sites (helper L111–112 and init L137–138):
   `p_eff_K = std::max(p, p_floor)`. Whether the floor argument should be
   `pc_char` (semiExpl precedent), a fraction of it, or a new value is a
   §1.1 material-parameter decision → Vinay (no candidate value proposed
   here). NOTE [D]: under D5(β) (K from total p incl. swelling compression)
   the floor may be redundant once swelling builds — but the IC step (s_sw ≈
   0, p0 = σ0) is exactly where it is needed; the two knobs are complements,
   not substitutes.

---

## 7. F — Which variants get the port [DECISION D3] + the surfaced D5

- **absP — FIRST and primary. Gains the most:** its two failure facets are
  precisely what the port addresses — (a) the step-1 eigenstrain shock
  through the soft K (gone by construction: Δσ_sw leaves `deto`), and (b)
  under D5(β) the swelling compression actively keeps the K-law pressure
  compressive, away from the p→0 singularity — **[PRED] internal compressive
  s_sw may cure the step-1 singularity outright**; with the E.2 floor as
  IC backstop. absP is also the physically preferred elasticity for the
  swelling regime (κ-line, pressure-dependent stiffness) — the variant the
  Task-13 campaign actually wants to run.
- **constE — SECOND.** Cheap cross-check (no K singularity, so it isolates
  the yield-traversal facet: if intSwell-constE traverses yield where
  external constE cascaded −1, the shock-removal mechanism is confirmed
  independently of the absP elasticity). Implementation cost note: the
  `@Brick StandardElasticity` (L36) owns the elastic prediction and stress;
  injecting the `s_sw·I` shift cleanly likely means dropping the brick and
  hand-writing the (constant-E) elasticity as in the other two parents —
  a larger diff than absP's. Acceptable; flagged.
- **plain semiExpl — NO port this iteration.** Not used by the Task-13
  campaign; keep as the untouched zero-s_sw reference parent for §8.1.

**[DECISION D5 — surfaced by the derivation, §2.3/§2.4]: which pressure does
the absP hypoelastic K see?**
(α) p_MCC (swelling-excluded; falls out of style S1) — conservative,
"stiffness belongs to the skeleton net of eigenstress"; keeps K soft at the
dry IC (relies on shock-removal + floor alone) [PRED];
(β) p_total (swelling-included; falls out of style S2) — "the disjoining
pressure confines and stiffens the skeleton on the κ-line"; carries the full
(b)-claim of the motivation [PRED].
This is a constitutive call (does recompression stiffness see the
disjoining-pressure confinement?) → Vinay. The yield shift (§2.2) is the same
under both; only the K argument differs.

---

## 8. G — Verification plan (anchors per CLAUDE.md §3 — structure only, no
invented expected values)

mtest precedent: sibling `.mtest` single-Gauss-point cases exist in the same
directory (e.g. `MohrCoulombAbboSloan_*.mtest`); mtest supports
`@ExternalStateVariable`.

1. **Zero-ESV reduction** [anchor (f): regression baseline — the parent
   behaviours, the shipped suite]: mtest, `s_sw ≡ 0` over an elastic+plastic
   strain path; intSwell variant must be **bit-identical** to its parent
   (same residuals evaluated, §2.3). Discriminating for any accidental
   formulation drift in the port.
2. **Free-strain s_sw ramp** [anchor (a): analytical limit]: mtest, total
   strain held at 0, ramp `s_sw` from 0 (elastic regime, yield never
   reached): the returned total stress must equal `s_sw·I` **exactly**
   (σ_MCC(eel=0) = σ0 = 0 ⇒ σ = s_sw·I, §2.2). Blind to every MCC parameter.
3. **FD check of the consistent tangent** [anchor: derived identity]:
   FD-vs-analytic dσ/dΔε at states on both sides of yield with `s_sw ≠ 0`;
   also FD-confirm ∂σ/∂s_sw = I numerically at the mtest level (the identity
   the global explicit coupling relies on). The `NumericalJacobian`
   diagnostic switch (§6.1) cross-checks the internal Jacobian.
4. **Task-13 stress-path traversal — THE discriminating test** [anchor (e):
   published/recorded benchmark state — the measured GP state of
   TASK13_MCC_BLOCKAGE.md (p′ = 1.6 MPa, q = 3.2 MPa, p_c = 6 MPa probe),
   cited, not invented]: mtest strain path + s_sw history reconstructing the
   stage-1 approach to yield on the dry side. PASS criterion (structural,
   no value asserted): the integrator returns status 1 through yield onset
   and `dlp > 0` engages — exactly where the external-mode constE probe
   cascaded −1 with ε_p ≡ 0. Quantitative endpoint assertions:
   `TODO(Vinay)`.
5. **OGS level — eigenstrain-mode regression** [anchor (f)]: the existing
   MCC RM cases (the Task-13 PRJs of `task13_absp_2026-06-11` /
   `task13_mcc_probe_2026-06-11`, and any registered MCC benchmarks) run in
   `eigenstrain` (default) mode on the new build: **bit-for-bit** unchanged
   (the flag defaults off; same isolation discipline as the strained-film
   off-mode regression, STRAINED_FILM_IMPLEMENTATION.md §10 — two binaries,
   one input).
6. **OGS level — Task-13 stage-1 internal-swelling run** vs the recorded
   1.459 external-eigenstrain baseline. **§2 discipline: the DD endpoint is
   VERIFICATION, not calibration** — no parameter (K, p_c, λ, κ, M, floor)
   may be tuned against the 1.4139 experimental value in the same exercise;
   the run answers Vinay's mechanism question (does dry-side plastic
   dilation walk DD toward 1.414?) as a prediction test [PRED until run].
7. **MS LE standard suite untouched:** ANCHORS_MS33 Model I/III/IV/VII LE
   PRJs do not use MCC; still re-run as release gate per CLAUDE.md §12.3
   (default-mode build regression only).
8. **Crawl policy applies:** slow-but-converging internal-mode runs get the
   tiered crawl-≠-stall diagnosis (Δt trace, reject count) before any knob is
   touched — per `project_tuller_hydraulic_stiffness_tiers` /
   `modelling-method.md`; no tolerance or algorithm change without the tier
   walk.

---

## 9. H — Implementation plan (numbered code sites)

Dependency order; 1–3 are MFront-side, 4–8 OGS-side, 9–10 tests/docs.

1. **New behaviour file(s)**:
   `MaterialLib/SolidModels/MFront/ModCamClay_semiExpl_absP_intSwell.mfront`
   (D3: absP first; constE variant second). Copy of parent + (i) scalar
   `@ExternalStateVariable stress s_sw;` with entry name per D4, (ii) the
   §2.2/§2.3 shift in the chosen style (S1 vs S2 per D5), (iii) header
   comment documenting the constitutive statement, conventions and this
   design doc. Zero-ESV reduction property is the review criterion of the
   diff.
2. **Build registration**: add the new name(s) to `_mfront_behaviours`
   (`MaterialLib/SolidModels/MFront/CMakeLists.txt` L22 ff., MCC block
   L30–32).
3. **mtest cases** (§8.1–8.4) next to the behaviour, following the
   `MohrCoulombAbboSloan_*.mtest` precedent.
4. **MPL variable**: new scalar member + enum entry + name string in
   `MaterialLib/MPL/VariableType.h` (member list L165–194; enum L20 ff.;
   names array L53–55 ff.) — e.g. `swelling_stress_mean` (name → D4).
5. **Adapter generalization** (`MaterialLib/SolidModels/MFront/`):
   - `Variable.h`: new metadata struct (MFront name ↔ `mpl_var`), after the
     `Temperature` pattern (L203–211);
   - `MFrontGeneric.h`: relax the static_assert (L368–370) and the runtime
     name check (L372–379) to accept the extended `ExtStateVars` list;
     replace the positional ESV copy (L516–532) with name-mapped
     `mp_for_each<ExtStateVars>` (mirroring Gradients, L542–547). ESVs
     remain *inputs without tangent blocks* — no interface change;
   - `MFront.h`: instantiate the solid-mechanics wrapper with
     `mp_list<Temperature, SwellingStressMean>` (L15–25) — NOTE this list is
     a superset check: behaviours declaring only Temperature must keep
     working (most do); the check loop must tolerate behaviours that declare
     a *prefix* of the list. (Smallest viable generalization; alternative —
     a second wrapper type only for intSwell — if the tolerant check turns
     out invasive.)
6. **Process flag**: parse `<swelling_application>` in
   `CreateRichardsMechanicsProcess.cpp` (pattern of L653–656), enum +
   default in `PotentialExchangeParameters.h` (pattern of L250); create-time
   exclusivity checks of §5.3 (behaviour ESV capability via
   `_behaviour.esvs` after load, CreateMFrontGeneric.cpp L290).
7. **FEM wiring** (`RichardsMechanicsFEM-impl.h`): in `internal` mode,
   (a) publish `s_sw`/`s_sw_prev` from `sigma_sw`/`sigma_sw_prev` into the
   VariableArrays in `updateSwellingState` (L2299–2317) and at the four
   `eps_m` sites' scope; (b) gate the `C_el⁻¹·σ_sw` addend OUT of the four
   `eps_m` sites (L2861–2864, L3181–3183, L3793–3798, L5245–5247);
   (c) D4b treatment of `volumetric_mechanical_strain` (L2312–2317) per
   Vinay's call.
8. **`IntegrationPointData.h`**: thread the prev ESV into the locally built
   `variable_array_prev` (L80–85) — extend `updateConstitutiveRelation`
   signature (or pre-fill a passed-in `variables_prev`); same for
   `computeElasticTangentStiffness` (L30–58) only if the intSwell behaviour
   is ever used for C_el extraction (it is — same `solid_material_`; the
   zero-ESV state is fine there since C_el is evaluated as ELASTIC operator,
   but the NaN-poisoning rule means the member must be SET — to the current
   s_sw — not left NaN).
9. **Tests**: ctest registration of the §8.5/8.6 OGS cases is a SUPPLEMENT
   to the existing matrix (CLAUDE.md §3.4) and waits for §12-compliant PRJ
   headers; unit tests for the adapter generalization (a loaded behaviour
   with 2 ESVs round-trips both).
10. **Docs**: this file's STATUS section; `DSM/AGENTS.md` entry;
    ProcessLib docs page for the new tag; memory-file pointer.

**Expected diff scope [PRED, not a promise]:** 2 new .mfront files + 1 CMake
line; ~4 focused OGS files (VariableType.h, Variable.h, MFrontGeneric.h/
MFront.h, CreateRichardsMechanicsProcess.cpp, PotentialExchangeParameters.h,
FEM-impl.h, IntegrationPointData.h); no change to PotentialExchange.h, the
micro solve, or any existing PRJ.

---

## 10. I — DECISIONS for Vinay

- **D1 — driver**: ESV = `s_sw` from the unchanged DSM eigenstress site
  (option (i), §3). Includes sign-off on the §2.4 semantic point: yield
  operates on σ_MCC = σ − s_sw·I (swelling removed from the yield argument —
  the formulation content of the port). Option (iv) process-side shift
  available as a cheap pre-port falsification probe — run it? (yes/no).
- **D2 — robustness knobs** (§6): (a) LM is already active in absP/constE —
  moot; keep `NumericalJacobian` as documented diagnostic only? (b) absP
  p-floor knob: adopt? floor source = `pc_char` (semiExpl L103 in-file
  precedent) or a Vinay-supplied value — §1.1 decision, no candidate proposed.
  Both knobs default OFF.
- **D3 — variant scope** (§7): absP first (primary), constE second
  (brick-drop cost flagged), plain semiExpl not ported. Approve order/scope.
- **D4 — naming + wiring** (§5): behaviour names
  (`ModCamClay_semiExpl_absP_intSwell` working name), ESV entry name
  (e.g. `SwellingStressMean`), PRJ tag
  `<swelling_application>eigenstrain|internal</swelling_application>`
  (default `eigenstrain`), exclusivity = PRJ flag cross-checked against
  declared-ESV capability (not behaviour-name string matching).
  **D4b (sub-decision):** in internal mode, does
  `volumetric_mechanical_strain` for the MPL porosity balance keep the
  C_el⁻¹·σ_sw shift (proposal: yes) or drop it?
- **D5 — absP hypoelastic pressure argument** (§7, surfaced by §2.3):
  K from p_MCC (α, style S1) or from p_total including the swelling
  compression (β, style S2). Constitutive call; β carries motivation claim
  (b) in full, α is the conservative reading.

---

## 11. STATUS

- **2026-06-11: DESIGN ONLY — NOTHING IMPLEMENTED.** This file is the
  deliverable. No code, build, PRJ, or test was changed; no run performed.
  All [PRED] tags above are unverified predictions awaiting the §8 program.
- Blocking: D1–D5 (Vinay). After D1/D5, items §9.1–9.3 (MFront side) can
  proceed independently of D4b.
- Surprises found during code reading, already folded in: (1) LM already
  active in absP/constE under the Task-13 binary (§1.1) — retires half of
  D2; (2) the adapter hard-pins Temperature as the ONLY ESV
  (MFrontGeneric.h L368–379) — the adapter generalization (§9.5) is a
  mandatory part of option (i), not an optional nicety; (3) absP has no
  `fpc`/`frpc` residual (explicit p_c update, L184/L243) — the port's
  residual change reduces to the `flp`/predictor shift there; (4) plain
  semiExpl already floors its K at `pc_char` (L103) — in-file precedent for
  the E.2 knob with an existing, cited parameter.

---

## D1 AMENDMENT (Vinay-approved, 2026-06-12) — under film_energy_route=exact, option (ii) supersedes option (i)

Context: this doc was drafted against the operational strained-film cut. The
branch now carries the EXACT one-Psi route (form (b), mu(Pi(n_l, eps_v));
dsm_native_Pi_fofnlev), under which sigma_sw depends on eps_v WITHIN the
mechanical iteration. That inverts the D1 cost-benefit:

- Option (i) (pass precomputed s_sw as ESV) now inserts a one-iteration LAG
  exactly in the regime where form (b) is stiff (compressed contact corners —
  the measured III gap-closure failure at t~15 d, eps_v~-0.73, hyperbolic
  dPi/deps_v; see ms33_formB_suite_2026-06-12 forensics: tolerance probe and
  1W-floor probe both refuted, stiffening+lag is the surviving diagnosis).
- Option (ii) (pass the Pi-law parameters + n_l as ESVs; MCC evaluates
  s_sw(n_l, eps_v) from ITS OWN strain inside the implicit solve) makes
  MFront's implicit DSL differentiate s_sw w.r.t. strain automatically — the
  consistent tangent comes FOR FREE, the lag vanishes, and the corner
  stiffening sits inside the Newton that must track it.

AMENDED RECOMMENDATION: with film_energy_route=exact, implement option (ii).
Option (i) remains the minimal choice only for form (a) (frozen-Pi) couplings.
Composition argument (2026-06-12 discussion): under form (b), sigma_sw
self-relaxes under the compression that plastic consolidation produces — a
negative feedback into the integrator increment, softening exactly the
eigenstress shock that killed constE at yield onset [PRED, testable today:
MS33 MCC ABSP suite III/IV under form (b) vs (a), before the port exists].
Extra interface cost vs (i): Pi-law parameters + n_l as ESVs (the
Temperature-only ESV unpinning, plan item 9.5, covers this anyway).

### D1-amendment [PRED] — TESTED 2026-06-12: NOT confirmed

MS33 MCC ABSP III/IV under forms (a) vs (b) (floored, 2x2; record:
ogs/ms33_mcc_AB_2026-06-12_failed/MCC_AB_RESULTS.md): form (b) did NOT
outlast form (a) — III died 2.75 vs 3.15 d, IV byte-identical step-#1 death.
ROOT CAUSE OF THE PREDICTION FAILURE: the MCC yield onset sits at ~3 d where
eps_v is still small, and there the truncation identity makes (a)=(b) — the
self-relaxing eigenstress cannot act in a regime where the forms coincide.
The prediction was applied outside its validity domain (an avoidable error:
the bitwise dd1600 control already implied it).

STANDING FINDINGS from the same test: (1) Model III is the FIRST MS33 run
documented to reach the MCC yield surface (both forms enter the plastic
corrector at ~3 d) — the integrator's "Negative plastic increment!" failure,
not the physics, blocks traversal. (2) The wall is ModCamClay_semiExpl
itself, form-independent — consistent with TASK13_MCC_BLOCKAGE.md. The
option-(ii) port + integrator robustness work is therefore THE path; no
form choice substitutes for it.
