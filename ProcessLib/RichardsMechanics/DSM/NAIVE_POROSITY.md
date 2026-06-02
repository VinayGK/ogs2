# Naive / Additive Porosity Split — Historical Record

**Era:** Naive / additive (superseded)  
**Branches:** `dsm_native`, `dsm_mfront`  
**Status:** Not the production path. Preserved for reference and audit.

---

## Split law

In the naive/additive convention, the micro water content `n_l` (state variable
`MicroWaterContent`, aggregate-referenced) is identified directly with the
REV-scale micro porosity:

```
phi_m := n_l
phi_M := phi - n_l       (additive remainder)
```

This makes the porosity split a simple scalar subtraction. It is an approximation
that ignores the reference volume scaling by `(1 - phi_M)` between aggregate and
REV frames.

---

## What was implemented

**Core file:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

- `computeTransportPorosityUpdate(...)` evaluated the additive split directly.
- `micro_porosity` output field was identical to `n_l` (no REV scaling).
- Config key `macro_porosity_update_mode="additive_macro_porosity_rate_mode"` was
  the active mode on these branches.
- Storage term was `phi_m * rho_lR` with `phi_m == n_l`.

**MFront path (`dsm_mfront` branch):**  
A parallel MFront bridge (`RichardsMechanicsNotebookBridge_MCC`) used the same
naive split on the notebook side, solving a local microstate `(n_l, rho_lR)` and
passing the result to the RM carrier. Full documentation in
`materialmodels/src/TPM/THMDSMRichardsRM_MFront_transition.tex`.

---

## Constitutive math

| File | Location | Role |
|------|----------|------|
| `THM_DSM_Richards.nb` | `materialmodels/src/TPM/` | Main constitutive notebook (AceGen/AceFEM) — naive split |
| `THM_DSM_Richards_driver.nb` | `materialmodels/src/TPM/` | Driver wrapper for `THM_DSM_Richards.nb` |
| `THMDSMRichardsVK.nb` | `materialmodels/src/TPM/` | VK-variant constitutive notebook |
| `THM_DSM_Richards_driver_VK.nb` | `materialmodels/src/TPM/` | Driver wrapper for VK variant |

Call-graph flowcharts for these notebooks are in
`materialmodels/src/TPM/DSM/CONSTITUTIVE_FLOWCHART.md`.

---

## Why it was superseded

The additive split `phi_m = n_l` is inconsistent with the PINION derivation
(Nagel et al.) where `n_l` is the aggregate-frame micro water content and the
REV-frame micro porosity is:

```
phi_m = (1 - phi_M) * n_l
```

Consequences of the naive split:
- Swelling driver `delta_phi_m` underestimated swelling by factor `≈ (1 - phi_M)
  ≈ 0.5` for typical bentonite dry densities.
- The effective calibrated swelling slope was therefore off by the same factor.
- `phi_m <= phi` was not guaranteed for all constitutive states.

These issues are resolved in the hierarchical split — see `HIERARCHICAL_POROSITY.md`.

---

## Archive pointers

- MFront transition state: `materialmodels/src/TPM/DSM/ARCHIVE/THMDSMRichardsRM_MFront_state_checkpoint.md`
- VK branch state: `materialmodels/src/TPM/DSM/ARCHIVE/THMDSMRichardsVK_state_checkpoint.md`
- Minimal resume guide (MFront era): `materialmodels/src/TPM/DSM/ARCHIVE/THMDSMRichardsRM_MFront_resume_min.md`
