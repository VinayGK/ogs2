# DSM Documentation Index — RichardsMechanics

This folder is the **single discoverable entry point** for all DSM (Dual-Structure
Model) documentation in the OGS repository. Start here.

---

## Implementation chronology

### Era 1 — Naive / additive porosity split

**Branches:** `dsm_native`, `dsm_mfront`  
**Status:** Superseded. Not the production path.

The micro water content `n_l` (state variable `MicroWaterContent`) was identified
directly with the REV-scale micro porosity:

```
phi_m := n_l          (additive / naive convention)
phi_M := phi - n_l
```

This is a simple additive split. It is physically consistent only when the macro
skeleton volume fraction `(1 - phi_M)` is unity — i.e., it ignores the reference
volume scaling between aggregate and REV.

Implemented in: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`  
Constitutive math: Mathematica notebooks `THM_DSM_Richards.nb`,
`THMDSMRichardsVK.nb`, `THM_DSM_Richards_driver.nb`,
`THM_DSM_Richards_driver_VK.nb` (see `materialmodels/src/TPM/DSM/README.md`).  
Detailed narrative: `NAIVE_POROSITY.md` in this folder.

---

### Era 2 — Hierarchical porosity split

**Branch:** `dsm_native_hierarchical` (current production branch)  
**Status:** Active. All benchmarks run from here.

The split is promoted to a proper two-level (hierarchical) relation:

```
phi   = phi_M + (1 - phi_M) * n_l        [REV-scale total porosity]
phi_m = (1 - phi_M) * n_l               [REV-scale micro porosity]
phi_M = (phi - n_l) / (1 - n_l)         [inversion for macro porosity]
```

`n_l` remains the local (aggregate-referenced) micro water content. `phi_m` is
now the REV-scale micro porosity — different from `n_l` by the factor
`(1 - phi_M)`.

Key commits on `dsm_native_hierarchical`:

| Step | Description | Commit |
|------|-------------|--------|
| 1 | REV-scale storage + split consistency | `0d7a9edd64` |
| 2 | Thermodynamic swelling stress + K recalibration | `88d42c98fd` |
| 3 | Pi-path Gibbs–Duhem consistency + flag cleanup | `c4888b6db4`, `ce9178fa96` |
| 5 | vdW dimensional fix (`/rho_lR`) + literature A lock | `0d579e8aeb` |
| 6 | DSM hardening (viscosity guards, micro-pressure density default) | `66b782afa1` |
| 7 | Dead-code removal (compatibility overload/unused flag) | `4d47efff55` |
| 8 | DSM micro-macro test refactor (13/13 passing) | `3ac6b7de1f` |

Swelling driver fix (2026-05-19): changed from `delta_phi_m` to `delta_n_l` as
the swelling slope argument — see `HIERARCHICAL_POROSITY.md` for full rationale.

Constitutive math: Jupyter notebook `Potentials.ipynb` in
`tex/cc2024/VK_B35_Pinion_May_2026/examples/DSM/` (canonical copy).  
Detailed narrative: `HIERARCHICAL_POROSITY.md` in this folder.

---

## Repository map

| Repo | Path | Role |
|------|------|------|
| `ogs` | `ProcessLib/RichardsMechanics/DSM/` | **This folder** — implementation docs, AGENTS, chronology |
| `ogs` | `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` | Core implementation file |
| `ogs` | `Tests/Data/RichardsMechanics/ANCHORS_MS33_Model*/` | Benchmark PRJ files and run data |
| `materialmodels` | `src/TPM/DSM/` | Constitutive math index, flowcharts, archive |
| `materialmodels` | `src/TPM/THM_DSM_Richards.nb` | Naive-era constitutive notebook (AceGen) |
| `materialmodels` | `src/TPM/THMDSMRichardsVK.nb` | Naive-era VK-variant notebook |
| `tex/cc2024` | `VK_B35_Pinion_May_2026/examples/DSM/Potentials.ipynb` | Hierarchical-era canonical notebook |
| `tex/cc2024` | `VK_SB_EURAD_DSM/` | ANCHORS MS33 presentation (beamer deck + data) |
| `tex/dsm-bgr-paper` | `draft/paper_DSM.tex` | BGR paper manuscript |

---

## How to pick up

1. Read `AGENTS.md` (this folder) — commit-ref roadmap and physics invariants.
2. Read `HIERARCHICAL_POROSITY.md` — current split law, swelling driver fix,
   validation checklist.
3. If you need naive-era context: read `NAIVE_POROSITY.md`.
4. For constitutive math entry points: `materialmodels/src/TPM/DSM/README.md`.
5. For benchmark status: `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/`
   (stub pointing here) and `tex/cc2024/VK_SB_EURAD_DSM/run_summary.md`.
