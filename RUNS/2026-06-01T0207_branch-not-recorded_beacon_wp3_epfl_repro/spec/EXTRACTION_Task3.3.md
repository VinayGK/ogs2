# BEACON Task 3.3 — EPFL granular MX-80 stress-path benchmark — provenance-anchored extraction

Three public sources, all on disk:

- **PRIMARY EXPERIMENT (peer-reviewed):** Ferrari, Bosch, Baryla & Rosone (2022),
  "Volume change response and fabric evolution of granular MX80 bentonite along
  different hydro-mechanical stress paths," **Acta Geotechnica**,
  DOI 10.1007/s11440-022-01481-0. Open access (CC). File:
  `refs/epfl_6777999f.txt` (from infoscience.epfl.ch bitstream 6777999f).
  → This is the authoritative measured dataset (P1-3, P2-1 + MIP/SEM).
- **TASK SPEC + INTERCOMPARISON:** BEACON D3.3 (Gens et al., 31/10/2021), §4
  "Task 3.3: verification" + Appendix-1 EPFL report (Bosch, Ferrari & Laloui,
  03/07/2019). File `spec/D3.3.txt` (genuine: `pdfinfo` → 338 pp, author
  antonio.gens, ModDate 2021-12-01). Public:
  https://www.beacon-h2020.eu/deliverables/
- **MATERIAL/WRC ANCHOR:** Seiphoori, Ferrari & Laloui (2014), Géotechnique
  64(9), 721-734 (cited throughout; numerical WRC params PENDING — to fetch).

> NOTE (2026-05-31): an earlier in-session Read of `D3.3.txt` returned spurious
> editorial strings (a bad tool result, NOT file contents). Re-extraction +
> `pdfinfo` confirmed `D3.3.txt` is the genuine 14 486-line deliverable, clean.
> No file was contaminated; nothing was renamed or discarded. The ICL/EPFL
> parameter sets below are read from the verified file.

---

## 1. Test definition (Acta Geotech §3; D3.3 App-1 §2)

Granular MX-80 bentonite, high-pressure oedometer (radially confined, drainage
top+bottom), **two HM stress paths from the SAME as-poured state to the SAME
~20 MPa axial stress**, differing only in how saturation is reached:

- **Path 1 (P1, samples P1-1/2/3):** hydrate at **constant low axial stress
  σ_a = 21 kPa** (free swelling) → then load in steps to 20 MPa.
- **Path 2 (P2, samples P2-1/2/3):** **isochoric** hydration (σ_a raised to hold
  volume, max Δε_a = ±0.6% ↔ Δρ_d ±0.01) → then load in steps to 20 MPa.

Back-pressure **u_w = 20 kPa constant in ALL stages**; σ'_a = σ_a − 20 kPa
(Acta:313). Apparatus: Ferrari, Favero & Laloui (2016) IJRMMS 88, 286-300.
Ring **h = 12.5 mm, d = 35.0 mm** (Acta:138). σ_a to 100 MPa (res 0.06 MPa),
LVDT 1 µm (Acta:147-150).

**THE TASK (D3.3:2324):** with ONE parameter set, reproduce e–vs–total-axial-
stress for P1-3 (A-B-C-D) and P2-1 (A-B'-C').

## 2. Material — granular MX-80 (Acta:101,129-141; D3.3 Table 1:2129)

| Property | Value | Source |
|---|---|---|
| Smectite content | 85 % | Acta:101 |
| Specific surface | 523 m²/g | Acta:101 |
| Specific gravity G_s | 2.74 → **ρ_s = 2740 kg/m³** | Acta:102 |
| Liquid limit / Plastic limit | 420 % / 65 % | Acta:129 |
| **Grain (aggregate) dry density** | **≈ 2.10 Mg/m³** | Acta:140 |
| Hygroscopic water content | 5–7 % | Acta:139 |
| GSD curvature / uniformity | C_C = 1.5 / C_U = 6 | Acta:133-134 |
| Pour (no-compaction) dry density | ≈ 1.5 Mg/m³ | Acta:135 |

Grain ρ_d 2.10 >> pour bulk ρ_d 1.5 ⇒ dense micro-porous aggregates loosely
packed → strong inter-aggregate (macro) porosity. **Bimodal as-poured PSD**
(macro ~10-20 µm + micro ~0.01 µm); single grain unimodal (Acta:148-151,
Seiphoori 2014).

## 3. WRC anchor (Seiphoori et al. 2014, via Acta:152-162)

As-poured (ρ_d 1.5): initial **S_r 20-25 %**, **air-entry ≈ 5 MPa total
suction**. **Adsorption-dominated**: unique water-content–suction trend to
saturation, ρ_d-independent over 1.5→1.8 Mg/m³, suction ~300→~1 MPa. (Mirrors
the adsorbed-film / free-water split — relevant to our Tuller-film work.)

## 4. Initial conditions — per sample (Acta Table 1:195-202)

| Sample | w (-) | ρ_d (Mg/m³) | e₀ (-) | S_r (-) | total suction (MPa) |
|---|---|---|---|---|---|
| P1-1 | 0.06 | 1.48 | 0.86 | 0.20 | 99.5 |
| P1-2 | 0.06 | 1.47 | 0.87 | 0.19 | 98.7 |
| **P1-3** | **0.07** | **1.50** | **0.83** | **0.23** | **90.5** |
| **P2-1** | **0.06** | **1.49** | **0.85** | **0.20** | **104.2** |
| P2-2 | 0.06 | 1.49 | 0.84 | 0.20 | 106.4 |
| P2-3 | 0.06 | 1.49 | 0.84 | 0.20 | 105.0 |

P2 pre-wetting seating stress σ_a = 190 kPa for contact (Acta:177; D3.3 main
text gives 0.12-0.31 MPa range).

## 5. MEASURED VALIDATION TARGETS

| Quantity | Measured | Source |
|---|---|---|
| Free-swell final axial strain (P1, A-B) | **≈ -80 %** (all samples) | Acta:263 |
| Saturated void ratio after free swell (B) | **2.25–2.35** | Acta:316 (D3.3 Table 4-1: 2.31-2.37, avg 2.34) |
| Swelling pressure (P2, A-B') | **3.2–3.6 MPa** | Acta:267 (D3.3: 3.12-3.55) |
| P1-3 yield stress on reloading | ≈ 0.65 MPa (stiff to 0.2-0.5; C_s' 0.05-0.08) | Acta:344,421 |
| Saturated **virgin compression index C_c** | **0.46** (Villar 0.45 oed; Tang-Cui 0.38 iso) | Acta:410-417 |
| Suction null at e | 1.62 (P1-3 after saturation) | Acta:425 |
| **VCL convergence** | curves merge at **σ_a ≈ 10 MPa, e ≈ 0.5** | Acta:372,441 |
| e at point C (path 1, σ_a 3.24 MPa) | ≈ 1.1 | Acta:446; D3.3:2293 |
| e at point B' (path 2) | ≈ 0.85 | Acta:447 |
| **STRESS-PATH-DEPENDENCY GAP at P_swell** | **≈ 0.26** (e_C − e_B') | Acta:448; D3.3 Table 4-1 |
| Δe between paths at end of saturation | ≈ 1.45 (2.3 vs 0.85) | Acta:443 |
| e at end (D / C', 20 MPa) | ≈ 0.5 (D3.3: 0.56 / 0.57) | Acta:441; D3.3 Table 4-1 |

**Unique virgin compression line / convergence (feature E, Acta §4.1):**
- Both paths merge onto **one VCL once σ_a > 12 MPa** (Acta:386); "the virgin
  compression curve at saturation is not affected by the initial fabric or stress
  path to saturation" (Acta:390). Saturated **C_c ≈ 0.46** (vs Villar 0.45 oed,
  Tang-Cui 0.38 iso) (Acta:410-417).

**Gap reconciliation (two authoritative numbers — keep both):**
- D3.3 Table 4-1 headline gap = **0.26**, defined e_C − e_B' at the swelling-
  pressure stress, derived from ρ_d 1.30 (C, e≈1.11) vs 1.50 (B', e≈0.83)
  (D3.3:1618).
- Acta refined per-sample: P1-2 (point C, σ_a 3.2 MPa) **e = 0.98**; P2-2 e =
  0.85 → direct gap **≈ 0.13** (Acta:376-378). The Acta curve is the more
  authoritative validation target; the 0.26 is what the eight teams were scored
  against. Report against both; expect the truth in 0.13-0.26.

**MIP intra/inter-assemblage boundary (for any micro/macro split calibration):**
pore Ø < **0.03 µm** = intra (micro) in all samples; as-poured bimodal peaks at
**0.01 µm (micro) / 16 µm (macro)** (Acta:421-429). P1-1 (free-swollen, e 2.35)
intra/inter boundary rises to 0.17 µm; P1-2 to 0.045 µm.

**MIP fabric story (the physical fingerprint a DSM must reproduce):**
- P1 saturate (P1-1): bimodal, dominant macro 10-20 µm + micro 0.01 µm.
- P1 saturate+compress (P1-2): macro peak ~vanishes, micro preserved → loading
  collapses MACRO only (Acta:397-403).
- P2 isochoric saturate (P2-2): partial macro closure + micro increase; vs P1-1
  → constant-volume gives higher micro, lower macro (Acta:411-424).
- P1-3 at 20 MPa: unimodal, micro only, e ≈ 0.5 = micro void ratio (Acta:425-431).
- **Micro ↔ degree of saturation (path-independent); macro ↔ stress path
  (carries the memory). Gap closes when macro is crushed out.**

## 6. Five reference features (D3.3:1790-1795): A large swelling strain · B sharp
yield · C P_swell in 2.5-4 MPa · **D path-dependency at P_swell (key)** · E VCL
convergence at high stress.

## 7. Eight-team intercomparison — measured vs models (D3.3 Table 4-1:1759)

| Team | analysis | e_A | e_B | P_swell B' (MPa) | gap | e_C' | e_D |
|---|---|---|---|---|---|---|---|
| **Experiment** | — | 0.83-0.85 | 2.31-2.37 | **3.12-3.55** | **0.26** | 0.57 | 0.56 |
| **BGR (OGS, single-struct BBM)** | BV | 0.85 | 2.39 | 2.71 | **1.33** | 0.52 | 1.62 |
| CU-CTU | BV/SE | 0.85 | 2.31 | 6.0 | 0.14 | 0.50 | 0.52 |
| Clay Tech | BV/SE | 0.85 | 2.40 | 3.72 | 0.09 | 0.58 | 0.57 |
| EPFL (ACMEG-TS, SE) | SE | 0.85 | 2.23 | 5.0 | 0.02 | 0.41 | 0.41 |
| ICL (IC-DSM, large-disp) | BV | 0.85 | 2.41 | 2.28 | 1.04 | 0.41 | 0.75 |
| LEI | BV | 0.83-0.84 | 2.27 | 3.15 | 0.90 | 0.55 | 0.61 |
| Quintessa | BV/SE | 0.83 | 2.11 | 3.29 | 0.01 | 0.34 | 0.41 |
| UPC | SE | 0.85 | 2.32 | 2.78 | 0.69 | 0.49 | 0.81 |
| ULg (Mohymar) | BV | 0.83 | 2.11 | 4.80 | 0.18 | 0.55 | 0.50 |

**BGR's own OGS attempt (single-structure BBM) overshot the gap 5× (1.33 vs
0.26) and missed convergence (D=1.62 vs 0.56)** — the documented motivation for
a double-structure model. That is the bar a native-DSM reproduction must clear.

## 8. KINEMATIC SCOPE FLAG (the central modelling fork)

Free-swell P1 A-B is e 0.83→2.3 (~80 % axial strain). EPFL: "hypothesis of
small strains… certainly beyond limits; large strains needed" (D3.3:6525). ICL
adopted a **large-displacement formulation** for P1-3 (D3.3:8190). OGS
RichardsMechanics is **small-strain** → Path-1 free-swell magnitude (feature A)
and the full B-C-D range are NOT representable without finite-strain kinematics.
Path-2 (isochoric P_swell + B'-C' compression, e 0.85→0.5) IS tractable and is
the native DSM's home turf. → User scientific call (§9), see response.

---

## §12-compliant provenance plan for an OGS native-DSM reproduction

| Group | Value | Source (family) | §12 status |
|---|---|---|---|
| ρ_s = 2740 kg/m³ (G_s 2.74) | Acta:102 / D3.3 Tab 1 | EPFL / Villar | ✓ |
| e₀, ρ_d, w, S_r, suction (per sample) | Acta Table 1 | EPFL | ✓ |
| geometry h=12.5, d=35.0 mm; u_w=20 kPa | Acta:138,167 | EPFL | ✓ |
| WRC (granular MX-80) | Seiphoori 2014 (PENDING numerical fit) | EPFL | ✓ family |
| swelling K (Pi-path) | P_swell 3.2-3.6 at ρ_d≈1.49 → Dixon MX-80 family | Dixon (§12.1) | ✓ |
| A_Hamaker = 2.2e-20 J | Israelachvili & Adams 1978 (fixed anchor, NOT a knob) | literature | ✓ |
| macro elastoplastic (λ,κ,M,p_c) | candidate cross-check vs EPFL Tab 4 / ICL Tab 4.1 | EPFL | **§2: identify on P1 NCL, VALIDATE on P2 gap — never assert the line we fit** |

CITABLE TEAM PARAMETER SETS (verified in D3.3.txt):
- **EPFL ACMEG-TS (D3.3 Table 4:6495):** κ=0.087, ν=0.25, λ_sat=0.20,
  φ'_c=φ'_e=12°, α=1.0 (MCC-equivalent), p'_r=1e-4 MPa, r=0.25, ζ=1.5, ξ=0.65;
  WRC a=0.9 MPa⁻¹, b=1.5, n=1.8, m=0.57, e_Cw,a=0.55.
- **ICL IC-DSM (D3.3 Tables 4.1-4.3:8232):** M_F=0.495, α_F=0.4, μ_F=0.9, p_c=80 kPa,
  λ(0)=0.5154, κ=0.0087, ν=0.3, κ_s=0.1415, κ_m=0.360, V_F=0.341, s₀=106 kPa;
  SWR α=1e-4 kPa⁻¹, m=0.47, n=1.90, ψ=2.0, S_r0=0.05; k_sat=3e-13 m/s.
  ICL result tables 6.1/6.2 (e, σ, suction, S_r at every load step) = digital
  reference curves for P1-3 and P2-1.

OPEN (user, §9): analysis type (SE vs BV); skeleton (LE vs MCC); large-strain
stance for Path 1. See response 2026-05-31.
