# Figure provenance — every plotted curve and its source

BGR curves are reductions of solver output (no fitting). Team curves are read
verbatim from the EURAD-2 data-collection `.xlsx` (Family-A master sheet, columns
located by header-string match per `TEAM_DATA_MAP.md`). Sample = (first, last) (x,y).

## Model I  (`modelI_*` figures)
BGR source: `reduced/{LE,MCC}_I_dd*_history.csv` (center probe; suction=-p/1e6, p=-(σxx+σyy+σzz)/3/1e6).
- BGR LE ρ_d=1.4: n=12, saturated P_s=5.030 MPa @ t=200 d (Sl=1.000)
- BGR LE ρ_d=1.6: n=14, saturated P_s=14.437 MPa @ t=22 d (Sl=1.000)
- BGR LE ρ_d=1.8: n=12, saturated P_s=40.863 MPa @ t=200 d (Sl=1.000)
- BGR MCC ρ_d=1.4: n=12, saturated P_s=4.922 MPa @ t=200 d
- BGR MCC ρ_d=1.6: n=12, saturated P_s=13.237 MPa @ t=200 d

Team Model I files read (saturated mean-stress endpoint + path + perm):
- **AMPOS21** (`Model_I_Amphos21.xlsx`): sat P_s by ρ_d {'1400': 4.97, '1600': 11.55, '1800': 35.06}; Dixon block: yes
- **BGE-CU-TUBAF-UFZ** (`FEM_Model_I_BGE_CU_TUBAF_UFZ.xlsx`): sat P_s by ρ_d {'1400': 6.37, '1600': 15.52, '1800': 34.08}; Dixon block: yes
- **CIMNE-UPC** (`Model_I_CIMNE_BBM_elastoplastic.xlsx`): sat P_s by ρ_d {'1400': 5.28, '1600': 13.84, '1800': 39.55}; Dixon block: yes
- **CTU-CU** (`FEM_Model_I_CTU_CU.xlsx`): sat P_s by ρ_d {'1400': 5.93, '1600': 14.79, '1800': 32.76}; Dixon block: yes
- **ClayTech** (`Model_I_Claytech.xlsx`): sat P_s by ρ_d {'1400': 4.19, '1600': 12.81, '1800': 28.3}; Dixon block: yes
- **EPFL** (`Model_I_EPFL.xlsx`): sat P_s by ρ_d {'1400': 4.8, '1600': 14.65, '1800': 40.83}; Dixon block: yes
- **ICL** (`Model_I_ICL.xlsx`): sat P_s by ρ_d {'1400': 8.14, '1600': 18.77, '1800': 36.65}; Dixon block: yes
- **IGN** (`Model_I_IGN_gamma_1e-4.xlsx`): sat P_s by ρ_d {'1400': 5.12, '1600': 12.59, '1800': 28.57}; Dixon block: yes
- **LEI** (`Model_I_LEI_updated.xlsx`): sat P_s by ρ_d {'1400': 4.89, '1600': 14.23, '1800': 40.2}; Dixon block: yes
- **PSI** (`Model_I_PSI.xlsx`): sat P_s by ρ_d {'1400': 5.04, '1600': 14.2, '1800': 40.58}; Dixon block: yes
- **TUDELFT** (`Model_I_TUDELFT.xlsx`): sat P_s by ρ_d {'1400': 4.89, '1600': 14.57, '1800': 41.84}; Dixon block: yes
- **UBERN** (`Model_I_UBERN.xlsx`): sat P_s by ρ_d {}; Dixon block: yes
- **UCLM** (`Model_I_UCLM.xlsx`): sat P_s by ρ_d {'1400': 4.21, '1600': 12.87, '1800': 37.3}; Dixon block: yes
- **ULIEGE** (`Model_I_ULiège.xlsx`): sat P_s by ρ_d {'1400': 5.04, '1600': 14.03, '1800': 39.36}; Dixon block: yes
- **VTT** (`Model_I_VTT.xlsx`): sat P_s by ρ_d {'1400': 3.5, '1600': 9.49, '1800': 43.6}; Dixon block: yes

## Model III  (`modelIII_interteam.png`)
- **CTU-CU** (`Model_III_CTU_CU.xlsx`):  |  mean stress: n=201 (0.002→2.57)  |  gap closure: n=201 (2→9.99e-05)
- **ICL** (`Model_III_ICL.xlsx`):  |  mean stress: n=401 (0.15→7.63)  |  gap closure: —
- **UCLM** (`Model_III_UCLM.xlsx`):  |  mean stress: n=41 (0.15→5.39)  |  gap closure: —
- **ULIEGE** (`Model_III_ULiège.xlsx`):  |  mean stress: n=131 (0.118→6.21)  |  gap closure: —

BGR III: mean stress from `reduced/{LE,MCC}_III_history.csv` (center); gap closure = **proxy** |u_r| at node nearest (r=0.023, z=0.040) m from VTU `displacement[0]`.

## Model IV  (`modelIV_interteam.png`)
- **AMPOS21** (`Model_IV_Amphos21.xlsx`):  |  mean stress: n=41 (0.15→8.3)  |  dry density: n=41 (0→0.226)
- **CTU-CU** (`Model_IV_CTU_CU.xlsx`):  |  mean stress: n=201 (0.002→1.69)  |  dry density: n=201 (0.9→1)
- **UCLM** (`Model_IV_UCLM.xlsx`):  |  mean stress: n=41 (0.15→1.4)  |  dry density: n=41 (0.9→1.11)

BGR IV: per-material pellet K (clay 85312.6 / pellet 13064 J/kg via <medium id=1>, Dixon 0.350 MPa). Mean stress center CSV (out_perK); dry density two zones from VTU `dry_density_solid` at the on-axis zone centroids (0,0.0525) clay [MaterialID 0, top] and (0,0.0175) pellet [MaterialID 1, bottom], ÷1000 kg/m³→g/cm³.

## Model VII  (`modelVII_interteam.png`)
- **BGE-CU-TUBAF-UFZ** (`Model_VII_BGE-CU-TUBAF-UFZ.xlsx`):  |  void ratio: n=274 (0.736→1.09)  |  axial stress: n=274 (0.2→0.262)
- **CTU-CU** (`Model_VII_CTU_CU.xlsx`):  |  void ratio: n=241 (0.724→1.09)  |  axial stress: n=241 (0.002→0.4)
- **EPFL** (`Model_VII_EPFL.xlsx`):  |  void ratio: —  |  axial stress: —

BGR VII: void ratio e=φ/(1−φ) and axial stress −σyy from `reduced/LE_VII_history.csv` (center).

### Data gaps (team quantity absent / non-standard sheet → silently omitted from plot)
- Model III gap closure: only CTU-CU exposes a `Gap closure (mm)` column; ICL/UCLM/ULIEGE omit it.
- Model VII: EPFL master sheet not in the standard time-series template (read returns none).
