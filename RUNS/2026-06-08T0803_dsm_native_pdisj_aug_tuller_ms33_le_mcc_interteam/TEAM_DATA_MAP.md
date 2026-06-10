# MS33 Team Data Map — Excel ingestion reference

Read-only structural survey of the EURAD-2 MS33 benchmark data-collection
spreadsheets, for building pandas ingestion of cross-team overlay plots.
Generated 2026-06-08. Nothing was modified.

Tooling used: `/opt/anaconda3/bin/python3` with `openpyxl` (read_only,
data_only). All header row/column indices below are **0-indexed for rows**
(as openpyxl `iter_rows` yields them) and given **both as a 0-indexed
integer and as the Excel column letter** for columns. Read with
`pandas.read_excel(..., sheet_name=<master>, header=None)` then slice by
row/column, because the headers are multi-row and merged — `header=` will
not work directly.

---

## 0. Two distinct schema families (CRITICAL)

There are **two completely different layouts**. A robust parser must branch
on which one a file uses (detect by master sheet name / first-cell content).

### Family A — "18-team collection template" (the other teams)
Location: `…/g_Support_Section_Data Collection/<TEAM>_DATA/Model_<N>/…xlsx`

- **All numeric data lives in ONE master worksheet** named `Model_I_data`
  (Model I) or `Model_III` / `Model_IV` / `Model_VII`.
- Every other sheet in the workbook (`Suction_stress_path`,
  `Gap_closure`, `Axial_stress_void_ratio`, `Mean_stress`, `Porosity`,
  `Density (comparison)`, …) is an **Excel CHARTSHEET** (a plot object
  referencing the master sheet), **NOT a data table**. `openpyxl` exposes
  them as `Chartsheet` objects with no `.max_row`. **Ignore them entirely.**
  → In pandas, only ever read the single master sheet per file.
- Headers span multiple rows with merged cells and blank spacer columns
  between column-groups. Locate columns by **matching header strings in the
  header row**, not by fixed offset — offsets drift between teams
  (see irregularities, §6).

### Family B — "official per-quantity template" (BGR's own files, and the
blank master template `0_Required Figures/ANCHORS-Benchmark-Data-Collection_MODEL_I.xlsx`)
Location: `…/VK_SB_EURAD_DSM/ms33_data_collection/MS33_Model_*_BGR.xlsx`
and `…/BGR_DATA/…` (the BGR submission under the collection folder also uses
a hybrid — see §5).

- **One worksheet per physical quantity**: `Stresses`, `Suction`,
  `Saturation`, `Permeability`, plus model-specific `GapAperture`
  (Model III), `DryDensity` (Model IV), `VoidRatio` (Model VII).
- Fixed 5-row metadata header on every sheet, then data:
  - row 0: quantity title (e.g. `Total stresses`, `Mean stress (positive…)`)
  - row 1: `Points` | `<point-name>` … (e.g. `center`, `center mean`)
  - row 2: `Coordinates` | `<coord or "N/A">` …
  - row 3: `Time` | `<quantity label>` … (repeated per data column)
  - row 4: `d` | `<unit>` … (unit row; `d` = time unit = **days**)
  - row 5+: data. Col A = time (days), cols B+ = the quantity at each point.
- Wide format: official blank template pre-allocates up to ~1000 point
  columns (mostly empty); BGR files use only the first 1–3.

→ **pandas recipe Family B:** `read_excel(f, sheet_name='Stresses',
header=None, skiprows=5)`, col 0 = time(d), col 1 = value (MPa); read the
unit from cell row4/col1 and point-name from row1/col1.

---

## 1. File inventory per model

Root A = `/Users/vinaykumar/tex/eurad2_MS34/MSXX/g_Support_Section_Data Collection/`
Root B = `/Users/vinaykumar/tex/cc2024/VK_SB_EURAD_DSM/ms33_data_collection/`

Files whose basename contains `teamname` are **blank templates** (not filled)
— skip them. Listed for completeness but marked `[TEMPLATE]`.

### Model I (18 teams; FILLED = 15)
```
A/AMPOS21_DATA/Model_I/Model_I_Amphos21.xlsx
A/BGE-CU-TUBAF-UFZ_DATA/Model_I/FEM_Model_I_BGE_CU_TUBAF_UFZ.xlsx
A/BGE-CU-TUBAF-UFZ_DATA/Model_I/Model_I_teamname.xlsx                [TEMPLATE]
A/BGR_DATA/Model_I/Model_I_teamname.xlsx                             [TEMPLATE — but Family-A layout]
A/CIMNE-UPC_DATA/Model_I/Model_I_CIMNE_BBM_elastoplastic.xlsx
A/CTU-CU_DATA/Model_I/FEM_Model_I_CTU_CU.xlsx
A/CTU-CU_DATA/Model_I/Model_I_teamname.xlsx                          [TEMPLATE]
A/ClayTech_DATA/Model_I/Model_I_Claytech.xlsx
A/EPFL_DATA/Model_I/Model_I_EPFL.xlsx
A/ICL_DATA/Model_I/Model_I_ICL.xlsx
A/IGN_DATA/Model_I/Model_I_IGN_gamma_1e-4.xlsx
A/LEI_DATA/Model_I/Model_I_LEI_updated.xlsx
A/PSI_DATA/Model_I/Model_I_PSI.xlsx
A/TUDELFT_DATA/Model_I/Model_I_TUDELFT.xlsx
A/UBERN_DATA/Model_I/Model_I_UBERN.xlsx
A/UCLM_DATA/Model_I/Model_I_UCLM.xlsx
A/ULIEGE_DATA/Model_I/Model_I_ULiège.xlsx                            (NOTE: non-ASCII 'è' in filename)
A/VTT_DATA/Model_I/Model_I_VTT.xlsx
B/MS33_Model_I_dd1400_BGR.xlsx                                       (BGR, Family B)
B/MS33_Model_I_dd1600_BGR.xlsx                                       (BGR, Family B)
B/MS33_Model_I_dd1800_BGR.xlsx                                       (BGR, Family B)
```
Blank master template (defines Family B schema, no team data):
`A/0_Required Figures/ANCHORS-Benchmark-Data-Collection_MODEL_I.xlsx`

### Model III (FILLED = 3)
```
A/CIMNE-UPC_DATA/Model_III/Model_III_teamname.xlsx                   [TEMPLATE]
A/CTU-CU_DATA/Model_III/Model_III_CTU_CU.xlsx
A/CTU-CU_DATA/Model_III/Model_III_teamname.xlsx                     [TEMPLATE]
A/ICL_DATA/Model_III/Model_III_ICL.xlsx
A/UCLM_DATA/Model_III/Model_III_UCLM.xlsx
A/ULIEGE_DATA/Model_III/Model_III_ULiège.xlsx
B/MS33_Model_III_BGR.xlsx                                            (BGR, Family B)
```

### Model IV (FILLED = 3)
```
A/AMPOS21_DATA/Model_IV/Model_IV_Amphos21.xlsx
A/CIMNE-UPC_DATA/Model_IV/Model_IV_teamname.xlsx                     [TEMPLATE]
A/CTU-CU_DATA/Model_IV/Model_IV_CTU_CU.xlsx
A/CTU-CU_DATA/Model_IV/Model_IV_teamname.xlsx                       [TEMPLATE]
A/UCLM_DATA/Model_IV/Model_IV_UCLM.xlsx
B/MS33_Model_IV_BGR.xlsx                                             (BGR, Family B)
```

### Model VII (FILLED = 3)
```
A/BGE-CU-TUBAF-UFZ_DATA/Model_VII/Model_VII_BGE-CU-TUBAF-UFZ.xlsx
A/CTU-CU_DATA/Model_VII/Model_VII_CTU_CU.xlsx
A/ClayTech_DATA/Model_VII/Model_VII_teamname.xlsx                   [TEMPLATE]
A/EPFL_DATA/Model_VII/Model_VII_EPFL.xlsx
A/TUDELFT_DATA/Model_VII/Model_VII_teamname.xlsx                    [TEMPLATE]
A/UCLM_DATA/Model_VII/Model_VII_teamname.xlsx                       [TEMPLATE]
B/MS33_Model_VII_BGR.xlsx                                            (BGR, Family B)
```

Team name is reliably derivable from the parent `<TEAM>_DATA/` folder name
(more reliable than the filename, which varies: `FEM_Model_I_…`,
`…_BBM_elastoplastic`, `…_gamma_1e-4`, `…_updated`, `…_v4`).

---

## 2. MODEL I — Family A master sheet `Model_I_data`

Model I is a **steady-state suction-ramp** (NOT a time series): each row is a
suction level. The deliverable is swelling/mean-stress and permeability **vs
suction**, at three target dry densities, plus the experimental calibration
target (Dixon 2023).

Layout: 4 column-blocks laid side-by-side with blank spacer columns between.

| block | dry density | std cols (letter) | header row | meaning |
|------|-------------|-------------------|-----------|---------|
| 1 | 1.4 g/cm³ | B, C, D | row 4 | Suction, Mean stress, Intrinsic permeability |
| 2 | 1.6 g/cm³ | H, I, J | row 4 | same |
| 3 | 1.8 g/cm³ | N, O, P | row 4 | same |
| 4 | Dixon (2023) experimental | T, U, V | row 4 | Suction, **Swelling pressure**, Dry density |

Row map (0-indexed):
- row 1: `Dry density` labels (cols C, I, O); col U = `Data`
- row 2: `1.4 g/cm3` / `1.6 g/cm3` / `1.8 g/cm3` (cols C, I, O); col U = `Dixon (2023)`
- row 4: column headers (the real header row)
- row 5+: data (one row per suction level). Typically 200–1000 rows;
  real data length = count of non-null in the block's suction column
  (trailing rows can be blank padding up to ~1005).

Header strings (vary slightly per team — match case-insensitively, strip):
- `Suction - MPa`
- `Mean effective stress (r') -MPa`  **OR**  `Mean total stress (r) -MPa`
  → **effective vs total differs by team** (EPFL & BGE-CU report *total*;
  AMPOS21, CIMNE, UCLM report *effective*). Record which.
- `Intrinsic permeability (m2)`
- Dixon block: `Suction (MPA)`, `Sweling pressure (MPa)` (sic — typo
  "Sweling"), `Dry density (g/cm3)`

Units: stress & suction in **MPa**, permeability in **m²**, dry density in
**g/cm³**. Suction reported as a **positive** magnitude, descending from 100
→ ~0 MPa down the rows.

First data rows, EPFL block 1 (1.4):
```
Suction  MeanTotalStress  Permeability
100      0.15             1.3066e-20
99.9     0.157404         1.3066e-20
99.8     0.164832         1.3066e-20
```

**Comparable quantity (Model I):** swelling pressure ≡ mean stress at full
saturation (suction → 0), i.e. the LAST data row of the `Mean stress` column
in each density block; overlay vs the Dixon block (cols U/V) as the
experimental target. Also: full mean-stress-vs-suction and
permeability-vs-suction curves per density.

### Time-series sub-quantities for Model I (per Family-B BGR files only)
Model I in **Family B** (BGR) additionally carries time series: sheets
`Stresses` (mean stress vs time, d/MPa), `Suction` (MPA), `Saturation`
(noDim), `Permeability` (m²). The blank Family-B template
(`ANCHORS-Benchmark-Data-Collection_MODEL_I.xlsx`) confirms these four
sheets with the standard 5-row header. The Family-A teams do **not** provide
Model I time series (only the suction-ramp endpoints).

---

## 3. MODEL III, IV, VII — Family A master sheet (shared template)

Models III/IV/VII share ONE master-sheet template (sheet named `Model_III`,
`Model_IV`, `Model_VII`). It is a **time series** (rows = time steps; up to
~16000 rows for ULIEGE Model III). Two regions:

### 3a. Left metadata + summary block (cols A–J)
- rows 5–9, cols A–E: `Locations` table — Top / Central / Bottom with
  X,Y,Z. Standard coords: Top (0,0,70), Central (0,0,40), Bottom (0,0,10)
  [mm along the column axis Z].
- col G label + cols H,I,J values: **per-location summary table**, one row
  each. Labels (string-match in col G, value in H=Top, I=Central, J=Bottom):
  - `Final mean stress`, `Final radial stress`, `Final axial stress`
  - `Final porosity`, `Final dry density`
  - `Time (d) to reach saturation`
  - `Initial permeability`, `Final permeability`
  (exact row offset of this block varies: ~rows 9–17. Locate by the col-G
  label, not a fixed row.)
  → For quick cross-team scalar comparisons (final stress / final dry
  density per location) this summary block is the cleanest source.

### 3b. Right time-series column-groups (cols O onward)
Each quantity = one `Time`/`(days)` column followed by **3 value columns**
(Top=coord 0,0,70 / Central=0,0,40 / Bottom=0,0,10). Header rows:
- row 4: `Top` / `Central` / `Bottom`
- row 5: `Time` (for the time col) or `Cordinates : 0, 0, 70|40|10` (sic,
  "Cordinates")
- row 6: the unit/quantity label (`(days)`, `Mean stress (MPa)`, …)
- row 7+: data.

**Standard 6 groups (same column letters across III/IV/VII):**
| group | time col | value cols (Top/Cen/Bot) | row-6 label | unit |
|-------|----------|--------------------------|-------------|------|
| Mean stress    | O  | P, Q, R    | `Mean stress (MPa)`   | MPa |
| Radial stress  | V  | W, X, Y    | `Radial stress (MPa)` | MPa |
| Axial stress   | AC | AD, AE, AF | `Axial stress (MPa)`  | MPa |
| Suction        | AJ | AK, AL, AM | `Suction (MPa)`       | MPa (positive) |
| Porosity       | AQ | AR, AS, AT | `Porosity`            | – |
| Permeability   | AX | AY, AZ, BA | `Permeability (m2)`   | m² |

Time is in **days**. (ULIEGE/ICL time columns can include sub-day values
like 1.157e-4 d — keep float.)

**Void ratio** is generally NOT a standalone column in III/IV — derive from
porosity: `e = n / (1 - n)`. (Model VII is the exception, see below.)

### Model III — extra groups (after col BA)
Appended after the standard 6, at ~cols BE–BL (offsets vary per team):
- `Radial stress (MPa)` at gap location (coord 23,0,40)
- `Radial displacements (mm)` at coord 23,0,40
- **`Gap closure (mm)`** at coord 23,0,40 — Model III key deliverable
  (gap aperture closure vs time). Time col is its own `Time (days)`.
- (ICL/UCLM) also `Mean Stress` at coord 0,0,40.
- The chartsheets `Gap_closure`, `Mean_stress_comparison`,
  `Radial_displacements` plot these but hold no data.
- CTU-CU Model III: gap-closure values confirmed present (e.g. col BI,
  header row 6 `Gap closure (mm)`, time col BK).
  r=23 mm is the canister/gap radius location.

**Comparable quantities (Model III):** gap closure (mm) vs time;
mean/radial/axial stress vs time at Top/Central/Bottom; final-stress summary.

### Model IV — extra groups (after col BA)
At ~cols BE–BT (offsets vary):
- `Radial stress`/`Mean stress` at coord 20,0,40
- `Radial displacements (mm)` at 20,0,40
- **`Dry density`** at coord 20,0,40 (header row 4 `Dry density `,
  col ~BH) and **`Dry density - clay block`** at coord 0,0,40 (col ~BJ;
  UCLM labels it `density (g/cm3)`) — Model IV key deliverable
  (dry-density evolution, two zones: pellet vs clay block).
- `Reference Configuration` mean stress at 0,0,40.
- UCLM appends team-specific `Time step`, `Liquid_Pressure` (cols BS, BT).
- Dedicated chartsheets: `Density (comparison)`,
  `Mean_stress (Ref - Model IV)`, `Mean_stress (clay-pellets)` (no data).

**Comparable quantities (Model IV):** dry-density evolution (two zones) vs
time; mean/radial/axial stress vs time; porosity (→ void ratio) vs time.

### Model VII — extra group (after col BA)
At cols BD–BF, location = center:
- time col BD `Time (days)`
- **col BE `Axial Stress (MPa)`**, **col BF `Void ratio`**
- group title row 4: `Model VII : Free Swelling` (free-swelling oedometer)
- This is the **void-ratio vs axial-stress stress path** (the Model VII key
  deliverable). Dedicated chartsheet `Axial_stress_void_ratio` (no data).

**Comparable quantities (Model VII):** void ratio vs axial stress
(stress path); void ratio / volumetric strain vs time; the suction-vs-mean-
stress stress path comes from cross-plotting the `Suction` (AK–AM) and
`Mean stress` (P–R) columns at a chosen location.

---

## 4. MODEL III / IV / VII — Family B (BGR) sheet layout

BGR files (`B/MS33_Model_<N>_BGR.xlsx`) use the per-quantity Family-B schema
(see §0 Family B). Sheets and the model-specific extra sheet:

| model | sheets |
|-------|--------|
| III | `Stresses`, `Suction`, `Saturation`, `Permeability`, **`GapAperture`** |
| IV  | `Stresses`, **`DryDensity`**, `Suction`, `Saturation`, `Permeability` |
| VII | `Stresses`, **`VoidRatio`**, `Suction`, `Saturation`, `Permeability` |

Each sheet: 5-row header, col A = `Time` in `d` (days), data cols labeled in
row 1 (`Points`, e.g. `center mean` / `center axial` / `center radial` for
Model III stresses — 3 stress invariants in one sheet), row 4 = unit (MPa).
Model III `Stresses` sheet has 3 columns (mean/axial/radial) at coord
(0.0125, 0.035) m. Model IV `Stresses` has 2 columns: `bentonite (clay)`
(clay-zone mean) and `pellets centre` (pellet-zone mean). Model VII
`Stresses` title is `Mean net stress`.

→ Family-B coordinate convention is in **metres** (e.g. 0.0125, 0.035),
whereas Family-A `Locations` are in **mm** (0,0,70 etc.). Watch units when
aligning BGR vs other-team locations.

---

## 5. BGR's two submissions (disambiguation)

- `A/BGR_DATA/Model_I/Model_I_teamname.xlsx` — Family-A layout (the
  collection template), but left as a **blank `teamname` template** here
  (the dd-block headers present, data sparse). Treat BGR's authoritative
  Model I data as Root B `MS33_Model_I_dd{1400,1600,1800}_BGR.xlsx`.
- Root B `MS33_Model_*_BGR.xlsx` — Family-B, fully populated, BGR's compiled
  OGS results. **Use these for the BGR curve in overlays.**

---

## 6. Irregularities / parser hazards

1. **Chartsheets masquerade as data sheets.** Most non-master sheets are
   `Chartsheet` objects (no `.max_row`). In openpyxl, guard with
   `getattr(ws,'max_row',None)`. In pandas, `pd.read_excel` will raise or
   return junk on a chartsheet — explicitly read only the master sheet name.
2. **Multi-row merged headers** (rows 0–6 depending on family). Always read
   with `header=None` and slice; never trust auto-header detection.
3. **Blank spacer columns** between Family-A blocks/groups (e.g. cols E,F,G
   between Model I blocks; col offsets are not contiguous).
4. **Column-offset drift between teams (Family A, Model I):**
   - Standard stride = 6 cols/block (B,C,D + 3 blanks).
   - **BGE-CU-TUBAF-UFZ uses a 7-col stride** with an *extra repeated*
     `Suction - MPa` column before permeability (cols B,C,D,E then I,J,K,L…).
     → DO NOT hard-code offsets; **find columns by header-string match
     within the header row** of each block.
5. **effective vs total mean stress** differs by team (Model I) — label
   accordingly; do not assume one.
   5b. **ULIEGE Model I breaks the 3-col block stride entirely.** It inserts
   two extra permeability-variant columns per block
   (`Perm=fct(Macro-Void ratio)`, `Perm=fct(Macro-Porosity)`), giving an
   ~8-col stride: block 1 = B,C,D(,E,F), block 2 starts at **J** (not H),
   block 3 at **R** (not N). Consequently ULIEGE's cols T/U/V are NOT the
   Dixon block — they hold permeability values (~1e-20). **ULIEGE has no
   Dixon experimental block.** Confirms: find the Dixon target by
   header-string `sweling pressure` (returns empty for ULIEGE), and find
   each density block's columns by header match, never by fixed letter.
   IGN, VTT, ClayTech, AMPOS21, CIMNE, EPFL, BGE-CU DO carry the Dixon
   block (medians 5 / 14.16 / 40.6 MPa at ρ_d 1.4 / 1.6 / 1.8 g/cm³).
6. **Header typos to match on:** `Cordinates` (not Coordinates),
   `Sweling pressure` (not Swelling), `Mean stress (positive…)` truncations,
   `(MPA)` vs `(MPa)` casing. Normalize: lower-case, strip, collapse spaces.
7. **Model III/IV/VII extra-group column letters drift** (BE–BT region)
   because teams add/remove team-specific columns (UCLM Liquid_Pressure;
   ICL extra radial-displacement at 23,0,40). Locate the Gap closure / Dry
   density / Void ratio groups by **row-4/row-6 header string**, not letter.
8. **Trailing blank rows**: master sheets pre-allocate up to ~1000–16000
   rows; real data length = count of non-null in the relevant time/suction
   column. Drop NaN rows after read.
9. **Non-ASCII filenames**: `Model_*_ULiège.xlsx` (è). Glob with care.
10. **Units differ between families**: Family-A locations in mm + days;
    Family-B (BGR) coords in m + days. Stresses MPa, suction positive MPa,
    permeability m² throughout. Time in **days** in both families (no
    seconds seen).
11. **Empty/near-empty filled files**: some teams' master sheets are thin;
    always verify a quantity's column actually has data before plotting.

---

## 7. Suggested parser skeleton (Family A)

```python
import openpyxl, pandas as pd
def load_master(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    # master sheet = first sheet that is a real Worksheet (has max_row) and
    # whose name starts with 'Model_'
    master = next(s for s in wb.sheetnames
                  if getattr(wb[s], 'max_row', None) and s.startswith('Model_'))
    return pd.read_excel(path, sheet_name=master, header=None)

def find_col(df, header_row, needle):
    row = df.iloc[header_row].astype(str).str.lower().str.strip()
    hits = row[row.str.contains(needle.lower(), na=False)].index.tolist()
    return hits  # list of column integer-indices; match per block
```
For Model I: header_row = 4, needles `suction - mpa`, `mean`, `permeability`,
`sweling pressure`. For III/IV/VII: header_row = 6 for `mean stress`,
`radial stress`, `axial stress`, `suction`, `porosity`, `permeability`,
`gap closure`, `void ratio`; header_row = 4 for `dry density`; the matching
time column is the nearest `(days)`/`time (days)` to the left.
For Family B: per-quantity sheet, `skiprows=5`, col0=time(d), col1+=value;
unit at original row index 4.
