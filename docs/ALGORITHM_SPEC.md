# MEGANv3.21: Verified Algorithm Specification (as implemented in *this* code)

**Purpose.** This document is transcribed directly from the `MEGANv3.21` Fortran source (Guenther, 2024b), not from papers or memory. It is the reference the Python implementation is built and validated against. Every constant and equation below cites the source file it came from; sections marked **[VERIFIED]** were read line by line.

Source location in workspace: `reference/src/…`

**Primary source (archived):** Guenther, A. (2024b). *MEGANv3.21.tar.gz emission calculator code* (3.21). Zenodo. https://doi.org/10.5281/zenodo.10526206. This port is transcribed from that exact code.

**Input files (archived):** Guenther, A. (2024a). *MEGAN 3.2 input files* [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.10939297. The emission-factor and speciation tables (`EFP240120`), the ecotype map (`EVT3d`) and the growth-form maps (`GFv7nc`) named below are files of this record.

---

## 0. Pipeline & data flow  [VERIFIED from work/step*.csh]

Eight programs run in sequence:

| Step | Program | Role | Key output |
|---|---|---|---|
| 1 | `txt2ioapi` | convert text input maps → IOAPI netCDF | EF map, LAI, canopy type, ecotype, LDF, … |
| 2 | `met2mgn` | convert WRF/MCIP meteorology → MEGAN met | MET.MEGAN…ncf |
| 3 | `daymet` | daily-average meteorology | D_TEMP, D_PPFD, MaxT, MinT, MaxWS |
| 4 | `megcan` | canopy radiation + leaf energy balance | SunT, ShaT, SunP, ShaP, SunF (per layer) |
| 5 | `megsea` | soil moisture activity + soil NO | GAMSM, GAMNO |
| 6 | `megvea` | **vegetation emission activity (all γ factors)** | ER per class (non-dimensional) |
| 7 | `mgn2mech` | apply EF, speciate, convert units | mol s⁻¹ or tonne hr⁻¹ |
| 8 | `ioapi2uam` | (optional) reformat for CAMx | |

For the Python port, steps 4–7 are the science; steps 1–2–8 are I/O plumbing we replace with the Earth-Engine + xarray layer.

---

## 1. Emission classes & coefficients  [VERIFIED from MEGVEA.EXT]

`NCLASS = 19`. Class index → name and per-class coefficients:

| # | MGN_SPC | β (light-indep. T) | LDF | CT1 | Cleo | Anew | Agro | Amat | Aold |
|---|---|---|---|---|---|---|---|---|---|
| 1 | ISOP (isoprene) | 0.13 | 1.0 | 95 | 2.00 | 0.05 | 0.6 | 1.0 | 0.9 |
| 2 | MBO | 0.13 | 1.0 | 95 | 2.00 | 0.05 | 0.6 | 1.0 | 0.9 |
| 3 | MT_PINE (α/β-pinene) | 0.10 | 0.6* | 80 | 1.83 | 2.0 | 1.8 | 1.0 | 1.05 |
| 4 | MT_ACYC (myrcene, ocimene) | 0.10 | 0.9* | 80 | 1.83 | 2.0 | 1.8 | 1.0 | 1.05 |
| 5 | MT_CAMP (carene, camphene) | 0.10 | 0.2* | 80 | 1.83 | 2.0 | 1.8 | 1.0 | 1.05 |
| 6 | MT_SABI (sabinene, limonene) | 0.10 | 0.4* | 80 | 1.83 | 2.0 | 1.8 | 1.0 | 1.05 |
| 7 | MT_AROM (cymenes + oxygenated) | 0.10 | 0.1 | 80 | 1.83 | 2.0 | 1.8 | 1.0 | 1.05 |
| 8 | NO (soil nitric oxide) | 0.10 | 0.0 | 80 | 1.83 | 1.0 | 1.0 | 1.0 | 1.0 |
| 9 | SQT_HR (e.g. caryophyllene) | 0.17 | 0.5 | 130 | 2.37 | 0.4 | 0.6 | 1.0 | 0.95 |
| 10 | SQT_LR (e.g. longifolene) | 0.17 | 0.4 | 130 | 2.37 | 0.4 | 0.6 | 1.0 | 0.95 |
| 11 | MEOH (methanol) | 0.08 | 0.8 | 60 | 1.60 | 3.5 | 3.0 | 1.0 | 1.2 |
| 12 | ACTO (acetone) | 0.10 | 0.2 | 80 | 1.83 | 1.0 | 1.0 | 1.0 | 1.0 |
| 13 | ETOH (acetaldehyde+ethanol) | 0.13 | 0.8 | 95 | 2.00 | 1.0 | 1.0 | 1.0 | 1.0 |
| 14 | ACID (formic/acetic/pyruvic) | 0.13 | 0.8 | 95 | 2.00 | 1.0 | 1.0 | 1.0 | 1.0 |
| 15 | LVOC (C2–C4 HC, e.g. ethene) | 0.10 | 0.2 | 80 | 1.83 | 1.0 | 1.0 | 1.0 | 1.0 |
| 16 | OXPROD (oxidation products) | 0.10 | 0.2 | 80 | 1.83 | 1.0 | 1.0 | 1.0 | 1.0 |
| 17 | STRESS (e.g. linalool) | 0.10 | 0.8 | 80 | 1.83 | 1.0 | 1.0 | 1.0 | 1.0 |
| 18 | OTHER (indole, pentane, …) | 0.10 | 0.1 | 80 | 1.83 | 1.0 | 1.0 | 1.0 | 1.0 |
| 19 | CO (carbon monoxide) | 0.08 | 1.0 | 60 | 1.60 | 1.0 | 1.0 | 1.0 | 1.0 |

\* For classes **3,4,5,6** the LDF is **not** the constant above but is read from a gridded `LDFILE` map (`megvea.f` lines 651–666). For all other classes LDF is the constant in this table.

Other constants (MEGVEA.EXT): `CO2 = 400`, `ISmax = 1.344`, `CO2h = 1.4614`, `Cstar = 585`, `PSTD = 200`, canopy-depth `CCD1 = -0.2`, `CCD2 = 1.3`. Stress coefficients (CAQ/CHW/CHT/CLT) and thresholds (TAQ/THW/THT/TLT + deltas) per class are tabulated in MEGVEA.EXT; **all stress factors default OFF** in the run scripts (see §4).

---

## 2. Activity-factor functions  [VERIFIED from megvea.f]

**Light-dependent temperature response** `GAMTLD(T1, T24, S)` (leaf temp T1; daily-avg T24). Note: this code sets `T240 = T24` (no separate 240-h history):
```
if T1 < 260 K:  GAMTLD = 0
Topt   = 312.5 + 0.6*(T240 - 297)
X      = (1/Topt - 1/T1) / 0.00831
Eopt   = Cleo(S) * exp(0.05*(T24-297)) * exp(0.05*(T240-297))
GAMTLD = Eopt * 230 * exp(CT1(S)*X) / (230 - CT1(S)*(1 - exp(230*X)))     # Ct2=230
```

**Light-independent temperature response** `GAMTLI(T, S)`:  `exp(β(S)·(T − 303.15))`

**Light response** `GAMP(PPFD1, PPFD24)` (leaf PPFD; daily-avg PPFD only used as a gate):
```
if PPFD24 < 0.01:  GAMP = 0
else: GAMP = (0.004 * 1.03 * PPFD1) / sqrt(1 + 0.004² · PPFD1²)     # Alpha=0.004, C1=1.03 (constant)
```
> Note: the PPFD24-dependent `C1` form is commented out in this code; C1 is fixed at 1.03. The port must replicate this exactly.

**Canopy-depth response** `GAMMA_CD` (per layer): `CDEA = -0.2·min(LAI·Cdepth_k, 3) + 1.3`, with 5-layer nodes `Cdepth = [0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899]`.

**Canopy integration** `GAMTP`: 5-layer Gaussian quadrature, weights `VPGWT = [0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635]`:
```
Ea1L_k = CDEA_k·GAMTLD(SunT_k)·GAMP(SunP_k)·SunF_k + GAMTLD(ShaT_k)·GAMP(ShaP_k)·(1-SunF_k)     # light-dependent (CDEA on the sunlit term only; see section 11)
Ea2L_k = GAMTLI(SunT_k)·SunF_k + GAMTLI(ShaT_k)·(1-SunF_k)                                       # light-independent
GAMTP  = Σ_k [ (Ea1L_k·LDF + Ea2L_k·(1-LDF)) · VPGWT_k ]
```

**Leaf-age response** `GAMMA_A` (uses current/previous LAI and daily temp; timestep `t`=8 days for 46-period LAI, 10 for 37, 30 for 12):
```
if LAIp < LAIc (growth):
    ti = 5.0 + 0.7·(300 - Tt) if Tt≤303 else 2.9 ;  tm = 2.3·ti
    Fnew = (1 - LAIp/LAIc) · (1 if ti≥t else ti/t)
    Fmat = (LAIp/LAIc)     · (1 if tm≥t else …) ;  Fgro = 1-Fnew-Fmat ;  Fold = 0
elif LAIp == LAIc:  Fnew=0, Fgro=0.1, Fmat=0.8, Fold=0.1
else (senescence):  Fnew=0, Fgro=0, Fold=(LAIp-LAIc)/LAIp, Fmat=1-Fold
GAMLA = Fnew·Anew(S) + Fgro·Agro(S) + Fmat·Amat(S) + Fold·Aold(S)
```

**CO₂ inhibition** `GAMMA_CO2` (isoprene only): with fixed CO2=400 → `GAMCO2 = 1`; otherwise `Ci = 0.7·CO2`, `GAMCO2 = ISmax − ISmax·Ci^1.4614 / (585^1.4614 + Ci^1.4614)`.

**Bidirectional LAI** `GAMMA_LAIbidir` (class 13 only): 0.5·LAI (LAI<2); 1−0.0625·(LAI−2) (2≤LAI≤6); 0.75 (LAI>6).

**Stress** (`GAMMA_HT/LT/HW/AQ`): piecewise-linear ramps from 1 to the class C-coefficient across the delta threshold. Default OFF.

---

## 3. Final emission-activity assembly  [VERIFIED from megvea.f lines 758–773]

Per class S, per grid cell (stress terms = 1 when disabled):
```
S = 1  (isoprene):  ER = LAIc · GAMTP · GAMCO2 · GAMLA · GAMHW·GAMAQ·GAMHT·GAMLT · GAMSM
S = 13 (etoh/acet): ER = LAIc · GAMTP · GAMBD  · GAMLA · GAMHW·GAMAQ·GAMHT·GAMLT · GAMSM
else:               ER = LAIc · GAMTP           · GAMLA · GAMHW·GAMAQ·GAMHT·GAMLT · GAMSM
ER = max(ER, 0)
```
This `ER` is the **non-dimensional emission activity** written by megvea (one field per class).

---

## 4. Run-time switches  [VERIFIED from work/step6.run.megvea]

`GAMVEA` reads YN flags (default `.FALSE.`): `GAMAQ_YN, GAMCO2_YN, GAMHW_YN, GAMHT_YN, GAMLT_YN, GAMSM_YN, GAMBD_YN`; and `Layers` (default 5), `NLAI` (46/37/12), `N_MaxT/N_MinT/N_MaxWS` (default 5). The as-shipped SE-US test enables `GAMSM_YN` (soil moisture) and `GAMBD_YN`; CO2 and the storm/temperature stresses are typically off. The settings of the Santa Catarina run are listed in section 11.

---

## 5. Emission factors & unit conversion  [VERIFIED from mgn2mech.F]

- EF is a gridded map per class, variable `EF_<MGN_SPC>`, in **nanomol m⁻² s⁻¹** (MEGAN3 changed from µg m⁻² h⁻¹; comment lines 55–56).
- Speciation: 19 classes → 201 detailed species (`SPCA`) via `mg20_map`/`spca_map` and fractional factors `effs_all`. Molar weights `spca_mwt` per species; mechanism weights `mech_mwt`.
- Constants: `nmol2mol = 1e-9`, `hr2sec = 3600`, `g2tonne = 1e-6`, `garea = xcell·ycell` (m²).

Per-species emission build (non-NO): `base = ER · EF · effs_all`  [nanomol m⁻² s⁻¹]. NO uses `GAMNO·EF·effs_all`.
```
mol_m2_hr   = base · 1e-9 · 3600
mol_s_cell  = mol_m2_hr · garea / 3600          # "native" per-cell rate (mol s⁻¹)
tonne_hr    = mol_s_cell · MW · 1e-6 · 3600     # tonnes per hour (MW = g/mol)
```
Monthly tonnes = Σ over the month of `tonne_hr · Δt_hours`, per species per cell. **Our geodesic cell area replaces `garea` (which assumes a constant-area projected grid).**

---

## 6. Inputs required by the engine

| Field | Produced by | Our source |
|---|---|---|
| SunT, ShaT, SunP, ShaP, SunF (per layer) | megcan | port megcan; drive with ERA5 T, PPFD, wind, humidity, pressure + solar geometry |
| D_TEMP, D_PPFD, MaxT, MinT, MaxWS | daymet | daily aggregation of ERA5 in xarray |
| LAIc, LAIp | LAI map | MODIS MCD15A3H via Earth Engine |
| GAMSM, GAMNO | megsea | port megsea; ERA5 soil water + wilting point |
| EF_* (19), LDF map, canopy type, ecotype | txt2ioapi from global data | provided `EFP240120` / `EVT3d` / `GFv7nc` (MEGAN 3.2 input files, Guenther, 2024a), clipped to SC |

---

## 7. Verification status

- **[VERIFIED + PORTED + TESTED]** MEGVEA.EXT + megvea.f: 19-class coefficients + every γ function → `engine/species.py`, `engine/gamma.py`; the ER assembly loop → `engine/megvea.py`; megcan.f + MEGCAN.EXT (canopy radiation + leaf energy balance) → `engine/canopy.py`; megsea.f soil-moisture block + MEGSEA.EXT → `engine/soil.py`; daymet.f → `engine/daymet.py`; mgn2mech.F + SPC_NOCONVER.EXT: 201-species speciation, molar weights, EF application, unit chain → `engine/speciation.py` (+ auto-generated `engine/speciation_data.py`), with unit tests that include an end-to-end met→canopy→ER→speciation→tonnes integration test (`tests/`).
- **[VERIFIED + PORTED + TESTED]** the soil-NO Yienger-Levy 1995 routine (`soilnox.F` + `soilnox_fx.F` + megsea GAMNO assembly) → `engine/soil.py` (`soilnox`, `growseason`, `fertilizer_adj`, `veg_adj`, `gamma_no`), wired into `run.run_timestep`, with its unit tests (`tests/test_soilno.py`).

**Soil-NO growing-season decision (SH):** the Fortran `GROWSEASON` compares day-of-year against month-day literals, so in the Southern Hemisphere the growing season degenerates and the fertilizer factor goes negative for much of the year (unphysical). We DEFAULT to the corrected SH season (Nov–May growing) and clamp the fertilizer factor ≥ 0; `corrected=False` reproduces the exact source behaviour for bit-comparison. This is the one place we deviate from the source by default, for physical robustness, and it is flagged here and in the code.

**Source quirk preserved (not corrected):** SPC_NOCONVER.EXT duplicates the `cadinene_d` map line, leaving `selinene_b` (species 67) unmapped → it emits 0, and class speciation fractions sum to OXPROD 0.995 / OTHER 1.005. This is MEGANv3.21's actual behaviour and is reproduced exactly.

---

## 8. Canopy model: MEGCAN  [VERIFIED, ported to engine/canopy.py]

Per grid cell and canopy type (6 types, area-weighted by canopy-type fractions), for 5 Gaussian layers:

- **Solar geometry**: `SinDelta = -sin(0.40907)·cos(6.28(Day+10)/365)`; solar elevation `Calcbeta`; eccentricity `1 + 0.033·cos(2·3.14(Day-10)/365)`. Solar hour = UTC hour + lon/15.
- **PPFD/Solar**: `PPFD[µmol m⁻² s⁻¹] = PAR[W m⁻²]·4.5`; `Solar[W m⁻²] = PPFD/2.25`; `Maxsolar = sinβ·1367·eccentricity`.
- **Solar fractions** (Spitters/Goudriaan): transmission `Solar/Maxsolar` → diffuse/beam split of visible & near-IR (`FracDiff = 0.156 + 0.86/(1+e^{11.1(τ-0.53)})`, `PPFDfrac = 0.55 - 0.12τ`).
- **Two-stream radiation** (`CanopyRad`): `Kb = cluster·0.5/sinβ`, `Kd = 0.8·cluster`; per layer `Sunfrac = e^{-Kb·LAIdepth}`; sun/shade PPFD from absorbed diffuse+scattered+beam, with `ConvertShadePPFD = 4.6`, `ConvertSunPPFD = 4.0`. Night (sinβ≤0.002 or no light): Sunfrac = 0.2, PPFD = 0.
- **Leaf energy balance** (`CanopyEB`/`LeafEB`): per-layer air-temp profile `TairK0 + Trate·Ldepth` (lapse rate from `Stability`), humidity profile, wind profile; leaf temperature solved by a fixed 10-iteration Newton balance of absorbed radiation, sensible heat (`LeafH`), latent heat (`LeafLE`, `SvdTk`, `LHV`), and IR (`LeafIR`, ε_atm = 0.642·(e/T)^{1/7}); `Tleaf = TairK + Δ`, Δ clamped to ±10 K. 6 canopy-type characteristics from the `Canopychar(17,6)` table in MEGCAN.EXT.

## 9. Soil moisture & daily met  [VERIFIED, ported]

- **GAMSM** (`engine/soil.py`): `0` if θ<wilt; `(θ-wilt)/0.04` if wilt≤θ<wilt+0.04; else `1`. Wilting point `WWLT(16)` by soil type (MEGSEA.EXT).
- **Daily met** (`engine/daymet.py`): `D_TEMP = mean(T)`, `D_PPFD = mean(PAR[W m⁻²])·4.5`, `MaxT/MinT = max/min(T)`, `MaxWS = max(wind)` over the 24 hourly steps.

---

## 10. MapBiomas integration (high-resolution vegetation)

A user-selectable option (`datasets.vegetation.growth_form_source: mapbiomas`) uses **MapBiomas Brazil** (30 m, annual, `projects/mapbiomas-public/assets/brazil/lulc/v1`) for the vegetation structure. It is wired in *coherently*, without displacing MEGAN's calibrated science:

**What MapBiomas provides** (`landcover.py` crosswalk → `gee/mapbiomas.py`):
- **Growth-form fractions** (tree/shrub/herb/crop) per analysis cell, from a documented class→growth-form crosswalk aggregated from 30 m. These replace the global GFv7nc growth-form file and are Brazil-specific.
- **Canopy-type fractions** for `megcan` (6 MEGAN types), from growth form + a forest-type assumption (temperate broadleaf for subtropical SC, overridable).
- **A 30 m vegetation mask** (water/urban/bare → 0), sharper than the global ecotype's coastline.

**What stays MEGAN's, to protect data quality:**
- The **ecotype (EVT) → species → emission-factor** composition is unchanged; MapBiomas growth-form fractions only *weight* it: `EF_cell = Σ_growthform frac_MapBiomas · EF_ecotype[growthform]` (`run.build_ef_grid_multiform`). This is exactly MEGAN's EF-processor structure with a better growth-form source.
- **LAI stays the observed MODIS field**: land cover is not a leaf-area measurement, so MapBiomas is used to *mask/sharpen* LAI (`lai_mask_with_mapbiomas`), never to invent LAI magnitudes. A per-class LAI-climatology option exists but is off by default (it would substitute assumption for observation).

**Documented assumptions** (`landcover.CROSSWALK_NOTES`, all overridable): Savanna→50/50 tree/shrub; Mosaic/Farming→50/50 herb/crop; Silviculture→broadleaf tree; Palm→crop; forest→temperate-broadleaf canopy for SC. Verified against the MapBiomas Collection 9/10 legend.

Coverage: `landcover.py` (crosswalk + aggregation), `gee/mapbiomas.py` (server-side remap + reduceResolution + export), `run.build_ef_grid_multiform`. 5 crosswalk unit tests.


## 11. Input and run procedures of `fixes.py` and the settings of the Santa Catarina run

**GAMTP, as megvea.f computes it.** In megvea.f (lines 725-731) the canopy-depth factor `CDEA_k` multiplies only the
sunlit term of `Ea1L_k` (Fortran operator precedence), as written in section 2:

```
Ea1L_k = CDEA_k·GAMTLD(SunT_k)·GAMP(SunP_k)·SunF_k + GAMTLD(ShaT_k)·GAMP(ShaP_k)·(1-SunF_k)
```

The engine module alone (`engine/gamma.py`, `gamtp`) multiplies both the sunlit and the shaded terms by `CDEA_k`;
`fixes.gamtp_megan321` implements the form above and replaces it at run time (FIX-6). Against the Fortran it agrees to a
maximum relative difference of 4 × 10⁻⁶ (`../VALIDATION.md`).

**The eight procedures** (`src/megan_sc/fixes.py`; the SHA-256 of the file, prefix a933b29602929a82, is stamped in
every output file and checked by C0 of notebook 11). FIX-1 to FIX-5 and FIX-8 prepare the inputs and drive the engine;
FIX-6 and FIX-7 are the two engine patches, and `fixes.apply_engine_patches()` installs them together with FIX-8:

| Id | Where | Procedure |
|---|---|---|
| FIX-1 | driver | ERA5-Land covers land only, so some coastal cells have no data, and the engine alone would turn the fill value into T = 340 K and PAR = 0; each such cell takes the meteorology of the nearest valid cell (all bands from the same donor) before the engine runs, and the run stops if a vegetated cell keeps invalid meteorology |
| FIX-2 | driver | D_TEMP and D_PPFD are the 24-h means of the simulation day (daymet.f), not T(hour) − 2 K and 0.5 · PPFD(hour), the defaults of the engine alone |
| FIX-3 | driver | leaf age with the previous period's LAI (NLAI = 12, TSTLEN = 30 d) |
| FIX-4 | data | monthly MODIS LAI gaps in vegetated cells filled from the adjacent months |
| FIX-5 | data | ecotype no-data pixels take the nearest valid ecotype |
| FIX-6 | engine | GAMTP as above |
| FIX-7 | engine | soil NO, no-soil path: air temperature capped at 303 K before both the grass factor and TSOI = 0.72·TAIR + 82.28 (soilnox.F) |
| FIX-8 | data | ERA5-Land radiation: hourly flux = difference with the previous hour at every hour except 01 UTC, where the accumulation restarts; the field valid at 00 UTC holds the 24 h ending at 00 UTC (ECMWF ERA5-Land documentation) |

**Settings of the Santa Catarina production run** (stated in the manuscript):
* ERA5-Land sampled every 3 h (00, 03, ..., 21 UTC); each step stands for 3 h in the monthly integration.
* Wilting point constant at 0.12 m³ m⁻³ for GAMSM (section 9 describes the soil-type table of the Fortran, which the
  run did not use); GAMSM on; GAMBD on for ethanol/acetaldehyde.
* Stress factors (GAMHT, GAMLT, GAMHW, GAMAQ) off, so MaxT, MinT and MaxWS do not affect the result; CO₂ inhibition off
  (γCO2 = 1, 400 ppm).
* Light-dependent fraction of MT_PINE, MT_ACYC, MT_CAMP and MT_SABI from the class constants (asterisked in section 1)
  instead of the gridded LDF map read by the Fortran.
* Soil NO: Yienger-Levy air-temperature path (no soil type), growing season November-May south of 23° S (section 7).
* Canopy loss and production factor: the MEGANv3.21 code applies none, and the engine writes isoprene as the code
  computes it; the analysis multiplies isoprene by ρ = 0.96 (Guenther et al., 2006, section 3.3 and Table 1) when it
  reads the run (`RHO_ISOP` in notebooks 11, 13, 14 and 15; `RHO` in `analysis/v10/common.py`).
* Canopy type for trees: temperate broadleaf (MEGCAN type 3); silviculture mapped to the tree growth form.
