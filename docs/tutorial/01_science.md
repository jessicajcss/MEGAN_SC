# 1. The science: MEGAN in one chapter

MEGAN (Model of Emissions of Gases and Aerosols from Nature) estimates the flux of volatile compounds from vegetation
and soil as an emission factor under standard conditions, scaled by activity factors that describe how light,
temperature, leaf area, leaf age and soil water move the flux away from that standard (Guenther et al., 2006, 2012).
MEGAN-SC implements in Python the code of MEGANv3.21 (Guenther, 2024b), and `docs/ALGORITHM_SPEC.md` transcribes
every equation and constant from the Fortran source with its file and line. This chapter explains the same content in
the order in which a value is computed, and ends with a worked example you can run.

## 1.1 What is computed

For each grid cell, each time step and each of the 19 emission classes c, the engine computes

```
flux_c   = ER_c × EF_c                                  nmol m⁻² s⁻¹
E_c      = Σ_k  flux_c × f_ck × 10⁻⁹ × A × Δt × M_k × 10⁻⁶      tonnes per cell and step
```

* `ER_c` is the **emission activity**, a dimensionless product of activity factors (section 1.2);
* `EF_c` is the **emission factor** of the cell (section 1.7);
* the class flux is split into the MEGAN compounds k with fractions `f_ck` (201 compounds), each converted to moles
  (10⁻⁹ mol per nmol), multiplied by the cell area `A` (m²) and the length of the time step `Δt` (s), and by the molar
  mass `M_k` (g mol⁻¹), then to tonnes (10⁻⁶ t per g);
* the compounds are summed back to their class, and the steps to months and years.

The MEGAN equation also holds a canopy loss and production factor ρ, E = ε · γ · ρ, where ε is the emission of a
compound into the canopy (Guenther et al., 2006). For isoprene, Guenther et al. (2006, section 3.3 and Table 1)
estimated ρ between 0.93 and 0.99 and recommended the constant value 0.96. The MEGANv3.21 code applies no such factor,
so the files of `data/run_v10` hold isoprene as the code computes it; every script that reads the run multiplies
isoprene by 0.96 (`RHO_ISOP` in notebooks 11, 13, 14 and 15, `RHO` in `analysis/v10/common.py`), and ρ = 1 for the
other classes. Every published number includes the factor.

In the Santa Catarina run, ERA5-Land is sampled at 00, 03, ..., 21 UTC, so each step stands for `Δt` = 3 h, and the
cell area is computed on the authalic sphere (R = 6,371,007.181 m), the sphere whose total surface equals that of the
WGS 84 ellipsoid; at the latitudes of Santa Catarina the spherical cell areas are 0.12 % to 0.19 % larger than the
ellipsoidal ones
(`results/v10/final_review/final_review_numbers.csv`, key `cell_area_sphere_minus_ellipsoid_pct`).

The 19 classes and their coefficients (MEGVEA.EXT; Table 1 of the article):

| # | Class | Representative compounds | β (K⁻¹) | LDF | CT1 | Cleo |
|---|---|---|---|---|---|---|
| 1 | ISOP | isoprene | 0.13 | 1.0 | 95 | 2.00 |
| 2 | MBO | 2-methyl-3-buten-2-ol | 0.13 | 1.0 | 95 | 2.00 |
| 3 | MT_PINE | α- and β-pinene | 0.10 | 0.6 | 80 | 1.83 |
| 4 | MT_ACYC | myrcene, ocimene | 0.10 | 0.9 | 80 | 1.83 |
| 5 | MT_CAMP | carene, camphene | 0.10 | 0.2 | 80 | 1.83 |
| 6 | MT_SABI | sabinene, limonene | 0.10 | 0.4 | 80 | 1.83 |
| 7 | MT_AROM | cymenes, oxygenated monoterpenes | 0.10 | 0.1 | 80 | 1.83 |
| 8 | NO | soil nitric oxide | 0.10 | 0.0 | 80 | 1.83 |
| 9 | SQT_HR | highly reactive sesquiterpenes (caryophyllene) | 0.17 | 0.5 | 130 | 2.37 |
| 10 | SQT_LR | less reactive sesquiterpenes (longifolene) | 0.17 | 0.4 | 130 | 2.37 |
| 11 | MEOH | methanol | 0.08 | 0.8 | 60 | 1.60 |
| 12 | ACTO | acetone | 0.10 | 0.2 | 80 | 1.83 |
| 13 | ETOH | acetaldehyde and ethanol | 0.13 | 0.8 | 95 | 2.00 |
| 14 | ACID | formic, acetic and pyruvic acids | 0.13 | 0.8 | 95 | 2.00 |
| 15 | LVOC | C2-C4 hydrocarbons (ethene and others) | 0.10 | 0.2 | 80 | 1.83 |
| 16 | OXPROD | oxidation products | 0.10 | 0.2 | 80 | 1.83 |
| 17 | STRESS | stress compounds (linalool and others) | 0.10 | 0.8 | 80 | 1.83 |
| 18 | OTHER | other compounds (indole, pentane and others) | 0.10 | 0.1 | 80 | 1.83 |
| 19 | CO | carbon monoxide | 0.08 | 1.0 | 60 | 1.60 |

The documents group them as isoprene (ISOP), monoterpenes (the five MT classes), sesquiterpenes (SQT_HR + SQT_LR),
oxygenated VOC (MEOH + ACTO + ETOH + ACID + OXPROD), light VOC (LVOC), stress and other (STRESS + OTHER) and MBO;
BVOC is the sum of these 17 classes, and soil NO and CO are reported apart.

## 1.2 The emission activity

MEGVEA (megvea.f, lines 758-773) assembles, for each class:

```
ISOP:       ER = LAI · γTP · γCO2 · γage · γSM
ETOH:       ER = LAI · γTP · γBD  · γage · γSM
others:     ER = LAI · γTP        · γage · γSM
```

The four stress factors (high and low temperature, high wind, air quality) multiply every class in the Fortran; they
default to 1 and are off in this run. Soil NO follows a separate algorithm (section 1.8).

## 1.3 Light and temperature: γTP

Each class is a mixture of a **light-dependent** part (fraction LDF), emitted as it is produced in the chloroplast,
and a **light-independent** part (1 − LDF), released from storage pools.

**Light-independent response** (an exponential in leaf temperature T, with β from the table):

```
γTLI(T) = exp[ β · (T − 303.15) ]
```

**Light-dependent temperature response** (T leaf temperature, T24 mean temperature of the past day; megvea.f sets the
240-h mean T240 equal to T24):

```
Topt  = 312.5 + 0.6 · (T240 − 297)
Eopt  = Cleo · exp[0.05 · (T24 − 297)] · exp[0.05 · (T240 − 297)]
x     = (1/Topt − 1/T) / 0.00831
γTLD  = Eopt · 230 · exp(CT1 · x) / [230 − CT1 · (1 − exp(230 · x))]          (0 below 260 K)
```

**Light response** (PPFD at the leaf, µmol m⁻² s⁻¹; the daily mean PPFD only acts as a night gate in MEGANv3.21):

```
γP = 0.004 · 1.03 · PPFD / sqrt(1 + 0.004² · PPFD²)          (0 if the daily mean PPFD < 0.01)
```

**Canopy depth.** The canopy is divided into five layers at the Gaussian nodes 0.0469, 0.2308, 0.5, 0.7692 and
0.9531 of the cumulative LAI, with weights 0.1185, 0.2393, 0.2844, 0.2393, 0.1185. A factor CDEA = −0.2 · min(LAI ·
depth, 3) + 1.3 raises the light-dependent emission of the upper layers and lowers it in the lower ones.

**Canopy integration.** With SunT, ShaT (leaf temperature of sunlit and shaded leaves), SunP, ShaP (their PPFD) and
SunF (sunlit fraction) for each layer k:

```
Ea1_k = CDEA_k · γTLD(SunT_k) · γP(SunP_k) · SunF_k  +  γTLD(ShaT_k) · γP(ShaP_k) · (1 − SunF_k)
Ea2_k = γTLI(SunT_k) · SunF_k  +  γTLI(ShaT_k) · (1 − SunF_k)
γTP   = Σ_k  [ Ea1_k · LDF + Ea2_k · (1 − LDF) ] · weight_k
```

CDEA multiplies only the sunlit term. This is how megvea.f evaluates its expression under Fortran operator precedence,
and MEGAN-SC evaluates it the same way through procedure FIX-6 (section 1.9).

## 1.4 The canopy model (MEGCAN)

MEGCAN turns the meteorology above the canopy into the leaf temperatures and light of section 1.3, for each of six
canopy types weighted by their fractions in the cell (`src/megan_sc/engine/canopy.py`, from megcan.f):

1. **Solar geometry**: solar elevation from the day of year, latitude and solar hour (UTC hour + longitude/15).
2. **Radiation units**: PAR = 0.45 × downward short-wave (W m⁻²); PPFD = 4.5 × PAR (µmol m⁻² s⁻¹); the total solar
   radiation entering the leaf energy balance is PPFD/2.25, that is 0.9 × short-wave.
3. **Beam and diffuse fractions** of visible and near-infrared light from the ratio of the observed to the top-of-
   atmosphere radiation.
4. **Radiation through the canopy**: extinction coefficients for beam and diffuse light, sunlit fraction
   exp(−Kb · cumulative LAI) at each layer, and the PPFD of sunlit and shaded leaves.
5. **Leaf energy balance**: for each layer, the leaf temperature that balances absorbed radiation, sensible heat,
   latent heat and long-wave exchange, solved with ten fixed iterations and limited to ±10 K from the air temperature.

At night (sun elevation sine ≤ 0.002, or no light) the sunlit fraction is set to 0.2 and the PPFD to zero, so only the
light-independent part emits.

## 1.5 Leaf age: γage

New leaves emit little isoprene and many monoterpenes and much methanol; mature leaves emit the reverse. MEGAN splits
the foliage into new, growing, mature and old fractions from the change of LAI between the previous period (LAIp)
and the current one (LAIc):

```
LAIp < LAIc (growth):   ti = 5 + 0.7·(300 − T24) if T24 ≤ 303 K, else 2.9;   tm = 2.3 · ti
                        Fnew = (1 − LAIp/LAIc) · min(1, ti/t)     Fmat from LAIp/LAIc and tm;   Fold = 0
LAIp = LAIc:            Fnew = 0, Fgro = 0.1, Fmat = 0.8, Fold = 0.1
LAIp > LAIc (loss):     Fnew = 0, Fgro = 0, Fold = (LAIp − LAIc)/LAIp, Fmat = 1 − Fold
γage = Fnew·Anew + Fgro·Agro + Fmat·Amat + Fold·Aold
```

with the class coefficients Anew, Agro, Amat, Aold of `docs/ALGORITHM_SPEC.md` (for isoprene 0.05, 0.6, 1.0, 0.9; for
monoterpenes 2.0, 1.8, 1.0, 1.05) and `t` the length of the LAI period (30 days for monthly LAI, NLAI = 12). With
LAIp = LAIc, γage of isoprene is 0.1 × 0.6 + 0.8 × 1.0 + 0.1 × 0.9 = 0.95. Passing the current LAI as LAIp would keep
γage at this steady value; as MEGANv3.21 reads the previous LAI period, the pipeline passes the previous month's LAI as
LAIp (FIX-3).

## 1.6 Soil moisture, CO₂ and bidirectional exchange

**Soil moisture** (megsea.f): with θ the volumetric water of the top soil layer and θw the wilting point,

```
γSM = 0                    if θ < θw
γSM = (θ − θw) / 0.04      if θw ≤ θ < θw + 0.04
γSM = 1                    otherwise
```

The Fortran takes θw from the soil type; this run used a constant θw = 0.12 m³ m⁻³, so γSM < 1 only below
0.16 m³ m⁻³. In the fix logs of the run, the monthly mean γSM was below 1 in 0.63 % of the cell-months (report key
`gamsm_lt1`), and the state mean stayed between 0.992 and 1.000 (test S3).

**CO₂ inhibition of isoprene** (megvea.f): with the intercellular concentration Ci = 0.7 × CO₂ (ppm),

```
γCO2 = 1.344 − 1.344 · Ci^1.4614 / (585^1.4614 + Ci^1.4614)
```

MEGANv3.21 fixes CO₂ at 400 ppm, where the code sets γCO2 = 1, and this run keeps that setting, so the rise of CO₂
over 2003-2025 does not enter the isoprene trend. For scale, the same formula gives 1.025 at 376 ppm and 0.979 at
425 ppm (computed with `megan_sc.engine.gamma.gamma_co2`; Heald et al., 2009, discuss this response).

**Bidirectional exchange** (class ETOH only): γBD = 0.5 · LAI for LAI < 2; 1 − 0.0625 · (LAI − 2) for 2 ≤ LAI ≤ 6;
0.75 above.

## 1.7 Emission factors and land cover

MEGANv3.21 builds the emission factor of a cell from its **ecotype** (which sets the plant species and their
emission factors) and its **growth forms** (tree, shrub, herb, crop):

```
EF_cell,c = Σ_g  fraction_g × EF_ecotype,g,c
```

In MEGAN-SC the growth-form fractions come from MapBiomas Collection 10 at 30 m (Souza Junior et al., 2020),
aggregated to the 0.025° grid through a crosswalk (`src/megan_sc/landcover.py`): forest formation, mangrove,
floodable forest and wooded restinga to tree; savanna 50 % tree and 50 % shrub; wetland, grassland, herbaceous restinga
and pasture to herb; farming and mosaic of uses 50 % herb and 50 % crop; temporary and perennial crops to crop;
silviculture to tree (no separation of *Pinus* and *Eucalyptus*); palm oil to crop; water, urban, beach, mining, rock
and other non-vegetated classes carry no vegetation. The ecotype map and the emission-factor and speciation tables come from the MEGAN 3.2 input files
(Guenther, 2024a): `data/sc_ecotype.nc` is a clip of `EVT3d.nc`, and `data/megan_inputs/EFP240120/` holds the tables. The canopy type of trees is temperate broadleaf (MEGCAN type 3). MapBiomas ends in
2024, so 2025 uses the 2024 land cover.

LAI is the monthly mean of the good-quality MODIS MCD15A3H retrievals (Myneni; Knyazikhin; Park, 2021) aggregated to
the grid; land cover is used to mask it (no LAI where there is no vegetation) but never to invent LAI.

## 1.8 Soil NO

Soil NO follows Yienger and Levy (1995) as coded in MEGANv3.21 (soilnox.F, megsea.f), on the path that uses air
temperature only (no soil temperature, soil moisture or rain pulses):

```
T      = min(T_air, 303 K)                        (FIX-7: the cap applies to both branches)
grass  CFNOG = exp(0.04686 · T − 14.30579)        for T > 268.869 K
Tsoil  = 0.72 · T + 82.28, limited to 273.16-303.16 K
dry    CFNODRY = (Tsoil − 273.16) / 90
wet    CFNOWET = (Tsoil − 273.16) · e^(−3.09) · 0.28      if Tsoil ≤ 283.16 K
               = exp(0.103 · (Tsoil − 273.16)) · e^(−3.09)  otherwise
crop   CFNO = 0.5 · (CFNOWET + CFNODRY) × fertiliser factor × canopy reduction
canopy reduction = [exp(−0.24 · LAI) + exp(−0.0525 · LAI)] / 2
```

During the growing season the soil-NO activity is the canopy-type-weighted mean of CFNO (crops) and CFNOG (all other
types); outside it, CFNOG everywhere. The fertiliser factor depends on the day of the growing season. South of 23° S,
MEGANv3.21 compares the day of the year with month-day constants (1101, 0531), which makes its season degenerate and
the fertiliser factor negative in winter; MEGAN-SC uses the intended November-May season and clips the factor at zero,
a documented deviation.

## 1.9 The input and run procedures FIX-1 to FIX-8

Eight procedures of `src/megan_sc/fixes.py` prepare the inputs and run the engine so that it follows MEGANv3.21 and the
ERA5-Land conventions. Notebook 10 writes this file, and its SHA-256 prefix a933b29602929a82 is stored in every output
file and checked by notebook 11 (check C0).

| Id | What the procedure does | Why it is needed |
|---|---|---|
| FIX-1 | every band of a cell without valid meteorology is taken from the nearest valid cell (all bands from the same cell) before the engine runs; the run stops if a vegetated cell keeps invalid meteorology | ERA5-Land covers land only, so some cells on the coast, on lagoons and on islands have no data; the engine alone would turn the fill value into T = 340 K and PAR = 0, which silences the light-dependent classes and multiplies the light-independent ones by exp[β(340 − 303.15)] (about 40 for monoterpenes, 525 for sesquiterpenes) |
| FIX-2 | the daily means of temperature and PPFD are the 24-h means of the simulation day | MEGANv3.21 computes daily means (daymet.f) and uses them in the light-dependent temperature response and as the night gate of the light response (section 1.3); given no daily means, the engine alone would approximate them by T − 2 K and 0.5 × PPFD of the current hour |
| FIX-3 | the previous month's LAI is passed as LAIp | MEGANv3.21 computes leaf age from the change of LAI between the previous period and the current one (section 1.5) |
| FIX-4 | a vegetated cell without a good-quality MODIS retrieval in a month takes its LAI from the adjacent months | some cell-months have no good-quality MCD15A3H retrieval, and LAI = 0 in a vegetated cell would give it zero emission |
| FIX-5 | ecotype no-data pixels take the nearest valid ecotype | the MEGAN ecotype map (`EVT3d.nc`) has no-data pixels on a coastal strip, which would give EF = 0 to the vegetated cells they cover |
| FIX-6 | γTP with CDEA on the sunlit term only | megvea.f applies the canopy-depth factor to the sunlit term only (section 1.3), while the engine module `gamma.gamtp` alone applies it to both terms |
| FIX-7 | soil NO with the 303 K cap on air temperature before both branches | soilnox.F caps the air temperature at 303 K before both the grass factor and the soil-temperature proxy (section 1.8), while the engine module `soil.soilnox` alone caps the grass branch only |
| FIX-8 | the hourly radiation flux is the accumulation minus that of the previous hour, at every hour except 01 UTC, where the accumulation restarts | ERA5-Land radiation is accumulated from 00 UTC and the field valid at 00 UTC holds the previous 24 h (ECMWF, 2025); read as an hourly value, the 00 UTC sample would carry a daily total (about 2,700 W m⁻² of PAR at night) |

Notebook 12 compares the engine with the official Fortran executables (`VALIDATION.md`): with FIX-6 and FIX-7, every
stage of the engine agrees with the Fortran to float32 precision, and controls run with the engine modules alone show
that the comparison detects a departure from megvea.f and soilnox.F.

## 1.10 Settings of the Santa Catarina run

* Grid 0.025° (about 2.8 km), 144 × 224 cells over the bounding box of the state; a cell is in the state when its
  centre lies inside the IBGE polygon (13,970 cells). Municipal and mesoregion totals split the emission of each cell
  among the polygons in proportion to the area of the cell inside each one.
* ERA5-Land every 3 h (eight steps per day), 2003-2025; MODIS LAI monthly; MapBiomas annual (2025 uses 2024).
* Wilting point 0.12 m³ m⁻³; γSM on; γBD on; stress factors and CO₂ inhibition off.
* Light-dependent fraction of MT_PINE, MT_ACYC, MT_CAMP and MT_SABI from the class constants, where the Fortran reads a
  gridded map.
* Soil NO on the air-temperature path, growing season November-May south of 23° S.
* Isoprene multiplied by ρ = 0.96 when the run is read (section 1.1).
* Trees as temperate broadleaf canopy; silviculture mapped to trees.

These choices are stated in both documents, and test S3 (chapter 5) measures the effect of γSM.

## 1.11 Worked example

The activity factors are plain functions of `src/megan_sc/engine/gamma.py`. From the repository root:

```python
import sys; sys.path.insert(0, 'src')
from megan_sc.engine import gamma as g, species as sp
isop, pine, sqt = sp.index_of('ISOP'), sp.index_of('MT_PINE'), sp.index_of('SQT_HR')
print(g.gamtld(303.15, 297.0, isop))     # light-dependent T response of isoprene, leaf at 30 °C, daily mean 24 °C
print(g.gamp(1500.0, 400.0))             # light response at PPFD 1500 µmol m-2 s-1
print(g.gamtli(340.0, pine), g.gamtli(340.0, sqt))   # what the engine would do with a fill value without FIX-1
print(g.gamma_a(3.0, 3.0, 297.0, isop, 12))          # leaf age with LAIp = LAIc
print(g.gamma_co2(376.0), g.gamma_co2(425.0))        # CO2 inhibition (off in the run)
```

Output, rounded:

```
1.0543        γTLD of isoprene for a leaf at 30 °C after a day at 24 °C
1.0160        γP close to its plateau at full sun
39.8  525.6   at 340 K the light-independent emission of α-pinene and caryophyllene is multiplied by 40 and 526
0.95          steady-state leaf age of isoprene
1.0250  0.9794
```

The third line shows what the engine would do with a fill value without FIX-1: a cell run at 340 K with no light would
emit no isoprene and forty times the monoterpenes it should have. This is why FIX-1 fills such cells and why notebook 11
looks for forested cells without isoprene (check C1).
