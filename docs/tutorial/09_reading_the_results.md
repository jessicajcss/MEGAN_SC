# 9. Reading the results

A biogenic inventory is a model estimate: the emission factors of the MEGAN tables, scaled by activity factors computed
from the meteorology, the leaf area and the land cover that the model read (chapter 1). This chapter explains what each
result of a run means, how to read it, and what it can and cannot support. The examples come from the Santa Catarina
run (`results/v10/`); notebook 25 writes the same kinds of files for any run.

## 9.1 Files, units and conventions

| File (notebook 25, `<run>/results/`) | Content | Unit |
|---|---|---|
| `monthly_totals_by_class.csv` | total of the region for every month, 19 classes | t month⁻¹ |
| `annual_totals_by_class.csv`, `annual_totals_by_group.csv` | totals of the region for every complete year | t yr⁻¹ |
| `seasonal_climatology.csv`, `fig_seasonal_cycle.png` | mean, minimum and maximum of each calendar month over the complete years | t month⁻¹ |
| `mean_intensity_<group>.tif`, `fig_maps_mean_intensity.png` | mean annual emission per area of each cell | kg ha⁻¹ yr⁻¹ |
| `drivers_annual.csv` | annual regional means of 2 m temperature, PAR, soil water and LAI as the model read them | °C, W m⁻², m³ m⁻³, m² m⁻² |
| `correlations_emissions_drivers.csv` | Pearson r of each group with each driver, raw and detrended | - |
| `trends.csv` | Theil-Sen slope, its 95 % interval, Mann-Kendall and Hamed-Rao p | t yr⁻² and % of the mean per year |
| `zones_annual_totals.csv` | annual totals of each zone of the configuration | t yr⁻¹ |

Conventions that hold for every file:

* **Mass of the compounds.** The tonnes are the mass of the compounds of each class (chapter 1, section 1.1), not the
  mass of carbon. Another inventory given in carbon, or in moles of a mechanism species, must be converted before a
  comparison.
* **Isoprene includes ρ = 0.96**, the canopy loss and production factor of Guenther et al. (2006), applied when the run
  is read; the files of notebook 10 hold isoprene without it (section 1.1).
* **Groups.** BVOC = the 17 classes other than soil NO and CO; MT = the five monoterpene classes; SQT = SQT_HR +
  SQT_LR; OVOC = MEOH + ACTO + ETOH + ACID + OXPROD; STRESS_OTHER = STRESS + OTHER; SOA = MT + SQT. "SOA" names the
  precursors of secondary organic aerosol, not aerosol formed: how much aerosol forms depends on the chemistry of the
  atmosphere, which an emission model does not compute.
* **Region.** A cell belongs to the region when its centre lies inside the boundary. Zone totals add the cells whose
  centre lies in each zone (notebook 10); notebook 13 of the example splits each cell among the municipalities by its
  area inside each one, so its municipal totals add up to 99.83 % of the state total of 2025 (chapter 7, section 7.10).
* **Time.** Months and years are calendar months and years in UTC, the time of ERA5-Land (chapter 8, section 8.4).

## 9.2 Annual totals

The level of an annual total is set mainly by the emission factors and by the area of vegetation, and its variation
from year to year by the activity factors. In MEGANv3.21 the emission activity of every class is the product of the
leaf area and the activity factors (chapter 1, section 1.2). The light-independent part of each class grows
exponentially with leaf temperature, with the coefficient β of the class, and the light-dependent part follows a
temperature response with an optimum and a light response (section 1.3). A warmer, sunnier year with denser foliage
therefore gives larger totals of the BVOC classes; soil NO follows its own algorithm (section 1.8).

In the example, BVOC averaged 657 kt yr⁻¹ over 2003-2025, with a coefficient of variation of 7.9 %; isoprene made
33.2 % of the mean BVOC, monoterpenes 25.8 % and oxygenated VOC 19.3 % (`notebook13_meso2024/key_numbers.csv`, keys
`bvoc_mean_kt`, `bvoc_cv_pct`, `mean_*`).

**Comparing with another inventory.** Check, before comparing numbers: the species definition (which classes make
"BVOC"), the mass convention (compound or carbon), the factor ρ, the area (cell-centre mask or area weighting, land or
land and water), the year, the land-cover map and its resolution, the emission factors and the soil-NO scheme. In the
example, the 2019 totals are 4 to 19 times those of the previous state inventory (Hoinaski et al., 2022): isoprene
252.4 Gg against 33.3 Gg, CO 27.8 Gg against 7.50 Gg, soil NO 10.8 Gg against 0.571 Gg (keys `cmp2019_*`). Such a ratio
says that the two inventories differ, and the checklist above says where to look; it does not say which one is right.
Notebooks 18 and 23 of the example go further, comparing the 2019 run class by class and cell by cell with the BRAIN
inventory of Hoinaski et al. (2024), which used the same MEGANv3.21 code with other inputs (`results/v10/brain_comparison`,
`results/v10/brain_differences`).

## 9.3 The seasonal cycle

Monthly totals follow the seasonal cycles of temperature, light and leaf area. In the example, BVOC peaks in January
(100 kt month⁻¹ on average) and is lowest in June and July (about 20 kt month⁻¹)
(`notebook13_meso2024/seasonal_climatology.csv`). Soil NO follows the soil-temperature factor and its growing season
(section 1.8), so its cycle can differ in shape from that of the BVOC. Where temperature and day length vary little
through the year, as near the equator, the activity factors of temperature and light vary little too, and the seasonal
cycle of a run follows mainly the leaf area and the soil water.

The minimum-maximum band of `fig_seasonal_cycle.png` shows the spread between years of each calendar month; a month
whose band is wide is one where the interannual variation of the drivers matters most.

## 9.4 Maps of intensity

`mean_intensity_<group>.tif` divides the mean annual emission of each cell by its area, so cells of different size and
regions of different area can be compared. The pattern follows first the emission factor of each cell, which weights
the factors of the ecotype by the fractions of trees, shrubs, herbs and crops in the cell (chapter 1, section 1.7);
`results/v10/posthoc_review/ecotype_tree_ef.csv` lists the tree emission factors of the ecotypes of Santa Catarina.
Leaf area and temperature then modulate the pattern; in the model, the cooler air of high ground lowers the emission
of a forest that would emit more at sea level.

Totals by polygon (`zones_annual_totals.csv`) favour large units; divide by the area of each unit before ranking them.
In the example, the SOA-precursor intensity of the six mesoregions follows their native-forest share (MapBiomas 2024):
Pearson r = 0.971 with a 95 % interval of 0.750 to 0.997, Spearman ρ = 0.943 with exact p = 0.017 (chapter 7,
section 7.8). With six units, this is a description of the six mesoregions, not a law that would hold elsewhere.

## 9.5 Drivers and correlations

The drivers of `drivers_annual.csv` are the inputs as the model read them (monthly maps written by notebook 10), not
independent observations. A correlation between an emission series and a driver therefore shows which input moved the
emissions of the model from year to year; it is a diagnostic of the model, not independent evidence of a physical
mechanism. Two precautions apply:

* **Detrending.** When an emission and a driver both have a trend, the raw correlation mixes the trends with the
  year-to-year variation. `r_detrended` removes a linear trend from both series first. In the example, soil NO and
  temperature correlate at r = 0.30 on the raw series and 0.82 on the detrended series, because their trends have
  opposite signs (`fig_detrend.png` of notebook 13; chapter 7, section 7.5).
* **Collinearity.** Drivers covary. In the example, PAR and soil water correlate at r = −0.90, so their separate
  regression coefficients are unstable, and temperature and LAI covary as well (chapter 7, sections 7.5 and 7.6).
  A partial correlation, or a regression with variance inflation factors, tells whether a driver adds information
  beyond another.

## 9.6 Trends

`trends.csv` gives the Theil-Sen slope (Sen, 1968), the median of the slopes between all pairs of years, with its 95 %
interval; the Mann-Kendall test (Mann, 1945; Kendall, 1975); and the Hamed and Rao (1998) correction of that test for
autocorrelation. Slopes in % of the mean per year compare groups of different size. In the example, isoprene increased
by +0.88 % yr⁻¹ (95 % interval +0.33 to +1.50) over 2003-2025, with Mann-Kendall p = 0.004 (chapter 5, section 5.4).

How to read a trend of a model inventory:

* **It is the trend of the inputs seen through the model.** A positive temperature trend raises the BVOC through the
  activity factors; a trend in leaf area or in the land-cover map does the same through the leaf area and the emission
  factors. The regression of section 7.7 separates these parts in the example.
* **Changes in the input products can create trends.** A new satellite sensor, an orbit drift or a new land-cover
  collection change the inputs without a change of the vegetation. In the example, the statewide MODIS LAI of 2025 was
  11 % below its 2003-2023 mean, a decrease larger than in any earlier year in cells dominated by trees or herbaceous
  vegetation, while the orbits of Terra and Aqua drifted; the run uses a VIIRS-corrected LAI for 2024 and 2025
  (chapter 6, section 6.7), and the trends are also given ending in 2023 and in 2024 (chapter 7, section 7.4).
* **Several tests at once.** Testing nine groups multiplies the chance of a false positive; the Benjamini-Hochberg
  control of the false discovery rate (chapter 7, section 7.3) keeps it at the chosen level.
* **Length.** With fewer than about 10 years a monotonic trend test has little power; notebook 25 does not compute
  trends below `MIN_YEARS_STATS = 10` complete years.

## 9.7 What the results can and cannot support

The results of a run **can** support: the spatial distribution of the emission potential given the land cover and
the MEGAN emission factors; the seasonal cycle and the year-to-year variation as driven by the reanalysis and the
satellite LAI; the effect of a change of input or setting, measured by running the model twice (as tests S1 and S3 and
the sensitivity runs of notebook 20 do for the example); and a comparison of regions computed in the same way.

They **cannot**, by themselves, establish the absolute magnitude of the emissions. Every class scales with its
emission factor, and the emission factors of the MEGAN tables come from measurements made elsewhere; the repository
checks that the code computes what MEGANv3.21 computes (notebook 12) and that the run is internally consistent
(notebook 11), not that the emission factors suit the vegetation of the region. Flux or concentration measurements in
the region are what would test the magnitude; the satellite test of chapter 10 tests the variation, and for the
example it supports the month-to-month anomalies of isoprene but not the spatial pattern, the seasonal cycle or the
magnitude (chapter 10, section 10.6).

Settings that bound what a run represents, all written in the configuration and in `provenance_v10.json`:

* constant wilting point (0.12 m³ m⁻³ by default) for the soil-moisture factor;
* light-dependent fractions of the monoterpene classes taken from the class constants, where the Fortran can read a
  gridded map;
* soil NO from the scheme of Yienger and Levy (1995) as coded in MEGANv3.21, on the path that uses air temperature
  only (no soil temperature, soil moisture or rain pulses), with the growing season of section 1.8;
* CO₂ inhibition of isoprene and the stress factors off;
* one canopy type for all trees of a run;
* coastal and island cells with the meteorology of the nearest land cell (FIX-1), listed in the fix logs.
