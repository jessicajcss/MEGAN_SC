# 7. The statistics

Every statistical result of the Santa Catarina example rests on one of the estimators below, and notebook 25 applies
the same ones to any run. For each estimator, this chapter gives the definition, the code that computes it and the
values of the example. All series are **regional totals**: the sum over the cells whose centre lies inside the region
(13,970 cells for Santa Catarina), by year (n = 23, 2003-2025) or by month.

Notebook 13 computes the descriptive statistics, the trends and the correlations with the drivers. The estimators added
as additional diagnostics, marked "post hoc" in the logs because they were defined after the first results were seen, are computed by
`analysis/v10/review/final_review_diagnostics.py` and `analysis/v10/nb14_diagnostics.py`.

## 7.1 Level and variability

* **Mean** of the annual totals, 2003-2025: BVOC 656.5 kt yr⁻¹ (`notebook13_meso2024/key_numbers.csv`, `bvoc_mean_kt`).
* **Coefficient of variation** with the sample standard deviation (ddof = 1): 7.9 % for BVOC (`bvoc_cv_pct`).
* **Shares** of each group in the mean BVOC: isoprene 33.2 %, monoterpenes 25.8 %, oxygenated VOC 19.3 %
  (`mean_*`). Isoprene includes the factor ρ = 0.96 (chapter 1, section 1.1).

## 7.2 Trends: Theil-Sen slope and Mann-Kendall test

**Slope.** The Theil-Sen estimator (Sen, 1968) is the median of the slopes between all pairs of years; its 95 %
interval comes from the distribution of those pairwise slopes (`scipy.stats.theilslopes`). The documents give it in
per cent of the 2003-2025 mean per year.

**Test.** The Mann-Kendall statistic (Mann, 1945; Kendall, 1975) counts the pairs of years in which the later value is
larger, minus those in which it is smaller; its p-value tests the absence of a monotonic trend
(`pymannkendall.original_test`).

**Autocorrelation.** Serial correlation inflates the significance of the test. The Hamed-Rao adjustment of the variance
(Hamed; Rao, 1998) rescales the variance of the statistic by the autocorrelation of the ranks of the detrended series,
using only the lags whose autocorrelation exceeds 1.96/√n (0.409 for n = 23). Here the largest absolute rank
autocorrelation of any group at any lag is 0.308, so the adjustment keeps every p-value unchanged
(`final_review_numbers.csv`: `hr_threshold`, `hr_max_abs_acf`, `hr_changes_p`). Forcing the lag-1 term into the
adjustment for soil NO, whose lag-1 rank autocorrelation is the largest (0.308), raises its p-value from 0.002 to 0.013
(`NO_p_with_lag1`, post hoc).

| Group | Theil-Sen (% yr⁻¹, 95 % CI) | Mann-Kendall p | Benjamini-Hochberg q |
|---|---|---|---|
| BVOC | +0.61 (0.17 to 1.09) | 0.008 | 0.019 |
| Isoprene | +0.88 (0.33 to 1.50) | 0.004 | 0.017 |
| Monoterpenes | +0.56 (0.13 to 1.02) | 0.011 | 0.020 |
| Sesquiterpenes | +0.76 (0.23 to 1.50) | 0.008 | 0.019 |
| Oxygenated VOC | +0.47 (0.04 to 0.89) | 0.015 | 0.023 |
| Light VOC | +0.36 (−0.03 to 0.79) | 0.064 | 0.081 |
| MBO | +0.52 (−0.09 to 1.04) | 0.081 | 0.081 |
| Soil NO | −0.37 (−0.58 to −0.12) | 0.002 | 0.017 |
| CO | +0.45 (−0.07 to 0.88) | 0.073 | 0.081 |

Sources: `notebook13_meso2024/trends.csv` (slopes, p) and `final_review/trends_fdr_hr.csv` (q).

## 7.3 Multiple testing: Benjamini-Hochberg

Nine trends are tested at once, so some small p-values are expected by chance. The procedure of Benjamini and Hochberg
(1995) sorts the nine p-values, p(1) ≤ ... ≤ p(9), and gives each the adjusted value
q(i) = min over j ≥ i of [9 · p(j) / j]; a trend is retained at a false discovery rate of 5 % when q < 0.05
(`statsmodels.stats.multitest.multipletests(..., method='fdr_bh')`). Six trends are retained (q = 0.017 to 0.023):
isoprene, soil NO, BVOC, monoterpenes, sesquiterpenes and oxygenated VOC; light VOC, MBO and CO are not (q = 0.081).
With the uncorrected MODIS LAI of 2024-2025, only isoprene and soil NO were retained
(`results/v10_modis_lai_2024_2025/`). The procedure holds for independent and positively dependent tests; section 7.4
adds the Benjamini-Yekutieli procedure (Benjamini; Yekutieli, 2001), which holds under any dependence.

## 7.4 Robustness of the trends to the period

Post hoc, in `final_review/trends_fdr_hr.csv`:

* starting in 2005 instead of 2003: p = 0.012 (isoprene), 0.020 (soil NO), 0.024 (BVOC), 0.027 (monoterpenes), 0.020
  (sesquiterpenes), 0.032 (oxygenated VOC);
* ending in 2024: isoprene +1.00 % yr⁻¹ (p = 0.0040), BVOC +0.67 % yr⁻¹ (p = 0.0080), monoterpenes +0.61 % yr⁻¹
  (p = 0.0080), sesquiterpenes +0.83 % yr⁻¹ (p = 0.0080), oxygenated VOC +0.51 % yr⁻¹ (p = 0.013), soil NO
  −0.36 % yr⁻¹ (p = 0.0057).

Post hoc, in `final_review/trend_robustness.csv` (`analysis/v10/review/trend_robustness.py`):

* 2003-2023, the years before the LAI correction applies: the same six trends, with p = 0.001 to 0.032 and
  Benjamini-Hochberg q ≤ 0.048 (isoprene +0.82 % yr⁻¹, soil NO −0.43 % yr⁻¹);
* trend-free pre-whitening (Yue et al., 2002; `pymannkendall.trend_free_pre_whitening_modification_test`, applied to
  every series without a prior test of the lag-1 autocorrelation): p ≤ 0.024 for the six groups over 2003-2025;
* Benjamini-Yekutieli over the nine groups: isoprene and soil NO q = 0.047, the other four 0.053 to 0.064;
* Benjamini-Hochberg over the eight groups without total BVOC (a sum of the others): q ≤ 0.024 (2003-2025) and
  q ≤ 0.051 (2003-2023).

The six trends hold in every variant of the period; under the Benjamini-Yekutieli procedure, only isoprene and soil NO
keep q < 0.05. A moving-block bootstrap of the Theil-Sen slope was tried and dropped, because with 23 years and
negative lag-1 autocorrelation its intervals were narrower than the Theil-Sen intervals.

## 7.5 Correlation with the drivers

The drivers are the state means of the annual means of the monthly maps of the run: 2 m temperature (°C), PAR,
soil water, and LAI averaged over cells with LAI > 0 (`notebook13_meso2024/drivers_annual_SC.csv`). The correlation of
each class with each driver is Pearson's r on the raw annual series and on the series after removing a linear trend
from both, so that a common trend does not create a correlation (`correlations_emissions_drivers.csv`). The figure
`fig_detrend.png` of notebook 13 shows the procedure with soil NO and temperature, whose linear trends have opposite
signs: r = 0.30 on the raw series and 0.82 on the detrended ones.

| Driver | r, detrended, over the 17 BVOC classes |
|---|---|
| temperature | 0.87 to 0.92 |
| LAI | 0.62 to 0.76 |
| PAR | 0.08 to 0.29 |
| soil water | −0.33 to −0.15 |

Soil NO correlates with temperature at r = 0.30 on the raw series and 0.82 on the detrended series: the downward trend
of soil NO runs against the upward tendency of temperature (Theil-Sen slope +0.033 K yr⁻¹, Mann-Kendall p = 0.09; key
`T_trend` of `final_review/final_review_numbers.csv`),
which lowers the raw correlation.

**Partial correlations** (post hoc, `final_review/partial_correlations_detrended.csv`): with temperature held fixed, the
correlation of each BVOC class with LAI ranges from 0.42 (MBO) to 0.75 (MT_CAMP); the detrended correlation between
temperature and LAI is 0.51. LAI therefore adds information beyond temperature, but the two drivers covary, and the
regression partition of the trends between them is unstable (section 7.7).

## 7.6 Standardised regression

`analysis/v10/report/10_report_numbers_v10.py` regresses each annual class series, in z-scores (population standard
deviation, ddof = 0), on the four drivers in z-scores, without intercept. The coefficients β are comparable across
drivers. Temperature has the largest β for every BVOC class (0.65 to 0.84), and R² is 0.88 to 0.94
(`report/std_regression_v10.csv`; keys `beta_*` and `R2_BVOC` of `report/report_numbers_v10.csv`).

Collinearity is measured by the variance inflation factor, VIF = 1/(1 − R²) of each driver regressed on the others:
temperature 1.61, PAR 5.38, LAI 1.68, soil water 5.38 (key `VIF`). PAR and soil water are strongly anticorrelated
(r = −0.90), so their separate β values are unstable, and a sign of one of them can flip without physical meaning.

## 7.7 Temperature sensitivity and the residual year term

**Sensitivity.** The slope b of ln(E) on the annual mean temperature gives 100 · (e^b − 1) % per kelvin: isoprene
+20.6 ± 2.2 % per K, or +17.5 % per K on detrended series (`report/temp_sensitivity_v10.csv`, key `isop_pct_per_K`).
A log-log "elasticity" on temperature in °C is not used, because it has no physical meaning.

**Residual year term** (post hoc). An ordinary least-squares fit of ln(annual emission) on temperature, ln(LAI) and the
year asks whether a trend remains once the two drivers are accounted for: isoprene +0.40 ± 0.11 % yr⁻¹ and soil NO
−0.45 ± 0.05 % yr⁻¹ (`residual_year_term_*`).

## 7.8 Mesoregions (n = 6)

The SOA-precursor intensity of the six mesoregions is compared with their native-forest share (MapBiomas 2024):

* Pearson r = 0.971, with the 95 % interval from Fisher's transformation, z = artanh(r) ± 1.96/√(n − 3), back-
  transformed: 0.750 to 0.997 (`meso_r_fisher`);
* Spearman ρ = 0.943, exact two-sided p = 0.017 by enumeration of the 720 permutations (the asymptotic p, 0.005, is
  not reliable at n = 6; `meso_spearman_exact`).

With six points the correlation is descriptive. The documents present it as such.

## 7.9 Monthly anomalies against TROPOMI

The monthly series of modelled isoprene and of ΔHCHO (2019-2025, n = 84) are deseasonalised by removing each calendar
month's mean, then correlated. Autocorrelation reduces the information in the series; the effective sample size of
Bretherton et al. (1999, their equation 31) is

```
N_eff = N · (1 − r1 · r2) / (1 + r1 · r2)
```

with r1 and r2 the lag-1 autocorrelations of the two series, and the p-value uses N_eff − 2 degrees of freedom. For
the anomaly correlation, r = 0.45, N_eff = 72, p < 0.001 (`notebook14/tier3_anomalies_*.csv`). Taking the 12 monthly
means removed by the deseasonalisation out of N as well gives N_eff = 61.5 and p = 2.4 × 10⁻⁴
(`anom_r_isop_dHCHO_neff_minus12`, post hoc). Partial correlations with the background anomaly or the temperature held
fixed (0.65 and 0.48) are discussed in chapter 10.

## 7.10 Municipalities

Municipal and mesoregion totals split the emission of each cell among the polygons it overlaps in proportion to the
area of the cell inside each one (area weighting, which assumes a uniform emission within the cell); intensities divide
by the official IBGE area. The fractions come from the intersection of each cell rectangle with the IBGE polygons in
longitude and latitude (`area_weights` in notebook 13 and in `analysis/v10/common.py`); within a 0.025° cell they differ
from equal-area fractions by less than 0.02 %. The area attributed in this way is 0.9988 to 1.0061 times the IBGE area
of each municipality (median 1.0017) and 1.0014 to 1.0019 times that of each mesoregion (`mun_area_ratio`,
`meso_area_ratio`); the cell areas are computed on a sphere, 0.12 % to 0.19 % above the ellipsoidal ones (chapter 1,
section 1.1). The IBGE areas of the municipalities equal 5.1 to 384 cells of 6.87 km² (median 34), and 17 to 485 cells
touch each one (median 54; `mun_cells`). State totals keep the cell-centre mask, so the sum of the municipalities is
99.83 % of the state BVOC of 2025 (`mun_area_check`) and that of the mesoregions 99.82 % of the state SOA precursors
(`meso_sum`): at the border and on the coast, a cell counted whole in the state total contributes to the municipalities
only the share of its area inside them. The maps use Jenks natural breaks of the municipal values (Fisher's exact
algorithm, five classes).

## 7.11 What is descriptive

* Correlations over the six mesoregions, and the annual TROPOMI correlation (n = 7, r = −0.30), are descriptive.
* Diagnostics defined after seeing the results are labelled post hoc in their logs (`final_review.log`, the
  `diagnostics/` logs and the folder `posthoc_review/`).
* Every p-value is two-sided.
