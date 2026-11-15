# 8. Number map: where each number comes from

This chapter maps the numbers of the article (English and Portuguese), of its Online Resource 1 and of the state
inventory report to the files of `results/v10/` that hold them. All values use the
VIIRS-corrected LAI of 2024 and 2025 (article Section 2.3.3); the values with the uncorrected MODIS LAI are in
`results/v10_modis_lai_2024_2025/`. Three files carry most of the numbers as key-value tables:

| File | Keys | Written by (`reproduce.py` step) |
|---|---|---|
| `notebook13_meso2024/manuscript_numbers.csv` | columns `key`, `description`, `value` (full precision), `en` and `pt` (as printed) | notebook 13 (`nb13`) |
| `report/report_numbers_v10.csv` | columns `key`, `value` (as printed, Portuguese), `source` | `analysis/v10/report/10_report_numbers_v10.py` (`report`) |
| `final_review/final_review_numbers.csv` | columns `key`, `value`, `note` | `analysis/v10/review/final_review_diagnostics.py` (`final`) |

"MN", "RN" and "FR" stand for the three files above. Section, figure and table numbers are those of the article
(Figures 1-11 and S4-S5; the figure files keep the names given in chapter 9). The report numbers its sections
differently; its section is given where it holds a number that the article does not.

## Abstract and conclusion

| Number | File and key |
|---|---|
| mean BVOC 657 Gg yr⁻¹; isoprene 33 %, monoterpenes 26 %, oxygenated VOC 19 % | MN `bvoc_mean_kt`, `mean_ISOP`, `mean_MT`, `mean_OVOC` |
| 72 % emitted from October to March | MN `season_BVOC` |
| isoprene +0.88 % yr⁻¹ (p = 0.004); soil NO −0.37 % yr⁻¹ (p = 0.002); significant after FDR control, also over 2003-2023 | MN `trend_ISOP`, `trend_NO`; FR `trend_*`; `final_review/trend_robustness.csv` |
| isoprene trend 0.65 % yr⁻¹ with CO₂ inhibition | FR `isop_trend_with_co2_inhibition` |
| 2019: isoprene, CO and soil NO 7.6, 3.7 and 18.9 times the previous inventory | MN `cmp2019_ISOP`, `cmp2019_CO`, `cmp2019_NO` |
| 2019 isoprene 21 % below BRAIN (ratio 0.79 on true cell areas) | `brain_comparison/comparison_2019_classes_t.csv` (`ratio_vs_true`) |
| settings tested one at a time: 2019 BVOC −2.3 % to +2.6 %; planted-forest emission factors: isoprene −19 % to +14 % | `sensitivity_2019/full/variants_summary.csv`; `notebook15/S3_fixlog_annual.csv`; `notebook15/s1/s1_groups.csv`; `ijest_review/plantation_ef_sensitivity_2019.csv` |
| TROPOMI anomalies r = 0.45 | `notebook14/tier3_anomalies_cf0.3_sza70_fire120.csv` |
| PAR fraction explains about one sixth of the isoprene gap to BRAIN (ratio 0.79 to 0.82) | `brain_differences/F_sensitivity_vs_brain_full.csv`; `brain_differences/summary.csv` (block F) |

## 2 Materials and methods

| Section | Number | File and key |
|---|---|---|
| 2.1, Figure 1 | MapBiomas 2019 shares of the classified area and growth forms on the grid | `figures/study_area_shares.csv` (step `extra`) |
| 2.2 | ρ = 0.96 for isoprene (Guenther et al., 2006), applied when the run is read | `RHO_ISOP` of notebooks 11, 13, 14 and 15; `RHO` of `analysis/v10/common.py` |
| 2.2 | Python-Fortran equivalence (99th percentile below 3 × 10⁻⁵; totals within 2 × 10⁻⁶; Online Resource 1, Section S4) | `benchmark/results*/`; `benchmark/colab_run/12_fortran_benchmark_ran.ipynb` |
| 2.3.1 | five ecotypes; three cover 98 % of the state; tree entries and factors (Online Resource 1, Table S8) | `ijest_review/ecotype_tree_ef.csv` (step `ijest`) |
| 2.3.2 | largest monthly cell LAI 6.67; 12 cell-months with LAI gaps | FR `lai_max_cell_all_months`, `lai_gap_cell_months` |
| 2.3.3 | MODIS LAI −4.8 % (2024) and −11.1 % (2025); hindcast RMSE 0.362 and 0.386 m² m⁻²; corrected LAI +5.0 % and +11.6 % | `lai_viirs_check/` (step `laicorr`); `ijest_review/lai_annual_by_growth_form.csv` |
| 2.3.3 | three VIIRS satellites, February-December change from 2023 and corrected LAI from each record | `lai_noaa20_check/three_satellites_deviation_feb_dec_vs_2023.csv`, `three_satellites_corrected_lai_feb_dec.csv`, `corrected_lai_by_viirs_record.csv` (notebook 22) |
| 2.4 | 190 to 191 cells with filled meteorology (1,301 km²) | MN `met_filled`; FR `filled_cells_area` |
| 2.4 | grid of 144 × 224 cells; 13,970 cells and 95,969 km² inside the state | MN `grid_facts` |
| 2.4 | months counted from 03 UTC instead of 00 UTC: monthly monoterpenes and soil NO at most 0.3 %, isoprene unchanged, annual totals at most 0.01 % | `final_review/month_boundary_effect.csv`, `.log` (step `final`) |
| 2.6 | global CO₂ 375 ppm (2003) to 426 ppm (2025); γCO₂ 1.026 to 0.979 | FR `isop_trend_with_co2_inhibition` (note); `data/external/co2_annmean_gl.txt` |
| 2.6 | cell areas of the authalic sphere 0.12 % to 0.19 % above the ellipsoidal ones | FR `cell_area_sphere_minus_ellipsoid_pct` |
| 2.6, Figure S4 | example of the detrended correlation: soil NO and temperature, r = 0.30 (raw) and 0.82 (detrended) | MN `r_NO_T`; file `fig_detrend.png` (notebook 13) |
| 2.6 | the five settings of the 2019 reruns and the station-fitted coefficient (10.78 K km⁻¹) | `sensitivity_2019/full/` (notebook 20); `era5land_vs_inmet/lapse_rate.json` (notebook 21) |
| 2.7 | 16:30 UTC = 12:55 to 13:17 local mean solar time | FR `overpass_local_solar_time` |
| 2.7, Figure 2 | the five verification steps and their results | drawn by `analysis/v10/schematics/verification_flowchart.py` (step `extra`) from the files of this chapter |
| Table 1 | β, LDF and the other class coefficients | `src/megan_sc/engine/species.py` (from MEGVEA.EXT; chapter 1) |

## 3.1 Totals, composition, seasonality and trends

| Number | File and key |
|---|---|
| Online Resource 1, Table S1 (annual totals by group) | `notebook13_meso2024/table2_annual_totals_SC_groups.csv` |
| 2025: BVOC 672.9 Gg, isoprene 228.6 Gg, monoterpenes 171.7 Gg | MN `y2025_text`, `y2025_*_t` |
| maximum 739 Gg in 2019 and 2024, minimum 552 Gg in 2008; coefficient of variation 7.9 % | MN `bvoc_max`, `bvoc_min`, `bvoc_cv_pct`; Table S1 |
| Theil-Sen slopes with 95 % intervals and Mann-Kendall p of every group | `notebook13_meso2024/trends.csv`; MN `trend_*` |
| Hamed-Rao threshold 0.41 and largest autocorrelation 0.31; soil NO p = 0.013 with the lag-1 term forced | FR `hr_threshold`, `hr_max_abs_acf`, `NO_p_with_lag1` |
| Benjamini-Hochberg q; p from 2005; trends ending in 2024 | FR `trend_*`; `final_review/trends_fdr_hr.csv` |
| 2003-2023 trends (isoprene +0.82 %, soil NO −0.43 % yr⁻¹; q ≤ 0.048); TFPW p ≤ 0.024; Benjamini-Yekutieli q (Online Resource 1, Table S9) | `final_review/trend_robustness.csv`, `trend_robustness.log` (step `final`) |
| MODIS LAI of 2025 11 % below the 2003-2023 mean (5.5 standard deviations) without the correction; 0.8 % below with it; 12 and 3 of 14 months below one standard deviation | `results/v10_modis_lai_2024_2025/final_review/final_review_numbers.csv` (`lai_2025`); FR `lai_2025`, `lai_months_below_1sd_nov2024_dec2025` |
| mean of each group and its share of BVOC (Figure 4a); isoprene:monoterpene ratio 1.29:1 | MN `mean_*`, `isop_mt_ratio` |
| seasonal cycle: January 100.4 Gg and July 20.1 Gg; isoprene 36.2 and 5.0 Gg; soil NO 674 to 1,060 t; isoprene share 25 % to 36 % (Figure 4b, 4c) | MN `season_*`, `season_isop_share`; `notebook13_meso2024/seasonal_cycle_monthly.csv`, `seasonal_climatology.csv` |

## 3.2 Interannual variability and drivers

| Number | File and key |
|---|---|
| detrended r of the BVOC classes with temperature 0.87 to 0.92, LAI 0.62 to 0.76, PAR 0.08 to 0.29, soil water −0.33 to −0.15 (Figure 5) | MN `rdet_bvoc_classes_*`; `notebook13_meso2024/correlations_emissions_drivers.csv` |
| 2008 (17.0 °C, 552 Gg); 2019 and 2024 (18.6 and 18.7 °C, 739 Gg each) | `notebook13_meso2024/drivers_annual_SC.csv`; Table S1 |
| driver series | `notebook13_meso2024/drivers_annual_SC.csv` |

## 3.3 SOA precursors and spatial distribution

| Number | File and key |
|---|---|
| Table 2: total, share, intensity, sesquiterpene share and native forest of each mesoregion (Figures 6 and 7) | `notebook13_meso2024/table3_mesoregions_soa.csv`; MN `meso_<name>` |
| intensity against native forest: r = 0.97 (95 % CI 0.73 to 0.997), ρ = 0.94, exact p = 0.017 (Figure S5) | MN `meso_forest_r`; FR `meso_r_fisher`, `meso_spearman_exact` |
| sesquiterpenes 7.2 % to 8.1 % of the precursors in every mesoregion | MN `meso_sqt_share_range` |
| sum of the six mesoregions 0.18 % below the state total; attributed area 1.0014 to 1.0019 times the IBGE area | MN `meso_sum`; FR `meso_area_ratio` |
| Table 3 (the twelve largest totals in 2025) and Figure 8 | `notebook13_meso2024/table4_municipalities_2025.csv`; MN `mun_top12_bvoc` |
| intensity 14 to 136 kg ha⁻¹ yr⁻¹ (median 70; maximum Itapoá); top Jenks class from 102.9 (37 municipalities) | MN `mun_intensity`, `mun_top_jenks`; `figures/maps_class_limits.csv` |
| attributed area 0.999 to 1.006 times the IBGE area; sum of the 295 municipalities 0.17 % below the state BVOC of 2025 | MN `mun_area_check`, `mun_cells`; FR `mun_area_ratio` |

## 3.4 Comparison with the previous state inventory and with BRAIN (2019)

| Number | File and key |
|---|---|
| Table 5: 2019 totals of this inventory and of Hoinaski et al. (2022), ratios 7.6, 3.7, 18.9 | MN `cmp2019_*`; the three totals of the previous inventory are typed in notebook 13 from its Table 12 |
| isoprene: mean 218 Gg, 6.5 times the previous value, 6.9 in 2025 | MN `isop_mean_ratio` |
| Figure 9 (municipal maps of 2019 and 2025) | `notebook13_meso2024/municipalities_2019.csv`, `table4_municipalities_2025.csv` |
| BRAIN rates as delivered 1.26-1.27 times the true-area totals; ratios 0.79 (isoprene), 0.85 (monoterpenes), 0.87 (summed VOC), 1.13 (soil NO); monthly r = 0.96; cells r = 0.63 and 83 % within a factor of 2 (Online Resource 1, Table S12) | `brain_comparison/comparison_2019_classes_t.csv`, `comparison_2019_monthly_t.csv`, `comparison_2019_cells.csv` (notebook 18) |
| BRAIN isoprene 9.96 times, CO 3.81 times and soil NO 16.8 times the previous inventory | `brain_comparison/comparison_2019_classes_t.csv` (`brain_true_over_previous`) |
| PAR 0.5: ratio 0.79 → 0.82; ethanol 0.75 → 0.87 without GAMBD; soil-moisture factor < 0.1 % | `brain_differences/A_switches.csv`, `summary.csv` (notebook 23) |
| cell-type and mesoregion medians of the isoprene ratio; 64 % fixed spatial pattern | `brain_comparison/comparison_2019_cells.csv`, `comparison_2019_mesoregions_t.csv`; `brain_differences/D_seasonal_decomposition.csv` |
| terrain slope 1.39 per km (1.20 to 1.60), predicted 1.46 (0.98 to 1.84), −0.03 after the terrain-adjusted run | `brain_differences/E_terrain_with_interval.csv`, `F_lapse_vs_terrain_full.csv` |
| TROPOMI blocks 2019: r = 0.63 (BRAIN) and 0.28 (here); pressure 0.82, 0.52, 0.03; partial 0.40 and 0.46 | `brain_differences/G_tropomi_2019.csv` |
| soil NO in western SC: 53.6 % (BRAIN) and 34.0 % (here) | `brain_differences/H_soil_no_mesoregions.csv` |

## 3.5 Run verification and TROPOMI evaluation

| Number | File and key |
|---|---|
| checks C0 to C5: 137 PASS and 1 FAIL in 138 | `notebook11/verify_v10_report.csv` |
| C3 in 2015: 4 of 71 cells; +0.003 % on the state methanol; LAI-adjusted test 0 of 69; pooled 42 of 4,785 (24 expected) | `diagnostics/c3/`; FR `c3_pooled_*`, `c3_distinct_cells` |
| eight steps per day in every month; annual = sum of monthly files | `diagnostics/check_run_integrity.log` |
| spatial r = 0.12 and 0.18 (95 % CI −0.07 to 0.29 and −0.03 to 0.37); kernels 0.16, 0.09, −0.10; climatology r = −0.26 (p = 0.41); anomalies r = 0.45 (N_eff = 72); model anomalies with temperature r = 0.79; annual r = −0.30; slope test r = 0.18 | `notebook14/tier1_*`, `tier2_*`, `tier3_*`, `tier4_*`; `diagnostics/nb14/nb14_diagnostics.log` (D7) |
| split-half r 0.76 and 0.86; ΔHCHO with surface pressure 0.82 and 0.88; isoprene with pressure ≤ 0.10; partial r 0.15 to 0.21 | `diagnostics/nb14/nb14_diagnostics.csv`, `nb14_D5_terrain.csv` |
| background seasonal cycle r = 0.67; raw column r = 0.37; screening r = 0.80; without August-October r = 0.41 (N_eff = 58, p = 0.001); Nov-Apr climatology −0.10 | `diagnostics/nb14/nb14_diagnostics.csv` (D3, D4) |
| background anomalies r = 0.43; partial r 0.65 (background held fixed) and 0.48 (temperature held fixed); seasons 0.54, 0.45, −0.21, 0.67 | FR `anom_*` |
| summer ΔHCHO 2.2 × 10¹⁵ molecules cm⁻²; isoprene 1.8 × 10¹² atoms C cm⁻² s⁻¹; lifetime 0.8 to 1.2 h (1.2 to 1.8 h with the bias correction); ratio 45 % to 80 % (71 % to 125 %) | FR `implied_tau_h_*`, `ratio_obs_expected_palmer_tau`; `final_review/tropomi_magnitude_grid.csv` |
| S3: γSM state mean 0.992 to 1.000, BVOC −0.05 % | `notebook15/S3_fixlog_monthly.csv`, `S3_fixlog_annual.csv` |
| S1: −0.06 % (BVOC) and −0.18 % (isoprene) over 28 days | `notebook15/s1/s1_groups.csv`, `s1_class_totals.csv`, `s1_cell_ratios.csv` |

## 3.6 Interpretation and 3.7 limitations

| Section | Number | File and key |
|---|---|---|
| 3.6 | trend of ln(emission) split into temperature, LAI and time term (isoprene +0.92 = +0.35 + 0.16 + 0.40 % yr⁻¹; soil NO temperature +0.15, time term −0.45 % yr⁻¹); temperature trend 0.027 K yr⁻¹ (OLS; MK p 0.09) | FR `ols_trend_decomposition_*` (the OLS trend of temperature, 0.0268 K/yr, is in the note of these keys), `residual_year_term_*`, `T_trend` (Theil-Sen 0.0328 K/yr and the Mann-Kendall p 0.091) |
| 3.6 | detrended r(T, LAI) = 0.51; log-LAI coefficient of isoprene 0.61 (MODIS LAI) and 1.47 (corrected LAI) | FR `r_detrended_T_LAI`, `ols_trend_decomposition_ISOP` (note); `results/v10_modis_lai_2024_2025/final_review/final_review_numbers.csv` |
| 3.6 | growth-form areas 2003-2024 (tree +9.5 %, herbaceous −18.8 %, crop −0.2 %) and emission potential (isoprene +0.38 % yr⁻¹, soil NO −0.45 % yr⁻¹) | FR `growthform_area_*`, `ef_potential_*`; `final_review/landcover_ef_potential_SC.csv` |
| 3.6 | CO₂ inhibition: isoprene +0.65 % yr⁻¹ (p = 0.018; q = 0.03), period mean +0.3 %; with the MODIS LAI, +0.82 % yr⁻¹ and two significant trends | FR `isop_trend_with_co2_inhibition`, `fdr_with_co2_inhibition`, `isop_mean_with_co2_inhibition`; `results/v10_modis_lai_2024_2025/notebook13_meso2024/trends.csv` |
| 3.6 | PAR coefficient of variation 3.2 %; temperature standard deviation 0.48 °C; r(PAR, soil water) = −0.89 (detrended) | FR `driver_interannual` |
| 3.7 | plantation share 9,960 km² (20.6 % of the tree area); Pinus −18.7 % and Eucalyptus +13.7 % isoprene | `ijest_review/ijest_review_numbers.csv`, `plantation_ef_sensitivity_2019.csv` |
| 3.7 | PAR 0.5 +4.2 %, tropical canopy −2.1 %, clay wilting point −2.2 % (isoprene); terrain +0.2 % (+2.4 % in Southern SC) | `sensitivity_2019/full/variants_summary.csv` (Online Resource 1, Table S10) |
| 3.7 | ERA5-Land minus INMET: daytime +0.17 K (−0.09 to +0.42); radiation ratio 1.028 (54 stations); 10.8 K km⁻¹ (8.3 without A845) | `era5land_vs_inmet/bias_summary_2019.csv`, `radiation_2019_daily_totals.csv`, `lapse_rate_fits.csv`, `lapse_rate_influence.csv` (Online Resource 1, Table S11) |
| 3.7 | NOAA-20 correction of 2025 2.8 % above Suomi NPP (2018-2023 reference) | `lai_noaa20_check/corrected_lai_by_viirs_record.csv` |
| 3.7 | isoprene areal flux 2.0 g C m⁻² yr⁻¹ (cell median 2.0, maximum 6.1) | MN `isop_flux` |
| 3.7 | LAI gaps filled in 0 to 7 cell-months per year | MN `lai_gaps` |

## Online Resource 1

| Item | File |
|---|---|
| Table S6 (2024-2025 with both LAI inputs) | `notebook13_meso2024/table2_annual_totals_SC_groups.csv` and the same file in `results/v10_modis_lai_2024_2025/` |
| Section S2.5.6 and Table S7 | `lai_noaa20_check/` (notebook 22); r = 0.999 between the Suomi NPP record read from the HDF5 files and the Earth Engine record: `lai_noaa20_check/notebook22_printed_checks.txt` (printed by the notebook, not written to a table) |
| Table S8 | `ijest_review/ecotype_tree_ef.csv` |
| Section S4 | `benchmark/colab_run/12_fortran_benchmark_ran.ipynb` |
| Section S5 | `diagnostics/c3/`; FR `c3_pooled_*` |
| Table S9 | `final_review/trend_robustness.csv` |
| Tables S10 and S11 | `sensitivity_2019/full/variants_summary.csv`; `era5land_vs_inmet/bias_summary_2019.csv`, `radiation_2019_daily_totals.csv` |
| Section S8, terrain run: Δz mean −0.8 m, SD 85.6 m; ΔT −2.69 to +2.65 K (1st and 99th percentiles) | `sensitivity_2019/full/notebook20_printed_terrain_stats.txt` (printed by notebook 20, not written to a table) |
| Section S9: about 6 K across classes; Section S2.5.3: statewide corrected LAI within 0.1 percentage points; NOAA-21 record from 10 February 2023 | recomputed from `brain_differences/B_fingerprint_properties.csv` (slope of the log ratio on the temperature sensitivity, −6.03); from `lai_viirs_check/lai_correction_2024_2025_monthly.csv`; first VJ215A2H granule listed by notebook 22 |
| Section S9 and Table S12 | `brain_comparison/` (notebook 18); `brain_differences/` (notebook 23) |

## Numbers of the report only

| Report section | Number | RN key |
|---|---|---|
| 3.2 | means with CO and NO, 2025 totals by group | `mean_*_kt`, `y2025_*_kt`, `cv_BVOC_pct`, `peak`, `trough` |
| 3.3 | monthly climatology of isoprene (January about 7 times June) | `isop_clim_2003_2025`, `isop_clim_2015_2025` |
| 3.4, 4.2 | raw and detrended correlations, driver intercorrelations, VIF | `rraw_*`, `rdet_*`, `driver_corr`, `VIF`, `r_T_SOILM` |
| 3.5, 4.1 | largest masses and intensities by mesoregion; 2025 intensities (BVOC 91, 88, 82 ... 55 kg ha⁻¹ yr⁻¹) | `meso_top_*`, `soa_*`, `meso2025_*` (`report/meso_2025_v10.csv`) |
| 3.6 | municipal intensities by intermediate region, bottom and top quintiles | `mun_*` (`report/muni_2025_ranked_v10.csv`) |
| 4.2, Figure 12, Table 7 | standardised regression: β of temperature 0.65 to 0.84, R² 0.88 to 0.94; NO R² 0.28 | `beta_*`, `R2_BVOC`, `SQT`, `MT`, `NO_R2`, `CO` (`report/std_regression_v10.csv`) |
| 4.2 | isoprene +20.6 % per K (17.5 % detrended) in the simple regression; 0.1355 per K in the terrain test | `isop_pct_per_K` (`report/temp_sensitivity_v10.csv`); `sensitivity_2019/full/lapse_temperature_sensitivity.json` |
| Figure 13 | annual temperature against soil water, r = −0.22 (p = 0.30) | `r_T_SOILM` |
| 3.7, 4.2 | soil water below the γSM threshold in 0.09 % of the cell-months; monthly γSM below 1 in 0.63 % | `soilm_below_016`, `gamsm_lt1` |

When a number is not in this map, search the three key-value files for it; every value the documents print comes from
a file of `results/v10/`.
