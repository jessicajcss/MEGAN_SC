# Results of the Santa Catarina example

The production run of the example (notebook 10 with `config/santa_catarina_2003_2025.yaml`, Google Colab, October 2026;
FIX-1 to FIX-8, `fixes.py` sha256[:16] a933b29602929a82; 3-hourly sampling; 2003-2025) is shipped in `data/run_v10/`.
Its years 2024 and 2025 were computed by notebook 17 with the VIIRS-corrected LAI (`src/megan_sc/laicorr.py`) and
installed by notebook 19; the results with the uncorrected MODIS LAI are in `../v10_modis_lai_2024_2025/`. The files of
this folder were computed from the run by notebooks 11, 13, 14, 15 and 18-23 and by the scripts of `analysis/v10/`.
They are the published results of the example and the reference against which `python reproduce.py compare` checks a
reproduction. Isoprene includes the canopy loss and production factor
ρ = 0.96 (Guenther et al., 2006), applied when the run is read (the files of `data/run_v10/` hold isoprene as the
MEGANv3.21 code computes it; the BRAIN comparison leaves it out on both sides), and municipal and mesoregion totals
split each cell among the polygons in proportion to its area inside each one.

| Folder | Content | Written by | `reproduce.py` step |
|---|---|---|---|
| `notebook11/` | `verify_v10_report.csv` (checks C0-C5 for each year: 137 PASS and 1 FAIL in 138), state totals by class (`annual_totals_SC_v10_*.csv`), `verify_v10_totals.png` | notebook 11 | `nb11` |
| `notebook13_meso2024/` | `key_numbers.csv` (76 keys, each with its value and its Portuguese and English wording), `table2_annual_totals_SC_groups.csv`, `table3_mesoregions_soa.csv`, `table4_municipalities_2025.csv`, `trends.csv`, `correlations_emissions_drivers.csv`, `drivers_annual_SC.csv`, `annual_totals_SC_by_class.csv`, `municipalities_2019.csv`, `seasonal_cycle_monthly.csv`, `seasonal_climatology.csv`, `figures/{pt,en}/fig1-8.png` and `fig_detrend.png`, `run_nb13.log` | notebook 13 (`13_santa_catarina_analysis.ipynb`) with `data/boundaries/sc_mesoregions_2024.gpkg` (mesoregions with the islands, MapBiomas 2024 forest share) | `nb13` |
| `notebook14/` | TROPOMI HCHO evaluation, screening `cf0.3_sza70_fire120`: the four pre-specified tests (`tier1` to `tier4`), the monthly series (`monthly_series_*.csv`) and the summary figure | notebook 14 | `nb14` |
| `notebook15/` | `S3_fixlog_{annual,monthly}.csv` (test S3); `s1/` (test S1: class totals, groups, cell ratios, diurnal check, figure, log) | notebook 15; `analysis/v10/s1_compare.py` | `nb15`, `s1` |
| `diagnostics/` | checks defined after inspection, labelled post hoc in their logs: `c3/` (the C3 failure of 2015 and its effect on the totals), `nb14/` (post hoc diagnostics of the TROPOMI comparison, D1-D7), `check_run_integrity.log` | `analysis/v10/` | `c3`, `nb14diag`, `integrity` |
| `report/` | statistics by class, mesoregion and municipality: `report_numbers_v10.csv` (104 rows) and its tables (standardised regression, temperature sensitivity, mesoregions in 2025, ranked municipalities, isoprene monthly climatology, driver correlations) | `analysis/v10/report/10_report_numbers_v10.py` | `report` |
| `figures/` | `pt/` and `en/` (Portuguese and English labels): `fig_study_area.png` (MapBiomas 2024 and the model growth forms; `analysis/v10/maps/study_area_map.py`, with `study_area_shares.csv`), `fig_verification.png` (the checks and the evaluation; `analysis/v10/schematics/verification_flowchart.py`), `fig_detrend.png` (detrended correlation, from notebook 13), `fig1-4.png` and `fig6.png` (copied from `notebook13_meso2024/figures/`), the maps `fig5.png`, `fig7.png` and `fig8.png` (`analysis/v10/maps/12_maps_wf_style.py`) and `fig9.png` (TROPOMI; `analysis/v10/fig9_tropomi.py`); `pt_report/` (Portuguese only): `fig9.png` (intensity by group of compounds), `fig12.png` (tree fraction and SOA-precursor intensity), `figA2-figA9.png` (monthly inputs of 2025 and annual intensity) and `fig10.png` (standardised regression); `maps_class_limits.csv` (class limits of every map) | as stated | `nb13`, `maps`, `report`, `extra` |
| `lai_viirs_check/` | post hoc check and correction of the 2024-2025 MODIS LAI: tables and figure of notebook 16 (`lai_viirs_modis_*`), hindcast of the candidate corrections 2013-2023 (`lai_correction_hindcast_*.csv`), corrected LAI of 2024-2025 (`lai_viirs_corrected_2024_2025.npz`, read by notebook 17), monthly and annual summaries, `lai_correction.png`, `lai_correction_map_2025.png`, settings | notebook 16; `analysis/v10/review/lai_viirs_correction.py` | `laicorr` |
| `final_review/` | additional diagnostics (Benjamini-Hochberg control, Hamed-Rao check, start-year and 2025 sensitivity of the trends, LAI of 2024-2025, partial correlations, residual year term, parts of the trends (temperature, LAI, land cover), growth-form areas and emission potential 2003-2025, effect of CO₂ inhibition, background control and magnitude of the TROPOMI comparison, municipal statistics, pooled check C3, filled coastal cells): `final_review_numbers.csv` (key, value, note), `final_review.log`, `trends_fdr_hr.csv`, `partial_correlations_detrended.csv`, `lai_monthly_SC.csv` and its z-scores, `tropomi_magnitude_grid.csv`, `filled_cells_area.csv`, `landcover_ef_potential_SC.csv`; `trend_robustness.csv` and `.log` (2003-2023, trend-free pre-whitening, Benjamini-Yekutieli); `month_boundary_effect.csv` and `.log` (monthly totals counted from 03 UTC, local midnight at UTC-3, instead of 00 UTC; `docs/TIME_CONVENTIONS.md`) | `analysis/v10/review/final_review_diagnostics.py`, `trend_robustness.py`, `month_boundary.py` | `final` |
| `posthoc_review/` | post hoc diagnostics: planted-forest emission factors (`plantation_ef_sensitivity_2019.csv`), LAI by growth form (`lai_annual_by_growth_form.csv`, `lai_monthly_state_z.csv`), `posthoc_review_numbers.csv` and `.log`, and the ecotypes of the state with their tree emission factors (`ecotype_tree_ef.csv`) | `analysis/v10/review/posthoc_review_diagnostics.py`, `ecotype_tree_ef.py` | `posthoc` |
| `notebook19/` | installation of the corrected 2024-2025: checks of the rerun, NOAA CO₂ file, comparison of both LAI inputs (`lai_sensitivity_summary.*`, `lai_sensitivity_table2_2024_2025.csv`), provenance | notebook 19 (Colab, 1 November 2026); file names in its tables are those of that date (`manuscript_numbers.csv` is today's `key_numbers.csv`, `ijest_review/` today's `posthoc_review/`) | none |
| `lai_noaa20_check/` | LAI of the VIIRS on NOAA-20 (2018-2025) and NOAA-21 (from February 2023) against Suomi NPP; corrected LAI from each record | notebook 22 (Colab, 8 November 2026) | none |
| `brain_comparison/` | 2019 against BRAIN (Hoinaski et al., 2024; data set of 2023): classes, CB6X species, months, mesoregions, cells, file checks; rates as delivered and on true cell areas | notebook 18 (Colab, 1 November 2026) | none |
| `brain_differences/` | tests of the explanations of the differences from BRAIN (blocks A-H: default switches, class fingerprint, cell regression, seasonal decomposition, terrain, sensitivity runs, TROPOMI 2019, soil NO) | notebook 23 (Colab, 3 November 2026); `E_terrain_with_interval.csv` and `F_lapse_vs_terrain_full.csv` from a local rerun of the updated notebook on the same inputs | none |
| `sensitivity_2019/` | one-at-a-time reruns of 2019 (`full`: 365 days; `design`: 28 days): PAR 0.5, terrain-adjusted temperature, tropical canopy, needleleaf plantations, clay wilting point | notebook 20 (Colab, 2 November 2026) | none |
| `era5land_vs_inmet/` | ERA5-Land against the INMET automatic stations, 2019: temperature bias, daily radiation, diurnal lag, station-fitted elevation coefficient and its influence diagnostics | notebook 21 (Colab, 2 November 2026); `lapse_rate_influence.csv` and `bias_summary_2019.csv` from a local rerun of its updated cell 8 on the same station table | none |

Two variables of the example have no figure of their own in the scripts: the annual temperature against soil moisture
(`notebook13_meso2024/drivers_annual_SC.csv`; key `r_T_SOILM` of `report/report_numbers_v10.csv`: r = −0.22, p = 0.30)
and the growth-form fractions of 2024 (`data/run_v10/2025/inputs/growthform_2024.tif`). Chapter 9 of the tutorial
(`docs/tutorial/09_reading_the_results.md`) explains how to read each kind of file of this folder.

## Summary of the checks

* **Notebook 11.** The C3 failure (2015, methanol, 4 of 71 comparable filled cells above the threshold where the
  binomial bound was 2) was accepted and reported: the deviations follow the LAI contrast that the 25 % LAI margin of
  the test still admits, none of the 69 LAI-adjusted tests fails, and replacing the four cells by their neighbours'
  median activity changes the 2015 state methanol by +0.003 % (`diagnostics/c3/`). Pooled over all years and classes,
  42 of 4,785 comparisons exceed the 99.5th percentile of the valid cells, where 24 are expected; methanol (23 against
  8) and MT_PINE (14 against 8) carry the excess, isoprene does not (5 against 8) (`final_review/`).
* **Test S1.** Over the 28 days, the 3-hourly step changes the state totals by −0.06 % (BVOC) and −0.18 % (isoprene);
  monthly isoprene differences range from −2.2 % (July) to +0.5 % (January), and cell ratios from 0.981 to 1.017 (1st
  and 99th percentiles). Emissions at the hours shared by both runs are identical. Test S2 (Earth Engine default LAI
  resampling) was not run.
