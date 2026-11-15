# Published results of the production run

The production run (notebook 10, Google Colab, October 2026; FIX-1 to FIX-8, `fixes.py` sha256[:16]
a933b29602929a82; 3-hourly sampling; 2003-2025) is shipped in `data/run_v10/`. Its years 2024 and 2025 were computed
by notebook 17 with the VIIRS-corrected LAI (`src/megan_sc/laicorr.py`; article Section 2.3.3) and installed by
notebook 19; the results with the uncorrected MODIS LAI are in `../v10_modis_lai_2024_2025/`. The files of this folder were computed from the run by notebooks 11, 13, 14,
15 and 18-23 and by the scripts of `analysis/v10/`. They are the source of every number, table and figure of the
article, of its Online Resource 1 and of the state inventory report, and the reference against
which `python reproduce.py compare` checks a reproduction. Isoprene includes the canopy loss and production factor
ρ = 0.96 (Guenther et al., 2006), applied when the run is read (the files of `data/run_v10/` hold isoprene as the
MEGANv3.21 code computes it; the BRAIN comparison leaves it out on both sides), and municipal and mesoregion totals
split each cell among the polygons in proportion to its area inside each one.

| Folder | Content | Written by | `reproduce.py` step |
|---|---|---|---|
| `notebook11/` | `verify_v10_report.csv` (checks C0-C5 for each year: 137 PASS and 1 FAIL in 138), state totals by class (`annual_totals_SC_v10_*.csv`), `verify_v10_totals.png` | notebook 11 | `nb11` |
| `notebook13_meso2024/` | **source of the documents**: `manuscript_numbers.csv` (76 keys), `table2_annual_totals_SC_groups.csv`, `table3_mesoregions_soa.csv`, `table4_municipalities_2025.csv`, `trends.csv`, `correlations_emissions_drivers.csv`, `drivers_annual_SC.csv`, `annual_totals_SC_by_class.csv`, `municipalities_2019.csv`, `seasonal_cycle_monthly.csv`, `seasonal_climatology.csv`, `figures/{pt,en}/fig1-8.png` and `fig_detrend.png`, `run_nb13.log` | notebook 13 with `data/boundaries/sc_mesoregions_2024.gpkg` (mesoregions with the islands, MapBiomas 2019 forest share), recomputed by notebook 19 after the LAI correction (`run_nb13.log`: local run of 8 November 2026) | `nb13` |
| `notebook14/` | TROPOMI HCHO evaluation, screening `cf0.3_sza70_fire120`: the four pre-specified tests (`tier1` to `tier4`), the monthly series (`monthly_series_*.csv`) and the summary figure | notebook 14 | `nb14` |
| `notebook15/` | `S3_fixlog_{annual,monthly}.csv` (test S3); `s1/` (test S1: class totals, groups, cell ratios, diurnal check, figure, log) | notebook 15; `analysis/v10/s1_compare.py` | `nb15`, `s1` |
| `diagnostics/` | checks defined after inspection, labelled as such in the documents: `c3/` (the C3 failure of 2015 and its effect on the totals), `nb14/` (post hoc diagnostics of the TROPOMI comparison, D1-D7), `check_run_integrity.log` | `analysis/v10/` | `c3`, `nb14diag`, `integrity` |
| `report/` | numbers of the inventory report: `report_numbers_v10.csv` (104 rows) and its tables (standardised regression, temperature sensitivity, mesoregions in 2025, ranked municipalities, isoprene monthly climatology, driver correlations) | `analysis/v10/report/10_report_numbers_v10.py` | `report` |
| `figures/` | `pt/` and `en/`: `fig_study_area.png` and `fig_verification.png` (article Figures 1 and 2, `analysis/v10/maps/study_area_map.py` and `analysis/v10/schematics/verification_flowchart.py`, with `study_area_shares.csv`), `fig_detrend.png` (article Figure S4, from notebook 13) and `fig1-9.png` (article Figures 3-11 and S5; files 1-4 and 6 copied from `notebook13_meso2024/figures/`; maps 5, 7 and 8 from `analysis/v10/maps/12_maps_wf_style.py`; 9 from `analysis/v10/fig9_tropomi.py`; chapter 9 of the tutorial gives the full map); `ijest/`: the journal versions `Fig1-Fig11.png` and `FigS4-FigS5.png` (`analysis/v10/ijest/ijest_figures.py`); `pt_report/`: report Figures 11 (`fig9.png`), 14 (`fig12.png`) and A.2-A.9 (maps) and 12 (`fig10.png`, standardised regression); `maps_class_limits.csv` (class limits of every map) | as stated | `nb13`, `maps`, `report`, `ijest` |
| `lai_viirs_check/` | post hoc check and correction of the 2024-2025 MODIS LAI (Online Resource 1, Section S2): tables and figure of notebook 16 (`lai_viirs_modis_*`), hindcast of the candidate corrections 2013-2023 (`lai_correction_hindcast_*.csv`), corrected LAI of 2024-2025 (`lai_viirs_corrected_2024_2025.npz`, read by notebook 17), monthly and annual summaries, `lai_correction.png`, `lai_correction_map_2025.png`, settings | notebook 16; `analysis/v10/review/lai_viirs_correction.py` | `laicorr` |
| `final_review/` | additional diagnostics (Benjamini-Hochberg control, Hamed-Rao check, start-year and 2025 sensitivity of the trends, LAI of 2024-2025, partial correlations, residual year term, parts of the trends (temperature, LAI, land cover), growth-form areas and emission potential 2003-2025, effect of CO₂ inhibition, background control and magnitude of the TROPOMI comparison, municipal statistics, pooled check C3, filled coastal cells): `final_review_numbers.csv` (key, value, note), `final_review.log`, `trends_fdr_hr.csv`, `partial_correlations_detrended.csv`, `lai_monthly_SC.csv` and its z-scores, `tropomi_magnitude_grid.csv`, `filled_cells_area.csv`, `landcover_ef_potential_SC.csv`; `trend_robustness.csv` and `.log` (2003-2023, trend-free pre-whitening, Benjamini-Yekutieli; Online Resource 1, Table S9); `month_boundary_effect.csv` and `.log` (monthly totals counted from 03 UTC, local midnight at UTC-3, instead of 00 UTC; `docs/TIME_CONVENTIONS.md`) | `analysis/v10/review/final_review_diagnostics.py`, `trend_robustness.py`, `month_boundary.py` | `final` |
| `ijest_review/` | diagnostics of the co-authors' and the external review: planted-forest emission factors (`plantation_ef_sensitivity_2019.csv`), LAI by growth form, ecotypes of the state and their tree emission factors (`ecotype_tree_ef.csv`, Online Resource 1, Table S8) | `analysis/v10/review/ijest_review_diagnostics.py`, `ecotype_tree_ef.py` | `ijest` |
| `notebook19/` | installation of the corrected 2024-2025: checks of the rerun, NOAA CO₂ file, comparison of both LAI inputs (`lai_sensitivity_summary.*`, `lai_sensitivity_table2_2024_2025.csv`), provenance | notebook 19 (Colab, 1 November 2026) | none |
| `lai_noaa20_check/` | LAI of the VIIRS on NOAA-20 (2018-2025) and NOAA-21 (from February 2023) against Suomi NPP; corrected LAI from each record (Online Resource 1, Section S2.5.6, Table S7) | notebook 22 (Colab, 8 November 2026) | none |
| `brain_comparison/` | 2019 against BRAIN (Hoinaski et al., 2023): classes, CB6X species, months, mesoregions, cells, file checks; rates as delivered and on true cell areas (Online Resource 1, Section S9, Table S12) | notebook 18 (Colab, 1 November 2026) | none |
| `brain_differences/` | tests of the explanations of the differences from BRAIN (blocks A-H: default switches, class fingerprint, cell regression, seasonal decomposition, terrain, sensitivity runs, TROPOMI 2019, soil NO) | notebook 23 (Colab, 3 November 2026); `E_terrain_with_interval.csv` and `F_lapse_vs_terrain_full.csv` from a local rerun of the updated notebook on the same inputs | none |
| `sensitivity_2019/` | one-at-a-time reruns of 2019 (`full`: 365 days; `design`: 28 days): PAR 0.5, terrain-adjusted temperature, tropical canopy, needleleaf plantations, clay wilting point (Online Resource 1, Table S10) | notebook 20 (Colab, 2 November 2026) | none |
| `era5land_vs_inmet/` | ERA5-Land against the INMET automatic stations, 2019: temperature bias, daily radiation, diurnal lag, station-fitted elevation coefficient and its influence diagnostics (Online Resource 1, Table S11) | notebook 21 (Colab, 2 November 2026); `lapse_rate_influence.csv` and `bias_summary_2019.csv` from a local rerun of its updated cell 8 on the same station table | none |

The report reuses article Figure 1 as its Figure 1, article Figures 2-9 and S5 as its Figures 2-10, the combined TROPOMI
figure (article Figures 10 and 11) as its Figure A.10 and article Figure S4 as its Figure A.11. Report Figure 13
(annual temperature against soil moisture) and Figure A.1 (growth-form fractions, 2024) are not redrawn by the scripts
of the repository; their data are in `notebook13_meso2024/drivers_annual_SC.csv` (key `r_T_SOILM` of
`report/report_numbers_v10.csv`: r = −0.22, p = 0.30) and in `data/run_v10/2025/inputs/growthform_2024.tif`.

Chapter 8 of the tutorial (`docs/tutorial/08_number_map.md`) maps each number of both documents to a file and key of
this folder, and chapter 9 each figure and table to its file and script.

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
