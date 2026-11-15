# Scripts of the analysis of the production run

These scripts compute the diagnostics, the report numbers, the maps and the additional diagnostics from the production
run. Each one reads its folders from `_paths.py` and writes under `reproduced/` (or `$MEGAN_SC_OUT`); the published
copies of their outputs are in `results/v10/`. `python reproduce.py <step>` runs them with the right environment;
they can also be run directly from any folder (`python analysis/v10/<script>.py`).

## Paths (`_paths.py`)

| Variable | Default | Content |
|---|---|---|
| `MEGAN_SC_RUN` | `data/run_v10` | production run of notebook 10, one folder per year |
| `MEGAN_SC_BND` | `data/boundaries` | IBGE boundaries |
| `MEGAN_SC_TROPOMI` | `data/tropomi_hcho_v10` | monthly TROPOMI fields cached by notebook 14 |
| `MEGAN_SC_S1` | `data/s1_test` | the two test runs of test S1 |
| `MEGAN_SC_LAIVIIRS` | `data/lai_viirs_v10` | monthly MODIS and VIIRS LAI on the grid, 2013-2025 (cache of notebook 16) |
| `MEGAN_SC_RUN_MODIS_ARCHIVE` | `data/run_v10_modis_lai_2024_2025` | 2024 and 2025 computed with the MODIS LAI |
| `MEGAN_SC_RUN_MODIS` | `data/run_v10_modis_view` | the MODIS-LAI record of every year (view built by `tools/lai_variant_view.py`) |
| `MEGAN_SC_RESULTS` | `results/v10` | published results (the reference of `reproduce.py compare`) |
| `MEGAN_SC_OUT` | `reproduced` | where the recomputed files go |
| `MEGAN_SC_NB13` | `<OUT>/notebook13_meso2024` if notebook 13 has been rerun there, else `results/v10/notebook13_meso2024` | tables and the annual grid cache (`annual_grids_cache.npz`) of notebook 13 |

`python reproduce.py paths` prints the values in use. The scripts of `report/`, `maps/` and `s1_compare.py` also import
`common.py` of this folder (GeoPackage reader, model grid read from `data/run_v10`, state mask, cell areas, polygon
rasterisation, area weights of the cells in the polygons, and the factor ρ of isoprene).

## Scripts

| Script | What it does | Writes (under `reproduced/`) | Step |
|---|---|---|---|
| `check_run_integrity.py` | eight steps per day in each of the 276 months; annual files equal to the sum of the monthly ones | log only | `integrity` |
| `s1_compare.py` | test S1: 3-hourly against hourly sampling on the two 7-day test runs | `notebook15/s1/` | `s1` |
| `c3/diag_c3.py`, `c3/diag_c3_lai.py` | check C3 of notebook 11 cell by cell, and its LAI-adjusted version (defined after inspection) | `diagnostics/c3/` | `c3` |
| `c3/c3_impact_on_totals.py` | change of the state totals if the flagged, or all comparable, filled cells took their neighbours' median activity | log | `c3` |
| `nb14_diagnostics.py` | post hoc diagnostics of the TROPOMI comparison (D1-D7: split halves, surface pressure, magnitude, background, screening) | `diagnostics/nb14/` | `nb14diag` |
| `report/10_report_numbers_v10.py` | every number of the inventory report | `report/` | `report` |
| `report/11_report_figures_v10.py` | report Figure 12 (standardised regression) | `figures/pt_report/fig10.png` | `report` |
| `maps/12_maps_wf_style.py` (with `maps/wf_style.py`) | every map of the article (Figures 7, 8 and 9; files fig5, fig7, fig8, redrawn in the journal format by `ijest/ijest_figures.py`) and of the report (Figures 7, 9, 10, 11 and 14, files fig5, fig7, fig8, fig9 and fig12; Figures A.2-A.9, files figA2-figA9), with the class limits | `figures/{pt,en,pt_report}/`, `figures/maps_class_limits.csv` | `maps` |
| `fig9_tropomi.py` | the TROPOMI figure (file fig9; report Figure A.10), Portuguese and English labels; its panels are article Figures 10 and 11, drawn in the journal format by `ijest/ijest_figures.py` | `figures/{pt,en}/fig9.png` | `maps` |
| `meso2024/make_sc_mesoregions_2024.py` | mesoregions as unions of the 2024 municipalities (coastal islands included), with the MapBiomas 2019 forest share of `data/boundaries/landcover2019_mesoregion.csv`, and the municipality-mesoregion table; writes to `reproduced/boundaries/` unless `MEGAN_SC_BND_OUT` is set | `boundaries/` | `meso2024` |
| `review/final_review_diagnostics.py` | additional diagnostics: Benjamini-Hochberg control, Hamed-Rao check, start year and 2025 sensitivity of the trends, LAI of 2024-2025, partial correlations, residual year term, background control and magnitude of the TROPOMI comparison, mesoregion and municipal statistics, pooled check C3, filled coastal cells | `final_review/` | `final` |
| `review/ijest_review_diagnostics.py` | two post hoc diagnostics of the IJEST manuscript: effect of the emission factors of *Pinus* or *Eucalyptus* on the 2019 plantation area (MapBiomas class 9) and LAI of 2024-2025 by dominant growth form | `ijest_review/` | `ijest` |
| `ijest/ijest_figures.py` | Figures 1-11 of the English manuscript for IJEST and Figures S4-S5 of its Online Resource 1, from the same tables and grids as the published figures: no titles or source notes inside, 174 mm (Fig. 5: 84 mm) wide, 600 dpi, PNG and LZW-compressed TIFF, TROPOMI figure split into Figs. 10 and 11 | `figures/ijest/` | `ijest` |
| `review/ecotype_tree_ef.py` | MEGAN v3 ecotypes of the state, their tree composition and the composed tree emission factors (Online Resource 1, Table S8) | `ijest_review/ecotype_tree_ef.csv` | `ijest` |
| `review/trend_robustness.py` | robustness of the trends of the nine groups: 2003-2023, trend-free pre-whitening, Benjamini-Yekutieli, the eight groups without total BVOC (Online Resource 1, Table S9; post hoc) | `final_review/trend_robustness.{csv,log}` | `final` |
| `maps/study_area_map.py` | study area: MapBiomas Collection 10 in 2019 and the model growth forms (article and report Figure 1) | `figures/{pt,en}/fig_study_area.png`, `figures/study_area_shares.csv` | `extra` |
| `schematics/verification_flowchart.py` | verification and observational evaluation (article and report Figure 2), every number read from the result files | `figures/{pt,en}/fig_verification.png` | `extra` |
| `review/lai_viirs_correction.py` | VIIRS check of the 2024-2025 MODIS LAI: hindcast of the candidate corrections and the corrected LAI of 2024-2025 (Online Resource 1, Section S2, Figures S2-S3, Tables S4-S5) | `lai_viirs_check/lai_correction*`, `lai_viirs_corrected_2024_2025.npz` | `laicorr` |

Order matters in two places: the scripts of `report/`, `maps/` and `review/` read the outputs of notebook 13 (run
`reproduce.py nb13` first, or they use the published copies), and `report/11_report_figures_v10.py` reads the table
written by `report/10_report_numbers_v10.py`.
