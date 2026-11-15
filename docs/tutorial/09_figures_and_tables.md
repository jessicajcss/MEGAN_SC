# 9. Figures and tables

Every figure of both documents is a PNG of `results/v10/figures/` (or of `notebook13_meso2024/figures/`), written by
one script or notebook; every table is a CSV of `results/v10/`. The `pt/` folders hold the Portuguese labels, `en/`
the English ones, and `pt_report/` the figures that only the report uses.

## Article

The numbers are those of the article. The journal versions of the figures (no titles or source
notes inside the figures, Gg units) are `figures/ijest/Fig1`-`Fig11` and `FigS4`-`FigS5` (PNG; the script also writes
LZW-compressed TIFF copies in `reproduced/figures/ijest/`, not kept in `results/`), drawn by `analysis/v10/ijest/ijest_figures.py` (step `ijest`) from the same data as the language versions below. The
language versions keep the order in which notebook 13 writes them, so the number in a file name is not the number of
the figure in the article.

| Item | Content | File (language versions) | Written by | Step |
|---|---|---|---|---|
| Figure 1 | study area: Santa Catarina in South America, MapBiomas land use and land cover in 2019, tree fraction on the grid and growth-form shares | `figures/{pt,en}/fig_study_area.png`; shares in `figures/study_area_shares.csv` | `analysis/v10/maps/study_area_map.py` | `extra` |
| Figure 2 | verification and evaluation steps as a timeline: how each one was tested, its result, what it means and its status | `figures/{pt,en}/fig_verification.png` | `analysis/v10/schematics/verification_flowchart.py` | `extra` |
| Figure 3 | annual series 2003-2025 of BVOC, isoprene, monoterpenes, sesquiterpenes and soil NO | `figures/{pt,en}/fig1.png` | notebook 13 | `nb13` |
| Figure 4 | BVOC by chemical family: a) mean composition, b) seasonal cycle of the monthly emissions, c) cycle of each family normalised by its monthly mean | `figures/{pt,en}/fig2.png` | notebook 13 | `nb13` |
| Figure 5 | detrended correlation of each class with the drivers | `figures/{pt,en}/fig3.png` | notebook 13 | `nb13` |
| Figure 6 | SOA precursors by mesoregion, share and intensity | `figures/{pt,en}/fig4.png` | notebook 13 | `nb13` |
| Figure 7 | map of SOA-precursor intensity by mesoregion | `figures/{pt,en}/fig5.png` | `analysis/v10/maps/12_maps_wf_style.py` | `maps` |
| Figure 8 | municipal BVOC in 2025, total and intensity (Jenks natural breaks) | `figures/{pt,en}/fig7.png` | `maps/12_maps_wf_style.py` | `maps` |
| Figure 9 | municipal isoprene, CO and soil NO in 2019 and 2025 (Jenks natural breaks, common to both years; legend on the left of each row) | `figures/{pt,en}/fig8.png` | `maps/12_maps_wf_style.py` | `maps` |
| Figures 10 and 11 | TROPOMI: spatial (summer maps of isoprene and ΔHCHO, surface pressure, ΔHCHO against pressure) and temporal (seasonal cycles, anomalies) comparison; one combined figure in the language versions | `figures/{pt,en}/fig9.png` | `analysis/v10/fig9_tropomi.py` | `maps` |
| Figure S4 (Online Resource 1) | example of the correlation between detrended series: soil NO and temperature (r = 0.30 raw, 0.82 detrended) | `figures/{pt,en}/fig_detrend.png` | notebook 13 | `nb13` |
| Figure S5 (Online Resource 1) | intensity against native-forest share (n = 6) | `figures/{pt,en}/fig6.png` | notebook 13 | `nb13` |
| Figure S1 (Online Resource 1) | MODIS and VIIRS LAI, 2013-2025 | `lai_viirs_check/lai_viirs_modis.png` | notebook 16 | none (published copy of notebook 16) |
| Figures S2 and S3 (Online Resource 1) | correction of the 2024-2025 LAI, and its 2025 map | `lai_viirs_check/lai_correction.png`, `lai_correction_map_2025.png` | `analysis/v10/review/lai_viirs_correction.py` | `laicorr` |
| Table 1 | MEGAN classes, LDF and β | | `src/megan_sc/engine/species.py` (constants of MEGVEA.EXT) | none |
| Table 2 | SOA precursors by mesoregion | `notebook13_meso2024/table3_mesoregions_soa.csv` | notebook 13 | `nb13` |
| Table 3 | the twelve largest municipal totals of 2025 | `notebook13_meso2024/table4_municipalities_2025.csv` | notebook 13 | `nb13` |
| Table 4 | configuration of this inventory and of the previous one | text: `docs/ALGORITHM_SPEC.md` (section 11) and Hoinaski et al. (2022) | | none |
| Table 5 | 2019 totals of both inventories | `manuscript_numbers.csv` (`cmp2019_*`) and the 2019 row of Table S1 | notebook 13 | `nb13` |
| Table S1 (Online Resource 1) | annual totals by group, 2003-2025 | `notebook13_meso2024/table2_annual_totals_SC_groups.csv` | notebook 13 | `nb13` |
| Tables S2 and S3 | LAI inputs of the inventories and of BRAIN; LAI products compared | text (chapter 2 and the product documentation) | | none |
| Table S4 | 2024 and 2025 LAI against the 2013-2023 mean, by cell group | `lai_viirs_check/lai_viirs_modis_summary.csv` (notebook 16), `lai_correction_annual_groups.csv` | notebook 16; `lai_viirs_correction.py` | `laicorr` |
| Table S5 | hindcast of the candidate corrections | `lai_viirs_check/lai_correction_hindcast_2013_2023.csv`, `lai_correction_hindcast_annual.csv` | `lai_viirs_correction.py` | `laicorr` |
| Tables S6-S12 | 2024-2025 with both LAI inputs, three VIIRS satellites, ecotypes, trend robustness, 2019 sensitivity, ERA5-Land against INMET, BRAIN | see chapter 8, section "Online Resource 1" | notebooks 18-23; `trend_robustness.py`; `ecotype_tree_ef.py` | `final`, `ijest` |

The files `fig1-4.png`, `fig6.png` and `fig_detrend.png` are copies of those in `notebook13_meso2024/figures/` (pixel-identical). Notebook 13 draws its
own versions of `fig5`, `fig7` and `fig8` there too; the documents use the redrawn maps of `figures/`.

## Report

| Item | Content | File | Written by | Step |
|---|---|---|---|---|
| Figure 1 | the article Figure 1 (Portuguese) | `figures/pt/fig_study_area.png` | as above | |
| Figures 2-10 | the article Figures 2-9 and S5 (Portuguese) | `figures/pt/fig_verification.png`, `fig1-8.png` | as above | |
| Figure 11 | emissions of 2025 by group of compounds (kg ha⁻¹ yr⁻¹) | `figures/pt_report/fig9.png` | `maps/12_maps_wf_style.py` | `maps` |
| Figure 12 | standardised regression coefficients and R² | `figures/pt_report/fig10.png` | `analysis/v10/report/11_report_figures_v10.py` | `report` |
| Figure 13 | annual temperature against soil water | not redrawn by the repository scripts (data: `notebook13_meso2024/drivers_annual_SC.csv`; r in `report_numbers_v10.csv`, `r_T_SOILM`) | | none |
| Figure 14 | tree fraction and SOA-precursor intensity, 2025 | `figures/pt_report/fig12.png` | `maps/12_maps_wf_style.py` | `maps` |
| Figure A.1 | growth-form fractions, MapBiomas 2024 | not redrawn by the repository scripts (data: `data/run_v10/2025/inputs/growthform_2024.tif`) | | none |
| Figures A.2-A.8 | monthly maps of the inputs in 2025: LAI (corrected with VIIRS), temperature, PAR, soil water, wind, humidity, pressure | `figures/pt_report/figA2.png` to `figA8.png` | `maps/12_maps_wf_style.py` | `maps` |
| Figure A.9 | BVOC intensity of every year 2003-2025 | `figures/pt_report/figA9.png` | `maps/12_maps_wf_style.py` | `maps` |
| Figure A.10 | the article Figures 10 and 11, combined | `figures/pt/fig9.png` | `fig9_tropomi.py` | `maps` |
| Figure A.11 | the article Figure S4 (example of the detrended correlation) | `figures/pt/fig_detrend.png` | notebook 13 | `nb13` |
| Tables 2, 3, 4 and 6 | article Tables S1, 2, 3 and 5 | | notebook 13 | `nb13` |
| Table 7 | drivers, their role in MEGAN and the observed effect | text with the numbers of `report_numbers_v10.csv` (chapter 8) | | |
| Quadros A.1-A.3 | ERA5-Land variables, LAI product, MapBiomas crosswalk | text (chapters 1 and 2) | | |

Figure 13 and Figure A.1 of the report are not redrawn by the repository scripts. Their data are the state means of
temperature and soil water (`notebook13_meso2024/drivers_annual_SC.csv`) and the MapBiomas growth forms of the run
(`data/run_v10/2025/inputs/growthform_2024.tif`).

## Map conventions

The maps share one layout (`analysis/v10/maps/wf_style.py`), that of the companion biomass-burning article: a
geographic frame or a grid of panels, north arrow, scale bar inside the map frame at the bottom right (next to the last map when a
figure has several maps and no frame), source note and class legends, with
the colour ramps of the biogenic figures. The municipal maps (article Figures 8 and 9) use Jenks natural breaks (Fisher's exact
algorithm, five classes; in article Figure 9 common to 2019 and 2025 for each compound); the grid maps of the report (Figures 11,
A.9 and the intensity panel of Figure 14) use quintiles of the state cells; the mesoregion map (article Figure 7), the tree
fraction of report Figure 14 and the input maps of the appendix use fixed intervals. All class limits are written to
`figures/maps_class_limits.csv` (the article maps appear twice, once per language), so a map can be redrawn with the
same classes in any GIS program.
