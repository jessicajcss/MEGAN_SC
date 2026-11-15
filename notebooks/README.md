# Notebooks

Fourteen notebooks produce, check and analyse the inventory. Each one is built from the percent-format source with the same
name in `src/` (`python src/build_ipynb.py src/<name>.py`), so the `.py` and the `.ipynb` hold the same code. They were
written for Google Colab with Google Drive; 11, 13, 14 (with its cache) and 15 also run on any computer through
`tools/run_notebook.py`, which is what `reproduce.py` does.

| Order | Notebook | Where | Needs | Produces | Published copy |
|---|---|---|---|---|---|
| 1 | `12_fortran_benchmark.ipynb` | Colab (Linux x86-64), 15-25 min | `MEGANv3.21.tar.gz` (Zenodo 10.5281/zenodo.10526206), this repository, the 2019 inputs for the Santa Catarina days | `megan_sc_benchmark_results/` (CSV + `benchmark_summary.csv`) | `benchmark/results_colab/` |
| 2 | `10_production_run.ipynb` | Colab + Earth Engine; hours per model year, resumable | this repository in Drive, an Earth Engine project | `megan_sc_outputs_v10/<year>/...` | `data/run_v10/` (release archive) |
| 3 | `11_verify_run.ipynb` | Colab or local | output of 10 | `verify_v10_report.csv` (checks C0-C5), annual totals by class, a summary figure | `results/v10/notebook11/` |
| 4 | `13_manuscript_numbers_and_figures.ipynb` | Colab or local | output of 10, `data/boundaries` | Tables 2-4 (report numbering; article Tables S1, 2 and 3), `manuscript_numbers.csv`, seasonal cycle, files `fig1-8.png` (article Figures 3-9 and S5; map in `docs/tutorial/09_figures_and_tables.md`) and `fig_detrend.png` (article Figure S4), pt and en, `annual_grids_cache.npz` | `results/v10/notebook13_meso2024/` |
| 5 | `14_tropomi_hcho_evaluation.ipynb` | Colab + Earth Engine (41.2 EECU-hours and about 6 h for the 84 months); local with the cache | output of 10 for 2019-2025 | the four tests, the monthly series, the monthly HCHO fields (`cache_<TAG>/`) | `results/v10/notebook14/`, `data/tropomi_hcho_v10/` |
| 6 | `15_sensitivity_2019.ipynb` | Colab or local | output of 10 (fix logs); tests S1 and S2 need extra 2019 runs of notebook 10 | test S3 (S1 and S2 when their runs are present) | `results/v10/notebook15/` |
| 7 | `16_lai_independent_check.ipynb` | Colab + Earth Engine; 312 requests (2 per month and product) | output of 10 (grid and growth forms), `data/boundaries`, an Earth Engine project | `lai_viirs_modis_monthly.csv`, `lai_viirs_modis_annual.csv`, `lai_viirs_modis_summary.csv` and `lai_viirs_modis.png` in `<ROOT>/analysis_v10/lai_viirs_check/` | `data/lai_viirs_v10/`, `results/v10/lai_viirs_check/` |
| 8 | `17_rerun_2024_2025_lai_viirs.ipynb` | Colab + Earth Engine; several hours per year, resumable | output of 10 (LAI 2013-2025, growth forms), the cache of 16, an Earth Engine project | `megan_sc_outputs_v10_lai_viirs/{2024,2025}/` (layout of notebook 10), `lai_viirs_rerun_comparison.csv`; optionally installs the two years in the production run | `data/run_v10/{2024,2025}` (installed by 19); results with the MODIS-LAI years in `results/v10_modis_lai_2024_2025/` |
| 9 | `18_brain_comparison_2019.ipynb` | Colab or local; downloads the 364 daily BRAIN MEGAN files of 2019 from the Science Data Bank (one zip per day, about 190 MB unzipped), resumable by month; no Earth Engine | output of 10 for 2019, `data/boundaries`, `MEGANv3.21.tar.gz`, the list of BRAIN download links (`BRAIN_output_urls.txt`) | `analysis_v10/brain_comparison/`: annual, monthly, mesoregion and cell comparisons, speciation checks, parameter table, figures, monthly BRAIN grids | `results/v10/brain_comparison/` |
| 10 | `19_install_and_reproduce.ipynb` | Colab, 30-60 min; no Earth Engine | outputs of 10, 16 and 17 in Drive, the release archives of the run, this repository | checks of the run of 17 (cell 7a included), installation of 2024-2025 with the corrected LAI (MODIS years to `megan_sc_outputs_v10_modis_lai_2024_2025`), the official NOAA CO₂ file, `reproduce.py` three times (MODIS control, CO₂ file only, corrected LAI), `lai_sensitivity_summary.{csv,md}`; optionally `results/v10` with the corrected LAI, checksums, manifest and the five release zips | `results/v10/notebook19/` |
| 12 | `20_sensitivity_2019_inputs.ipynb` | Colab + Earth Engine; about 40 min (28 days) or about 15 h (full year; 914 min in the authors' run), resumable | output of 10 for 2019, an Earth Engine project, optionally the lapse rate of 21 | one-at-a-time runs of 2019: `base`, `lapse` (temperature at the terrain of each 0.025° cell), `par050`, `canopy_tropical`, `canopy_needle_silv`, `wilt_clay` (`analysis_v10/sensitivity_2019/`) | `results/v10/sensitivity_2019/` |
| 11 | `21_era5land_vs_inmet.ipynb` | Colab + Earth Engine; about 30 min | the INMET station lists (1981-2010, 1961-1990); downloads the INMET 1991-2020 normals and the 2019 hourly station archive | ERA5-Land minus station temperature (normals, 2019 hourly, by month and hour), daily radiation, empirical lapse rate `lapse_rate.json` read by 20 (`analysis_v10/era5land_vs_inmet/`) | `results/v10/era5land_vs_inmet/` |
| 13 | `22_lai_check_noaa20_viirs.ipynb` | Colab, 1-2 h; Earthdata account; no Earth Engine | output of 10 (MODIS LAI, growth forms), the cache of 16 | NOAA-20 VIIRS LAI (VJ115A2H v2) and Suomi NPP (VNP15A2H v2) read from the HDF5 files, 2018-2025, and NOAA-21 (VJ215A2H v2) from February 2023; the corrected 2024-2025 LAI recomputed from each VIIRS record (`analysis_v10/lai_noaa20_check/`) | `results/v10/lai_noaa20_check/` |
| 14 | `23_brain_differences_2019.ipynb` | Colab or local, about 2 min; no Earth Engine | the state of 18 (re-run from its cache), optionally the outputs of 14, 20 and 21 | tests of the explanations of the differences from BRAIN: MEGANv3.21 default switches, class fingerprint, cell regression, seasonal decomposition, terrain, sensitivity runs, TROPOMI 2019, soil NO (`analysis_v10/brain_differences/`) | `results/v10/brain_differences/` |

Notebook 12 does not depend on the others. Notebook 13 should run after notebook 11 reports zero FAIL rows, or after
every FAIL has been explained (the production run has one, C3 in 2015 for methanol, analysed in
`results/v10/diagnostics/c3/`).

Notebook 16 compares the monthly VIIRS LAI of Suomi NPP (`NASA/VIIRS/002/VNP15A2H`, main-algorithm retrievals) with
the MODIS LAI of the production run on the inventory grid, 2013-2025. Run by the authors on 31 October 2026, it shows
that VIIRS does not show the 2024-2025 decrease of the MODIS LAI (2025 statewide: MODIS −12.1 %, VIIRS −1.8 % of the
2013-2023 mean), so the decrease is treated as a defect of the recent MODIS record. Notebook 17 computes 2024 and 2025
with a corrected LAI (2013-2023 MODIS mean of each cell and month times the VIIRS anomaly ratio of its growth-form
zone; `src/megan_sc/laicorr.py`, chosen by the hindcast of `analysis/v10/review/lai_viirs_correction.py`). Everything
else is identical to notebook 10, and its cell 11 checks that only the LAI differs.

Notebook 18 compares the 2019 emissions with the MEGAN emissions of BRAIN, domain d02 (Hoinaski et al. 2024;
doi:10.57760/sciencedb.09886): same region (the IBGE polygon,
area-weighted on both grids), same species (the 19 classes converted to the 34 CB6X species of BRAIN with the
MEGANv3.21 tables `SPC_NOCONVER.EXT` and `MAP_CV2CB6X.EXT`), same units (moles and tonnes) and same period (2019, UTC
months), with this inventory remapped conservatively onto the 4-km Mercator grid of BRAIN for the maps and cell
statistics. It checks the BRAIN files (time steps, grid, units, species never produced by the CB6X map, speciation
fingerprints) and the co-registration of the two grids, and reports the BRAIN state totals both as delivered and for the
true cell areas, because MGN2MECH multiplies the fluxes by the nominal cell area (`garea = xcell3d * ycell3d`) of a
Mercator grid true at the equator. Tested on 15 January 2020, the BRAIN file of the same series kept in the project, and
run by the authors on 1 November 2026 on the 364 daily files of 2019 that the BRAIN list holds (31 December is missing;
`results/v10/brain_comparison/provenance_brain_comparison.json`).

Notebooks 19-23 were run by the authors in Colab between 1 and 8 November 2026 and checked against their outputs. Notebook
19 checks the run of 17, installs it, takes the official CO₂ file from NOAA GML and runs `reproduce.py` with the MODIS
LAI (which gives back `results/v10_modis_lai_2024_2025`) and with the corrected LAI, side by side. Notebook 20 runs 2019 again
with one input changed at a time, notebook 21 measures the ERA5-Land temperature and radiation against the INMET
stations, notebook 22 repeats the LAI check with the VIIRS of NOAA-20 and NOAA-21, and notebook 23 tests each explanation of the
differences from BRAIN with the data at hand (it uses 14, 20 and 21 when their outputs exist).

Notebook 21 runs before notebook 20: the `lapse` run of 20 reads the lapse rate of 21 (`lapse_rate.json`, 10.78 K km⁻¹
in the published run) and falls back to 6.5 K km⁻¹ when that file is missing. Notebook 19 was run once by the authors:
it starts from the run folders and archives as they were before the installation, which are not distributed, and its
last cells write the checksum lists of the installed run (`data/RELEASE_ASSETS_SHA256SUMS.txt`,
`data/DATA_SHA256SUMS.txt`). Its comparisons are repeated with `reproduce.py` (REPRODUCE.md, section B).

Each notebook opens in Google Colab at
`https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/<notebook>.ipynb` (repository
https://github.com/jessicajcss/MEGAN_SC, branch `main`).

## Settings

Each notebook has one settings cell (the first code cell); notebooks 10, 17 and 20 set the Earth Engine project in their
cell "3 · Earth Engine" and notebook 21 in its cell "2 · Earth Engine". The defaults assume this Google Drive layout:
`/content/drive/MyDrive/megan-sc` (this repository) and `/content/drive/MyDrive/megan_sc_outputs_v10` (the output of
notebook 10). Set your own Earth Engine project in notebooks 10, 14, 16, 17, 20 and 21 (`EE_PROJECT = 'your-ee-project'`). Outside
Colab, `tools/run_notebook.py` replaces the settings after the first cell:

```bash
python tools/run_notebook.py notebooks/13_manuscript_numbers_and_figures.ipynb \
    "ROOT='data/run_v10';BND='data/boundaries';OUT_DIR='reproduced/notebook13_meso2024'"
```

Notebook 10 writes `src/megan_sc/fixes.py` (the input and run procedures) into the repository folder before importing
it; the file it writes is the one committed here (sha256 prefix a933b29602929a82), and every output file carries that
prefix. If `OLD_ROOT` points to a folder that already holds the MapBiomas growth forms of a year
(`<year>/inputs/growthform_<year>.npz`), notebook 10 reuses them instead of requesting them again; otherwise it computes
them in Earth Engine.

## Variants of notebook 10

`VARIANT = 'base'` is the production run. `VARIANT = 'hourly'` (hourly sampling) and `VARIANT = 'lai_eedefault'`
(Earth Engine default LAI resampling) are the sensitivity runs of notebook 15; `TEST_DAYS = 7` limits each month to its
first seven days. Test S1 used `TEST_DAYS = 7` with `base` and `hourly` on January, April, July and October 2019
(`data/s1_test/`); test S2 was not run.
