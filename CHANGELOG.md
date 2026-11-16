# Changelog

## 1.0.0 (first release)

* Engine: Python implementation of the MEGANv3.21 emission code (canopy, activity factors, soil NO, speciation) with
  the input and run procedures FIX-1 to FIX-8 of `src/megan_sc/fixes.py` (SHA-256 prefix a933b29602929a82, stamped in
  every output file), compared stage by stage with the official Fortran executables (`VALIDATION.md`, notebook 12).
* Any region and period: a run is defined by a configuration file (`config/*.yaml`). Notebook 10 runs the model for
  it (MapBiomas or MODIS MCD12Q1 land cover, MapBiomas collection fixed by the configuration), notebook 11 checks the
  run and notebook 25 writes its standard results (totals, seasonal cycle, maps of intensity, drivers, correlations,
  trends, zone totals). The region wizard (`run_your_region/`) writes the configuration of a new region from its
  boundary file and tests it. Every notebook finds the repository by itself (`MyDrive/MEGAN_SC` in Colab, the
  repository folder elsewhere), and Earth Engine works in Colab, on any computer and with a service account
  (`megan_sc.gee.initialize`).
* Worked example: Santa Catarina, Brazil, 2003-2025, 0.025°, every 3 h (`config/santa_catarina_2003_2025.yaml`), with
  2024 and 2025 computed with the VIIRS-corrected LAI (notebooks 16, 17 and 19); the run and the MODIS-LAI years are
  distributed as release archives (`DATA_AVAILABILITY.md`). Its results are in `results/v10`, and the sensitivity case
  with the MODIS LAI of 2024-2025 in `results/v10_modis_lai_2024_2025`; checksums in `results/MANIFEST.sha256`. The
  land cover of reference is MapBiomas Collection 10 of 2024 (study-area figure, native-forest share of the
  mesoregions); the analyses of 2019 (comparisons with BRAIN and with the previous inventory, notebooks 20 and 23, the
  plantation diagnostic) keep MapBiomas 2019, the year they analyse.
* Checks and evaluations of the example: notebook 11 (checks C0-C5), TROPOMI formaldehyde (notebook 14), time step and
  soil moisture (notebook 15), 2019 sensitivity runs (notebook 20), ERA5-Land against INMET (notebook 21), VIIRS LAI of
  three satellites (notebooks 16 and 22), comparison with BRAIN 2019 (notebooks 18 and 23), and the time conventions of
  the inputs and of the model (`docs/TIME_CONVENTIONS.md`).
* `reproduce.py`: every number, table and figure of the example recomputed from the release archives and compared with
  the published copies (`docs/tutorial/reproduction_log.txt`).
* Documentation: illustrated tutorials in Portuguese and English (`TUTORIAL_MEGAN_SC_pt.pdf`,
  `TUTORIAL_MEGAN_SC_en.pdf`) and 12 technical chapters (`docs/tutorial/`), from the MEGAN equations to Google Drive,
  Earth Engine, your own region and the reading of the results; unit tests, including an offline run of notebooks 10,
  11 and 25 without Earth Engine (`tools/offline_run.py`).
* Double-blind review: author names, e-mail, ORCID and affiliations are withheld in `CITATION.cff`, and the field
  `ee_project` of the provenance files of the release archives reads `"withheld"` (`tools/anonymize_release.py`).

### Changes since the first push of the repository (commit 9df851d)

For a clone made before these changes:

* renamed: `notebooks/13_manuscript_numbers_and_figures.ipynb` to `13_santa_catarina_analysis.ipynb`; the table
  `manuscript_numbers.csv` of notebook 13 to `key_numbers.csv`; the step `ijest` of `reproduce.py` to `posthoc`, with
  `analysis/v10/review/ijest_review_diagnostics.py` to `posthoc_review_diagnostics.py` and the folder
  `results/*/ijest_review` to `posthoc_review` (files `posthoc_review.log`, `posthoc_review_numbers.csv`); the tutorial
  chapters 5, 6, 8, 11 and 12 (`05_reproduce_the_example.md`, `06_run_the_model.md`, `08_your_region.md`,
  `11_troubleshooting.md`, `12_provenance.md`);
* new: `config/santa_catarina_2003_2025.yaml`, notebook 25, `notebooks/src/_setup_repo.py`, `tools/offline_run.py`,
  `tests/test_offline_run.py`, tutorial chapter 9 (reading the results), the setting `MONTHS` of notebook 10;
* removed: the figures in journal format (`analysis/v10/ijest/`, `results/*/figures/ijest/`) and the tutorial chapters
  that mapped the results to an external document (number map, figures and tables);
* changed: notebooks 10 and 11 read a run configuration (for the example the results are unchanged); the region wizard
  ends by writing the configuration and leaves the full run to notebooks 10, 11 and 25; `.gitattributes` keeps the
  line endings of the MEGAN tables (`data/megan_inputs/** -text`), so their checksums match on every platform;
  `.gitignore` adds `data/regions/`, `megan_sc_outputs_*/` and `.Rproj.user/`; MapBiomas 2024 replaces 2019 as the
  land cover of reference (study-area figure, `forest_pct` of the mesoregions: Pearson r 0.968 to 0.971, 95 % interval
  0.73-0.997 to 0.75-0.997; Spearman ρ 0.94 and exact p = 0.017 unchanged);
* changed in `results/` (labels only, numbers unchanged): the notes of `figures/pt_report/figA2.png` to `figA8.png` cite
  FIX-1 of `src/megan_sc/fixes.py` and notebooks 16 and 17 instead of sections of an external document (redrawn in both
  results folders; the LAI note of the MODIS-LAI case now says that its LAI is uncorrected), the sources of
  `report/report_numbers_v10.csv` read `key_numbers.csv`, the sources of `brain_comparison/parametrizations.csv` cite
  notebooks and Hoinaski et al. (2024), and `notebook22_printed_checks.txt` no longer cites an external document;
  `results/MANIFEST.sha256` lists 337 files;
* fixed: the scripts that run cells of notebooks 11 and 14 (`c3/diag_c3.py`, `c3/diag_c3_lai.py`, `nb14_diagnostics.py`,
  `fig9_tropomi.py`) read the notebook source through `_paths.notebook_source`, which replaces the block that finds the
  repository; the reproduction of 16 November 2026 is in `docs/tutorial/reproduction_log.txt`.
