# MEGAN_SC: biogenic emissions of Santa Catarina, Brazil, 2003-2025

Code, model output and analysis behind the article *Biogenic emissions in Santa Catarina, Brazil,
2003-2025: a MEGANv3.21 inventory driven by MapBiomas, MODIS and ERA5-Land* and the state inventory report
(*Inventário das emissões biogênicas na qualidade do ar em Santa Catarina*, base year 2025 with the 2003-2025 series).
Both documents cite this repository at https://github.com/jessicajcss/MEGAN_SC.

> [Leia-me em português](README.pt-BR.md) · **Illustrated PDF tutorial: [`TUTORIAL_MEGAN_SC_en.pdf`](TUTORIAL_MEGAN_SC_en.pdf)** ([português](TUTORIAL_MEGAN_SC_pt.pdf)) · [Technical chapters](docs/tutorial/README.md) · [Reproduce in four commands](#quick-start-reproduce-every-number-in-a-few-minutes)

MEGAN_SC is a Python implementation of the MEGANv3.21 emission code (Guenther, 2024b) with a Google Earth Engine data
layer. For Santa Catarina it computes, on a 0.025° grid (144 × 224 cells, 13,970 of them inside the state) and every
3 h from 2003 to 2025, the emission of 19 MEGAN classes: isoprene, monoterpenes, sesquiterpenes, oxygenated and light
VOCs, soil NO and CO. The input and run procedures (FIX-1 to FIX-8) are in `src/megan_sc/fixes.py`, whose SHA-256 is
stamped in every output file, and the engine, as used in production, agrees stage by stage with the official Fortran
executables to float32 precision. When the analyses read the run, they multiply isoprene by the canopy loss and
production factor ρ = 0.96 (Guenther et al., 2006), which the MEGANv3.21 code does not apply.

## What this repository lets you do

| Level | Goal | Needs | Time | Start at |
|---|---|---|---|---|
| **A. Reproduce** | recompute every number, table and figure of the article and of the report from the shipped model output, and compare them with the published copies | Python 3.10+, about 1.2 GB of disk, no accounts | about 5 min of computing | [Quick start](#quick-start-reproduce-every-number-in-a-few-minutes), [tutorial chapter 5](docs/tutorial/05_level_A_reproduce.md) |
| **B. Replicate** | run the model again from the raw inputs (ERA5-Land, MODIS, MapBiomas), the checks, the Fortran benchmark and the TROPOMI evaluation | Google account with Earth Engine, Google Colab | hours of Earth Engine requests per model year; 41 EECU-hours for the TROPOMI notebook | [tutorial chapter 6](docs/tutorial/06_level_B_rerun_the_model.md) |
| **C. Extend** | build the same inventory for another state, province or country | as B, plus a boundary file | as B | [tutorial chapter 11](docs/tutorial/11_your_region.md), [`run_your_region/`](run_your_region/START_HERE.md) |

## Quick start: reproduce every number in a few minutes

```bash
git clone https://github.com/jessicajcss/MEGAN_SC megan-sc && cd megan-sc
pip install -r requirements-reproduce.txt          # or: conda env create -f environment.yml && conda activate megan-sc
python tools/data_release.py fetch https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0   # the data archives (DATA_AVAILABILITY.md)
python reproduce.py check-data                     # SHA-256 of every data file
python reproduce.py all                            # 16 steps; the last one compares every file with results/v10
```

The address above is that of the GitHub release (tag v1.0.0) that will hold the data archives once the repository is published. The archived copy on Zenodo, once deposited,
or a local folder holding the zip files works the same way (DATA_AVAILABILITY.md). Everything is written under `reproduced/`; nothing in
`data/` or `results/` changes. The last step prints the comparison with the published files. A run of 15 November 2026 on
a Linux machine with newer library versions than `requirements-lock.txt` (Python 3.13.18, numpy 2.5.3, pandas 3.0.5, matplotlib 3.11.2, geopandas 1.2.0), in one pass into an empty folder (5 min 22 s), ended with:

```
CSV files: 58 compared, 26 identical, 32 equal within tolerance, 0 different, 0 without a published copy
PNG files: 59 compared, 17 pixel-identical
```

Four of the identical CSV files (`report/report_numbers_v10.csv`, `final_review/trend_robustness.csv`,
`final_review/month_boundary_effect.csv`, `ijest_review/ecotype_tree_ef.csv`) and the 17 pixel-identical PNG files were
written to `results/v10` on that machine, so the comparison does not test them independently. The other PNG files were
drawn in Google Colab and differ in rendering (fonts and library versions), while the data behind them are equal. Of the
32 CSV files equal within tolerance, 28 differ from their published copies only in the last digits of double-precision
numbers and 4 (TROPOMI tests of notebook 14) by at most 9.5 × 10⁻⁷ in relative terms.
`docs/tutorial/reproduction_log.txt` holds the full comparison, file by file, and the result of the MODIS control
(REPRODUCE.md, section B). Running notebook 13 a second time into the same folder reads its single-precision grid cache
instead of the NetCDF files, and more files then come out as equal within tolerance.
PNG files depend on fonts and library versions, so on another machine some figures may differ by a few pixels while
every number stays the same ([chapter 5](docs/tutorial/05_level_A_reproduce.md) explains how to read the comparison).

## Run the notebooks in Google Colab

The fourteen notebooks of `notebooks/` run in Google Colab with Google Drive (levels B and C). Run them in this order; [notebooks/README.md](notebooks/README.md) gives, for each one, what it needs, what it writes, its settings cell and the time it takes.

| Order | Notebook | Open in Colab |
|---|---|---|
| 1 | `12_fortran_benchmark.ipynb`: comparison with the official Fortran executables | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/12_fortran_benchmark.ipynb) |
| 2 | `10_production_run.ipynb`: the production run (Earth Engine) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/10_production_run.ipynb) |
| 3 | `11_verify_run.ipynb`: checks C0-C5 of every year | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/11_verify_run.ipynb) |
| 4 | `13_manuscript_numbers_and_figures.ipynb`: numbers, tables and figures | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/13_manuscript_numbers_and_figures.ipynb) |
| 5 | `14_tropomi_hcho_evaluation.ipynb`: evaluation with TROPOMI formaldehyde (Earth Engine) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/14_tropomi_hcho_evaluation.ipynb) |
| 6 | `15_sensitivity_2019.ipynb`: sensitivity tests | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/15_sensitivity_2019.ipynb) |
| 7 | `16_lai_independent_check.ipynb`: VIIRS check of the 2024-2025 decrease of the MODIS LAI (Earth Engine) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/16_lai_independent_check.ipynb) |
| 8 | `17_rerun_2024_2025_lai_viirs.ipynb`: 2024 and 2025 with the VIIRS-corrected LAI (Earth Engine) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/17_rerun_2024_2025_lai_viirs.ipynb) |
| 9 | `18_brain_comparison_2019.ipynb`: comparison with the BRAIN MEGAN emissions of 2019 (d02, 4 km) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/18_brain_comparison_2019.ipynb) |
| 10 | `19_install_and_reproduce.ipynb`: checks and installation of 2024-2025 with the corrected LAI, NOAA CO₂ file, reproduction with both LAI inputs | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/19_install_and_reproduce.ipynb) |
| 11 | `21_era5land_vs_inmet.ipynb`: ERA5-Land temperature and radiation against INMET stations, and the lapse rate read by notebook 20 (Earth Engine) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/21_era5land_vs_inmet.ipynb) |
| 12 | `20_sensitivity_2019_inputs.ipynb`: 2019 reruns with one setting changed at a time (Earth Engine) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/20_sensitivity_2019_inputs.ipynb) |
| 13 | `22_lai_check_noaa20_viirs.ipynb`: LAI of the VIIRS on NOAA-20 and NOAA-21 against Suomi NPP (Earthdata login) | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/22_lai_check_noaa20_viirs.ipynb) |
| 14 | `23_brain_differences_2019.ipynb`: tests of the explanations of the differences from BRAIN | [open](https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/23_brain_differences_2019.ipynb) |

The links open the notebooks of the `main` branch. Set your Earth Engine project (`EE_PROJECT`) in the cell "3 · Earth Engine" of notebooks 10, 17 and 20, in the cell "2 · Earth Engine" of notebook 21 and in the settings cell of notebooks 14 and 16; notebook 22 needs an Earthdata login. Notebook 19 installed 2024-2025 with the corrected LAI in the run and was run once by the authors; `reproduce.py` repeats its comparisons ([notebooks/README.md](notebooks/README.md)).

## Main results and where each one comes from

| Result (2003-2025, Santa Catarina) | Value | File in `results/v10/` |
|---|---|---|
| Mean BVOC emission | 657 Gg yr⁻¹; isoprene 33.2 %, monoterpenes 25.8 %, oxygenated VOC 19.3 % | `notebook13_meso2024/manuscript_numbers.csv` (`bvoc_mean_kt`, `mean_*`) |
| Interannual coefficient of variation of BVOC | 7.9 % (peaks 2019 and 2024, 739 Gg; minimum 2008, 552 Gg) | same file (`bvoc_cv_pct`, `bvoc_max`, `bvoc_min`) |
| Trend of isoprene and of soil NO (Theil-Sen, Mann-Kendall) | +0.88 % yr⁻¹ (p = 0.004) and −0.37 % yr⁻¹ (p = 0.002); with total BVOC, monoterpenes, sesquiterpenes and oxygenated VOC, six trends survive Benjamini-Hochberg control (q ≤ 0.023), also over 2003-2023 (q ≤ 0.048) | `notebook13_meso2024/trends.csv`; `final_review/trends_fdr_hr.csv`, `trend_robustness.csv` |
| SOA precursors (monoterpenes + sesquiterpenes) | 183 Gg yr⁻¹; intensity against native-forest share of the six mesoregions: Spearman ρ = 0.94 (exact p = 0.017), Pearson r = 0.97 (95 % CI 0.73 to 1.00) | `notebook13_meso2024/table3_mesoregions_soa.csv`, `manuscript_numbers.csv` (`meso_forest_r`), `final_review/final_review_numbers.csv` (`meso_spearman_exact`) |
| 2019 against the previous state inventory | isoprene 7.6, CO 3.7 and soil NO 18.9 times larger | `manuscript_numbers.csv` (`cmp2019_*`) |
| 2019 against BRAIN (true cell areas) | isoprene 0.79, monoterpenes 0.85, summed VOC 0.87, soil NO 1.13 times BRAIN | `brain_comparison/comparison_2019_classes_t.csv`; `brain_differences/` |
| Parts of the trends (descriptive, unstable) | isoprene: linear trend of the logarithm +0.92 % yr⁻¹ = temperature +0.35, LAI +0.16, time term +0.40 (tree area +9.5 % from 2003 to 2024; isoprene emission potential +0.38 % yr⁻¹); soil NO time term −0.45 % yr⁻¹ (emission potential −0.45 % yr⁻¹) | `final_review/final_review_numbers.csv` (`ols_trend_decomposition_*`, `growthform_area_*`, `ef_potential_*`) |
| Isoprene trend with CO₂ inhibition (NOAA global CO₂, 375 to 426 ppm) | +0.65 % yr⁻¹ (p = 0.018; BH q = 0.03); period mean +0.3 % | same file (`isop_trend_with_co2_inhibition`, `fdr_with_co2_inhibition`, `isop_mean_with_co2_inhibition`) |
| TROPOMI formaldehyde, 2019-2025 | deseasonalised monthly anomalies r = 0.45 (N_eff = 72, p < 0.001); spatial pattern and seasonal cycle not reproduced; observed/expected column 45 % to 125 % | `notebook14/`, `diagnostics/nb14/`, `final_review/` |

All values use the VIIRS-corrected LAI of 2024 and 2025 (the values with the uncorrected MODIS LAI are in `results/v10_modis_lai_2024_2025/`). Isoprene includes ρ = 0.96 (except in the BRAIN comparison, where neither side applies it); municipal and mesoregion totals split each cell among the polygons in proportion to its
area inside each one.

[Chapter 8 of the tutorial](docs/tutorial/08_number_map.md) maps every number of both documents to its file, column
or key, and the script that writes it.

## How the work was checked

* **Against the official code.** `notebooks/12_fortran_benchmark.ipynb` runs the MEGANv3.21 Fortran executables and
  the Python engine on the same inputs. With the engine patches of `fixes.py`, as in production, every stage passes:
  for the test case shipped with the code, the emission activity differs by at most 3.8 × 10⁻⁶ (relative) when both
  use the Fortran canopy ([VALIDATION.md](VALIDATION.md)).
* **The production run itself.** Notebook 11 applies six automatic checks to every year (138 results: 137 PASS and
  1 FAIL, explained in `results/v10/diagnostics/c3/`).
* **Against observations.** Notebook 14 compares modelled isoprene with TROPOMI formaldehyde columns under tests fixed
  before the comparison ([chapter 10](docs/tutorial/10_tropomi.md)).
* **Sensitivity.** Time step (3-hourly against hourly: −0.06 % on BVOC over 28 days of 2019), soil-moisture factor
  (about −0.05 % on BVOC), five settings rerun for 2019 (notebook 20: −2.3 % to +2.6 % on BVOC) and planted-forest
  emission factors (−19 % to +14 % on isoprene).
* **Inputs.** The 2024-2025 decrease of the MODIS LAI was checked against the VIIRS of Suomi NPP, NOAA-20 and NOAA-21
  (notebooks 16 and 22) and corrected (notebooks 17 and 19); ERA5-Land was compared with INMET stations (notebook 21).
* **Against another inventory.** The 2019 emissions were compared with BRAIN (notebooks 18 and 23).
* **Reproduction.** `reproduce.py` recomputes everything from the shipped output and compares it with the published
  files, as above.

## Repository map

```
README.md, README.pt-BR.md this page, in English and Portuguese
REPRODUCE.md               the short version of tutorial chapters 5 and 6
DATA_AVAILABILITY.md       every input and output, its source, licence and checksum
VALIDATION.md              comparison with the Fortran code and with observations; checks of the run
CHANGELOG.md, CITATION.cff, LICENSE (MIT, code), LICENSE-docs-data.md (CC-BY-4.0, documents and results)
reproduce.py               level A in one command (python reproduce.py list)
TUTORIAL_MEGAN_SC_en.pdf   illustrated tutorial (also TUTORIAL_MEGAN_SC_pt.pdf); HTML sources in docs/tutorial/
docs/tutorial/             the technical chapters: 13 chapters, from the equations to the number map and provenance
docs/                      ALGORITHM_SPEC (equations transcribed from the Fortran), TIME_CONVENTIONS, INDEX
src/megan_sc/              the engine (canopy, activity factors, soil NO, speciation) and fixes.py (procedures FIX-1 to FIX-8)
tests/                     127 unit tests (pytest -q)
notebooks/                 10-23: production run, checks, Fortran benchmark, numbers and figures, TROPOMI, sensitivity,
                           VIIRS check of the LAI, rerun of 2024-2025, comparison with BRAIN 2019, installation and
                           reproduction, 2019 sensitivity runs, ERA5-Land against INMET, NOAA-20/21 LAI, BRAIN differences
analysis/v10/              scripts run on the output of notebook 10 (diagnostics, report numbers, maps, study-area map,
                           verification flowchart)
results/v10/               the published results: tables, numbers, figures, logs (corrected LAI of 2024-2025)
results/v10_modis_lai_2024_2025/  the same results with the uncorrected MODIS LAI of 2024-2025
data/                      small inputs in git (boundaries, MEGAN tables, data/external: MapBiomas 2019 window, IBGE
                           outlines, NOAA global CO2); the model output and TROPOMI fields arrive with tools/data_release.py
tools/                     notebook runner, result comparison, data release, test fixtures
benchmark/                 scripts and results of the Fortran comparison
run_your_region/           setup wizard for another region
scripts/                   generator of the region wizard notebook
```

## Data and licences

Inputs come from public products read through Google Earth Engine: ERA5-Land (Muñoz-Sabater et al., 2021), MODIS
MCD15A3H leaf area index (Myneni; Knyazikhin; Park, 2021), MapBiomas Collection 10 land cover (Souza Junior et al.,
2020) and TROPOMI formaldehyde (De Smedt et al., 2018); the global annual mean CO₂ comes from NOAA GML (Lan; Tans;
Thoning, 2026) and the outlines from IBGE. The emission-factor and speciation tables and the ecotype map come from the
MEGAN 3.2 input files (Guenther, 2024a, Zenodo 10.5281/zenodo.10939297, CC BY 4.0): the tables, unchanged, are in
`data/megan_inputs/` and the ecotype map, clipped to Santa Catarina, in `data/sc_ecotype.nc`. The Fortran code of
MEGANv3.21 (Guenther, 2024b, Zenodo 10.5281/zenodo.10526206) is not redistributed here. Code is under the MIT licence;
documents, tables and figures under CC-BY-4.0; third-party data keep their own licences
([DATA_AVAILABILITY.md](DATA_AVAILABILITY.md), [LICENSE-docs-data.md](LICENSE-docs-data.md)).

## Citation

Cite the article and this repository; `CITATION.cff` holds the metadata (GitHub shows it under "Cite this
repository"). The DOI of the archived version will be added there after deposit. Cite also the MEGANv3.21 code
(Guenther, 2024b) and the MEGAN 3.2 input files (Guenther, 2024a); `CITATION.cff` lists both.
