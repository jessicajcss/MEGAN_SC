# How to reproduce and replicate the results

This page is the short version. The tutorial explains every step: [chapter 3](docs/tutorial/03_setup.md) (setup,
Google Drive and Earth Engine), [chapter 5](docs/tutorial/05_reproduce_the_example.md) (reproduce the Santa Catarina
example), [chapter 6](docs/tutorial/06_run_the_model.md) (run the model) and
[chapter 8](docs/tutorial/08_your_region.md) (your own region and period).

## 0. Environment and data

```bash
pip install -r requirements-reproduce.txt        # reproduction of the example only
# or the full environment (Earth Engine, notebooks, maps):
conda env create -f environment.yml && conda activate megan-sc
pytest -q                                        # unit tests of the engine, fixes.py, the region wizard and an offline run

python tools/data_release.py fetch https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0
python reproduce.py check-data                   # every data file against data/DATA_SHA256SUMS.txt
```

`fetch` downloads `data/run_v10`, `tropomi_hcho_v10`, `s1_test`, `lai_viirs_v10` and `run_v10_modis_lai_2024_2025` from
the release v1.0.0; a folder holding the five zips also works.

## A. Reproduce every number, table and figure of the example from the shipped model output

```bash
python reproduce.py list        # the steps and what each one writes
python reproduce.py all         # everything, then the comparison with results/v10
python reproduce.py nb13 maps   # or any subset, in the order given
```

| Step | What it recomputes |
|---|---|
| `integrity` | eight steps per day in every month; annual files equal the sum of the monthly ones |
| `nb11` | notebook 11: checks C0-C5 for every year |
| `nb13` | notebook 13: annual totals by class and group, trends, correlations, drivers, seasonal cycle, mesoregions, municipalities, comparison of 2019 with the previous state inventory, `key_numbers.csv`, the files `fig1-8.png` and `fig_detrend.png`, the grid cache |
| `nb14` | notebook 14: the four TROPOMI tests, from the cached monthly HCHO fields (no Earth Engine) |
| `nb15` | notebook 15: test S3 (soil-moisture factor) from the fix logs |
| `s1` | test S1, 3-hourly against hourly sampling |
| `c3` | check C3 in detail and its effect on the totals |
| `nb14diag` | post hoc diagnostics of the TROPOMI comparison |
| `report` | statistics by class, mesoregion and municipality, standardised regression, temperature sensitivity and its figure (`pt_report/fig10.png`) |
| `maps` | every map (`fig5`, `fig7`, `fig8`, `pt_report/fig9`, `fig12`, `figA2`-`figA9`) and the TROPOMI figure (`fig9`) |
| `meso2024` | MapBiomas 2024 land-cover shares by mesoregion and the mesoregions as unions of the 2024 municipalities |
| `final` | additional diagnostics, including the parts of the trends, the effect of CO₂ inhibition, the robustness of the trends and the effect of UTC month boundaries |
| `extra` | study-area map and verification flowchart |
| `posthoc` | post hoc diagnostics (planted-forest emission factors, LAI of 2024-2025) and the ecotype emission factors of the state |
| `laicorr` | VIIRS check of the 2024-2025 MODIS LAI (cache of notebook 16): hindcast of the candidate corrections and the corrected LAI of 2024-2025 read by notebook 17 |
| `compare` | every reproduced file against its published copy |

Outputs go to `reproduced/` (or `$MEGAN_SC_OUT`), logs to `reproduced/logs/<step>.log`. On 16 November 2026 the 16 steps took 3 min 34 s on a Linux
machine with the versions of `requirements-lock.txt`: CSV files 59 compared, 54 identical, 5 equal within tolerance, 0
different; PNG files 46 compared, 34 pixel-identical (`docs/tutorial/reproduction_log.txt`). Fonts and library versions
can change a few pixels of a figure, or the position of its panels, on another machine; numbers do not change.

## B. Replicate: run the model and the evaluation of the example again (Google Colab and Earth Engine)

Put the repository in Google Drive as `MyDrive/MEGAN_SC` and register an Earth Engine project
([chapter 3](docs/tutorial/03_setup.md)). Each notebook then opens in Colab from Drive, or at
`https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/<notebook>.ipynb`.

| Order | Notebook | Needs | Produces |
|---|---|---|---|
| 1 | `notebooks/12_fortran_benchmark.ipynb` | Linux x86-64 (Colab), `MEGANv3.21.tar.gz` from Zenodo 10.5281/zenodo.10526206 | the table of `VALIDATION.md` |
| 2 | `notebooks/10_production_run.ipynb` with `config/santa_catarina_2003_2025.yaml` | an Earth Engine project | `megan_sc_outputs_v10/<year>/` (same layout as `data/run_v10`); resumable |
| 3 | `notebooks/11_verify_run.ipynb` | output of 10 | `verify_v10_report.csv` |
| 4 | `notebooks/13_santa_catarina_analysis.ipynb` | output of 10, `data/boundaries` | the tables, numbers and figures of `notebook13_meso2024/` |
| 5 | `notebooks/25_results_any_region.ipynb` | output of 10 | the standard results of the run |
| 6 | `notebooks/14_tropomi_hcho_evaluation.ipynb` | output of 10 for 2019-2025, an Earth Engine project (about 41 EECU-hours) | the TROPOMI tests and the monthly fields of `data/tropomi_hcho_v10` |
| 7 | `notebooks/15_sensitivity_2019.ipynb` | output of 10, plus the two test runs of notebook 10 for S1 | tests S1 and S3 |
| 8 | `notebooks/16_lai_independent_check.ipynb` | output of 10 (grid and growth forms), an Earth Engine project | monthly VIIRS and MODIS LAI on the inventory grid, 2013-2025 (cache in `data/lai_viirs_v10`, read by step `laicorr`) |
| 9 | `notebooks/17_rerun_2024_2025_lai_viirs.ipynb` | output of 10, the cache of 16, an Earth Engine project | 2024 and 2025 with the VIIRS-corrected LAI (`megan_sc_outputs_v10_lai_viirs`) |
| 10 | `notebooks/18_brain_comparison_2019.ipynb` | output of 10 for 2019, `data/boundaries`, `MEGANv3.21.tar.gz`, the list of BRAIN download links; no Earth Engine | the comparison with the BRAIN MEGAN emissions of 2019 |
| 11 | `notebooks/19_install_and_reproduce.ipynb` | outputs of 10, 16 and 17, the release archives of the run | installation of 2024-2025 with the corrected LAI, the NOAA CO₂ file, `reproduce.py` with the MODIS and the corrected LAI, release archives (run once to prepare the example: see below) |
| 12 | `notebooks/21_era5land_vs_inmet.ipynb` | INMET station lists, an Earth Engine project | ERA5-Land against the INMET stations; empirical lapse rate (`lapse_rate.json`, read by 20) |
| 13 | `notebooks/20_sensitivity_2019_inputs.ipynb` | output of 10 for 2019, an Earth Engine project, the lapse rate of 21 | one-at-a-time sensitivity runs of 2019 |
| 14 | `notebooks/22_lai_check_noaa20_viirs.ipynb` | output of 10, the cache of 16, an Earthdata account | the LAI check repeated with the VIIRS of NOAA-20 and NOAA-21 |
| 15 | `notebooks/23_brain_differences_2019.ipynb` | the outputs of 18 (and of 14, 20, 21 when present) | tests of the explanations of the differences from BRAIN |

The step `laicorr` reads the MODIS-LAI record of every year through the view that `tools/lai_variant_view.py` builds
from `data/run_v10` and `data/run_v10_modis_lai_2024_2025`. Notebook 19 was run once by the authors: it installed the
2024-2025 run of notebook 17 in the run folder, kept the MODIS-LAI years apart and packed the release archives, starting
from folders that are not distributed. Its MODIS control, which recomputes `results/v10_modis_lai_2024_2025` from the
MODIS-LAI record, is repeated with the steps that this folder holds:

```bash
python tools/lai_variant_view.py
MEGAN_SC_RUN=data/run_v10_modis_view MEGAN_SC_RESULTS=results/v10_modis_lai_2024_2025 MEGAN_SC_OUT=reproduced_modis \
    python reproduce.py integrity nb11 nb13 nb14 nb15 s1 c3 nb14diag report maps meso2024 final laicorr compare
```

The steps `extra` and `posthoc` are left out, as in the run quoted below: `extra` draws the verification figure, which
reads the 2019 sensitivity runs of notebook 20 (`sensitivity_2019/`) from `MEGAN_SC_RESULTS`, and that folder has none. On 16 November 2026 this
command ended with `CSV files: 53 compared, 53 identical, 0 equal within tolerance, 0 different, 0 without a published
copy` and `PNG files: 42 compared, 41 pixel-identical`.

Then point `reproduce.py` at the new run (`MEGAN_SC_RUN=/path/to/megan_sc_outputs_v10 python reproduce.py all`) to
recompute everything else from it. A rerun reads the Earth Engine collections as they are on that day, so its numbers
can differ slightly from the published ones; chapter 6 lists the sources of difference and how to limit them.

## C. Another region and period

1. `run_your_region/region_setup.ipynb` (wizard; `run_your_region/START_HERE.md`) writes `config/<region>_run.yaml`.
2. `notebooks/10_production_run.ipynb`, `11_verify_run.ipynb` and `25_results_any_region.ipynb`, in that order, with
   `RUN_CONFIG = 'config/<region>_run.yaml'` ([chapter 8](docs/tutorial/08_your_region.md)).

Notebook 10 and the wizard run the engine with the procedures of `fixes.py` (`megan_sc.fixes.run_timestep_fixed` after
`fixes.apply_engine_patches()`), never `run.run_timestep` alone: a land-only reanalysis such as ERA5-Land has no data
in some coastal cells, and the engine alone would turn the fill value into T = 340 K and PAR = 0; FIX-1 gives each such
cell the meteorology of the nearest valid land cell.
