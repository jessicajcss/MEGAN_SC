# How to reproduce and replicate the results

This page is the short version. The tutorial explains every step: [chapter 3](docs/tutorial/03_setup.md) (setup),
[chapter 5](docs/tutorial/05_level_A_reproduce.md) (level A) and [chapter 6](docs/tutorial/06_level_B_rerun_the_model.md)
(level B).

## 0. Environment and data

```bash
pip install -r requirements-reproduce.txt        # level A only
# or the full environment (Earth Engine, notebooks, maps):
conda env create -f environment.yml && conda activate megan-sc
pytest -q                                        # 127 unit tests of the engine, fixes.py and the region wizard

python tools/data_release.py fetch <base-url>    # data/run_v10, tropomi_hcho_v10, s1_test, lai_viirs_v10, run_v10_modis_lai_2024_2025
python reproduce.py check-data                   # every data file against data/DATA_SHA256SUMS.txt
```

`<base-url>` will be `https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0` (GitHub release) or
`https://zenodo.org/records/<record>/files` (archived copy) once the repository is published. A folder holding the zips also works.

## A. Reproduce every number, table and figure from the shipped model output

```bash
python reproduce.py list        # the 16 steps and what each one writes
python reproduce.py all         # everything, then the comparison with results/v10
python reproduce.py nb13 maps   # or any subset, in the order given
```

| Step | What it recomputes |
|---|---|
| `integrity` | eight steps per day in every month; annual files equal the sum of the monthly ones |
| `nb11` | notebook 11: checks C0-C5 for every year |
| `nb13` | notebook 13: Tables 2-4, trends, correlations, drivers, seasonal cycle, `manuscript_numbers.csv`, the files `fig1-8.png` and `fig_detrend.png` (article Figures 3-6, S4 and S5; drafts of the maps), the grid cache |
| `nb14` | notebook 14: the four TROPOMI tests, from the cached monthly HCHO fields (no Earth Engine) |
| `nb15` | notebook 15: test S3 (soil-moisture factor) from the fix logs |
| `s1` | test S1, 3-hourly against hourly sampling |
| `c3` | check C3 in detail and its effect on the totals |
| `nb14diag` | post hoc diagnostics of the TROPOMI comparison |
| `report` | the numbers and Figure 12 of the inventory report (file `pt_report/fig10.png`) |
| `maps` | every map of both documents (article Figures 7-9) and the TROPOMI figure (article Figures 10 and 11, report Figure A.10) |
| `meso2024` | the mesoregions as unions of the 2024 municipalities |
| `final` | additional diagnostics, including the parts of the trends, the effect of CO₂ inhibition, the robustness of the trends (Online Resource 1, Table S9) and the effect of UTC month boundaries |
| `extra` | study-area map and verification flowchart (article Figures 1 and 2) |
| `ijest` | English manuscript for IJEST: Figures 1-11 and S4-S5 in the journal format and its post hoc diagnostics (planted-forest emission factors, LAI of 2024-2025, ecotype emission factors) |
| `laicorr` | VIIRS check of the 2024-2025 MODIS LAI (cache of notebook 16): hindcast of the candidate corrections and the corrected LAI of 2024-2025 read by notebook 17 |
| `compare` | every reproduced file against its published copy |

Outputs go to `reproduced/` (or `$MEGAN_SC_OUT`), logs to `reproduced/logs/<step>.log`. Expected result of `compare`
(run of 15 November 2026, 16 steps in one pass, a Linux machine with newer library versions than `requirements-lock.txt` (Python 3.13.18, numpy 2.5.3, pandas 3.0.5, matplotlib 3.11.2, geopandas 1.2.0);
`docs/tutorial/reproduction_log.txt`): 58 CSV files compared, 26 identical and 32 equal within tolerance, none different,
and 17 of 59 PNG files pixel-identical. Four of the identical CSV files and the 17 identical PNG files were written to
`results/v10` on that machine, so they are not tested independently (the header of the log lists them). Fonts and library
versions can change a few pixels of a figure on another machine; numbers do not change.

## B. Replicate: run the model and the evaluation again (Google Colab and Earth Engine)

Each notebook opens in Colab at `https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/<notebook>.ipynb`.

| Order | Notebook | Needs | Produces |
|---|---|---|---|
| 1 | `notebooks/12_fortran_benchmark.ipynb` | Linux x86-64 (Colab), `MEGANv3.21.tar.gz` from Zenodo 10.5281/zenodo.10526206 | the table of `VALIDATION.md` |
| 2 | `notebooks/10_production_run.ipynb` | an Earth Engine project | `megan_sc_outputs_v10/<year>/` (same layout as `data/run_v10`); resumable |
| 3 | `notebooks/11_verify_run.ipynb` | output of 10 | `verify_v10_report.csv` |
| 4 | `notebooks/13_manuscript_numbers_and_figures.ipynb` | output of 10, `data/boundaries` | every number, table and figure of the article |
| 5 | `notebooks/14_tropomi_hcho_evaluation.ipynb` | output of 10 for 2019-2025, an Earth Engine project (about 41 EECU-hours) | the TROPOMI tests and the monthly fields of `data/tropomi_hcho_v10` |
| 6 | `notebooks/15_sensitivity_2019.ipynb` | output of 10, plus the two test runs of notebook 10 for S1 | tests S1 and S3 |
| 7 | `notebooks/16_lai_independent_check.ipynb` | output of 10 (grid and growth forms), an Earth Engine project | monthly VIIRS and MODIS LAI on the inventory grid, 2013-2025 (run on 31 October 2026; cache in `data/lai_viirs_v10`, read by step `laicorr`) |
| 8 | `notebooks/17_rerun_2024_2025_lai_viirs.ipynb` | output of 10, the cache of 16, an Earth Engine project | 2024 and 2025 with the VIIRS-corrected LAI (`megan_sc_outputs_v10_lai_viirs`) |
| 9 | `notebooks/18_brain_comparison_2019.ipynb` | output of 10 for 2019, `data/boundaries`, `MEGANv3.21.tar.gz`, the list of BRAIN download links; no Earth Engine | the comparison with the BRAIN MEGAN emissions of 2019 (Online Resource 1, Section S9) |
| 10 | `notebooks/19_install_and_reproduce.ipynb` | outputs of 10, 16 and 17, the release archives of the run | installation of 2024-2025 with the corrected LAI, the NOAA CO₂ file, `reproduce.py` with the MODIS and the corrected LAI, release archives (run once by the authors: see below) |
| 11 | `notebooks/21_era5land_vs_inmet.ipynb` | INMET station lists, an Earth Engine project | ERA5-Land against the INMET stations; empirical lapse rate (`lapse_rate.json`, read by 20) |
| 12 | `notebooks/20_sensitivity_2019_inputs.ipynb` | output of 10 for 2019, an Earth Engine project, the lapse rate of 21 | one-at-a-time sensitivity runs of 2019 |
| 13 | `notebooks/22_lai_check_noaa20_viirs.ipynb` | output of 10, the cache of 16, an Earthdata account | the LAI check repeated with the VIIRS of NOAA-20 and NOAA-21 |
| 14 | `notebooks/23_brain_differences_2019.ipynb` | the outputs of 18 (and of 14, 20, 21 when present) | tests of the explanations of the differences from BRAIN |

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

The steps `extra` and `ijest` are left out because they draw the verification figure, which reads the 2019 sensitivity
runs of notebook 20 (`sensitivity_2019/`) from `MEGAN_SC_RESULTS`, and that folder has none. On 15 November 2026 this
command ended with `CSV files: 52 compared, 21 identical, 31 equal within tolerance, 0 different, 0 without a published
copy`.

Then point `reproduce.py` at the new run (`MEGAN_SC_RUN=/path/to/megan_sc_outputs_v10 python reproduce.py all`) to
recompute everything else from it. A rerun reads the Earth Engine collections as they are on that day, so its numbers
can differ slightly from the published ones; chapter 6 lists the sources of difference and how to limit them.

## C. Another region

`run_your_region/START_HERE.md` and `region_setup.ipynb` (wizard), or notebook 10 with the edits of
[chapter 11](docs/tutorial/11_your_region.md). Use the engine with the procedures of `fixes.py`
(`megan_sc.fixes.run_timestep_fixed` after `fixes.apply_engine_patches()`), never `run.run_timestep` alone: a
land-only reanalysis such as ERA5-Land has no data in some coastal cells, and the engine alone would turn the fill value
into T = 340 K and PAR = 0; FIX-1 gives each such cell the meteorology of the nearest valid land cell.
