# 5. Level A: reproduce every result

Level A recomputes, from the model output shipped in `data/`, every table, number and figure of `results/v10/`, the
folder from which both documents take their values, and then compares the recomputed files with the published ones.
It needs no Earth Engine and no internet once the data are in place (chapter 3).

## 5.1 One command

```bash
python reproduce.py check-data
python reproduce.py all
```

In a run on 15 November 2026 (a Linux machine with newer library versions than `requirements-lock.txt` (Python 3.13.18, numpy 2.5.3, pandas 3.0.5, matplotlib 3.11.2, geopandas 1.2.0)), the 16 steps took 5 min 22 s in one pass into an
empty folder, and the comparison ended with

```
CSV files: 58 compared, 26 identical, 32 equal within tolerance, 0 different, 0 without a published copy
PNG files: 59 compared, 17 pixel-identical
```

`docs/tutorial/reproduction_log.txt` keeps that run: timings, and the status of every file (section 5.3). Outputs go to
`reproduced/` (about 160 MB with all 16 steps), one log per step to `reproduced/logs/`. A step that fails prints the end of its log;
the other steps still run.

Any subset can be run alone, in the order given, for example `python reproduce.py nb13 report maps compare`. Steps
that read the outputs of notebook 13 (`report`, `maps`, `final`) use the reproduced copy when `nb13` has run, and the
published copy otherwise.

## 5.2 The steps

Each paragraph gives what the step does, what it writes under `reproduced/`, and what to look for in its log.

**`integrity`** (14 s). Reads the fix log of each of the 276 months and checks that every month accumulated eight
steps per day, then that each annual file equals the sum of its twelve monthly files. Log:

```
(a) months with steps != ndays*8: none (276 months checked)
(b) max relative |annual - sum of monthly| over years and classes: 6.2e-10
```

**`nb11`** (30 s). Notebook 11 applies the six checks C0 to C5 to every year. `notebook11/verify_v10_report.csv` has
138 rows: 137 PASS and 1 FAIL, C3 in 2015 for methanol (4 of 71 comparable filled cells above the threshold, where the
bound was 2). The log ends with `FAILS: 1`. `results/v10/diagnostics/c3/` explains why the failure was accepted.

| Check | What it tests |
|---|---|
| C0 provenance | each year carries `provenance_v10.json` with the three engine patches and the same `fixes.py` hash |
| C1 forested cells without isoprene | no state cell with tree fraction > 0.1, BVOC > 0 and isoprene = 0 (the signature of a cell computed without valid meteorology) |
| C2 outliers | no state cell above 5 × the 99th percentile of the cells with valid meteorology |
| C3 filled cells | the emission activity E/(EF · area) of the cells filled by FIX-1 against the median of valid neighbours in a 7 × 7 window, for LAI-matched cells; at most the binomial bound may exceed |
| C4 completeness | 12 months, 12 fix logs, no checkpoint left |
| C5 PAR | monthly state-mean PAR between 20 and 200 W m⁻² (a test of the de-accumulation of the ERA5-Land radiation, FIX-8) |

**`nb13`** (29 s). Notebook 13 computes Tables 2 to 4, the trends, the driver correlations, the 76 entries of
`manuscript_numbers.csv`, the files `fig1.png` to `fig8.png` and `fig_detrend.png` in Portuguese and English (the article figures they correspond to are listed in chapter 9), and `annual_grids_cache.npz`, the annual grids of
all years in one file, which the later steps read. Every value marked in the article comes from this notebook or from
the scripts below (chapter 8).

**`nb14`** (7 s). Notebook 14 runs the four TROPOMI tests from the shipped monthly fields. The log prints, among
others, the anomaly test (`HCHO vs model isoprene (monthly anomalies) 0.452 84.0 ... 71.704 0.000`) and the
climatology (`r_climatology -0.262, p 0.411`). Chapter 10 explains each test.

**`nb15`** (12 s). Notebook 15 computes test S3 from the fix logs: the monthly state mean of the soil-moisture factor
γSM and its approximate effect on BVOC, year by year (0.02 % to 0.10 %).

**`s1`** (2 s). Test S1 compares the two 28-day runs of 2019 (3-hourly and hourly). State totals change by −0.06 %
(BVOC) and −0.18 % (isoprene); emissions at the shared hours are identical.

**`c3`** (55 s). Check C3 in detail: every filled cell of every year with its neighbours, the LAI-adjusted version,
and the effect on the totals. Log line to look for:
`2015 MEOH 4 cells; emission t [2.3 1.9 1.7 0.6] share 0.01476% change if neighbour-median activity +0.00273%`.

**`nb14diag`** (8 s). The diagnostics defined after the TROPOMI tests (D1 to D7): split halves, surface pressure,
magnitude, background, screening.

**`report`** (9 s). Every number of the inventory report (`report/report_numbers_v10.csv`, 104 rows) and report Figure 12 (file `pt_report/fig10.png`).

**`maps`** (40 s). Every map of both documents (article Figures 7 to 9; report Figures 7, 9 to 11, 14 and A.2 to A.9) and the TROPOMI figure (article Figures 10 and 11 in one file, report Figure A.10), with `figures/maps_class_limits.csv`.

**`meso2024`** (2 s). Rebuilds `sc_mesoregions_2024.gpkg` and the municipality table from the IBGE layers into
`reproduced/boundaries/` (the copies in `data/boundaries` stay untouched); the forest share of each mesoregion
(attribute `forest_pct`, used by article Table 2 and Figure S5 and by report Table 3 and Figure 8) is taken from `data/boundaries/landcover2019_mesoregion.csv`.

**`final`** (9 s). Additional diagnostics: Benjamini-Hochberg control, Hamed-Rao check, start-year and 2025 sensitivity
of the trends, LAI of 2024 and 2025, partial correlations, residual year term, background control and magnitude of the
TROPOMI comparison, mesoregion and municipal statistics, pooled C3, filled coastal cells, interannual statistics of the
drivers, the parts of the trends (regression on temperature, LAI and year), the growth-form areas and emission potential
of every year (`landcover_ef_potential_SC.csv`), the LAI gaps and largest LAI, and the effect of CO₂ inhibition on the
isoprene trend with the NOAA global annual means of `data/external/co2_annmean_gl.txt`
(`final_review/final_review_numbers.csv`, `final_review.log`), then the robustness of the trends (`analysis/v10/review/trend_robustness.py`: 2003-2023, trend-free pre-whitening, Benjamini-Yekutieli and the eight groups without total BVOC; `final_review/trend_robustness.csv`, Online Resource 1, Table S9).

**`ijest`** (about 30 s). The figures of the English manuscript in the journal format (`figures/ijest/Fig1` to `Fig11`, `FigS4`, `FigS5`; TIFF copies in `reproduced/figures/ijest/`), its post hoc diagnostics (`ijest_review/`) and the ecotype emission factors of the state. Its log ends with a check line (`check: r_clim -0.262 | anomalies r 0.452 N_eff 71.7 ...`) and the size of each figure.

**`laicorr`**. The VIIRS check of the 2024-2025 MODIS LAI from the cache of notebook 16: hindcast of the candidate corrections over 2013-2023 and the corrected LAI of 2024-2025 (`lai_viirs_check/`). It reads the MODIS-LAI record of every year through the view of `tools/lai_variant_view.py`.

**`extra`** (21 s). The study-area map (article Figure 1: MapBiomas Collection 10 in 2019 from `data/external/`, the
2019 growth forms of the run and their state shares, `figures/study_area_shares.csv`) and the verification flowchart
(article Figure 2, with numbers read from the files of `results/v10`), in Portuguese and English.

**`compare`** (2 s). The comparison of section 5.3.

## 5.3 Reading the comparison

`tools/compare_results.py` walks `reproduced/` and looks for the file with the same relative path in `results/v10/`
(and in `data/boundaries` for the `boundaries/` folder). Each file gets one status:

| Status | Meaning |
|---|---|
| IDENTICAL | same bytes |
| EQUAL | CSV: same columns and rows, every number equal within a relative tolerance of 10⁻⁶ or an absolute tolerance of 10⁻⁹, text equal after removing the `np.float64(...)` wrappers that NumPy 2 prints; GPKG: same attributes and geometries |
| DIFFERENT | a CSV that fails the test above; the line shows the first differing column and values |
| PIXEL-IDENTICAL | PNG with the same pixels (the file bytes can differ by metadata) |
| PNG DIFF | PNG with differing pixels; the number is the share of pixels that differ |
| PRESENT | a log, which carries timings and paths and is not compared |
| NEW | no published copy: the grid cache of notebook 13 and the TIFF copies of the IJEST figures |

The exit code is 1 if any CSV is DIFFERENT. A PNG DIFF alone is not a failure: figures depend on the fonts and on the
versions of matplotlib and of its libraries, so a figure drawn on another machine can differ by a few pixels while its
numbers do not. `pip install -r requirements-lock.txt` installs a tested set of versions.

In the run of 15 November 2026, 26 of the 58 CSV files are identical and 32 equal within tolerance, none different.
Four of the identical CSV files (`report/report_numbers_v10.csv`, `final_review/trend_robustness.csv`,
`final_review/month_boundary_effect.csv`, `ijest_review/ecotype_tree_ef.csv`) and the 17 pixel-identical PNG files were
written to `results/v10` on that machine, so the comparison does not test them independently. Of the 32 files equal
within tolerance, 28 differ only in the last digits of double-precision numbers and 4 (notebook 14) by at most
9.5 × 10⁻⁷; that machine ran newer library versions than `requirements-lock.txt`. A second run
of notebook 13 into the same folder reads the annual grids from the cache it wrote, stored in single precision, instead
of the NetCDF files; the last digits of some numbers then change, and the comparison reports those files as EQUAL.

## 5.4 Checking a number by hand

The comparison shows that the code reproduces the published files. To see how a number is made, recompute it from a
table in a few lines. The isoprene trend of the abstract (+0.88 % per year, p = 0.004):

```python
import numpy as np, pandas as pd, pymannkendall as pmk
from scipy import stats

t = pd.read_csv('results/v10/notebook13_meso2024/table2_annual_totals_SC_groups.csv')   # state totals, t/yr
y = t['ISOP'].values; x = np.arange(len(y))
slope, intercept, lo, hi = stats.theilslopes(y, x, 0.95)
print(100 * slope / y.mean(), 100 * lo / y.mean(), 100 * hi / y.mean())    # 0.88  0.33  1.50  (% of the mean per year)
print(pmk.original_test(y).p, pmk.hamed_rao_modification_test(y).p)       # 0.0037  0.0037
```

The control of the false discovery rate over the nine groups (q = 0.017 for isoprene and soil NO, 0.023 or less for four other groups):

```python
from statsmodels.stats.multitest import multipletests
tr = pd.read_csv('results/v10/final_review/trends_fdr_hr.csv')
print(multipletests(tr['MK_p'], method='fdr_bh')[1].round(4))
# [0.0186 0.0165 0.0202 0.0186 0.0227 0.0813 0.0813 0.0165 0.0813]   (BVOC, ISOP, MT, SQT, OVOC, LVOC, MBO, NO, CO)
```

Chapter 7 explains these estimators, and chapter 8 lists where every other number comes from.

## 5.5 Reproducing from a different run

All steps read the run through `MEGAN_SC_RUN`. To recompute the results from a rerun of notebook 10 (chapter 6), or
from a copy of the run kept elsewhere:

```bash
MEGAN_SC_RUN=/path/to/megan_sc_outputs_v10 MEGAN_SC_OUT=reproduced_rerun python reproduce.py all
```

The comparison then measures how far the rerun is from the published run. Chapter 6 lists the differences to expect.

A rerun of notebook 10 computes 2024 and 2025 with the MODIS LAI. Until those two years are computed with the
corrected LAI (chapter 6, section 6.7), compare it with `results/v10_modis_lai_2024_2025`, using the steps that this
folder holds:

```bash
MEGAN_SC_RUN=/path/to/megan_sc_outputs_v10 MEGAN_SC_RESULTS=results/v10_modis_lai_2024_2025 MEGAN_SC_OUT=reproduced_rerun \
    python reproduce.py integrity nb11 nb13 nb14 nb15 s1 c3 nb14diag report maps meso2024 final laicorr compare
```

`REPRODUCE.md` (section B) gives the same command for the shipped MODIS-LAI record (the MODIS control) and its result.
