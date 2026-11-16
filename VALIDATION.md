# Validation of the Python engine against the official MEGANv3.21 Fortran

The engine in `src/megan_sc/` is a Python implementation of the MEGANv3.21 code (Guenther, 2024). It is compared,
stage by stage and under identical inputs, with the official Fortran executables built from `MEGANv3.21.tar.gz`
(Zenodo `10.5281/zenodo.10526206`), as used in production, that is with the engine patches of `src/megan_sc/fixes.py`.
Scripts are in `benchmark/`; notebook `notebooks/12_fortran_benchmark.ipynb` runs the whole comparison in Google Colab
(Linux x86-64).

## What is compared

* **S0, harness check.** Our runs of the official binaries on the test case shipped with the code (south-east United
  States, 2 May 2018, 163 × 163 cells, 24 h) against the reference output in the tarball: DAYMET, CANMET and the 19
  MGNERS classes are bit-identical (`benchmark/results/S0_harness_vs_shipped_reference.csv`), so the Fortran settings
  are right.
* **S1-S4.** DAYMET (daily means), MEGCAN (canopy), MEGVEA with the Fortran canopy as input, and end to end, for the
  engine as used in production, with `fixes.apply_engine_patches()` (FIX-6, FIX-7, FIX-8; rows "patched"), and, as a
  control, for the engine modules alone, without the two engine patches FIX-6 and FIX-7 (rows "original port"), which
  shows that the comparison detects a departure from megvea.f and soilnox.F.
* **Soil NO** (Yienger-Levy as coded in MEGSEA; BDSNP off): soil-data path and no-soil path; GAMSM.
* **Speciation and units** (MGN2MECH): 201 compounds in t h⁻¹.
* **Santa Catarina grid** (0.025°, Southern Hemisphere, production canopy setting): 15 January and 15 July 2019, with
  idealised diurnal cycles scaled to the 2019 monthly means. These days test code equivalence, not weather.

## Results (Colab run of notebook 12, `benchmark/results_colab/benchmark_summary.csv`)

The relative differences are cell-hour values unless stated; "domain-sum" compares the totals over the grid and
"day-tonnes" the daily totals in tonnes with the production emission-factor grid.

| Case | Stage | max rel. diff. | p99 rel. diff. | domain-sum rel. diff. | day-tonnes rel. diff. | Verdict |
|---|---|---|---|---|---|---|
| SE_US, original port | S1 daymet | 2.7e-07 | 1.8e-07 | 1.8e-08 |  | PASS |
| SE_US, original port | S2 megcan | 2.8e-02 | 8.1e-05 | 4.2e-06 |  | PASS |
| SE_US, original port | S3 megvea (Fortran canopy in) | 2.5e-01 | 1.5e-01 | 1.8e-02 |  | expected FAIL (port bug, FIX-6) |
| SE_US, original port | S4 end-to-end ER | 2.5e-01 | 1.5e-01 | 1.8e-02 |  | expected FAIL (port bug, FIX-6) |
| SE_US, patched | S1 daymet | 2.7e-07 | 1.8e-07 | 1.8e-08 |  | PASS |
| SE_US, patched | S2 megcan | 2.8e-02 | 8.1e-05 | 4.2e-06 |  | PASS |
| SE_US, patched | S3 megvea (Fortran canopy in) | 3.8e-06 | 2.9e-06 | 7.3e-07 |  | PASS |
| SE_US, patched | S4 end-to-end ER | 1.0e-01 | 1.5e-05 | 7.5e-07 |  | PASS |
| SE_US, patched, GAMSM and GAMBD on | S1 daymet | 2.7e-07 | 1.8e-07 | 1.8e-08 |  | PASS |
| SE_US, patched, GAMSM and GAMBD on | S2 megcan | 2.8e-02 | 8.1e-05 | 4.2e-06 |  | PASS |
| SE_US, patched, GAMSM and GAMBD on | S3 megvea (Fortran canopy in) | 3.8e-06 | 2.9e-06 | 7.7e-07 |  | PASS |
| SE_US, patched, GAMSM and GAMBD on | S4 end-to-end ER | 1.0e-01 | 2.2e-05 | 1.1e-06 |  | PASS |
| SE_US soil, original port | GAMNO YL nosoil-path, corrected=False | 2.6e-01 | 7.8e-02 | 2.3e-03 |  | expected FAIL (port bug, FIX-7) |
| SE_US soil, original port | GAMNO YL nosoil-path, corrected=True | 2.6e-01 | 7.8e-02 | 2.3e-03 |  | expected FAIL (port bug, FIX-7) |
| SE_US soil, original port | GAMNO YL soil-path, corrected=False | 6.7e-07 | 5.6e-07 | 4.3e-07 |  | PASS |
| SE_US soil, original port | GAMNO YL soil-path, corrected=True | 6.7e-07 | 5.6e-07 | 4.3e-07 |  | PASS |
| SE_US soil, original port | GAMSM (WWLT by soil type) | 4.0e-02 | 2.1e-06 | 8.5e-09 |  | PASS |
| SE_US soil, patched | GAMNO YL nosoil-path, corrected=False | 1.8e-06 | 9.0e-07 | 3.8e-07 |  | PASS |
| SE_US soil, patched | GAMNO YL nosoil-path, corrected=True | 1.8e-06 | 9.0e-07 | 3.8e-07 |  | PASS |
| SE_US soil, patched | GAMNO YL soil-path, corrected=False | 6.7e-07 | 5.6e-07 | 4.3e-07 |  | PASS |
| SE_US soil, patched | GAMNO YL soil-path, corrected=True | 6.7e-07 | 5.6e-07 | 4.3e-07 |  | PASS |
| SE_US soil, patched | GAMSM (WWLT by soil type) | 4.0e-02 | 2.1e-06 | 8.5e-09 |  | PASS |
| SC grid, 15 Jan 2019 | S1 daymet | 3.1e-07 | 2.1e-07 | 7.0e-08 |  | PASS |
| SC grid, 15 Jan 2019 | S2 megcan | 1.7e-02 | 1.6e-04 | 4.5e-06 |  | PASS |
| SC grid, 15 Jan 2019 | S3+S4 ER, patched engine (FIX-6) | 1.4e-01 | 2.5e-05 | 3.7e-07 | 4.0e-07 | PASS |
| SC grid, 15 Jan 2019 | S4 ER, original port (no FIX-6) | 2.4e-01 | 1.1e-01 | 3.0e-02 | 3.6e-02 | expected FAIL (port bug, FIX-6) |
| SC grid, 15 Jan 2019 | soil NO GAMNO (FIX-7, SH season as Fortran) | 1.3e-06 | 1.2e-06 | 4.6e-07 |  | PASS |
| SC grid, 15 Jan 2019 | soil NO GAMNO (FIX-7, SH season corrected) | 1.0e+00 | 1.0e+00 | 3.3e-01 |  | deviation by design (SH season) |
| SC grid, 15 Jan 2019 | soil NO GAMNO (original port, SH season as Fortran) | 1.3e-06 | 1.2e-06 | 4.6e-07 |  | PASS |
| SC grid, 15 Jan 2019 | soil NO GAMNO (original port, SH season corrected) | 1.0e+00 | 1.0e+00 | 3.3e-01 |  | deviation by design (SH season) |
| SC grid, 15 Jul 2019 | S1 daymet | 3.2e-07 | 2.1e-07 | 7.4e-08 |  | PASS |
| SC grid, 15 Jul 2019 | S2 megcan | 2.8e-02 | 3.8e-05 | 1.5e-06 |  | PASS |
| SC grid, 15 Jul 2019 | S3+S4 ER, patched engine (FIX-6) | 1.6e-01 | 9.5e-06 | 6.3e-07 | 1.2e-06 | PASS |
| SC grid, 15 Jul 2019 | S4 ER, original port (no FIX-6) | 2.1e-01 | 7.2e-02 | 1.0e-02 | 1.0e-02 | expected FAIL (port bug, FIX-6) |
| SC grid, 15 Jul 2019 | soil NO GAMNO (FIX-7, SH season as Fortran) | 1.5e+03 | 4.2e-01 | 1.5e-02 |  | see note (SH crop cells) |
| SC grid, 15 Jul 2019 | soil NO GAMNO (FIX-7, SH season corrected) | 4.3e+04 | 1.5e+01 | 7.4e-01 |  | deviation by design (SH season) |
| SC grid, 15 Jul 2019 | soil NO GAMNO (original port, SH season as Fortran) | 1.5e+03 | 4.2e-01 | 1.5e-02 |  | see note (SH crop cells) |
| SC grid, 15 Jul 2019 | soil NO GAMNO (original port, SH season corrected) | 4.3e+04 | 1.5e+01 | 7.4e-01 |  | deviation by design (SH season) |
| SE_US speciation | mgn2mech 201 species (t/h) | 3.8e-07 |  | 1.1e-07 |  | PASS |

fixes.py sha256[:16] = a933b29602929a82; rows: 39; PASS: 27; expected FAIL (port bug, FIX-6): 4; deviation by design (SH season): 4; expected FAIL (port bug, FIX-7): 2; see note (SH crop cells): 2

Reading the table:

* The control rows ("original port", engine modules without FIX-6 and FIX-7) fail at S3/S4, where the engine module
  alone applies the canopy-depth factor to the shaded term too (FIX-6 applies it to the sunlit term only, as megvea.f
  does), and in the no-soil NO path, where the engine module alone caps the air temperature in the grass branch only
  (FIX-7 caps it before both branches, as soilnox.F does). Both are expected failures and show that the comparison
  detects such a departure. With the patches, as in production, every stage passes.
* The SC rows marked "deviation by design" compare the Southern-Hemisphere growing season used in production
  (November-May) with the Fortran season, which in that branch compares the day of year with month-day constants;
  "see note (SH crop cells)" marks the July rows where the Fortran fertiliser factor is negative in crop cells and the
  port clips it at zero. Non-crop cells agree.
* The local run (`benchmark/results/`) and the Colab run (`benchmark/results_colab/`) produced identical CSV files.

## Unit tests

`pytest -q` from the repository root runs the engine tests and `tests/test_fixes.py` (12 tests, FIX-1 to
FIX-8): 127 tests in all, all passing (`pytest -q`). Unit tests check the transcription and the procedures of
`fixes.py`; they do not validate the model against observations.

## Checks of the production run (notebook 11)

Six automatic checks per year on the 2003-2025 output (`results/v10/notebook11/verify_v10_report.csv`): C0 provenance
(all procedures active), C1 forested cells without isoprene (the signature of a cell computed without valid
meteorology), C2 no outlier cells, C3 consistency of the cells with filled meteorology with their valid neighbours, C4
completeness, C5 monthly state-mean PAR between 20 and 200 W m⁻². Of 138 results, 137 PASS and 1 FAIL: C3 in 2015 for
methanol, 4 of 71 comparable filled cells above the threshold where the binomial bound was 2. The deviations (-20% to
+25% in emission activity) follow the LAI contrast that the 25% margin of the test still admits; after an LAI
adjustment defined after inspection, none of the 69 tests fails, and replacing the four cells by their neighbours'
median activity changes the 2015 state methanol by +0.003% (`results/v10/diagnostics/c3/`). The failure was accepted
and reported. Every one of the 276 months holds eight steps per day and the annual files equal the sum of the monthly
ones (`results/v10/diagnostics/check_run_integrity.log`).

## Evaluation against TROPOMI formaldehyde (notebook 14)

Tropospheric HCHO columns of TROPOMI (L3 OFFL product in Earth Engine, 2019-2025; cloud fraction < 0.3, solar zenith
angle < 70°, pixels within 20 km of FIRMS fire detections of the same or previous day removed; Atlantic background
subtracted) against modelled isoprene at the overpass time. Pre-specified tests (`results/v10/notebook14/`):

| Test | Result |
|---|---|
| summer spatial pattern, 0.1° and 0.2° blocks | r = 0.12 and 0.18 (95% CI -0.07 to 0.29 and -0.03 to 0.37); r ≤ 0.16 with 10-40 km smoothing |
| monthly climatology | r = -0.26 (p = 0.41); model maximum in January, ΔHCHO maximum in September |
| deseasonalised monthly anomalies | r = 0.45 (N_eff = 72, p < 0.001); ΔHCHO with temperature r = 0.21 |
| column-emission slope, 0.2° blocks | r = 0.18, too weak to interpret |

Diagnostics defined after these tests (`results/v10/diagnostics/nb14/`, `results/v10/final_review/`):

* **Spatial pattern.** The summer ΔHCHO pattern is reproducible between halves of the period (r = 0.76 and 0.86 on
  0.1° and 0.2° blocks) and follows surface pressure (r = 0.82 and 0.88), with larger columns on the coast and in the
  lowlands; modelled isoprene is nearly unrelated to pressure (r ≤ 0.10). Holding pressure fixed, the partial
  correlation between ΔHCHO and modelled isoprene is 0.15 to 0.21 for smoothing of 0 to 40 km.
* **Seasonal cycle and background.** The Atlantic background has a seasonal cycle close to that of modelled isoprene
  (r = 0.67 between the monthly means), so its subtraction removes most of the summer rise; ΔHCHO peaks in August and
  September together with the fraction of observations removed by the screening (r = 0.80).
* **Anomalies.** The monthly anomalies of the background also follow those of modelled isoprene (r = 0.43), a regional
  signal common to both areas. With the background anomaly held fixed, the partial correlation between the ΔHCHO and
  the isoprene anomalies is 0.65; with the 2 m temperature held fixed, 0.48. By season the anomaly correlation is 0.54 (December-February), 0.45 (March-May), −0.21 (June-August) and 0.67 (September-November); without August-October it is 0.41 (N_eff = 58, p = 0.001).
* **Magnitude.** The summer state means, ΔHCHO = 2.2 × 10¹⁵ molecules cm⁻² and modelled isoprene at the overpass
  time = 1.8 × 10¹² atoms C cm⁻² s⁻¹, imply an HCHO lifetime of 0.8 to 1.2 h for yields of 0.43 to 0.29 per carbon, or
  1.2 to 1.8 h with the bias correction of Vigouroux et al. (2020), TROPOMI = 1.10 × 10¹⁵ + 0.64 × FTIR. With the
  lifetimes of 1.5 to 1.8 h that Palmer et al. (2003) simulated for summer columns over North America, the observed
  ΔHCHO is 45 % to 80 % of the value expected from the modelled emissions without the bias correction and 71 % to 125 % with it. The ratio spans 1, so it does not tell whether modelled isoprene is too high or too low.

## Limits of the validation

The model has not been compared with measured fluxes or concentrations in Santa Catarina. The TROPOMI comparison
supports the temporal anomalies of modelled isoprene; it does not constrain the spatial pattern, which in the
observations follows the terrain, and, given the range of lifetimes, yields and bias corrections in the literature, it
leaves the magnitude open.

## Time step (test S1 of notebook 15)

Two test runs of notebook 10 on the first seven days of January, April, July and October 2019, identical except for the
sampling (3-hourly, as in production, and hourly), give identical emissions at the hours they share; over the 28 days the
3-hourly step changes the state totals by −0.06 % for BVOC and −0.18 % for isoprene (monthly isoprene −2.2 % in July to
+0.5 % in January; cell ratios 0.981 to 1.017 between the 1st and 99th percentiles). Files: `results/v10/notebook15/s1/`.

## Soil-moisture factor (test S3 of notebook 15)

From the fix logs of the production run (`results/v10/notebook15/S3_fixlog_*.csv`): the monthly state mean of the
soil-moisture factor γSM stayed between 0.992 and 1.000, so γSM lowers BVOC by about 0.05 % (0.02 % to 0.10 % by
year).

## Reproduction of the published results

`python reproduce.py all` recomputes every table, number and figure in `results/v10/` from the shipped model output and
compares each file with its published copy (`tools/compare_results.py`). The run of 16 November 2026 is kept in
`docs/tutorial/reproduction_log.txt`: 59 CSV files compared, 54 identical and 5 equal within tolerance, none different;
34 of 46 PNG files pixel-identical. Chapter 5 of the tutorial explains how to read it.
