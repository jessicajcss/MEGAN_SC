# 6. Level B: run the model again

Level B repeats the first half of the pipeline (chapter 4): the production run from the raw inputs, the Fortran
benchmark, the TROPOMI extraction and the sensitivity runs. It needs Google Colab, Google Drive and an Earth Engine
project (chapter 3, section 3.6). The order below is the order used in October 2026 for notebooks 10 to 15 and between
31 October and 8 November 2026 for notebooks 16 to 23.

## 6.1 Before starting

* Copy the repository to `MyDrive/megan-sc`.
* Decide where the output goes. Notebook 10 writes to `MyDrive/megan_sc_outputs_v10`; the analysis notebooks write
  into folders inside it.
* Expect the cost. A production year takes hours of Earth Engine requests (eight 3-hourly fields and one daily mean
  per day, plus the monthly LAI). The notebook checkpoints after every step, so Colab disconnections cost little:
  run all cells again and it resumes. The TROPOMI notebook used 41.2 EECU-hours and about 6 h for 84 months.

## 6.2 Notebook 12: the Fortran benchmark (independent of the rest)

1. Download `MEGANv3.21.tar.gz` from Zenodo (10.5281/zenodo.10526206) to `MyDrive`.
2. Open `notebooks/12_fortran_benchmark.ipynb` in Colab and set, in its first code cell, `MGN_TGZ` (the tarball),
   `MEGAN_DIR` (the repository) and `RAW` (the run folder whose 2019 inputs build the Santa Catarina test case:
   `<RAW>/2019/inputs/lai_2019_MM.tif`, `met_2019_MM.tif`, `growthform_2019.npz`).
3. Run all cells (15 to 25 min). The results go to `MyDrive/megan_sc_benchmark_results/`.

For the published results, `RAW` was a run folder whose 2019 inputs (2019 monthly LAI with Earth Engine default
resampling, and growth forms) build the Santa Catarina test case. The repository ships the production run
(`data/run_v10`, same layout); with it as `RAW`, the Santa Catarina test days start from the 2019 inputs of the
production run, so their numbers differ from `benchmark/results_colab/` while the test itself, the same inputs given to
the Fortran and to the Python engine, is unchanged. The test case shipped with the Fortran (south-east United States)
does not depend on `RAW`, and its rows should match `benchmark/results_colab/benchmark_summary.csv`.
`python benchmark/make_table.py` prints the table of `VALIDATION.md` from that summary.

What to look for: S0 (our Fortran runs against the reference output in the tarball) bit-identical, and every row of
the patched engine marked PASS.

## 6.3 Notebook 10: the production run

Open `notebooks/10_production_run.ipynb` in Colab. In the settings:

| Cell | Setting | Value for the published run |
|---|---|---|
| 1 | `MEGAN_DIR`, `OUT_BASE`, `OLD_ROOT` | the repository in Drive; `MyDrive`; optional: a folder whose cached MapBiomas growth forms are reused (skipped if absent) |
| 3 | `EE_PROJECT` | your Earth Engine project |
| 4 | `YEARS` | `list(range(2003, 2026))`, or a few years at a time |
| 4 | `SAMPLE_HOURS`, `VARIANT`, `TEST_DAYS`, `LAI_AGG` | 3, `'base'`, `None`, `'mean'` |

Then run all cells. For each year the notebook writes `megan_sc_outputs_v10/<year>/` with the layout of
`data/run_v10/<year>/` (chapter 2), and at the end a short check that must report `signature_cells = 0` and 12 fix
logs per year.

**Keeping the land cover of the published run.** The MapBiomas asset in Earth Engine points to the newest collection,
so a rerun some time after October 2026 may read a newer collection than Collection 10. Notebook 10 reuses the growth
forms of a year when `<OUT>/<year>/inputs/growthform_<vyear>.npz` exists. To keep the published land cover, copy the
shipped files before running:

```bash
# from the repository root, with the Drive folder mounted or synchronised at $DRIVE
for y in $(seq 2003 2025); do
  v=$(( y < 2024 ? y : 2024 ))
  mkdir -p "$DRIVE/megan_sc_outputs_v10/$y/inputs"
  cp "data/run_v10/$y/inputs/growthform_$v.npz" "$DRIVE/megan_sc_outputs_v10/$y/inputs/"
done
```

**What can differ in a rerun.** ERA5-Land and the MODIS collection 6.1 are archived products, and Earth Engine samples
them on the model grid in the same way, so a rerun with the same land cover is expected to come close to the published
run. Two sources of difference remain: a reprocessing of an input collection by its producer, and the fact that the
published run was computed once, so this repository has no measured spread between two production runs. Recompute
everything from the rerun with `MEGAN_SC_RUN` (chapter 5, section 5.5) and look at the comparison.

## 6.4 Notebooks 11, 13 and 15 on the rerun

These read files only. Run them in Colab with their default Drive paths, or locally:

```bash
MEGAN_SC_RUN=/path/to/megan_sc_outputs_v10 MEGAN_SC_OUT=reproduced_rerun python reproduce.py integrity nb11 nb13 nb15
```

Notebook 11 must report zero FAIL rows, or every FAIL must be explained before the numbers are used; the published run
had one (C3, 2015, methanol), explained in `results/v10/diagnostics/c3/`.

## 6.5 Notebook 14: the TROPOMI extraction

Open `notebooks/14_tropomi_hcho_evaluation.ipynb` in Colab, set `EE_PROJECT` and `MODEL_ROOT` (the run) and run all
cells. Cell 3 submits one Earth Engine batch export per month to an asset folder of your project
(`projects/<EE_PROJECT>/assets/megan_sc_hcho_cf0p3_sza70_fire120`), waits for the tasks, reads each month back on the
model grid and caches it in `<OUT_DIR>/cache_cf0.3_sza70_fire120/`. Re-running the cell resumes: months already cached
or already running are not submitted again. The cell prints the months done and the EECU-hours used.

The cache is what `data/tropomi_hcho_v10/` holds. With the cache in place, the notebook runs without Earth Engine,
which is how `reproduce.py nb14` runs it.

Settings that define the published test: cloud fraction < 0.3, solar zenith angle < 70°, fire mask of 20 km on the same
and the previous UTC day (`TAG = cf0.3_sza70_fire120`), Atlantic box 45° W to 41° W and 30° S to 25° S, overpass time
16:30 UTC, `BIAS_CORRECT = False`. Changing any of them changes `TAG`, so a new cache is built next to the old one.

## 6.6 The sensitivity runs of notebook 15

Test S3 needs only the fix logs of the production run. Test S1 needs two short extra runs of notebook 10 for
2019, identical except for the sampling. The published ones used, in the settings of notebook 10:

```python
OUT_BASE  = '/content/drive/MyDrive/s1_test'     # keeps the test runs apart from the production run
YEARS     = [2019]
MONTHS    = [1, 4, 7, 10]
TEST_DAYS = 7                                    # first seven days of each month; files get the suffix _test7d
VARIANT   = 'base'                               # then again with VARIANT = 'hourly' (24 steps per day)
```

`analysis/v10/s1_compare.py` compares the two runs (`MEGAN_SC_S1` = the `s1_test` folder). Test S2
(`VARIANT = 'lai_eedefault'`, the Earth Engine default LAI resampling) was planned and not run.

## 6.7 The LAI of 2024 and 2025 (notebooks 16, 17 and 19)

A rerun of notebook 10 computes 2024 and 2025 with the MODIS LAI, so its results compare with
`results/v10_modis_lai_2024_2025` until the two years are computed with the corrected LAI:

1. Notebook 16 (`EE_PROJECT` in its settings cell) reads the monthly VIIRS LAI of Suomi NPP (`NASA/VIIRS/002/VNP15A2H`)
   and the MODIS LAI on the inventory grid, 2013-2025, in 312 Earth Engine requests, and writes the cache that
   `data/lai_viirs_v10` holds.
2. The step `laicorr` of `reproduce.py` repeats, from that cache, the hindcast that chose the correction
   (`analysis/v10/review/lai_viirs_correction.py`, Online Resource 1, Table S5).
3. Notebook 17 (`EE_PROJECT` in its cell "3 · Earth Engine") computes 2024 and 2025 with the corrected LAI of
   `src/megan_sc/laicorr.py` into `megan_sc_outputs_v10_lai_viirs`; its cell 11 must show that only the LAI differs from
   the production run. Each year takes about as long as a year of notebook 10, and the notebook resumes after a
   disconnection.
4. Cell 12 of notebook 17 (`INSTALL = True`) moves the MODIS years to `megan_sc_outputs_v10_modis_lai_2024_2025` and
   copies the corrected ones into the run; then run notebooks 11, 13 and 14 again. For the published run, notebook 19
   made this installation; it was run once by the authors (`notebooks/README.md`).

## 6.8 Notebooks 18 and 20 to 23

These add the comparison with BRAIN and the tests of Online Resource 1, Sections S8 and S9.

| Notebook | Before running | Time |
|---|---|---|
| 18 | `BRAIN_output_urls.txt` (the BRAIN download links) and `MEGANv3.21.tar.gz` in Drive; no Earth Engine | download of the 364 daily files of 2019, resumable by month |
| 21 | `EE_PROJECT` in its cell "2 · Earth Engine"; the INMET station lists | about 30 min |
| 20 | `EE_PROJECT` in its cell "3 · Earth Engine"; notebook 21 first, because the `lapse` run reads its `lapse_rate.json` (else it uses 6.5 K km⁻¹) | about 40 min for 28 days, about 15 h for the year |
| 22 | an Earthdata account (HDF5 files of VNP15A2H, VJ115A2H and VJ215A2H); no Earth Engine | 1-2 h |
| 23 | the outputs of 18, and those of 14, 20 and 21 when present | about 2 min |

Executed copies of the Earth Engine notebooks show the project ID in their output (the authors' copies of 17, 20 and 21
do); check the outputs before sharing an executed notebook.

## 6.9 After the rerun

1. `python reproduce.py check-data` does not apply to a new run (its checksums differ); check provenance instead:
   every `<year>/provenance_v10.json` and every emission file must carry `megan_sc_fixes_sha256 = a933b29602929a82`,
   and notebook 11 check C0 verifies it.
2. Recompute everything from the new run (chapter 5, section 5.5). Compare it with `results/v10` once 2024 and 2025
   carry the corrected LAI (section 6.7), and with `results/v10_modis_lai_2024_2025` before that.
3. If you publish results from a rerun, ship its run folder and its `results` folder together, as this repository
   does.
