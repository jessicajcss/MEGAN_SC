# 6. Run the model

This chapter runs the first half of the pipeline (chapter 4): notebook 10 computes an inventory from the raw inputs in
Google Earth Engine, notebook 11 checks it, and notebook 25 turns it into tables, maps and trends. The same three
notebooks serve the Santa Catarina example and any other region: a run is defined by a configuration file, and
`config/santa_catarina_2003_2025.yaml` is the configuration of the example. Chapter 8 writes a configuration for another
region. The last sections cover the other notebooks of the example (TROPOMI, sensitivity tests, LAI check, benchmark).

## 6.1 Before starting

* The repository is in Google Drive at `MyDrive/MEGAN_SC` (chapter 3, section 3.2), or on a computer with Jupyter and
  Earth Engine credentials (section 3.7.3).
* You have the id of a Cloud project registered for Earth Engine (section 3.7.1).
* Outputs go to `OUT_BASE/<output.dir of the configuration>`. In Colab `OUT_BASE` is `MyDrive`, so the example writes to
  `MyDrive/megan_sc_outputs_v10`; elsewhere it is the folder that contains the repository.
* Expect the cost. Each simulated day takes about nine Earth Engine requests (section 3.7.5); a year of the example
  takes hours. The notebook saves a checkpoint after every step, so a Colab disconnection costs little: run all cells
  again and it resumes. Start with one month or with `TEST_DAYS = 2` to see the whole chain work in minutes.

## 6.2 The run configuration

A run configuration is a YAML file with seven blocks. The example:

| Block | Key | Example | Meaning |
|---|---|---|---|
| `run` | `name` | `santa_catarina_2003_2025` | name of the run |
| `period` | `start`, `end` | `"2003-01"`, `"2025-12"` | first and last month, both included (UTC calendar months) |
| | `sample_hours` | 3 | ERA5-Land read at 00, 03, ..., 21 UTC: eight model steps per day |
| `domain` | `mode` | `bbox` | `bbox`: the grid covers the box `bbox`; `shapefile`: the grid covers the boundary file `shapefile` |
| | `resolution_deg` | 0.025 | cell size in degrees (144 × 224 cells for the example) |
| `region` | `name`, `abbr` | `Santa Catarina`, `SC` | names used in the outputs (`annual_tonnes_SC`) |
| | `boundary` | `data/boundaries/sc_state_ibge.gpkg` | the area reported: a cell counts when its centre lies inside |
| | `zones`, `zones_id` | municipalities, `CD_MUN` | optional: one table row per polygon and month |
| `datasets` | `vegetation.growth_form_source` | `mapbiomas` | `mapbiomas` (Brazil, 30 m) or `modis` (MCD12Q1, worldwide, 500 m) |
| | `vegetation.collection`, `first_year`, `last_year` | 10, 1985, 2024 | land-cover product and its years; a year outside them uses the nearest one (2025 uses 2024) |
| | `vegetation.tree_canopy_type` | 3 | canopy type of trees in MEGCAN: 3 temperate broadleaf, 2 tropical broadleaf, 1 needleleaf |
| | `soil.constant_wilt` | 0.12 | wilting point (m³ m⁻³) of the soil-moisture factor |
| | `emission_factors` | `data/megan_inputs/...`, `data/sc_ecotype.nc` | the MEGAN 3.2 tables and the ecotype map |
| `earth_engine` | `project` | `your-ee-project` | your Cloud project id (or set `EE_PROJECT` in the notebook) |
| `species` | `mechanism` | `CB6X` | mechanism of the speciated output |
| `output` | `dir`, `tag` | `megan_sc_outputs_v10`, `SC_v10` | output folder and the tag of its summary files |

Paths inside the repository are relative to it; files elsewhere (a boundary in your Drive) can be absolute.

## 6.3 Notebook 10: the production run

Open `notebooks/10_production_run.ipynb` in Colab (double-click it in Drive) or in Jupyter. Its settings cell (cell 3):

| Setting | Default | Meaning |
|---|---|---|
| `RUN_CONFIG` | `config/santa_catarina_2003_2025.yaml` | the run configuration |
| `EE_PROJECT` | `None` | your project id; `None` takes `earth_engine.project` of the configuration |
| `OUT_BASE` | `MyDrive` in Colab, else the folder of the repository's parent | where the output folder is created |
| `YEARS` | `None` | every year of the period, or a few years, for example `[2019]` |
| `MONTHS` | `None` | every month of the period, or some calendar months, for example `[1, 4, 7, 10]` |
| `TEST_DAYS` | `None` | the first N days of each month, for a preview; files get the suffix `_testNd` and are not read by notebooks 11 and 25 |
| `VARIANT` | `'base'` | `'hourly'` (24 steps per day) or `'lai_eedefault'` (Earth Engine default LAI resampling): the sampling tests of the example |
| `LAI_AGG` | `'mean'` | area mean of the good-quality 500 m MODIS pixels in each cell |

Then **Runtime > Run all**. Cell 4 prints the region, the years, the steps per day, the land-cover source and the output
folder, and stops with a message if the project id is missing or a file of the configuration is not found. For each
year the notebook writes `<output>/<year>/` with the layout described in `data/README.md` (section `run_v10/`):
monthly and annual emission files in tonnes per cell, the inputs as the model saw them, a fix log per month, the
diurnal profiles, the monthly statistics, the zone tables and `provenance_v10.json`. At the end it prints, for every
year, `signature_cells` (must be 0: no forested cell of the region without isoprene) and the number of fix logs
(must be 12 for a complete year).

**Rerunning the Santa Catarina example.** Keep `RUN_CONFIG` as it is and set `EE_PROJECT`. The configuration fixes
MapBiomas Collection 10 (`collection: 10`), so a later Collection 11 in Earth Engine does not change the land cover;
when Earth Engine holds several releases of Collection 10, notebook 10 takes the newest one
(`mapbiomas.resolve_collection_id`) and, when it reads the land cover from Earth Engine, writes the id it used into
`provenance_v10.json`. To reuse exactly the growth
forms of the shipped run, copy them into the output folder before running; notebook 10 reads
`<output>/<year>/inputs/growthform_<land-cover year>.npz` when it exists:

```bash
# from the repository root, with the release archives unpacked and Drive mounted or synchronised at $DRIVE
for y in $(seq 2003 2025); do
  v=$(( y < 2024 ? y : 2024 ))
  mkdir -p "$DRIVE/megan_sc_outputs_v10/$y/inputs"
  cp "data/run_v10/$y/inputs/growthform_$v.npz" "$DRIVE/megan_sc_outputs_v10/$y/inputs/"
done
```

**What can differ in a rerun.** ERA5-Land and MODIS collection 6.1 are archived products, and Earth Engine samples them
on the model grid in the same way, so a rerun with the same land cover is expected to come close to the shipped run.
Two sources of difference remain: a reprocessing of an input collection by its producer, and the fact that the shipped
run was computed once, so the repository has no measured spread between two production runs. The years 2024 and 2025
of the shipped run carry a VIIRS-corrected LAI (section 6.7); a rerun of notebook 10 computes them with the MODIS LAI,
so compare it with `results/v10_modis_lai_2024_2025` for those years (chapter 5, section 5.5).

## 6.4 Notebook 11: the checks

Notebook 11 reads the run and its configuration and applies six checks to every year. Set the same `RUN_CONFIG` (and
`OUT_BASE` if you changed it) and run all cells. It writes `verify_v10_report.csv` and the annual totals of the region
(`annual_totals_<tag>_all_classes.csv`, isoprene times ρ = 0.96) into the run folder.

| Check | What it tests | Pass criterion |
|---|---|---|
| C0 provenance | every year was produced by notebook 10 with the same code | `provenance_v10.json` present, engine patches FIX-6, FIX-7 and FIX-8 listed, the same `fixes.py` hash in every year |
| C1 cells without meteorology | a forested cell that emits BVOC but no isoprene, the signature of a cell computed with zero light | no such cell in the region |
| C2 outliers | cells far outside the intensity distribution | no cell above 5 × the 99th percentile of the cells with valid meteorology |
| C3 filled cells | cells whose meteorology was filled (FIX-1) behave like their valid neighbours | the emission activity E/(EF · area) of filled cells against the median of valid neighbours in a 7 × 7 window, for LAI-matched cells; exceedances within a binomial 99 % bound |
| C4 completeness | 12 months, 12 fix logs, no checkpoint left | all present (a partial year fails C4 only) |
| C5 PAR | the de-accumulation of the ERA5-Land radiation (FIX-8) | monthly regional mean PAR between 20 and 200 W m⁻² |

Every FAIL must be understood before the numbers are used. The shipped run has one, C3 in 2015 for methanol (four
filled cells above the threshold where the bound was two); `results/v10/diagnostics/c3/` shows why it was accepted.

## 6.5 Notebook 14: the TROPOMI extraction (Santa Catarina example)

Open `notebooks/14_tropomi_hcho_evaluation.ipynb`, set `EE_PROJECT` and `MODEL_ROOT` (the run) and run all cells.
Cell 3 submits one Earth Engine batch export per month to an asset folder of your project
(`projects/<EE_PROJECT>/assets/megan_sc_hcho_cf0p3_sza70_fire120`), waits for the tasks, reads each month back on the
model grid and caches it in `<OUT_DIR>/cache_cf0.3_sza70_fire120/`. Running the cell again resumes: months already cached
or already running are not submitted again. The cell prints the months done and the EECU-hours used.

The cache is what `data/tropomi_hcho_v10/` holds. With the cache in place, the notebook runs without Earth Engine,
which is how `reproduce.py nb14` runs it.

Settings that define the test of the example: cloud fraction < 0.3, solar zenith angle < 70°, fire mask of 20 km on the
same and the previous UTC day (`TAG = cf0.3_sza70_fire120`), Atlantic background box 45° W to 41° W and 30° S to
25° S, overpass time 16:30 UTC, `BIAS_CORRECT = False`. Changing any of them changes `TAG`, so a new cache is built
next to the old one. For another region, the background box must be moved to a clean area near it (chapter 10).

## 6.6 Notebook 25: the results of any run

Notebook 25 reads the run of a configuration and writes the standard results of an inventory into `<run>/results/`.
Set `RUN_CONFIG` (and `OUT_BASE` or `RUN_ROOT` if needed) and run all cells; it needs no Earth Engine. Chapter 9
explains each file and how to read it. The drivers, the correlations and the trends need at least 10 complete years
(`MIN_YEARS_STATS`); with fewer complete years the notebook writes the totals, the seasonal cycle and the maps, and with
no complete year (a run of a few months) the monthly totals and the seasonal cycle of the months present.

For the Santa Catarina example, notebook 13 computes, in addition, what needs the boundaries of the state (mesoregions,
municipalities) and the previous state inventory; both notebooks use the same group definitions, the same factor ρ and
the same estimators.

## 6.7 The LAI of 2024 and 2025 of the example (notebooks 16, 17 and 19)

The statewide mean MODIS LAI of 2025 was 11 % below its 2003-2023 mean (`results/v10/posthoc_review`), and the orbits
of the MODIS satellites have drifted since their orbit maintenance ended. The shipped run therefore computes 2024 and
2025 with a LAI corrected against VIIRS:

1. Notebook 16 (`EE_PROJECT` in its settings cell) reads the monthly VIIRS LAI of Suomi NPP (`NASA/VIIRS/002/VNP15A2H`)
   and the MODIS LAI on the inventory grid, 2013-2025, in 312 Earth Engine requests, and writes the cache that
   `data/lai_viirs_v10` holds.
2. The step `laicorr` of `reproduce.py` repeats, from that cache, the hindcast that chose the correction
   (`analysis/v10/review/lai_viirs_correction.py`, `results/v10/lai_viirs_check/`).
3. Notebook 17 (`EE_PROJECT` in its cell "3 · Earth Engine") computes 2024 and 2025 with the corrected LAI of
   `src/megan_sc/laicorr.py` into `megan_sc_outputs_v10_lai_viirs`; its cell 11 must show that only the LAI differs from
   the production run. Each year takes about as long as a year of notebook 10, and the notebook resumes after a
   disconnection.
4. Cell 12 of notebook 17 (`INSTALL = True`) moves the MODIS years to `megan_sc_outputs_v10_modis_lai_2024_2025` and
   copies the corrected ones into the run; then run notebooks 11, 13 and 14 again. For the shipped run, notebook 19
   made this installation once (`notebooks/README.md`).

Notebook 22 repeats the check with the VIIRS of NOAA-20 and NOAA-21 (an Earthdata account, section 3.8).

## 6.8 The other notebooks of the example

| Notebook | What it does | Before running | Time |
|---|---|---|---|
| 12 | benchmark of the Python engine against the official MEGANv3.21 Fortran executables (`VALIDATION.md`) | `MEGANv3.21.tar.gz` in `MyDrive` (section 3.8); `RAW` = a run folder whose 2019 inputs build the test case (`data/run_v10` works) | 15 to 25 min |
| 15 | test S3 (soil-moisture factor) from the fix logs; test S1 from two short runs of notebook 10 | for S1: notebook 10 with `OUT_BASE = MyDrive/s1_test`, `YEARS = [2019]`, `MONTHS = [1, 4, 7, 10]`, `TEST_DAYS = 7`, once with `VARIANT = 'base'` and once with `'hourly'` | minutes |
| 18 | comparison with the BRAIN MEGAN emissions of 2019 (Hoinaski et al., 2024) | `BRAIN_output_urls.txt` (the BRAIN download links) and `MEGANv3.21.tar.gz` in Drive; no Earth Engine | download of the 364 daily files of 2019, resumable by month |
| 21 | ERA5-Land against the INMET stations, empirical lapse rate | `EE_PROJECT`; the INMET station lists | about 30 min |
| 20 | one-at-a-time sensitivity runs of 2019 | `EE_PROJECT`; notebook 21 first (its `lapse_rate.json`, else 6.5 K km⁻¹) | about 40 min for 28 days, about 15 h for the year |
| 23 | tests of the explanations of the differences from BRAIN | the outputs of 18, and those of 14, 20 and 21 when present | about 2 min |

In notebook 12, the test case shipped with the Fortran (south-east United States) does not depend on `RAW`, and its
rows should match `benchmark/results_colab/benchmark_summary.csv`; what to look for is S0 (the Fortran against the
reference output in the tarball) bit-identical, and every row of the patched engine marked PASS.
`python benchmark/make_table.py` prints the table of `VALIDATION.md` from the summary.

Executed copies of the Earth Engine notebooks show the project id in their output; check the outputs before sharing an
executed notebook.

## 6.9 After a run

1. `python reproduce.py check-data` does not apply to a new run (its checksums differ). Check provenance instead: every
   `<year>/provenance_v10.json` and every emission file must carry `megan_sc_fixes_sha256 = a933b29602929a82`, and
   check C0 of notebook 11 verifies it.
2. For a rerun of the example, recompute everything from it (`MEGAN_SC_RUN=<run folder> python reproduce.py all`,
   chapter 5, section 5.5) and look at the comparison.
3. When you publish results from a run, publish its configuration, its run folder and its results folder together,
   as this repository does for the example (chapter 12).
