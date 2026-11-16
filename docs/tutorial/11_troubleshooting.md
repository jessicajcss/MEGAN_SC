# 11. Troubleshooting

## 11.1 Reproducing the example (`reproduce.py`)

| What you see | Why | What to do |
|---|---|---|
| `missing run_v10/...` from `check-data` | the release archives were not unpacked | `python tools/data_release.py fetch <base-url>` (chapter 3) |
| `checksum differs ...` | an interrupted download, or a file changed after unpacking | delete the folder under `data/` and fetch again |
| `... annual_grids_cache.npz not found: run notebook 13 first` | a step that reads notebook 13 ran before it | `python reproduce.py nb13`, then the step again |
| `ModuleNotFoundError: pymannkendall` (or `statsmodels`, `PIL`) | package missing | `pip install -r requirements-reproduce.txt` |
| `AttributeError: ... has no attribute 'union_all'` | geopandas older than 1.0 | upgrade geopandas, or use `requirements-lock.txt` |
| a CSV marked `DIFFERENT` by `compare` | different library versions, or a changed input | read the first differing column printed on its line; install `requirements-lock.txt`; check `python reproduce.py check-data` and `python reproduce.py paths` |
| PNG marked `PNG DIFF` | fonts or plotting libraries differ from those that drew the published figure | harmless for the numbers |
| `NEW` lines in the comparison | files without a published copy: the grid cache of notebook 13 | expected |
| `OSError: [Errno 28] No space left on device` | the outputs need about 160 MB, plus copies of the TROPOMI cache | point `MEGAN_SC_OUT` to a disk with space; the step can be rerun alone |
| notebook 14 starts contacting Earth Engine during `reproduce.py nb14` | the cache is incomplete (it needs the 84 monthly files of `data/tropomi_hcho_v10/`) | check the folder and `check-data` |
| a script cannot find a folder | paths changed | `python reproduce.py paths`; set the `MEGAN_SC_*` variables (chapter 3, section 3.6) |

## 11.2 Running notebooks (Colab, Jupyter, Earth Engine)

| What you see | Why | What to do |
|---|---|---|
| `The MEGAN_SC folder was not found` | the setup block found no folder with `src/megan_sc/fixes.py` | in Colab, put the repository at `MyDrive/MEGAN_SC` (chapter 3, section 3.2); elsewhere, open the notebook from the `notebooks/` folder of the repository, or set `os.environ["MEGAN_SC_REPO"]` before the first cell |
| `Earth Engine needs the id of a Google Cloud project ...` | `EE_PROJECT` is empty, or the configuration still holds `your-ee-project` | set `EE_PROJECT` in the settings cell, or `earth_engine.project` in the configuration (chapter 3, section 3.7.1) |
| `ee.Initialize` fails with a permission or "not registered" message | the project is not registered for Earth Engine, or the noncommercial questions were not answered, or a service account lacks the roles of section 3.7.4 | open https://console.cloud.google.com/earth-engine with the account that owns the project |
| `region boundary not found` (notebook 10), or a file-not-found error when notebook 11 or 25 reads the boundary | the path of `region.boundary` (or `domain.shapefile`) in the configuration is wrong | paths inside the repository are relative to it; a file in Drive needs its full path `/content/drive/MyDrive/...` |
| notebook 11 prints `years found: []`, or notebook 25 stops with `no monthly emission files` | `RUN_CONFIG` or `OUT_BASE` differs from those of the notebook-10 run, or only `TEST_DAYS` files exist | set the same `RUN_CONFIG` and `OUT_BASE`, or `NEW_ROOT` / `RUN_ROOT` to the run folder |
| `/content/drive holds local files written while Google Drive was not mounted` | an earlier cell wrote to `/content/drive` while the mount had failed | Runtime > Disconnect and delete runtime, then run again from the first cell |
| notebook 13: `expected annual_emis or 12 monthly files, found 0` | `ROOT` does not point to the run, or Drive is not mounted | check `ROOT` in the settings cell; the notebook lists what it found |
| `EE retry ...` lines in notebook 10 | transient Earth Engine errors; each request is retried four times | let it continue; if a request fails four times, run all cells again (the notebook resumes) |
| `Computation timed out` in notebook 14 | an interactive request was too large; the notebook uses Earth Engine batch exports instead | use the notebook shipped here (cell 1 prints "notebook 14 (Earth Engine batch exports)") |
| `FIX-1: no valid meteorology anywhere on the grid` | Earth Engine returned no valid ERA5-Land value for the step, for example for a date outside the archive | check the period of the configuration, `YEARS`, `MONTHS` and the availability of the requested dates in Earth Engine (chapter 8, section 8.4) |
| notebook 12: `libiomp5.so not found` | the Intel OpenMP runtime the Fortran binaries need | cell 2 installs it with `pip install intel-openmp`; run on Linux x86-64 (Colab) |
| notebook 10 is slow | each day needs nine Earth Engine requests (eight 3-hourly fields and the daily mean temperature), plus the monthly LAI | run a few years per session; the notebook resumes after interruptions |

## 11.3 Region wizard

`run_your_region/START_HERE.md` ends with a table of its messages (boundary file not found, repository folder not
found, Earth Engine sign-in, empty grid, ecotype file needed outside Santa Catarina, empty test map).

## 11.4 Reporting a problem

Open an issue on the repository with the command, the end of the log (`reproduced/logs/<step>.log`), the output of
`python reproduce.py paths`, and the versions of Python and of the packages (`pip freeze`).
