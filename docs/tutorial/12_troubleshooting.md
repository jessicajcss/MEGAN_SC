# 12. Troubleshooting

## Level A (`reproduce.py`)

| What you see | Why | What to do |
|---|---|---|
| `missing run_v10/...` from `check-data` | the release archives were not unpacked | `python tools/data_release.py fetch <base-url>` (chapter 3) |
| `checksum differs ...` | an interrupted download, or a file changed after unpacking | delete the folder under `data/` and fetch again |
| `... annual_grids_cache.npz not found: run notebook 13 first` | a step that reads notebook 13 ran before it | `python reproduce.py nb13`, then the step again |
| `ModuleNotFoundError: pymannkendall` (or `statsmodels`, `PIL`) | package missing | `pip install -r requirements-reproduce.txt` |
| `AttributeError: ... has no attribute 'union_all'` | geopandas older than 1.0 | upgrade geopandas, or use `requirements-lock.txt` |
| a CSV marked `DIFFERENT` by `compare` | different library versions, or a changed input | read the first differing column printed on its line; install `requirements-lock.txt`; check `python reproduce.py check-data` and `python reproduce.py paths` |
| PNG marked `PNG DIFF` | fonts or plotting libraries differ from those that drew the published figure | harmless for the numbers |
| `NEW` lines in the comparison | files without a published copy: the grid cache of notebook 13 and the TIFF copies of the IJEST figures | expected |
| `OSError: [Errno 28] No space left on device` | the outputs need about 160 MB, plus copies of the TROPOMI cache | point `MEGAN_SC_OUT` to a disk with space; the step can be rerun alone |
| notebook 14 starts contacting Earth Engine during `reproduce.py nb14` | the cache is incomplete (it needs the 84 monthly files of `data/tropomi_hcho_v10/`) | check the folder and `check-data` |
| a script cannot find a folder | paths changed | `python reproduce.py paths`; set the `MEGAN_SC_*` variables (chapter 3, section 3.5) |

## Level B (Colab and Earth Engine)

| What you see | Why | What to do |
|---|---|---|
| `/content/drive holds local files written while Google Drive was not mounted` | an earlier cell wrote to `/content/drive` while the mount had failed | Runtime > Disconnect and delete runtime, then run again from the first cell |
| notebook 13: `expected annual_emis or 12 monthly files, found 0` | `ROOT` does not point to the run, or Drive is not mounted | check `ROOT` in the settings cell; the notebook lists what it found |
| `ee.Initialize` fails, or "project not found" | Earth Engine not enabled for the project, or a wrong id | check `EE_PROJECT` and the project's Earth Engine registration |
| `EE retry ...` lines in notebook 10 | transient Earth Engine errors; each request is retried four times | let it continue; if a request fails four times, run all cells again (the notebook resumes) |
| `Computation timed out` in notebook 14 | an interactive request was too large; the notebook uses Earth Engine batch exports instead | use the notebook shipped here (cell 1 prints "notebook 14 (Earth Engine batch exports)") |
| `FIX-1: no valid meteorology anywhere on the grid` | Earth Engine returned no valid ERA5-Land value for the step, for example for a date outside the archive | check `YEARS`, `MONTHS` and the availability of the requested dates in Earth Engine |
| notebook 12: `libiomp5.so not found` | the Intel OpenMP runtime the Fortran binaries need | cell 2 installs it with `pip install intel-openmp`; run on Linux x86-64 (Colab) |
| notebook 10 is slow | each day needs nine Earth Engine requests (eight 3-hourly fields and the daily mean temperature), plus the monthly LAI | run a few years per session; the notebook resumes after interruptions |

## Region wizard

`run_your_region/START_HERE.md` ends with a table of its messages (boundary file not found, `megan-sc` folder not
found, Earth Engine sign-in, empty grid, ecotype file needed outside Santa Catarina, empty test map).

## Reporting a problem

Open an issue on the repository with the command, the end of the log (`reproduced/logs/<step>.log`), the output of
`python reproduce.py paths`, and the versions of Python and of the packages (`pip freeze`).
