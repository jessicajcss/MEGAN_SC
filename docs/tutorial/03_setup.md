# 3. Setting up

Level A (reproduce) needs Python and about 1 GB of disk. Level B (rerun the model) needs, in addition, a Google account
with Earth Engine and Google Colab. This chapter prepares both.

## 3.1 Get the repository

```bash
git clone <repository-url> megan-sc
cd megan-sc
```

or, once the release is published, download its source archive from GitHub and unpack it. Every command of this tutorial is run from
this folder.

## 3.2 Python environment

Three ways, from the lightest to the most complete:

**pip, level A only** (Python 3.10 or later):

```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements-reproduce.txt
```

**pip, pinned versions.** `requirements-lock.txt` pins a tested set of versions (Python 3.10.12 on Linux). The numbers
are the same with it and with newer versions; figures can differ by a few pixels, because they depend on fonts and
library versions:

```bash
pip install -r requirements-lock.txt
```

**conda, everything** (Earth Engine client, notebooks, map projections), Python 3.11:

```bash
conda env create -f environment.yml
conda activate megan-sc
```

A `Dockerfile` builds the conda environment in a container (`docker build -t megan-sc .`); it has not been exercised yet.

The reproduction was made on Linux. The packages are available for Windows and macOS as well; on Windows, conda is the
simplest route for `rasterio`, `geopandas` and `netCDF4`.

## 3.3 Unit tests

```bash
pytest -q
```

Expected: `127 passed` (a few deprecation warnings from `shapely` are harmless). The tests cover the transcription of
the Fortran (activity factors, canopy, soil, daily means, speciation and units, input and output) and the eight
procedures FIX-1 to FIX-8 (`tests/test_fixes.py`). They need no download. They check the code, and say nothing about
the model against observations (chapter 10 does that).

## 3.4 Data

The model output, the TROPOMI fields, the test runs of test S1, the LAI cache of notebook 16 and the MODIS-LAI years 2024-2025 will be distributed as five zip files attached to the release (0.90 GB to download; see `DATA_AVAILABILITY.md`):

```bash
python tools/data_release.py fetch https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0
# or from the archived copy:      python tools/data_release.py fetch https://zenodo.org/records/<record>/files
# or from a folder with the zips:  python tools/data_release.py fetch /path/to/zips
python reproduce.py check-data
```

`fetch` downloads the five zips, checks each one against `SHA256SUMS_assets.txt` (or the copy in
`data/RELEASE_ASSETS_SHA256SUMS.txt`), and unpacks it into `data/`; on a slow disk this takes a few minutes.
`check-data` then checks every file under `data/` against `data/DATA_SHA256SUMS.txt` and should end with
`... data files checked, 0 missing or different`.

Disk space: about 1 GB for the unpacked data and about 0.2 GB for the outputs of `reproduce.py all`.

## 3.5 Paths

Every script finds its folders through `analysis/v10/_paths.py`. The defaults point inside the repository
(`data/run_v10`, `results/v10`, `reproduced`); environment variables override them, for example to keep the data on
another disk:

```bash
export MEGAN_SC_RUN=/big/disk/run_v10
export MEGAN_SC_OUT=/tmp/megan_reproduced
python reproduce.py paths       # prints the folders in use and whether they exist
```

| Variable | Default |
|---|---|
| `MEGAN_SC_RUN` | `data/run_v10` |
| `MEGAN_SC_BND` | `data/boundaries` |
| `MEGAN_SC_TROPOMI` | `data/tropomi_hcho_v10` |
| `MEGAN_SC_S1` | `data/s1_test` |
| `MEGAN_SC_LAIVIIRS` | `data/lai_viirs_v10` (monthly MODIS and VIIRS LAI of notebook 16) |
| `MEGAN_SC_RUN_MODIS_ARCHIVE` | `data/run_v10_modis_lai_2024_2025` (2024 and 2025 with the MODIS LAI) |
| `MEGAN_SC_RUN_MODIS` | `data/run_v10_modis_view` (MODIS-LAI record of every year, built by `tools/lai_variant_view.py`) |
| `MEGAN_SC_RESULTS` | `results/v10` |
| `MEGAN_SC_OUT` | `reproduced` |

## 3.6 For level B: Earth Engine and Colab

1. **Earth Engine.** Register a Google Cloud project for Earth Engine (free for research and non-commercial use) and
   note its id (for example `ee-yourname`). Notebooks 10, 14, 16, 17, 20 and 21 ask for it (`EE_PROJECT`); notebook 22 needs a NASA Earthdata account instead.
2. **Google Drive.** Copy the repository to `MyDrive/megan-sc`. The notebooks write to `MyDrive/megan_sc_outputs_v10`
   (notebook 10) and to folders inside it (notebooks 11 to 15).
3. **Colab.** Open the notebook from Drive with Google Colab, run the first cell (it mounts Drive), authorise Earth
   Engine when asked, then run the rest. A failed Drive mount stops the notebook on purpose: a notebook that went on
   would read an empty local folder.
4. **MEGANv3.21 code** (notebook 12 only): download `MEGANv3.21.tar.gz` from Zenodo (10.5281/zenodo.10526206) to
   `MyDrive`.

Chapter 6 runs these notebooks in order.
