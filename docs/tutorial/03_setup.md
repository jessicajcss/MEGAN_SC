# 3. Setting up

The repository can be used in two places, and most users need both:

* **on a computer** with Python, to reproduce the Santa Catarina example from the shipped model output (chapter 5), to
  run the unit tests, and to analyse any run already computed;
* **in Google Colab**, with the repository in Google Drive, to run the notebooks that read Google Earth Engine:
  notebook 10 (the model itself, for any region and period), the region wizard, and the notebooks that read satellite
  data (14, 16, 17, 20, 21).

Earth Engine also works from a local Jupyter session (section 3.7.3) and from a server with a service account
(section 3.7.4). This chapter prepares each case.

## 3.1 Get the repository

With git:

```bash
git clone https://github.com/jessicajcss/MEGAN_SC.git MEGAN_SC
cd MEGAN_SC
```

Without git: on the GitHub page of the repository choose **Code > Download ZIP**, unpack it, and rename the folder
`MEGAN_SC-main` to `MEGAN_SC`. Every command of this tutorial is run from this folder, and every path is relative to it.

## 3.2 Put the repository in Google Drive (for Colab)

Every notebook opens with the same setup block (`notebooks/src/_setup_repo.py`). In Colab it mounts Google Drive and
looks for the repository in this order: the folder named by the environment variable `MEGAN_SC_REPO`;
`MyDrive/MEGAN_SC` (also `MyDrive/megan_sc` or `MyDrive/megan-sc`); a `MEGAN_SC` folder one level down in `MyDrive`;
`/content/MEGAN_SC`. A folder counts as the repository when it holds `src/megan_sc/fixes.py`. Outside Colab, the block
looks in the current folder and in its parent, so a notebook opened from `notebooks/` finds the repository by itself.

There are two ways to place the repository at `MyDrive/MEGAN_SC`.

**Upload the folder.** In a browser, open Google Drive, choose **New > Folder upload** and select the `MEGAN_SC` folder
of section 3.1. Google Drive for desktop (the synchronisation client) does the same by copying the folder into
`My Drive`.

**Clone it from Colab.** Open a new Colab notebook and run one cell:

```python
from google.colab import drive
drive.mount('/content/drive')
!git clone https://github.com/jessicajcss/MEGAN_SC.git /content/drive/MyDrive/MEGAN_SC
```

The notebooks that read the shipped model output (11, 13, 14 from its cache, 15 and 25 on the example) also need the
release archives inside the Drive copy (section 3.5). In Colab:

```python
%cd /content/drive/MyDrive/MEGAN_SC
!python tools/data_release.py fetch https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0
!python reproduce.py check-data
```

Notebook 10 and the region wizard need no release archive: they read the MEGAN tables, the ecotype tile and the
boundaries kept in git, and everything else from Earth Engine.

To open a notebook, double-click it in Google Drive and choose **Open with > Google Colaboratory** (the first time,
**Connect more apps > Colaboratory**). The first code cell asks for permission to mount Drive; a failed mount stops the
notebook on purpose, because a notebook that went on would read and write an empty local folder.

## 3.3 Python environment (local)

Three ways, from the lightest to the most complete.

**pip, reproduction only** (Python 3.10 or later):

```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements-reproduce.txt
```

**pip, pinned versions.** `requirements-lock.txt` pins a tested set of versions (Python 3.10.12 on Linux). The numbers
are the same with it and with newer versions; figures can differ by a few pixels or in the position of their panels, because they depend on fonts and
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

## 3.4 Unit tests

```bash
pytest -q
```

Expected: every test passes (a few deprecation warnings from `shapely` are harmless). The tests cover the
transcription of the Fortran (activity factors, canopy, soil, daily means, speciation and units, input and output),
the eight procedures FIX-1 to FIX-8 (`tests/test_fixes.py`), the region wizard, and a run of notebooks 10, 11 and 25
for another region with Earth Engine replaced by synthetic fields (`tests/test_offline_run.py`, through
`tools/offline_run.py`). They need no download and no Earth Engine account. They check the code, and say nothing about
the model against observations (chapter 10 does that).

## 3.5 The release archives

The model output of the Santa Catarina example, the TROPOMI fields, the test runs of test S1, the LAI cache of
notebook 16 and the MODIS-LAI years 2024-2025 are distributed as five zip files attached to the release v1.0.0 of the
repository (0.90 GB to download; `DATA_AVAILABILITY.md` lists them):

```bash
python tools/data_release.py fetch https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0
# or from a folder that holds the five zips:  python tools/data_release.py fetch /path/to/zips
python reproduce.py check-data
```

`fetch` downloads the five zips, checks each one against `SHA256SUMS_assets.txt` (or the copy in
`data/RELEASE_ASSETS_SHA256SUMS.txt`), and unpacks it into `data/`; on a slow disk this takes a few minutes.
`check-data` then checks every file under `data/` against `data/DATA_SHA256SUMS.txt` and should end with
`... data files checked, 0 missing or different`.

Disk space: about 1 GB for the unpacked data and about 0.2 GB for the outputs of `reproduce.py all`.

## 3.6 Paths

Every script of `analysis/v10/` and `reproduce.py` find their folders through `analysis/v10/_paths.py`. The defaults
point inside the repository (`data/run_v10`, `results/v10`, `reproduced`); environment variables override them, for
example to keep the data on another disk:

```bash
export MEGAN_SC_RUN=/big/disk/run_v10
export MEGAN_SC_OUT=/tmp/megan_reproduced
python reproduce.py paths       # prints the folders in use
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

The notebooks take their folders from their settings cell instead (chapter 6).

## 3.7 Google Earth Engine

Notebook 10 reads ERA5-Land, the MODIS LAI and the land cover through Earth Engine (Gorelick et al., 2017), and so do
the wizard and notebooks 14, 16, 17, 20 and 21. Every Earth Engine request is made on behalf of a Google Cloud project
registered for Earth Engine. The steps below follow the Earth Engine guides "Authentication and Initialization",
"Access" and "Access control" (https://developers.google.com/earth-engine/guides/auth,
https://developers.google.com/earth-engine/guides/access and https://developers.google.com/earth-engine/guides/access_control).
A person who signs in with the account that owns the project already has these permissions.

### 3.7.1 Register a Cloud project (once)

1. Sign in with a Google account and open https://console.cloud.google.com/earth-engine (the Earth Engine page of the
   Google Cloud console). Create a project or choose an existing one, and register it for Earth Engine.
2. For noncommercial use (research, education, non-profit work), answer the noncommercial-eligibility questions on
   https://console.cloud.google.com/earth-engine/configuration. Noncommercial projects need no billing account. Since
   26 September 2025 every noncommercial project must have answered these questions to keep using Earth Engine.
3. Note the **project id** (for example `my-ee-project`). It goes into the setting `EE_PROJECT` of a notebook, or into
   `earth_engine.project` of a run configuration (chapter 6). The placeholder `your-ee-project` of the shipped files
   stops the notebook with a message until it is replaced.

Every notebook calls `megan_sc.gee.initialize(project)`, which first tries the credentials stored by an earlier sign-in
and asks for a sign-in only when there are none.

### 3.7.2 In Google Colab

Nothing to install. The first time a notebook reaches Earth Engine, Colab asks to allow access with your Google account;
accept with the account that owns the project. The credentials last for the Colab session.

### 3.7.3 On your own computer (Windows, macOS or Linux)

```bash
pip install earthengine-api          # included in environment.yml
earthengine authenticate             # opens a browser; stores the credentials in ~/.config/earthengine/credentials
earthengine set_project my-ee-project    # optional: default project of the command-line tool
```

`earthengine authenticate --force` signs in again with another account. Inside a Jupyter notebook,
`ee.Authenticate()` prints a link to open in a browser and asks for the code it returns; on a remote machine without a
browser pass `auth_mode='notebook'` (`megan_sc.gee.initialize(project, auth_mode='notebook')`). With the credentials
stored, open the notebooks of `notebooks/` in Jupyter from the repository: the setup block finds the repository from
the notebook folder, and the outputs go next to the repository (setting `OUT_BASE`).

### 3.7.4 With a service account (servers and scheduled runs)

For runs without a person to sign in, use a service account of the same Cloud project:

1. In the Cloud console, **IAM & Admin > Service accounts > Create service account**, then **Keys > Add key > JSON**,
   and keep the JSON file outside the repository.
2. In **IAM**, give the service account the role **Earth Engine Resource Writer**, which allows interactive
   computations and export tasks (notebook 10 makes the first, notebook 14 the second), and **Service Usage Consumer**:
   without the permission `serviceusage.services.use` on the project, `ee.Initialize` fails (Earth Engine guide
   "Access control", https://developers.google.com/earth-engine/guides/access_control).
3. Write the key path into the configuration: `earth_engine: {project: my-ee-project, service_account_key: /secure/path/key.json}`.
   `megan_sc.gee.initialize` then builds `ee.ServiceAccountCredentials` from the key and initialises Earth Engine with
   it; no browser is involved.

### 3.7.5 What Earth Engine costs in time

A month of notebook 10 makes about nine requests per simulated day (eight 3-hourly ERA5-Land fields and the 24 hourly
values of the daily mean temperature, sampled as one daily image), plus the monthly LAI and the land cover of the year.
A year of the Santa Catarina example (144 × 224 cells) takes hours of wall-clock time in Colab; the notebook saves a
checkpoint after every step and resumes after a disconnection (chapter 6). The TROPOMI extraction of notebook 14 used
41.2 EECU-hours (Earth Engine compute units) for 84 months.

## 3.8 Other accounts and files, for some notebooks only

* **MEGANv3.21 Fortran code** (notebook 12, the benchmark against the official executables): download
  `MEGANv3.21.tar.gz` from Zenodo (Guenther, 2024b, https://doi.org/10.5281/zenodo.10526206) to `MyDrive`.
* **Global ecotype map** (the region wizard, for a region outside the Santa Catarina tile): `EVT3d.zip` of the MEGAN 3.2
  input files (Guenther, 2024a, https://doi.org/10.5281/zenodo.10939297), chapter 8.
* **NASA Earthdata account** (notebook 22 only, the VIIRS LAI of NOAA-20 and NOAA-21 read as HDF5 files).
