# 🌎 Run MEGAN for your own region: start here

This folder lets you run the biogenic-emissions model (**MEGAN**) for **any state,
province, or country in the world**, Santa Catarina included, **without writing any
code**. You fill in a few boxes, upload the outline of your region, and press one button.

There are two files in this folder:

| File | What it is |
|------|------------|
| **`START_HERE.md`** | This guide. Read it once, top to bottom. |
| **`region_setup.ipynb`** | The wizard you actually run. Open it in **Google Colab**. |

> **Nothing here touches the original Santa Catarina setup.** The wizard only *reads* the
> shared model code and *writes* new files named after your region. Your existing outputs
> stay exactly as they are.

> The wizard runs the engine of notebook 10 with the same input and run procedures
> (FIX-1 to FIX-8 of `src/megan_sc/fixes.py`), active in its test (STEP 7) and in notebook 10.

---

## What you'll end up with

- A short **test map** of your region's isoprene emissions (proof it works), in minutes.
- A saved settings file, `config/<your-region>_run.yaml`.
- Then, with notebooks 10, 11 and 25 of `notebooks/`, the **full monthly inventory** of your region in
  `megan_sc_outputs_<your-region>/`, its checks and its results (this part takes hours, so start small).

---

## Before you start: three things (all free)

### 1. A boundary file of your region

This is just the **outline** of the area you want: a state, a province, a whole country.
It comes as either:

- a **shapefile**: a set of files that share a name (`myregion.shp`, `myregion.dbf`,
  `myregion.shx`, `myregion.prj`), all of them together, **or**
- a **GeoPackage**: a single `myregion.gpkg` file.

**Where to get one (free):**

- **Any country:** [geoBoundaries.org](https://www.geoboundaries.org/) or
  [Natural Earth](https://www.naturalearthdata.com/) → download the admin outline.
- **Brazilian states/municipalities:** IBGE
  ([ibge.gov.br](https://www.ibge.gov.br/en/geosciences/territorial-organization/territorial-meshes.html)),
  the same source the Santa Catarina files came from.
- **US states / counties:** the US Census TIGER/Line files.

Put the file(s) in your Google Drive inside the repository folder, for example in
`MyDrive/MEGAN_SC/data/regions/` (git ignores that folder, so your files stay apart from the repository). Note the full
path, for example `/content/drive/MyDrive/MEGAN_SC/data/regions/myregion.shp`: you'll paste it into the wizard.

> **Tip:** if you downloaded a whole-country file with many provinces and you only want
> one, that's fine: the outline can contain several pieces; the model uses the whole
> shape you give it. If you want just one province, open it in the free program **QGIS**,
> select that province, and *Export → Save Selected Features* as a new file.

### 2. A Google account with Earth Engine turned on

The model pulls weather, leaf-area, and land-cover data live from **Google Earth Engine**
(free for research/non-commercial use).

1. Open [console.cloud.google.com/earth-engine](https://console.cloud.google.com/earth-engine) with your Google
   account, create a Google Cloud project (or choose one) and register it for Earth Engine.
2. For research, teaching or other noncommercial use, answer the noncommercial-eligibility questions on
   [console.cloud.google.com/earth-engine/configuration](https://console.cloud.google.com/earth-engine/configuration);
   no billing account is needed for noncommercial use.
3. Note the **project id** (for example `my-ee-project`), shown in the project selector at the top of the Cloud
   console. You'll type it into the wizard as `EE_PROJECT`. Chapter 3 of the tutorial (`docs/tutorial/03_setup.md`)
   has the details, including Earth Engine on your own computer.

### 3. The `MEGAN_SC` folder in your Google Drive

The whole repository must be in your Drive as `MyDrive/MEGAN_SC` (a `MEGAN_SC` folder one level down in `MyDrive` also
works). Two ways to put it there: download the repository as a ZIP from GitHub (**Code > Download ZIP**), unpack it,
rename `MEGAN_SC-main` to `MEGAN_SC` and upload the folder (**New > Folder upload** in Google Drive); or, in a Colab
cell, mount Drive and run
`!git clone https://github.com/jessicajcss/MEGAN_SC.git /content/drive/MyDrive/MEGAN_SC`.
If you can open `MyDrive/MEGAN_SC/run_your_region/region_setup.ipynb` in Drive, you're ready.

---

## Running the wizard (step by step)

1. **Open the wizard in Colab.** In Google Drive, double-click `region_setup.ipynb`. If it
   asks how to open it, choose **Google Colaboratory**. (First time only: *Connect more
   apps → Colaboratory → Install*.)

2. **Fill in STEP 1.** This is the only cell you edit. Each line has a plain-English
   comment. The important ones:
   - `REGION_NAME`: a short nickname, no spaces (e.g. `parana`, `kenya`).
   - `SHAPEFILE_PATH`: the full path to your boundary file in Drive.
   - `COUNTRY_IS_BRAZIL`: `True` for Brazil, `False` for anywhere else (this picks the
     land-cover source; see *Brazil vs the rest of the world* below).
   - `START_YEAR`, `START_MONTH`, `END_YEAR`, `END_MONTH`: the first and the last month you want
     (both included; months are numbers, 1 = January). One month: the same year and month twice. One year:
     month 1 to month 12 of the same year. The cell stops with a message if the end comes before the start.
   - `EE_PROJECT`: your Earth Engine project id from thing #2.
   - `RESOLUTION_DEG`: map detail. `0.05` (≈5.5 km) is a good, fast default.

3. **Run everything.** Menu: **Runtime ▸ Run all**. Then just watch. When Earth Engine
   asks you to sign in (STEP 4), do it with the account from thing #2.

4. **Check the preview and the test map.** STEP 3 draws your region's outline so you can
   confirm it's the right shape. STEP 7 runs a short simulation with the full engine and draws an
   isoprene map. **If that map looks sensible, your region works.**

5. **Your settings are saved** automatically in STEP 8 as `config/<region>_run.yaml`, and STEP 9 prints what to type
   in notebooks 10, 11 and 25 for the full run.

If a step needs something from you, it stops and prints a plain-language message telling
you exactly what to fix. Change STEP 1 and run again.

---

## Brazil vs the rest of the world (the one difference that matters)

Everything the model needs is worldwide **except** the land-cover map:

- **In Brazil** (`COUNTRY_IS_BRAZIL = True`): uses **MapBiomas**, a 30-metre, Brazil-only
  land-cover map. Highest quality: it recognises Brazil's forests correctly.
- **Anywhere else** (`COUNTRY_IS_BRAZIL = False`): uses the **global MODIS** land-cover
  map (500 metres). It covers the entire planet, but it's coarser and can be less accurate
  in some forests. Everything still runs; just know the land-cover input is blunter.

A local land-cover map would need code that turns its classes into the four growth forms of MEGAN (trees, shrubs,
herbs, crops), as `src/megan_sc/gee/modis_lc.py` does for MODIS; the built-in global option works everywhere with no
extra steps.

---

## STEP 5 in detail: the "ecotype" file (only if you're outside Santa Catarina)

MEGAN uses an **ecotype map**: think of it as *which emission recipe applies where*. The
pipeline ships an ecotype tile that covers **only Santa Catarina**.

- **If your region is inside the Santa Catarina box**, the wizard reuses that tile
  automatically, with nothing to do.
- **If your region is anywhere else**, you need MEGAN's **global** ecotype file **once**.
  It's the same global dataset (called **EVT / EVT3d**) that the Santa Catarina tile was
  cut from. Get it, put it in your Drive, and set `GLOBAL_ECOTYPE_PATH` in STEP 1; the
  wizard then **clips it to your region for you** and remembers it.

**How to get the global ecotype (EVT) file:** download `EVT3d.zip` (11 MB) from the MEGAN 3.2
input files on Zenodo (Guenther, 2024, https://doi.org/10.5281/zenodo.10939297, licence CC BY 4.0)
and unzip it. It holds `EVT3d.nc` (1.5 GB), the global 1 km² map of MEGAN emission types. The
emission-factor and speciation tables already in `data/megan_inputs/` come from `EFP240120.zip` of
the same record. The wizard accepts the map as either a **NetCDF** (`.nc`) or a **GeoTIFF** (`.tif`).

> If you are unsure which file it is, point `GLOBAL_ECOTYPE_PATH` at it anyway: STEP 5 checks that it
> covers your region and stops with a message if it does not.

---

## The full inventory (after the test works)

The test in STEP 7 already used the **full** model engine. The complete monthly inventory for all your months and years
is computed by three notebooks of the `notebooks/` folder, opened in Colab the same way as the wizard:

1. **`10_production_run.ipynb`**: in its settings cell, set `RUN_CONFIG = 'config/<your-region>_run.yaml'` (the value
   printed by STEP 9) and your `EE_PROJECT`, then **Runtime ▸ Run all**. It writes one folder per year in
   `MyDrive/megan_sc_outputs_<your-region>/`, saves a checkpoint after every step and resumes after a disconnection
   (just run all cells again).
2. **`11_verify_run.ipynb`**: same `RUN_CONFIG`; six automatic checks of every year.
3. **`25_results_any_region.ipynb`**: same `RUN_CONFIG`; tables, maps, seasonal cycle and trends in
   `megan_sc_outputs_<your-region>/results/`.

Start with **one month** and a **coarse resolution** the first time: a full year can take hours. Chapter 8 of the
tutorial (`docs/tutorial/08_your_region.md`) lists the choices to review for your region, and chapter 9 explains how to
read the results.

---

## What this does *not* do (so there are no surprises)

- It **needs the internet and Earth Engine**: the data are pulled live, not stored locally.
- A **full run takes time** (hours per year at fine resolution). The *test* is quick.
- Outside Santa Catarina it **needs the global ecotype file once** (see STEP 5 above).
- Tables by unit (municipalities, counties) need a zones file (`ZONES_PATH` in STEP 1); without it, notebook 25 gives
  the totals of the whole region. The maps of `notebooks/13_santa_catarina_analysis.ipynb` are specific to Santa
  Catarina.
- It does not change or delete any existing Santa Catarina files.

---

## Troubleshooting

| Message you see | What to do |
|-----------------|-----------|
| `I can't find your boundary file` | The `SHAPEFILE_PATH` in STEP 1 is wrong, or the file isn't in Drive yet. For a shapefile, make sure the `.shp`, `.dbf`, `.shx`, `.prj` are all in the same folder. |
| `The MEGAN_SC folder was not found` | Put the whole repository at `MyDrive/MEGAN_SC` (see thing #3 above) and run all cells again. |
| `Earth Engine sign-in failed` | Check `EE_PROJECT` in STEP 1, and that you enabled Earth Engine for that Google account. |
| `no cells fell inside it` | Your file loaded but isn't a filled outline (polygon). Use an administrative *boundary/area* file, not a point or line file. |
| `Your region is OUTSIDE Santa Catarina …` | You need the global ecotype file: see **STEP 5** above. |
| The test map is all one colour / empty | Try a warmer `TEST_MONTH`, confirm your region has vegetation, and that land cover pulled sensible growth-form means (printed in STEP 5). |

---

*This wizard reuses the same tested model code as notebook 10, with the same procedures, and implements none of the
science anew. If the built-in test runs, your region is ready for notebook 10.*
