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
> (FIX-1 to FIX-8 of `src/megan_sc/fixes.py`), active in STEP 7 and STEP 9.

---

## What you'll end up with

- A short **test map** of your region's isoprene emissions (proof it works), in minutes.
- A saved settings file, `config/<your-region>_run.yaml`.
- Optionally, a **full monthly inventory** for your region in
  `megan_sc_outputs_<your-region>/` (this part takes hours, so start small).

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

Put the file(s) in your Google Drive inside the pipeline folder, for example
`MyDrive/megan-sc/data/boundaries/`. Note the full path: you'll paste it into the wizard.

> **Tip:** if you downloaded a whole-country file with many provinces and you only want
> one, that's fine: the outline can contain several pieces; the model uses the whole
> shape you give it. If you want just one province, open it in the free program **QGIS**,
> select that province, and *Export → Save Selected Features* as a new file.

### 2. A Google account with Earth Engine turned on

The model pulls weather, leaf-area, and land-cover data live from **Google Earth Engine**
(free for research/non-commercial use).

1. Go to [earthengine.google.com](https://earthengine.google.com/) and click **Get
   Started** / **Sign Up**. Approve the one-time access with your Google account.
2. You'll get a **project id** that looks like `ee-yourname`. You can see it at
   [console.cloud.google.com](https://console.cloud.google.com/) (top project dropdown),
   or in the Earth Engine Code Editor under *Assets*. Keep it handy: you'll type it into
   the wizard as `EE_PROJECT`.

### 3. The `megan-sc` folder in your Google Drive

The whole pipeline folder must be at `MyDrive/megan-sc` (or
`MyDrive/MEGAN_072026/megan-sc`). If you can open it in Drive, you're ready.

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

5. **Your settings are saved** automatically in STEP 8 as `config/<region>_run.yaml`.

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

If you work outside Brazil and your country has an excellent local land-cover map, that can
be plugged in later (ask whoever maintains the pipeline), and the built-in global option
works everywhere with no extra steps.

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

The test in STEP 7 already used the **full** model engine. To make the complete monthly
inventory for all your months and years, you have two options:

**Easy (recommended):** in the wizard, go to **STEP 9**, change `RUN_FULL = False` to
`RUN_FULL = True`, and run that cell. It writes a NetCDF of tonnes and a statistics CSV per
month into `megan_sc_outputs_<region>/<year>/`. Start with **one year** and a **coarse
resolution** the first time: a full year can take hours.

**Advanced (full-featured):** `notebooks/10_production_run.ipynb`, the notebook
that produced the Santa Catarina inventory, adds checkpointing (it resumes after a
disconnect), per-municipality tables, input GeoTIFFs and the fix logs. The edits that
point it at your region are listed in `docs/tutorial/11_your_region.md`.

Because the wizard wrote your box into the config, the run will target your region either
way.

---

## What this does *not* do (so there are no surprises)

- It **needs the internet and Earth Engine**: the data are pulled live, not stored locally.
- A **full run takes time** (hours per year at fine resolution). The *test* is quick.
- Outside Santa Catarina it **needs the global ecotype file once** (see STEP 5 above).
- The per-municipality **choropleth figures** in the analysis notebooks assume Brazilian
  IBGE units. For another country, supply your own administrative outline for those
  figures, or skip them; the emissions themselves don't depend on them.
- It does not change or delete any existing Santa Catarina files.

---

## Troubleshooting

| Message you see | What to do |
|-----------------|-----------|
| `I can't find your boundary file` | The `SHAPEFILE_PATH` in STEP 1 is wrong, or the file isn't in Drive yet. For a shapefile, make sure the `.shp`, `.dbf`, `.shx`, `.prj` are all in the same folder. |
| `Could not find the 'megan-sc' folder` | Put the whole `megan-sc` folder at `MyDrive/megan-sc` and re-run STEP 2. |
| `Earth Engine sign-in failed` | Check `EE_PROJECT` in STEP 1, and that you enabled Earth Engine for that Google account. |
| `no cells fell inside it` | Your file loaded but isn't a filled outline (polygon). Use an administrative *boundary/area* file, not a point or line file. |
| `Your region is OUTSIDE Santa Catarina …` | You need the global ecotype file: see **STEP 5** above. |
| The test map is all one colour / empty | Try a warmer `TEST_MONTH`, confirm your region has vegetation, and that land cover pulled sensible growth-form means (printed in STEP 5). |

---

*This wizard reuses the same tested model code as the main pipeline (notebook 10, with the
same procedures) and implements none of the science anew. If the built-in test runs, your region is
ready.*
