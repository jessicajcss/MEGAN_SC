# 8. Your own region and period

The engine applies to any land region of the world and to any period covered by its inputs: the grid, the inputs and
the reporting mask all follow from a boundary file and a run configuration. The route has four steps, all in Google
Colab with the repository at `MyDrive/MEGAN_SC` (chapter 3):

1. `run_your_region/region_setup.ipynb`, the wizard, writes the configuration `config/<region>_run.yaml` and tests it;
2. `notebooks/10_production_run.ipynb` runs the model with that configuration;
3. `notebooks/11_verify_run.ipynb` checks the run;
4. `notebooks/25_results_any_region.ipynb` writes the tables, maps and trends.

`run_your_region/START_HERE.md` is the plain-language guide to the wizard, with the sources of free boundary files.

## 8.1 What you need

* **A boundary file** of the region: a shapefile (`.shp` with its `.dbf`, `.shx` and `.prj`) or a GeoPackage
  (`.gpkg`), with one or several polygons, in any coordinate system. Put it in your Drive, for example in
  `MyDrive/MEGAN_SC/data/regions/` (a folder that git ignores, so your files never mix with the repository).
* **Optionally, a zones file**: smaller units inside the region (municipalities, counties) for tables by unit.
* **An Earth Engine project** (chapter 3, section 3.7).
* **Outside the Santa Catarina tile, the global ecotype map** of the MEGAN 3.2 input files: `EVT3d.zip` (11 MB, with
  `EVT3d.nc`, 1.5 GB, the global 1 km map of MEGAN emission types) from Zenodo (Guenther, 2024a,
  https://doi.org/10.5281/zenodo.10939297, licence CC BY 4.0). The emission-factor and speciation tables of
  `data/megan_inputs/` come from `EFP240120.zip` of the same record and serve every region.

## 8.2 The wizard, step by step

Only STEP 1 is edited:

| Setting | Meaning |
|---|---|
| `REGION_NAME` | a short name without spaces, used in file names |
| `SHAPEFILE_PATH` | the boundary file |
| `COUNTRY_IS_BRAZIL` | `True`: MapBiomas Collection 10 (30 m, Brazil only); `False`: MODIS MCD12Q1 (500 m, worldwide) |
| `START_YEAR`, `START_MONTH`, `END_YEAR`, `END_MONTH` | the first and the last month, both included |
| `EE_PROJECT` | your Earth Engine project id |
| `RESOLUTION_DEG` | cell size: 0.05° (about 5.5 km) is a fast default, 0.025° about four times slower, 0.1° for a quick look |
| `GLOBAL_ECOTYPE_PATH` | the global ecotype map (`.nc` or `.tif`), needed outside the Santa Catarina tile |
| `ZONES_PATH`, `ZONES_ID` | optional zones file and its identifier column |
| `CANOPY_TYPE`, `TEST_MONTH`, `TEST_DAYS` | advanced; `auto` picks tropical broadleaf (type 2) within 23.5° of the equator and temperate broadleaf (type 3) elsewhere, and a warm month of the period for the test |

Then **Runtime > Run all**. The wizard:

1. checks the period and lists its months (STEP 1);
2. mounts Drive and finds the repository (STEP 2);
3. reads the boundary in any coordinate system, builds a grid over it and the mask of cells whose centre lies inside,
   and draws a preview (STEP 3);
4. signs in to Earth Engine and checks that ERA5-Land and the MODIS LAI exist for the first and last day of the period
   (STEP 4);
5. reads the growth-form fractions of the land cover and selects the ecotype map: the Santa Catarina tile when the
   region lies inside it, otherwise a clip of the global map saved in `data/regions/<region>/`; then fills the
   no-data pixels of the ecotype map (FIX-5) (STEP 5);
6. builds the emission factors and the canopy types (STEP 6);
7. runs the engine of notebook 10, with the procedures FIX-1 to FIX-8, over `TEST_DAYS` days of the test month every
   6 hours and draws an isoprene map (STEP 7);
8. writes `config/<region>_run.yaml`, in the same form as the configuration of the example (STEP 8);
9. prints the value of `RUN_CONFIG` to set in notebooks 10, 11 and 25 (STEP 9).

If a step needs something from you, it stops with a message that names what to fix in STEP 1.

## 8.3 The full run

In notebooks 10, 11 and 25 set `RUN_CONFIG = 'config/<region>_run.yaml'` and `EE_PROJECT` (notebook 10), and run them in
that order (chapter 6, sections 6.3, 6.4 and 6.6). The output goes to `MyDrive/megan_sc_outputs_<region>/`, one
folder per year, and the results of notebook 25 to `MyDrive/megan_sc_outputs_<region>/results/`.

The time grows with the number of days and of grid cells. The first time, run one month at a coarse resolution; then
one year; then the whole period, a few years per session if needed (the notebook resumes).

## 8.4 Choices to review for another region

* **Canopy type of trees.** Santa Catarina used temperate broadleaf (type 3). The wizard picks tropical broadleaf
  (type 2) within 23.5° of the equator; check `datasets.vegetation.tree_canopy_type` in the configuration. A region
  dominated by conifers may need type 1 (needleleaf).
* **Land cover.** MapBiomas (Brazil) gives 30 m growth forms from its own classes; MODIS MCD12Q1 version 6.1 gives
  500 m classes of the LC_Type5 (plant functional type) legend, which `megan_sc.gee.modis_lc` maps to trees (classes
  1 to 4), shrubs (5), herbaceous vegetation (6) and crops (7 and 8); water, urban, snow and barren cells (0, 9, 10, 11)
  carry no vegetation. Coarse classes can misplace forest edges and mosaics of crops and trees.
* **Land-cover years.** In Earth Engine, MapBiomas Collection 10 covers 1985-2024 and MCD12Q1 version 6.1 covers
  2001-2024. A year outside the range uses the nearest year of the product, as Santa Catarina did for 2025, and the
  year used is written in `provenance_v10.json` (`land_cover_year`).
* **Wilting point.** The run used a constant 0.12 m³ m⁻³. `src/megan_sc/gee/soil_texture.py` and
  `soil.wilting_point_from_texture` can derive it from soil texture, which matters where soils dry out (section 1.6).
* **Soil-NO season.** `soil.growseason` treats the tropics (year-round), the southern temperate zone (November-May, the
  season adjusted to the Southern Hemisphere) and the northern temperate zone; check that the season suits your
  agriculture (section 1.8).
* **Canopy loss and production factor.** Notebook 10 writes isoprene as the MEGANv3.21 code computes it, without ρ.
  Notebooks 11, 13 and 25 multiply it by ρ = 0.96 (Guenther et al., 2006) when they read the run (section 1.1); state
  which convention your numbers follow.
* **Coastlines and islands.** FIX-1 fills cells without ERA5-Land data from the nearest valid land cell. Islands far from
  land take the meteorology of a distant coast; look at the fix logs (`filled_met_steps`) and decide whether to report
  those cells.
* **Resolution.** ERA5-Land is 0.1°; a grid much finer than 0.025° adds land-cover detail but no meteorological detail.
* **Time of day.** Inputs and outputs are in UTC, the time of ERA5-Land. Solar elevation is computed for each cell in
  local mean solar time (UTC hour + longitude/15, as in MEGANv3.21), so the daily cycle of light is right in any time
  zone. The daily mean temperature (00-23 UTC, as `daymet.f`), the days and the calendar months of the output stay in
  UTC: in Santa Catarina the UTC day starts at about 20:30 local solar time, at 120° E it starts at 08:00. Radiation
  de-accumulated from ERA5-Land is the mean of the hour before each step, while temperature is the value at the step.
  `docs/TIME_CONVENTIONS.md` gives the sources and what each of these does to the results of the example.
* **Period.** In the Earth Engine catalog, ERA5-Land hourly (`ECMWF/ERA5_LAND/HOURLY`) starts on 1 January 1950 and
  the MODIS LAI (`MODIS/061/MCD15A3H`) on 4 July 2002, so the earliest complete year of a run is 2003. Both collections
  are extended with a delay of weeks to months; STEP 4 of the wizard checks that the first and the last day of the
  period are available.

## 8.5 What to report

Report the configuration file, the engine version (`fixes.py` hash), the choices of section 8.4, the result of the
checks of notebook 11, and the limits of chapter 9 (section 9.7) that apply to your region. Without an
observational comparison, present the inventory as a model estimate; the TROPOMI test of chapter 10 can be run for
another region by moving the background box and the grid.
