# 11. Your own region

The engine applies to any region: the grid, the inputs and the reporting mask all follow from a boundary file. Two routes lead to an inventory for another region: the setup wizard, which is quick and needs no
programming, and notebook 10 adapted to the region, which keeps every feature of the production run.

## 11.1 The setup wizard

`run_your_region/START_HERE.md` walks through it; `run_your_region/region_setup.ipynb` is the notebook (Colab). In
short: fill in the first cell (region name, boundary file, first and last month of the period, Earth Engine project,
resolution), run all cells,
and the wizard

1. builds the grid over the boundary and the mask of cells inside it;
2. reads the land cover (MapBiomas in Brazil, MODIS MCD12Q1 elsewhere) and turns it into growth-form fractions;
3. selects the ecotype map and fills its no-data pixels (FIX-5): inside the state, the Santa Catarina tile; elsewhere,
   a clip of `EVT3d.nc`, the global map of the MEGAN 3.2 input files (Guenther, 2024a), downloaded once;
4. builds the emission factors and the canopy types;
5. runs a short test (two days of a warm month) with the engine of notebook 10 and draws an isoprene map;
6. writes `config/<region>_run.yaml`;
7. optionally (`RUN_FULL = True`) runs every month of the chosen period, 3-hourly, and writes monthly NetCDF files and
   statistics.

The engine of steps 5 and 7 is the one of notebook 10, with the procedures FIX-1 to FIX-8, except that the daily mean
temperature is the mean of the sampled steps of the day, while notebook 10 averages the 24 hourly values (one more
Earth Engine request per day).

## 11.2 Notebook 10 for another region

For checkpoints, input GeoTIFFs, fix logs and municipal tables, adapt notebook 10. The edits, by cell:

| Cell | Edit |
|---|---|
| 1 | `OUT_BASE` to your Drive folder; `OLD_ROOT = ''` (no folder of cached growth forms to reuse) |
| 4 | `YEARS`; keep `SAMPLE_HOURS = 3` and `VARIANT = 'base'` |
| 5 | load the configuration written by the wizard: replace the path in `_cfgp = MEGAN_DIR + '/config/example_config.yaml'` by `MEGAN_DIR + '/config/<region>_run.yaml'`. The notebook sets `cfg['domain']['mode'] = 'bbox'`, and the wizard has written the bounding box of the region into the configuration, so the grid covers the region. The ecotype path of that configuration is the filled ecotype of the wizard; the cell copies it to a file named `sc_ecotype_filled.nc` in the same folder (variable `_eco_f`), so give that file a name of your region to avoid mixing regions |
| 6 | `SC_GPKG` to your boundary file (the reporting mask) and `SHP` to your administrative units, if any (the identifier column is `CD_MUN` when present, otherwise the first column) |
| 8 | outside Brazil, replace the MapBiomas call in `setup_year` (`mapbiomas.growth_form_fractions_gridded`) by the MODIS land-cover call used in STEP 9 of the wizard (`modis_lc.growth_form_fractions` and `mapbiomas.sample_image_tiled`), and adapt `VEG_YEAR` to the years of that product |

Then run notebook 11 (checks C0 to C5) on the new run (`NEW_ROOT` = your output, `SC_GPKG` = your boundary). Notebook
13 is specific to the Santa Catarina documents (mesoregions, municipalities, previous inventory); for another region,
start from the monthly statistics of the run (`<year>/stats/stats_<year>_<mm>.csv`) or adapt its cells.

## 11.3 Choices to revisit for another region

* **Canopy type of trees.** Santa Catarina used temperate broadleaf (type 3). The wizard picks tropical broadleaf
  (type 2) within 23.5° of the equator; check `datasets.vegetation.tree_canopy_type` in the configuration.
* **Wilting point.** The run used a constant 0.12 m³ m⁻³. `src/megan_sc/gee/soil_texture.py` and
  `soil.wilting_point_from_texture` can derive it from soil texture, which matters where soils dry out.
* **Soil-NO season.** `soil.growseason` treats the tropics (year-round), the southern temperate zone (November-May,
  the season adjusted to the Southern Hemisphere) and the northern temperate zone; check that the season suits your
  agriculture.
* **Canopy loss and production factor.** Notebook 10 and the wizard write isoprene as the MEGANv3.21 code computes it,
  without ρ. The Santa Catarina analysis multiplies it by ρ = 0.96 (Guenther et al., 2006) when it reads the run
  (chapter 1, section 1.1); apply the same factor when you report isoprene, or state that you did not.
* **Land-cover years.** MapBiomas covers 1985-2024 in Collection 10; MODIS MCD12Q1 has its own range. The run reuses
  the last available year for later years, as Santa Catarina did for 2025.
* **Coastlines.** FIX-1 fills cells without ERA5-Land data from the nearest valid cell. Islands far from land take the
  meteorology of a distant coast; look at the fix logs (`filled_met_steps`) and decide whether to report those cells.
* **Resolution.** ERA5-Land is 0.1°; a grid much finer than 0.025° adds land-cover detail but no meteorological detail.
* **Time of day.** Inputs and outputs are in UTC, the time of ERA5-Land. Solar elevation is computed for each cell in
  local mean solar time (UTC hour + longitude/15, as in MEGANv3.21), so the daily cycle of light is right in any time
  zone. The daily mean temperature (00-23 UTC, as `daymet.f`), the days and the calendar months of the output stay in
  UTC: in Santa Catarina the UTC day starts at about 20:30 local solar time, at 120° E it starts at 08:00. Radiation
  de-accumulated from ERA5-Land is the mean of the hour before each step, while temperature is the value at the step.
  `docs/TIME_CONVENTIONS.md` gives the sources and what each of these does to the Santa Catarina results.

## 11.4 What to report

Report the engine version (`fixes.py` hash), the settings above, the checks of notebook 11, and the limits of section
1.10 that apply to your region. Without an observational comparison, present the inventory as a model estimate; the
TROPOMI test of chapter 10 can be run for any region by changing the background box and the grid.
