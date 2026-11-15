# `data/`: inputs and model output

The small inputs are kept in git. The five large folders (`run_v10/` 2,075 files, `tropomi_hcho_v10/` 84, `s1_test/` 72, `lai_viirs_v10/` 5, `run_v10_modis_lai_2024_2025/` 180; 2,416 files, distributed as five zip files of 0.90 GB together) are unpacked here by
`python tools/data_release.py fetch <base-url>`; `git` ignores them (`.gitignore`). `DATA_SHA256SUMS.txt` lists the
SHA-256 of every file below (2,435 files: all except this README and the two checksum lists), and
`python reproduce.py check-data` verifies them; `RELEASE_ASSETS_SHA256SUMS.txt` holds the checksums of the five zip files. Sources and licences are in `../DATA_AVAILABILITY.md`.

## Conventions shared by every grid

* **Grid**: 0.025° in latitude and longitude, 144 rows × 224 columns, cell centres from 29.4875° S to 25.9125° S and
  from 53.8875° W to 48.3125° W (the bounding box of Santa Catarina). Geographic coordinates, WGS 84 (EPSG:4326).
* **Row order**: NetCDF files and NumPy arrays (`.nc`, `.npz`) store row 0 at the **south** (latitude ascending), as the
  engine does. GeoTIFF files are north-up (row 0 at the north), as GIS software expects; flip them (`numpy.flipud`)
  before combining them with the NetCDF fields.
* **State mask**: a cell belongs to Santa Catarina when its centre lies inside `boundaries/sc_state_ibge.gpkg`; 13,970
  cells (43.3 % of the box, 95,969 km² of grid area). The engine computes cell areas on the authalic sphere
  (R = 6,371,007.181 m, `src/megan_sc/domain.py`), and the analysis scripts on the sphere of mean radius
  (R = 6,371,008.8 m); the two areas differ by less than 10⁻⁶.
* **Emission classes**: the 19 MEGAN classes, in this order: ISOP, MBO, MT_PINE, MT_ACYC, MT_CAMP, MT_SABI, MT_AROM,
  NO, SQT_HR, SQT_LR, MEOH, ACTO, ETOH, ACID, LVOC, OXPROD, STRESS, OTHER, CO. BVOC = the 17 classes other than NO and
  CO.

## In git

| Path | Content |
|---|---|
| `boundaries/sc_state_ibge.gpkg` | state polygon (IBGE) |
| `boundaries/SC_Municipios_2024.gpkg` | 295 municipalities, 2024 mesh, simplified; fields `CD_MUN`, `NM_MUN`, `NM_RGINT` (intermediate geographic region), `NM_UF`, `SIGLA_UF`, `AREA_KM2` |
| `boundaries/sc_mesoregions_2024.gpkg` | six mesoregions with the islands; `name_meso`, `n_mun`, `area_ibge_km2`, `forest_pct` (MapBiomas 2019) |
| `boundaries/sc_municipios_2024_mesorregiao.csv` | municipality to mesoregion |
| `boundaries/sc_mesoregions.gpkg` | IBGE mesoregion polygons (without the coastal islands), used only to attribute each 2024 municipality to its mesoregion (`analysis/v10/meso2024/make_sc_mesoregions_2024.py`); `name_meso`, `forest_pct` (documented for three mesoregions only; not read by any script), `reflor_pct` |
| `boundaries/landcover2019_mesoregion.csv` | share of land-cover groups by mesoregion, MapBiomas Collection 10, 2019, from the biomass-burning pipeline (`tools/study_area_map.py` of that pipeline); fields `group`, `mesorregiao`, `km2`, `pct_of_group`, `pct_of_mesoregion`; input of `make_sc_mesoregions_2024.py`, which takes the `forest` rows as `forest_pct` of `sc_mesoregions_2024.gpkg` |
| `boundaries/oni.csv` | Oceanic Niño Index |
| `megan_inputs/EFP240120/` | `EFv240120.csv` (emission factors by ecotype and growth form) and `Speciation{Tree,Shrub,Herb,Crop}*.csv`, unchanged from `EFP240120.zip` of the MEGAN 3.2 input files (Guenther, 2024, Zenodo 10.5281/zenodo.10939297) |
| `sc_ecotype.nc` | MEGAN ecotype codes on the Santa Catarina tile, clipped from `EVT3d.nc` of the same record (variable `EVT`; 65535 = no data, filled by FIX-5 at run time) |
| `external/mapbiomas_c10_lulc2019_SC_120m.npz` | MapBiomas Collection 10, 2019, Santa Catarina window at about 120 m: `a` (class codes, uint8), `inside` (state mask), `transform` (affine of the window); study-area map |
| `external/mapbiomas_c10_lulc2019_SC_class_areas.csv` | `class_id`, `class_name`, `area_km2`, `pct_classified`: MapBiomas classes inside the state in 2019 |
| `external/ibge_br_uf_2024_simplified.geojson`, `external/ibge_america_do_sul_2021_simplified.geojson` | IBGE outlines (states 2024, South America 2021), simplified |
| `external/co2_annmean_gl.txt` | NOAA GML global annual mean CO₂ (`year mean unc`, ppm); read by `analysis/v10/review/final_review_diagnostics.py` |

## `run_v10/`: the production run (notebook 10, 2003-2025)

One folder per year with 90 files, plus five tables at the top.

| File | Format | Content |
|---|---|---|
| `<year>/emissions/emis_<year>_<mm>.nc` | NetCDF (CF-1.8), 19 float32 variables (lat, lon) | emission of each class in **tonnes per cell per month**; fill value −9999; global attributes `year`, `month`, `mechanism`, `megan_sc_fixes_sha256` (a933b29602929a82 in every file) and `variant` (`base`; `lai_viirs` in 2024 and 2025, see below) |
| `<year>/emissions/annual_emis_<year>.nc` | same | **tonnes per cell per year**, the sum of the 12 monthly files |
| `<year>/inputs/met_<year>_<mm>.tif` | GeoTIFF, 6 float32 bands | monthly mean over all 3-hourly steps of the forcing as sampled from ERA5-Land: `TEMP2` (K), `PAR` (W m⁻², 0.45 × short-wave, de-accumulated as in FIX-8), `WINDSPD` (m s⁻¹), `PRES` (Pa), `QV` (kg kg⁻¹), `SOILM` (m³ m⁻³); land cells only, NaN where ERA5-Land has no data (the cells filled by FIX-1 at run time) |
| `<year>/inputs/lai_<year>_<mm>.tif` | GeoTIFF, 1 band | monthly MODIS LAI (in 2024 and 2025 the VIIRS-corrected LAI of notebook 17) used as the current LAI, after the gap filling of FIX-4 (m² m⁻²) |
| `<year>/inputs/growthform_<vyear>.tif`, `.npz` | GeoTIFF with bands `tree`, `shrub`, `herb`, `crop`; NumPy archive with the same arrays and `veg` | fractions of each MEGAN growth form from MapBiomas (`vyear` = min(year, 2024)) |
| `<year>/inputs/ef_<year>.npz` | `ef` float32 (19, 144, 224), `classes` | emission factor of each class after weighting the ecotype factors by the growth forms (nmol m⁻² s⁻¹) |
| `<year>/inputs/fixlog_<year>_<mm>.npz` | NumPy archive | `filled_met_steps` (per cell: number of steps whose meteorology came from a neighbour, FIX-1), `lai_gap_cells` (cells filled by FIX-4 that month), `gamsm_mean` (mean soil-moisture factor γSM), `no_tonnes_mgn321_season` (soil NO in tonnes with the growing season as coded in MEGANv3.21, written by notebook 10 and not used by the analyses), `steps` (number of 3-hourly steps) |
| `<year>/hourly/diurnal_<year>_<mm>.csv` | CSV, classes × UTC hours | diurnal profile: mean over the days of the month of the class flux summed over all cells of the box (nmol m⁻² s⁻¹ summed over cells; a shape diagnostic, used by test S1) |
| `<year>/stats/stats_<year>_<mm>.csv` | CSV | state total (t), number of cells, maximum and mean per cell for each class and month |
| `<year>/stats/municipality_<year>_<mm>.csv` | CSV | tonnes by IBGE municipality code (`CD_MUN`), class and month (notebook 10 zonal sums; the published municipal numbers come from notebook 13) |
| `<year>/stats/annual_by_family_<year>.csv` | CSV | state total (t) of each class and its chemical family |
| `<year>/provenance_v10.json` | JSON | settings of the run of that year: sampling (3 h), LAI aggregation, wilting point (0.12), `fixes.py` hash, patches, ecotype file, engine folder, Earth Engine project (`"withheld"` in this release, for double-blind review; `tools/anonymize_release.py`), start time |
| `annual_totals_SC_v10_*.csv` | CSV | state totals by class and year, as printed by notebook 10 (one file per batch of years) and by notebook 11 (`all_classes`) |
| `verify_v10_report.csv` | CSV | the result of checks C0 to C5 of notebook 11 for each year |

## `tropomi_hcho_v10/`: TROPOMI formaldehyde on the model grid (notebook 14)

`hcho_<year>_<mm>.npz`, 84 files (2019-2025), one per month, written by cell 3 of notebook 14 with the screening
`cf0.3_sza70_fire120` (cloud fraction < 0.3, solar zenith angle < 70°, fire mask of 20 km). Arrays use the engine row
order (south first).

| Key | Shape | Content |
|---|---|---|
| `H` | (144, 224) float32 | mean tropospheric HCHO column of the month (molecules cm⁻²) |
| `N` | (144, 224) float32 | mean number of kept observations per native pixel |
| `bg` | scalar | Atlantic background: monthly median over 45° W to 41° W, 30° S to 25° S (molecules cm⁻²) |
| `removed`, `kept` | scalars | screened-out and kept observation counts summed on the grid (only their ratio is used) |
| `orbits` | scalar | number of TROPOMI L3 images of the month |

## `lai_viirs_v10/`: MODIS and VIIRS LAI on the model grid (notebook 16)

Written by notebook 16 in Colab on 31 October 2026 (Earth Engine project of the authors) and copied here unchanged.

| File | Content |
|---|---|
| `lai_viirs_modis_cache.npz` | keys `<sensor>_<year>_<mm>_LAI` and `<sensor>_<year>_<mm>_FRAC`, sensor `MODIS` or `VIIRS`, 2013-2025 (624 arrays, (144, 224) float32, south first): monthly LAI, area mean of the 500 m main-algorithm pixels of each cell (MODIS MCD15A3H C6.1 with `FparLai_QC` bit 0 = 0, as notebook 10; VIIRS VNP15A2H C2 with SCF_QC 0 or 1), and the share of retrievals made by the main algorithm |
| `lai_viirs_modis_monthly.csv`, `lai_viirs_modis_annual.csv`, `lai_viirs_modis_summary.csv`, `lai_viirs_modis.png` | the tables and figure of the run (also in `results/v10/lai_viirs_check`) |

The MODIS arrays equal the MODIS LAI inputs of the production run wherever both are positive (`run_v10` for
2013-2023, `run_v10_modis_lai_2024_2025` for 2024 and 2025) (maximum difference 5 × 10⁻⁷; the run also masks non-vegetated cells and fills 12 cell-months). **The VIIRS `LAI` arrays of this cache are LAI divided by 10**: the run that wrote it multiplied the Earth Engine band, which already holds LAI, by 0.1. Every use of them (deviations, z-scores, the ratios of `src/megan_sc/laicorr.py`) is a ratio and is unaffected. Notebook 16 as distributed reads the band without that factor. This folder is released as its own zip file (`lai_viirs_v10_2013-2025.zip`).

## `run_v10_modis_lai_2024_2025/` and `run_v10_modis_view/`: the MODIS-LAI years

`run_v10/2024` and `run_v10/2025` hold the years computed by notebook 17 with the VIIRS-corrected LAI (provenance variant
`lai_viirs`), and the same years computed with the MODIS LAI (variant `base`) are in `run_v10_modis_lai_2024_2025/` (same
layout, released as its own zip). `python tools/lai_variant_view.py` builds
`run_v10_modis_view/`, a folder of links holding the MODIS record of every year; the step `laicorr` of `reproduce.py`
builds and reads it. The view is never summed or packed.

## `s1_test/`: test S1 (time step)

Two short runs of notebook 10 for 2019 with `TEST_DAYS = 7` (first seven days of January, April, July and October),
identical except for the sampling: `megan_sc_outputs_v10/` (3-hourly, as production) and
`megan_sc_outputs_v10_hourly/` (hourly). Same layout as `run_v10/`, with the suffix `_test7d` in the file names.
Compared by `analysis/v10/s1_compare.py`.
