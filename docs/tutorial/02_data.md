# 2. The data

The inventory rests on five public inputs, the tables of the MEGAN code, and the IBGE boundaries. This chapter says
what each one is, how the model uses it, what is shipped with the repository, and how to open every file.
`DATA_AVAILABILITY.md` lists identifiers and licences, and `data/README.md` the format of every file.

## 2.1 Inputs

| Input | Product | Resolution | How the model uses it |
|---|---|---|---|
| Meteorology | ERA5-Land hourly reanalysis (Muñoz-Sabater et al., 2021), `ECMWF/ERA5_LAND/HOURLY` | 0.1° (about 9 km), hourly | 2 m temperature, short-wave radiation (PAR = 0.45 × short-wave), 10 m wind, surface pressure, humidity (from the dewpoint) and layer-1 soil water drive the canopy model, the activity factors and soil NO; sampled every 3 h |
| Leaf area index | MODIS MCD15A3H version 6.1 (Myneni; Knyazikhin; Park, 2021), `MODIS/061/MCD15A3H` | 500 m, 4 days | monthly mean of the good-quality retrievals in each 0.025° cell; LAI multiplies every class and sets the leaf age |
| Land cover | MapBiomas Collection 10 (Souza Junior et al., 2020) | 30 m, annual | growth-form fractions (tree, shrub, herb, crop) that weight the emission factors, the canopy types, and the vegetation mask |
| Emission factors and speciation | MEGAN 3.2 input files: tables `EFP240120` and ecotype map `EVT3d` (Guenther, 2024a) | ecotype tile at the resolution of the MEGAN map | emission factor of each class by ecotype and growth form; split of each class into 201 compounds |
| Boundaries | IBGE state, 2024 municipal mesh (IBGE, 2024) and the mesoregions built from it | vector | reporting mask (cell centre inside the state); municipal and mesoregion sums, with each cell split among the polygons by its area inside each one |
| HCHO columns (evaluation only) | Sentinel-5P TROPOMI L3 OFFL (De Smedt et al., 2018), `COPERNICUS/S5P/OFFL/L3_HCHO` | about 1.1 km in Earth Engine | monthly means on the model grid, compared with modelled isoprene (chapter 10) |

All inputs except the MEGAN tables are read through Google Earth Engine (Gorelick et al., 2017) by
`src/megan_sc/gee/` and notebooks 10 and 14.

Two facts about the inputs shape the method:

* **ERA5-Land covers land only.** Cells on the coast, on lagoons and on islands have no data in some or all bands. In
  the run, 190 vegetated cells of the state (1,301 km²; 191 cells in 2021 and 2022) take their meteorology from the
  nearest valid cell (FIX-1).
* **The radiation field is accumulated.** ERA5-Land stores short-wave radiation as an accumulation from 00 UTC, and the
  field valid at 00 UTC holds the previous 24 h (ECMWF, 2025). The hourly flux is the difference between consecutive
  hours, except at 01 UTC (FIX-8).

## 2.2 What is shipped

| Folder | In git? | Size | Content |
|---|---|---|---|
| `data/boundaries/` | yes | 3.9 MB | IBGE state, municipalities and mesoregions; MapBiomas 2019 land-cover shares by mesoregion (`landcover2019_mesoregion.csv`); ONI |
| `data/megan_inputs/`, `data/sc_ecotype.nc` | yes | 2.6 MB | MEGAN emission-factor and speciation tables, ecotype tile |
| `data/run_v10/` | release archive | 840 MB | the production run of notebook 10, 2003-2025, with 2024 and 2025 from notebook 17 (corrected LAI; 2,075 files) |
| `data/tropomi_hcho_v10/` | release archive | 18 MB | 84 monthly TROPOMI fields on the model grid |
| `data/s1_test/` | release archive | 44 MB | the two 28-day test runs of test S1 |
| `data/lai_viirs_v10/` | release archive | 59 MB | monthly VIIRS and MODIS LAI on the model grid, 2013-2025, cached by notebook 16 (5 files) |
| `data/run_v10_modis_lai_2024_2025/` | release archive | 73 MB | 2024 and 2025 computed with the MODIS LAI (layout of `run_v10`; 180 files) |

Chapter 3 shows how to download and check the release archives (five zip files, 0.90 GB). Sizes above are unpacked sizes. Each year of
`data/run_v10/` holds the monthly and annual emission grids, the inputs as the model saw them (monthly means of the
meteorology, the gap-filled LAI, the growth forms, the emission factors), a fix log per month, the diurnal profiles,
the monthly statistics and a provenance file. The ERA5-Land hourly fields themselves are not shipped: they would take
hundreds of gigabytes and can be read again from Earth Engine (chapter 6).

## 2.3 Orientation of the grids

The grid has 144 rows and 224 columns of 0.025°, with cell centres from 29.4875° S to 25.9125° S and from 53.8875° W to
48.3125° W. Two conventions coexist and are the most common source of mistakes:

* NetCDF files and NumPy arrays have row 0 at the **south** (latitude ascending), as the engine computes;
* GeoTIFF files have row 0 at the **north**, as GIS programs expect.

Flip a GeoTIFF array with `numpy.flipud` before combining it with a NetCDF field. GIS programs read both correctly,
because the NetCDF files carry their latitude coordinate.

## 2.4 Opening the files in Python

**State total of a class from an annual grid.** A cell belongs to the state when its centre lies inside the IBGE
polygon:

```python
import numpy as np, xarray as xr, geopandas as gpd, shapely

ds = xr.open_dataset('data/run_v10/2019/emissions/annual_emis_2019.nc')     # tonnes per cell per year
lon2, lat2 = np.meshgrid(ds.lon.values, ds.lat.values)
state = gpd.read_file('data/boundaries/sc_state_ibge.gpkg').to_crs(4326).geometry.union_all()
inside = shapely.contains_xy(state, lon2, lat2)

RHO_ISOP = 0.96                                           # canopy loss and production factor of isoprene (section 1.1)
print(int(inside.sum()))                                  # 13970 cells
isop = float(ds['ISOP'].where(inside).sum())
print(isop, RHO_ISOP * isop)                              # 262900.6 t as computed by the code; 252384.6 t with rho
bvoc = sum(ds[c] * (RHO_ISOP if c == 'ISOP' else 1.0) for c in ds.data_vars if c not in ('NO', 'CO'))
print(float(bvoc.where(inside).sum()) / 1e3)              # 739.4 kt of BVOC: the 2019 peak of article Figure 3 (2024: 739.3 kt)
```

The isoprene total without ρ equals the value written by notebook 10 in
`data/run_v10/2019/stats/annual_by_family_2019.csv` (262,900.6 t); with ρ it is the 2019 isoprene of Table S1 of the article (Table 2 of the
report; 252,385 t), and the BVOC total is the 2019 value of the same table. (`union_all` needs geopandas 1.0 or later; older
versions use `unary_union`.)

**Intensity per area.** Divide by the cell area. The analysis scripts use a sphere of radius 6,371,008.8 m:

```python
R, d = 6371008.8, np.deg2rad(0.025)
lat = ds.lat.values
area = (R**2 * d * np.abs(np.sin(np.deg2rad(lat + 0.0125)) - np.sin(np.deg2rad(lat - 0.0125))))[:, None] * np.ones((1, ds.lon.size))
print(area[inside].sum() / 1e6)                           # 95969 km² of grid area inside the state
kg_ha = bvoc * 1e3 / (area / 1e4)                         # kg per hectare per year
```

**Meteorology and LAI as the model saw them** (monthly means of all 3-hourly steps, land cells only):

```python
import rasterio
with rasterio.open('data/run_v10/2019/inputs/met_2019_01.tif') as src:
    print(src.descriptions)                               # ('TEMP2', 'PAR', 'WINDSPD', 'PRES', 'QV', 'SOILM')
    t2 = np.flipud(src.read(1))                           # K, now south-first like the NetCDF grids
print(np.nanmean(np.where(inside, t2, np.nan)) - 273.15)  # 23.6 °C, state mean of January 2019
```

**Fix logs** (what the procedures FIX-1 to FIX-8 did in a month):

```python
z = np.load('data/run_v10/2019/inputs/fixlog_2019_01.npz')
filled = z['filled_met_steps'] > 0                        # cells whose meteorology came from a neighbour (FIX-1)
print(int((filled & inside).sum()), int(z['lai_gap_cells']), int(z['steps']))
```

**TROPOMI fields** (chapter 10):

```python
z = np.load('data/tropomi_hcho_v10/hcho_2019_01.npz')
print(sorted(z.files))                                    # ['H', 'N', 'bg', 'kept', 'orbits', 'removed']
dhcho = z['H'] - z['bg']                                  # ΔHCHO, molecules cm-2, south-first
```

## 2.5 Opening the files in QGIS

Drag a `.nc` file into QGIS and choose the variable (one layer per class); drag a `.tif` file to see its bands; the
`.gpkg` files open as vector layers. Values are tonnes per cell per month (or year); to map intensities, use the
raster calculator with the cell area, or the maps of `results/v10/figures/`.

## 2.6 Units at a glance

| Quantity | Unit in the files | Conversion used in the documents |
|---|---|---|
| Emission per cell | t per month (`emis_*`) or t per year (`annual_emis_*`) | state totals in t yr⁻¹ or kt yr⁻¹ (10³ t), Gg (10³ t) in the comparison with the previous inventory |
| Intensity | not stored | kg ha⁻¹ yr⁻¹ = t × 1000 / (area in m² / 10⁴) |
| Isoprene areal flux | not stored | g C m⁻² yr⁻¹ = t × 10⁶ × (5 × 12.011 / 68.117) / area in m² |
| Emission factor | nmol m⁻² s⁻¹ (`ef_*.npz`) | |
| Temperature | K (`TEMP2`) | °C in the driver tables |
| PAR | W m⁻² | PPFD = 4.5 × PAR in µmol m⁻² s⁻¹ |
| HCHO column | molecules cm⁻² | |

## 2.7 Checksums

`data/DATA_SHA256SUMS.txt` holds the SHA-256 of every file under `data/`, and `python reproduce.py check-data` checks
them. Every emission file also records the hash prefix of `fixes.py` (`megan_sc_fixes_sha256` = a933b29602929a82),
which identifies the code of the procedures that produced it.
