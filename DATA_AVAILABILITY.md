# Data availability and provenance

Every file that the Santa Catarina example depends on is either in this repository, in one of the data archives
attached to the release v1.0.0, or in a public archive read at run time; a run of another region reads the same public
archives. This page lists each item with its source, identifier,
licence and the way to obtain or regenerate it. `data/README.md` describes the content and format of every data
file, and chapter 2 of the tutorial shows how to open them.

## 1. In the repository (git)

| Path | Content | Source |
|---|---|---|
| `data/boundaries/sc_state_ibge.gpkg` | Santa Catarina state polygon: the reporting mask (a cell belongs to the state when its centre lies inside) | IBGE, public |
| `data/boundaries/SC_Municipios_2024.gpkg` | the 295 municipalities, 2024 territorial mesh, simplified copy (about 82 thousand vertices); municipal sums split each cell among the polygons in proportion to its area inside each one | IBGE, public |
| `data/boundaries/sc_mesoregions_2024.gpkg` | the six IBGE geographic mesoregions as unions of the 2024 municipalities, coastal islands included; `forest_pct` = MapBiomas Collection 10, 2024, classes 3, 5 and 49 (% of the classified area; 2019 until 16 November 2026); reporting units of the mesoregion results (`table3_mesoregions_soa.csv`, `fig5`, `fig6`) | built by `analysis/v10/meso2024/make_sc_mesoregions_2024.py` from the IBGE layers |
| `data/boundaries/sc_municipios_2024_mesorregiao.csv` | mesoregion of each municipality | same script |
| `data/boundaries/sc_mesoregions.gpkg` | IBGE mesoregion polygons (Santa Catarina and São Francisco do Sul islands left out), used only to attribute each 2024 municipality to its mesoregion in `analysis/v10/meso2024/make_sc_mesoregions_2024.py` | IBGE |
| `data/boundaries/landcover2024_mesoregion.csv` | share of land-cover groups by mesoregion, MapBiomas Collection 10, 2024; input of `make_sc_mesoregions_2024.py`, which writes the `forest_pct` attribute of `sc_mesoregions_2024.gpkg` | MapBiomas (Souza Junior et al., 2020), computed by `analysis/v10/maps/landcover_mesoregion.py` from `data/external/mapbiomas_c10_lulc2024_SC_120m.npz`; with the full-resolution IBGE 2024 municipal mesh instead of the simplified mesh of `data/boundaries`, the shares differ by at most 0.02 percentage point, and the same script reproduces the 2019 table exactly |
| `data/boundaries/landcover2019_mesoregion.csv` | the same for 2019 (sha256 1c4f8d02...); source of `forest_pct` until 16 November 2026, kept as the record of the earlier values | MapBiomas (Souza Junior et al., 2020), computed with the method of `analysis/v10/maps/landcover_mesoregion.py` (before that script was added) with the full-resolution IBGE 2024 municipal mesh |
| `data/boundaries/oni.csv` | Oceanic Niño Index, annual and DJF, 2003-2025 | NOAA CPC ONI v5, public |
| `data/megan_inputs/EFP240120/` | MEGAN 3.2 emission-factor table (`EFv240120.csv`) and the four speciation tables by growth form, unchanged | `EFP240120.zip` of the MEGAN 3.2 input files (Guenther, 2024a, Zenodo 10.5281/zenodo.10939297), CC BY 4.0 |
| `data/sc_ecotype.nc` | MEGAN ecotype tile of the Santa Catarina bounding box (variable `EVT`) | clipped from `EVT3d.nc` (`EVT3d.zip`) of the MEGAN 3.2 input files (Guenther, 2024a), CC BY 4.0 |
| `results/v10/` | every published table, number, figure and log of the production run | notebooks 11, 13-16 and 18-23 and `analysis/v10/` |
| `data/external/mapbiomas_c10_lulc2019_SC_120m.npz` | MapBiomas Collection 10 land use and land cover of 2019, Santa Catarina window read at about 120 m (4 × the native 30 m) from the internal overviews of the authors' copy of `brazil_coverage_2019.tif`; arrays `a` (class codes), `inside` (state mask) and `transform`; read by `analysis/v10/review/posthoc_review_diagnostics.py`, notebooks 20 and 23 and, until 16 November 2026, by the study-area map | MapBiomas (Souza Junior et al., 2020); extracted with the method of `analysis/v10/maps/extract_mapbiomas_sc.py` (before that script was added); to regenerate it, export 2019 with `analysis/v10/maps/mapbiomas_lulc_SC_export.ipynb` (`YEAR = 2019`) and run `extract_mapbiomas_sc.py <tif> 2019` |
| `data/external/mapbiomas_c10_lulc2019_SC_class_areas.csv` | area of each MapBiomas class inside the state (IBGE 2024 limit) in 2019, from the 60 m overview with pixel areas on the authalic sphere (km² and % of the classified area; class 0, coastal water inside the limit, excluded) | same |
| `data/external/mapbiomas_c10_lulc2024_SC_120m.npz` | MapBiomas Collection 10 land use and land cover of 2024 (the map also used for 2025), Santa Catarina window at about 120 m (nearest-neighbour reading at 1/4 of the native 30 m); same arrays; read by `analysis/v10/maps/study_area_map.py` (study-area figure) and `analysis/v10/maps/landcover_mesoregion.py` (land cover of the mesoregions) | MapBiomas (Souza Junior et al., 2020); exported on 16 November 2026 with `analysis/v10/maps/mapbiomas_lulc_SC_export.ipynb` and extracted with `analysis/v10/maps/extract_mapbiomas_sc.py`; aggregated as in the run it reproduces `growthform_2024.npz` of the run (tree fraction r = 0.99999) |
| `data/external/mapbiomas_c10_lulc2024_SC_class_areas.csv` | area of each MapBiomas class inside the state (IBGE 2024 municipal mesh) in 2024, at about 60 m with pixel areas on the authalic sphere (km² and % of the classified area; class 0 excluded) | same |
| `data/external/ibge_br_uf_2024_simplified.geojson`, `ibge_america_do_sul_2021_simplified.geojson` | outlines of the Brazilian states (IBGE BR_UF_2024, simplified at 0.01°) and of South America (IBGE Lim_america_do_sul_2021, simplified at 0.02°), for the study-area map | IBGE, public |
| `data/external/co2_annmean_gl.txt` | globally averaged marine surface annual mean CO₂, 1979-2025 (ppm), used for the effect of CO₂ inhibition on the isoprene trend | NOAA GML (Lan; Tans; Thoning, 2026), public domain; the official file downloaded from NOAA GML by notebook 19 on 1 November 2026 (file creation 9 October 2026); `co2_annmean_gl.SOURCE.txt` records the address, retrieval time, creation line and SHA-256 (9872655642690b06...) |
| `data/DATA_SHA256SUMS.txt` | SHA-256 of every file under `data/` (2,438 files: all except `README.md` and the two checksum lists), including the five unpacked archives below | `python tools/data_release.py sums` |

## 2. Data archives of the release (outside git)

Five zip files hold the data that is too large for git (0.90 GB to download; 2,416 files once unpacked). Notebook 19
(`PACK = True`) packed them. They will be attached to the GitHub release and deposited with the archived version (DOI in
`CITATION.cff` after deposit). `SHA256SUMS_assets.txt`, attached with them, gives the checksum of each zip; its copy is
`data/RELEASE_ASSETS_SHA256SUMS.txt`. Unpacked on a machine other than the one that packed them, the five zips pass
`python reproduce.py check-data`: the 2,438 data files listed in `data/DATA_SHA256SUMS.txt`, 0 missing or different.

| Archive (size, SHA-256 prefix) | Unpacks to | Content | Used by |
|---|---|---|---|
| `megan_sc_run_v10_2003-2025.zip` (723 MB, 206b0aa4) | `data/run_v10/` | the production run of notebook 10, 2003-2025, with 2024 and 2025 computed with the VIIRS-corrected LAI (notebook 17): monthly and annual emission grids (NetCDF), the inputs sampled on the model grid (monthly meteorology and LAI GeoTIFFs, growth forms, emission factors), the fix logs, the diurnal profiles, the monthly statistics and the provenance file of each year; the state totals and the report of notebook 11 | every step of `reproduce.py` |
| `tropomi_hcho_v10_2019-2025.zip` (18 MB, 32bfc1ed) | `data/tropomi_hcho_v10/` | 84 monthly TROPOMI HCHO fields on the model grid (mean column, observation count, Atlantic background, screened and kept counts), as cached by notebook 14 with the screening `cf0.3_sza70_fire120` | `reproduce.py nb14`, `nb14diag`, `maps`, `final` |
| `megan_sc_s1_test_2019.zip` (39 MB, 0e79ec0a) | `data/s1_test/` | the two 28-day test runs of notebook 10 for test S1 (first seven days of January, April, July and October 2019, 3-hourly and hourly sampling) | `reproduce.py s1` |
| `lai_viirs_v10_2013-2025.zip` (59 MB, 32ef9cee) | `data/lai_viirs_v10/` | the monthly VIIRS and MODIS LAI on the model grid, 2013-2025, cached by notebook 16 (`data/README.md` describes the VIIRS scale of this cache) | `reproduce.py laicorr`, notebook 17 |
| `megan_sc_run_v10_modis_lai_2024_2025.zip` (63 MB, d2e819c9) | `data/run_v10_modis_lai_2024_2025/` | 2024 and 2025 computed with the MODIS LAI (layout of `run_v10`) | `reproduce.py laicorr` (MODIS view), the MODIS control (REPRODUCE.md, section B) |

Download and unpack with `python tools/data_release.py fetch <base-url>`, then check with
`python reproduce.py check-data`. The base URL will be `https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0` for the
GitHub release or `https://zenodo.org/records/<record>/files` for the archived copy; a local folder holding the zips
also works.

The `provenance_v10.json` file of every model year (23 in `run_v10`, 2 in `run_v10_modis_lai_2024_2025`, 2 in `s1_test`)
records the Earth Engine project of the run in the field `ee_project`. For double-blind review, that value reads
`"withheld"` in the archives of this release (`tools/anonymize_release.py`, which copies every other member of the
archives unchanged, with the same CRC-32 and size); `data/DATA_SHA256SUMS.txt` lists those files as released.

## 3. Read at run time through Google Earth Engine (notebook 10, the region wizard and the satellite notebooks)

| Variable | Product | Earth Engine collection | Notes |
|---|---|---|---|
| 2 m temperature and dewpoint, surface pressure, 10 m wind, downward short-wave radiation, layer-1 soil water | ERA5-Land hourly (Muñoz-Sabater et al., 2021) | `ECMWF/ERA5_LAND/HOURLY` | land only, so coastal cells need FIX-1; radiation is accumulated from 00 UTC and the field valid at 00 UTC holds the previous 24 h (FIX-8); PAR = 0.45 × short-wave |
| Leaf area index | MODIS MCD15A3H version 6.1, 4-day, 500 m (Myneni; Knyazikhin; Park, 2021) | `MODIS/061/MCD15A3H` | good-quality retrievals only (FparLai_QC bit 0); monthly mean; the record starts in mid-2002, so 2003 is the first full year |
| Land cover, Brazil | MapBiomas Collection 10, 30 m, annual (Souza Junior et al., 2020) | `projects/mapbiomas-public/assets/brazil/lulc/v1` | the configuration fixes the collection (`collection: 10`) and `src/megan_sc/gee/mapbiomas.py` takes the newest release within it; 1985-2024, so 2025 reuses 2024 |
| Land cover, any country | MODIS MCD12Q1 version 6.1, 500 m, annual, legend LC_Type5 (Sulla-Menashe et al., 2019) | `MODIS/061/MCD12Q1` via `src/megan_sc/gee/modis_lc.py` | growth-form fractions when `growth_form_source: modis` (notebook 10, region wizard); 2001-2024 |
| Tropospheric HCHO column | Sentinel-5P TROPOMI, L3 OFFL (De Smedt et al., 2018) | `COPERNICUS/S5P/OFFL/L3_HCHO` | Earth Engine keeps pixels with qa_value ≥ 0.5; notebook 14 adds cloud fraction < 0.3 and solar zenith angle < 70° |
| Active fires, for the HCHO screening | MODIS active fires (FIRMS) | `FIRMS` | pixels within 20 km of a detection on the same or previous UTC day are removed |

A rerun reads these collections as they are on the day of the rerun. ERA5-Land and the MODIS collection 6.1 are
fixed archives; the MapBiomas asset receives new collections and new releases, so the configuration fixes the
collection, and notebook 10 reuses the cached growth forms when `growthform_<year>.npz` is present (tutorial chapter 6
shows how to seed them from `data/run_v10`).

## 4. Obtained separately

* **MEGANv3.21 code** (`MEGANv3.21.tar.gz`, Guenther, 2024b, Zenodo 10.5281/zenodo.10526206): needed only by the
  Fortran benchmark (notebook 12). The tables the engine reads are not part of it (section 1).
* **MEGAN global ecotype map** (`EVT3d.zip` of the MEGAN 3.2 input files, Guenther, 2024a, Zenodo
  10.5281/zenodo.10939297; 11 MB, unpacks to the 1.5 GB `EVT3d.nc`): needed only to run a region outside the Santa
  Catarina tile (region wizard, STEP 5).
* **NOAA GML global annual CO₂** (`co2_annmean_gl.txt`): shipped in `data/external/` (section 1); input of
  `analysis/v10/review/final_review_diagnostics.py` for the isoprene trend with CO₂ inhibition, which the production run
  leaves off. The published values with CO₂ inhibition (isoprene +0.65 % yr⁻¹, q = 0.03) depend on it; the section is
  skipped without the file.

## 5. Licences

IBGE boundaries (public), NOAA CPC ONI (public), ERA5-Land (Copernicus licence), MODIS products (NASA, public),
MapBiomas (licence stated by the MapBiomas project), TROPOMI (Copernicus Sentinel data), MEGAN 3.2 input files
(CC BY 4.0: the tables of `data/megan_inputs/` are redistributed unchanged and `data/sc_ecotype.nc` is a clip of
`EVT3d.nc`), VIIRS LAI products (NASA, public). The INMET station records and the BRAIN datasets
(Science Data Bank) are not redistributed; `results/v10/era5land_vs_inmet/`, `brain_comparison/` and
`brain_differences/` hold tables derived from them, and notebooks 21 and 18 download the originals from their sources. The model output in the release archives and the results in `results/` are released under
CC-BY-4.0 (`LICENSE-docs-data.md`); the input fields sampled from ERA5-Land, MODIS and MapBiomas keep the attribution
terms of those products.
