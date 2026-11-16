# %% [markdown]
# # 17 · Rerun of 2024 and 2025 with the VIIRS-corrected LAI
#
# **What this notebook does.** It recomputes the inventory of 2024 and 2025 exactly as notebook 10 does, with the same
# engine, the same input and run procedures (`fixes.py`, written by cell 2 with the same SHA-256 prefix a933b29602929a82),
# the same three-hourly ERA5-Land sampling, the same MapBiomas growth forms and emission factors (reused from the
# production run) and the same output layout. One input changes: the monthly leaf area index (LAI) of the 24 months of
# 2024 and 2025, which is replaced by a corrected LAI built from the MODIS reference climatology and the VIIRS anomaly
# of each month (cell 7a). The outputs go to a separate folder, `megan_sc_outputs_v10_lai_viirs`, so the two versions
# can be compared (cell 11) before the corrected years replace the original ones (cell 12, optional and off by default).
#
# Code shared with notebook 10 is copied line for line; the changes are listed in the table at the end of this cell.
#
# ## 1. Why the two years are recomputed
#
# LAI enters MEGAN in three places (Guenther et al. 2006, 2012): it multiplies the emission activity
# (γ = LAI × γTP × γage × γSM × γCO2), it sets how much of the foliage is sunlit or shaded in the five-layer canopy
# model, and its change from the previous month sets the leaf-age factor γage. It also enters the canopy reduction of
# soil NO and the bidirectional factor of ethanol.
#
# The production run (notebook 10) uses the monthly MODIS LAI of Terra and Aqua (MCD15A3H Collection 6.1; Myneni et al.
# 2021), main radiative-transfer algorithm only (Knyazikhin et al. 1998), averaged over the 500 m pixels of each 0.025°
# cell. Its statewide value was 5.8% below the 2013-2023 mean in 2024 and 12.1% below in 2025, more than any earlier
# year (lowest earlier year: −2.0%). In the cells dominated by trees, herbaceous vegetation or crops the 2025 deficits
# were 11.4%, 12.0% and 8.7%.
#
# Notebook 16 compared it, on the same grid and with the same aggregation, with the LAI of the VIIRS instrument on Suomi
# NPP (VNP15A2H Collection 2; Myneni 2023), produced with the MODIS radiative-transfer algorithm adapted to VIIRS to
# continue the MODIS record (Yan et al. 2018; Román et al. 2024):
#
# * in 2013-2023 the two sensors agreed: annual deviations r = 0.84 to 0.91 across the four cell groups, deseasonalised
#   monthly statewide anomalies r = 0.79 (n = 132), annual MODIS-minus-VIIRS difference with a standard deviation of
#   0.75 percentage points; Xu et al. (2018) found the same agreement globally for 2012-2016;
# * in 2024 and 2025 the difference was −4.7 and −10.2 points, and MODIS fell below VIIRS by more than two standard
#   deviations of the reference difference in 11 of the 13 months from December 2024 (4 of 132 reference months);
# * the VIIRS deviations of 2025 (−1.8% statewide; trees −2.3%, herbaceous −0.2%, crops +1.9%) lie inside the
#   2013-2023 range of every group;
# * the share of MODIS retrievals made by the main algorithm in 2024-2025 was within its 2013-2023 range.
#
# Terra's orbit has drifted toward earlier overpass times since its final inclination adjustment of 27 February 2020,
# and that drift has biased at least one other MODIS land-product time series (MOD10A1 albedo; Feng et al. 2024). NASA
# reports the drift of both Terra and Aqua and advises users to move to VIIRS (NASA MODIS Land Team 2025; NASA Earthdata
# 2025). The decrease is therefore treated as a defect of the recent MODIS record. The cause is not established by
# these data; the correction rests on the comparison with VIIRS, not on a model of the drift.
#
# ## 2. The correction (`src/megan_sc/laicorr.py`, written into your `megan-sc` folder by cell 2b)
#
# For cell c of month m of year Y (2024 or 2025):
#
#     LAI*(c, Y, m) = M(c, m) × V_z(Y, m) / V_z(m)
#
# * **M(c, m)**: mean MODIS LAI of cell c in calendar month m over 2013-2023, read from the LAI files the production run
#   used (`<year>/inputs/lai_<year>_<mm>.tif`), over the years in which the cell has a value (at least 6 of 11).
# * **V_z(Y, m) / V_z(m)**: sum of the VIIRS LAI of month m of year Y over the cells of zone z, divided by the sum of
#   their 2013-2023 VIIRS means of month m, over the cells valid in both (VIIRS values from the notebook 16 cache: same
#   grid, area mean of the 500 m main-algorithm pixels). The quotient of two VIIRS values does not depend on the scale
#   of the VIIRS band.
# * **Zones z**: cells dominated by trees, by herbaceous vegetation or by crops (MapBiomas growth-form fraction above
#   0.6, growth forms of the run) and the other cells of Santa Catarina; zones with fewer than 25 valid cells, and cells
#   outside the state, take the statewide ratio.
# * Cells without M(c, m) keep their MODIS value multiplied by the statewide correction factor of the month (no such
#   cell in Santa Catarina in 2024-2025).
#
# The corrected field keeps the level, seasonal cycle and spatial pattern of the MODIS record of 2003-2023 and takes the
# year-to-year change of each zone from VIIRS.
#
# **How this form was chosen.** Each year of 2013-2023 was predicted from the other ten years with eleven candidate forms
# and compared with the LAI the run actually used (`analysis/v10/review/lai_viirs_correction.py`, outputs in
# `results/v10/lai_viirs_check`). The rule, fixed before the comparison, was the lowest root-mean-square error of the
# cell LAI over all Santa Catarina cell-months. Seven candidates were defined first (MODIS mean alone, one statewide
# ratio, one ratio per cell with and without bounds, 3 × 3 window with and without bounds, 5 × 5 window); because the
# per-cell and 3 × 3 ratios did worse than the MODIS mean alone, the list was extended once with intermediate scales
# (11 × 11 and 21 × 21 windows, mesoregions, growth-form zones). Result: growth-form zones 0.362 m² m⁻², statewide ratio
# 0.364, 21 × 21 window 0.364, mesoregions 0.365, MODIS mean alone 0.386, one ratio per cell 0.439. With the selected
# form, the statewide annual LAI of each hindcast year was predicted within 0.6% on average and 1.9% at most.
#
# **Effect on the LAI.** The corrected statewide annual LAI is 3.29 in 2024 (+5.0% over the MODIS value used, 3.14) and
# 3.27 in 2025 (+11.6% over 2.93), that is, 1.1% and 1.8% below the 2013-2023 mean. December 2023, the previous month of
# January 2024 for the leaf-age factor, keeps its MODIS value, as do all months before 2024. January 2026, used only to
# fill gaps of December 2025, is sampled from MODIS as in notebook 10.
#
# ## 3. How to run
#
# 1. Run notebook 16 first (it writes the VIIRS cache this notebook reads) and keep the production run of notebook 10 in
#    `ROOT_V10` (cell 1).
# 2. Set your Earth Engine project (cell 3). Earth Engine is needed for ERA5-Land, as in notebook 10; the MapBiomas
#    growth forms are copied from `ROOT_V10`, so they are not requested again.
# 3. Run all cells. Each year takes about as long as one year of notebook 10 (several hours in Colab). Finished months
#    are skipped and interrupted months resume from their checkpoint, so after a disconnection run all cells again.
# 4. Read cell 11: only the LAI must differ from the production run (the meteorology check must print 0), and the
#    emission changes are printed by class.
# 5. To make the corrected years part of the run, set `INSTALL = True` in cell 12 and run it: the original 2024 and 2025
#    folders are moved to `megan_sc_outputs_v10_modis_lai_2024_2025` and the corrected ones are copied into `ROOT_V10`.
#    Then run notebooks 11, 13 and 14 again (notebook 11 accepts the variant `lai_viirs` for 2024 and 2025 only).
#
# ## 4. Changes with respect to notebook 10
#
# | Cell | Change |
# |---|---|
# | 1 | `ROOT_V10` and `VIIRS_CACHE`; `OLD_ROOT = ROOT_V10`, so the growth forms of the production run are reused |
# | 2b (new) | writes `megan_sc/laicorr.py` (the correction) and prints its SHA-256 prefix |
# | 4 | `YEARS = [2024, 2025]`, `VARIANT = 'lai_viirs'` (outputs in `megan_sc_outputs_v10_lai_viirs`) |
# | 7 | `lai_raw` returns the corrected field for the months of cell 7a, and the MODIS sample otherwise |
# | 7a (new) | corrected LAI of the 24 months; check against the published fields |
# | 8 | provenance and NetCDF attributes record the LAI source and the `laicorr.py` hash |
# | 11, 12 (new) | comparison with the production run; optional installation |
#
# ## References
#
# * Feng S, Wehrlé A, Cook JM, Anesio AM, Box JE, Benning LG, Tranter M (2024) The apparent effect of orbital drift on
#   time series of MODIS MOD10A1 albedo on the Greenland ice sheet. Science of Remote Sensing 9:100116.
#   https://doi.org/10.1016/j.srs.2023.100116
# * Guenther A, Karl T, Harley P, Wiedinmyer C, Palmer PI, Geron C (2006) Estimates of global terrestrial isoprene
#   emissions using MEGAN (Model of Emissions of Gases and Aerosols from Nature). Atmospheric Chemistry and Physics
#   6(11):3181-3210. https://doi.org/10.5194/acp-6-3181-2006
# * Guenther AB et al (2012) The Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2.1): an extended
#   and updated framework for modeling biogenic emissions. Geoscientific Model Development 5(6):1471-1492.
#   https://doi.org/10.5194/gmd-5-1471-2012
# * Knyazikhin Y, Martonchik JV, Myneni RB, Diner DJ, Running SW (1998) Synergistic algorithm for estimating vegetation
#   canopy leaf area index and fraction of absorbed photosynthetically active radiation from MODIS and MISR data. Journal
#   of Geophysical Research: Atmospheres 103(D24):32257-32275. https://doi.org/10.1029/98JD02462
# * Myneni R (2023) VIIRS/NPP Leaf Area Index/FPAR 8-Day L4 Global 500m SIN Grid V002. NASA EOSDIS Land Processes
#   Distributed Active Archive Center. https://doi.org/10.5067/VIIRS/VNP15A2H.002
# * Myneni R, Knyazikhin Y, Park T (2021) MODIS/Terra+Aqua Leaf Area Index/FPAR 4-Day L4 Global 500m SIN Grid V061.
#   NASA EOSDIS Land Processes Distributed Active Archive Center. https://doi.org/10.5067/MODIS/MCD15A3H.061
# * NASA Earthdata (2025) VIIRS instruments become more essential as Terra and Aqua drift from their traditional
#   orbits. https://www.earthdata.nasa.gov/news/feature-articles/viirs-instruments-become-more-essential-terra-aqua-drift-traditional-orbits
# * NASA MODIS Land Team (2025) MODIS Land news. https://modis-land.gsfc.nasa.gov/news.html
# * Román MO et al (2024) Continuity between NASA MODIS Collection 6.1 and VIIRS Collection 2 land products. Remote
#   Sensing of Environment 302:113963. https://doi.org/10.1016/j.rse.2023.113963
# * Xu B, Park T, Yan K, Chen C, Zeng Y, Song W, Yin G, Li J, Liu Q, Myneni RB (2018) Analysis of global LAI/FPAR products
#   from VIIRS and MODIS sensors for spatio-temporal consistency and uncertainty from 2012-2016. Forests 9(2):73.
#   https://doi.org/10.3390/f9020073
# * Yan K, Park T, Chen C, Xu B, Song W, Yang B, Zeng Y, Liu Z, Yan G, Knyazikhin Y, Myneni RB (2018) Generating global
#   products of LAI and FPAR from SNPP-VIIRS data: theoretical background and implementation. IEEE Transactions on
#   Geoscience and Remote Sensing 56(4):2119-2137. https://doi.org/10.1109/TGRS.2017.2775247

# %%
# 1 · Setup ---------------------------------------------------------------------------------------------
__FILE:notebooks/src/_setup_repo.py__
%pip -q install earthengine-api numpy scipy xarray netCDF4 matplotlib pyyaml pandas rasterio geopandas
from google.colab import drive; drive.mount('/content/drive')
MEGAN_DIR = REPO_DIR                      # the repository folder (found by the setup block)
OUT_BASE  = '/content/drive/MyDrive'                   # <-- the outputs go to OUT_BASE/megan_sc_outputs_v10_lai_viirs
ROOT_V10  = '/content/drive/MyDrive/megan_sc_outputs_v10'   # <-- the production run of notebook 10 (LAI 2013-2025, growth forms)
VIIRS_CACHE = ROOT_V10 + '/analysis_v10/lai_viirs_check/lai_viirs_modis_cache.npz'   # <-- monthly VIIRS LAI, notebook 16
OLD_ROOT  = ROOT_V10                                   # growth forms of the production run are reused (same inputs)
import sys, os, json, time, datetime, calendar, logging, hashlib, glob, re, shutil
import numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, MEGAN_DIR + '/src')
for _lg in ('google_auth_httplib2', 'googleapiclient.http', 'google.auth'):
    logging.getLogger(_lg).setLevel(logging.ERROR)
print('engine:', MEGAN_DIR)

# %%
# 2 · Write megan_sc/fixes.py (the input and run procedures; src/megan_sc/fixes.py of the repository, 12 unit tests) --
FIXES_PY = r'''__FIXES_PY__'''
FIXES_SHA = hashlib.sha256(FIXES_PY.encode('utf-8')).hexdigest()[:16]
with open(os.path.join(MEGAN_DIR, 'src', 'megan_sc', 'fixes.py'), 'w', encoding='utf-8') as _f:
    _f.write(FIXES_PY)
import importlib
for _m in [m for m in list(sys.modules) if m.startswith('megan_sc')]:
    del sys.modules[_m]                                  # fresh import of the whole package
from megan_sc import fixes
print('fixes.py written | sha256[:16] =', FIXES_SHA)

# %%
# 2b · Write megan_sc/laicorr.py (the VIIRS-anomaly correction of the LAI; src/megan_sc/laicorr.py of the repository) --
LAICORR_PY = r'''__FILE:src/megan_sc/laicorr.py__'''
LAICORR_SHA = hashlib.sha256(LAICORR_PY.encode('utf-8')).hexdigest()[:16]
with open(os.path.join(MEGAN_DIR, 'src', 'megan_sc', 'laicorr.py'), 'w', encoding='utf-8') as _f:
    _f.write(LAICORR_PY)
from megan_sc import laicorr
print('laicorr.py written | sha256[:16] =', LAICORR_SHA, '| settings', laicorr.DEFAULT, '| reference', laicorr.REF_YEARS[0], '-', laicorr.REF_YEARS[-1])

# %%
# 3 · Earth Engine ---------------------------------------------------------------------------------------
EE_PROJECT = 'your-ee-project'         # <-- your Earth Engine (Google Cloud) project id
import ee
ee.Authenticate(); ee.Initialize(project=EE_PROJECT)
print('Earth Engine ready:', EE_PROJECT)

# %%
# 4 · Parameters -----------------------------------------------------------------------------------------
YEARS        = [2024, 2025]             # the two years whose LAI is corrected (do not change; see the header)
MONTHS       = list(range(1, 13))
SAMPLE_HOURS = 3                        # 8 steps/day, production setting (VARIANT 'hourly' sets 1)
TEST_DAYS    = None                     # e.g. 2 for a quick preview (outputs get a _testNd suffix, not used downstream)
VARIANT      = 'lai_viirs'              # production settings with the VIIRS-corrected LAI of 2024-2025 (this notebook)
LAI_AGG      = 'mean'                   # 'mean': area mean of QC-valid 500 m MODIS pixels per 0.025° cell (explicit,
                                        #  reduceResolution) | 'ee_default': Earth Engine default resampling (test S2)
MAKE_MAPS    = False                    # monthly/annual PNG maps (slow on Drive); the analysis notebooks make the figures
DIAG         = True                     # cheap diagnostics in the fix log (soil-NO season alternative, mean GAMSM)
MECH         = 'CB6X'
if VARIANT == 'hourly':
    SAMPLE_HOURS = 1
if VARIANT == 'lai_eedefault':
    LAI_AGG = 'ee_default'
assert VARIANT == 'lai_viirs' and SAMPLE_HOURS == 3 and LAI_AGG == 'mean' and YEARS == [2024, 2025], VARIANT
LAI_SOURCE = 'MODIS MCD15A3H C6.1 2013-2023 monthly mean x VIIRS VNP15A2H C2 anomaly ratio by growth-form zone (laicorr.py)'
OUT_ROOT = os.path.join(OUT_BASE, 'megan_sc_outputs_v10' + ('' if VARIANT == 'base' else f'_{VARIANT}'))
os.makedirs(OUT_ROOT, exist_ok=True)
print(f'years {YEARS[0]}-{YEARS[-1]} ({len(YEARS)}) | {24 // SAMPLE_HOURS} steps/day | LAI_AGG={LAI_AGG} | -> {OUT_ROOT}')

# %%
# 5 · Engine, configuration, patches and FIX-5 ------------------------------------------------------------
from megan_sc import domain as dom, run, landcover, stats as st, fixes
from megan_sc.io import geotiff, netcdf as ncio
from megan_sc.engine import speciation as spc, species as sp, soil
from megan_sc.gee import era5, lai as laimod, mapbiomas, region_from_bounds
PATCHES = fixes.apply_engine_patches()                  # FIX-6, FIX-7, FIX-8
assert fixes.patches_active(), 'engine patches not active'
_cfgp = MEGAN_DIR + '/config/example_config.yaml'
cfg = run.load_config(_cfgp if os.path.exists(_cfgp) else MEGAN_DIR + '/example_config.yaml')
cfg['domain']['mode'] = 'bbox'
for k in ('ef_table', 'speciation_dir', 'ecotype'):
    v = cfg['datasets']['emission_factors'].get(k)
    if v and not os.path.isabs(v):
        cfg['datasets']['emission_factors'][k] = os.path.join(MEGAN_DIR, v)
_eco = cfg['datasets']['emission_factors']['ecotype']
_eco_f = os.path.join(os.path.dirname(_eco), 'sc_ecotype_filled.nc')            # FIX-5
if not os.path.exists(_eco_f):
    print('FIX-5: ecotype no-data pixels filled:', fixes.fill_ecotype_nearest(_eco, _eco_f))
cfg['datasets']['emission_factors']['ecotype'] = _eco_f
domain = dom.build_domain(cfg); region = region_from_bounds(domain['bounds'])
res = cfg['domain']['resolution_deg']; veg = cfg['datasets']['vegetation']
lon, lat, area = domain['lon'], domain['lat'], domain['area']
WILT = float(cfg['datasets'].get('soil', {}).get('constant_wilt', 0.12))
HOURS = list(range(0, 24, SAMPLE_HOURS)); SEC_PER_STEP = SAMPLE_HOURS * 3600.0
CLS = list(sp.NAMES); MET_BANDS = ['TEMP2', 'PAR', 'WINDSPD', 'PRES', 'QV', 'SOILM']
# soil-NO compound(s) for the DIAG alternative (EF_NO x GAMNO x speciation fraction)
_NO_MAP = [(int(spc.MAP_SPCA[s]), float(spc.MAP_EFFS[s])) for s in range(spc.MAP_SPCA.size) if spc.MAP_MG20[s] == sp.IDX_NO]
print('patches:', PATCHES, '| grid', lat.size, 'x', lon.size, '| wilt', WILT, '| NO compounds', _NO_MAP)

# %%
# 6 · Reporting setup: SC mask, palettes, map helper, municipality zones -------------------------------------
import geopandas as gpd
from rasterio.features import rasterize as _rasterize
from mpl_toolkits.axes_grid1 import make_axes_locatable
DPI = 300
ASPECT = 1.0 / np.cos(np.deg2rad(float(np.mean(lat))))
SC_GPKG = MEGAN_DIR + '/data/boundaries/sc_state_ibge.gpkg'
SHP = MEGAN_DIR + '/data/boundaries/SC_Municipios_2024.shp'
if not os.path.exists(SHP):                      # the repository ships the 2024 mesh as a (simplified) GeoPackage
    SHP = MEGAN_DIR + '/data/boundaries/SC_Municipios_2024.gpkg'
def _load_sc():
    for p in (SC_GPKG, SHP):
        if os.path.exists(p):
            gg = gpd.read_file(p).to_crs('EPSG:4326')
            return gg.geometry.union_all() if hasattr(gg.geometry, 'union_all') else gg.geometry.unary_union
    return None
_BND = _load_sc()
assert _BND is not None, 'SC boundary not found in megan-sc/data/boundaries'
_tr = geotiff._transform(lon, lat)
SC_MASK = np.flipud(_rasterize([(_BND, 1)], out_shape=(lat.size, lon.size), transform=_tr, fill=0, dtype='uint8')).astype(bool)
domain['mask'] = SC_MASK
def sc_clip(field): return np.where(SC_MASK, field, np.nan)
def draw_boundary(ax, lw=0.5):
    for gm in getattr(_BND, 'geoms', [_BND]):
        try: xs, ys = gm.exterior.xy; ax.plot(xs, ys, color='0.25', lw=lw)
        except Exception: pass
FAMILIES = [('Isoprene-type', ['ISOP', 'MBO'], 'Greens'), ('Monoterpenes', ['MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM'], 'Blues'),
            ('Sesquiterpenes', ['SQT_HR', 'SQT_LR'], 'Purples'), ('Oxygenated VOC', ['MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD'], 'Oranges'),
            ('Light hydrocarbons', ['LVOC'], 'YlGn'), ('Other / stress VOC', ['STRESS', 'OTHER'], 'YlOrBr'),
            ('Soil NO', ['NO'], 'Reds'), ('CO', ['CO'], 'Greys')]
FAM_OF = {s: (nm, cm) for nm, ss, cm in FAMILIES for s in ss}
def save_map(field, title, path_png, cmap='viridis', vmin=None, vmax=None):
    if not MAKE_MAPS: return
    cm_ = plt.get_cmap(cmap).copy(); cm_.set_bad('white')
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    im = ax.imshow(np.ma.masked_invalid(sc_clip(field)), origin='lower', extent=[lon[0], lon[-1], lat[0], lat[-1]],
                   cmap=cm_, vmin=vmin, vmax=vmax, aspect=ASPECT)
    draw_boundary(ax); ax.set_title(title); ax.axis('off')
    fig.colorbar(im, cax=make_axes_locatable(ax).append_axes('right', size='4.5%', pad=0.05))
    fig.tight_layout(); fig.savefig(path_png, dpi=DPI, bbox_inches='tight'); plt.close(fig)
_ZONE = _IDS = None
try:
    _g = gpd.read_file(SHP).to_crs('EPSG:4326'); _idc = 'CD_MUN' if 'CD_MUN' in _g.columns else _g.columns[0]
    _ZONE, _IDS = st.rasterize_municipalities(_g, lon, lat, _idc)
    print(f'municipalities: {len(_IDS)} zones')
except Exception as e:
    print('no municipality shapefile (per-municipality CSVs skipped):', e)
print(f'SC mask: {int(SC_MASK.sum())} of {SC_MASK.size} cells')

# %%
# 7 · Earth Engine samplers: meteorology, 24-h mean temperature, LAI (FIX-4 gap fill; explicit aggregation) --
def sample_multi(img, bands):
    grid = {'dimensions': {'width': int(lon.size), 'height': int(lat.size)},
            'affineTransform': {'scaleX': res, 'shearX': 0.0, 'translateX': float(lon[0] - res / 2),
                                'shearY': 0.0, 'scaleY': -res, 'translateY': float(lat[-1] + res / 2)},
            'crsCode': 'EPSG:4326'}
    for attempt in range(4):                              # retry transient EE errors
        try:
            st_ = ee.data.computePixels({'expression': img.select(bands), 'grid': grid, 'fileFormat': 'NUMPY_NDARRAY'})
            return {b: np.flipud(np.asarray(st_[b], dtype='float64')) for b in bands}
        except Exception as e:
            if attempt == 3: raise
            print('   EE retry', attempt + 1, repr(e)[:90]); time.sleep(10 * (attempt + 1))

def daily_tmean_image(day_iso):
    """FIX-2: mean of the 24 hourly ERA5-Land 2 m temperatures of the UTC day (daymet.f: D_TEMP)."""
    d = ee.Date(day_iso)
    return (ee.ImageCollection(era5.COLLECTION).filterBounds(region).filterDate(d, d.advance(1, 'day'))
            .select('temperature_2m').mean().rename('D_TEMP'))

def lai_image(y, m):
    img = laimod.monthly_lai(ee, region, y, m).rename('LAI')
    if LAI_AGG == 'mean':                                 # area mean of the QC-valid 500 m pixels in each grid cell
        proj = ee.ImageCollection(laimod.COLLECTION).filterBounds(region).first().select('Lai').projection()
        img = img.setDefaultProjection(proj).reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)
    return img

_lai_raw, _lai_cache, LAI_GAP_LOG = {}, {}, {}
def _shift(y, m, k):
    t = y * 12 + (m - 1) + k
    return t // 12, t % 12 + 1
def lai_raw(y, m):
    if (y, m) not in _lai_raw:
        if (y, m) in LAI_CORR:                                # 2024-2025: corrected LAI of cell 7a
            a = np.array(LAI_CORR[(y, m)], dtype='float64')
        else:                                                 # other months (December 2023, January 2026): MODIS, as notebook 10
            a = np.nan_to_num(sample_multi(lai_image(y, m), ['LAI'])['LAI'], nan=0.0, posinf=0.0, neginf=0.0)
        if veg.get('lai_mask', True):
            a = a * (VEG_FRAC > 0)
        _lai_raw[(y, m)] = a
    return _lai_raw[(y, m)]
def month_lai(m, y=None):
    """FIX-4: gap-filled monthly LAI (used as LAIc, and for month-1 as LAIp)."""
    y = YEAR if y is None else y
    if (y, m) not in _lai_cache:
        yp, mp = _shift(y, m, -1); yn, mn = _shift(y, m, +1)
        filled, gap = fixes.fill_lai_gaps(lai_raw(y, m), lai_raw(yp, mp), lai_raw(yn, mn), VEG_FRAC)
        _lai_cache[(y, m)] = filled; LAI_GAP_LOG[(y, m)] = int(gap.sum())
    return _lai_cache[(y, m)]

def no_tonnes(gamno, ef_no):
    """Soil-NO tonnes for one step from a GAMNO field (same unit chain as speciation.py / mgn2mech.F)."""
    t = np.zeros_like(gamno)
    for ci, f in _NO_MAP:
        t += gamno * ef_no * f * spc.NANO2MOL * area * spc.SPCA_MWT[ci] * spc.G2TONNE * SEC_PER_STEP
    return t

# %%
# 7a · Corrected LAI of the 24 months of 2024 and 2025 (laicorr.corrected_run_lai) --------------------------------------
assert os.path.exists(VIIRS_CACHE), 'run notebook 16 first: ' + VIIRS_CACHE
LAI_CORR, LAI_CORR_LOG = laicorr.corrected_run_lai(ROOT_V10, VIIRS_CACHE, SC_MASK, years=YEARS)
_pub = os.path.join(MEGAN_DIR, 'results', 'v10', 'lai_viirs_check', 'lai_viirs_corrected_2024_2025.npz')
if os.path.exists(_pub):                                  # the fields published with the repository (same function)
    _P = np.load(_pub)
    _d = max(float(np.abs(_P[f'LAI_{y}_{m:02d}'] - LAI_CORR[(y, m)]).max()) for (y, m) in LAI_CORR)
    print(f'published corrected fields: max abs difference {_d:.2e}'); assert _d < 1e-4, _d
_rows = []
for (y, m), a in sorted(LAI_CORR.items()):
    o = laicorr.read_south_first(os.path.join(ROOT_V10, str(y), 'inputs', f'lai_{y}_{m:02d}.tif'))
    sel = SC_MASK & (o > 0) & (a > 0)
    inf = [i for i in LAI_CORR_LOG if i['year'] == y and i['month'] == m][0]
    _rows.append(dict(year=y, month=m, modis_used=round(float(o[sel].mean()), 3), corrected=round(float(a[sel].mean()), 3),
                      change_pct=round(100 * (a[sel].mean() / o[sel].mean() - 1), 1), statewide_viirs_ratio=round(inf['state_ratio'], 3),
                      zones_on_state_ratio=inf['ratio_fallback_cells'], cells_without_modis_mean=inf['modis_fallback_cells']))
LAI_CORR_TABLE = pd.DataFrame(_rows); print(LAI_CORR_TABLE.to_string(index=False))
LAI_CORR_INFO = dict(source=LAI_SOURCE, laicorr_sha256=LAICORR_SHA, settings=str(laicorr.DEFAULT), reference_years=list(laicorr.REF_YEARS),
                     viirs_cache=VIIRS_CACHE, modis_reference_root=ROOT_V10)

# %%
# 8 · One month and one year ------------------------------------------------------------------------------
def _finalize_month(m, tag, tonnes, meanflux, meanmet, diurnal):
    ncio.write_emissions(os.path.join(OUT, 'emissions', f'emis_{tag}.nc'), tonnes, CLS, lon, lat,
                         units='tonnes/month', global_attrs={'year': YEAR, 'month': m, 'mechanism': MECH,
                                                             'megan_sc_fixes_sha256': FIXES_SHA, 'variant': VARIANT,
                                                             'lai_source': LAI_SOURCE, 'megan_sc_laicorr_sha256': LAICORR_SHA})
    geotiff.write_species_stack(os.path.join(OUT, 'emissions', f'emis_{tag}.tif'), tonnes, lon, lat, CLS)
    for nm in ['ISOP', 'MT_PINE', 'SQT_HR', 'NO', 'CO']:
        save_map(tonnes[sp.index_of(nm)], f'{nm} {tag} (t/month, SC)', os.path.join(OUT, 'maps', f'{nm}_{tag}.png'),
                 cmap=FAM_OF.get(nm, (None, 'YlOrRd'))[1], vmin=0)
    geotiff.write_geotiff(os.path.join(OUT, 'inputs', f'met_{tag}.tif'), meanmet, lon, lat, band_names=MET_BANDS, nodata=np.nan)
    geotiff.write_geotiff(os.path.join(OUT, 'inputs', f'lai_{tag}.tif'), month_lai(m), lon, lat, band_names=['LAI'])
    pd.DataFrame(st.species_stats_rows(tonnes, CLS, YEAR, m, mask=domain.get('mask'))).to_csv(
        os.path.join(OUT, 'stats', f'stats_{tag}.csv'), index=False)
    if _ZONE is not None:
        pd.DataFrame(st.municipality_table(tonnes, CLS, _ZONE, _IDS, YEAR, m)).to_csv(
            os.path.join(OUT, 'stats', f'municipality_{tag}.csv'), index=False)
    pd.DataFrame(diurnal, index=CLS, columns=[f'{h:02d}h' for h in HOURS]).to_csv(os.path.join(OUT, 'hourly', f'diurnal_{tag}.csv'))

def run_month(m):
    ndays_full = calendar.monthrange(YEAR, m)[1]
    ndays = ndays_full if TEST_DAYS is None else min(int(TEST_DAYS), ndays_full)
    tag = f'{YEAR}_{m:02d}' + ('' if TEST_DAYS is None else f'_test{ndays}d')
    if os.path.exists(os.path.join(OUT, 'emissions', f'emis_{tag}.nc')):
        print(f'  {tag}: already done, skipping'); return
    ck = os.path.join(OUT, 'checkpoints', f'{tag}.npz')
    steps = [(d, h) for d in range(1, ndays + 1) for h in HOURS]
    ny, nx = lat.size, lon.size
    keys = ('t', 'f', 'met', 'metn', 'fill', 'di', 'gsm', 'noalt')
    loaded = st.load_checkpoint(ck)
    if loaded is not None:
        A_, done = loaded
        acc = {k: A_[k] for k in keys}
        print(f'  {tag}: resuming ({len(done)}/{len(steps)} steps done)')
    else:
        acc = dict(t=np.zeros((len(CLS), ny, nx)), f=np.zeros((len(CLS), ny, nx)), met=np.zeros((len(MET_BANDS), ny, nx)),
                   metn=np.zeros((ny, nx)), fill=np.zeros((ny, nx)), di=np.zeros((len(CLS), len(HOURS))),
                   gsm=np.zeros((ny, nx)), noalt=np.zeros((ny, nx)))
        done = set()
    lai_c = month_lai(m)                                       # FIX-4
    yp, mp = _shift(YEAR, m, -1); lai_p = month_lai(mp, yp)     # FIX-3
    ef_no = ef[sp.IDX_NO]
    for d in range(1, ndays + 1):
        idx_day = [i for i, (dd, _) in enumerate(steps) if dd == d]
        if all(i in done for i in idx_day):
            continue
        doy = datetime.date(YEAR, m, d).timetuple().tm_yday; day_iso = f'{YEAR}-{m:02d}-{d:02d}'; date = YEAR * 1000 + doy
        mets = {h: sample_multi(era5.megan_met_image(ee, region, ee.Date(f'{day_iso}T{h:02d}:00')), MET_BANDS) for h in HOURS}
        D_TEMP = fixes.fill_field_nearest(sample_multi(daily_tmean_image(day_iso), ['D_TEMP'])['D_TEMP'], *fixes.VALID_RANGE['temp'])
        par_f = []
        for h in HOURS:
            mm_, _ = fixes.fill_met_nearest({'temp': mets[h]['TEMP2'], 'par': mets[h]['PAR'], 'wind': mets[h]['WINDSPD'],
                                             'pres': mets[h]['PRES'], 'qv': mets[h]['QV'], 'soilm': mets[h]['SOILM']})
            par_f.append(mm_['par'])
        D_PPFD = np.mean(par_f, axis=0) * fixes.PAR_WM2_TO_PPFD                                 # FIX-2
        for i in idx_day:
            if i in done:
                continue
            h = steps[i][1]; met = mets[h]
            metd = {'temp': met['TEMP2'], 'par': met['PAR'], 'wind': met['WINDSPD'], 'pres': met['PRES'],
                    'qv': met['QV'], 'soilm': met['SOILM'], 'lai': lai_c}
            r = fixes.run_timestep_fixed(domain, ef, metd, day=doy, hour_utc=h, canopy_type_fractions=ctf, d_temp=D_TEMP,
                                         d_ppfd=D_PPFD, laip=lai_p, nlai=12, wilt=WILT, date_yyyyddd=date, strict=True)
            tph = spc.compound_tonnes(spc.flux_to_mol_s(r['compounds'], area[None, :, :]), SEC_PER_STEP)
            acc['t'] += spc.compounds_to_classes(tph); acc['f'] += r['classes']; acc['fill'] += r['filled']
            land = np.isfinite(met['TEMP2']) & (met['TEMP2'] > 150) & (met['TEMP2'] < 340)   # input maps stay land-only
            acc['met'] += np.stack([np.where(land, np.nan_to_num(metd[k]), 0.0) for k in ['temp', 'par', 'wind', 'pres', 'qv', 'soilm']])
            acc['metn'] += land
            acc['di'][:, HOURS.index(h)] += np.nansum(r['classes'], axis=(1, 2))
            if DIAG:
                mf, _ = fixes.fill_met_nearest({k: metd[k] for k in ('temp', 'par', 'wind', 'pres', 'qv', 'soilm')})
                acc['gsm'] += soil.gamma_sm(mf['soilm'], WILT)
                cfno, cfnog = soil.soilnox(mf['temp'], lai_c, domain['LAT'], date, corrected=False)   # MEGANv3.21 SH season
                g_alt = soil.gamma_no(cfno, cfnog, ctf, soil.growseason(date, domain['LAT'], corrected=False)[0])
                acc['noalt'] += no_tonnes(g_alt, ef_no)
            done.add(i)
            st.save_checkpoint(ck, acc, done)
        if d % 5 == 0 or d == ndays:
            print(f'    {tag} day {d}/{ndays}')
    n = len(steps)
    meanmet = acc['met'] / np.where(acc['metn'] > 0, acc['metn'], np.nan)
    _finalize_month(m, tag, acc['t'], acc['f'] / n, meanmet, acc['di'] / ndays)
    np.savez_compressed(os.path.join(OUT, 'inputs', f'fixlog_{tag}.npz'), filled_met_steps=acc['fill'],
                        lai_gap_cells=np.array(LAI_GAP_LOG.get((YEAR, m), 0)), gamsm_mean=acc['gsm'] / n,
                        no_tonnes_mgn321_season=acc['noalt'], steps=np.array(n))
    os.remove(ck)
    print(f'  {tag}: done | vegetated cells with filled met: {int((acc["fill"] > 0).sum())} | LAI gap cells filled: {LAI_GAP_LOG.get((YEAR, m), 0)}')

def setup_year(y):
    """Per-year globals: output folders, MapBiomas growth forms (cached), EF grid, canopy fractions."""
    global YEAR, VEG_YEAR, OUT, gf_np, ef, ctf, VEG_FRAC
    YEAR = y; VEG_YEAR = min(y, 2024); OUT = os.path.join(OUT_ROOT, str(y))
    for sub in ('inputs', 'emissions', 'maps', 'stats', 'hourly', 'checkpoints'):
        os.makedirs(os.path.join(OUT, sub), exist_ok=True)
    gf_cache = os.path.join(OUT, 'inputs', f'growthform_{VEG_YEAR}.npz')
    old = os.path.join(OLD_ROOT, str(y), 'inputs', f'growthform_{VEG_YEAR}.npz')
    if not os.path.exists(gf_cache) and os.path.exists(old):
        shutil.copyfile(old, gf_cache); print('  growth forms: reused', old)
    if os.path.exists(gf_cache):
        gf_np = {k: v for k, v in np.load(gf_cache).items()}
    else:
        print('  MapBiomas growth forms (tiled in-memory pull, ~1-2 min)...')
        gf_np = mapbiomas.growth_form_fractions_gridded(ee, region, domain, VEG_YEAR); np.savez(gf_cache, **gf_np)
    ef = run.build_ef_grid_multiform(cfg, domain, gf_np)
    np.savez_compressed(os.path.join(OUT, 'inputs', f'ef_{y}.npz'), ef=ef.astype('float32'), classes=np.array(CLS))   # for notebook 11 (C3)
    ctf = landcover.canopy_type_fractions(gf_np, veg['tree_canopy_type'])
    VEG_FRAC = sum(gf_np[g] for g in landcover.GROWTH_FORMS)
    geotiff.write_geotiff(os.path.join(OUT, 'inputs', f'growthform_{VEG_YEAR}.tif'),
                          np.stack([gf_np[g] for g in landcover.GROWTH_FORMS]), lon, lat, band_names=list(landcover.GROWTH_FORMS))
    json.dump(dict(year=y, variant=VARIANT, sample_hours=SAMPLE_HOURS, lai_aggregation=LAI_AGG, wilt=WILT,
                   fixes_sha256=FIXES_SHA, lai_correction=LAI_CORR_INFO, patches=PATCHES, ecotype=_eco_f, engine=MEGAN_DIR, ee_project=EE_PROJECT,
                   started=datetime.datetime.utcnow().isoformat() + 'Z'),
              open(os.path.join(OUT, 'provenance_v10.json'), 'w'), indent=1)

def annual_rollup(y):
    from netCDF4 import Dataset
    mfiles = [f for f in sorted(glob.glob(os.path.join(OUT, 'emissions', f'emis_{y}_*.nc'))) if 'test' not in os.path.basename(f)]
    if len(mfiles) < 12:
        print(f'  {y}: {len(mfiles)} months finished, annual roll-up skipped'); return None
    annual = {c: np.zeros((lat.size, lon.size)) for c in CLS}
    for f in mfiles:
        with Dataset(f) as ds:
            for c in CLS:
                a = np.array(ds.variables[c][:], dtype='float64'); annual[c] += np.where(np.isfinite(a) & (a > -9990), a, 0.0)
    stack = np.stack([annual[c] for c in CLS])
    ncio.write_emissions(os.path.join(OUT, 'emissions', f'annual_emis_{y}.nc'), stack, CLS, lon, lat, units='tonnes/year',
                         global_attrs={'year': y, 'months': '1-12', 'mechanism': MECH, 'megan_sc_fixes_sha256': FIXES_SHA, 'variant': VARIANT,
                                                             'lai_source': LAI_SOURCE, 'megan_sc_laicorr_sha256': LAICORR_SHA})
    geotiff.write_species_stack(os.path.join(OUT, 'emissions', f'annual_emis_{y}.tif'), stack, lon, lat, CLS)
    tot = pd.DataFrame([{'species': c, 'annual_tonnes_SC': float(np.nansum(sc_clip(annual[c]))), 'family': FAM_OF.get(c, ('', ''))[0]} for c in CLS])
    tot.to_csv(os.path.join(OUT, 'stats', f'annual_by_family_{y}.csv'), index=False)
    return tot.set_index('species')['annual_tonnes_SC']

# %%
# 9 · Run (safe to re-run after any interruption) -----------------------------------------------------------
summary = {}
for y in YEARS:
    t0 = time.time(); print(f'=== {y} ===')
    setup_year(y); _lai_cache.clear(); LAI_GAP_LOG.clear()
    for m in MONTHS:
        run_month(m)
    tot = annual_rollup(y)
    if tot is not None:
        summary[y] = tot
        print(f'  {y} SC totals (kt): BVOC {tot.drop(["NO", "CO"]).sum() / 1e3:.1f} | ISOP {tot["ISOP"] / 1e3:.1f} | '
              f'SQT {(tot["SQT_HR"] + tot["SQT_LR"]) / 1e3:.1f} | NO {tot["NO"] / 1e3:.2f} | {(time.time() - t0) / 60:.0f} min')
    _lai_raw.clear()
if summary:
    pd.DataFrame(summary).T.to_csv(os.path.join(OUT_ROOT, f'annual_totals_SC_v10_{min(summary)}_{max(summary)}.csv'))

# %%
# 10 · Quick self-check of the years just run (the full check is notebook 11) -----------------------------------
rows = []
for y in summary:
    from netCDF4 import Dataset
    with Dataset(os.path.join(OUT_ROOT, str(y), 'emissions', f'annual_emis_{y}.nc')) as ds:
        E = {c: np.nan_to_num(np.array(ds.variables[c][:], dtype='float64')) for c in CLS}
    bv = sum(E[c] for c in CLS if c not in ('NO', 'CO'))
    tree = np.load(os.path.join(OUT_ROOT, str(y), 'inputs', f'growthform_{min(y, 2024)}.npz'))['tree']
    sig = SC_MASK & (tree > 0.1) & (bv > 0) & (E['ISOP'] == 0)               # the 340 K / PAR = 0 signature
    logs = sorted(glob.glob(os.path.join(OUT_ROOT, str(y), 'inputs', f'fixlog_{y}_*.npz')))
    filled = sum(int((np.load(f)['filled_met_steps'] > 0).sum()) for f in logs)
    rows.append(dict(year=y, signature_cells=int(sig.sum()), fixlogs=len(logs), filled_cell_months=filled))
print(pd.DataFrame(rows).to_string(index=False))
print('signature_cells must be 0 in every year; fixlogs must be 12.')

# %%
# 11 · Comparison with the production run: only the LAI may differ -----------------------------------------------
import rasterio
def _bands(path):                                         # all bands, south first
    with rasterio.open(path) as s:
        return np.flip(s.read().astype('float64'), axis=1)
def _lai_mean(root, y, m):
    a = laicorr.read_south_first(os.path.join(root, str(y), 'inputs', f'lai_{y}_{m:02d}.tif'))
    return float(a[SC_MASK & (a > 0)].mean())
cmp_rows, met_max = [], 0.0
for y in YEARS:
    for m in MONTHS:
        a = _bands(os.path.join(ROOT_V10, str(y), 'inputs', f'met_{y}_{m:02d}.tif'))
        b = _bands(os.path.join(OUT_ROOT, str(y), 'inputs', f'met_{y}_{m:02d}.tif'))
        ok = SC_MASK[None] & np.isfinite(a) & np.isfinite(b); met_max = max(met_max, float(np.abs(a[ok] - b[ok]).max()))
print(f'meteorology (monthly means of the six input bands, SC cells): max abs difference {met_max:.3g} (must be 0)')
GROUPS = {'BVOC': [c for c in CLS if c not in ('NO', 'CO')], 'ISOP': ['ISOP'], 'MT': [c for c in CLS if c.startswith('MT_')],
          'SQT': [c for c in CLS if c.startswith('SQT')], 'OVOC': ['MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD'], 'LVOC': ['LVOC'],
          'MBO': ['MBO'], 'NO': ['NO'], 'CO': ['CO']}
for y in YEARS:
    t0 = pd.read_csv(os.path.join(ROOT_V10, str(y), 'stats', f'annual_by_family_{y}.csv')).set_index('species')['annual_tonnes_SC']
    t1 = pd.read_csv(os.path.join(OUT_ROOT, str(y), 'stats', f'annual_by_family_{y}.csv')).set_index('species')['annual_tonnes_SC']
    l0 = np.mean([_lai_mean(ROOT_V10, y, m) for m in MONTHS]); l1 = np.mean([_lai_mean(OUT_ROOT, y, m) for m in MONTHS])
    cmp_rows.append(dict(year=y, item='LAI (state mean of the monthly means)', modis_run=round(l0, 3), lai_viirs_run=round(l1, 3), change_pct=round(100 * (l1 / l0 - 1), 1)))
    for g, cl in GROUPS.items():
        a0, a1 = float(t0[cl].sum()), float(t1[cl].sum())
        cmp_rows.append(dict(year=y, item=f'{g} (t/yr, SC, before rho)', modis_run=round(a0), lai_viirs_run=round(a1), change_pct=round(100 * (a1 / a0 - 1), 1)))
CMP = pd.DataFrame(cmp_rows); CMP.to_csv(os.path.join(OUT_ROOT, 'lai_viirs_rerun_comparison.csv'), index=False)
print(CMP.to_string(index=False))
print('Isoprene in the texts is multiplied by rho = 0.96 (notebooks 11 and 13); the percentage changes are not affected.')

# %%
# 12 · Optional: make the corrected years part of the production run --------------------------------------------------
INSTALL = False                         # <-- set True only after reading cell 11
ARCHIVE = os.path.join(OUT_BASE, 'megan_sc_outputs_v10_modis_lai_2024_2025')   # where the original 2024-2025 folders go
if INSTALL:
    os.makedirs(ARCHIVE, exist_ok=True)
    for y in YEARS:
        src_new, dst = os.path.join(OUT_ROOT, str(y)), os.path.join(ROOT_V10, str(y))
        assert os.path.exists(os.path.join(src_new, 'emissions', f'annual_emis_{y}.nc')), f'{y} not finished'
        assert not os.path.exists(os.path.join(ARCHIVE, str(y))), f'{ARCHIVE}/{y} exists: the original was already moved'
        shutil.move(dst, os.path.join(ARCHIVE, str(y)))      # original MODIS-LAI year kept in the archive folder
        shutil.copytree(src_new, dst)
        print(f'{y}: original moved to {ARCHIVE}/{y}; corrected year copied into {ROOT_V10}')
    print('Next: run notebooks 11, 13 and 14 again on', ROOT_V10)
else:
    print('INSTALL = False: the production run in', ROOT_V10, 'is unchanged.')
