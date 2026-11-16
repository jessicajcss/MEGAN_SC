# %% [markdown]
# # 10 · Production run: biogenic emissions for any region and period
#
# This notebook computes the emission inventory described by a run configuration (`config/*.yaml`): the grid, the
# region, the period, the land-cover source and the MEGAN input tables. For every month of the period it samples
# ERA5-Land, MODIS LAI and the land cover in Google Earth Engine on the model grid, runs the MEGANv3.21 engine of
# `megan_sc` (canopy, activity factors, soil NO, speciation) every `sample_hours` hours, and writes the monthly and
# annual emission grids, the inputs as sampled and a log of each procedure. The input and run procedures live in
# `src/megan_sc/fixes.py`, which this notebook writes into the repository folder (cell 2); its SHA-256 is stamped in
# every output file, so each output can be traced to the code that produced it.
#
# | Procedure (label in `fixes.py`) | What it does |
# |---|---|
# | FIX-1 | cells without ERA5-Land data (the reanalysis covers land only) take all meteorological variables from the nearest valid land cell (one donor for all bands); the run stops if a vegetated cell keeps invalid meteorology |
# | FIX-2 | D_TEMP = 24-h mean 2 m temperature of the UTC day (24 hourly ERA5-Land records, as daymet.f); D_PPFD from the day's PAR |
# | FIX-3 | leaf age with the previous month's LAI as LAIp (NLAI = 12) |
# | FIX-4 | monthly MODIS LAI gaps filled from the adjacent valid months |
# | FIX-5 | no-data pixels of the MEGAN ecotype map take the nearest valid ecotype |
# | FIX-6 | GAMTP applies the canopy-depth factor to the sunlit term only, as megvea.f does |
# | FIX-7 | air temperature capped at 303 K before the soil-temperature proxy, as soilnox.F does |
# | FIX-8 | ERA5-Land radiation de-accumulated: previous hour subtracted at every hour except 01 UTC (ECMWF ERA5-Land convention; the 00 UTC field holds the previous day's 24-h total) |
#
# Settings that follow the MEGANv3.21 code or are fixed by the configuration: the soil-NO growing season of each
# latitude band (`soil.growseason`), constant light-dependent fractions for the monoterpene classes, a constant wilting
# point (`datasets.soil.constant_wilt`), the canopy type of trees (`datasets.vegetation.tree_canopy_type`), and the
# soil-moisture and bidirectional-exchange factors on.
#
# **Configurations.** `config/santa_catarina_2003_2025.yaml` (default) reproduces the Santa Catarina example run shipped
# in `data/run_v10`. For another region, `run_your_region/region_setup.ipynb` writes `config/<region>_run.yaml`; set
# `RUN_CONFIG` to that file. The tutorial (`docs/tutorial/06_run_the_model.md`) explains every setting.
#
# **How to run.** Set `RUN_CONFIG` and `EE_PROJECT` (cell 3), optionally restrict `YEARS`, and run all cells.
# Finished months are skipped and an interrupted month resumes from its checkpoint, so after a Colab disconnect just
# run all cells again. Outputs go to `OUT_BASE/<output.dir of the configuration>`. When the run finishes, run
# `11_verify_run.ipynb` (checks) and `25_results_any_region.ipynb` (tables, maps, trends) with the same `RUN_CONFIG`.

# %%
# 1 · Setup ---------------------------------------------------------------------------------------------
__FILE:notebooks/src/_setup_repo.py__
%pip -q install earthengine-api numpy scipy xarray netCDF4 matplotlib pyyaml pandas rasterio geopandas
import json, time, datetime, calendar, logging, hashlib, re, shutil
import numpy as np, pandas as pd, matplotlib.pyplot as plt
MEGAN_DIR = REPO_DIR                                    # the repository folder found above
for _lg in ('google_auth_httplib2', 'googleapiclient.http', 'google.auth'):
    logging.getLogger(_lg).setLevel(logging.ERROR)

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
# 3 · Settings -------------------------------------------------------------------------------------------
RUN_CONFIG   = 'config/santa_catarina_2003_2025.yaml'   # <-- the run configuration (path relative to the repository)
EE_PROJECT   = None                     # <-- your Earth Engine project id; None = earth_engine.project of RUN_CONFIG
OUT_BASE     = '/content/drive/MyDrive' if IN_COLAB else os.path.dirname(REPO_DIR)   # <-- the output folder goes here
YEARS        = None                     # <-- None = every year of the period; or a few, e.g. [2019], and come back later
MONTHS       = None                     # None = every month of the period; or calendar months, e.g. [1, 4, 7, 10]
TEST_DAYS    = None                     # e.g. 2 for a quick preview (outputs get a _testNd suffix, not used downstream)
VARIANT      = 'base'                   # 'base' | 'hourly' or 'lai_eedefault' (sampling tests of the Santa Catarina example)
LAI_AGG      = 'mean'                   # 'mean': area mean of QC-valid 500 m MODIS pixels per grid cell (explicit
                                        #  reduceResolution) | 'ee_default': Earth Engine default resampling
MAKE_MAPS    = False                    # monthly and annual PNG maps (slow on Drive); notebook 25 draws the result maps
DIAG         = True                     # cheap diagnostics in the fix log (soil-NO season alternative, mean GAMSM)

# %%
# 4 · Configuration, Earth Engine, engine patches and FIX-5 ------------------------------------------------
from megan_sc import domain as dom, run, landcover, stats as st, fixes
from megan_sc import gee
from megan_sc.io import geotiff, netcdf as ncio
from megan_sc.engine import speciation as spc, species as sp, soil
from megan_sc.gee import era5, lai as laimod, mapbiomas, modis_lc, region_from_bounds
def in_repo(p):
    """Absolute path of a configuration entry (relative entries are relative to the repository)."""
    return p if (not p or os.path.isabs(p)) else os.path.join(MEGAN_DIR, p)
cfg = run.load_config(in_repo(RUN_CONFIG))
for k in ('ef_table', 'speciation_dir', 'ecotype'):
    cfg['datasets']['emission_factors'][k] = in_repo(cfg['datasets']['emission_factors'][k])
if cfg['domain'].get('mode') == 'shapefile':
    cfg['domain']['shapefile'] = in_repo(cfg['domain']['shapefile'])
REG = cfg.get('region', {})
REGION_NAME = REG.get('name', cfg['run']['name']); REGION_ABBR = REG.get('abbr', cfg['run']['name'])
BOUNDARY = in_repo(REG.get('boundary') or cfg['domain'].get('shapefile'))
ZONES, ZONES_ID = in_repo(REG.get('zones')), REG.get('zones_id')
veg = cfg['datasets']['vegetation']
GF_SOURCE = veg.get('growth_form_source', 'mapbiomas')
assert GF_SOURCE in ('mapbiomas', 'modis'), f'growth_form_source {GF_SOURCE!r}: this notebook reads mapbiomas or modis'
LC_FIRST, LC_LAST = {'mapbiomas': (1985, 2024), 'modis': (2001, 2024)}[GF_SOURCE]   # years of the products in Earth Engine
LC_FIRST, LC_LAST = int(veg.get('first_year', LC_FIRST)), int(veg.get('last_year', LC_LAST))
def veg_year(y):
    """Land-cover year used for year y: y itself, or the nearest year of the land-cover product."""
    return min(max(int(y), LC_FIRST), LC_LAST)
_p0, _p1 = cfg['period']['start'], cfg['period']['end']
PERIOD_MONTHS = [(y, m) for y in range(int(_p0[:4]), int(_p1[:4]) + 1) for m in range(1, 13)
                 if _p0 <= f'{y}-{m:02d}' <= _p1]
YEARS = sorted({y for y, _ in PERIOD_MONTHS}) if YEARS is None else [y for y in YEARS if any(yy == y for yy, _ in PERIOD_MONTHS)]
SAMPLE_HOURS = int(cfg['period'].get('sample_hours', 3))   # production setting 3: 8 steps per day
if VARIANT == 'hourly':
    SAMPLE_HOURS = 1
if VARIANT == 'lai_eedefault':
    LAI_AGG = 'ee_default'
assert VARIANT in ('base', 'hourly', 'lai_eedefault'), VARIANT
MECH = cfg.get('species', {}).get('mechanism', 'CB6X')
OUT_DIR_NAME = cfg.get('output', {}).get('dir', f"megan_sc_outputs_{cfg['run']['name']}")
OUT_TAG = cfg.get('output', {}).get('tag', cfg['run']['name'])
OUT_ROOT = os.path.join(OUT_BASE, OUT_DIR_NAME + ('' if VARIANT == 'base' else f'_{VARIANT}'))
os.makedirs(OUT_ROOT, exist_ok=True)
EE_PROJECT = EE_PROJECT or cfg.get('earth_engine', {}).get('project')
ee = gee.initialize(EE_PROJECT, service_account_key=cfg.get('earth_engine', {}).get('service_account_key'))
print('Earth Engine ready:', EE_PROJECT)
PATCHES = fixes.apply_engine_patches()                  # FIX-6, FIX-7, FIX-8
assert fixes.patches_active(), 'engine patches not active'
_eco = cfg['datasets']['emission_factors']['ecotype']
_stem = os.path.splitext(os.path.basename(_eco))[0]
_eco_f = _eco if _stem.endswith('_filled') else os.path.join(os.path.dirname(_eco), _stem + '_filled.nc')   # FIX-5
if not os.path.exists(_eco_f):
    print('FIX-5: ecotype no-data pixels filled:', fixes.fill_ecotype_nearest(_eco, _eco_f))
cfg['datasets']['emission_factors']['ecotype'] = _eco_f
domain = dom.build_domain(cfg); region = region_from_bounds(domain['bounds'])
res = cfg['domain']['resolution_deg']
lon, lat, area = domain['lon'], domain['lat'], domain['area']
WILT = float(cfg['datasets'].get('soil', {}).get('constant_wilt', 0.12))
HOURS = list(range(0, 24, SAMPLE_HOURS)); SEC_PER_STEP = SAMPLE_HOURS * 3600.0
CLS = list(sp.NAMES); MET_BANDS = ['TEMP2', 'PAR', 'WINDSPD', 'PRES', 'QV', 'SOILM']
# soil-NO compound(s) for the DIAG alternative (EF_NO x GAMNO x speciation fraction)
_NO_MAP = [(int(spc.MAP_SPCA[s]), float(spc.MAP_EFFS[s])) for s in range(spc.MAP_SPCA.size) if spc.MAP_MG20[s] == sp.IDX_NO]
print(f'{REGION_NAME}: years {YEARS[0]}-{YEARS[-1]} ({len(PERIOD_MONTHS)} months in the period) | {24 // SAMPLE_HOURS} steps/day '
      f'| land cover {GF_SOURCE} {LC_FIRST}-{LC_LAST} | LAI_AGG={LAI_AGG} | -> {OUT_ROOT}')
print('patches:', PATCHES, '| grid', lat.size, 'x', lon.size, '| wilt', WILT, '| NO compounds', _NO_MAP)

# %%
# 5 · Reporting setup: region mask, palettes, map helper, zones ----------------------------------------------
import geopandas as gpd
from rasterio.features import rasterize as _rasterize
from mpl_toolkits.axes_grid1 import make_axes_locatable
DPI = 300
ASPECT = 1.0 / np.cos(np.deg2rad(float(np.mean(lat))))
assert BOUNDARY and os.path.exists(BOUNDARY), f'region boundary not found: {BOUNDARY}'
_gg = gpd.read_file(BOUNDARY).to_crs('EPSG:4326')
_BND = _gg.geometry.union_all() if hasattr(_gg.geometry, 'union_all') else _gg.geometry.unary_union
_tr = geotiff._transform(lon, lat)
MASK = np.flipud(_rasterize([(_BND, 1)], out_shape=(lat.size, lon.size), transform=_tr, fill=0, dtype='uint8')).astype(bool)
domain['mask'] = MASK                                   # a cell belongs to the region when its centre lies inside
def region_clip(field): return np.where(MASK, field, np.nan)
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
    im = ax.imshow(np.ma.masked_invalid(region_clip(field)), origin='lower', extent=[lon[0], lon[-1], lat[0], lat[-1]],
                   cmap=cm_, vmin=vmin, vmax=vmax, aspect=ASPECT)
    draw_boundary(ax); ax.set_title(title); ax.axis('off')
    fig.colorbar(im, cax=make_axes_locatable(ax).append_axes('right', size='4.5%', pad=0.05))
    fig.tight_layout(); fig.savefig(path_png, dpi=DPI, bbox_inches='tight'); plt.close(fig)
_ZONE = _IDS = None
if ZONES:
    try:
        _g = gpd.read_file(ZONES).to_crs('EPSG:4326')
        _idc = ZONES_ID if ZONES_ID in _g.columns else ('CD_MUN' if 'CD_MUN' in _g.columns else _g.columns[0])
        _ZONE, _IDS = st.rasterize_municipalities(_g, lon, lat, _idc)
        print(f'zones: {len(_IDS)} polygons ({os.path.basename(ZONES)}, id {_idc})')
    except Exception as e:
        print('zones file not read (per-zone tables skipped):', e)
print(f'{REGION_ABBR} mask: {int(MASK.sum())} of {MASK.size} cells')

# %%
# 6 · Earth Engine samplers: meteorology, 24-h mean temperature, LAI (FIX-4 gap fill; explicit aggregation) --
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
# 7 · One month and one year ------------------------------------------------------------------------------
def _finalize_month(m, tag, tonnes, meanflux, meanmet, diurnal):
    ncio.write_emissions(os.path.join(OUT, 'emissions', f'emis_{tag}.nc'), tonnes, CLS, lon, lat,
                         units='tonnes/month', global_attrs={'year': YEAR, 'month': m, 'mechanism': MECH,
                                                             'megan_sc_fixes_sha256': FIXES_SHA, 'variant': VARIANT})
    geotiff.write_species_stack(os.path.join(OUT, 'emissions', f'emis_{tag}.tif'), tonnes, lon, lat, CLS)
    for nm in ['ISOP', 'MT_PINE', 'SQT_HR', 'NO', 'CO']:
        save_map(tonnes[sp.index_of(nm)], f'{nm} {tag} (t/month, {REGION_ABBR})', os.path.join(OUT, 'maps', f'{nm}_{tag}.png'),
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

def growth_forms(year):
    """Growth-form fractions (tree, shrub, herb, crop, veg) of a land-cover year on the grid, from Earth Engine."""
    if GF_SOURCE == 'mapbiomas':
        cid = mapbiomas.resolve_collection_id(ee, year, veg.get('collection'))
        return mapbiomas.growth_form_fractions_gridded(ee, region, domain, year, collection_id=cid), cid
    img = modis_lc.growth_form_fractions(ee, region, year, out_scale_m=res * 111000.0)
    bands = list(landcover.GROWTH_FORMS) + ['veg']
    g = mapbiomas.sample_image_tiled(ee, img, domain, bands, ntiles=4)
    return {b: np.nan_to_num(g[b], nan=0.0) for b in bands}, None

def setup_year(y):
    """Per-year globals: output folders, growth forms (cached), EF grid, canopy fractions."""
    global YEAR, VEG_YEAR, OUT, gf_np, ef, ctf, VEG_FRAC
    YEAR = y; VEG_YEAR = veg_year(y); OUT = os.path.join(OUT_ROOT, str(y))
    for sub in ('inputs', 'emissions', 'maps', 'stats', 'hourly', 'checkpoints'):
        os.makedirs(os.path.join(OUT, sub), exist_ok=True)
    gf_cache = os.path.join(OUT, 'inputs', f'growthform_{VEG_YEAR}.npz')
    cid = None
    if os.path.exists(gf_cache):
        gf_np = {k: v for k, v in np.load(gf_cache).items()}
    else:
        print(f'  {GF_SOURCE} growth forms of {VEG_YEAR} (tiled in-memory pull, ~1-2 min)...')
        gf_np, cid = growth_forms(VEG_YEAR); np.savez(gf_cache, **gf_np)
    ef = run.build_ef_grid_multiform(cfg, domain, gf_np)
    np.savez_compressed(os.path.join(OUT, 'inputs', f'ef_{y}.npz'), ef=ef.astype('float32'), classes=np.array(CLS))   # for notebook 11 (C3)
    ctf = landcover.canopy_type_fractions(gf_np, veg['tree_canopy_type'])
    VEG_FRAC = sum(gf_np[g] for g in landcover.GROWTH_FORMS)
    geotiff.write_geotiff(os.path.join(OUT, 'inputs', f'growthform_{VEG_YEAR}.tif'),
                          np.stack([gf_np[g] for g in landcover.GROWTH_FORMS]), lon, lat, band_names=list(landcover.GROWTH_FORMS))
    json.dump(dict(year=y, variant=VARIANT, sample_hours=SAMPLE_HOURS, lai_aggregation=LAI_AGG, wilt=WILT,
                   fixes_sha256=FIXES_SHA, patches=PATCHES, ecotype=_eco_f, engine=MEGAN_DIR, ee_project=EE_PROJECT,
                   config=RUN_CONFIG, region=REGION_NAME, growth_form_source=GF_SOURCE, land_cover_year=VEG_YEAR,
                   mapbiomas_collection_id=cid, started=datetime.datetime.utcnow().isoformat() + 'Z'),
              open(os.path.join(OUT, 'provenance_v10.json'), 'w'), indent=1)

def annual_rollup(y):
    """Annual sum of the 12 monthly files and the totals of the region (skipped when a month is missing)."""
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
                         global_attrs={'year': y, 'months': '1-12', 'mechanism': MECH, 'megan_sc_fixes_sha256': FIXES_SHA, 'variant': VARIANT})
    geotiff.write_species_stack(os.path.join(OUT, 'emissions', f'annual_emis_{y}.tif'), stack, lon, lat, CLS)
    tot = pd.DataFrame([{'species': c, f'annual_tonnes_{REGION_ABBR}': float(np.nansum(region_clip(annual[c]))), 'family': FAM_OF.get(c, ('', ''))[0]} for c in CLS])
    tot.to_csv(os.path.join(OUT, 'stats', f'annual_by_family_{y}.csv'), index=False)
    return tot.set_index('species')[f'annual_tonnes_{REGION_ABBR}']

# %%
# 8 · Run (safe to re-run after any interruption) -----------------------------------------------------------
summary = {}
for y in YEARS:
    t0 = time.time(); print(f'=== {y} ===')
    setup_year(y); _lai_cache.clear(); LAI_GAP_LOG.clear()
    for m in [mm for yy, mm in PERIOD_MONTHS if yy == y and (MONTHS is None or mm in MONTHS)]:
        run_month(m)
    tot = annual_rollup(y)
    if tot is not None:
        summary[y] = tot
        print(f'  {y} {REGION_ABBR} totals (kt): BVOC {tot.drop(["NO", "CO"]).sum() / 1e3:.1f} | ISOP {tot["ISOP"] / 1e3:.1f} | '
              f'SQT {(tot["SQT_HR"] + tot["SQT_LR"]) / 1e3:.1f} | NO {tot["NO"] / 1e3:.2f} | {(time.time() - t0) / 60:.0f} min')
    _lai_raw.clear()
if summary:
    pd.DataFrame(summary).T.to_csv(os.path.join(OUT_ROOT, f'annual_totals_{OUT_TAG}_{min(summary)}_{max(summary)}.csv'))

# %%
# 9 · Quick self-check of the years just run (the full check is notebook 11) -----------------------------------
rows = []
for y in summary:
    from netCDF4 import Dataset
    with Dataset(os.path.join(OUT_ROOT, str(y), 'emissions', f'annual_emis_{y}.nc')) as ds:
        E = {c: np.nan_to_num(np.array(ds.variables[c][:], dtype='float64')) for c in CLS}
    bv = sum(E[c] for c in CLS if c not in ('NO', 'CO'))
    tree = np.load(os.path.join(OUT_ROOT, str(y), 'inputs', f'growthform_{veg_year(y)}.npz'))['tree']
    sig = MASK & (tree > 0.1) & (bv > 0) & (E['ISOP'] == 0)                  # the 340 K / PAR = 0 signature
    logs = sorted(glob.glob(os.path.join(OUT_ROOT, str(y), 'inputs', f'fixlog_{y}_*.npz')))
    filled = sum(int((np.load(f)['filled_met_steps'] > 0).sum()) for f in logs)
    rows.append(dict(year=y, signature_cells=int(sig.sum()), fixlogs=len(logs), filled_cell_months=filled))
print(pd.DataFrame(rows).to_string(index=False))
print('signature_cells must be 0 in every year; fixlogs must be 12.')
