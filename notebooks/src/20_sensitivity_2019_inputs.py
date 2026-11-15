# %% [markdown]
# # 20 · One-at-a-time sensitivity runs of 2019: temperature at the terrain of each cell, light, canopy type, wilting point
#
# **Why.** The comparison with BRAIN (notebook 18) left isoprene 21 % lower here than in BRAIN over the state, with the
# largest gaps on the coast (Grande Florianópolis and Sul, ratio 0.57) and near agreement in the west and on the plateau
# (Oeste and Serrana, 1.01-1.06). Notebook 23 tests explanations with the data at hand; this notebook supplies the
# model experiments. Each variant changes **one** input or setting of the production run of 2019 and keeps everything
# else, so its effect is the effect of that choice alone.
#
# | Variant | Change from the production run | Why it is tested |
# |---|---|---|
# | `base` | none | reproduces the production run (check) and is the reference of every other variant |
# | `lapse` | 2 m temperature (hourly and daily mean) corrected from the mean elevation of the ERA5-Land cell to the mean elevation of each 0.025° cell: ΔT = −Γ (z_cell − z_ERA5) | ERA5-Land (≈9 km) cannot resolve the valleys and ridges of the Serra Geral and the coastal escarpment; BRAIN's WRF runs at 4 km |
# | `par050` | PAR = 0.50 × shortwave instead of 0.45 × shortwave (hourly and daily mean) | 0.5 is the default of MEGANv3.21 when PAR comes from MCIP (`SRAD_TO_PAR` in `src/MET2MGN/readmcip.f`); the distributed run script does not change it |
# | `canopy_tropical` | tree canopy type 2 (tropical forest) instead of 3 (temperate broadleaf) | canopy structure of the Atlantic Forest; MEGAN's own canopy-type files use several types |
# | `canopy_needle_silv` | the silviculture share of each cell (MapBiomas class 9, 2019) as needleleaf trees (type 1), the rest of the trees as type 3 | Pinus plantations are needleleaf; the run treats all silviculture as broadleaf (`landcover.py`) |
# | `wilt_clay` | wilting point 0.286 m³ m⁻³ (WWLT of clay, MEGSEA.EXT) instead of 0.12 | upper bound of the soil-moisture factor; MEGANv3.21 takes the wilting point from the soil type of the meteorological model |
# | `wilt_texture` (optional) | wilting point from MapBiomas Soil sand and clay through the USDA texture classes and WWLT | needs the MapBiomas Soil assets (`gee/soil_texture.py`; check the asset names first) |
#
# All variants share the same Earth Engine samples of every time step, so the differences contain no sampling noise.
#
# **Lapse rate.** Γ comes from notebook 21 (`lapse_rate.json`: ERA5-Land minus INMET station temperature against the
# station-minus-cell elevation) when that file exists, else 6.5 K km⁻¹ (standard atmosphere). The elevation of each
# 0.025° cell is the mean of SRTM (CGIAR SRTM 90 m v4) over the cell. The elevation of the ERA5-Land cell is the mean of
# SRTM over the 0.1° ERA5-Land pixel (Earth Engine), or, if you set `ZOROG_NC`, the ERA5-Land orography itself
# (geopotential file of the ERA5-Land documentation divided by g = 9.80665 m s⁻²). Pressure and humidity are not
# changed. The meteorological fields are filled (FIX-1) before ΔT is added, so coastal cells filled from a neighbour get
# their own correction.
#
# **Modes.** `MODE = 'design'`: the 28 days of test S1 (days 1-7 of January, April, July and October), about 40 min
# for six variants. `MODE = 'full'`: every day of 2019, several hours (run all cells again after a disconnection;
# finished days are skipped). Monthly class grids of every variant go to `OUT_ROOT/<variant>/` for notebook 23.
#
# **Inputs reused from the production run** (`ROOT_V10/2019/inputs`): MapBiomas growth forms, and the monthly LAI the
# run used (GeoTIFF, float32; December 2018 for the leaf age of January). Only ERA5-Land and the DEM come from Earth
# Engine. The engine is the same as notebooks 10 and 17 (`fixes.py` with SHA-256 prefix a933b29602929a82).

# %%
# 1 · Settings -------------------------------------------------------------------------------------------------------------
%pip -q install earthengine-api numpy scipy xarray netCDF4 matplotlib pyyaml pandas rasterio geopandas
from google.colab import drive; drive.mount('/content/drive')
MEGAN_DIR = '/content/drive/MyDrive/megan-sc'                      # <-- engine repo in Drive (as notebooks 10 and 17)
REPO_DIR  = '/content/drive/MyDrive/MEGAN_SC'                      # <-- this repository (MapBiomas 2019 at 120 m, boundaries)
ROOT_V10  = '/content/drive/MyDrive/megan_sc_outputs_v10'          # <-- production run (2019 inputs and outputs)
OUT_ROOT  = ROOT_V10 + '/analysis_v10/sensitivity_2019'            # outputs of this notebook (read by notebook 23)
N21_JSON  = ROOT_V10 + '/analysis_v10/era5land_vs_inmet/lapse_rate.json'   # notebook 21 (optional)
ZOROG_NC  = None                       # optional: ERA5-Land orography (geopotential netCDF, variable z); None = SRTM mean
YEAR      = 2019
MODE      = 'design'                   # 'design' (28 days of test S1) or 'full' (every day of 2019)
VARIANTS  = ['base', 'lapse', 'par050', 'canopy_tropical', 'canopy_needle_silv', 'wilt_clay']   # add 'wilt_texture' if wanted
GAMMA_K_PER_KM = None                  # None: notebook 21 if available, else 6.5
PAR_FRAC_RUN, PAR_FRAC_ALT = 0.45, 0.50
WILT_ALT  = 0.286                      # MEGSEA.EXT WWLT, clay
MET_CACHE = True                       # keep the Earth Engine samples of each day (float32, about 7 MB per day)
SAMPLE_HOURS = 3; MECH = 'CB6X'
import sys, os, json, time, datetime, calendar, logging, hashlib, glob, re, shutil
import numpy as np, pandas as pd
sys.path.insert(0, MEGAN_DIR + '/src')
for _lg in ('google_auth_httplib2', 'googleapiclient.http', 'google.auth'):
    logging.getLogger(_lg).setLevel(logging.ERROR)
DAYS = ([(m, d) for m in (1, 4, 7, 10) for d in range(1, 8)] if MODE == 'design'
        else [(m, d) for m in range(1, 13) for d in range(1, calendar.monthrange(YEAR, m)[1] + 1)])
OUT = os.path.join(OUT_ROOT, MODE); os.makedirs(OUT, exist_ok=True)
print(f'{MODE}: {len(DAYS)} days, {len(VARIANTS)} variants -> {OUT}')

# %%
# 2 · Write megan_sc/fixes.py (same file as notebooks 10 and 17) --------------------------------------------------------------
FIXES_PY = r'''__FIXES_PY__'''
FIXES_SHA = hashlib.sha256(FIXES_PY.encode('utf-8')).hexdigest()[:16]
assert FIXES_SHA == 'a933b29602929a82', FIXES_SHA
with open(os.path.join(MEGAN_DIR, 'src', 'megan_sc', 'fixes.py'), 'w', encoding='utf-8') as _f:
    _f.write(FIXES_PY)
for _m in [m for m in list(sys.modules) if m.startswith('megan_sc')]:
    del sys.modules[_m]
from megan_sc import fixes
print('fixes.py written | sha256[:16] =', FIXES_SHA)

# %%
# 3 · Earth Engine ------------------------------------------------------------------------------------------------------------
EE_PROJECT = 'your-ee-project'         # <-- your Earth Engine (Google Cloud) project id
import ee
ee.Authenticate(); ee.Initialize(project=EE_PROJECT)
print('Earth Engine ready:', EE_PROJECT)

# %%
# 4 · Engine and grid (as cell 5 of notebook 17) ----------------------------------------------------------------------------------
from megan_sc import domain as dom, run, landcover, stats as st
from megan_sc.engine import speciation as spc, species as sp, soil
from megan_sc.gee import era5, mapbiomas, region_from_bounds
PATCHES = fixes.apply_engine_patches(); assert fixes.patches_active()
_cfgp = MEGAN_DIR + '/config/example_config.yaml'
cfg = run.load_config(_cfgp if os.path.exists(_cfgp) else MEGAN_DIR + '/example_config.yaml')
cfg['domain']['mode'] = 'bbox'
for k in ('ef_table', 'speciation_dir', 'ecotype'):
    v = cfg['datasets']['emission_factors'].get(k)
    if v and not os.path.isabs(v):
        cfg['datasets']['emission_factors'][k] = os.path.join(MEGAN_DIR, v)
_eco = cfg['datasets']['emission_factors']['ecotype']
_eco_f = os.path.join(os.path.dirname(_eco), 'sc_ecotype_filled.nc')
if not os.path.exists(_eco_f):
    print('FIX-5: ecotype no-data pixels filled:', fixes.fill_ecotype_nearest(_eco, _eco_f))
cfg['datasets']['emission_factors']['ecotype'] = _eco_f
domain = dom.build_domain(cfg); region = region_from_bounds(domain['bounds'])
res = cfg['domain']['resolution_deg']; veg = cfg['datasets']['vegetation']
lon, lat, area = domain['lon'], domain['lat'], domain['area']
WILT = float(cfg['datasets'].get('soil', {}).get('constant_wilt', 0.12))
HOURS = list(range(0, 24, SAMPLE_HOURS)); SEC_PER_STEP = SAMPLE_HOURS * 3600.0
CLS = list(sp.NAMES); MET_BANDS = ['TEMP2', 'PAR', 'WINDSPD', 'PRES', 'QV', 'SOILM']
assert WILT == 0.12 and veg['tree_canopy_type'] == 3, (WILT, veg['tree_canopy_type'])
import geopandas as gpd, rasterio
from rasterio.features import rasterize
from megan_sc.io import geotiff
BND = REPO_DIR + '/data/boundaries'
_sc = gpd.read_file(BND + '/sc_state_ibge.gpkg').to_crs('EPSG:4326')
_tr = geotiff._transform(lon, lat)
SC_MASK = np.flipud(rasterize([(_sc.geometry.union_all() if hasattr(_sc.geometry, 'union_all') else _sc.geometry.unary_union, 1)],
                              out_shape=(lat.size, lon.size), transform=_tr, fill=0, dtype='uint8')).astype(bool)
_meso = gpd.read_file(BND + '/sc_mesoregions_2024.gpkg').to_crs('EPSG:4326')
MESO = {r['name_meso']: np.flipud(rasterize([(r.geometry, 1)], out_shape=(lat.size, lon.size), transform=_tr, fill=0, dtype='uint8')).astype(bool)
        for _, r in _meso.iterrows()}
print('patches:', PATCHES, '| grid', lat.size, 'x', lon.size, '| SC cells', int(SC_MASK.sum()), '| mesoregions', list(MESO))

# %% [markdown]
# ## Cell 5 · Static inputs of the variants
#
# Growth forms, emission factors and canopy fractions as in the production run; the LAI the run used; the elevation
# fields of the lapse correction; the silviculture share; and the optional texture-based wilting point. Everything is
# saved in `OUT_ROOT/static_2019.npz` (notebook 23 reads the elevations and the land-cover shares).

# %%
# 5 · Static inputs -------------------------------------------------------------------------------------------------------------
gf_np = {k: v for k, v in np.load(os.path.join(ROOT_V10, str(YEAR), 'inputs', f'growthform_{YEAR}.npz')).items()}
ef = run.build_ef_grid_multiform(cfg, domain, gf_np)
VEG_FRAC = sum(gf_np[g] for g in landcover.GROWTH_FORMS)
CTF = {'base': landcover.canopy_type_fractions(gf_np, 3), 'tropical': landcover.canopy_type_fractions(gf_np, 2)}
def read_sf(p):
    with rasterio.open(p) as s:
        return np.flipud(s.read(1).astype('float64'))
LAI = {m: read_sf(os.path.join(ROOT_V10, str(YEAR), 'inputs', f'lai_{YEAR}_{m:02d}.tif')) for m in range(1, 13)}
LAI[0] = read_sf(os.path.join(ROOT_V10, str(YEAR - 1), 'inputs', f'lai_{YEAR - 1}_12.tif'))     # December of the year before
ef_pub = np.load(os.path.join(ROOT_V10, str(YEAR), 'inputs', f'ef_{YEAR}.npz'))['ef']
print(f'EF rebuilt vs ef_{YEAR}.npz of the run: max rel dev {np.nanmax(np.abs(ef - ef_pub) / np.maximum(np.abs(ef_pub), 1e-12)):.2e} (float32 storage)')

# silviculture (class 9) and forest formation (class 3) shares of each 0.025° cell, MapBiomas Collection 10, 2019, ~120 m
Z = np.load(REPO_DIR + '/data/external/mapbiomas_c10_lulc2019_SC_120m.npz')
def class_shares(Z, lon, lat, res, groups):
    """Share of the classified ~120 m pixels of each 0.025° cell (pixel centres) that fall in each group of class codes."""
    A, T6 = Z['a'], Z['transform']
    ix = np.floor((T6[2] + T6[0] * (np.arange(A.shape[1]) + 0.5) - (lon[0] - res / 2)) / res).astype(np.int64)
    iy = np.floor((T6[5] + T6[4] * (np.arange(A.shape[0]) + 0.5) - (lat[0] - res / 2)) / res).astype(np.int64)
    ok = (iy[:, None] >= 0) & (iy[:, None] < lat.size) & (ix[None, :] >= 0) & (ix[None, :] < lon.size) & (A > 0)
    flat_ = (iy[:, None] * lon.size + ix[None, :])[ok]; a_ok = A[ok]; n = lat.size * lon.size
    n_all = np.bincount(flat_, minlength=n).reshape(lat.size, lon.size)
    return {g: np.bincount(flat_[np.isin(a_ok, codes)], minlength=n).reshape(lat.size, lon.size) / np.maximum(n_all, 1)
            for g, codes in groups.items()}, n_all
SH, N120 = class_shares(Z, lon, lat, res, {'silv': [9], 'forest': [3]})
SILV, FOREST = SH['silv'], SH['forest']
needle = np.minimum(SILV, gf_np['tree'])
CTF['needle_silv'] = CTF['base'].copy(); CTF['needle_silv'][0] += needle * 100.0; CTF['needle_silv'][2] -= needle * 100.0
print(f"state shares: silviculture {100 * (SILV * area)[SC_MASK].sum() / area[SC_MASK].sum():.1f} %, forest formation "
      f"{100 * (FOREST * area)[SC_MASK].sum() / area[SC_MASK].sum():.1f} % (MapBiomas class areas: 10.5 % and 39.0 % of the classified area)")

# elevations: mean SRTM over each 0.025° cell, and over each ERA5-Land pixel sampled at the 0.025° cells
dem = ee.Image('CGIAR/SRTM90_V4').select('elevation')
era_proj = ee.ImageCollection(era5.COLLECTION).first().select('temperature_2m').projection()
z_fine = mapbiomas.sample_image_tiled(ee, dem.reduceResolution(ee.Reducer.mean(), maxPixels=4096).rename('z'), domain, ['z'], ntiles=4)['z']
z_era_dem = mapbiomas.sample_image_tiled(ee, dem.reduceResolution(ee.Reducer.mean(), maxPixels=65536).reproject(era_proj).rename('z'),
                                         domain, ['z'], ntiles=4)['z']
z_era = z_era_dem
if ZOROG_NC:
    import xarray as xr
    zz = xr.open_dataset(ZOROG_NC); zv = zz['z'].squeeze() / 9.80665
    la_n = [c for c in zv.dims if 'lat' in c][0]; lo_n = [c for c in zv.dims if 'lon' in c][0]
    LON2, LAT2 = np.meshgrid(lon, lat)
    zl = zv.sel({la_n: xr.DataArray(np.round(LAT2, 1), dims=('y', 'x')),
                 lo_n: xr.DataArray(np.round(np.mod(LON2, 360.0) if float(zv[lo_n].max()) > 180 else LON2, 1), dims=('y', 'x'))}, method='nearest')
    z_era = np.asarray(zl.values, dtype='float64')
    d_ = (z_era - z_era_dem)[SC_MASK & np.isfinite(z_era_dem)]
    print(f'ERA5-Land orography minus SRTM mean of the ERA5-Land pixel (state cells): mean {d_.mean():.1f} m, SD {d_.std():.1f} m')
if GAMMA_K_PER_KM is None:
    GAMMA_K_PER_KM = json.load(open(N21_JSON)).get('gamma_for_notebook20_K_per_km', 6.5) if os.path.exists(N21_JSON) else 6.5
GAMMA_SRC = 'notebook 21' if os.path.exists(N21_JSON) else 'standard atmosphere'
DZ = np.where(np.isfinite(z_fine) & np.isfinite(z_era), z_fine - z_era, 0.0)
DT = -GAMMA_K_PER_KM * DZ / 1000.0
q = np.percentile(DT[SC_MASK], [1, 5, 25, 50, 75, 95, 99])
print(f'lapse rate {GAMMA_K_PER_KM:.2f} K/km ({GAMMA_SRC}); state cells: dz mean {DZ[SC_MASK].mean():+.1f} m, SD {DZ[SC_MASK].std():.1f} m; '
      f'dT percentiles 1-5-25-50-75-95-99: {np.round(q, 2)} K')
WILT_TEX = None
if 'wilt_texture' in VARIANTS:
    from megan_sc.gee import soil_texture
    sc_ = mapbiomas.sample_image_tiled(ee, soil_texture.sand_clay_image(ee, region).reduceResolution(ee.Reducer.mean(), maxPixels=4096),
                                       domain, ['sand', 'clay'], ntiles=4)
    WILT_TEX = np.where(np.isfinite(sc_['sand']) & np.isfinite(sc_['clay']), soil.wilting_point_from_texture(sc_['sand'], sc_['clay']), WILT)
    print('texture wilting point, state cells: percentiles 5-50-95', np.round(np.percentile(WILT_TEX[SC_MASK], [5, 50, 95]), 3))
os.makedirs(OUT_ROOT, exist_ok=True)
np.savez_compressed(os.path.join(OUT_ROOT, f'static_{YEAR}.npz'), lat=lat, lon=lon, area=area, sc_mask=SC_MASK, z_fine=z_fine, z_era5=z_era,
                    z_era5_srtm=z_era_dem, dz=DZ, dT=DT, gamma_K_per_km=GAMMA_K_PER_KM, silviculture=SILV, forest_formation=FOREST,
                    **{f'gf_{k}': v for k, v in gf_np.items()}, wilt_texture=WILT_TEX if WILT_TEX is not None else np.nan)

# %% [markdown]
# ## Cell 6 · Samplers (as cell 7 of notebook 17) and the variant definitions

# %%
# 6 · Samplers and variants -------------------------------------------------------------------------------------------------------
def sample_multi(img, bands):
    grid = {'dimensions': {'width': int(lon.size), 'height': int(lat.size)},
            'affineTransform': {'scaleX': res, 'shearX': 0.0, 'translateX': float(lon[0] - res / 2),
                                'shearY': 0.0, 'scaleY': -res, 'translateY': float(lat[-1] + res / 2)},
            'crsCode': 'EPSG:4326'}
    for attempt in range(4):
        try:
            st_ = ee.data.computePixels({'expression': img.select(bands), 'grid': grid, 'fileFormat': 'NUMPY_NDARRAY'})
            return {b: np.flipud(np.asarray(st_[b], dtype='float64')) for b in bands}
        except Exception as e:
            if attempt == 3: raise
            print('   EE retry', attempt + 1, repr(e)[:90]); time.sleep(10 * (attempt + 1))
def daily_tmean_image(day_iso):
    d = ee.Date(day_iso)
    return (ee.ImageCollection(era5.COLLECTION).filterBounds(region).filterDate(d, d.advance(1, 'day'))
            .select('temperature_2m').mean().rename('D_TEMP'))
CACHE = os.path.join(OUT_ROOT, 'met_cache'); os.makedirs(CACHE, exist_ok=True)
def day_met(m, d):
    """The 8 three-hourly samples and the 24-h mean temperature of one day (Earth Engine, or the cache)."""
    f = os.path.join(CACHE, f'met_{YEAR}_{m:02d}_{d:02d}.npz')
    if os.path.exists(f):
        z = np.load(f)
        return {h: {b: z[f'{b}_{h:02d}'].astype('float64') for b in MET_BANDS} for h in HOURS}, z['D_TEMP'].astype('float64')
    day_iso = f'{YEAR}-{m:02d}-{d:02d}'
    mets = {h: sample_multi(era5.megan_met_image(ee, region, ee.Date(f'{day_iso}T{h:02d}:00')), MET_BANDS) for h in HOURS}
    dtm = sample_multi(daily_tmean_image(day_iso), ['D_TEMP'])['D_TEMP']
    if MET_CACHE:
        np.savez_compressed(f, D_TEMP=dtm.astype('float32'), **{f'{b}_{h:02d}': mets[h][b].astype('float32') for h in HOURS for b in MET_BANDS})
    return mets, dtm
KF = {'temp': 'TEMP2', 'par': 'PAR', 'wind': 'WINDSPD', 'pres': 'PRES', 'qv': 'QV', 'soilm': 'SOILM'}
def variant_inputs(v, metf, d_temp, d_ppfd):
    """(met, d_temp, d_ppfd, canopy fractions, wilting point) of variant v from the filled base fields."""
    met, ctf, wilt = dict(metf), CTF['base'], WILT
    if v == 'lapse':
        met['temp'] = metf['temp'] + DT; d_temp = d_temp + DT
    elif v == 'par050':
        met['par'] = metf['par'] * (PAR_FRAC_ALT / PAR_FRAC_RUN); d_ppfd = d_ppfd * (PAR_FRAC_ALT / PAR_FRAC_RUN)
    elif v == 'canopy_tropical':
        ctf = CTF['tropical']
    elif v == 'canopy_needle_silv':
        ctf = CTF['needle_silv']
    elif v == 'wilt_clay':
        wilt = WILT_ALT
    elif v == 'wilt_texture':
        wilt = WILT_TEX
    elif v != 'base':
        raise ValueError(v)
    return met, d_temp, d_ppfd, ctf, wilt
print('variants:', VARIANTS)

# %% [markdown]
# ## Cell 7 · Run
#
# One checkpoint per month (after each finished day). A finished month is written as `OUT/<variant>/class_tonnes_<year>_<mm>.npz`
# (19 classes × grid, tonnes over the days run in that month) and `diurnal_<year>_<mm>.csv` (state tonnes per class and
# hour). Safe to run again.

# %%
# 7 · Run --------------------------------------------------------------------------------------------------------------------------------
def run_month(m):
    days = [d for (mm, d) in DAYS if mm == m]
    if not days:
        return
    outs = {v: os.path.join(OUT, v, f'class_tonnes_{YEAR}_{m:02d}.npz') for v in VARIANTS}
    if all(os.path.exists(p) for p in outs.values()):
        print(f'  {YEAR}-{m:02d}: done'); return
    ck = os.path.join(OUT, f'checkpoint_{YEAR}_{m:02d}.npz')
    ny, nx = lat.size, lon.size
    loaded = st.load_checkpoint(ck)
    if loaded is not None and set(k.split('__')[0] for k in loaded[0]) == set(VARIANTS):
        acc, done = loaded; print(f'  {YEAR}-{m:02d}: resuming, {len(done)} of {len(days)} days done')
    else:
        acc = {**{f'{v}__t': np.zeros((len(CLS), ny, nx)) for v in VARIANTS}, **{f'{v}__di': np.zeros((len(CLS), len(HOURS))) for v in VARIANTS}}
        done = set()
    lai_c, lai_p = LAI[m], LAI[m - 1]
    for d in days:
        if d in done:
            continue
        t0 = time.time()
        doy = datetime.date(YEAR, m, d).timetuple().tm_yday; date = YEAR * 1000 + doy
        mets, dtm = day_met(m, d)
        D_TEMP = fixes.fill_field_nearest(dtm, *fixes.VALID_RANGE['temp'])
        filled = {}
        for h in HOURS:
            mf, _ = fixes.fill_met_nearest({k: mets[h][b] for k, b in KF.items()})          # FIX-1 before any change
            filled[h] = mf
        D_PPFD = np.mean([filled[h]['par'] for h in HOURS], axis=0) * fixes.PAR_WM2_TO_PPFD   # FIX-2
        for hi, h in enumerate(HOURS):
            for v in VARIANTS:
                met, dt_, dp_, ctf, wilt = variant_inputs(v, filled[h], D_TEMP, D_PPFD)
                met = {**met, 'lai': lai_c}
                r = fixes.run_timestep_fixed(domain, ef, met, day=doy, hour_utc=h, canopy_type_fractions=ctf, d_temp=dt_,
                                             d_ppfd=dp_, laip=lai_p, nlai=12, wilt=wilt, date_yyyyddd=date, strict=True)
                t = spc.compounds_to_classes(spc.compound_tonnes(spc.flux_to_mol_s(r['compounds'], area[None, :, :]), SEC_PER_STEP))
                acc[f'{v}__t'] += t; acc[f'{v}__di'][:, hi] += (t * SC_MASK).sum(axis=(1, 2))
        done.add(d); st.save_checkpoint(ck, acc, done)
        print(f'    {YEAR}-{m:02d}-{d:02d} ({time.time() - t0:.0f} s)', flush=True)
    for v in VARIANTS:
        os.makedirs(os.path.join(OUT, v), exist_ok=True)
        np.savez_compressed(outs[v], tonnes=acc[f'{v}__t'].astype('float32'), days=np.array(days), classes=np.array(CLS),
                            lat=lat, lon=lon, variant=v)
        pd.DataFrame(acc[f'{v}__di'], index=CLS, columns=[f'{h:02d}h' for h in HOURS]).to_csv(os.path.join(OUT, v, f'diurnal_{YEAR}_{m:02d}.csv'))
    os.remove(ck); print(f'  {YEAR}-{m:02d}: written for {len(VARIANTS)} variants')
json.dump(dict(notebook='20_sensitivity_2019_inputs', year=YEAR, mode=MODE, days=DAYS, variants=VARIANTS, fixes_sha256=FIXES_SHA,
               patches=PATCHES, gamma_K_per_km=GAMMA_K_PER_KM, gamma_source=GAMMA_SRC, zorog=ZOROG_NC or 'SRTM mean over the ERA5-Land pixel',
               par_frac=[PAR_FRAC_RUN, PAR_FRAC_ALT], wilt=[WILT, WILT_ALT], ee_project='(set in cell 3)',
               started=datetime.datetime.utcnow().isoformat() + 'Z'), open(os.path.join(OUT, 'provenance_notebook20.json'), 'w'), indent=1)
T_START = time.time()
for m in range(1, 13):
    run_month(m)
print(f'finished in {(time.time() - T_START) / 60:.0f} min')

# %% [markdown]
# ## Cell 8 · Checks and summary
#
# * `base` against the production run: the state totals of every class over the same days must agree to better than
#   1e-4 (relative; the LAI is read from its float32 GeoTIFF). In `design` mode the reference is the 3-hourly S1 test run
#   (`emis_2019_<mm>_test7d.nc`), in `full` mode the monthly files of the production run.
# * Each variant against `base`: state totals by group, isoprene by mesoregion, and, for `lapse`, the cell-level
#   response of isoprene to ΔT (slope of Δln E on ΔT, an empirical temperature sensitivity used by notebook 23).
# * If notebook 18 has run, the ratio of this inventory to BRAIN (true areas) implied for each variant:
#   ratio(base, notebook 18) × variant / base. In `design` mode this is indicative only (28 days).

# %%
# 8 · Checks and summary ----------------------------------------------------------------------------------------------------------------
from netCDF4 import Dataset
MONTHS_RUN = sorted({m for m, _ in DAYS})
V = {v: {m: np.load(os.path.join(OUT, v, f'class_tonnes_{YEAR}_{m:02d}.npz'))['tonnes'].astype('float64') for m in MONTHS_RUN} for v in VARIANTS}
ANN = {v: sum(V[v].values()) for v in VARIANTS}
ref_rows = []
for m in MONTHS_RUN:
    f = os.path.join(ROOT_V10, str(YEAR), 'emissions', f'emis_{YEAR}_{m:02d}' + ('_test7d' if MODE == 'design' else '') + '.nc')
    if not os.path.exists(f):
        print('reference file not found:', f); continue
    with Dataset(f) as ds:
        R = np.stack([np.nan_to_num(np.asarray(ds.variables[c][:], 'float64')) for c in CLS])
    a, b = (V['base'][m] * SC_MASK).sum((1, 2)), (R * SC_MASK).sum((1, 2))
    ref_rows += [dict(month=m, cls=c, base=a[i], production=b[i], rel_diff=a[i] / b[i] - 1 if b[i] else np.nan) for i, c in enumerate(CLS)]
REF = pd.DataFrame(ref_rows); REF.to_csv(os.path.join(OUT, 'check_base_vs_production.csv'), index=False)
if len(REF):
    print(f'base vs production ({"S1 test run" if MODE == "design" else "monthly files"}): max |rel diff| {REF.rel_diff.abs().max():.2e} '
          f'({"OK" if REF.rel_diff.abs().max() < 1e-4 else "LARGER THAN 1e-4: check before using the variants"})')
GROUPS = {'BVOC': [c for c in CLS if c not in ('NO', 'CO')], 'ISOP': ['ISOP'], 'MT': [c for c in CLS if c.startswith('MT_')],
          'SQT': [c for c in CLS if c.startswith('SQT')], 'OVOC': ['MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD'], 'LVOC': ['LVOC'],
          'MBO': ['MBO'], 'STRESS_OTHER': ['STRESS', 'OTHER'], 'NO': ['NO'], 'CO': ['CO']}
IX = {c: i for i, c in enumerate(CLS)}
gsum = lambda a, cl, msk=SC_MASK: float(sum((a[IX[c]] * msk).sum() for c in cl))
rows = []
for v in VARIANTS:
    for g, cl in GROUPS.items():
        rows.append(dict(variant=v, region='Santa Catarina', group=g, tonnes=gsum(ANN[v], cl), change_vs_base_pct=100 * (gsum(ANN[v], cl) / gsum(ANN['base'], cl) - 1)))
    for nm, msk in MESO.items():
        rows.append(dict(variant=v, region=nm, group='ISOP', tonnes=gsum(ANN[v], ['ISOP'], msk),
                         change_vs_base_pct=100 * (gsum(ANN[v], ['ISOP'], msk) / gsum(ANN['base'], ['ISOP'], msk) - 1)))
SUM = pd.DataFrame(rows); SUM.to_csv(os.path.join(OUT, 'variants_summary.csv'), index=False)
pd.set_option('display.width', 200)
print(SUM[SUM.region == 'Santa Catarina'].pivot(index='group', columns='variant', values='change_vs_base_pct').round(2).to_string())
print(SUM[SUM.group == 'ISOP'].pivot(index='region', columns='variant', values='change_vs_base_pct').round(2).to_string())
if 'lapse' in VARIANTS:
    e0, e1 = ANN['base'][IX['ISOP']], ANN['lapse'][IX['ISOP']]
    sel = SC_MASK & (e0 > 0) & (e1 > 0) & (np.abs(DT) > 0.05)
    y_, x_ = np.log(e1[sel] / e0[sel]), DT[sel]
    slope = float(np.sum(x_ * y_) / np.sum(x_ * x_))                      # through the origin: no change where dT = 0
    print(f'lapse: isoprene d ln E / d T over {sel.sum()} state cells = {slope:.4f} per K (through the origin); '
          f'state change {100 * (e1[SC_MASK].sum() / e0[SC_MASK].sum() - 1):+.2f} %')
    json.dump(dict(isop_dlnE_dT_per_K=slope, cells=int(sel.sum()), gamma_K_per_km=GAMMA_K_PER_KM, mode=MODE),
              open(os.path.join(OUT, 'lapse_temperature_sensitivity.json'), 'w'), indent=1)
B18 = os.path.join(ROOT_V10, 'analysis_v10', 'brain_comparison', f'comparison_{YEAR}_classes_t.csv')
if os.path.exists(B18):
    tb = pd.read_csv(B18).set_index('item')
    MAP18 = {'ISOP': ['ISOP'], 'MT (as TERP)': GROUPS['MT'], 'SQT (as SQT)': GROUPS['SQT'], 'MEOH': ['MEOH'], 'ETOH': ['ETOH'], 'NO': ['NO'], 'CO': ['CO']}
    imp = []
    for it, cl in MAP18.items():
        if it in tb.index:
            for v in VARIANTS:
                imp.append(dict(item=it, variant=v, ratio_vs_brain_true=tb.loc[it, 'ratio_vs_true'] * gsum(ANN[v], cl) / gsum(ANN['base'], cl)))
    IMP = pd.DataFrame(imp).pivot(index='item', columns='variant', values='ratio_vs_brain_true')
    IMP.to_csv(os.path.join(OUT, 'implied_ratio_vs_brain.csv'))
    print('implied ratio this inventory / BRAIN (true areas)' + (' - indicative, 28 days' if MODE == 'design' else '') + ':')
    print(IMP.round(3).to_string())
print('outputs in', OUT)
