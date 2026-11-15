# %% [markdown]
# # 14 · Observational check: TROPOMI formaldehyde (HCHO) against the modelled isoprene
#
# Independent evaluation of the isoprene backbone of the inventory with Sentinel-5P TROPOMI tropospheric HCHO
# columns (Earth Engine `COPERNICUS/S5P/OFFL/L3_HCHO`, from 2018-12-05; Earth Engine removes pixels with qa_value < 0.5).
# Four tests (Tiers 1-4, fixed in this notebook) for 2019-2025, the full TROPOMI years.
#
# Physical basis (Palmer et al., 2003, eq. 1): without transport the HCHO column is Ω = (1/k_HCHO) Σ Y_i E_i, where E_i is
# the emission in atoms C cm⁻² s⁻¹ and Y_i the HCHO yield per C atom. Isoprene lives about one hour and HCHO a few hours,
# so the column follows isoprene emission with some downwind smearing. Palmer et al. give ultimate isoprene yields of
# 0.3-0.45 per C in their mechanism (0.43 per C at 1 ppb NOx, 0.29 per C at 0.1 ppb NOx, where only 55 % of it forms within
# 2 h). HCHO is therefore an isoprene-weighted test: it cannot check monoterpenes, sesquiterpenes or soil NO separately.
#
# | Tier | Question | Estimator |
# |---|---|---|
# | 1 spatial | does HCHO sit where the model puts isoprene? | Pearson/Spearman r of 2019-2025 summer (DJF) means on 0.1° and 0.2° blocks, after smoothing the model with Gaussian kernels of 0-40 km (transport smearing); 95 % CI by spatial block bootstrap (0.5° blocks) |
# | 2 seasonal | same annual cycle? | 12-month climatologies of SC-mean ΔHCHO and modelled isoprene: r, peak month |
# | 3 interannual | same anomalies? | deseasonalised monthly anomalies (n = 84): r with the effective sample size N(1 − r₁r₂)/(1 + r₁r₂) of Bretherton et al. (1999); the same against SC-mean 2 m temperature; annual means (n = 7) are shown only descriptively |
# | 4 magnitude | is the column-emission slope physically plausible? | slope S of ΔΩ (molec cm⁻²) on E (atoms C cm⁻² s⁻¹) at the overpass time, summer, 0.2° blocks; S = Y/k_HCHO, so S/Y is an implied HCHO lifetime to compare with "a few hours" |
#
# Choices that the reader must see (each has a sensitivity switch below): orbit-level cloud fraction < 0.3 and solar zenith
# < 70°; fire screening that removes HCHO pixels within 20 km of a FIRMS detection on the same or previous UTC day (MODIS active
# fires; FIRMS is not "science quality" according to its catalogue page); background = monthly median over an Atlantic box
# east of SC (45° to 41° W, 30° to 25° S); its anomalies share a regional signal with modelled isoprene (additional diagnostic); TROPOMI bias: Vigouroux et
# al. (2020) found TROPOMI = 1.10×10¹⁵ + 0.64 × FTIR (overestimation of 26 % at very low columns, underestimation of 31 %
# above 8×10¹⁵ molec cm⁻²); `BIAS_CORRECT = True` inverts that regression as a sensitivity case.
#
# The modelled isoprene is the output of notebook 10 times the canopy loss and production factor rho = 0.96 of Guenther et
# al. (2006, Sect. 3.3) (RHO_ISOP, cell 1), as in notebook 13.
#
# A result here is evidence, not proof: a mismatch can come from the retrieval (noise, air-mass factor, clouds), the
# NOx-dependent yield, transport, or fires as much as from the inventory.
#
# Earth Engine side: one batch export per month to an asset folder of your project (`megan_sc_hcho_<TAG>`), read back
# on the model grid and cached; re-running cell 3 resumes. The tasks use Earth Engine compute (EECU) of your project; cell 3
# prints the EECU-hours per month as they finish.
#
# **Send back:** the folder `OUT_DIR`.

# %%
# 1 · Settings --------------------------------------------------------------------------------------------
import os, glob, json, datetime, calendar, warnings
try:
    from google.colab import drive
    IN_COLAB = True
except ImportError:
    IN_COLAB = False
if IN_COLAB:
    # A failed mount must stop the notebook: if it is hidden, the notebook reads and writes a local /content/drive
    # folder that is not Google Drive (notebook 13 then stops with "found 0").
    if os.path.isdir('/content/drive') and not os.path.ismount('/content/drive') and os.listdir('/content/drive'):
        raise RuntimeError('/content/drive holds local files written while Google Drive was not mounted, so Drive cannot be '
                           'mounted there. Use Runtime > Disconnect and delete runtime, then run the notebook again from cell 1.')
    drive.mount('/content/drive')
%pip -q install earthengine-api rasterio geopandas netCDF4 scipy pandas matplotlib
EE_PROJECT = 'your-ee-project'                                     # <-- your Earth Engine (Google Cloud) project id
MODEL_ROOT = '/content/drive/MyDrive/megan_sc_outputs_v10'          # <-- production run (notebook 10)
BND        = '/content/drive/MyDrive/megan-sc/data/boundaries'
OUT_DIR    = os.path.join(MODEL_ROOT, 'analysis_v10', 'tropomi_hcho')
YEARS      = list(range(2019, 2026))
CF_MAX, SZA_MAX = 0.3, 70.0              # orbit-level screening (sensitivity: CF_MAX = 0.5)
FIRE_MASK, FIRE_BUF_KM, FIRE_DAYS = True, 20, 1                    # sensitivity: FIRE_MASK = False
BG_BOX     = [-45.0, -30.0, -41.0, -25.0]                          # lon_min, lat_min, lon_max, lat_max (Atlantic, east of the state)
BIAS_CORRECT = False                     # True: invert TROPOMI = 1.10e15 + 0.64 FTIR (Vigouroux et al., 2020)
OVERPASS_UTC = 16.5                      # TROPOMI overpass 13:30 LT; 16:30 UTC = 12:55-13:17 local mean solar time in SC
RHO_ISOP = 0.96                          # isoprene canopy loss and production factor (Guenther et al., 2006); notebook 10 writes rho = 1
TAG = f'cf{CF_MAX:g}_sza{SZA_MAX:g}_fire{int(FIRE_MASK)}{FIRE_BUF_KM}'
print('notebook 14 (Earth Engine batch exports) | TAG', TAG)

# %%
# 2 · Model side: grid, SC mask, monthly isoprene, diurnal shape, monthly temperature -------------------------
os.makedirs(os.path.join(OUT_DIR, 'cache_' + TAG), exist_ok=True)
import numpy as np, pandas as pd, rasterio, geopandas as gpd, matplotlib.pyplot as plt
from netCDF4 import Dataset
from rasterio.features import rasterize
from rasterio.transform import from_origin
from scipy import stats, ndimage
warnings.filterwarnings('ignore', category=RuntimeWarning)
with Dataset(sorted(glob.glob(os.path.join(MODEL_ROOT, str(YEARS[0]), 'emissions', f'emis_{YEARS[0]}_01.nc')))[0]) as ds:
    lat = np.array(ds.variables['lat'][:]); lon = np.array(ds.variables['lon'][:])
NY, NX = lat.size, lon.size; RES = float(lon[1] - lon[0])
TR = from_origin(lon[0] - RES / 2, lat[-1] + RES / 2, RES, RES)
g = gpd.read_file(os.path.join(BND, 'sc_state_ibge.gpkg')).to_crs('EPSG:4326')
SC_GEOM = g.geometry.union_all() if hasattr(g.geometry, 'union_all') else g.geometry.unary_union
SC = np.flipud(rasterize([(SC_GEOM, 1)], out_shape=(NY, NX), transform=TR, fill=0, dtype='uint8')).astype(bool)
A_M2 = (6371008.8 ** 2 * np.deg2rad(RES) * np.abs(np.sin(np.deg2rad(lat + RES / 2)) - np.sin(np.deg2rad(lat - RES / 2))))[:, None] * np.ones((1, NX))
NA = 6.02214076e23
MONTHS = [(y, m) for y in YEARS for m in range(1, 13)]
def model_month(y, m):
    """Isoprene for one month: per-cell mean flux at the overpass (atoms C cm-2 s-1) and SC total (t)."""
    with Dataset(os.path.join(MODEL_ROOT, str(y), 'emissions', f'emis_{y}_{m:02d}.nc')) as ds:
        t = np.nan_to_num(np.array(ds.variables['ISOP'][:], dtype='float64')) * RHO_ISOP   # t/month, south-first
    sec = calendar.monthrange(y, m)[1] * 86400.0
    mean_flux = t * 1e6 / 68.117 * 5 * NA / (A_M2 * 1e4) / sec                          # atoms C cm-2 s-1, monthly mean
    di = pd.read_csv(os.path.join(MODEL_ROOT, str(y), 'hourly', f'diurnal_{y}_{m:02d}.csv'), index_col=0).loc['ISOP']
    hrs = np.array([int(c[:2]) for c in di.index]); v = di.values.astype(float)
    shape_at_overpass = np.interp(OVERPASS_UTC, np.r_[hrs, 24], np.r_[v, v[0]]) / v.mean()   # SC diurnal shape, 3-hourly samples
    return mean_flux * shape_at_overpass, float((t * SC).sum())
def sc_mean_tif(path, band):
    with rasterio.open(path) as s: b = np.flipud(s.read(band).astype('float64'))
    return float(np.nanmean(np.where(SC, b, np.nan)))
MOD = {}; ROWS = []
for y, m in MONTHS:
    E, tot = model_month(y, m); MOD[(y, m)] = E
    ROWS.append(dict(year=y, month=m, ISOP_t=tot, T2_C=sc_mean_tif(os.path.join(MODEL_ROOT, str(y), 'inputs', f'met_{y}_{m:02d}.tif'), 1) - 273.15))
SER = pd.DataFrame(ROWS)
print('model months:', len(MOD), '| SC isoprene 2019-2025 mean (t/month):', round(SER.ISOP_t.mean()))

# %%
# 3 · Satellite side (Earth Engine): monthly screened HCHO on the model grid, background, fire statistics --------
# Earth Engine work runs as one batch export per month into an asset folder of your project. An interactive request per
# month times out ("Computation timed out"), because the 20 km fire buffer is built for every orbit on the native 1.1 km
# TROPOMI L3 grid. Batch tasks have no such limit and keep running if Colab disconnects:
# re-run this cell to resume. Each finished asset (5 bands on the model grid) is read once into OUT_DIR/cache_<TAG>;
# everything after this cell uses only that cache. Tasks also show in the Tasks tab of the Earth Engine Code Editor.
CACHE = os.path.join(OUT_DIR, 'cache_' + TAG)
SAFE = TAG.replace('.', 'p')                                      # asset IDs allow letters, digits, '_' and '-'
ASSET_DIR = f'projects/{EE_PROJECT}/assets/megan_sc_hcho_{SAFE}'
WAIT_HOURS, POLL_S = 8, 60
def cache_file(y, m): return os.path.join(CACHE, f'hcho_{y}_{m:02d}.npz')
def desc(y, m): return f'hcho_{SAFE}_{y}_{m:02d}'
todo = [(y, m) for (y, m) in MONTHS if not os.path.exists(cache_file(y, m))]
print(f'months already cached: {len(MONTHS) - len(todo)} of {len(MONTHS)}')
if todo:
    import ee, time
    try:
        ee.Initialize(project=EE_PROJECT)
    except Exception:
        ee.Authenticate(); ee.Initialize(project=EE_PROJECT)
    S5P = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_HCHO'); FIRMS = ee.ImageCollection('FIRMS')
    NATIVE = ee.Projection('EPSG:4326').atScale(1113.2)
    REGION = ee.Geometry.Rectangle([float(lon[0] - RES), float(lat[0] - RES), float(lon[-1] + RES), float(lat[-1] + RES)])
    BG = ee.Geometry.Rectangle(BG_BOX)
    X0, Y0 = float(lon[0] - RES / 2), float(lat[-1] + RES / 2)
    GRID = {'dimensions': {'width': NX, 'height': NY}, 'crsCode': 'EPSG:4326',
            'affineTransform': {'scaleX': RES, 'shearX': 0, 'translateX': X0, 'shearY': 0, 'scaleY': -RES, 'translateY': Y0}}
    def screened(img):
        d = ee.Date(img.get('system:time_start'))
        ok = img.select('cloud_fraction').lt(CF_MAX).And(img.select('solar_zenith_angle').lt(SZA_MAX))
        if FIRE_MASK:
            day = ee.Date.fromYMD(d.get('year'), d.get('month'), d.get('day'))      # UTC calendar day of the orbit
            fc = FIRMS.filterDate(day.advance(-FIRE_DAYS, 'day'), day.advance(1, 'day'))   # same day and FIRE_DAYS before
            fire = ee.Image(ee.Algorithms.If(fc.size().gt(0), fc.select('T21').max().unmask(0).gt(0), ee.Image(0)))
            fire = fire.focal_max(radius=FIRE_BUF_KM * 1000, kernelType='circle', units='meters')
            ok = ok.And(fire.Not())
        h = img.select('tropospheric_HCHO_column_number_density')
        return h.updateMask(ok).rename('HCHO').addBands(h.mask().And(ok.Not()).rename('removed'))
    def month_image(y, m):
        """Bands on the model grid: HCHO (mean, mol m-2), N (mean count of kept observations per native pixel), RM (same for
        observations removed by the screening), BG (median of the monthly mean over the Atlantic box, constant), ORB (orbits)."""
        d0 = ee.Date.fromYMD(y, m, 1); d1 = d0.advance(1, 'month')
        col = S5P.filterDate(d0, d1).filterBounds(REGION.union(BG, 1)).map(screened)
        hm = col.select('HCHO').mean()
        grid = (hm.addBands(col.select('HCHO').count().rename('N')).addBands(col.select('removed').sum().rename('RM'))
                  .setDefaultProjection(NATIVE).reduceResolution(reducer=ee.Reducer.mean(), maxPixels=256).unmask(-9999))
        bg = hm.reduceRegion(ee.Reducer.median(), BG, 5566, maxPixels=1e9).get('HCHO')
        return (grid.addBands(ee.Image.constant(ee.Number(ee.Algorithms.If(bg, bg, -9999))).rename('BG'))
                    .addBands(ee.Image.constant(col.size()).rename('ORB')).toFloat())
    def start_export(y, m):
        t = ee.batch.Export.image.toAsset(image=month_image(y, m), description=desc(y, m), assetId=f'{ASSET_DIR}/hcho_{y}_{m:02d}',
                                          crs='EPSG:4326', crsTransform=[RES, 0, X0, 0, -RES, Y0], dimensions=f'{NX}x{NY}',
                                          maxPixels=int(1e9))
        t.start(); return t
    def assets_done():
        try:
            return {a['name'].split('/')[-1] for a in ee.data.listAssets({'parent': ASSET_DIR}).get('assets', [])}
        except Exception:
            return set()
    def operations():
        """Latest task per description of this TAG: state (PENDING, RUNNING, SUCCEEDED, FAILED, CANCELLED), error, EECU-s."""
        out = {}
        for o in ee.data.listOperations():
            md = o.get('metadata', {}); k = md.get('description', '')
            if k.startswith(f'hcho_{SAFE}_') and (k not in out or md.get('updateTime', '') > out[k][3]):
                out[k] = (md.get('state'), (o.get('error') or {}).get('message'), float(md.get('batchEecuUsageSeconds') or 0), md.get('updateTime', ''))
        return out
    def read_asset(y, m):
        a = ee.data.computePixels({'expression': ee.Image(f'{ASSET_DIR}/hcho_{y}_{m:02d}'), 'fileFormat': 'NUMPY_NDARRAY', 'grid': GRID})
        H = np.flipud(np.where(a['HCHO'] == -9999, np.nan, a['HCHO']).astype('float64')) * NA / 1e4       # molec cm-2
        N = np.flipud(np.where(a['N'] == -9999, 0, a['N']).astype('float64'))
        RM = np.flipud(np.where(a['RM'] == -9999, 0, a['RM']).astype('float64'))
        bg = float(a['BG'][0, 0])
        # removed / kept: sums over the grid of the per-cell mean counts; only their ratio is used below
        np.savez_compressed(cache_file(y, m), H=H.astype('float32'), N=N.astype('float32'), bg=np.array(np.nan if bg == -9999 else bg * NA / 1e4),
                            removed=np.array(RM.sum()), kept=np.array(N.sum()), orbits=np.array(int(round(float(a['ORB'][0, 0])))))
    try:
        ee.data.getAsset(ASSET_DIR)
    except Exception:
        ee.data.createAsset({'type': 'FOLDER'}, ASSET_DIR)
    tasks, T0, last = {}, time.time(), -1e9                      # tasks: the exports started by this run of the cell
    def state_of(k, ops):
        if k in tasks:
            s_ = tasks[k].status(); return {'READY': 'PENDING', 'COMPLETED': 'SUCCEEDED'}.get(s_['state'], s_['state']), s_.get('error_message')
        v = ops.get(desc(*k)); return (v[0], v[1]) if v else (None, None)
    while True:
        ops = operations(); done = assets_done()         # in this order: a task reported SUCCEEDED has its asset listed
        for k in todo:
            if f'hcho_{k[0]}_{k[1]:02d}' in done and not os.path.exists(cache_file(*k)):
                read_asset(*k)
        left = [k for k in todo if not os.path.exists(cache_file(*k))]
        if not left:
            break
        st = {k: state_of(k, ops) for k in left}
        for k in left:                     # start what is neither finished nor queued/running elsewhere (once per run)
            if k not in tasks and st[k][0] not in ('PENDING', 'RUNNING') and f'hcho_{k[0]}_{k[1]:02d}' not in done:
                tasks[k] = start_export(*k); st[k] = ('PENDING', None)
        active = [k for k in left if st[k][0] in ('PENDING', 'RUNNING', 'UNSUBMITTED')]
        failed = [(k, st[k][1]) for k in left if k in tasks and st[k][0] in ('FAILED', 'CANCELLED')]
        if failed and not active:
            raise RuntimeError(f'{len(failed)} export task(s) ended without output, e.g. {failed[0][0]}: {failed[0][1]}. '
                               'Re-run this cell to retry them.')
        if time.time() - last > 300:
            used = [v[2] for v in ops.values() if v[0] == 'SUCCEEDED' and v[2] > 0]
            print(f"{time.strftime('%H:%M')}  cached {len(MONTHS) - len(left)}/{len(MONTHS)} | queued or running {len(active)} | "
                  f"failed {len(failed)} | EECU-h per finished month {np.mean(used) / 3600 if used else float('nan'):.2f}", flush=True)
            last = time.time()
        if time.time() - T0 > WAIT_HOURS * 3600:
            raise TimeoutError(f'still waiting after {WAIT_HOURS} h; the tasks keep running on Earth Engine. Re-run this cell later.')
        time.sleep(POLL_S)
    used = [v[2] for v in operations().values() if v[0] == 'SUCCEEDED']
    print(f'all {len(MONTHS)} months cached | batch compute of this TAG: {sum(used) / 3600:.1f} EECU-hours')
SAT = {(y, m): {k: v for k, v in np.load(cache_file(y, m)).items()} for (y, m) in MONTHS}
def corrected(H):
    return (H - 1.10e15) / 0.64 if BIAS_CORRECT else H
for i, (y, m) in enumerate(MONTHS):
    s = SAT[(y, m)]; H = corrected(s['H']); bg = float(corrected(np.array(s['bg'])))
    SER.loc[i, 'HCHO_SC'] = float(np.nanmean(np.where(SC, H, np.nan))); SER.loc[i, 'BG'] = bg
    SER.loc[i, 'dHCHO_SC'] = SER.loc[i, 'HCHO_SC'] - bg
    SER.loc[i, 'frac_removed_screening'] = float(s['removed'] / max(s['removed'] + s['kept'], 1))
    SER.loc[i, 'orbits'] = int(s['orbits'])
SER.to_csv(os.path.join(OUT_DIR, f'monthly_series_{TAG}{"_biascorr" if BIAS_CORRECT else ""}.csv'), index=False)
print(SER.head(13).round(3).to_string())

# %%
# 4 · Tier 1, spatial pattern (summer DJF, 2019-2025) -----------------------------------------------------------
def block(a, k, w=None):
    ny, nx = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    a = a[:ny, :nx]; w = np.ones_like(a) if w is None else w[:ny, :nx]
    num = np.nansum((a * w).reshape(ny // k, k, nx // k, k), axis=(1, 3)); den = np.nansum((np.isfinite(a) * w).reshape(ny // k, k, nx // k, k), axis=(1, 3))
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
DJF = [(y, m) for (y, m) in MONTHS if m in (12, 1, 2)]
def dmean(keys, field):
    return np.nanmean(np.stack([field(k) for k in keys]), axis=0)
H_djf = dmean(DJF, lambda k: corrected(SAT[k]['H']) - float(corrected(np.array(SAT[k]['bg']))))
E_djf = dmean(DJF, lambda k: MOD[k])
scf = SC.astype(float)
def block_boot_r(x, y, k_blocks, n=2000, seed=1):
    """95 % CI of Pearson r by resampling spatial blocks of k_blocks x k_blocks cells (keeps spatial autocorrelation inside blocks)."""
    rng = np.random.default_rng(seed); ok = np.isfinite(x) & np.isfinite(y)
    by, bx = np.indices(x.shape); bid = ((by // k_blocks) * 1000 + (bx // k_blocks))[ok]; xv, yv = x[ok], y[ok]
    ids = np.unique(bid); members = [np.where(bid == b)[0] for b in ids]; rs = []
    for _ in range(n):
        sel = np.concatenate([members[i] for i in rng.integers(0, ids.size, ids.size)])
        if sel.size > 3 and np.std(xv[sel]) > 0 and np.std(yv[sel]) > 0: rs.append(np.corrcoef(xv[sel], yv[sel])[0, 1])
    return np.percentile(rs, [2.5, 97.5])
T1 = []
for k, kb in ((4, 5), (8, 3)):                                   # 0.1° blocks (bootstrap blocks 0.5°), 0.2° blocks (~0.6°)
    frac = block(scf, k)
    for sig_km in (0, 10, 20, 40):
        Es = ndimage.gaussian_filter(np.nan_to_num(E_djf), sig_km / (RES * 111.0)) if sig_km else E_djf
        x = np.where(frac >= 0.5, block(np.where(SC, Es, np.nan), k), np.nan); y_ = np.where(frac >= 0.5, block(np.where(SC, H_djf, np.nan), k), np.nan)
        ok = np.isfinite(x) & np.isfinite(y_)
        r = np.corrcoef(x[ok], y_[ok])[0, 1]; rho = stats.spearmanr(x[ok], y_[ok])[0]; lo, hi = block_boot_r(x, y_, kb)
        T1.append(dict(block_deg=k * RES, smoothing_km=sig_km, n_blocks=int(ok.sum()), pearson_r=r, r_lo95=lo, r_hi95=hi, spearman_rho=rho))
T1 = pd.DataFrame(T1); T1.to_csv(os.path.join(OUT_DIR, f'tier1_spatial_{TAG}.csv'), index=False); print(T1.round(3).to_string(index=False))

# %%
# 5 · Tiers 2 and 3, seasonal cycle and interannual anomalies -----------------------------------------------
clim = SER.groupby('month')[['ISOP_t', 'dHCHO_SC', 'T2_C']].mean()
r2, p2 = stats.pearsonr(clim.ISOP_t, clim.dHCHO_SC)
T2 = dict(r_climatology=r2, p=p2, peak_month_model=int(clim.ISOP_t.idxmax()), peak_month_hcho=int(clim.dHCHO_SC.idxmax()),
          min_month_model=int(clim.ISOP_t.idxmin()), min_month_hcho=int(clim.dHCHO_SC.idxmin()))
an = SER.copy()
for c in ('ISOP_t', 'dHCHO_SC', 'T2_C'):
    an[c] = SER[c] - SER.groupby('month')[c].transform('mean')
def r_neff(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float); ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
    r = np.corrcoef(a, b)[0, 1]; r1 = np.corrcoef(a[:-1], a[1:])[0, 1]; r2_ = np.corrcoef(b[:-1], b[1:])[0, 1]
    neff = a.size * (1 - r1 * r2_) / (1 + r1 * r2_); t = r * np.sqrt(max(neff - 2, 1) / (1 - r ** 2))
    return dict(r=r, n=int(a.size), lag1_x=r1, lag1_y=r2_, n_eff=neff, p_neff=2 * stats.t.sf(abs(t), max(neff - 2, 1)))
T3 = pd.DataFrame({'HCHO vs model isoprene (monthly anomalies)': r_neff(an.ISOP_t, an.dHCHO_SC),
                   'HCHO vs SC 2 m temperature (monthly anomalies)': r_neff(an.T2_C, an.dHCHO_SC),
                   'model isoprene vs SC 2 m temperature (monthly anomalies)': r_neff(an.T2_C, an.ISOP_t)}).T
ann = SER.groupby('year')[['ISOP_t', 'dHCHO_SC', 'T2_C']].mean()
T3.loc['annual means (n = 7, descriptive)'] = dict(r=np.corrcoef(ann.ISOP_t, ann.dHCHO_SC)[0, 1], n=len(ann))
T3.to_csv(os.path.join(OUT_DIR, f'tier3_anomalies_{TAG}.csv')); pd.Series(T2).to_csv(os.path.join(OUT_DIR, f'tier2_seasonal_{TAG}.csv'))
print('Tier 2:', {k: (round(v, 3) if isinstance(v, float) else v) for k, v in T2.items()}); print(T3.round(3).to_string())

# %%
# 6 · Tier 4, column-emission slope (Palmer et al., 2003) -------------------------------------------------------
# HCHO is smeared downwind, so an OLS slope against the unsmoothed emission is diluted. The slope is computed against the
# model field unsmoothed and smoothed at the kernel width that maximised the 0.2° spatial correlation in Tier 1.
frac8 = block(scf, 8)
best = int(T1[T1.block_deg == T1.block_deg.max()].sort_values('pearson_r').smoothing_km.iloc[-1])
T4 = []
for sig_km in sorted({0, best}):
    Es = ndimage.gaussian_filter(np.nan_to_num(E_djf), sig_km / (RES * 111.0)) if sig_km else E_djf
    x = np.where(frac8 >= 0.5, block(np.where(SC, Es, np.nan), 8), np.nan).ravel(); y_ = np.where(frac8 >= 0.5, block(np.where(SC, H_djf, np.nan), 8), np.nan).ravel()
    ok = np.isfinite(x) & np.isfinite(y_); x, y_ = x[ok], y_[ok]
    ols = stats.linregress(x, y_); rma = np.sign(ols.rvalue) * np.std(y_) / np.std(x)
    row = dict(smoothing_km=sig_km, n_blocks=int(ok.sum()), slope_ols_s=ols.slope, slope_ols_se=ols.stderr, intercept_molec_cm2=ols.intercept, r=ols.rvalue, slope_rma_s=rma)
    for Y in (0.29, 0.43):
        row[f'tau_h_Y{Y}_ols'] = ols.slope / Y / 3600; row[f'tau_h_Y{Y}_rma'] = rma / Y / 3600
    T4.append(row)
T4 = pd.DataFrame(T4); T4.to_csv(os.path.join(OUT_DIR, f'tier4_slope_{TAG}{"_biascorr" if BIAS_CORRECT else ""}.csv'), index=False)
print(T4.round(3).to_string(index=False))
print('Reading: S = Y/k_HCHO, so S/Y is an implied HCHO lifetime; "a few hours" (Palmer et al., 2003) is consistent. Much longer, '
      'much shorter or negative values point to the background, the retrieval bias or the emission magnitude; the slope alone cannot '
      'tell which. A non-zero intercept is expected from the residual background and non-isoprene HCHO sources.')

# %%
# 7 · Figures ---------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
ext = [lon[0] - RES / 2, lon[-1] + RES / 2, lat[0] - RES / 2, lat[-1] + RES / 2]
for ax, fld, ttl, cm in ((axs[0, 0], np.where(SC, E_djf, np.nan) / 1e12, 'Model isoprene at overpass, DJF 2019-2025 (10¹² atoms C cm⁻² s⁻¹)', 'Greens'),
                         (axs[0, 1], np.where(SC, H_djf, np.nan) / 1e15, 'TROPOMI ΔHCHO, DJF 2019-2025 (10¹⁵ molec cm⁻²)', 'YlOrBr')):
    im = ax.imshow(fld, origin='lower', extent=ext, cmap=cm); fig.colorbar(im, ax=ax, shrink=0.75); ax.set_title(ttl, fontsize=10)
    ax.set_aspect(1 / np.cos(np.deg2rad(-27.3)))
ax = axs[1, 0]; mm = np.arange(1, 13)
ax.plot(mm, clim.ISOP_t / clim.ISOP_t.max(), 'o-', color='#2e8b57', label='model isoprene (normalised)')
ax.plot(mm, clim.dHCHO_SC / clim.dHCHO_SC.max(), 's-', color='#b8860b', label='TROPOMI ΔHCHO (normalised)')
ax.set_xticks(mm); ax.set_title(f'Tier 2: climatology, r = {r2:.2f}', fontsize=10); ax.legend(frameon=False)
ax = axs[1, 1]; t = pd.to_datetime(dict(year=SER.year, month=SER.month, day=15))
z = lambda s: (s - s.mean()) / s.std()
ax.plot(t, z(an.ISOP_t), color='#2e8b57', label='model isoprene'); ax.plot(t, z(an.dHCHO_SC), color='#b8860b', label='TROPOMI ΔHCHO')
rr = T3.iloc[0]; ax.set_title(f"Tier 3: deseasonalised anomalies, r = {rr.r:.2f} (N_eff = {rr.n_eff:.0f}, p = {rr.p_neff:.3f})", fontsize=10)
ax.legend(frameon=False); ax.axhline(0, color='0.6', lw=0.8)
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, f'hcho_evaluation_{TAG}{"_biascorr" if BIAS_CORRECT else ""}.png'), dpi=250); plt.show()
