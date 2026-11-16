# %% [markdown]
# # 22 · Second and third LAI checks: VIIRS on NOAA-20 (VJ115A2H) and NOAA-21 (VJ215A2H) against VIIRS on Suomi NPP and MODIS
#
# **Why.** The correction of the 2024-2025 MODIS LAI (notebooks 16 and 17, `src/megan_sc/laicorr.py`) takes the
# year-to-year change of each growth-form zone from the VIIRS LAI of Suomi NPP (VNP15A2H, Collection 2). If the Suomi NPP
# record itself drifted, the correction would inherit it. The VIIRS on NOAA-20 is a second instrument, on another
# satellite, processed with the same algorithm (VJ115A2H, version 2; Myneni 2023, NASA LP DAAC,
# https://doi.org/10.5067/VIIRS/VJ115A2H.002; from 1 January 2018). Earth Engine does not hold it, so this notebook reads
# the HDF5 files from NASA (Earthdata login, `earthaccess`).
#
# **What it does.**
# 1. Downloads the 8-day composites of tile h13v11 (Santa Catarina lies inside it) for 2018-2025, for NOAA-20 and, by
#    default, also for Suomi NPP (`VNP15A2H` version 2), so both sensors pass through the same code. Each file is read
#    and deleted.
# 2. Applies the quality rule of notebook 16 (LAI ≤ 10, i.e. raw ≤ 100; `FparLai_QC` bits 0-2 equal to 0 or 1: main
#    radiative-transfer method), averages the kept 500 m values of each month per pixel, then averages the pixels whose
#    centres fall in each 0.025° cell of the inventory grid (notebook 16 uses Earth Engine's area-weighted mean; the
#    pixel-centre mean is checked against it for Suomi NPP).
# 3. Compares the three records by growth-form group with 2018-2023 as reference: annual deviations, NOAA-20/Suomi NPP
#    ratio by year (a drift of either sensor shows here), and deseasonalised monthly anomalies.
# 4. Recomputes the corrected LAI of 2024-2025 with `laicorr.corrected_run_lai` from each VIIRS record with the same
#    reference years (2018-2023): Suomi NPP (Earth Engine cache of notebook 16), Suomi NPP (this pipeline) and NOAA-20.
#    The published correction uses 2013-2023; the 2018-2023 runs isolate the effect of the sensor.
#
# 5. (cell 8) The three VIIRS satellites together. NOAA-21 (VJ215A2H, version 2, https://doi.org/10.5067/VIIRS/VJ215A2H.002)
#    starts on 10 February 2023, so it cannot share the 2018-2023 reference. Cell 8 compares MODIS and the three VIIRS
#    records over February-December of 2023, 2024 and 2025, on the cells valid in all four records of each month: the
#    deviation of 2024 and 2025 from 2023, the ratios between the satellites (paired by calendar month, which shows
#    which record drifts), and the corrected LAI of February-December 2024-2025 from each satellite with 2023 as the
#    only reference year.
#
# **Reading the result.** The statewide corrected LAI from NOAA-20 within ±1.9 % of the one from Suomi NPP (the largest
# error of the leave-one-year-out hindcast of the correction, `results/v10/lai_viirs_check`) confirms the correction.
# A larger gap means the choice of VIIRS instrument matters, and a rerun of notebook 17 with the NOAA-20 field (cell 7)
# is the way to measure its effect on the emissions; cell 7 gives only a linearised estimate.
#
# **Needs** an Earthdata account (free, https://urs.earthdata.nasa.gov) and about 1-2 h in Colab on a first run (about 880
# files of 5-20 MB each are downloaded and deleted one batch at a time); with the cache of an earlier run only the NOAA-21
# files (about 140) are downloaded. No Earth Engine.

# %%
# 1 · Settings ----------------------------------------------------------------------------------------------------------------
__FILE:notebooks/src/_setup_repo.py__
%pip -q install earthaccess h5py numpy pandas scipy rasterio geopandas matplotlib netCDF4
import os, sys, re, glob, json, time, datetime, shutil
try:
    from google.colab import drive; drive.mount('/content/drive')
except Exception:
    pass
# REPO_DIR: this repository, found by the setup block above
ROOT_V10  = '/content/drive/MyDrive/megan_sc_outputs_v10'           # <-- production run (MODIS LAI 2013-2023, growth forms)
ARCHIVE   = '/content/drive/MyDrive/megan_sc_outputs_v10_modis_lai_2024_2025'   # MODIS-LAI 2024-2025 once installed (notebook 19)
NB16_CACHE = ROOT_V10 + '/analysis_v10/lai_viirs_check/lai_viirs_modis_cache.npz'   # notebook 16 (Earth Engine, Suomi NPP and MODIS)
OUT_DIR   = ROOT_V10 + '/analysis_v10/lai_noaa20_check'
WORK      = '/content/lai_h5'
YEARS     = list(range(2018, 2026)); REF = (2018, 2023)
PRODUCTS  = {'VJ1': ('VJ115A2H', '002'), 'VNP': ('VNP15A2H', '002'), 'VJ2': ('VJ215A2H', '002')}   # VJ2 = NOAA-21, from 2023-02-10
THREE_YEARS = (2023, 2024, 2025); THREE_MONTHS = tuple(range(2, 13))      # cell 8: February-December, the months NOAA-21 covers in 2023
TILE      = 'h13v11'
DOMINANT  = 0.6
BATCH     = 23                                                      # files downloaded per batch (deleted after reading)
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(WORK, exist_ok=True)
sys.path.insert(0, REPO_DIR + '/src'); sys.path.insert(0, REPO_DIR + '/tools')
import numpy as np, pandas as pd, h5py
from megan_sc import laicorr
import lai_variant_view

# %%
# 2 · Earthdata login and granule search ----------------------------------------------------------------------------------------
import earthaccess
auth = earthaccess.login(strategy='interactive', persist=True)      # asks for your Earthdata user name and password once
GRAN = {}
for key, (sn, ver) in PRODUCTS.items():
    res_ = earthaccess.search_data(short_name=sn, version=ver, temporal=(f'{YEARS[0]}-01-01', f'{YEARS[-1]}-12-31'),
                                   bounding_box=(-53.9, -29.5, -48.3, -25.9))
    g = [r for r in res_ if f'.{TILE}.' in r.data_links()[0]]
    GRAN[key] = sorted(g, key=lambda r: r.data_links()[0])
    assert GRAN[key], f'no granule of {sn} v{ver} for tile {TILE}'
    print(f'{sn} v{ver}: {len(res_)} granules found, {len(GRAN[key])} of tile {TILE} '
          f'({os.path.basename(GRAN[key][0].data_links()[0])} ... {os.path.basename(GRAN[key][-1].data_links()[0])})')

# %% [markdown]
# ## Cell 3 · Grid, mask, growth-form groups; the 500 m sinusoidal pixels of tile h13v11
#
# MODIS/VIIRS sinusoidal grid: sphere R = 6,371,007.181 m, tiles of 1,111,950.5197665 m, 2400 × 2400 pixels per tile at
# 500 m. Pixel (i, j) of tile (h, v) has its centre at x = (h − 18) T + (j + 0.5) w, y = (9 − v) T − (i + 0.5) w,
# latitude y / R, longitude x / (R cos latitude).

# %%
# 3 · Grid and pixel binning ---------------------------------------------------------------------------------------------------------
import geopandas as gpd, rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from netCDF4 import Dataset
with Dataset(os.path.join(ROOT_V10, '2019', 'emissions', 'annual_emis_2019.nc')) as ds:
    lat = np.asarray(ds.variables['lat'][:], 'float64'); lon = np.asarray(ds.variables['lon'][:], 'float64')
res = float(lon[1] - lon[0])
_tr = Affine(res, 0.0, float(lon[0]) - res / 2, 0.0, -res, float(lat[-1]) + res / 2)
_sc = gpd.read_file(REPO_DIR + '/data/boundaries/sc_state_ibge.gpkg').to_crs('EPSG:4326')
SC = np.flipud(rasterize([(g, 1) for g in _sc.geometry], out_shape=(lat.size, lon.size), transform=_tr, fill=0, dtype='uint8')).astype(bool)
gf = dict(np.load(os.path.join(ROOT_V10, '2019', 'inputs', 'growthform_2019.npz')))
GROUPS = {'all': SC, 'tree': SC & (gf['tree'] > DOMINANT), 'herb': SC & (gf['herb'] > DOMINANT), 'crop': SC & (gf['crop'] > DOMINANT)}
R_SIN, T_TILE, NPIX = 6371007.181, 1111950.5197665233, 2400
W_PIX = T_TILE / NPIX
def tile_lonlat(h, v):
    j = np.arange(NPIX); i = np.arange(NPIX)
    x = (h - 18) * T_TILE + (j + 0.5) * W_PIX; y = (9 - v) * T_TILE - (i + 0.5) * W_PIX
    la = np.degrees(y / R_SIN)
    lo = np.degrees(x[None, :] / (R_SIN * np.cos(np.radians(la))[:, None]))
    return lo, np.broadcast_to(la[:, None], lo.shape)
h_, v_ = int(TILE[1:3]), int(TILE[4:6])
PLON, PLAT = tile_lonlat(h_, v_)
IX = np.floor((PLON - (lon[0] - res / 2)) / res).astype(np.int64); IY = np.floor((PLAT - (lat[0] - res / 2)) / res).astype(np.int64)
INGRID = (IX >= 0) & (IX < lon.size) & (IY >= 0) & (IY < lat.size)
CELL = np.where(INGRID, IY * lon.size + IX, -1)
print(f'tile {TILE}: {int(INGRID.sum()):,} of {NPIX * NPIX:,} pixels fall in the inventory grid; SC cells {int(SC.sum())}; '
      f'pixels per SC cell, median {np.median(np.bincount(CELL[INGRID], minlength=lat.size * lon.size).reshape(lat.size, lon.size)[SC]):.0f}')
assert np.bincount(CELL[INGRID], minlength=lat.size * lon.size).reshape(lat.size, lon.size)[SC].min() > 0, 'some SC cells get no pixel'

# %% [markdown]
# ## Cell 4 · Reading the HDF5 files
#
# The dataset names are found by walking the file (the HDF-EOS5 group names differ between products and versions). The
# scale factor, fill value and valid range are read from the attributes of the `Lai` dataset and printed for the first
# file of each product; raw values above 100 (LAI 10) are fill or special codes in both products and are dropped, as in
# notebook 16.

# %%
# 4 · HDF5 reader ----------------------------------------------------------------------------------------------------------------------------
def find(h5, name):
    hits = []
    h5.visititems(lambda k, o: hits.append(k) if isinstance(o, h5py.Dataset) and k.split('/')[-1] == name else None)
    assert len(hits) == 1, (name, hits)
    return h5[hits[0]]
SHOWN = set()
def read_lai(path, key):
    with h5py.File(path, 'r') as h5:
        d = find(h5, 'Lai'); qc = find(h5, 'FparLai_QC')
        raw = d[()].astype('float64'); q = qc[()].astype(np.int64)
        att = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in d.attrs.items()}
        if key not in SHOWN:
            print(f'{key}: {d.name} {d.shape} {d.dtype} attributes {att}'); SHOWN.add(key)
    sf = float(np.ravel(att.get('scale_factor', [0.1]))[0])
    ok = (raw <= 100) & ((q & 7) <= 1)                                 # notebook 16: raw <= 100 and SCF_QC (bits 0-2) in {0, 1}
    return np.where(ok, raw * sf, np.nan)
def start_date(path):
    m = re.search(r'\.A(\d{4})(\d{3})\.', os.path.basename(path))
    return datetime.date(int(m[1]), 1, 1) + datetime.timedelta(days=int(m[2]) - 1)

# %% [markdown]
# ## Cell 5 · Monthly fields on the inventory grid (resumable)
#
# For each product and month: per-pixel mean of the kept values of the composites that start in that month (as
# Earth Engine's `filterDate` on the start time, notebook 16), then the mean over the pixels of each 0.025° cell.
# Months already in the cache are skipped.

# %%
# 5 · Download, read, bin --------------------------------------------------------------------------------------------------------------------
CACHE = os.path.join(OUT_DIR, 'lai_viirs_h5_cache.npz')
cache = dict(np.load(CACHE)) if os.path.exists(CACHE) else {}
NC = lat.size * lon.size
def bin_month(stack_sum, stack_n):
    pix = np.where(stack_n > 0, stack_sum / np.maximum(stack_n, 1), np.nan)
    good = INGRID & np.isfinite(pix)
    s = np.bincount(CELL[good], weights=pix[good], minlength=NC); n = np.bincount(CELL[good], minlength=NC)
    return np.where(n > 0, s / np.maximum(n, 1), np.nan).reshape(lat.size, lon.size), n.reshape(lat.size, lon.size)
for key, grans in GRAN.items():
    by_month = {}
    for g in grans:
        d0 = start_date(g.data_links()[0])
        if d0.year in YEARS:
            by_month.setdefault((d0.year, d0.month), []).append(g)
    todo = [ym for ym in sorted(by_month) if f'{key}_{ym[0]}_{ym[1]:02d}_LAI' not in cache]
    print(f'{key}: {len(by_month)} months, {len(todo)} to do')
    for ym in todo:
        files = earthaccess.download(by_month[ym], local_path=WORK, threads=8)
        S = np.zeros((NPIX, NPIX)); N = np.zeros((NPIX, NPIX))
        for f in files:
            a = read_lai(f, key); ok = np.isfinite(a); S[ok] += a[ok]; N[ok] += 1; os.remove(f)
        field, npx = bin_month(S, N)
        cache[f'{key}_{ym[0]}_{ym[1]:02d}_LAI'] = field.astype('float32'); cache[f'{key}_{ym[0]}_{ym[1]:02d}_NPIX'] = npx.astype('int16')
        np.savez_compressed(CACHE, **cache)
        print(f'  {key} {ym[0]}-{ym[1]:02d}: {len(files)} composites, SC mean LAI {np.nanmean(field[SC]):.3f}, '
              f'SC cells with a value {int(np.isfinite(field[SC]).sum())}', flush=True)
for key in GRAN:                                                    # caches in the format laicorr reads (keys VIIRS_<y>_<mm>_LAI)
    np.savez_compressed(os.path.join(OUT_DIR, f'lai_{key.lower()}_as_viirs_cache.npz'),
                        **{f'VIIRS_{y}_{m:02d}_LAI': cache[f'{key}_{y}_{m:02d}_LAI'] for y in YEARS for m in range(1, 13)
                           if f'{key}_{y}_{m:02d}_LAI' in cache})

# %% [markdown]
# ## Cell 6 · The three records, 2018-2023 reference
#
# * pipeline check: Suomi NPP from this notebook against the Earth Engine cache of notebook 16 (cell values, every
#   month; correlation and mean ratio; the Earth Engine VIIRS values of that cache are ten times too small, a constant
#   that cancels in every ratio);
# * annual deviation from the 2018-2023 mean by group, for MODIS (the LAI the production run used; 2024-2025 from the
#   MODIS-LAI years), Suomi NPP and NOAA-20;
# * the NOAA-20/Suomi NPP ratio of the annual state means, by year;
# * deseasonalised monthly state anomalies (deviation from the 2018-2023 mean of the calendar month), correlations.

# %%
# 6 · Comparison ----------------------------------------------------------------------------------------------------------------------------------
C16 = np.load(NB16_CACHE)
MOD_ROOT = lai_variant_view.modis_root(ROOT_V10, ARCHIVE, '/content/run_v10_modis_view')
def read_sf(p):
    with rasterio.open(p) as s:
        return np.flipud(s.read(1).astype('float64'))
REC = {'MODIS': lambda y, m: read_sf(os.path.join(MOD_ROOT, str(y), 'inputs', f'lai_{y}_{m:02d}.tif')),
       'VNP_EE': lambda y, m: C16[f'VIIRS_{y}_{m:02d}_LAI'].astype('float64')}
for key in GRAN:
    REC[key] = (lambda k: (lambda y, m: cache[f'{k}_{y}_{m:02d}_LAI'].astype('float64')))(key)
FULL = [k for k in GRAN if all(f'{k}_{y}_{m:02d}_LAI' in cache for y in YEARS for m in range(1, 13))]   # cells 6-7: whole 2018-2025
print('records with every month of 2018-2025 (cells 6 and 7):', FULL, '| shorter records (cell 8 only):', [k for k in GRAN if k not in FULL])
REC6 = {k: REC[k] for k in ['MODIS', 'VNP_EE'] + FULL}
pos = lambda a: np.where(np.isfinite(a) & (a > 0), a, np.nan)
if 'VNP' in GRAN:
    x, yv = [], []
    for y in YEARS:
        for m in range(1, 13):
            a, b = pos(REC['VNP'](y, m))[SC], pos(REC['VNP_EE'](y, m))[SC]; ok = np.isfinite(a) & np.isfinite(b)
            x.append(a[ok]); yv.append(b[ok])
    x, yv = np.concatenate(x), np.concatenate(yv)
    print(f'Suomi NPP, this pipeline vs Earth Engine (notebook 16), all SC cell-months: r = {np.corrcoef(x, yv)[0, 1]:.3f}, '
          f'median ratio {np.median(x / yv):.3f} (about 10 expected: the Earth Engine cache is ten times too small)')
rows = []
for rec, fn in REC6.items():
    for y in YEARS:
        L = np.stack([pos(fn(y, m)) for m in range(1, 13)])
        for m in range(12):
            rows.append(dict(record=rec, year=y, month=m + 1, group='all', lai=float(np.nanmean(L[m][SC]))))
        ann = np.nanmean(L, 0)
        for gname, msk in GROUPS.items():
            rows.append(dict(record=rec, year=y, month=0, group=gname, lai=float(np.nanmean(ann[msk]))))
D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT_DIR, 'lai_three_records_monthly_annual.csv'), index=False)
A = D[D.month == 0].copy()
A['ref'] = A.groupby(['record', 'group']).lai.transform(lambda s: s[A.loc[s.index, 'year'].between(*REF)].mean())
A['dev_pct'] = 100 * (A.lai / A.ref - 1)
DEV = A.pivot_table(index=['group', 'year'], columns='record', values='dev_pct').round(2)
DEV.to_csv(os.path.join(OUT_DIR, 'lai_annual_deviation_2018_2023ref.csv'))
print('annual deviation from the 2018-2023 mean (%):'); print(DEV.to_string())
if 'VNP' in GRAN:
    al = A[A.group == 'all'].pivot(index='year', columns='record', values='lai')
    rat = (al['VJ1'] / al['VNP']).rename('VJ1_over_VNP'); rat.to_csv(os.path.join(OUT_DIR, 'lai_noaa20_over_snpp_by_year.csv'))
    from scipy import stats
    tr = stats.theilslopes(rat.values, rat.index.values)
    print('NOAA-20 / Suomi NPP, state annual mean, by year:', rat.round(4).to_dict())
    print(f'trend of the ratio: {100 * tr[0] / rat.mean():+.2f} % per year (Theil-Sen 95 % interval {100 * tr[2] / rat.mean():+.2f} to {100 * tr[3] / rat.mean():+.2f})')
M = D[(D.month > 0) & (D.group == 'all')].copy()
M['clim'] = M.groupby(['record', 'month']).lai.transform(lambda s: s[M.loc[s.index, 'year'].between(*REF)].mean())
M['anom_pct'] = 100 * (M.lai / M.clim - 1)
AN = M.pivot_table(index=['year', 'month'], columns='record', values='anom_pct')
AN.to_csv(os.path.join(OUT_DIR, 'lai_monthly_anomalies_2018_2023ref.csv'))
print('correlation of the deseasonalised monthly state anomalies, 2018-2023:'); print(AN.loc[REF[0]:REF[1]].corr().round(2).to_string())
print('mean monthly anomaly 2024 and 2025 (%):'); print(AN.groupby(level=0).mean().loc[[2024, 2025]].round(2).to_string())

# %% [markdown]
# ## Cell 7 · The corrected LAI of 2024-2025 from each VIIRS record (same function, reference 2018-2023)
#
# `laicorr.corrected_run_lai(MODIS record, VIIRS cache, SC mask, years, ref_years)`: published settings (one ratio per
# growth-form zone and month). Statewide annual means are compared with the published correction (reference
# 2013-2023). The emission change is a linear estimate from the class changes of notebook 17 (cell 11) per percent of
# LAI change; it is only indicative.

# %%
# 7 · Corrected LAI from each record ----------------------------------------------------------------------------------------------------------------
def state_annual(fields, y):          # mean over the state cells with LAI > 0 (water and urban cells excluded), as in cell 6
    return float(np.mean([np.nanmean(pos(fields[(y, m)])[SC]) for m in range(1, 13)]))
REF_Y = list(range(REF[0], REF[1] + 1))
CORR = {'published (VNP, Earth Engine, 2013-2023)': laicorr.corrected_run_lai(MOD_ROOT, NB16_CACHE, SC, years=(2024, 2025))[0],
        'VNP, Earth Engine, 2018-2023': laicorr.corrected_run_lai(MOD_ROOT, NB16_CACHE, SC, years=(2024, 2025), ref_years=REF_Y)[0]}
for key in FULL:
    CORR[f'{key}, this notebook, 2018-2023'] = laicorr.corrected_run_lai(MOD_ROOT, os.path.join(OUT_DIR, f'lai_{key.lower()}_as_viirs_cache.npz'),
                                                                         SC, years=(2024, 2025), ref_years=REF_Y)[0]
modis_used = {(y, m): REC['MODIS'](y, m) for y in (2024, 2025) for m in range(1, 13)}
rows = []
for y in (2024, 2025):
    base = state_annual(modis_used, y)
    rows.append(dict(year=y, correction='MODIS (no correction)', state_lai=round(base, 3), vs_modis_pct=0.0))
    for lab, F in CORR.items():
        v = state_annual(F, y); rows.append(dict(year=y, correction=lab, state_lai=round(v, 3), vs_modis_pct=round(100 * (v / base - 1), 2)))
CT = pd.DataFrame(rows)
ref_lab = 'VNP, Earth Engine, 2018-2023'
CT['vs_vnp_2018_2023_pct'] = [round(100 * (r.state_lai / CT[(CT.year == r.year) & (CT.correction == ref_lab)].state_lai.iloc[0] - 1), 2) for r in CT.itertuples()]
CT.to_csv(os.path.join(OUT_DIR, 'corrected_lai_by_viirs_record.csv'), index=False)
print(CT.to_string(index=False))
NB17 = {2024: dict(lai=5.0, BVOC=3.2, ISOP=2.3, MT=3.4, SQT=4.4), 2025: dict(lai=11.6, BVOC=7.0, ISOP=4.3, MT=7.9, SQT=9.9)}   # notebook 17, cell 11
for y in (2024, 2025):
    r = CT[(CT.year == y) & CT.correction.str.startswith('VJ1')]
    if len(r):
        d = float(r.vs_vnp_2018_2023_pct.iloc[0])
        assert np.isfinite(d), f'{y}: the state LAI of the NOAA-20 correction is not finite'
        est = {k: round(NB17[y][k] / NB17[y]['lai'] * d, 2) for k in ('BVOC', 'ISOP', 'MT', 'SQT')}
        verdict = 'within' if abs(d) <= 1.9 else 'OUTSIDE'
        print(f'{y}: NOAA-20 correction differs from the Suomi NPP one (same reference) by {d:+.2f} % of state LAI, {verdict} the '
              f'hindcast bound of 1.9 %; linearised emission change (%): {est}')

# %% [markdown]
# ## Cell 8 · The three VIIRS satellites, February-December 2023-2025, 2023 as reference
#
# * Common cells: for each month, the state cells with LAI > 0 in MODIS and in every VIIRS record, so that differences in
#   cloud gaps do not enter the comparison.
# * Deviation of the February-December mean of 2024 and 2025 from that of 2023, by record and growth-form group.
# * Ratios between the satellites: for each calendar month, ln(A/B) in 2024 or 2025 minus ln(A/B) in 2023 (paired by
#   month, 11 pairs; mean with a t interval, Wilcoxon signed-rank p). A ratio that changes says that one of the two
#   records drifts relative to the other; with three satellites the drifting one is the one whose two ratios change.
# * Corrected LAI of February-December 2024-2025 from each satellite with 2023 as the only reference year
#   (`laicorr.correct_month`, published settings), against the published correction (Suomi NPP, 2013-2023) over the same
#   months. A single reference year makes these corrections noisier than the published one, so they compare the
#   satellites, not the reference periods; the line "VNP, Earth Engine, 2023" shows the effect of the reference alone.

# %%
# 8 · Three VIIRS satellites -------------------------------------------------------------------------------------------------------------------
from scipy import stats
R8 = ['MODIS', 'VNP', 'VJ1', 'VJ2']
miss = [(k, y, m) for k in R8[1:] for y in THREE_YEARS for m in THREE_MONTHS if f'{k}_{y}_{m:02d}_LAI' not in cache]
assert not miss, f'months missing in the cache for cell 8: {miss[:6]}'
F8 = {(k, y, m): pos(REC[k](y, m)) for k in R8 for y in THREE_YEARS for m in THREE_MONTHS}
COMMON = {(y, m): SC & np.logical_and.reduce([np.isfinite(F8[(k, y, m)]) for k in R8]) for y in THREE_YEARS for m in THREE_MONTHS}
rows = []
for (y, m), cm in COMMON.items():
    for k in R8:
        for gname, msk in GROUPS.items():
            c = cm & msk
            rows.append(dict(record=k, year=y, month=m, group=gname, cells=int(c.sum()), lai=float(np.mean(F8[(k, y, m)][c]))))
T8 = pd.DataFrame(rows); T8.to_csv(os.path.join(OUT_DIR, 'three_satellites_monthly_2023_2025.csv'), index=False)
print('common state cells per month: min', int(T8[(T8.group == 'all') & (T8.record == 'MODIS')].cells.min()),
      'median', int(T8[(T8.group == 'all') & (T8.record == 'MODIS')].cells.median()), 'of', int(SC.sum()))
AN8 = T8.groupby(['group', 'record', 'year']).lai.mean().unstack('year')
DV8 = pd.DataFrame({f'dev_{y}_vs_2023_pct': 100 * (AN8[y] / AN8[2023] - 1) for y in (2024, 2025)}).round(2)
DV8.to_csv(os.path.join(OUT_DIR, 'three_satellites_deviation_feb_dec_vs_2023.csv'))
print('\nFebruary-December mean, deviation from 2023 (%), common cells:'); print(DV8.unstack('record').to_string())
PAIRS = [('VJ1', 'VNP'), ('VJ2', 'VNP'), ('VJ2', 'VJ1')]
al8 = T8[T8.group == 'all'].pivot_table(index=['year', 'month'], columns='record', values='lai')
rows = []
for a, b in PAIRS:
    lr = np.log(al8[a] / al8[b])
    for y in (2024, 2025):
        d = (lr.loc[y] - lr.loc[2023]).values; h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size)
        rows.append(dict(ratio=f'{a}/{b}', year=y, ratio_2023=float(np.exp(lr.loc[2023].mean())), ratio_year=float(np.exp(lr.loc[y].mean())),
                         change_pct=100 * (np.exp(d.mean()) - 1), lo_pct=100 * (np.exp(d.mean() - h) - 1), hi_pct=100 * (np.exp(d.mean() + h) - 1),
                         wilcoxon_p=float(stats.wilcoxon(d).pvalue), months=int(d.size)))
RT8 = pd.DataFrame(rows); RT8.to_csv(os.path.join(OUT_DIR, 'three_satellites_ratio_change_vs_2023.csv'), index=False)
print('\nratios between the satellites, change from 2023 (paired by calendar month):'); print(RT8.round(4).to_string(index=False))
# corrected LAI, February-December, 2023 as the only reference year
zones8 = {y: laicorr.growth_form_zones(dict(np.load(sorted(glob.glob(os.path.join(MOD_ROOT, str(y), 'inputs', 'growthform_*.npz')))[0])), SC, laicorr.DOMINANT)
          for y in (2024, 2025)}
V8 = dict(REC); V8['VNP_EE'] = REC['VNP_EE']
pub = CORR['published (VNP, Earth Engine, 2013-2023)']
rows = []
for y in (2024, 2025):
    feb_dec = lambda F: float(np.mean([np.nanmean(pos(F[(y, m)])[SC]) for m in THREE_MONTHS]))
    base = feb_dec(modis_used); p_ = feb_dec(pub)
    rows.append(dict(year=y, correction='MODIS (no correction)', state_lai_feb_dec=base))
    rows.append(dict(year=y, correction='published (VNP, Earth Engine, 2013-2023)', state_lai_feb_dec=p_))
    for k in ['VNP_EE', 'VNP', 'VJ1', 'VJ2']:
        F = {(y, m): laicorr.correct_month([REC['MODIS'](2023, m)], [V8[k](2023, m)], REC['MODIS'](y, m), V8[k](y, m),
                                           ref_mask=SC, zones=zones8[y])[0] for m in THREE_MONTHS}
        lab = {'VNP_EE': 'VNP, Earth Engine, 2023', 'VNP': 'VNP, this notebook, 2023', 'VJ1': 'VJ1 (NOAA-20), 2023', 'VJ2': 'VJ2 (NOAA-21), 2023'}[k]
        rows.append(dict(year=y, correction=lab, state_lai_feb_dec=feb_dec(F)))
C8 = pd.DataFrame(rows)
for y in (2024, 2025):
    s_ = C8.year == y
    C8.loc[s_, 'vs_modis_pct'] = 100 * (C8.loc[s_, 'state_lai_feb_dec'] / C8.loc[s_ & (C8.correction == 'MODIS (no correction)'), 'state_lai_feb_dec'].iloc[0] - 1)
    C8.loc[s_, 'vs_published_pct'] = 100 * (C8.loc[s_, 'state_lai_feb_dec'] / C8.loc[s_ & C8.correction.str.startswith('published'), 'state_lai_feb_dec'].iloc[0] - 1)
    C8.loc[s_, 'vs_vnp_2023_pct'] = 100 * (C8.loc[s_, 'state_lai_feb_dec'] / C8.loc[s_ & (C8.correction == 'VNP, this notebook, 2023'), 'state_lai_feb_dec'].iloc[0] - 1)
C8 = C8.round(3); C8.to_csv(os.path.join(OUT_DIR, 'three_satellites_corrected_lai_feb_dec.csv'), index=False)
print('\ncorrected LAI, February-December state mean, 2023 as reference:'); print(C8.to_string(index=False))
sp = C8[C8.correction.isin(['VNP, this notebook, 2023', 'VJ1 (NOAA-20), 2023', 'VJ2 (NOAA-21), 2023'])].groupby('year').state_lai_feb_dec
print('spread of the three satellites (max/min - 1, %):', (100 * (sp.max() / sp.min() - 1)).round(2).to_dict())
json.dump(dict(notebook='22_lai_check_noaa20_viirs', products={k: list(v) for k, v in PRODUCTS.items()}, tile=TILE, years=YEARS, ref=REF,
               three_satellites=dict(years=list(THREE_YEARS), months=list(THREE_MONTHS), reference=2023),
               granules={k: len(v) for k, v in GRAN.items()}, run_utc=datetime.datetime.utcnow().isoformat(timespec='seconds') + 'Z'),
          open(os.path.join(OUT_DIR, 'provenance_notebook22.json'), 'w'), indent=1)
print('outputs in', OUT_DIR)
