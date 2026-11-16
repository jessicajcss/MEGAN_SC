# %% [markdown]
# # 15 · Sensitivity tests (2019) and diagnostics from the fix logs of the production run
#
# Quantifies the choices that a reviewer is most likely to question, using runs of notebook 10 for 2019:
#
# | Test | Runs needed (notebook 10, `YEARS = [2019]`) | What is reported |
# |---|---|---|
# | S1 temporal sampling | `VARIANT = 'base'` (3-hourly, production) and `VARIANT = 'hourly'` | class totals, monthly totals and the isoprene map, 3-hourly / hourly |
# | S2 LAI aggregation | `VARIANT = 'lai_eedefault'` (Earth Engine default resampling) | class totals and SC-mean LAI, explicit 500 m area mean vs default |
# | S3 soil-moisture factor (GAMSM on in production; off in the official script) | fix logs of the base run, all years | SC-mean GAMSM by month and the approximate reduction of the non-soil classes |
#
# The 2019 runs cost about one production year each (the hourly run about three). S3 needs no extra run. Isoprene is read
# times the canopy loss and production factor rho = 0.96 of Guenther et al. (2006, Sect. 3.3) (RHO_ISOP), as in notebook 13;
# the ratios of S1 and S2 do not depend on it.
#
# **Send back:** the folder `OUT_DIR`.

# %%
# 1 · Paths -----------------------------------------------------------------------------------------------
__FILE:notebooks/src/_setup_repo.py__
import os, glob, warnings
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
%pip -q install rasterio geopandas netCDF4 pandas matplotlib
BASE_ROOT   = '/content/drive/MyDrive/megan_sc_outputs_v10'
HOURLY_ROOT = '/content/drive/MyDrive/megan_sc_outputs_v10_hourly'
LAIDEF_ROOT = '/content/drive/MyDrive/megan_sc_outputs_v10_lai_eedefault'
BND         = os.path.join(REPO_DIR, 'data/boundaries')
OUT_DIR     = os.path.join(BASE_ROOT, 'analysis_v10', 'sensitivity')
YEAR = 2019; YEARS_ALL = list(range(2003, 2026))
RHO_ISOP = 0.96   # isoprene canopy loss and production factor (Guenther et al., 2006); notebook 10 writes rho = 1

# %%
# 2 · Helpers -----------------------------------------------------------------------------------------------
os.makedirs(OUT_DIR, exist_ok=True)
import numpy as np, pandas as pd, rasterio, geopandas as gpd, matplotlib.pyplot as plt
from netCDF4 import Dataset
from rasterio.features import rasterize
from rasterio.transform import from_origin
warnings.filterwarnings('ignore', category=RuntimeWarning)
CLASSES = ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM', 'NO', 'SQT_HR', 'SQT_LR',
           'MEOH', 'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO']
BVOC = [c for c in CLASSES if c not in ('NO', 'CO')]
def month_grid(root, y, m):
    f = os.path.join(root, str(y), 'emissions', f'emis_{y}_{m:02d}.nc')
    if not os.path.exists(f): return None
    with Dataset(f) as ds:
        return ({c: np.nan_to_num(np.array(ds.variables[c][:], dtype='float64')) * (RHO_ISOP if c == 'ISOP' else 1.0) for c in CLASSES},
                np.array(ds.variables['lat'][:]), np.array(ds.variables['lon'][:]))
first = month_grid(BASE_ROOT, YEAR, 1); assert first is not None, 'base run for the year not found'
_, lat, lon = first; NY, NX = lat.size, lon.size; RES = float(lon[1] - lon[0])
g = gpd.read_file(os.path.join(BND, 'sc_state_ibge.gpkg')).to_crs('EPSG:4326')
geom = g.geometry.union_all() if hasattr(g.geometry, 'union_all') else g.geometry.unary_union
SC = np.flipud(rasterize([(geom, 1)], out_shape=(NY, NX), transform=from_origin(lon[0] - RES / 2, lat[-1] + RES / 2, RES, RES), fill=0, dtype='uint8')).astype(bool)
def monthly_totals(root, y):
    rows = {}
    for m in range(1, 13):
        r = month_grid(root, y, m)
        if r is None: return None
        rows[m] = {c: float((r[0][c] * SC).sum()) for c in CLASSES}
    T = pd.DataFrame(rows).T; T['BVOC'] = T[BVOC].sum(1); T.index.name = 'month'; return T
def annual_grid(root, y, c):
    return sum(month_grid(root, y, m)[0][c] for m in range(1, 13))
def compare(root_alt, label):
    base = monthly_totals(BASE_ROOT, YEAR); alt = monthly_totals(root_alt, YEAR) if os.path.isdir(root_alt) else None
    if alt is None:
        print(f'{label}: run not found (or incomplete) under {root_alt}; skipped'); return None
    ann = pd.DataFrame({'base_t': base.sum(), f'{label}_t': alt.sum()}); ann[f'base/{label}'] = ann.base_t / ann[f'{label}_t']
    mon = (base['BVOC'] / alt['BVOC']).rename(f'BVOC base/{label}').to_frame(); mon[f'ISOP base/{label}'] = base['ISOP'] / alt['ISOP']
    ann.to_csv(os.path.join(OUT_DIR, f'{label}_{YEAR}_annual.csv')); mon.to_csv(os.path.join(OUT_DIR, f'{label}_{YEAR}_monthly.csv'))
    print(f'--- {label}: annual ratio base / {label}'); print(ann.round(4).to_string()); print(mon.round(4).to_string())
    return ann, mon

# %%
# 3 · S1 temporal sampling: 3-hourly (production) vs hourly ----------------------------------------------------
S1 = compare(HOURLY_ROOT, 'hourly')
if S1 is not None:
    ib = annual_grid(BASE_ROOT, YEAR, 'ISOP'); ih = annual_grid(HOURLY_ROOT, YEAR, 'ISOP')
    ratio = np.where(SC & (ih > 0), ib / np.where(ih > 0, ih, 1), np.nan)
    v = ratio[np.isfinite(ratio)]
    print(f'isoprene cell ratio 3-hourly/hourly: median {np.median(v):.4f}, 1st-99th percentile {np.percentile(v, 1):.4f}-{np.percentile(v, 99):.4f}')
    fig, ax = plt.subplots(figsize=(8, 5.5)); im = ax.imshow(ratio, origin='lower', cmap='RdBu_r', vmin=0.9, vmax=1.1,
                                                          extent=[lon[0], lon[-1], lat[0], lat[-1]])
    ax.set_aspect(1 / np.cos(np.deg2rad(-27.3))); fig.colorbar(im, ax=ax, shrink=0.8, label='3-hourly / hourly')
    ax.set_title(f'Isoprene {YEAR}: 3-hourly sampling relative to hourly'); fig.savefig(os.path.join(OUT_DIR, f'S1_isoprene_ratio_{YEAR}.png'), dpi=200); plt.show()

# %%
# 4 · S2 LAI aggregation: explicit 500 m area mean (production) vs Earth Engine default resampling ------------------
S2 = compare(LAIDEF_ROOT, 'lai_eedefault')
if S2 is not None:
    def lai_mean(root):
        vals = []
        for m in range(1, 13):
            with rasterio.open(os.path.join(root, str(YEAR), 'inputs', f'lai_{YEAR}_{m:02d}.tif')) as s: L = np.flipud(s.read(1).astype('float64'))
            vals.append(np.nanmean(np.where(SC & (L > 0), L, np.nan)))
        return np.array(vals)
    lb, ld = lai_mean(BASE_ROOT), lai_mean(LAIDEF_ROOT)
    print('SC-mean LAI by month, base vs EE default:', np.round(lb, 3), np.round(ld, 3), '| annual ratio', round(lb.mean() / ld.mean(), 4))

# %%
# 5 · S3 soil-moisture factor, from the fix logs (all years of the base run) ---------------------------------------
rows = []
for y in YEARS_ALL:
    for m in range(1, 13):
        f = os.path.join(BASE_ROOT, str(y), 'inputs', f'fixlog_{y}_{m:02d}.npz')
        r = month_grid(BASE_ROOT, y, m)
        if not os.path.exists(f) or r is None: continue
        z = np.load(f); E = r[0]; bv = sum(E[c] for c in BVOC); gsm = z['gamsm_mean']
        ok = SC & (bv > 0) & np.isfinite(gsm) & (gsm > 0)
        rows.append(dict(year=y, month=m, gamsm_sc_mean=float(np.nanmean(np.where(SC, gsm, np.nan))),
                         frac_cells_gamsm_lt_0_99=float(np.mean(gsm[SC] < 0.99)),
                         bvoc_reduction_approx=float(1 - bv[ok].sum() / (bv[ok] / gsm[ok]).sum()),
                         BVOC_t=float((bv * SC).sum())))
D = pd.DataFrame(rows)
if len(D):
    D.to_csv(os.path.join(OUT_DIR, 'S3_fixlog_monthly.csv'), index=False)
    A = D.groupby('year')[['BVOC_t']].sum()
    A['bvoc_reduction_by_gamsm_approx'] = D.groupby('year').apply(lambda d: np.average(d.bvoc_reduction_approx, weights=d.BVOC_t))
    A.to_csv(os.path.join(OUT_DIR, 'S3_fixlog_annual.csv'))
    print('S3 GAMSM: SC-mean factor by month (climatology):'); print(D.groupby('month').gamsm_sc_mean.mean().round(4).to_string())
    print('approximate BVOC reduction by GAMSM (emission-weighted over all months): %.2f %%' % (100 * np.average(D.bvoc_reduction_approx, weights=D.BVOC_t)))
    print('(approximate because the fix log keeps the monthly mean GAMSM per cell, not its hourly covariance with emission)')
else:
    print('no fix logs found: S3 needs the notebook-10 base run')
