# %% [markdown]
# # 16 · Independent check of the 2024-2025 LAI decrease: VIIRS LAI and the MODIS quality flags
#
# In the production run, the statewide mean MODIS LAI of 2025 was 11% below its 2003-2023 mean, and the decrease was
# larger than in any earlier year in cells dominated by trees or herbaceous vegetation (`results/v10/posthoc_review`,
# from `analysis/v10/review/posthoc_review_diagnostics.py`). The orbits of Terra and Aqua have drifted since their orbit
# maintenance ended, so part of the decrease could come from the MODIS record itself. This notebook asks two questions:
#
# 1. Does an independent sensor show a decrease of similar size in 2024-2025? It uses the VIIRS LAI of Suomi NPP
#    (VNP15A2H, Collection 2), the product developed for the continuity of the MODIS LAI record.
# 2. Did the share of MODIS retrievals made by the main algorithm change in 2024-2025?
#
# | Product | Earth Engine collection | 500 m pixels kept for LAI | Main-algorithm flag |
# |---|---|---|---|
# | MODIS Terra + Aqua (as in notebook 10) | `MODIS/061/MCD15A3H` (4-day) | `FparLai_QC` bit 0 = 0 | `FparLai_QC` bit 0 = 0 |
# | VIIRS Suomi NPP | `NASA/VIIRS/002/VNP15A2H` (8-day) | `FparLai_QC` bits 0-2 (SCF_QC) = 0 or 1 | SCF_QC = 0 or 1 (main method, with or without saturation) |
#
# For each month of 2013-2025 the notebook averages the kept 500 m values of the month, takes the area mean of these
# pixels in each 0.025° cell of the inventory grid (as `LAI_AGG = 'mean'` in notebook 10, without the gap filling of
# FIX-4), and samples the monthly fraction of retrievals made by the main algorithm. The Santa Catarina mask and the
# groups of cells dominated by trees, herbaceous vegetation or crops (growth-form fraction above 0.6, from the
# MapBiomas growth forms of the run) are those of `posthoc_review_diagnostics.py`. The reference period is 2013-2023
# because VIIRS starts in January 2012, so the MODIS deviations here differ slightly from those of the inventory
# diagnostics (reference 2003-2023).
#
# **Reading the result.** A 2025 VIIRS deficit of similar size in the same groups indicates that the decrease is not
# specific to MODIS. A VIIRS LAI close to its 2013-2023 mean while MODIS is not points to the MODIS record, and then the
# 2025 emissions inherit a sensor artifact. Only relative deviations are compared, so the result does not depend on the
# scale factor of either product.
#
# **Status.** Run by the authors on 31 October 2026 (outputs in `data/lai_viirs_v10` and `results/v10/lai_viirs_check`).
# That run multiplied the VIIRS band by 0.1, as for MODIS, but the Earth Engine band of VNP15A2H already holds LAI values
# (VIIRS medians printed by cell 5: 0.38-0.40, MODIS 3.90-4.03), so the VIIRS `lai` columns and cache of that run are LAI
# divided by 10. This notebook reads the VIIRS band without that factor. Every deviation, z-score and main-algorithm
# fraction is a ratio or a flag and is unaffected, and so is the correction of notebook 17, which uses VIIRS ratios only. The number
# of Earth Engine requests is 2 per month and product (312 in all).
#
# **Needs:** an Earth Engine project, the repository in Drive (`MEGAN_DIR`) and the output of notebook 10 (`ROOT`, for the
# grid and the growth forms). **Send back:** the folder `OUT_DIR`.

# %%
# 1 · Settings -------------------------------------------------------------------------------------------
__FILE:notebooks/src/_setup_repo.py__
import os, sys, glob, time
try:
    from google.colab import drive
    IN_COLAB = True
except ImportError:
    IN_COLAB = False
if IN_COLAB:
    if os.path.isdir('/content/drive') and not os.path.ismount('/content/drive') and os.listdir('/content/drive'):
        raise RuntimeError('/content/drive holds local files written while Google Drive was not mounted, so Drive cannot be '
                           'mounted there. Use Runtime > Disconnect and delete runtime, then run the notebook again from cell 1.')
    drive.mount('/content/drive')
%pip -q install earthengine-api rasterio geopandas xarray netCDF4 pandas matplotlib
MEGAN_DIR  = REPO_DIR                      # the repository folder (found by the setup block)
ROOT       = '/content/drive/MyDrive/megan_sc_outputs_v10'      # <-- output of notebook 10 (grid and growth forms)
BND        = os.path.join(MEGAN_DIR, 'data', 'boundaries')
OUT_DIR    = os.path.join(ROOT, 'analysis_v10', 'lai_viirs_check')
EE_PROJECT = 'your-ee-project'                                   # <-- your Earth Engine (Google Cloud) project id
YEARS      = list(range(2013, 2026))
REF        = (2013, 2023)                                         # reference years for the deviations
DOMINANT   = 0.6                                                  # growth-form fraction that defines a dominated cell
os.makedirs(OUT_DIR, exist_ok=True)
print('outputs ->', OUT_DIR)

# %%
# 2 · Earth Engine ---------------------------------------------------------------------------------------
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import ee
ee.Authenticate(); ee.Initialize(project=EE_PROJECT)
sys.path.insert(0, os.path.join(MEGAN_DIR, 'src'))
from megan_sc.gee import lai as laimod                          # MODIS: same collection and quality mask as notebook 10
print('Earth Engine ready:', EE_PROJECT, '| MODIS collection', laimod.COLLECTION)

# %%
# 3 · Grid, Santa Catarina mask and growth-form groups (as analysis/v10/review/posthoc_review_diagnostics.py) ----------
import xarray as xr, geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import Affine
with xr.open_dataset(os.path.join(ROOT, '2025', 'emissions', 'annual_emis_2025.nc')) as ds:
    lat, lon = ds.lat.values, ds.lon.values                      # cell centres, south first
res = float(lon[1] - lon[0])
assert abs(res - 0.025) < 1e-9 and lat[0] < lat[-1]
_tr = Affine(res, 0.0, float(lon[0]) - res / 2, 0.0, -res, float(lat[-1]) + res / 2)
_sc = gpd.read_file(os.path.join(BND, 'sc_state_ibge.gpkg')).to_crs('EPSG:4326')
SC = np.flipud(rasterize([(g, 1) for g in _sc.geometry], out_shape=(lat.size, lon.size), transform=_tr, fill=0,
                         dtype='uint8')).astype(bool)                # cell-centre rule, south first
def groups(y):
    f = sorted(glob.glob(os.path.join(ROOT, str(y), 'inputs', 'growthform_*.npz')))[0]
    g = dict(np.load(f))
    return {'all': SC, 'tree': SC & (g['tree'] > DOMINANT), 'herb': SC & (g['herb'] > DOMINANT), 'crop': SC & (g['crop'] > DOMINANT)}
G0 = groups(2019)
print(f'grid {lat.size} x {lon.size} | SC cells {int(SC.sum())} | 2019 cells: ' +
      ', '.join(f'{k} {int(v.sum())}' for k, v in G0.items()))

# %%
# 4 · Monthly images on the inventory grid -------------------------------------------------------------------
region = ee.Geometry.Rectangle([float(lon[0] - res / 2), float(lat[0] - res / 2), float(lon[-1] + res / 2), float(lat[-1] + res / 2)],
                               proj='EPSG:4326', geodesic=False)
VIIRS = 'NASA/VIIRS/002/VNP15A2H'
GRID = {'dimensions': {'width': int(lon.size), 'height': int(lat.size)},
        'affineTransform': {'scaleX': res, 'shearX': 0.0, 'translateX': float(lon[0] - res / 2),
                            'shearY': 0.0, 'scaleY': -res, 'translateY': float(lat[-1] + res / 2)},
        'crsCode': 'EPSG:4326'}

def sample(img, bands):
    """Bands of img on the inventory grid (south first), as sample_multi in notebook 10."""
    for attempt in range(4):
        try:
            a = ee.data.computePixels({'expression': img.select(bands), 'grid': GRID, 'fileFormat': 'NUMPY_NDARRAY'})
            return {b: np.flipud(np.asarray(a[b], dtype='float64')) for b in bands}
        except Exception as e:
            if attempt == 3: raise
            print('   EE retry', attempt + 1, repr(e)[:90]); time.sleep(10 * (attempt + 1))

def to_grid(img, collection):
    """Area mean of the 500 m pixels in each 0.025° cell (explicit aggregation, as LAI_AGG = 'mean' in notebook 10)."""
    proj = ee.ImageCollection(collection).filterBounds(region).first().select('Lai').projection()
    return img.setDefaultProjection(proj).reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)

def modis_month(y, m):
    start = ee.Date.fromYMD(y, m, 1); end = start.advance(1, 'month')
    col = ee.ImageCollection(laimod.COLLECTION).filterBounds(region).filterDate(start, end)
    def flags(img):
        qc = img.select('FparLai_QC'); lai = img.select('Lai')
        made = lai.lte(100).And(qc.rightShift(5).bitwiseAnd(7).neq(4))     # a retrieval exists (SCF_QC != 4) and is a LAI value
        main = made.And(qc.bitwiseAnd(1).eq(0))                           # MODLAND_QC bit 0 = 0: main algorithm
        return made.rename('made').addBands(main.rename('main'))
    s = col.map(flags).sum()
    frac = s.select('main').divide(s.select('made')).updateMask(s.select('made').gt(0)).rename('FRAC')
    lai = laimod.monthly_lai(ee, region, y, m).rename('LAI')               # bit 0 = 0, scale 0.1 (notebook 10)
    return to_grid(lai.addBands(frac), laimod.COLLECTION)

def viirs_month(y, m):
    start = ee.Date.fromYMD(y, m, 1); end = start.advance(1, 'month')
    col = ee.ImageCollection(VIIRS).filterBounds(region).filterDate(start, end)
    def prep(img):
        scf = img.select('FparLai_QC').bitwiseAnd(7); raw = img.select('Lai')
        made = raw.lte(100).And(scf.neq(4))
        main = made.And(scf.lte(1))                                       # SCF_QC 0 or 1: main (RT) method
        lai = raw.updateMask(main).rename('LAI')                          # the Earth Engine band already holds LAI (m2 m-2)
        return lai.addBands(made.rename('made')).addBands(main.rename('main'))
    c = col.map(prep)
    s = c.select(['made', 'main']).sum()
    frac = s.select('main').divide(s.select('made')).updateMask(s.select('made').gt(0)).rename('FRAC')
    return to_grid(c.select('LAI').mean().rename('LAI').addBands(frac), VIIRS)

print('VIIRS images in 2025:', ee.ImageCollection(VIIRS).filterBounds(region).filterDate('2025-01-01', '2026-01-01').size().getInfo())

# %%
# 5 · Sampling (resumable: months already in the cache file are skipped) --------------------------------------
CACHE = os.path.join(OUT_DIR, 'lai_viirs_modis_cache.npz')
cache = dict(np.load(CACHE)) if os.path.exists(CACHE) else {}
for y in YEARS:
    for m in range(1, 13):
        for name, fn in (('MODIS', modis_month), ('VIIRS', viirs_month)):
            k = f'{name}_{y}_{m:02d}'
            if k + '_LAI' in cache: continue
            a = sample(fn(y, m), ['LAI', 'FRAC'])
            cache[k + '_LAI'] = a['LAI'].astype('float32'); cache[k + '_FRAC'] = a['FRAC'].astype('float32')
            if y == YEARS[0] and m <= 2:
                v = a['LAI'][SC]
                print(f'{k}: LAI in SC cells min {np.nanmin(v):.2f} median {np.nanmedian(v):.2f} max {np.nanmax(v):.2f} | '
                      f'main-algorithm fraction median {np.nanmedian(a["FRAC"][SC]):.2f}')
        np.savez_compressed(CACHE, **cache)
    print(y, 'done')

# %%
# 6 · Monthly and annual means by group; deviations from 2013-2023 -------------------------------------------------
rows = []
for y in YEARS:
    G = groups(y)
    for name in ('MODIS', 'VIIRS'):
        L = np.stack([cache[f'{name}_{y}_{m:02d}_LAI'] for m in range(1, 13)]).astype('float64')
        F = np.stack([cache[f'{name}_{y}_{m:02d}_FRAC'] for m in range(1, 13)]).astype('float64')
        L = np.where(L > 0, L, np.nan)
        for m in range(12):
            rows.append(dict(sensor=name, year=y, month=m + 1, lai=float(np.nanmean(L[m][SC])), frac_main=float(np.nanmean(F[m][SC]))))
        ann = np.nanmean(L, 0)
        for gname, msk in G.items():
            rows.append(dict(sensor=name, year=y, month=0, group=gname, lai=float(np.nanmean(ann[msk])), n_cells=int((msk & np.isfinite(ann)).sum())))
D = pd.DataFrame(rows)
MON = D[D.month > 0].drop(columns=['group', 'n_cells'], errors='ignore')
ANN = D[D.month == 0].drop(columns=['month', 'frac_main'], errors='ignore')
MON.to_csv(os.path.join(OUT_DIR, 'lai_viirs_modis_monthly.csv'), index=False)
out = []
for (name, gname), d in ANN.groupby(['sensor', 'group']):
    d = d.set_index('year'); ref = d.loc[REF[0]:REF[1], 'lai']
    for y in d.index:
        out.append(dict(sensor=name, group=gname, year=y, lai=d.loc[y, 'lai'], n_cells=d.loc[y, 'n_cells'],
                        dev_pct=100 * (d.loc[y, 'lai'] / ref.mean() - 1), z=(d.loc[y, 'lai'] - ref.mean()) / ref.std(ddof=1)))
A = pd.DataFrame(out)
A.to_csv(os.path.join(OUT_DIR, 'lai_viirs_modis_annual.csv'), index=False)
summ = []
for (name, gname), d in A.groupby(['sensor', 'group']):
    d = d.set_index('year')
    summ.append(dict(sensor=name, group=gname, dev_2024_pct=round(d.loc[2024, 'dev_pct'], 1), dev_2025_pct=round(d.loc[2025, 'dev_pct'], 1),
                     z_2025=round(d.loc[2025, 'z'], 1), lowest_ref_year_pct=round(d.loc[REF[0]:REF[1], 'dev_pct'].min(), 1)))
S = pd.DataFrame(summ).sort_values(['group', 'sensor'])
S.to_csv(os.path.join(OUT_DIR, 'lai_viirs_modis_summary.csv'), index=False)
clim = MON[MON.year <= REF[1]].groupby(['sensor', 'month']).frac_main.mean()
Q = MON.assign(frac_ref=[clim.loc[(r.sensor, r.month)] for r in MON.itertuples()])
Q = Q.assign(frac_diff=Q.frac_main - Q.frac_ref)
print('Annual LAI deviation from the 2013-2023 mean (%):'); print(S.to_string(index=False))
print('\nMain-algorithm fraction, 2024-2025 minus the 2013-2023 mean of the same month (state mean):')
print(Q[Q.year >= 2024].pivot_table(index=['year', 'month'], columns='sensor', values='frac_diff').round(3).to_string())

# %%
# 7 · Figure ------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3, 1, figsize=(7.2, 8.2), constrained_layout=True)
t = lambda d: pd.to_datetime(dict(year=d.year, month=d.month, day=15))
for name, col in (('MODIS', '#1b9e77'), ('VIIRS', '#d95f02')):
    d = MON[MON.sensor == name]
    mref = d[d.year <= REF[1]].groupby('month').lai.mean()           # 2013-2023 mean of each calendar month
    rel = d.lai.values / d.month.map(mref).values
    axs[0].plot(t(d), 100 * (rel - 1), color=col, lw=0.9, label=name)
    axs[2].plot(t(d), d.frac_main, color=col, lw=0.9, label=name)
axs[0].axhline(0, color='0.6', lw=0.6); axs[0].set_ylabel('% of the 2013-2023\nmean of the month'); axs[0].legend(frameon=False)
axs[0].set_title('(a) Statewide monthly LAI, deviation from the monthly mean', loc='left', fontsize=9, fontweight='bold')
for i, gname in enumerate(('all', 'tree', 'herb', 'crop')):
    for j, (name, col) in enumerate((('MODIS', '#1b9e77'), ('VIIRS', '#d95f02'))):
        d = A[(A.sensor == name) & (A.group == gname)].set_index('year')
        axs[1].bar(i * 3 + j, d.loc[2025, 'dev_pct'], color=col, label=name if i == 0 else None)
        axs[1].plot([i * 3 + j - 0.4, i * 3 + j + 0.4], [d.loc[REF[0]:REF[1], 'dev_pct'].min()] * 2, color='k', lw=1.0)
axs[1].set_xticks([0.5, 3.5, 6.5, 9.5]); axs[1].set_xticklabels(['all cells', 'tree-dominated', 'herbaceous-dominated', 'crop-dominated'])
axs[1].axhline(0, color='0.6', lw=0.6); axs[1].set_ylabel('2025 deviation (%)'); axs[1].legend(frameon=False)
axs[1].set_title('(b) Annual LAI of 2025 versus 2013-2023 (black line: lowest year of 2013-2023)', loc='left', fontsize=9, fontweight='bold')
axs[2].set_ylabel('main-algorithm\nfraction'); axs[2].legend(frameon=False)
axs[2].set_title('(c) Monthly fraction of retrievals by the main algorithm (state mean)', loc='left', fontsize=9, fontweight='bold')
fig.savefig(os.path.join(OUT_DIR, 'lai_viirs_modis.png'), dpi=200); plt.show()
print('written:', sorted(os.listdir(OUT_DIR)))
