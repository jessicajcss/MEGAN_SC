# %% [markdown]
# # 25 · Results of a run, for any region and period
#
# Reads the output folder of `10_production_run.ipynb` (after `11_verify_run.ipynb` passes) for the configuration
# `RUN_CONFIG`, and writes the standard results of an inventory into `OUT_DIR`:
#
# | File | Content |
# |---|---|
# | `monthly_totals_by_class.csv` | totals of the region for every month of the run, t month⁻¹, 19 classes |
# | `annual_totals_by_class.csv`, `annual_totals_by_group.csv` | totals of the region for every complete year, t yr⁻¹ |
# | `seasonal_climatology.csv`, `fig_seasonal_cycle.png` | mean, minimum and maximum of each calendar month over the complete years |
# | `fig_annual_series.png` | annual totals of the groups |
# | `mean_intensity_<group>.tif`, `fig_maps_mean_intensity.png` | mean annual emission per area of each cell, kg ha⁻¹ yr⁻¹ |
# | `drivers_annual.csv`, `correlations_emissions_drivers.csv` | annual regional means of temperature, PAR, soil water and LAI as the model saw them; Pearson r of the raw and linearly detrended series (10 complete years or more) |
# | `trends.csv` | Theil-Sen slope with its 95 % interval, Mann-Kendall p and Hamed-Rao p (10 complete years or more) |
# | `zones_annual_totals.csv` | annual totals of each zone of the configuration (when notebook 10 wrote per-zone tables) |
#
# Definitions (tutorial, chapters 1 and 7):
# * a cell belongs to the region when its centre lies inside the boundary of the configuration;
# * isoprene is multiplied by the canopy loss and production factor ρ = 0.96 of Guenther et al. (2006, Sect. 3.3)
#   when it is read (`RHO_ISOP`); notebook 10 writes isoprene as the MEGANv3.21 code computes it, without ρ;
# * groups: BVOC = the 17 classes other than soil NO and CO; MT = the five monoterpene classes; SQT = SQT_HR + SQT_LR;
#   OVOC = MEOH + ACTO + ETOH + ACID + OXPROD; STRESS_OTHER = STRESS + OTHER; SOA = MT + SQT (precursors of secondary
#   organic aerosol, not aerosol formed);
# * intensity per area divides by the cell area on the authalic sphere (`megan_sc.domain`, the areas of the engine);
# * trends: Theil-Sen slope (Sen, 1968) with its 95 % interval, original Mann-Kendall test and the Hamed and Rao (1998)
#   variance correction for autocorrelation (package `pymannkendall`); slopes also in % of the mean per year;
# * per-zone tables add the cells whose centre lies in each zone (notebook 10, `stats/municipality_<year>_<mm>.csv`).

# %%
__FILE:notebooks/src/_setup_repo.py__
%pip -q install pymannkendall rasterio geopandas netCDF4 scipy pandas matplotlib pyyaml
RUN_CONFIG = 'config/santa_catarina_2003_2025.yaml'           # <-- the configuration of the run (as in notebooks 10 and 11)
OUT_BASE   = '/content/drive/MyDrive' if IN_COLAB else os.path.dirname(REPO_DIR)   # <-- as in notebook 10
RUN_ROOT   = None                                             # None = OUT_BASE/<output.dir of RUN_CONFIG>; or the run folder
OUT_DIR    = None                                             # None = RUN_ROOT/results
RHO_ISOP   = 0.96                                             # isoprene canopy loss and production factor (Guenther et al., 2006)
MIN_YEARS_STATS = 10                                          # trends and correlations need at least this many complete years

# %%
# 1 · Configuration, grid, region mask and cell areas ---------------------------------------------------------------
import json, re, warnings, yaml
import numpy as np, pandas as pd, geopandas as gpd, rasterio, matplotlib.pyplot as plt
from netCDF4 import Dataset
from rasterio.features import rasterize
from rasterio.transform import from_origin
from scipy import stats
from megan_sc import domain as dom
from megan_sc.io import geotiff
warnings.filterwarnings('ignore', category=RuntimeWarning)
def in_repo(p): return p if (not p or os.path.isabs(p)) else os.path.join(REPO_DIR, p)
CFG = yaml.safe_load(open(in_repo(RUN_CONFIG), encoding='utf-8'))
REG = CFG.get('region', {}); REGION_NAME = REG.get('name', CFG['run']['name'])
RUN_ROOT = RUN_ROOT or os.path.join(OUT_BASE, CFG.get('output', {}).get('dir', f"megan_sc_outputs_{CFG['run']['name']}"))
OUT_DIR = OUT_DIR or os.path.join(RUN_ROOT, 'results'); os.makedirs(OUT_DIR, exist_ok=True)
BOUNDARY = in_repo(REG.get('boundary') or CFG['domain'].get('shapefile'))
CLASSES = ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM', 'NO', 'SQT_HR', 'SQT_LR',
           'MEOH', 'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO']
BVOC = [c for c in CLASSES if c not in ('NO', 'CO')]
MT = [c for c in CLASSES if c.startswith('MT_')]
def groups(g):
    o = {'ISOP': g['ISOP'], 'MBO': g['MBO'], 'MT': sum(g[k] for k in MT), 'SQT': g['SQT_HR'] + g['SQT_LR'],
         'OVOC': sum(g[k] for k in ('MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD')), 'LVOC': g['LVOC'],
         'STRESS_OTHER': g['STRESS'] + g['OTHER'], 'NO': g['NO'], 'CO': g['CO']}
    o['BVOC'] = sum(g[k] for k in BVOC); o['SOA'] = o['MT'] + o['SQT']
    return o
def monthly_files(y):
    """The monthly emission files of year y (test runs excluded)."""
    return sorted(glob.glob(os.path.join(RUN_ROOT, str(y), 'emissions', f'emis_{y}_[01][0-9].nc')))
def read_nc(fn):
    """Per-class grids (t per cell) of one file, south-first; isoprene times RHO_ISOP."""
    with Dataset(fn) as ds:
        lat = np.array(ds.variables['lat'][:]); lon = np.array(ds.variables['lon'][:])
        E = {}
        for c in CLASSES:
            a = np.array(ds.variables[c][:], dtype='float64'); E[c] = np.where(np.isfinite(a) & (a > -9990), a, 0.0)
    E['ISOP'] = E['ISOP'] * RHO_ISOP
    return E, lat, lon
YEARS_FOUND = sorted(int(os.path.basename(p)) for p in glob.glob(os.path.join(RUN_ROOT, '[12][0-9][0-9][0-9]')) if monthly_files(int(os.path.basename(p))))
assert YEARS_FOUND, f'no monthly emission files (emis_<year>_<mm>.nc) under {RUN_ROOT}'
COMPLETE = [y for y in YEARS_FOUND if len(monthly_files(y)) == 12]
_, lat, lon = read_nc(monthly_files(YEARS_FOUND[0])[0])
NY, NX = lat.size, lon.size; RES = float(lon[1] - lon[0])
TR = from_origin(lon[0] - RES / 2, lat[-1] + RES / 2, RES, RES)
_g = gpd.read_file(BOUNDARY).to_crs('EPSG:4326')
_geom = _g.geometry.union_all() if hasattr(_g.geometry, 'union_all') else _g.geometry.unary_union
MASK = np.flipud(rasterize([(_geom, 1)], out_shape=(NY, NX), transform=TR, fill=0, dtype='uint8')).astype(bool)
AREA = dom.cell_area_grid(lat, RES, NX)                    # m2, as in the engine
print(f'{REGION_NAME}: {int(MASK.sum())} cells, {AREA[MASK].sum() / 1e6:,.0f} km2 of grid | years with output {YEARS_FOUND[0]}-'
      f'{YEARS_FOUND[-1]} | complete years {len(COMPLETE)} | results -> {OUT_DIR}')

# %%
# 2 · Monthly and annual totals of the region --------------------------------------------------------------------
rows, ANN = [], {}
for y in YEARS_FOUND:
    acc = {c: np.zeros((NY, NX)) for c in CLASSES}
    for fn in monthly_files(y):
        m = int(re.search(r'emis_\d{4}_(\d{2})\.nc$', fn).group(1))
        E, _, _ = read_nc(fn)
        rows.append({'year': y, 'month': m, **{c: float((E[c] * MASK).sum()) for c in CLASSES}})
        for c in CLASSES: acc[c] += E[c]
    if y in COMPLETE:
        ANN[y] = acc
MON = pd.DataFrame(rows).set_index(['year', 'month']).sort_index()
MON.to_csv(os.path.join(OUT_DIR, 'monthly_totals_by_class.csv'))
CL = pd.DataFrame({y: {c: float((ANN[y][c] * MASK).sum()) for c in CLASSES} for y in COMPLETE}).T; CL.index.name = 'year'
G = pd.DataFrame({y: {k: float((v * MASK).sum()) for k, v in groups(ANN[y]).items()} for y in COMPLETE}).T; G.index.name = 'year'
CL.to_csv(os.path.join(OUT_DIR, 'annual_totals_by_class.csv')); G.to_csv(os.path.join(OUT_DIR, 'annual_totals_by_group.csv'))
SHOW = ['BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'NO', 'CO']
if len(G):
    print('annual totals of the region (t/yr):'); print(G[SHOW].round(0).astype(int).to_string())
else:
    print('no complete year: annual tables, maps, trends and correlations are skipped; the monthly table is written')

# %%
# 3 · Seasonal cycle -----------------------------------------------------------------------------------------------
MG = pd.DataFrame([{**{'year': y, 'month': m}, **groups(MON.loc[(y, m)].to_dict())} for y, m in MON.index]).set_index(['year', 'month'])
use = MG.loc[COMPLETE] if COMPLETE else MG
CLIM = use.groupby(level='month').agg(['mean', 'min', 'max'])
CLIM.to_csv(os.path.join(OUT_DIR, 'seasonal_climatology.csv'))
fig, axs = plt.subplots(1, 4, figsize=(15, 3.4))
for ax, k in zip(axs, ['ISOP', 'MT', 'SQT', 'NO']):
    ax.fill_between(CLIM.index, CLIM[(k, 'min')] / 1e3, CLIM[(k, 'max')] / 1e3, color='#9ecae1', alpha=.6, label='min-max')
    ax.plot(CLIM.index, CLIM[(k, 'mean')] / 1e3, 'o-', color='#08519c', ms=3, label='mean')
    ax.set_title(k); ax.set_xticks(range(1, 13)); ax.set_xlabel('month (UTC)'); ax.set_ylabel('kt month⁻¹'); ax.grid(alpha=.3)
axs[0].legend(fontsize=8); fig.suptitle(f'{REGION_NAME}: seasonal cycle ({len(COMPLETE) or "no complete"} years)'); fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'fig_seasonal_cycle.png'), dpi=200); plt.show()

# %%
# 4 · Annual series and mean maps of intensity ----------------------------------------------------------------------
if len(G):
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.4))
    for ax, k in zip(axs, ['BVOC', 'ISOP', 'MT', 'NO']):
        ax.plot(G.index, G[k] / 1e3, 'o-', ms=3, color='#1b6ca8'); ax.set_title(k); ax.set_ylabel('kt yr⁻¹'); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, 'fig_annual_series.png'), dpi=200); plt.show()
    MEAN = {k: np.mean([groups(ANN[y])[k] for y in COMPLETE], axis=0) for k in ('BVOC', 'ISOP', 'MT', 'SQT', 'NO')}
    fig, axs = plt.subplots(1, 5, figsize=(19, 4.2))
    ext = [lon[0] - RES / 2, lon[-1] + RES / 2, lat[0] - RES / 2, lat[-1] + RES / 2]
    for ax, (k, v) in zip(axs, MEAN.items()):
        inten = np.where(MASK, v * 1e3 / (AREA / 1e4), np.nan)          # t -> kg; m2 -> ha
        geotiff.write_geotiff(os.path.join(OUT_DIR, f'mean_intensity_{k}.tif'), inten, lon, lat, band_names=[k], nodata=np.nan)
        im = ax.imshow(np.ma.masked_invalid(inten), origin='lower', extent=ext, cmap='YlGn' if k != 'NO' else 'OrRd',
                       aspect=1.0 / np.cos(np.deg2rad(float(np.mean(lat)))))
        for gm in getattr(_geom, 'geoms', [_geom]):
            try: xs, ys = gm.exterior.xy; ax.plot(xs, ys, color='0.3', lw=0.4)
            except Exception: pass
        ax.set_title(f'{k}, kg ha⁻¹ yr⁻¹'); ax.set_xticks([]); ax.set_yticks([]); fig.colorbar(im, ax=ax, shrink=0.75)
    fig.suptitle(f'{REGION_NAME}: mean annual intensity, {COMPLETE[0]}-{COMPLETE[-1]}'); fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'fig_maps_mean_intensity.png'), dpi=200); plt.show()

# %%
# 5 · Drivers as the model saw them, correlations and trends (10 complete years or more) ----------------------------
def regional_mean(y, kind, band, positive=False):
    """Regional mean of the annual mean of the monthly input maps written by notebook 10 (inputs/<kind>_<y>_<mm>.tif)."""
    files = sorted(glob.glob(os.path.join(RUN_ROOT, str(y), 'inputs', f'{kind}_{y}_[01][0-9].tif')))
    if not files: return np.nan
    acc = np.zeros((NY, NX)); cnt = np.zeros((NY, NX))
    for f in files:
        with rasterio.open(f) as s: b = np.flipud(s.read(band).astype('float64'))
        ok = np.isfinite(b); acc[ok] += b[ok]; cnt[ok] += 1
    m = np.where(cnt > 0, acc / np.where(cnt > 0, cnt, 1), np.nan)
    if positive: m = np.where(m > 0, m, np.nan)
    return float(np.nanmean(np.where(MASK, m, np.nan)))
def detr(v):
    x = np.arange(len(v)); v = np.asarray(v, float); return v - np.polyval(np.polyfit(x, v, 1), x)
if len(COMPLETE) >= MIN_YEARS_STATS:
    import pymannkendall as pmk
    DRV = pd.DataFrame({y: {'TEMP2_C': regional_mean(y, 'met', 1) - 273.15, 'PAR': regional_mean(y, 'met', 2),
                            'SOILM': regional_mean(y, 'met', 6), 'LAI': regional_mean(y, 'lai', 1, positive=True)}
                        for y in COMPLETE}).T
    DRV.index.name = 'year'; DRV.to_csv(os.path.join(OUT_DIR, 'drivers_annual.csv'))
    rows = []
    for k in ('BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'NO', 'CO'):
        for d in DRV.columns:
            if DRV[d].notna().all():
                r, p = stats.pearsonr(G[k].values, DRV[d].values); rd, pdd = stats.pearsonr(detr(G[k].values), detr(DRV[d].values))
                rows.append(dict(series=k, driver=d, r=r, p=p, r_detrended=rd, p_detrended=pdd))
    COR = pd.DataFrame(rows); COR.to_csv(os.path.join(OUT_DIR, 'correlations_emissions_drivers.csv'), index=False)
    def trend(s):
        yv = np.asarray(s, float); x = np.arange(len(yv)); o = pmk.original_test(yv); hr = pmk.hamed_rao_modification_test(yv)
        sl, ic, lo, hi = stats.theilslopes(yv, x, 0.95)
        return dict(mean=yv.mean(), sen=sl, pct_yr=sl / yv.mean() * 100, pct_lo95=lo / yv.mean() * 100,
                    pct_hi95=hi / yv.mean() * 100, MK_Z=o.z, MK_p=o.p, HR_p=hr.p)
    TRD = pd.DataFrame({k: trend(G[k]) for k in ('BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'MBO', 'NO', 'CO')}).T
    TRD.index.name = 'series'; TRD.to_csv(os.path.join(OUT_DIR, 'trends.csv'))
    print('trends (% of the mean per year; Theil-Sen, 95 % interval; Mann-Kendall and Hamed-Rao p):')
    print(TRD[['pct_yr', 'pct_lo95', 'pct_hi95', 'MK_p', 'HR_p']].round(3).to_string())
    print('\ndetrended correlations with the drivers:'); print(COR.pivot(index='series', columns='driver', values='r_detrended').round(2).to_string())
else:
    print(f'{len(COMPLETE)} complete year(s): trends and correlations need {MIN_YEARS_STATS} (MIN_YEARS_STATS); skipped')

# %%
# 6 · Zones (when notebook 10 wrote per-zone tables) --------------------------------------------------------------
zf = sorted(f for y in COMPLETE for f in glob.glob(os.path.join(RUN_ROOT, str(y), 'stats', f'municipality_{y}_[01][0-9].csv')))
if zf:
    Z = pd.concat([pd.read_csv(f) for f in zf], ignore_index=True)
    Z.loc[Z.species == 'ISOP', 'tonnes'] *= RHO_ISOP
    ZA = Z.groupby(['year', 'municipality', 'species']).tonnes.sum().unstack('species')
    ZA = ZA.assign(BVOC=ZA[BVOC].sum(axis=1)).rename_axis(index={'municipality': 'zone'})
    ZA.to_csv(os.path.join(OUT_DIR, 'zones_annual_totals.csv'))
    print(f'zones_annual_totals.csv: {ZA.index.get_level_values("zone").nunique()} zones x {len(COMPLETE)} years')
else:
    print('no per-zone tables in the run (region.zones not set in the configuration)')
print('\nwritten to', OUT_DIR, ':', sorted(os.listdir(OUT_DIR)))
