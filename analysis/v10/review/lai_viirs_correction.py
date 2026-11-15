"""Design, validation and result of the VIIRS-anomaly correction of the 2024-2025 MODIS LAI (post hoc).

Inputs
  ROOT/<year>/inputs/lai_<year>_<mm>.tif   monthly LAI used by the production run (notebook 10), 2013-2025
  ROOT/<year>/inputs/growthform_*.npz      MapBiomas growth-form fractions of the run (cell groups, vegetated cells)
  ROOT/2025/emissions/annual_emis_2025.nc  grid (cell centres, south first)
  CACHE  lai_viirs_modis_cache.npz         monthly VIIRS LAI on the same grid, written by notebook 16
  BND/sc_state_ibge.gpkg                   Santa Catarina polygon (cell-centre rule, as notebooks 10, 11 and 16)

Steps
  1. Reference distribution of the cell ratios r = V(c, y, m) / V(c, m) in 2013-2023 (clip limits).
  2. Leave-one-year-out hindcast: each year Y of 2013-2023 is predicted from the other ten years with every candidate
     form of the correction (src/megan_sc/laicorr.py) and compared with the MODIS LAI used by the run in Y.
     Selection rule, fixed before the comparison: the candidate with the lowest root-mean-square error over all
     Santa Catarina cell-months; the baseline 'clim' (MODIS reference mean only) shows what VIIRS adds.
  3. Corrected LAI of every month of 2024 and 2025 with the selected candidate and the eleven reference years.
  4. Tables and figures in OUT.

Usage: python analysis/v10/review/lai_viirs_correction.py [ROOT] [CACHE] [OUT]
Defaults (analysis/v10/_paths.py): MEGAN_SC_RUN (data/run_v10), MEGAN_SC_LAIVIIRS/lai_viirs_modis_cache.npz
(data/lai_viirs_v10) and <MEGAN_SC_OUT>/lai_viirs_check; published copy in results/v10/lai_viirs_check
(`python reproduce.py laicorr`).
"""
import os, sys, glob, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'src')); sys.path.insert(0, os.path.join(REPO, 'analysis', 'v10'))
import _paths as PTH
from megan_sc import laicorr as lc
import rasterio, geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import Affine
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = sys.argv[1] if len(sys.argv) > 1 else PTH.RUN
CACHE = sys.argv[2] if len(sys.argv) > 2 else os.path.join(PTH.LAIVIIRS, 'lai_viirs_modis_cache.npz')
OUT = sys.argv[3] if len(sys.argv) > 3 else PTH.out('lai_viirs_check')
BND = os.path.join(PTH.BND, 'sc_state_ibge.gpkg')
os.makedirs(OUT, exist_ok=True)
REF = list(lc.REF_YEARS); TARGET = [2024, 2025]; YEARS = REF + TARGET; MONTHS = range(1, 13); DOMINANT = 0.6

# ---------------------------------------------------------------- grid, mask, inputs
with Dataset(os.path.join(ROOT, '2025', 'emissions', 'annual_emis_2025.nc')) as ds:
    lat = np.array(ds.variables['lat'][:]); lon = np.array(ds.variables['lon'][:])
res = float(lon[1] - lon[0]); assert lat[0] < lat[-1] and abs(res - 0.025) < 1e-9
tr = Affine(res, 0.0, float(lon[0]) - res / 2, 0.0, -res, float(lat[-1]) + res / 2)
SC = np.flipud(rasterize([(g, 1) for g in gpd.read_file(BND).to_crs('EPSG:4326').geometry], out_shape=(lat.size, lon.size),
                         transform=tr, fill=0, dtype='uint8')).astype(bool)

def tif(path):
    with rasterio.open(path) as s:
        return np.flipud(s.read(1).astype('float64'))                     # south first, as the run arrays

M = {(y, m): tif(os.path.join(ROOT, str(y), 'inputs', f'lai_{y}_{m:02d}.tif')) for y in YEARS for m in MONTHS}
C = np.load(CACHE)
V = {(y, m): C[f'VIIRS_{y}_{m:02d}_LAI'].astype('float64') for y in YEARS for m in MONTHS}
M16 = {(y, m): C[f'MODIS_{y}_{m:02d}_LAI'].astype('float64') for y in YEARS for m in MONTHS}

def gf(y):
    f = sorted(glob.glob(os.path.join(ROOT, str(y), 'inputs', 'growthform_*.npz')))[0]
    return dict(np.load(f))
GF = {y: gf(y) for y in YEARS}
def groups(y):
    g = GF[y]
    return {'all': SC, 'tree': SC & (g['tree'] > DOMINANT), 'herb': SC & (g['herb'] > DOMINANT), 'crop': SC & (g['crop'] > DOMINANT)}
VEG = {y: sum(GF[y][k] for k in ('tree', 'shrub', 'herb', 'crop') if k in GF[y]) > 0 for y in YEARS}

# the notebook-16 MODIS samples reproduce the run's LAI (same product, mask and aggregation; the run also applies
# the vegetated-cell mask and the FIX-4 gap fill)
both = [(k, (M[k] > 0) & (M16[k] > 0)) for k in M]
dmax = max(float(np.abs(M[k][b] - M16[k][b]).max()) for k, b in both if b.any())
assert dmax < 1e-5, dmax

# ---------------------------------------------------------------- 1. reference distribution of the cell ratios
rr = []
for m in MONTHS:
    vref, _ = lc.monthly_reference([V[(y, m)] for y in REF])
    for y in REF:
        v = lc.valid(V[(y, m)]); ok = SC & np.isfinite(v) & np.isfinite(vref)
        rr.append(v[ok] / vref[ok])
rr = np.concatenate(rr)
P = {q: float(np.percentile(rr, q)) for q in (0.5, 1, 5, 50, 95, 99, 99.5)}
CLIP = (P[0.5], P[99.5])
print('reference cell ratios: n', rr.size, '| percentiles', {k: round(v, 3) for k, v in P.items()})

# ---------------------------------------------------------------- 2. leave-one-year-out hindcast
MESO_GPKG = os.path.join(PTH.BND, 'sc_mesoregions_2024.gpkg')
_meso = gpd.read_file(MESO_GPKG).to_crs('EPSG:4326')
MESO = np.flipud(rasterize([(g, i + 1) for i, g in enumerate(_meso.geometry)], out_shape=(lat.size, lon.size), transform=tr,
                           fill=0, dtype='int32')) * SC                 # 1-6, cell-centre rule
def gform_zones(y):
    return lc.growth_form_zones(GF[y], SC, DOMINANT)                   # 1 tree, 2 herb, 3 crop, 4 other SC cells
ZONES = {'meso': lambda y: MESO, 'gform': gform_zones}
# Candidate list. The first seven were fixed before any comparison. After the first comparison showed that the cell
# and 3 x 3 ratios add noise (they did worse than the MODIS reference mean alone) and that the statewide ratio did
# best, the list was extended once with intermediate spatial scales (11 x 11 and 21 x 21 windows, mesoregions, and
# the dominant growth form); the same selection rule applies to the extended list.
CANDIDATES = [dict(name='clim', method='clim', k=0, clip=None), dict(name='state', method='state', k=0, clip=None),
              dict(name='cell', method='cell', k=0, clip=None), dict(name='cell_clip', method='cell', k=0, clip=CLIP),
              dict(name='smooth3', method='smooth', k=1, clip=None), dict(name='smooth3_clip', method='smooth', k=1, clip=CLIP),
              dict(name='smooth5', method='smooth', k=2, clip=None),
              dict(name='smooth11', method='smooth', k=5, clip=None), dict(name='smooth21', method='smooth', k=10, clip=None),
              dict(name='mesoregion', method='zone', k=0, clip=None, zones='meso'),
              dict(name='growth_form', method='zone', k=0, clip=None, zones='gform')]
rows, arows = [], []
for cnd in CANDIDATES:
    sse = sae = se = n = 0.0; an_p, an_t = [], []; nfb = ncl = nmf = 0; ncells = 0
    for Y in REF:
        others = [y for y in REF if y != Y]
        P_, T_ = np.full((12,) + SC.shape, np.nan), np.full((12,) + SC.shape, np.nan)
        for m in MONTHS:
            mref, _ = lc.monthly_reference([M[(y, m)] for y in others])
            lai, info = lc.correct_month([M[(y, m)] for y in others], [V[(y, m)] for y in others], M[(Y, m)], V[(Y, m)],
                                         ref_mask=SC, method=cnd['method'], k=cnd['k'], clip=cnd['clip'],
                                         zones=ZONES[cnd['zones']](Y) if 'zones' in cnd else None)
            t = lc.valid(M[(Y, m)]); E = SC & np.isfinite(t) & np.isfinite(mref)    # cells with a MODIS reference mean
            d = lai[E] - t[E]
            sse += float((d ** 2).sum()); sae += float(np.abs(d).sum()); se += float(d.sum()); n += d.size
            an_p.append(lai[E] / mref[E] - 1); an_t.append(t[E] / mref[E] - 1)
            nfb += info['ratio_fallback_cells']; ncl += info['clipped_cells']; nmf += info['modis_fallback_cells']; ncells += int(SC.sum())
            P_[m - 1][E] = lai[E]; T_[m - 1][E] = t[E]
        pa, ta = np.nanmean(P_, 0), np.nanmean(T_, 0)                       # annual mean per cell (same cell-months)
        for gname, msk in groups(Y).items():
            ok = msk & np.isfinite(pa) & np.isfinite(ta)
            arows.append(dict(candidate=cnd['name'], year=Y, group=gname, true_lai=float(ta[ok].mean()),
                              pred_lai=float(pa[ok].mean()), error_pct=100 * (pa[ok].mean() / ta[ok].mean() - 1)))
    ap, at = np.concatenate(an_p), np.concatenate(an_t)
    rows.append(dict(candidate=cnd['name'], method=cnd['method'], window=(2 * cnd['k'] + 1) if cnd['method'] == 'smooth' else '',
                     zones=cnd.get('zones', ''),
                     clip='' if cnd['clip'] is None else f'{cnd["clip"][0]:.3f}-{cnd["clip"][1]:.3f}',
                     cell_months=int(n), rmse=float(np.sqrt(sse / n)), mae=float(sae / n), bias=float(se / n),
                     anomaly_r=float(np.corrcoef(ap, at)[0, 1]), ratio_fallback_pct=100 * nfb / ncells,
                     clipped_pct=100 * ncl / ncells, modis_fallback_pct=100 * nmf / ncells))
H = pd.DataFrame(rows); HA = pd.DataFrame(arows)
base = float(H.loc[H.candidate == 'clim', 'rmse'].iloc[0])
H['skill_vs_clim'] = 1 - H.rmse / base
for c in ('all', 'tree', 'herb', 'crop'):
    s = HA[HA.group == c].groupby('candidate').error_pct
    H[f'annual_abs_err_{c}_mean_pct'] = H.candidate.map(s.apply(lambda x: float(np.abs(x).mean())))
    H[f'annual_abs_err_{c}_max_pct'] = H.candidate.map(s.apply(lambda x: float(np.abs(x).max())))
H.to_csv(os.path.join(OUT, 'lai_correction_hindcast_2013_2023.csv'), index=False)
HA.to_csv(os.path.join(OUT, 'lai_correction_hindcast_annual.csv'), index=False)
best = H.sort_values('rmse').iloc[0]; CHOICE = [c for c in CANDIDATES if c['name'] == best.candidate][0]
print(H[['candidate', 'rmse', 'mae', 'bias', 'anomaly_r', 'skill_vs_clim', 'annual_abs_err_all_mean_pct',
         'annual_abs_err_all_max_pct', 'ratio_fallback_pct', 'clipped_pct']].round(4).to_string(index=False))
print('selected:', CHOICE)

# ---------------------------------------------------------------- 3. corrected LAI of 2024 and 2025
# the module's DEFAULT must be the selected candidate; the fields come from laicorr.corrected_run_lai, the function
# that notebook 17 calls, so the published fields and the rerun share one code path
assert (CHOICE['method'], CHOICE['k'], CHOICE['clip'], CHOICE.get('zones')) == \
       (lc.DEFAULT['method'], lc.DEFAULT['k'], lc.DEFAULT['clip'], 'gform'), ('update laicorr.DEFAULT', CHOICE)
CORR, INFOS = lc.corrected_run_lai(ROOT, CACHE, SC, years=TARGET)
FIELDS, mrows = {}, []
for Y in TARGET:
    for m in MONTHS:
        lai = CORR[(Y, m)]; info = [i for i in INFOS if i['year'] == Y and i['month'] == m][0]
        FIELDS[f'LAI_{Y}_{m:02d}'] = lai.astype('float32')
        run = lai * VEG[Y]                                                     # as lai_raw in notebook 10 (vegetated cells)
        o, c = lc.valid(M[(Y, m)]), lc.valid(run)
        mref, _ = lc.monthly_reference([M[(y, m)] for y in REF])
        sel = SC & np.isfinite(o) & np.isfinite(c) & np.isfinite(mref)
        mrows.append(dict(year=Y, month=m, modis_lai=float(o[sel].mean()), corrected_lai=float(c[sel].mean()),
                          modis_dev_pct=100 * (o[sel].mean() / mref[sel].mean() - 1), corrected_dev_pct=100 * (c[sel].mean() / mref[sel].mean() - 1),
                          viirs_state_ratio_dev_pct=100 * (info['state_ratio'] - 1), cells=int(sel.sum()),
                          ratio_fallback_cells=info['ratio_fallback_cells'], clipped_cells=info['clipped_cells'],
                          modis_fallback_cells=info['modis_fallback_cells'],
                          corrected_over_modis_p05=float(np.percentile(c[sel] / o[sel], 5)),
                          corrected_over_modis_p50=float(np.percentile(c[sel] / o[sel], 50)),
                          corrected_over_modis_p95=float(np.percentile(c[sel] / o[sel], 95))))
MT = pd.DataFrame(mrows); MT.to_csv(os.path.join(OUT, 'lai_correction_2024_2025_monthly.csv'), index=False)
np.savez_compressed(os.path.join(OUT, 'lai_viirs_corrected_2024_2025.npz'), lat=lat, lon=lon,
                    method=np.array(CHOICE['method']), window=np.array(2 * CHOICE['k'] + 1),
                    clip=np.array([np.nan, np.nan] if CHOICE['clip'] is None else list(CHOICE['clip'])),
                    ref_years=np.array(REF), **FIELDS)

# annual deviation of each group from its 2013-2023 mean: MODIS as used, VIIRS, corrected
def annual(stack_fn, y):
    A = np.stack([lc.valid(stack_fn(y, m)) for m in MONTHS]); return np.nanmean(A, 0)
ann = []
for y in YEARS:
    for gname, msk in groups(y).items():
        for sensor, fn in (('MODIS', lambda y_, m_: M[(y_, m_)]), ('VIIRS', lambda y_, m_: V[(y_, m_)]),
                           ('corrected', lambda y_, m_: FIELDS[f'LAI_{y_}_{m_:02d}'] * VEG[y_] if y_ in TARGET else M[(y_, m_)])):
            a = annual(fn, y); ok = msk & np.isfinite(a)
            ann.append(dict(sensor=sensor, group=gname, year=y, lai=float(a[ok].mean())))
AN = pd.DataFrame(ann)
out = []
for (s, g), d in AN.groupby(['sensor', 'group']):
    d = d.set_index('year'); ref = d.loc[REF, 'lai']
    for y in YEARS:
        out.append(dict(sensor=s, group=g, year=y, lai=d.loc[y, 'lai'], dev_pct=100 * (d.loc[y, 'lai'] / ref.mean() - 1),
                        z=(d.loc[y, 'lai'] - ref.mean()) / ref.std(ddof=1)))
AD = pd.DataFrame(out); AD.to_csv(os.path.join(OUT, 'lai_correction_annual_groups.csv'), index=False)
SUM = AD[AD.year.isin(TARGET)].pivot_table(index=['group', 'year'], columns='sensor', values='dev_pct').round(1)
print('annual deviation from 2013-2023 (%):'); print(SUM.to_string())

# ---------------------------------------------------------------- 4. figures
plt.rcParams.update({'font.size': 8})
fig, ax = plt.subplots(3, 1, figsize=(7.0, 8.4))
mon = pd.date_range('2013-01-01', '2025-12-01', freq='MS')
def state_anom(get, ref_get):
    out = []
    for t in mon:
        y, m = t.year, t.month
        a = lc.valid(get(y, m)); r_, _ = lc.monthly_reference([ref_get(yy, m) for yy in REF])
        ok = SC & np.isfinite(a) & np.isfinite(r_); out.append(100 * (a[ok].mean() / r_[ok].mean() - 1))
    return np.array(out)
mo = state_anom(lambda y, m: M[(y, m)], lambda y, m: M[(y, m)])
vi = state_anom(lambda y, m: V[(y, m)], lambda y, m: V[(y, m)])
co = state_anom(lambda y, m: FIELDS[f'LAI_{y}_{m:02d}'] * VEG[y] if y in TARGET else M[(y, m)], lambda y, m: M[(y, m)])
sel = mon >= '2022-01-01'
ax[0].plot(mon[sel], mo[sel], color='#1b9e77', lw=1.4, label='MODIS, as used by the run')
ax[0].plot(mon[sel], vi[sel], color='#d95f02', lw=1.4, label='VIIRS')
ax[0].plot(mon[mon >= '2024-01-01'], co[mon >= '2024-01-01'], color='k', lw=1.6, ls='--', label='corrected (2024-2025)')
ax[0].axhline(0, color='0.6', lw=0.7); ax[0].set_ylabel('% of the 2013-2023 mean\nof the month'); ax[0].legend(frameon=False, loc='lower left')
ax[0].set_title('(a) Statewide monthly LAI anomaly, 2022-2025', loc='left', fontweight='bold')
G = ['all', 'tree', 'herb', 'crop']; lab = ['all cells', 'tree-dominated', 'herbaceous-dominated', 'crop-dominated']
x = np.arange(len(G)); w = 0.13
for i, (yy, s, col, h) in enumerate([(2024, 'MODIS', '#1b9e77', ''), (2024, 'VIIRS', '#d95f02', ''), (2024, 'corrected', '0.35', ''),
                                     (2025, 'MODIS', '#1b9e77', '//'), (2025, 'VIIRS', '#d95f02', '//'), (2025, 'corrected', '0.35', '//')]):
    vals = [AD[(AD.sensor == s) & (AD.group == g) & (AD.year == yy)].dev_pct.iloc[0] for g in G]
    ax[1].bar(x + (i - 2.5) * w, vals, w, color=col, hatch=h, edgecolor='white' if not h else 'k', lw=0.3, label=f'{s} {yy}')
ax[1].axhline(0, color='0.6', lw=0.7); ax[1].set_ylim(-17.5, 3.5); ax[1].set_xticks(x); ax[1].set_xticklabels(lab); ax[1].set_ylabel('annual deviation (%)')
ax[1].legend(frameon=False, ncol=2, fontsize=7, loc='lower left')
ax[1].set_title('(b) Annual LAI of 2024 and 2025 versus 2013-2023', loc='left', fontweight='bold')
order = H.sort_values('rmse')
ax[2].barh(order.candidate, order.rmse, color=['k' if c == CHOICE['name'] else '0.7' for c in order.candidate])
for i, (r_, s_) in enumerate(zip(order.rmse, order.skill_vs_clim)):
    ax[2].text(r_ + 0.005, i, f'{r_:.3f} (skill {s_:+.2f})', va='center', fontsize=7)
ax[2].set_xlabel('root-mean-square error of the cell LAI (m² m⁻²), leave-one-year-out, 2013-2023')
ax[2].set_xlim(0, float(order.rmse.max()) * 1.35)
ax[2].set_title('(c) Hindcast of the correction (selected: black)', loc='left', fontweight='bold')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'lai_correction.png'), dpi=300); plt.close(fig)

a25o = annual(lambda y, m: M[(y, m)], 2025); a25c = annual(lambda y, m: FIELDS[f'LAI_{y}_{m:02d}'] * VEG[y], 2025)
dmap = np.where(SC & np.isfinite(a25o) & np.isfinite(a25c), 100 * (a25c / a25o - 1), np.nan)
fig, ax = plt.subplots(figsize=(6.2, 4.4))
cm = plt.get_cmap('BrBG').copy(); cm.set_bad('white')
lim = float(np.nanpercentile(np.abs(dmap), 99))
im = ax.imshow(np.ma.masked_invalid(dmap), origin='lower', extent=[lon[0] - res / 2, lon[-1] + res / 2, lat[0] - res / 2, lat[-1] + res / 2],
               cmap=cm, vmin=-lim, vmax=lim, aspect=1 / np.cos(np.deg2rad(lat.mean())))
for g in gpd.read_file(BND).to_crs('EPSG:4326').geometry:
    for p in getattr(g, 'geoms', [g]):
        xs, ys = p.exterior.xy; ax.plot(xs, ys, color='0.3', lw=0.5)
ax.set_xlim(-54.0, -48.2); ax.set_ylim(-29.5, -25.8); ax.set_xlabel('longitude'); ax.set_ylabel('latitude')
cb = fig.colorbar(im, ax=ax, shrink=0.8); cb.set_label('corrected / MODIS - 1 (%)')
ax.set_title('Annual mean LAI of 2025: corrected relative to MODIS', loc='left', fontweight='bold', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'lai_correction_map_2025.png'), dpi=300); plt.close(fig)

json.dump(dict(reference_years=REF, target_years=TARGET, min_frac=lc.MIN_FRAC, clip_percentiles=[0.5, 99.5],
               clip_limits=list(CLIP), ratio_percentiles=P, selected=dict(name=CHOICE['name'], method=CHOICE['method'], zones=CHOICE.get('zones'),
               window=2 * CHOICE['k'] + 1, clip=None if CHOICE['clip'] is None else list(CHOICE['clip'])),
               modis_cache_vs_run_max_abs_diff=dmax),
          open(os.path.join(OUT, 'lai_correction_settings.json'), 'w'), indent=1)
print('written to', OUT)
