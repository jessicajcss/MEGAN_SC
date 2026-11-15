# %% [markdown]
# # 13 · Every number and figure of the manuscript, from the production run
#
# Computes, from the outputs of `10_production_run.ipynb` (after `11_verify_run.ipynb` passes), the numbers of the
# manuscript, Tables 2-4 and 6, Figures 1-8 and the example of the detrended correlation (`fig_detrend.png`, article
# Figure 2) in Portuguese and English.
# Nothing is typed in by hand: the only external numbers are the three 2019 totals of the previous state inventory
# (Hoinaski *et al.*, 2022, LCQAr report, Table 12).
#
# Estimators:
# * isoprene: the canopy loss and production factor of Guenther et al. (2006, Sect. 3.3), rho = 0.96 (RHO_ISOP, cell 1),
#   multiplies the isoprene of notebook 10 when it is read; the MEGANv3.21 code applies no such factor (rho = 1);
# * series statistics on state totals (SC mask, cell-centre rule); CV with ddof = 1;
# * trends: Theil-Sen slope with its 95 % CI (scipy), Mann-Kendall p and the Hamed-Rao autocorrelation-corrected p
#   (pymannkendall); slopes in % of the 2003-2025 mean per year;
# * emission-driver correlations: Pearson r of the raw and of the linearly detrended annual series (least-squares line
#   against the year removed from each series; `fig_detrend.png` shows the procedure for soil NO and temperature); drivers are the SC means
#   of the annual-mean monthly maps written by notebook 10 (TEMP2 in °C, PAR, SOILM; LAI averaged over cells with
#   LAI > 0);
# * mesoregions and municipalities: the emission of every grid cell is split among the polygons in proportion to the area
#   of the cell inside each one (area weighting), so the area attributed to a municipality is its polygon area; intensity
#   per IBGE area (for a mesoregion, the sum of the IBGE areas of its municipalities);
# * mesoregions: 2003-2025 mean of MT + SQT (SOA precursors); Pearson r with Fisher 95 % CI and Spearman ρ against the
#   native-forest attribute of the mesoregion file (n = 6, descriptive only);
# * municipalities: 2025 totals; Jenks natural breaks (Fisher's exact algorithm) in the maps;
# * seasonal cycle: monthly state totals by class from the stats_<year>_<mm>.csv files of notebook 10 (Figure 2b, 2c);
# * isoprene areal flux in g C m⁻² yr⁻¹ (5 C atoms, 68.117 g mol⁻¹).
#
# **Send back:** the folder `OUT_DIR` (CSV tables, `manuscript_numbers.csv`, figures `pt/` and `en/`).

# %%
# 1 · Paths and settings ---------------------------------------------------------------------------------
import os, glob, json, warnings
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
%pip -q install pymannkendall rasterio geopandas netCDF4 scipy pandas matplotlib
ROOT    = '/content/drive/MyDrive/megan_sc_outputs_v10'          # <-- production run (notebook 10); cell 2 searches MyDrive if absent
BND     = '/content/drive/MyDrive/megan-sc/data/boundaries'      # <-- sc_state_ibge.gpkg, sc_mesoregions_2024.gpkg, SC_Municipios_2024.gpkg
OUT_DIR = None                                                   # None = <ROOT>/analysis_v10
YEARS   = list(range(2003, 2026))
PREV_2019_KG = dict(ISOP=33_280_584.0, CO=7_495_807.0, NO=571_106.94)   # Hoinaski et al. (2022), LCQAr report, Table 12
RHO_ISOP = 0.96   # isoprene canopy loss and production factor, Guenther et al. (2006), Sect. 3.3; notebook 10 writes rho = 1
MAKE_FIGS = True

# %%
# 2 · Load the annual grids, the SC / mesoregion / municipality rasters, drivers and fix logs ----------------
def year_files(root, y):
    """annual_emis_<y>.nc (notebook 10) or, failing that, the 12 monthly files; [] if neither is complete."""
    d = os.path.join(root, str(y), 'emissions'); f = os.path.join(d, f'annual_emis_{y}.nc')
    if os.path.exists(f): return [f]
    m = sorted(glob.glob(os.path.join(d, f'emis_{y}_[01][0-9].nc')))
    return m if len(m) == 12 else []
def resolve_root(root):
    """ROOT of cell 1 if it holds every year; otherwise the only folder named megan_sc_outputs_v10* (MyDrive, up to
    two levels down) that does. Anything else stops with a listing of what was found."""
    if all(year_files(root, y) for y in YEARS): return root
    cands = []
    if IN_COLAB:
        for pat in ('megan_sc_outputs_v10*', '*/megan_sc_outputs_v10*', '*/*/megan_sc_outputs_v10*'):
            cands += [p for p in glob.glob(os.path.join('/content/drive/MyDrive', pat)) if os.path.isdir(p) and p != root]
    full = [p for p in cands if all(year_files(p, y) for y in YEARS)]
    if len(full) == 1:
        print(f'ROOT {root} lacks some years; using {full[0]}'); return full[0]
    have = sorted(os.listdir(root)) if os.path.isdir(root) else 'folder not found'
    raise FileNotFoundError(f'ROOT = {root}\n  contents: {have}\n'
                            f'  years without annual_emis_<y>.nc or 12 monthly emis_<y>_MM.nc: {[y for y in YEARS if not year_files(root, y)]}\n'
                            f'  other megan_sc_outputs_v10* folders under MyDrive: {cands or "none"} (complete: {full or "none"})\n'
                            'Set ROOT in cell 1 to the folder written by notebook 10 (one sub-folder per year).')
ROOT = resolve_root(ROOT)
OUT_DIR = OUT_DIR or os.path.join(ROOT, 'analysis_v10')
_miss = [_f for _f in ('sc_state_ibge.gpkg', 'sc_mesoregions_2024.gpkg', 'SC_Municipios_2024.gpkg') if not os.path.exists(os.path.join(BND, _f))]
if _miss:
    raise FileNotFoundError(f'not found in BND = {BND}: {_miss}. Upload the three files of data/boundaries of the repository '
                            '(MEGAN_SC) to a Drive folder and set BND in cell 1 to that folder.')
os.makedirs(OUT_DIR, exist_ok=True)
print('ROOT', ROOT, '| BND', BND, '| OUT_DIR', OUT_DIR)
import numpy as np, pandas as pd, rasterio, geopandas as gpd, matplotlib as mpl, matplotlib.pyplot as plt
from netCDF4 import Dataset
from rasterio.features import rasterize
from rasterio.transform import from_origin
from scipy import stats, ndimage
import pymannkendall as pmk
warnings.filterwarnings('ignore', category=RuntimeWarning)
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

def read_year(y):
    """Annual totals per class (t/yr), south-first. annual_emis_<y>.nc (notebook 10) or the sum of the 12 monthly files;
    isoprene times RHO_ISOP."""
    files = year_files(ROOT, y)
    assert files, f'{y}: no annual_emis_{y}.nc and not 12 monthly files in {os.path.join(ROOT, str(y), "emissions")}'
    E = {c: 0.0 for c in CLASSES}
    for fn in files:
        with Dataset(fn) as ds:
            lat = np.array(ds.variables['lat'][:]); lon = np.array(ds.variables['lon'][:])
            for c in CLASSES:
                a = np.array(ds.variables[c][:], dtype='float64'); E[c] = E[c] + np.where(np.isfinite(a) & (a > -9990), a, 0.0)
    E['ISOP'] = E['ISOP'] * RHO_ISOP                    # canopy loss and production factor (cell 1)
    return E, lat, lon
cache = os.path.join(OUT_DIR, 'annual_grids_cache.npz')
SIG = json.dumps([os.path.realpath(ROOT), RHO_ISOP] + [[os.path.basename(p), os.path.getsize(p)] for y in YEARS for p in year_files(ROOT, y)])
def _cache_ok():
    if not os.path.exists(cache): return False
    Z = np.load(cache)
    return 'sig' in Z.files and str(Z['sig']) == SIG and json.loads(str(Z['years'])) == YEARS
if _cache_ok():
    Z = np.load(cache); lat, lon = Z['lat'], Z['lon']; GR = {y: {c: Z[f'{c}_{y}'].astype('float64') for c in CLASSES} for y in YEARS}
else:
    GR = {}
    for y in YEARS:
        GR[y], lat, lon = read_year(y)
    np.savez_compressed(cache, lat=lat, lon=lon, years=json.dumps(YEARS), sig=SIG, **{f'{c}_{y}': GR[y][c].astype('float32') for y in YEARS for c in CLASSES})
NY, NX = lat.size, lon.size; RES = float(lon[1] - lon[0])
TR = from_origin(lon[0] - RES / 2, lat[-1] + RES / 2, RES, RES)
def rast(geoms_ids):
    return np.flipud(rasterize(geoms_ids, out_shape=(NY, NX), transform=TR, fill=0, dtype='int32'))
R_E = 6371008.8
A = (R_E ** 2 * np.deg2rad(RES) * np.abs(np.sin(np.deg2rad(lat + RES / 2)) - np.sin(np.deg2rad(lat - RES / 2))))[:, None] * np.ones((1, NX))
def _union(g): return g.geometry.union_all() if hasattr(g.geometry, 'union_all') else g.geometry.unary_union
sc_g = gpd.read_file(os.path.join(BND, 'sc_state_ibge.gpkg')).to_crs('EPSG:4326')
SC = rast([(_union(sc_g), 1)]) == 1
meso_g = gpd.read_file(os.path.join(BND, 'sc_mesoregions_2024.gpkg')).to_crs('EPSG:4326')   # unions of the 2024 municipalities (islands included)
mun_g = gpd.read_file(os.path.join(BND, 'SC_Municipios_2024.gpkg')).to_crs('EPSG:4326')
MESO_NAMES = list(meso_g['name_meso']); FOREST = dict(zip(meso_g['name_meso'], meso_g['forest_pct']))   # MapBiomas Col. 10, 2019
AREA_MESO = dict(zip(meso_g['name_meso'], meso_g['area_ibge_km2']))   # km2, sum of the IBGE areas of the municipalities
import shapely
from scipy import sparse
_x0 = np.broadcast_to(lon[None, :] - RES / 2, (NY, NX)); _y0 = np.broadcast_to(lat[:, None] - RES / 2, (NY, NX))
CELLS = shapely.box(_x0.ravel(), _y0.ravel(), (_x0 + RES).ravel(), (_y0 + RES).ravel())   # south-first, row by row
_TREE = shapely.STRtree(CELLS)
def area_weights(geoms):
    """Sparse matrix (polygons x cells) of the fraction of each grid cell inside each polygon (area weighting). The fraction
    is taken in longitude-latitude; inside a 0.025° cell it differs from the equal-area fraction by less than 0.02 %."""
    r_, c_, v_ = [], [], []
    for k, g in enumerate(geoms):
        idx = _TREE.query(g, predicate='intersects'); f = shapely.area(shapely.intersection(CELLS[idx], g)) / RES ** 2
        keep = f > 0; r_ += [k] * int(keep.sum()); c_ += list(idx[keep]); v_ += list(f[keep])
    return sparse.csr_matrix((v_, (r_, c_)), shape=(len(geoms), NY * NX))
def wsum(W, field):
    """Area-weighted sum of a grid field over each polygon of W."""
    return np.asarray(W @ np.asarray(field, dtype='float64').ravel()).ravel()
WMESO = area_weights(list(meso_g.geometry)); WMUN = area_weights(list(mun_g.geometry))
SC_AREA = float((A * SC).sum())
_ra = wsum(WMUN, A) / 1e6 / mun_g['AREA_KM2'].values
print(f'grid {NY} x {NX}; SC cells {int(SC.sum())} ({SC_AREA / 1e6:,.0f} km2); area attributed to each municipality / IBGE area: '
      f'{_ra.min():.4f} to {_ra.max():.4f}')

def sc_mean_monthly(y, kind, band, positive=False):
    """SC mean of the annual mean of the monthly maps (notebook 03 method)."""
    files = sorted(glob.glob(os.path.join(ROOT, str(y), 'inputs', f'{kind}_{y}_[01][0-9].tif')))
    acc = np.zeros((NY, NX)); cnt = np.zeros((NY, NX))
    for f in files:
        with rasterio.open(f) as s: b = np.flipud(s.read(band).astype('float64'))
        ok = np.isfinite(b); acc[ok] += b[ok]; cnt[ok] += 1
    m = np.where(cnt > 0, acc / np.where(cnt > 0, cnt, 1), np.nan)
    if positive: m = np.where(m > 0, m, np.nan)
    return float(np.nanmean(np.where(SC, m, np.nan))), len(files)
DRV = pd.DataFrame(index=YEARS, dtype=float)
for y in YEARS:
    DRV.loc[y, 'TEMP2_C'] = sc_mean_monthly(y, 'met', 1)[0] - 273.15
    DRV.loc[y, 'PAR'] = sc_mean_monthly(y, 'met', 2)[0]
    DRV.loc[y, 'SOILM'] = sc_mean_monthly(y, 'met', 6)[0]
    DRV.loc[y, 'LAI'] = sc_mean_monthly(y, 'lai', 1, positive=True)[0]
DRV.index.name = 'year'; DRV.to_csv(os.path.join(OUT_DIR, 'drivers_annual_SC.csv'))
LOGS = {}
for y in YEARS:
    fl = sorted(glob.glob(os.path.join(ROOT, str(y), 'inputs', f'fixlog_{y}_[01][0-9].npz')))
    if fl:
        z = [np.load(f) for f in fl]
        LOGS[y] = dict(months=len(fl), filled_met_cells=int(((sum(q['filled_met_steps'] for q in z) > 0) & SC).sum()),
                       lai_gap_cells=int(sum(int(q['lai_gap_cells']) for q in z)))
LOGS = pd.DataFrame(LOGS).T
print(DRV.round(3).to_string()); print(LOGS.to_string() if len(LOGS) else 'no fix logs (not a notebook-10 run)')
if DRV.PAR.max() > 200:
    print('WARNING: SC-mean PAR above 200 W m-2: these met maps carry the ERA5-Land 00 UTC accumulation error (FIX-8 not applied)')

# %%
# 3 · State totals by year (Table 2) and by class -----------------------------------------------------------
CL = pd.DataFrame({y: {c: float((GR[y][c] * SC).sum()) for c in CLASSES} for y in YEARS}).T; CL.index.name = 'year'
G = pd.DataFrame({y: {k: float((v * SC).sum()) for k, v in groups(GR[y]).items()} for y in YEARS}).T; G.index.name = 'year'
CL.to_csv(os.path.join(OUT_DIR, 'annual_totals_SC_by_class.csv')); G.to_csv(os.path.join(OUT_DIR, 'table2_annual_totals_SC_groups.csv'))
print(G[['BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'NO', 'CO']].round(0).astype(int).to_string())

# %%
# 4 · Numbers for the text: series, trends, composition --------------------------------------------------------
NUM = []
def fmt(v, nd=1, lang='en'):
    s = f'{v:,.{nd}f}'
    return s.replace(',', '§').replace('.', ',').replace('§', '.') if lang == 'pt' else s
def put(key, desc, value, en=None, pt=None):
    en = en if en is not None else str(value)
    NUM.append(dict(key=key, description=desc, value=value, en=en, pt=pt if pt is not None else en))
bv = G['BVOC']
put('bvoc_mean_kt', 'mean BVOC 2003-2025 (kt/yr)', bv.mean() / 1e3, fmt(bv.mean() / 1e3, 0), fmt(bv.mean() / 1e3, 0, 'pt'))
cv = bv.std(ddof=1) / bv.mean() * 100
put('bvoc_cv_pct', 'interannual CV of BVOC (%, ddof=1)', cv, fmt(cv, 1), fmt(cv, 1, 'pt'))
put('bvoc_max', 'BVOC maximum: year / thousand t', bv.max(), f'{bv.idxmax()} / {fmt(bv.max() / 1e3, 0)}', f"{bv.idxmax()} / {fmt(bv.max() / 1e3, 0, 'pt')}")
put('bvoc_min', 'BVOC minimum: year / thousand t', bv.min(), f'{bv.idxmin()} / {fmt(bv.min() / 1e3, 0)}', f"{bv.idxmin()} / {fmt(bv.min() / 1e3, 0, 'pt')}")
for k in ('BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'NO', 'CO'):
    put(f'y2025_{k}_t', f'2025 {k} (t)', round(G.loc[2025, k]), fmt(G.loc[2025, k], 0), fmt(G.loc[2025, k], 0, 'pt'))
put('y2025_text', '2025 BVOC / ISOP / MT (thousand t)', G.loc[2025, 'BVOC'],
    f"{fmt(G.loc[2025, 'BVOC'] / 1e3)} / {fmt(G.loc[2025, 'ISOP'] / 1e3)} / {fmt(G.loc[2025, 'MT'] / 1e3)}",
    f"{fmt(G.loc[2025, 'BVOC'] / 1e3, 1, 'pt')} / {fmt(G.loc[2025, 'ISOP'] / 1e3, 1, 'pt')} / {fmt(G.loc[2025, 'MT'] / 1e3, 1, 'pt')}")
def trend(s):
    y = np.asarray(s, float); x = np.arange(len(y)); o = pmk.original_test(y); hr = pmk.hamed_rao_modification_test(y)
    sl, ic, lo, hi = stats.theilslopes(y, x, 0.95)
    res = y - np.polyval(np.polyfit(x, y, 1), x)
    return dict(mean=y.mean(), sen=sl, pct_yr=sl / y.mean() * 100, pct_lo95=lo / y.mean() * 100, pct_hi95=hi / y.mean() * 100,
                MK_Z=o.z, MK_p=o.p, HR_p=hr.p, ar1_detrended=np.corrcoef(res[:-1], res[1:])[0, 1])
TRD = pd.DataFrame({k: trend(G[k]) for k in ('BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'MBO', 'NO', 'CO')}).T
TRD.index.name = 'series'; TRD.to_csv(os.path.join(OUT_DIR, 'trends.csv'))
for k, r in TRD.iterrows():
    put(f'trend_{k}', f'trend {k}: Sen %/yr (95% CI); Hamed-Rao p; original MK p', r.pct_yr,
        f'{r.pct_yr:+.2f}% yr-1 ({r.pct_lo95:.2f} to {r.pct_hi95:.2f}; p = {r.HR_p:.3f}; MK p = {r.MK_p:.3f})',
        f"{fmt(r.pct_yr, 2, 'pt')}% ano-1 ({fmt(r.pct_lo95, 2, 'pt')} a {fmt(r.pct_hi95, 2, 'pt')}; p = {fmt(r.HR_p, 3, 'pt')}; MK p = {fmt(r.MK_p, 3, 'pt')})")
for k in ('ISOP', 'MT', 'OVOC', 'LVOC', 'STRESS_OTHER', 'SQT', 'MBO'):
    m = G[k].mean() / 1e3; sh = 100 * G[k].mean() / bv.mean(); nd = 0 if m >= 1 else 1
    put(f'mean_{k}', f'mean 2003-2025 {k} kt/yr (share of BVOC)', m, f'{fmt(m, nd)} kt yr-1 ({fmt(sh)}%)', f"{fmt(m, nd, 'pt')} kt ano-1 ({fmt(sh, 1, 'pt')}%)")
ratio = G['ISOP'].mean() / G['MT'].mean()
put('isop_mt_ratio', 'isoprene:monoterpene ratio of the 2003-2025 means', ratio, f'{ratio:.2f}:1', f"{fmt(ratio, 2, 'pt')}:1")
print(pd.DataFrame(NUM)[['key', 'en']].to_string(index=False))

# %%
# 5 · Emission-driver correlations (raw and linearly detrended annual series) ---------------------------------
def detr(v):
    x = np.arange(len(v)); v = np.asarray(v, float); return v - np.polyval(np.polyfit(x, v, 1), x)
rows = []
SER = {**{c: CL[c] for c in CLASSES}, **{k: G[k] for k in ('BVOC', 'MT', 'SQT', 'OVOC')}}
for k, s in SER.items():
    for d in ('TEMP2_C', 'PAR', 'LAI', 'SOILM'):
        r, p = stats.pearsonr(s.values, DRV[d].values); rd, pd_ = stats.pearsonr(detr(s.values), detr(DRV[d].values))
        rows.append(dict(series=k, driver=d, r=r, p=p, r_detrended=rd, p_detrended=pd_))
COR = pd.DataFrame(rows); COR.to_csv(os.path.join(OUT_DIR, 'correlations_emissions_drivers.csv'), index=False)
cb = COR[COR.series.isin(BVOC)]
for d in ('TEMP2_C', 'LAI', 'PAR', 'SOILM'):
    v = cb[cb.driver == d].r_detrended
    put(f'rdet_bvoc_classes_{d}', f'detrended r, 17 BVOC classes vs {d}: min to max', (round(v.min(), 3), round(v.max(), 3)),
        f'{v.min():.2f} to {v.max():.2f}', f"{fmt(v.min(), 2, 'pt')} a {fmt(v.max(), 2, 'pt')}")
no = COR[(COR.series == 'NO') & (COR.driver == 'TEMP2_C')].iloc[0]
put('r_NO_T', 'soil NO vs temperature: raw r ; detrended r', (round(no.r, 3), round(no.r_detrended, 3)),
    f'{no.r:.2f} ; {no.r_detrended:.2f}', f"{fmt(no.r, 2, 'pt')} ; {fmt(no.r_detrended, 2, 'pt')}")
print(COR.pivot(index='series', columns='driver', values='r_detrended').round(2).to_string())

# %%
# 6 · Mesoregions: SOA precursors (MT + SQT), 2003-2025 mean ------------------------------------------------
def meso_table(W):
    acc = {k: np.zeros(len(MESO_NAMES)) for k in ('MT', 'SQT', 'ISOP', 'NO', 'BVOC')}
    for y in YEARS:
        gg = groups(GR[y])
        for k in acc: acc[k] += wsum(W, gg[k]) / len(YEARS)
    T = pd.DataFrame(acc, index=MESO_NAMES); T['SOA'] = T.MT + T.SQT; T['pct_SC'] = 100 * T.SOA / T.SOA.sum()
    T['area_ibge_ha'] = [AREA_MESO[n] * 100 for n in T.index]; T['area_grid_ha'] = wsum(W, A) / 1e4
    T['intensity_kg_ha'] = T.SOA * 1e3 / T.area_ibge_ha; T['SQT_share_pct'] = 100 * T.SQT / T.SOA
    T['forest_pct'] = [FOREST[n] for n in T.index]; T.index.name = 'mesoregion'
    return T.sort_values('SOA', ascending=False)
MESO = meso_table(WMESO)
MESO.to_csv(os.path.join(OUT_DIR, 'table3_mesoregions_soa.csv'))
state_soa = np.mean([groups(GR[y])['SOA'][SC].sum() for y in YEARS])
put('soa_total_kt', 'SOA precursors (MT+SQT) 2003-2025 mean, state (kt/yr)', state_soa / 1e3, fmt(state_soa / 1e3, 0), fmt(state_soa / 1e3, 0, 'pt'))
put('meso_sum', 'sum of the six mesoregions (area weighting) / state total (cell-centre rule), SOA 2003-2025 mean (%)',
    100 * MESO.SOA.sum() / state_soa, f'{100 * MESO.SOA.sum() / state_soa:.2f}%', f"{fmt(100 * MESO.SOA.sum() / state_soa, 2, 'pt')}%")
for n, r in MESO.iterrows():
    put(f'meso_{n}', f'{n}: SOA t/yr; % of SC; kg/ha/yr; SQT %', r.SOA,
        f'{r.SOA:,.0f}; {r.pct_SC:.1f}%; {r.intensity_kg_ha:.1f}; {r.SQT_share_pct:.1f}%',
        f"{fmt(r.SOA, 0, 'pt')}; {fmt(r.pct_SC, 1, 'pt')}%; {fmt(r.intensity_kg_ha, 1, 'pt')}; {fmt(r.SQT_share_pct, 1, 'pt')}%")
x, yv = MESO.forest_pct.values.astype(float), MESO.intensity_kg_ha.values
r, p = stats.pearsonr(x, yv); rho, prho = stats.spearmanr(x, yv); zz = np.arctanh(r); ci = np.tanh([zz - 1.96 / np.sqrt(len(x) - 3), zz + 1.96 / np.sqrt(len(x) - 3)])
put('meso_forest_r', 'SOA intensity vs native forest: Pearson r (Fisher 95% CI); Spearman rho; n', r,
    f'r = {r:.2f} (95% CI {ci[0]:.2f} to {ci[1]:.2f}; p = {p:.3f}); rho = {rho:.2f} (p = {prho:.3f}); n = {len(x)}',
    f"r = {fmt(r, 2, 'pt')} (IC 95% {fmt(ci[0], 2, 'pt')} a {fmt(ci[1], 2, 'pt')}; p = {fmt(p, 3, 'pt')}); rho = {fmt(rho, 2, 'pt')} (p = {fmt(prho, 3, 'pt')}); n = {len(x)}")
put('meso_sqt_share_range', 'SQT share of SOA precursors, min to max over mesoregions (%)', (round(MESO.SQT_share_pct.min(), 2), round(MESO.SQT_share_pct.max(), 2)),
    f'{MESO.SQT_share_pct.min():.1f}% to {MESO.SQT_share_pct.max():.1f}%', f"{fmt(MESO.SQT_share_pct.min(), 1, 'pt')}% a {fmt(MESO.SQT_share_pct.max(), 1, 'pt')}%")
put('meso_top_isop', 'mesoregion with the largest isoprene total', MESO.ISOP.idxmax())
put('meso_top_no', 'mesoregions with the largest soil NO (t/yr)', '', '; '.join(f'{n} {v:,.0f}' for n, v in MESO.NO.sort_values(ascending=False).head(3).items()))
print(MESO.round(1).to_string())

# %%
# 7 · Municipalities, 2025 (Table 4, Figure 7) and 2019 (Figure 8) --------------------------------------------
def fisher_jenks(values, k=5):
    """Jenks natural breaks by Fisher's exact dynamic programme (same limits as mapclassify.FisherJenks): k + 1 limits, the
    first value of each class after the first, with the minimum and the maximum at the ends (analysis/v10/maps/wf_style.py)."""
    x = np.sort(np.asarray(values, float)); x = x[np.isfinite(x)]; n = x.size
    if n <= k: return np.unique(x)
    c1 = np.concatenate([[0.0], np.cumsum(x)]); c2 = np.concatenate([[0.0], np.cumsum(x * x)])
    D = np.full((k + 1, n + 1), np.inf); B = np.zeros((k + 1, n + 1), int); D[0, 0] = 0.0
    for m in range(1, k + 1):
        for j in range(m, n + 1):
            i = np.arange(m - 1, j); s1 = c1[j] - c1[i]; cost = D[m - 1, i] + (c2[j] - c2[i]) - s1 * s1 / (j - i)
            a = int(np.argmin(cost)); D[m, j], B[m, j] = cost[a], i[a]
    ends, j = [], n
    for m in range(k, 0, -1):
        ends.append(j); j = B[m, j]
    ends = sorted(ends)
    return np.array([x[0]] + [x[e] for e in ends[:-1]] + [x[-1]])
NCELL_MUN = np.asarray((WMUN > 0).sum(1)).ravel(); CELL_KM2 = float(A[SC].mean() / 1e6)
def muni_table(y, W=WMUN):
    """Municipal totals of year y by area weighting; n_cells = grid cells that intersect the municipality."""
    gg = groups(GR[y])
    T = pd.DataFrame({'mun': list(mun_g['NM_MUN']), 'area_ibge_km2': mun_g['AREA_KM2'].astype(float).values, 'n_cells': NCELL_MUN,
                      'area_grid_km2': wsum(W, A) / 1e6, **{k: wsum(W, gg[k]) for k in ('BVOC', 'ISOP', 'MT', 'SQT', 'NO', 'CO')}}).set_index('mun')
    T['intensity_kg_ha'] = T.BVOC * 1e3 / (T.area_ibge_km2 * 100)
    return T
MU = muni_table(2025); MU.to_csv(os.path.join(OUT_DIR, 'table4_municipalities_2025.csv'))
MU19 = muni_table(2019); MU19.to_csv(os.path.join(OUT_DIR, 'municipalities_2019.csv'))
top = MU.sort_values('BVOC', ascending=False).head(12)
put('mun_top12_bvoc', 'top-12 municipalities 2025: BVOC t (kg/ha/yr)', '', '; '.join(f'{n} {r_.BVOC:,.0f} ({r_.intensity_kg_ha:.0f})' for n, r_ in top.iterrows()))
put('mun_top5_isop', 'top-5 isoprene municipalities 2025 (t)', '', '; '.join(f'{n} {v:,.0f}' for n, v in MU.ISOP.sort_values(ascending=False).head(5).items()))
put('mun_top5_no', 'top-5 soil-NO municipalities 2025 (t)', '', '; '.join(f'{n} {v:,.0f}' for n, v in MU.NO.sort_values(ascending=False).head(5).items()))
it = MU.intensity_kg_ha
put('mun_intensity', 'municipal BVOC intensity 2025 (kg/ha/yr, IBGE area): min to max; median', it.median(),
    f'{it.min():.0f} to {it.max():.0f} (median {it.median():.0f}; max {it.idxmax()})', f"{it.min():.0f} a {it.max():.0f} (mediana {it.median():.0f}; máximo {it.idxmax()})")
JB = fisher_jenks(it.values, 5); topc = it[it >= JB[-2]].sort_values(ascending=False)
put('mun_top_jenks', 'top Jenks class of intensity: lower limit (kg/ha/yr), number of municipalities, the first ones', JB[-2],
    f'>= {JB[-2]:.1f} kg ha-1 yr-1 ({len(topc)} municipalities; limits ' + ', '.join(f'{v:.1f}' for v in JB) + '): ' + ', '.join(topc.index[:8]),
    f"≥ {fmt(JB[-2], 1, 'pt')} kg ha-1 ano-1 ({len(topc)} municípios; limites " + '; '.join(fmt(v, 1, 'pt') for v in JB) + '): ' + ', '.join(topc.index[:8]))
eq = MU.area_ibge_km2 / CELL_KM2
put('mun_cells', 'municipalities: grid cells that intersect each one (min, median, max); IBGE area in cell areas (min, median, max)', float(eq.median()),
    f'{MU.n_cells.min()}, {MU.n_cells.median():.0f}, {MU.n_cells.max()}; {eq.min():.1f}, {eq.median():.0f}, {eq.max():.0f} (cell {CELL_KM2:.2f} km2)',
    f"{MU.n_cells.min()}, {MU.n_cells.median():.0f}, {MU.n_cells.max()}; {fmt(eq.min(), 1, 'pt')}, {eq.median():.0f}, {eq.max():.0f} (célula de {fmt(CELL_KM2, 2, 'pt')} km²)")
ra = MU.area_grid_km2 / MU.area_ibge_km2; sm = 100 * MU.BVOC.sum() / G.loc[2025, 'BVOC']
put('mun_area_check', 'area attributed to each municipality / IBGE area (min to max); sum of the 295 municipalities / state total, BVOC 2025 (%)', sm,
    f'{ra.min():.4f} to {ra.max():.4f}; {sm:.2f}%', f"{fmt(ra.min(), 3, 'pt')} a {fmt(ra.max(), 3, 'pt')}; {fmt(sm, 2, 'pt')}%")

# %%
# 8 · Isoprene areal flux, 2019 comparison with the previous inventory, fix-log statistics --------------------
CF = 5 * 12.011 / 68.117
mean_iso = np.mean([GR[y]['ISOP'] for y in YEARS], axis=0)
state_flux = (mean_iso * SC).sum() * 1e6 * CF / SC_AREA
cell = (mean_iso * 1e6 * CF / A)[SC & (mean_iso > 0)]
put('isop_flux', 'isoprene areal flux g C m-2 yr-1: state mean (cell median; p99; max)', state_flux,
    f'{state_flux:.1f} (cell median {np.median(cell):.1f}; p99 {np.percentile(cell, 99):.1f}; max {cell.max():.1f})',
    f"{fmt(state_flux, 1, 'pt')} (mediana por célula {fmt(np.median(cell), 1, 'pt')}; p99 {fmt(np.percentile(cell, 99), 1, 'pt')}; máximo {fmt(cell.max(), 1, 'pt')})")
rat = {}
for k in ('ISOP', 'CO', 'NO'):
    now = G.loc[2019, k] / 1e3; prev = PREV_2019_KG[k] / 1e6; rat[k] = now / prev
    nd = 3 if prev < 1 else (2 if prev < 10 else 1)
    put(f'cmp2019_{k}', f'2019 {k}: this inventory Gg vs previous Gg; ratio', rat[k], f'{now:.1f} Gg vs {prev:.{nd}f} Gg; ratio {rat[k]:.1f}',
        f"{fmt(now, 1, 'pt')} Gg contra {fmt(prev, nd, 'pt')} Gg; razão {fmt(rat[k], 1, 'pt')}")
put('cmp2019_range', 'range of the three 2019 ratios', (round(min(rat.values()), 2), round(max(rat.values()), 2)),
    f'{min(rat.values()):.0f} to {max(rat.values()):.0f} times', f'{min(rat.values()):.0f} a {max(rat.values()):.0f} vezes')
iso = G['ISOP'] / 1e3; pv = PREV_2019_KG['ISOP'] / 1e6
put('isop_mean_ratio', 'isoprene: 2003-2025 mean Gg; mean/previous; 2025/previous; peak year', iso.mean() / pv,
    f'{iso.mean():.0f} Gg; {iso.mean() / pv:.1f}; {iso.loc[2025] / pv:.1f}; peak {iso.idxmax()}',
    f"{iso.mean():.0f} Gg; {fmt(iso.mean() / pv, 1, 'pt')}; {fmt(iso.loc[2025] / pv, 1, 'pt')}; pico {iso.idxmax()}")
if len(LOGS):
    put('lai_gaps', 'LAI gap cell-months filled (FIX-4) per year: min to max', '', f"{int(LOGS.lai_gap_cells.min())} to {int(LOGS.lai_gap_cells.max())}")
    put('met_filled', 'SC cells with filled ERA5-Land meteorology (FIX-1) per year: min to max', '', f"{int(LOGS.filled_met_cells.min())} to {int(LOGS.filled_met_cells.max())}")
put('grid_facts', 'grid; SC cells (% of box); SC area on grid (km2)', int(SC.sum()), f'{NY} x {NX}; {int(SC.sum())} ({100 * SC.mean():.1f}%); {SC_AREA / 1e6:,.0f}')

# %%
# 8b · Seasonal cycle: monthly state totals by class (stats_<year>_<mm>.csv of notebook 10; isoprene times RHO_ISOP) ----
MROWS = []
for y in YEARS:
    for m_ in range(1, 13):
        f = os.path.join(ROOT, str(y), 'stats', f'stats_{y}_{m_:02d}.csv')
        d = pd.read_csv(f); t = dict(zip(d.species, d.total_tonnes)); t['ISOP'] = t['ISOP'] * RHO_ISOP
        MROWS.append(dict(year=y, month=m_, **groups(t)))
MON = pd.DataFrame(MROWS); MON['nonISOP'] = MON.BVOC - MON.ISOP
_ann = MON.groupby('year')[['BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'NO', 'CO']].sum()
_dev = float(((_ann - G[_ann.columns]) / G[_ann.columns]).abs().max().max())
assert _dev < 1e-6, f'monthly stats do not add up to the annual grids (max relative difference {_dev:.2e})'
MON.to_csv(os.path.join(OUT_DIR, 'seasonal_cycle_monthly.csv'), index=False)
SEAS_K = ['BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'STRESS_OTHER', 'MBO', 'nonISOP', 'NO', 'CO']
CLIM = MON.groupby('month')[SEAS_K].mean(); CMIN = MON.groupby('month')[SEAS_K].min(); CMAX = MON.groupby('month')[SEAS_K].max()
pd.concat({'mean': CLIM, 'min': CMIN, 'max': CMAX}, axis=1).to_csv(os.path.join(OUT_DIR, 'seasonal_climatology.csv'))
MES_EN = ['', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
MES_PT = ['', 'jan', 'fev', 'mar', 'abr', 'mai', 'jun', 'jul', 'ago', 'set', 'out', 'nov', 'dez']
for k in SEAS_K:
    c_ = CLIM[k]; oc = c_.loc[[10, 11, 12, 1, 2, 3]].sum() / c_.sum() * 100; djf = c_.loc[[12, 1, 2]].sum() / c_.sum() * 100
    put(f'season_{k}', f'{k} monthly climatology 2003-2025: max month (t); min month (t); max/min; Oct-Mar and DJF share of the year (%)',
        c_.max() / c_.min(),
        f'{MES_EN[c_.idxmax()]} {c_.max():,.0f} t; {MES_EN[c_.idxmin()]} {c_.min():,.0f} t; ratio {c_.max() / c_.min():.2f}; Oct-Mar {oc:.1f}%; DJF {djf:.1f}%',
        f"{MES_PT[c_.idxmax()]} {fmt(c_.max(), 0, 'pt')} t; {MES_PT[c_.idxmin()]} {fmt(c_.min(), 0, 'pt')} t; razão {fmt(c_.max() / c_.min(), 2, 'pt')}; "
        f"out-mar {fmt(oc, 1, 'pt')}%; dez-fev {fmt(djf, 1, 'pt')}%")
mxm = MON.loc[MON.groupby('year').BVOC.idxmax(), 'month'].value_counts().sort_index()
mnm = MON.loc[MON.groupby('year').BVOC.idxmin(), 'month'].value_counts().sort_index()
put('season_bvoc_max_min_months', 'BVOC: number of years with the monthly maximum / minimum in each month', '',
    'max: ' + ', '.join(f'{MES_EN[m_]} {n_}' for m_, n_ in mxm.items()) + ' | min: ' + ', '.join(f'{MES_EN[m_]} {n_}' for m_, n_ in mnm.items()),
    'máx.: ' + ', '.join(f'{MES_PT[m_]} {n_}' for m_, n_ in mxm.items()) + ' | mín.: ' + ', '.join(f'{MES_PT[m_]} {n_}' for m_, n_ in mnm.items()))
shr = 100 * CLIM.ISOP / CLIM.BVOC
put('season_isop_share', 'isoprene share of monthly BVOC (climatology): Jun, Jul, Dec, Jan (%)', '',
    '; '.join(f'{MES_EN[m_]} {shr[m_]:.1f}%' for m_ in (6, 7, 12, 1)), '; '.join(f"{MES_PT[m_]} {fmt(shr[m_], 1, 'pt')}%" for m_ in (6, 7, 12, 1)))
print((CLIM / 1e3).round(2).to_string())

# %%
# 9 · The table of manuscript numbers ------------------------------------------------------------------------
N = pd.DataFrame(NUM); N.to_csv(os.path.join(OUT_DIR, 'manuscript_numbers.csv'), index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 170)
print(N[['key', 'en']].to_string(index=False))

# %%
# 10 · Figures 1-8 (PT and EN), same layout as the revised manuscript ----------------------------------------
from matplotlib.ticker import MultipleLocator
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPoly, Patch
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": "#4d4d4d", "axes.linewidth": 0.9,
                     "text.color": "#1a1a1a", "axes.labelcolor": "#1a1a1a", "xtick.color": "#4d4d4d", "ytick.color": "#4d4d4d"})
PALC = dict(tot="#3d4a57", iso="#0072B2", mt="#D55E00", sqt="#E69F00", ov="#009E73", lv="#CC79A7", oth="#999999", mbo="#56B4E9", no="#C0392B")
LBL = {'pt': dict(f1t="Emissões biogênicas anuais, Santa Catarina, 2003–2025", f1s="painéis b e c em eixo comum",
                  tot="Total de COVB", iso="Isopreno", mt="Monoterpenos", sqt="Sesquiterpenos", no="NO do solo", ov="COV oxigenados",
                  lv="Hidrocarbonetos leves (LVOC)", oth="Outros/estresse", mbo="MBO", unit="(kt ano⁻¹)",
                  f2t="Composição média das emissões de COVB (2003–2025), total ≈ {:.0f} kt ano⁻¹", f2x="Emissão média 2003–2025 (kt ano⁻¹)",
                  f3t="Correlação entre emissões anuais e forçantes\n(séries 2003–2025 sem tendência linear)", drv=["Temperatura (2 m)", "PAR", "IAF", "Umidade do solo"], cb="r de Pearson",
                  f4t="Precursores de AOS por mesorregião (média 2003–2025)", f4x="Emissão de precursores de AOS (kt ano⁻¹)",
                  f5t="Intensidade de precursores de AOS por mesorregião (kg ha⁻¹ ano⁻¹), média 2003–2025",
                  f6t="Intensidade de AOS vs. floresta nativa (n = 6)", f6x="Fração de floresta nativa, MapBiomas 2019 (%)", f6y="Intensidade de AOS (kg ha⁻¹ ano⁻¹)",
                  f7a="(a) Emissão total de COVB (t ano⁻¹)", f7b="(b) Intensidade de COVB (kg ha⁻¹ ano⁻¹)", f7t="Emissões de COVB por município, 2025",
                  f8t="Isopreno, CO e NO do solo por município (t ano⁻¹): 2019 e 2025", co="CO", q="quebras naturais de Jenks", yr='ano⁻¹',
                  f2b="b) Ciclo sazonal (média 2003–2025)", f2c="c) Ciclo normalizado pela média mensal", f2y="Emissão mensal (kt mês⁻¹)",
                  f2yn="Emissão / média mensal", band="COVB total: menor e maior valor", ref="média mensal (= 1)", mes=MES_PT[1:]),
       'en': dict(f1t="Annual biogenic emissions, Santa Catarina, 2003–2025", f1s="panels b and c share a common axis",
                  tot="Total BVOC", iso="Isoprene", mt="Monoterpenes", sqt="Sesquiterpenes", no="Soil NO", ov="Oxygenated VOC",
                  lv="Light hydrocarbons (LVOC)", oth="Other/stress", mbo="MBO", unit="(kt yr⁻¹)",
                  f2t="Mean BVOC composition (2003–2025), total ≈ {:.0f} kt yr⁻¹", f2x="Mean emission 2003–2025 (kt yr⁻¹)",
                  f3t="Correlation between annual emissions and drivers\n(linearly detrended 2003–2025 series)", drv=["Temperature (2 m)", "PAR", "LAI", "Soil moisture"], cb="Pearson r",
                  f4t="SOA precursors by mesoregion (2003–2025 mean)", f4x="SOA-precursor emission (kt yr⁻¹)",
                  f5t="SOA-precursor intensity by mesoregion (kg ha⁻¹ yr⁻¹), 2003–2025 mean",
                  f6t="SOA intensity vs. native forest (n = 6)", f6x="Native-forest fraction, MapBiomas 2019 (%)", f6y="SOA intensity (kg ha⁻¹ yr⁻¹)",
                  f7a="(a) Total BVOC emission (t yr⁻¹)", f7b="(b) BVOC intensity (kg ha⁻¹ yr⁻¹)", f7t="BVOC emissions by municipality, 2025",
                  f8t="Isoprene, CO and soil NO by municipality (t yr⁻¹): 2019 and 2025", co="CO", q="Jenks natural breaks", yr='yr⁻¹',
                  f2b="b) Seasonal cycle (2003–2025 mean)", f2c="c) Cycle normalised by the monthly mean", f2y="Monthly emission (kt month⁻¹)",
                  f2yn="Emission / monthly mean", band="Total BVOC: lowest and highest value", ref="monthly mean (= 1)", mes=MES_EN[1:])}
DLBL = {'pt': dict(t="Correlação com e sem a tendência linear: temperatura do ar e NO do solo, 2003–2025",
                   a="a) Temperatura média do ar a 2 m", b="b) NO do solo (total estadual)", ty="Temperatura (°C)", ny="NO do solo (kt ano⁻¹)",
                   fit="reta de mínimos quadrados", res="desvio em relação à reta", pts="valor anual (cor: ano)", c="c) Séries originais",
                   d="d) Desvios em relação às retas", tx="Temperatura (°C)", dtx="Desvio da temperatura (°C)",
                   dny="Desvio do NO do solo (kt ano⁻¹)", yr="Ano"),
        'en': dict(t="Correlation with and without the linear trend: air temperature and soil NO, 2003–2025",
                   a="a) Mean 2 m air temperature", b="b) Soil NO (state total)", ty="Temperature (°C)", ny="Soil NO (kt yr⁻¹)",
                   fit="least-squares line", res="departure from the line", pts="annual value (colour: year)", c="c) Original series",
                   d="d) Departures from the lines", tx="Temperature (°C)", dtx="Temperature departure (°C)",
                   dny="Soil NO departure (kt yr⁻¹)", yr="Year")}
ASPECT = 1 / np.cos(np.deg2rad(-27.3))
def clean(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
def nice(v, step): return step * np.ceil(v * 1.08 / step)
def dec(s, lang): return s.replace('.', ',') if lang == 'pt' else s
def patches(geom): return [MplPoly(np.asarray(p.exterior.coords), closed=True) for p in getattr(geom, 'geoms', [geom])]
def choropleth(ax, names, geoms, values, cmap, bins=None, vmin=None, vmax=None, lw=0.15):
    P, V = [], []
    for k, g in zip(names, geoms):
        for p in patches(g): P.append(p); V.append(values.get(k, np.nan))
    V = np.asarray(V, float)
    if bins is not None:
        idx = np.clip(np.digitize(V, bins[1:-1]), 0, len(bins) - 2); cm = plt.get_cmap(cmap, len(bins) - 1)
        pc = PatchCollection(P, facecolor=[cm(i) for i in idx], edgecolor='#606060', linewidth=lw)
    else:
        pc = PatchCollection(P, cmap=cmap, edgecolor='#606060', linewidth=lw); pc.set_array(V); pc.set_clim(vmin, vmax)
    ax.add_collection(pc); ax.set_xlim(-54.0, -48.2); ax.set_ylim(-29.5, -25.8); ax.set_aspect(ASPECT); ax.axis('off'); return pc
def legend_bins(ax, bins, cmap, title, lang):
    cm = plt.get_cmap(cmap, len(bins) - 1)
    lab = lambda a, b: (f"{a:,.0f}–{b:,.0f}".replace(',', '.') if lang == 'pt' else f"{a:,.0f}–{b:,.0f}")
    h = [Patch(facecolor=cm(i), edgecolor='#606060', label=lab(bins[i], bins[i + 1])) for i in range(len(bins) - 1)]
    ax.legend(handles=h, title=title, loc='lower left', fontsize=8, title_fontsize=8.5, frameon=True)
def fig_detrend(lang, od):
    """Example of the detrended correlation (article Figure 2): soil NO and temperature, whose linear trends have opposite signs."""
    yr = np.array(YEARS)
    from matplotlib.lines import Line2D
    D0 = DLBL[lang]; x0 = np.arange(len(yr)); tt = DRV['TEMP2_C'].values.astype(float); nn = (G.NO / 1e3).values.astype(float)
    fit_t = np.polyval(np.polyfit(x0, tt, 1), x0); fit_n = np.polyval(np.polyfit(x0, nn, 1), x0); rt, rn = tt - fit_t, nn - fit_n
    r_raw, r_det = stats.pearsonr(tt, nn)[0], stats.pearsonr(rt, rn)[0]
    c0 = COR[(COR.series == 'NO') & (COR.driver == 'TEMP2_C')].iloc[0]
    assert abs(r_raw - c0.r) < 1e-9 and abs(r_det - c0.r_detrended) < 1e-9          # the same r as correlations_emissions_drivers.csv
    cm0 = plt.get_cmap('viridis'); nr0 = mpl.colors.Normalize(yr.min(), yr.max()); col0 = cm0(nr0(yr))
    fig, axs = plt.subplots(2, 2, figsize=(9.6, 7.8)); fig.subplots_adjust(left=0.085, right=0.875, top=0.895, bottom=0.08, hspace=0.42, wspace=0.32)
    for ax, s, f, lab, ttl in ((axs[0, 0], tt, fit_t, D0['ty'], D0['a']), (axs[0, 1], nn, fit_n, D0['ny'], D0['b'])):
        ax.vlines(yr, f, s, color='#9a9a9a', lw=0.9, zorder=2); ax.plot(yr, f, color='#3d4a57', lw=1.4, ls='--', zorder=3)
        ax.scatter(yr, s, c=col0, s=34, edgecolor='white', linewidth=0.6, zorder=4)
        ax.set_ylabel(lab); ax.set_xticks(np.arange(2003, 2026, 4)); clean(ax); ax.set_title(ttl, fontsize=11, fontweight='bold', loc='left')
        ax.yaxis.grid(True, color="#e8e8e6", lw=0.8, zorder=0); ax.set_axisbelow(True)
    hh = [Line2D([], [], color='#3d4a57', lw=1.4, ls='--'), Line2D([], [], color='#9a9a9a', lw=0.9),
          Line2D([], [], marker='o', ls='', color=cm0(0.5), ms=5.5, markeredgecolor='white')]
    axs[0, 0].legend(hh, [D0['fit'], D0['res'], D0['pts']], fontsize=8, loc='upper left', frameon=False)
    for ax, xv, yv_, xl, yl, ttl, r in ((axs[1, 0], tt, nn, D0['tx'], D0['ny'], D0['c'], r_raw), (axs[1, 1], rt, rn, D0['dtx'], D0['dny'], D0['d'], r_det)):
        b1, b0 = np.polyfit(xv, yv_, 1); xx = np.linspace(xv.min(), xv.max(), 50)
        ax.plot(xx, b1 * xx + b0, color='#C0392B', lw=1.2, zorder=3); ax.scatter(xv, yv_, c=col0, s=34, edgecolor='white', linewidth=0.6, zorder=4)
        ax.set_xlabel(xl); ax.set_ylabel(yl); clean(ax); ax.set_title(ttl, fontsize=11, fontweight='bold', loc='left')
        ax.text(0.03, 0.97, dec(f"r = {r:.2f}", lang), transform=ax.transAxes, ha='left', va='top', fontsize=11.5, fontweight='bold')
    axs[1, 1].axhline(0, color='0.78', lw=0.8, zorder=1); axs[1, 1].axvline(0, color='0.78', lw=0.8, zorder=1)
    cax = fig.add_axes([0.905, 0.12, 0.016, 0.72]); cb0 = fig.colorbar(mpl.cm.ScalarMappable(norm=nr0, cmap=cm0), cax=cax)
    cb0.set_ticks(np.arange(2004, 2026, 4)); cb0.ax.set_yticklabels([str(v) for v in np.arange(2004, 2026, 4)]); cb0.set_label(D0['yr'])
    if lang == 'pt':                                                 # decimal comma on the axes of the Portuguese figure
        class _Comma(mpl.ticker.ScalarFormatter):
            def __call__(self, v, pos=None): return super().__call__(v, pos).replace('.', ',')
        for ax in axs.flat: ax.xaxis.set_major_formatter(_Comma()); ax.yaxis.set_major_formatter(_Comma())
    fig.suptitle(D0['t'], x=0.085, y=0.975, ha='left', fontsize=12.5, fontweight='bold')
    fig.savefig(f"{od}/fig_detrend.png", dpi=300, facecolor="white"); plt.close(fig)
mnames, mgeoms = list(meso_g['name_meso']), list(meso_g.geometry); unames, ugeoms = list(mun_g['NM_MUN']), list(mun_g.geometry)
for lang, L in (LBL.items() if MAKE_FIGS else []):
    od = os.path.join(OUT_DIR, 'figures', lang); os.makedirs(od, exist_ok=True); yr = np.array(YEARS)
    # Fig 1 - annual series
    common_ax = nice(max(G.ISOP.max(), G.MT.max()) / 1e3, 50)     # panels b and c share an axis; d has its own (authors' comment)
    panels = [("a", L['tot'], G.BVOC / 1e3, PALC['tot'], nice(G.BVOC.max() / 1e3, 100), 200), ("b", L['iso'], G.ISOP / 1e3, PALC['iso'], common_ax, 50),
              ("c", L['mt'], G.MT / 1e3, PALC['mt'], common_ax, 50), ("d", L['sqt'], G.SQT / 1e3, PALC['sqt'], nice(G.SQT.max() / 1e3, 5), 5),
              ("e", L['no'], G.NO / 1e3, PALC['no'], nice(G.NO.max() / 1e3, 1), 3)]
    fig, axes = plt.subplots(5, 1, figsize=(8.0, 12.4), sharex=True); fig.subplots_adjust(left=0.14, right=0.965, top=0.895, bottom=0.07, hspace=0.42)
    for ax, (lab, name, s, col, ymax, ymaj) in zip(axes, panels):
        ax.bar(yr, s.values, width=0.72, color=col, edgecolor="white", linewidth=0.3, zorder=3)
        ax.set_ylim(0, ymax); ax.yaxis.set_major_locator(MultipleLocator(ymaj)); ax.set_ylabel(f"{name}\n{L['unit']}", fontsize=10.5)
        clean(ax); ax.yaxis.grid(True, color="#e8e8e6", lw=0.8, zorder=0); ax.set_axisbelow(True)
        ax.text(1.0, 1.05, f"({lab})", transform=ax.transAxes, ha="right", va="bottom", fontsize=12.5, fontweight="bold")
    axes[-1].set_xlim(yr.min() - 0.7, yr.max() + 0.7); axes[-1].set_xticks(np.arange(YEARS[0], YEARS[-1] + 1, 2))
    fig.suptitle(L['f1t'], x=0.14, y=0.975, ha="left", fontsize=12.5, fontweight="bold"); fig.text(0.14, 0.945, L['f1s'], ha="left", fontsize=9, color="#6b6b6b")
    fig.savefig(f"{od}/fig1.png", dpi=300, facecolor="white"); plt.close(fig)
    # Fig 2 - composition (a) and seasonal cycle (b, c)
    fam = [(L['iso'], G.ISOP.mean(), PALC['iso']), (L['mt'], G.MT.mean(), PALC['mt']), (L['ov'], G.OVOC.mean(), PALC['ov']), (L['lv'], G.LVOC.mean(), PALC['lv']),
           (L['sqt'], G.SQT.mean(), PALC['sqt']), (L['oth'], G.STRESS_OTHER.mean(), PALC['oth']), (L['mbo'], G.MBO.mean(), PALC['mbo'])]
    fam = sorted([(n, v / 1e3, cc) for n, v, cc in fam], key=lambda t: t[1], reverse=True); total = sum(v for _, v, _ in fam)
    fig = plt.figure(figsize=(9.2, 10.4))
    gsA = fig.add_gridspec(1, 1, left=0.30, right=0.965, top=0.96, bottom=0.585)          # a) as the one-panel figure
    gsB = fig.add_gridspec(1, 2, left=0.085, right=0.975, top=0.495, bottom=0.145, wspace=0.28)   # one legend below b and c
    ax = fig.add_subplot(gsA[0, 0]); y_ = np.arange(len(fam))[::-1]
    for yi, (n, v, cc) in zip(y_, fam):
        ax.barh(yi, v, color=cc, edgecolor="white", height=0.72, zorder=3)
        ax.text(v + total * 0.006, yi, dec(f"{v:.0f} kt  ({100 * v / total:.1f}%)" if v >= 1 else f"{v:.1f} kt  ({100 * v / total:.1f}%)", lang), va="center", ha="left", fontsize=10)
    ax.set_yticks(y_); ax.set_yticklabels([n for n, _, _ in fam]); ax.set_xlim(0, nice(fam[0][1] * 1.25, 50)); ax.set_xlabel(L['f2x'])
    ax.xaxis.grid(True, color="#e8e8e6", lw=0.8, zorder=0); ax.set_axisbelow(True); clean(ax)
    fig.text(0.085, 0.975, 'a) ' + L['f2t'].format(total), fontsize=11.5, fontweight="bold", ha='left', va='bottom')
    mm_ = np.arange(1, 13)
    ax = fig.add_subplot(gsB[0, 0])
    ax.fill_between(mm_, CMIN.BVOC.values / 1e3, CMAX.BVOC.values / 1e3, color=PALC['tot'], alpha=0.15, lw=0, label=L['band'])
    for k, col, nm in (('BVOC', PALC['tot'], L['tot']), ('ISOP', PALC['iso'], L['iso']), ('MT', PALC['mt'], L['mt']),
                       ('OVOC', PALC['ov'], L['ov']), ('LVOC', PALC['lv'], L['lv'])):
        ax.plot(mm_, CLIM[k].values / 1e3, 'o-', ms=3.5, lw=1.6, color=col, label=nm)
    ax.set_xticks(mm_); ax.set_xticklabels(L['mes'], fontsize=8.5); ax.set_ylabel(L['f2y']); ax.set_ylim(0, None); clean(ax)
    ax.yaxis.grid(True, color="#e8e8e6", lw=0.8, zorder=0); ax.set_axisbelow(True)
    ax.set_title(L['f2b'], fontsize=11, fontweight='bold', loc='left'); HL = dict(zip(*ax.get_legend_handles_labels()[::-1]))
    ax = fig.add_subplot(gsB[0, 1])
    for k, col, nm in (('ISOP', PALC['iso'], L['iso']), ('SQT', PALC['sqt'], L['sqt']), ('CO', '#8c6d31', 'CO'), ('MT', PALC['mt'], L['mt']),
                       ('LVOC', PALC['lv'], L['lv']), ('NO', PALC['no'], L['no'])):          # CO in brown: not confused with total BVOC
        ax.plot(mm_, (CLIM[k] / CLIM[k].mean()).values, marker='s' if k == 'CO' else 'o', ms=3.2, lw=1.5, color=col, label=nm,
                ls='--' if k == 'NO' else '-')
    ax.axhline(1, color='0.6', lw=0.9, label=L['ref']); ax.set_xticks(mm_); ax.set_xticklabels(L['mes'], fontsize=8.5); ax.set_ylabel(L['f2yn']); clean(ax)
    ax.yaxis.grid(True, color="#e8e8e6", lw=0.8, zorder=0); ax.set_axisbelow(True)
    for h_, l_ in zip(*ax.get_legend_handles_labels()): HL.setdefault(l_, h_)          # one legend for b and c
    order = [L['tot'], L['band'], L['iso'], L['mt'], L['sqt'], L['ov'], L['lv'], 'CO', L['no'], L['ref']]
    fig.legend([HL[k] for k in order], order, loc='upper center', bbox_to_anchor=(0.53, 0.088), ncol=5, fontsize=8.2, frameon=False,
               columnspacing=1.1, handlelength=2.0, handletextpad=0.5)
    ax.yaxis.set_major_locator(MultipleLocator(0.5)); ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: dec(f'{v:.1f}', lang)))
    ax.set_title(L['f2c'], fontsize=11, fontweight='bold', loc='left')
    fig.savefig(f"{od}/fig2.png", dpi=300, facecolor="white", bbox_inches='tight', pad_inches=0.12); plt.close(fig)
    # Fig 3 - detrended correlations
    Rm = np.array([[COR[(COR.series == k) & (COR.driver == d)].r_detrended.iloc[0] for d in ('TEMP2_C', 'PAR', 'LAI', 'SOILM')] for k in CLASSES])
    fig, ax = plt.subplots(figsize=(6.4, 9.6)); im = ax.imshow(Rm, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i in range(Rm.shape[0]):
        for j in range(Rm.shape[1]):
            ax.text(j, i, dec(f"{Rm[i, j]:.2f}", lang), ha='center', va='center', fontsize=8.5, color='white' if abs(Rm[i, j]) > 0.6 else 'black')
    ax.set_xticks(range(4)); ax.set_xticklabels(L['drv'], rotation=30, ha='right'); ax.set_yticks(range(len(CLASSES))); ax.set_yticklabels(CLASSES)
    cbar = fig.colorbar(im, ax=ax, shrink=0.5); cbar.set_label(L['cb']); ax.set_title(L['f3t'], fontsize=11.5, fontweight='bold')
    fig.tight_layout(); fig.savefig(f"{od}/fig3.png", dpi=300, facecolor="white"); plt.close(fig)
    fig_detrend(lang, od)                                          # article Figure 2 (Methods)
    # Fig 4 - SOA precursors by mesoregion
    S = MESO.sort_values('SOA', ascending=False)
    fig, ax = plt.subplots(figsize=(9.4, 5.6)); fig.subplots_adjust(left=0.22, right=0.965, top=0.90, bottom=0.12)
    yy = np.arange(len(S))[::-1]; mtv = S.MT.values / 1e3; sqv = S.SQT.values / 1e3
    ax.barh(yy, mtv, color=PALC['mt'], edgecolor="white", height=0.68, zorder=3, label=L['mt'])
    ax.barh(yy, sqv, left=mtv, color=PALC['sqt'], edgecolor="white", height=0.68, zorder=3, label=L['sqt'])
    for yi, (idx, r_) in zip(yy, S.iterrows()):
        ax.text(r_.SOA / 1e3 + S.SOA.max() / 1e3 * 0.015, yi, dec(f"{r_.pct_SC:.0f}%  |  {r_.intensity_kg_ha:.0f} kg ha⁻¹", lang), va="center", ha="left", fontsize=10)
    ax.set_yticks(yy); ax.set_yticklabels(list(S.index)); ax.set_xlim(0, nice(S.SOA.max() / 1e3 * 1.35, 10)); ax.set_xlabel(L['f4x'])
    ax.xaxis.grid(True, color="#e8e8e6", lw=0.8, zorder=0); ax.set_axisbelow(True); clean(ax)
    ax.set_title(L['f4t'], fontsize=12.5, fontweight="bold", pad=10); ax.legend(loc="lower right", fontsize=10)
    fig.savefig(f"{od}/fig4.png", dpi=300, facecolor="white"); plt.close(fig)
    # Fig 5 - intensity map by mesoregion
    vals = MESO.intensity_kg_ha.to_dict(); fig, ax = plt.subplots(figsize=(10, 7))
    pc = choropleth(ax, mnames, mgeoms, vals, 'YlOrRd', vmin=0, vmax=np.ceil(max(vals.values()) / 5) * 5, lw=0.6)
    for k, g in zip(mnames, mgeoms):
        pp = g.representative_point(); ax.text(pp.x, pp.y, f"{k}\n{vals[k]:.0f}", ha='center', va='center', fontsize=8, fontweight='bold',
                                              bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.6))
    cbar = fig.colorbar(pc, ax=ax, shrink=0.6); cbar.set_label('kg ha⁻¹ ' + L['yr'])
    ax.set_title(L['f5t'], fontsize=12, fontweight='bold'); fig.savefig(f"{od}/fig5.png", dpi=300, facecolor="white", bbox_inches='tight'); plt.close(fig)
    # Fig 6 - intensity vs native forest
    sl, icp = np.polyfit(x, yv, 1)
    fig, ax = plt.subplots(figsize=(8.6, 5.6)); fig.subplots_adjust(left=0.11, right=0.965, top=0.86, bottom=0.12)
    xs = np.array([x.min() - 2, x.max() + 2]); ax.plot(xs, icp + sl * xs, ls="--", color="#9a9a9a", lw=1.6, zorder=2)
    ax.scatter(x, yv, s=140, color=PALC['mt'], edgecolor="white", lw=0.8, zorder=3)
    for xi, yi, idx in zip(x, yv, MESO.index): ax.annotate(idx, (xi, yi), xytext=(7, 4), textcoords="offset points", fontsize=9)
    ax.set_xlim(x.min() - 4, x.max() + 14); ax.set_xlabel(L['f6x']); ax.set_ylabel(L['f6y']); clean(ax)
    ax.yaxis.grid(True, color="#e8e8e6", lw=0.8, zorder=0); ax.set_axisbelow(True)
    sub = (f"Pearson r = {r:.2f} (IC 95% {ci[0]:.2f}–{ci[1]:.2f}); Spearman ρ = {rho:.2f}" if lang == 'pt'
           else f"Pearson r = {r:.2f} (95% CI {ci[0]:.2f}–{ci[1]:.2f}); Spearman ρ = {rho:.2f}")
    ax.set_title(L['f6t'] + "\n" + dec(sub, lang), fontsize=11.5, fontweight="bold", pad=8)
    fig.savefig(f"{od}/fig6.png", dpi=300, facecolor="white"); plt.close(fig)
    # Fig 7 - municipalities 2025: totals and intensity, Jenks natural breaks
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.2))
    for ax, vals, ttl, cmap in ((axs[0], MU.BVOC.to_dict(), L['f7a'], 'Greens'), (axs[1], MU.intensity_kg_ha.to_dict(), L['f7b'], 'YlGn')):
        b = fisher_jenks(list(vals.values()), 5); choropleth(ax, unames, ugeoms, vals, cmap, bins=b)
        legend_bins(ax, b, cmap, ttl.split(') ')[1], lang); ax.set_title(ttl, fontsize=11)
    fig.suptitle(L['f7t'] + f" · {L['q']}", fontsize=12.5, fontweight='bold')
    fig.savefig(f"{od}/fig7.png", dpi=250, facecolor="white", bbox_inches='tight'); plt.close(fig)
    # Fig 8 - isoprene, CO and soil NO by municipality, 2019 and 2025
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))
    for j, (k, cmap, nm) in enumerate((('ISOP', 'Greens', L['iso']), ('CO', 'Blues', L['co']), ('NO', 'Reds', L['no']))):
        b = fisher_jenks(list(MU19[k].values) + list(MU[k].values), 5)
        for i, (yy_, T_) in enumerate(((2019, MU19), (2025, MU))):
            ax = axs[i, j]; choropleth(ax, unames, ugeoms, T_[k].to_dict(), cmap, bins=b); legend_bins(ax, b, cmap, 't ' + L['yr'], lang)
            ax.set_title(f"{nm}, {yy_}", fontsize=11, fontweight='bold')
    fig.suptitle(L['f8t'] + f" · {L['q']}", fontsize=13, fontweight='bold')
    fig.savefig(f"{od}/fig8.png", dpi=220, facecolor="white", bbox_inches='tight'); plt.close(fig)
    print(lang, 'figures 1-8 ->', od)
