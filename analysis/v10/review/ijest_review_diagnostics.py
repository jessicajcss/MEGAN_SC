# -*- coding: utf-8 -*-
"""ijest_review_diagnostics.py - two post hoc diagnostics of the IJEST manuscript (planted forests; LAI of 2024-2025).

Both sections were defined after the published results were seen (POST HOC); they do not change any published number.

1. Planted forests (potential effect of their emission factors). The production run gives the MapBiomas class
   "forest plantation" (Collection 10, class 9) to the tree growth form, with the emission factors (EF) of the local
   MEGAN v3 ecotype. Here the tree EF on the plantation share of each cell is replaced by the EF of the genus Pinus or
   Eucalyptus from the same MEGAN v3 table (EFv240120.csv, VegID "Pinus" and "Eucalyptus"), for 2019:
       EF'_k = EF_k + f_plant * (EF_k(genus) - EF_k(ecotype, tree))
   In MEGANv3.21 the emission of class k in a cell is EF_k times an activity factor that does not depend on EF (canopy
   type, LAI, meteorology and the constant LDFs used here are unchanged), so the 2019 emission of each cell is scaled by
   EF'_k / EF_k. f_plant is the share of class 9 among the 120-m MapBiomas pixels of each model cell
   (data/external/mapbiomas_c10_lulc2019_SC_120m.npz), capped at the tree fraction of the run. The two cases bound the
   effect for a plantation area of Pinus only or Eucalyptus only; the canopy of both stays the broadleaf canopy type.
2. LAI decrease of 2024-2025. Annual mean MODIS LAI of the run (inputs/lai_<y>_<mm>.tif) in state cells grouped by the
   growth form that covers more than 60 % of the cell (2024 growth forms for 2024-2025, as in the run), and the monthly
   state mean as a z-score against the 2003-2023 mean and standard deviation of the same calendar month.

Usage: python analysis/v10/review/ijest_review_diagnostics.py   (writes <MEGAN_SC_OUT>/ijest_review/; published copy in
results/v10/ijest_review; also run by `python reproduce.py ijest`)
"""
import glob, os, sys
import numpy as np
import pandas as pd
import rasterio
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
import _paths as P
import common as C
REPO = P.REPO; sys.path.insert(0, os.path.join(REPO, 'src'))
from megan_sc import efproc, run, domain as dom
from netCDF4 import Dataset

OUT = os.path.join(P.OUT, 'ijest_review'); os.makedirs(OUT, exist_ok=True)
ROWS, LOG = [], []
def put(k, v, note=''):
    ROWS.append(dict(key=k, value=v, note=note)); LOG.append(f'{k:45s} {v}  {note}')

lat, lon = C.grid(); NY, NX = lat.size, lon.size
SC = C.sc_mask(lat, lon); A = C.cell_area_m2(lat, lon)
CL = C.CLASSES

# ------------------------------------------------------------------ 1. planted forests, 2019
Y = 2019
z = np.load(os.path.join(P.RUN, str(Y), 'inputs', f'ef_{Y}.npz')); EF = z['ef'].astype('float64'); assert list(z['classes']) == CL
gf = dict(np.load(os.path.join(P.RUN, str(Y), 'inputs', f'growthform_{Y}.npz')))
cfg = run.load_config(os.path.join(REPO, 'example_config.yaml')); cfg['domain']['mode'] = 'bbox'
for k in ('ef_table', 'speciation_dir', 'ecotype'):
    v = cfg['datasets']['emission_factors'].get(k)
    if v and not os.path.isabs(v): cfg['datasets']['emission_factors'][k] = os.path.join(REPO, v)
D = dom.build_domain(cfg)
assert np.allclose(D['lat'], lat) and np.allclose(D['lon'], lon), 'domain grid differs from the run grid'
efc = cfg['datasets']['emission_factors']
from megan_sc import fixes                                       # FIX-5 of notebook 10: ecotype no-data pixels take the nearest valid ecotype
_eco_f = P.out('_work', 'ijest', 'sc_ecotype_filled.nc')
if not os.path.exists(_eco_f): fixes.fill_ecotype_nearest(efc['ecotype'], _eco_f)
efc['ecotype'] = _eco_f
evt, flon, flat = run._load_ecotype(efc['ecotype']); tab = efproc.load_ef_table(efc['ef_table'])
EF_TREE = np.nan_to_num(dom.regrid_block_mean(efproc.map_ecotype_ef(evt, efproc.ef_by_ecotype(tab, efc['speciation_dir'], 'tree')),
                                              flon, flat, lon, lat, D['resolution_deg']), nan=0.0)
EF_CHK = run.build_ef_grid_multiform(cfg, D, {k: gf[k] for k in ('tree', 'shrub', 'herb', 'crop')})
dmax = float(np.max(np.abs(EF_CHK - EF)) / np.max(np.abs(EF)))
put('ef_rebuild_max_rel_dev', f'{dmax:.2e}', 'EF grid rebuilt from the MEGAN tables vs inputs/ef_2019.npz (must be ~0)')
assert dmax < 1e-5
GEN = {g: tab.loc[g, efproc.EF_COLS].to_numpy(dtype='float64') for g in ('Pinus', 'Eucalyptus')}

# plantation share of each model cell from the 120-m MapBiomas 2019 map
m = np.load(os.path.join(REPO, 'data', 'external', 'mapbiomas_c10_lulc2019_SC_120m.npz'))
cls, tr = m['a'], m['transform']
ny, nx = cls.shape
px_lon = tr[2] + tr[0] * (np.arange(nx) + 0.5); px_lat = tr[5] + tr[4] * (np.arange(ny) + 0.5)
valid = cls > 0
plant = dom.regrid_block_mean(np.where(valid, (cls == 9).astype('float64'), np.nan), px_lon, px_lat, lon, lat, D['resolution_deg'])
TREE_CLASSES = (3, 5, 6, 49, 9)   # forest formation, mangrove, floodable forest, wooded sandbank vegetation, forest plantation
tree120 = dom.regrid_block_mean(np.where(valid, np.isin(cls, TREE_CLASSES).astype('float64'), np.nan), px_lon, px_lat, lon, lat, D['resolution_deg'])
# fractions of the whole cell (pixels outside Brazil / no data count as not plantation)
cntall = dom.regrid_block_mean(np.ones_like(cls, dtype='float64'), px_lon, px_lat, lon, lat, D['resolution_deg'])
plant = np.nan_to_num(plant, nan=0.0); tree120 = np.nan_to_num(tree120, nan=0.0)
ok = SC & (gf['tree'] > 0)
r_tree = np.corrcoef(tree120[SC], gf['tree'][SC])[0, 1]
put('tree_frac_120m_vs_run_r', f'{r_tree:.3f}', 'tree fraction from the 120-m map vs the run (30 m, Earth Engine), state cells')
put('tree_area_120m_over_run', f'{(tree120 * A * SC).sum() / (gf["tree"] * A * SC).sum():.3f}', 'state tree area ratio')
fp = np.minimum(plant, gf['tree'])
put('plantation_area_km2_2019', f'{(fp * A * SC).sum() / 1e6:,.0f}', 'class 9 area in state cells (capped at the tree fraction)')
put('plantation_share_of_tree_area', f'{(fp * A * SC).sum() / (gf["tree"] * A * SC).sum() * 100:.1f}', '%')

def read_year(y):
    with Dataset(os.path.join(P.RUN, str(y), 'emissions', f'annual_emis_{y}.nc')) as ds:
        la = np.array(ds.variables['lat'][:]); assert np.allclose(la, lat)
        E = {}
        for c in CL:
            a = np.array(ds.variables[c][:], dtype='float64'); E[c] = np.where(np.isfinite(a) & (a > -9990), a, 0.0) * C.rho(c)
    return E
E = read_year(Y)
def grp(E):
    MT = sum(E[k] for k in CL if k.startswith('MT_'))
    o = dict(ISOP=E['ISOP'], MT=MT, SQT=E['SQT_HR'] + E['SQT_LR'], MBO=E['MBO'],
             OVOC=sum(E[k] for k in ('MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD')), NO=E['NO'], CO=E['CO'])
    o['BVOC'] = sum(E[k] for k in C.BVOC); return {k: float((v * SC).sum()) for k, v in o.items()}
base = grp(E)
# share of the 2019 emission of each group that the model draws from the plantation share of the cells
share = {}
Ep = {}
for i, c in enumerate(CL):
    with np.errstate(invalid='ignore', divide='ignore'):
        Ep[c] = np.where(EF[i] > 0, E[c] * fp * EF_TREE[i] / EF[i], 0.0)
gp = grp(Ep)
rows = []
for g in base: rows.append(dict(case='baseline', group=g, t=base[g], plantation_t=gp[g], plantation_pct=100 * gp[g] / base[g] if base[g] else np.nan))
for gname, efg in GEN.items():
    Es = {}
    for i, c in enumerate(CL):
        new = EF[i] + fp * (efg[i] - EF_TREE[i])
        with np.errstate(invalid='ignore', divide='ignore'):
            Es[c] = np.where(EF[i] > 0, E[c] * new / EF[i], 0.0)
    s = grp(Es)
    for g in base: rows.append(dict(case=f'plantations as {gname}', group=g, t=s[g], change_pct=100 * (s[g] / base[g] - 1) if base[g] else np.nan))
T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, 'plantation_ef_sensitivity_2019.csv'), index=False)
for _, r in T.iterrows():
    if r.case == 'baseline':
        put(f'base_{r.group}_2019_t', f'{r.t:,.0f}', f'from plantation share: {r.plantation_t:,.0f} t ({r.plantation_pct:.1f} %)')
    else:
        put(f'{r.case.split()[-1]}_{r.group}_2019_change_pct', f'{r.change_pct:+.1f}', f'{r.t:,.0f} t')
for gname, efg in GEN.items():
    put(f'EF_{gname}_ISOP_MT', f'{efg[0]} / {efg[2:7].sum():.2f}', 'MEGAN v3 table, genus row: EF1 and EF3-EF7')
i_eco = (EF_TREE[0] * fp * A * SC).sum() / (fp * A * SC).sum(); m_eco = (EF_TREE[2:7].sum(0) * fp * A * SC).sum() / (fp * A * SC).sum()
put('EF_ecotype_tree_on_plantations_ISOP_MT', f'{i_eco:.2f} / {m_eco:.2f}', 'area-weighted mean over the plantation share')
put('ISOP_MT_ratio_2019_baseline', f'{base["ISOP"] / base["MT"]:.2f}', '')

# ------------------------------------------------------------------ 2. LAI 2024-2025 by dominant growth form
YEARS = list(range(2003, 2026))
def lai_month(y, mth):
    with rasterio.open(os.path.join(P.RUN, str(y), 'inputs', f'lai_{y}_{mth:02d}.tif')) as s:
        return np.flipud(s.read(1).astype('float64'))
def gform(y):
    f = sorted(glob.glob(os.path.join(P.RUN, str(y), 'inputs', 'growthform_*.npz')))[0]; return dict(np.load(f))
rowsL, rowsM = [], []
for y in YEARS:
    g = gform(y)
    dom_ = {'tree': g['tree'] > 0.6, 'herb': g['herb'] > 0.6, 'crop': g['crop'] > 0.6}
    L = np.stack([lai_month(y, k) for k in range(1, 13)])
    L = np.where(L > 0, L, np.nan); ann = np.nanmean(L, 0)
    r = dict(year=y, all=float(np.nanmean(np.where(SC, ann, np.nan))))
    for k, msk in dom_.items(): r[k] = float(np.nanmean(np.where(SC & msk, ann, np.nan))); r[f'n_{k}'] = int((SC & msk & np.isfinite(ann)).sum())
    rowsL.append(r)
    for k in range(12): rowsM.append(dict(year=y, month=k + 1, lai=float(np.nanmean(np.where(SC, L[k], np.nan)))))
LA = pd.DataFrame(rowsL).set_index('year'); LM = pd.DataFrame(rowsM)
LA.to_csv(os.path.join(OUT, 'lai_annual_by_growth_form.csv'))
ref = LA.loc[2003:2023]
for k in ('all', 'tree', 'herb', 'crop'):
    for y in (2024, 2025):
        put(f'lai_{k}_{y}_vs_2003_2023_pct', f'{100 * (LA.loc[y, k] / ref[k].mean() - 1):+.1f}',
            f'z = {(LA.loc[y, k] - ref[k].mean()) / ref[k].std(ddof=1):+.1f}; n cells {int(LA.loc[y, "n_" + k]) if k != "all" else int(SC.sum())}')
    put(f'lai_{k}_min_2003_2023_pct', f'{100 * (ref[k].min() / ref[k].mean() - 1):+.1f}', 'lowest year of 2003-2023 vs its mean')
clim = LM[LM.year <= 2023].groupby('month').lai.agg(['mean', 'std'])
LM['z'] = [(r.lai - clim.loc[r.month, 'mean']) / clim.loc[r.month, 'std'] for r in LM.itertuples()]
LM.to_csv(os.path.join(OUT, 'lai_monthly_state_z.csv'), index=False)
zz = LM[(LM.year >= 2024)].copy()
put('lai_months_z_below_-1_2024_2025', ' '.join(f'{int(r.year)}-{int(r.month):02d}:{r.z:+.1f}' for r in zz.itertuples() if r.z < -1), '')
put('lai_months_z_below_-1_2003_2023', str(int((LM[LM.year <= 2023].z < -1).sum())), 'of 252 months')
pd.DataFrame(ROWS).to_csv(os.path.join(OUT, 'ijest_review_numbers.csv'), index=False)
open(os.path.join(OUT, 'ijest_review.log'), 'w').write('\n'.join(LOG) + '\n')
print('\n'.join(LOG)); print(LA.round(3).to_string())
