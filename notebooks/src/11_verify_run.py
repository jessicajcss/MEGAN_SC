# %% [markdown]
# # 11 · Verification of a production run
#
# Automatic checks on the outputs of `10_production_run.ipynb`, year by year, for the region and run of a configuration
# file (`RUN_CONFIG`, the same file notebook 10 used). It only reads files, and runs in Colab or on any computer.
#
# | Check | What it tests | Pass criterion |
# |---|---|---|
# | C0 provenance | the year was produced by notebook 10 with the final code (Santa Catarina example: 2024 and 2025 by notebook 17, VIIRS-corrected LAI) | `provenance_v10.json` present, engine patches FIX-6/7/8 listed, same `fixes.py` hash in every year, variant `base` (Santa Catarina 2024-2025: `base` or `lai_viirs`) |
# | C1 cells without meteorology | the signature of a cell computed without valid meteorology (zero light): a forested cell that emits BVOC but no isoprene | no cell of the region with tree fraction > 0.1, BVOC > 0 and isoprene = 0 |
# | C2 outliers | any cell far outside the intensity distribution | no cell of the region above 5 × the 99th percentile of valid-meteorology cells |
# | C3 filled cells | cells whose meteorology was filled (FIX-1) behave like their valid neighbours | emission activity γ = E/(EF·area) of filled cells vs the median of valid cells in a 7×7 window, only where the cell LAI is within ±25 % of the neighbours' LAI; exceedances of the valid-cell 99.5th percentile of \|Δlog γ\| not more than a binomial(n, 0.005) 99 % bound |
# | C4 completeness | 12 months, 12 fix logs, no checkpoint left | all present |
# | C5 PAR | the de-accumulation of the ERA5-Land radiation (FIX-8; the 00 UTC field holds the 24 h of the previous day) | regional mean of the monthly PAR between 20 and 200 W m⁻² |
#
# All checks must PASS before the results are used. C4 needs whole years: a run of part of a year shows FAIL in C4
# for that year, and the other checks still apply.
#
# The table of annual totals multiplies isoprene by the canopy loss and production factor rho = 0.96 of Guenther et al.
# (2006, Sect. 3.3), as the result notebooks do (notebook 10 writes the MEGANv3.21 output, which applies no such factor);
# the checks are ratios or zero tests and do not change with this constant factor.
#
# Why C3 is built this way: the emission factor differs from cell to cell (ecotype and growth-form mix), so comparing
# raw emissions or the isoprene share of filled coastal cells with inland cells mixes land cover with meteorology (in
# the Santa Catarina example, the cells of Babitonga Bay have a low isoprene share because of their EF, not because of
# the fill). Dividing by the cell EF leaves γ, which depends only on meteorology, LAI and canopy; the LAI match removes
# the remaining contrast of coastal LAI pixels.

# %%
__FILE:notebooks/src/_setup_repo.py__
%pip -q install rasterio netCDF4 geopandas pandas matplotlib pyyaml
RUN_CONFIG = 'config/santa_catarina_2003_2025.yaml'           # <-- the configuration of the run (as in notebook 10)
OUT_BASE   = '/content/drive/MyDrive' if IN_COLAB else os.path.dirname(REPO_DIR)   # <-- as in notebook 10
NEW_ROOT   = None                                             # None = OUT_BASE/<output.dir of RUN_CONFIG>; or the run folder
BOUNDARY   = None                                             # None = region boundary of RUN_CONFIG; or a .gpkg/.shp
REPORT_DIR = None                                             # None = NEW_ROOT; where the report CSV/PNG are written
YEARS = None                                                  # None = every year found; or e.g. [2019]
RHO_ISOP = 0.96                                               # isoprene canopy loss and production factor (Guenther et al., 2006)

# %%
import json, yaml
import numpy as np, pandas as pd, rasterio, geopandas as gpd, matplotlib.pyplot as plt
def in_repo(p): return p if (not p or os.path.isabs(p)) else os.path.join(REPO_DIR, p)
CFG = yaml.safe_load(open(in_repo(RUN_CONFIG), encoding='utf-8'))
OUT_TAG = CFG.get('output', {}).get('tag', CFG['run']['name'])
REG = CFG.get('region', {}); REGION_ABBR = REG.get('abbr', CFG['run']['name'])
NEW_ROOT = NEW_ROOT or os.path.join(OUT_BASE, CFG.get('output', {}).get('dir', f"megan_sc_outputs_{CFG['run']['name']}"))
BOUNDARY = in_repo(BOUNDARY or REG.get('boundary') or CFG['domain'].get('shapefile'))
REPORT_DIR = REPORT_DIR or NEW_ROOT
print('run:', NEW_ROOT, '| region:', BOUNDARY)
from scipy.ndimage import generic_filter
from scipy import stats
from netCDF4 import Dataset
from rasterio.features import rasterize
from rasterio.transform import from_origin
CLASSES = ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM', 'NO', 'SQT_HR', 'SQT_LR',
           'MEOH', 'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO']
BVOC = [c for c in CLASSES if c not in ('NO', 'CO')]

# %%
def read_annual(root, y):
    """Annual totals per class (t/yr) on the south-first grid. Uses annual_emis_<y>.nc when present
    (notebook 10), otherwise sums the 12 monthly NetCDFs; isoprene times RHO_ISOP."""
    d = os.path.join(root, str(y), 'emissions'); f = os.path.join(d, f'annual_emis_{y}.nc')
    files = [f] if os.path.exists(f) else sorted(glob.glob(os.path.join(d, f'emis_{y}_[01][0-9].nc')))
    if not files: return None, None, None
    E = {c: 0.0 for c in CLASSES}
    for fn in files:
        with Dataset(fn) as ds:
            lat = np.array(ds.variables['lat'][:]); lon = np.array(ds.variables['lon'][:])
            for c in CLASSES: E[c] = E[c] + np.nan_to_num(np.array(ds.variables[c][:], dtype='float64'))
            attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
    E['ISOP'] = E['ISOP'] * RHO_ISOP
    return E, (lat, lon), attrs
def south_first_tif(path, band=1):
    with rasterio.open(path) as s: return np.flipud(s.read(band).astype('float64'))
_geom = None
def region_mask(lat, lon):
    global _geom
    if _geom is None:
        g = gpd.read_file(BOUNDARY).to_crs('EPSG:4326'); _geom = g.geometry.union_all() if hasattr(g.geometry, 'union_all') else g.geometry.unary_union
    r = lon[1] - lon[0]; tr = from_origin(lon[0] - r / 2, lat[-1] + r / 2, r, r)
    return np.flipud(rasterize([(_geom, 1)], out_shape=(lat.size, lon.size), transform=tr, fill=0, dtype='uint8')).astype(bool)
def nbr_median(field, valid, size=7):
    """Median of the valid cells in a size x size window around each cell, excluding the cell itself (>= 3 needed)."""
    v = np.where(valid, field, np.nan); c = size * size // 2
    def f(w):
        w = np.delete(w, c); w = w[np.isfinite(w)]
        return np.median(w) if w.size >= 3 else np.nan
    return generic_filter(v, f, size=size, mode='constant', cval=np.nan)
def cell_area(lat, lon):
    R = 6371008.8; d = np.deg2rad(lon[1] - lon[0]); h = (lat[1] - lat[0]) / 2
    return (R ** 2 * d * np.abs(np.sin(np.deg2rad(lat + h)) - np.sin(np.deg2rad(lat - h))))[:, None] * np.ones((1, lon.size))

# %%
years = sorted(int(os.path.basename(p)) for p in glob.glob(os.path.join(NEW_ROOT, '[12][09][0-9][0-9]')))
if YEARS: years = [y for y in years if y in YEARS]
print('years found:', years)
rep, tot_rows, shas = [], [], set()
def add(y, check, ok, detail): rep.append(dict(year=y, check=check, result=('PASS' if ok else 'FAIL') if ok is not None else 'INFO', detail=detail))
import time; T0 = time.time()
print('progress (year, minutes since start; files are read from Drive, so a fresh session is slow):', flush=True)
for y in years:
    print(f'{y} ({(time.time() - T0) / 60:.1f} min)', end='  ', flush=True)
    D = os.path.join(NEW_ROOT, str(y))
    pv = os.path.join(D, 'provenance_v10.json')
    if os.path.exists(pv):
        P = json.load(open(pv)); shas.add(P.get('fixes_sha256'))
        ok = all(any(f in p for p in P.get('patches', [])) for f in ('FIX-6', 'FIX-7', 'FIX-8')) and (P.get('variant') == 'base' or (P.get('variant') == 'lai_viirs' and y in (2024, 2025)))
        add(y, 'C0 provenance', ok, f"sha {P.get('fixes_sha256')} | {P.get('patches')} | variant {P.get('variant')} | LAI {P.get('lai_aggregation')}")
    else:
        add(y, 'C0 provenance', False, 'provenance_v10.json missing (not produced by notebook 10)')
    E, ll, attrs = read_annual(NEW_ROOT, y)
    if E is None:
        add(y, 'C4 completeness', False, 'annual_emis missing'); continue
    lat, lon = ll; INSIDE = region_mask(lat, lon); A = cell_area(lat, lon); bv = sum(E[c] for c in BVOC)
    gfp = glob.glob(os.path.join(D, 'inputs', 'growthform_*.npz')); tree = np.load(gfp[0])['tree'] if gfp else np.ones_like(bv)
    c1 = INSIDE & (tree > 0.1) & (bv > 0) & (E['ISOP'] == 0)
    add(y, 'C1 forested cells without isoprene', not c1.any(), f'{int(c1.sum())} cells')
    filled = np.zeros_like(INSIDE); logs = sorted(glob.glob(os.path.join(D, 'inputs', f'fixlog_{y}_??.npz')))
    for f in logs: filled |= np.load(f)['filled_met_steps'] > 0
    inten = bv * 1e3 / (A / 1e4); v = inten[INSIDE & ~filled & (bv > 0)]
    thr = 5 * np.percentile(v, 99); c2 = INSIDE & (inten > thr)
    add(y, 'C2 intensity outliers', not c2.any(), f'{int(c2.sum())} cells > {thr:.0f} kg/ha/yr')
    efp = os.path.join(D, 'inputs', f'ef_{y}.npz')
    laif = sorted(glob.glob(os.path.join(D, 'inputs', f'lai_{y}_[01][0-9].tif')))
    if logs and filled.any() and os.path.exists(efp) and laif:
        Z = np.load(efp); EF = Z['ef'].astype('float64'); ecl = [str(c) for c in Z['classes']]
        lai = np.mean([south_first_tif(f) for f in laif], axis=0)
        laim = nbr_median(lai, INSIDE & ~filled & (lai > 0))
        match = (lai > 0) & (laim > 0) & (np.abs(np.log(np.where(lai > 0, lai, 1) / np.where(laim > 0, laim, 1))) < np.log(1.25))
        det, ok3 = [], True
        for c in ('ISOP', 'MT_PINE', 'MEOH'):
            ef_c = EF[ecl.index(c)]; okc = INSIDE & (ef_c > 0) & (E[c] > 0)
            lg = np.log(np.where(okc, E[c], 1.0) / np.where(okc, ef_c * A, 1.0)); lg[~okc] = np.nan
            d = lg - nbr_median(lg, okc & ~filled)
            v = okc & ~filled & match & np.isfinite(d); fsel = okc & filled & match & np.isfinite(d)
            thr = np.percentile(np.abs(d[v]), 99.5); k = int((np.abs(d[fsel]) > thr).sum()); n = int(fsel.sum())
            lim = int(stats.binom.ppf(0.99, n, 0.005)) if n else 0
            ok3 &= k <= lim; det.append(f'{c} {k}/{n} > {thr:.2f} (limit {lim})')
        nf = int((INSIDE & filled).sum()); nm = int((INSIDE & filled & match).sum())
        add(y, 'C3 filled cells vs valid neighbours', ok3, ' | '.join(det) + f' | {nm} of {nf} filled cells LAI-matched')
    else:
        add(y, 'C3 filled cells vs valid neighbours', None, 'not applicable (no fix logs, no filled cells, or no ef_<year>.npz / LAI maps)')
    nmon = len([f for f in glob.glob(os.path.join(D, 'emissions', f'emis_{y}_??.nc'))])
    ck = glob.glob(os.path.join(D, 'checkpoints', '*.npz'))
    add(y, 'C4 completeness', nmon == 12 and len(logs) == 12 and not ck, f'{nmon} months, {len(logs)} fix logs, {len(ck)} checkpoints')
    par = [np.nanmean(np.where(INSIDE, south_first_tif(os.path.join(D, 'inputs', f'met_{y}_{m:02d}.tif'), 2), np.nan)) for m in range(1, 13)
           if os.path.exists(os.path.join(D, 'inputs', f'met_{y}_{m:02d}.tif'))]
    add(y, f'C5 PAR monthly {REGION_ABBR} mean (W m-2)', bool(par) and 20 < min(par) and max(par) < 200, f'{min(par):.0f}-{max(par):.0f}' if par else 'no met maps')
    row = {'year': y, **{c: float((E[c] * INSIDE).sum()) for c in CLASSES}}
    tot_rows.append(row)
if len(shas) > 1:
    add('all', 'C0 same fixes.py in every year', False, str(shas))
R = pd.DataFrame(rep); R.to_csv(os.path.join(REPORT_DIR, 'verify_v10_report.csv'), index=False)
TOT = pd.DataFrame(tot_rows).set_index('year'); TOT.to_csv(os.path.join(REPORT_DIR, f'annual_totals_{OUT_TAG}_all_classes.csv'))
pd.set_option('display.width', 220); pd.set_option('display.max_colwidth', 120)
print(f'\ndone in {(time.time() - T0) / 60:.1f} min\n'); print(R.to_string(index=False)); print('\nFAILS:', int((R.result == 'FAIL').sum()))

# %%
# Totals of the region by year (kt/yr)
if len(TOT):
    grp = {'BVOC': BVOC, 'Isoprene': ['ISOP'], 'Monoterpenes': ['MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM'], 'Sesquiterpenes': ['SQT_HR', 'SQT_LR'], 'Soil NO': ['NO']}
    fig, axs = plt.subplots(1, len(grp), figsize=(4 * len(grp), 3.4))
    for ax, (g, cl) in zip(axs, grp.items()):
        ax.plot(TOT.index, TOT[cl].sum(1) / 1e3, 'o-', ms=3, color='#1b6ca8')
        ax.set_title(g); ax.set_ylabel('kt yr⁻¹'); ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(REPORT_DIR, 'verify_v10_totals.png'), dpi=200); plt.show()
