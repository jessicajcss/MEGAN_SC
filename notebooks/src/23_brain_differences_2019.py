# %% [markdown]
# # 23 · Why this inventory and BRAIN differ in 2019: tests of each explanation with the data at hand
#
# Notebook 18 found, for Santa Catarina in 2019 (true cell areas for BRAIN): isoprene 0.79 of BRAIN, monoterpenes 0.85,
# sesquiterpenes 0.77, ethanol/acetaldehyde 0.75, CO 0.97, soil NO 1.13; monthly cycles alike (isoprene r = 0.96);
# cells correlated (isoprene r = 0.63) with ratios from 0.57 on the coast (Grande Florianópolis, Sul) to 1.01-1.06 in the
# west and on the plateau (Oeste, Serrana). This notebook tests the candidate explanations one by one and states, for each,
# how much of the difference it can account for and how firmly.
#
# | Block | Question | Data | Kind of evidence |
# |---|---|---|---|
# | A | Settings that differ between this inventory and the distributed MEGANv3.21 scripts: bidirectional factor of ethanol/acetaldehyde (`GAMBD_YN`), soil-moisture factor (`GAMSM_YN`), PAR fraction (`SRAD_TO_PAR` = 0.5) | MEGANv3.21 source and run scripts; this inventory's LAI and fix logs; notebook 20 (`par050`) | exact (GAMBD), approximate (GAMSM), model experiment (PAR) |
# | B | Does the class-by-class pattern of the ratios follow the vegetation (tree share of the emission factor), the temperature sensitivity, or the light dependence of each class? | MEGAN tables of this repository; ratios of notebook 18 | exploratory, n = 12 classes, permutation tests |
# | C | Where in space do the ratios depart, and with what do they co-vary? | BRAIN grid, ≈3,000 cells inside the state | regression with spatial block bootstrap, relative importance (LMG), Moran's I |
# | D | Is the difference fixed in time (land cover, emission factors) or seasonal (meteorology, LAI phenology)? | monthly fields | two-way variance decomposition |
# | E | Does the ratio rise with the elevation of the 4 km cell above the ERA5-Land cell, as unresolved terrain would make it? | notebook 20 elevations, notebook 21 lapse rate | slope against the slope predicted by MEGAN's temperature response |
# | F | How far does each sensitivity run move the ratios? | notebook 20 | model experiments |
# | G | Which spatial pattern does the satellite formaldehyde of 2019 follow better? | TROPOMI (notebook 14 cache) | independent observation, correlations with block bootstrap |
# | H | Soil NO: BDSNP (BRAIN) against the MEGANv3.21 default (here) | notebook 18 | description |
#
# **What is known from the sources** (checked in the MEGANv3.21 distribution, `MEGANv3.21.tar.gz`, Zenodo
# 10.5281/zenodo.10526206, and in Hoinaski et al. 2024, https://doi.org/10.5194/essd-16-2385-2024):
# * `megvea.f` reads `GAMBD_YN` and `GAMSM_YN` with default false, and the distributed `step6.run.megvea.v32.csh` sets
#   both to `N`; this inventory applies both (notebook 12; `run_timestep_fixed`).
# * `readmcip.f` (MET2MGN) multiplies the MCIP shortwave `RGRND` by `SRAD_TO_PAR`, default 0.5; the distributed
#   `step2.run.met2mgn.v32.csh` does not set it. This inventory uses 0.45 (`gee/era5.py`).
# * `megvea.f` reads the light-dependent fraction of the four monoterpene classes from a gridded file (`LDFILE`); this
#   inventory uses the constants of `MEGVEA.EXT` for every class. Not quantified here (no LDF file for this grid).
# * `megsea.f` takes the wilting point from the soil type of the meteorological model (`SLTYP`) when the soil-moisture
#   factor is on.
# * BRAIN: MEGAN 3.2 driven by WRF 4.0 through MCIP (WRF: Noah land surface, Dudhia shortwave, YSU boundary layer,
#   Kain-Fritsch cumulus, MODIS land use; namelist in the supplement of Hoinaski et al. 2024); the article points to the
#   MEGAN distribution site for the input data and does not name the vegetation, LAI or emission-factor files; soil NO by
#   BDSNP (file names). Which switches the BRAIN run used is not stated: block A gives the size of each, conditional on it.
#
# **Run after** notebook 18 (its monthly BRAIN cache is reused; nothing is downloaded again). Notebooks 20, 21 and 14 are
# optional inputs: the blocks that need them are skipped with a message when their outputs are absent. No Earth Engine.

# %%
# 1 · Settings -------------------------------------------------------------------------------------------------------------------
%pip -q install netCDF4 pyproj geopandas shapely rasterio requests statsmodels scipy
import os, sys, re, json, glob, math, itertools, datetime, warnings, hashlib
try:
    from google.colab import drive; drive.mount('/content/drive')
except Exception:
    pass
NB18      = '/content/drive/MyDrive/MEGAN_SC/notebooks/18_brain_comparison_2019.ipynb'   # <-- this repository's notebook 18
REPO_DIR  = '/content/drive/MyDrive/MEGAN_SC'                       # <-- this repository (boundaries, MapBiomas 120 m, outlines)
MEGAN_DIR = '/content/drive/MyDrive/megan-sc'                       # <-- engine repo in Drive (emission factors, MEGAN functions)
ROOT      = '/content/drive/MyDrive/megan_sc_outputs_v10'           # <-- production run (as notebook 18)
OVERRIDES = dict(ROOT=ROOT, BND=REPO_DIR + '/data/boundaries', URL_FILE='/content/drive/MyDrive/BRAIN/BRAIN_output_urls.txt',
                 BRAIN_DIR='/content/drive/MyDrive/BRAIN/MEGAN_SC_4X4', MGN_TGZ='/content/drive/MyDrive/MEGANv3.21.tar.gz', MAKE_FIGS=False)
N20_DIR   = ROOT + '/analysis_v10/sensitivity_2019'                 # notebook 20 (optional)
N20_MODE  = 'full'                                                  # 'full' or 'design' (whichever notebook 20 ran; 'full' preferred)
N21_JSON  = ROOT + '/analysis_v10/era5land_vs_inmet/lapse_rate.json'   # notebook 21 (optional)
TROPOMI   = ROOT + '/analysis_v10/tropomi_hcho/cache_cf0.3_sza70_fire120'   # notebook 14 cache (optional)
OUT23     = ROOT + '/analysis_v10/brain_differences'
N_BOOT, N_PERM, BLOCK = 1000, 100_000, 10                           # bootstrap replicates, permutations, block size (BRAIN cells, 4 km)
SEED = 23
os.makedirs(OUT23, exist_ok=True)

# %% [markdown]
# ## Cell 2 · The state of notebook 18
#
# The code cells 1-13 of notebook 18 (this repository's version) run here in the same namespace, with the paths above;
# the BRAIN monthly fields come from its cache (`analysis_v10/brain_comparison/brain_monthly`), so nothing is downloaded.
# Cell 14 (figures) and cell 15 (provenance) are skipped.

# %%
# 2 · Rebuild the state of notebook 18 -------------------------------------------------------------------------------------------------
os.makedirs(OUT23, exist_ok=True)
NBJ = json.load(open(NB18, encoding='utf-8'))
CELLS18 = [''.join(c['source']) for c in NBJ['cells'] if c['cell_type'] == 'code']
def cell_no(code):
    m = re.search(r'^#\s*(\d+)\s*·', code, re.M)
    return int(m.group(1)) if m else None
G18 = {'__name__': '__nb18__'}
for code in CELLS18:
    k = cell_no(code)
    if k is None or k > 13:
        continue
    exec(compile('\n'.join(l for l in code.split('\n') if not l.lstrip().startswith(('%', '!'))), f'nb18:cell{k}', 'exec'), G18)
    if k == 1:
        G18.update(OVERRIDES); G18['OUT_DIR'] = os.path.join(ROOT, 'analysis_v10', 'brain_comparison')
    print(f'notebook 18, cell {k}: done')
globals().update({k: v for k, v in G18.items() if not k.startswith('__')})
import numpy as np, pandas as pd
from scipy import stats, ndimage
warnings.filterwarnings('ignore', category=RuntimeWarning)
pd.set_option('display.width', 220)
print('months compared:', [int(m) for m in MONTHS.month[MSEL]])
print(TB[['item', 'brain_true_t', 'this_inventory_t', 'ratio_vs_true']].round(3).to_string(index=False))

# %% [markdown]
# ## Cell 3 · Shared tools: fields on the BRAIN grid and statistics
#
# * `to_b(field)`: area-weighted mean of a 0.025° field over each BRAIN cell (the same overlap weights as notebook 18).
# * `lmg(X, y)`: the relative importance of each predictor as its average contribution to R² over all orders of entry
#   (Lindeman, Merenda and Gold; Grömping 2006, https://doi.org/10.18637/jss.v017.i01), exact enumeration of the subsets.
# * `block_boot`: resampling of 40 km × 40 km blocks of BRAIN cells (10 × 10), which keeps the spatial correlation within
#   blocks in every replicate.
# * `moran(z)`: Moran's I with queen contiguity on the BRAIN grid and a permutation p-value.

# %%
# 3 · Tools ---------------------------------------------------------------------------------------------------------------------------------
def to_b(field, valid=None):
    f = np.asarray(field, dtype='float64'); v = np.isfinite(f) if valid is None else (valid & np.isfinite(f))
    num = remap(np.where(v, f, 0.0) * A_O); den = remap(v * A_O)
    return np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
def r2_of(X, y):
    if X.shape[1] == 0:
        return 0.0
    Xc = np.column_stack([np.ones(len(y)), X]); b, *_ = np.linalg.lstsq(Xc, y, rcond=None)
    res = y - Xc @ b
    return 1.0 - res @ res / ((y - y.mean()) @ (y - y.mean()))
def lmg(X, y, groups=None):
    """Shapley decomposition of R2 over predictors (or over groups of columns). Returns (shares dict, R2 full)."""
    p = X.shape[1]; groups = groups or {f'x{j}': [j] for j in range(p)}
    names = list(groups); q = len(names); R = {}
    for mask in range(1 << q):
        cols = [c for i, g in enumerate(names) if mask >> i & 1 for c in groups[g]]
        R[mask] = r2_of(X[:, cols], y)
    out = {}
    for i, g in enumerate(names):
        tot = 0.0
        for mask in range(1 << q):
            if mask >> i & 1:
                continue
            s = bin(mask).count('1'); w = math.factorial(s) * math.factorial(q - s - 1) / math.factorial(q)
            tot += w * (R[mask | (1 << i)] - R[mask])
        out[g] = tot
    return out, R[(1 << q) - 1]
def block_ids(rows, cols, k=BLOCK):
    return (rows // k) * 100000 + (cols // k)
def block_boot(stat, n_items, bid, n=N_BOOT, seed=SEED):
    """Replicates of stat(index array) over block resamples."""
    rng = np.random.default_rng(seed); ids = np.unique(bid); mem = [np.nonzero(bid == b)[0] for b in ids]; out = []
    for _ in range(n):
        sel = np.concatenate([mem[i] for i in rng.integers(0, ids.size, ids.size)])
        out.append(stat(sel))
    return out
def moran(z, rows, cols, n_perm=999, seed=SEED):
    """Moran's I (queen contiguity among the given cells) and its permutation p-value (two-sided)."""
    pos = {(r, c): i for i, (r, c) in enumerate(zip(rows, cols))}; I_, J_ = [], []
    for i, (r, c) in enumerate(zip(rows, cols)):
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                j = pos.get((r + dr, c + dc))
                if (dr or dc) and j is not None:
                    I_.append(i); J_.append(j)
    I_, J_ = np.array(I_), np.array(J_); n = len(z); W = len(I_)
    def mi(v):
        v = v - v.mean(); return n / W * np.sum(v[I_] * v[J_]) / np.sum(v * v)
    obs = mi(np.asarray(z, float)); rng = np.random.default_rng(seed)
    perm = np.array([mi(rng.permutation(z)) for _ in range(n_perm)])
    return obs, (1 + np.sum(np.abs(perm) >= abs(obs))) / (n_perm + 1)
def perm_corr(x, y, n=N_PERM, seed=SEED, method='pearson'):
    f = (lambda a, b: np.corrcoef(a, b)[0, 1]) if method == 'pearson' else (lambda a, b: stats.spearmanr(a, b)[0])
    obs = f(x, y); rng = np.random.default_rng(seed)
    cnt = sum(abs(f(x, rng.permutation(y))) >= abs(obs) for _ in range(n))
    return obs, (cnt + 1) / (n + 1)
import rasterio
def read_sf(p, band=1):
    with rasterio.open(p) as s:
        return np.flipud(s.read(band).astype('float64'))
def met_band(name):
    with rasterio.open(os.path.join(ROOT, str(YEAR), 'inputs', f'met_{YEAR}_01.tif')) as s:
        names = list(s.descriptions)
    return names.index(name) + 1 if name in names else ['TEMP2', 'PAR', 'WINDSPD', 'PRES', 'QV', 'SOILM'].index(name) + 1
LAI_M = np.stack([read_sf(os.path.join(ROOT, str(YEAR), 'inputs', f'lai_{YEAR}_{m:02d}.tif')) for m in range(1, 13)])
MET = {b: np.stack([read_sf(os.path.join(ROOT, str(YEAR), 'inputs', f'met_{YEAR}_{m:02d}.tif'), met_band(b)) for m in range(1, 13)])
       for b in ('TEMP2', 'PAR', 'PRES', 'SOILM')}
FIXL = [np.load(os.path.join(ROOT, str(YEAR), 'inputs', f'fixlog_{YEAR}_{m:02d}.npz')) for m in range(1, 13)]
def state_t(arr_m_y_x, months=None):
    """State tonnes (area-weighted SC fraction of each 0.025° cell) of a (12, NY, NX) field over the compared months."""
    ms = MSEL if months is None else months
    return float((np.asarray(arr_m_y_x)[ms].sum(0) * F_SC_O).sum())
RATIO = TB.set_index('item').ratio_vs_true
print('tools ready; LAI and meteorology of', YEAR, 'loaded')

# %% [markdown]
# ## Cell 4 · A. Switches that differ from the distributed MEGANv3.21 scripts
#
# * **GAMBD** multiplies ethanol/acetaldehyde (class ETOH) by a function of the month's LAI only (`GAMMA_LAIbidir`:
#   0.5 LAI below 2, 1 − 0.0625 (LAI − 2) from 2 to 6, 0.75 above), so removing it is exact: each cell-month of ETOH is
#   divided by the factor of its LAI.
# * **GAMSM** multiplies every class at every step; the fix logs hold its monthly mean per cell, and dividing by it is
#   an approximation (exact if the factor is constant within the month).
# * **PAR fraction 0.5**: from the `par050` run of notebook 20, class by class.
#
# The ratio to BRAIN is recomputed with each change and with all of them together.

# %%
# 4 · A. MEGANv3.21 default switches ---------------------------------------------------------------------------------------------------------
def gambd(lai):
    lai = np.asarray(lai, float)
    return np.where(lai < 2.0, 0.5 * lai, np.where(lai <= 6.0, 1.0 - 0.0625 * (lai - 2.0), 0.75))
GBD = gambd(LAI_M)
ETOH = OT[:, I_C['ETOH']]
ETOH_noBD = np.where(GBD > 0, ETOH / np.where(GBD > 0, GBD, 1.0), 0.0)
f_bd = state_t(ETOH_noBD) / state_t(ETOH)
GSM = np.stack([f['gamsm_mean'] for f in FIXL])
w_isop = OT[:, I_C['ISOP']]
gsm_w = float((w_isop[MSEL] * GSM[MSEL] * F_SC_O).sum() / (w_isop[MSEL] * F_SC_O).sum())
f_sm = {c: state_t(np.where(GSM > 0, OT[:, I_C[c]] / np.where(GSM > 0, GSM, 1.0), 0.0)) / state_t(OT[:, I_C[c]]) for c in CLASSES}
low = float((GSM[MSEL][:, F_SC_O > 0] < 0.99).mean())
print(f'GAMBD: state ETOH without the factor = {f_bd:.3f} x the inventory (emission-weighted mean factor {1 / f_bd:.3f})')
print(f'GAMSM: isoprene-weighted mean monthly factor {gsm_w:.4f}; cell-months below 0.99: {100 * low:.2f} %; '
      f'isoprene without it x{f_sm["ISOP"]:.4f}')
ITEMS = {'ISOP': ['ISOP'], 'MT (as TERP)': ['MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM'], 'SQT (as SQT)': ['SQT_HR', 'SQT_LR'],
         'MEOH': ['MEOH'], 'ACTO': ['ACTO'], 'ETOH': ['ETOH'], 'ACID': ['ACID'], 'LVOC': ['LVOC'], 'OXPROD': ['OXPROD'],
         'STRESS': ['STRESS'], 'OTHER': ['OTHER'], 'CO': ['CO'], 'NO': ['NO']}
N20_OUT = os.path.join(N20_DIR, N20_MODE)
HAVE_N20 = os.path.exists(os.path.join(N20_OUT, 'base'))
def n20_annual(v):
    fs = sorted(glob.glob(os.path.join(N20_OUT, v, f'class_tonnes_{YEAR}_*.npz')))
    return sum(np.load(f)['tonnes'].astype('float64') for f in fs) if fs else None
f_par = {}
if HAVE_N20 and os.path.exists(os.path.join(N20_OUT, 'par050')):
    b0, b1 = n20_annual('base'), n20_annual('par050')
    for it, cl in ITEMS.items():
        f_par[it] = float(sum((b1[I_C[c]] * F_SC_O).sum() for c in cl) / sum((b0[I_C[c]] * F_SC_O).sum() for c in cl))
else:
    print('notebook 20 output (par050) not found in', N20_OUT, ': the PAR column is left empty')
rows = []
for it, cl in ITEMS.items():
    if it not in RATIO.index:
        continue
    r0 = float(RATIO[it]); bd = f_bd if it == 'ETOH' else 1.0
    sm_ = float(np.mean([f_sm[c] for c in cl])) if it != 'NO' else 1.0
    pr = f_par.get(it, np.nan)
    rows.append(dict(item=it, ratio_nb18=r0, x_no_gambd=r0 * bd, x_no_gamsm=r0 * sm_, x_par050=r0 * pr,
                     all_three=r0 * bd * sm_ * (pr if np.isfinite(pr) else 1.0), par_available=np.isfinite(pr)))
TA_ = pd.DataFrame(rows); TA_.to_csv(os.path.join(OUT23, 'A_switches.csv'), index=False)
print('ratio this inventory / BRAIN (true areas) with each MEGANv3.21 default applied here:')
print(TA_.round(3).to_string(index=False))

# %% [markdown]
# ## Cell 5 · B. Class fingerprint
#
# For each of the 12 VOC items (soil NO excluded: a different process), properties computed from the MEGAN tables and
# the inventory's own inputs:
# * `tree_share`: share of the state emission-factor potential Σ_g f_g EF_g that comes from trees (f_g: area share of
#   growth form g in the state, from the MapBiomas growth forms of 2019; EF_g: emission factor of the class for that
#   growth form, area-weighted over the cells where it grows);
# * `gsens`: d ln(potential) / d f_tree when tree cover replaces the other growth forms in proportion (with no shrubs,
#   as here, a linear function of `tree_share`, so its correlations are the same and only `tree_share` is tested);
# * `sT`: d ln γ / dT at 25 °C for a uniform warming (hourly and daily temperatures together), LDF-weighted mix of the
#   light-dependent (`gamtld`) and light-independent (`gamtli`) responses of MEGANv3.21;
# * `LDF`: light-dependent fraction (constants of `MEGVEA.EXT`, as the run).
# Items with several classes are weighted by this inventory's 2019 emissions. Correlations of ln(ratio) with each
# property get permutation p-values (100,000 permutations) because n is 12 and the properties are skewed. The analysis
# is repeated on the ratios after block A (all three switches), to see what remains.

# %%
# 5 · B. Fingerprint ---------------------------------------------------------------------------------------------------------------------------
for _m in [m for m in list(sys.modules) if m.startswith('megan_sc')]:
    del sys.modules[_m]
sys.path.insert(0, MEGAN_DIR + '/src')
from megan_sc import run as mrun, domain as mdom, landcover
from megan_sc.engine import gamma as mg, species as msp
cfg = mrun.load_config(MEGAN_DIR + '/config/example_config.yaml' if os.path.exists(MEGAN_DIR + '/config/example_config.yaml') else MEGAN_DIR + '/example_config.yaml')
cfg['domain']['mode'] = 'bbox'
for k in ('ef_table', 'speciation_dir', 'ecotype'):
    v = cfg['datasets']['emission_factors'].get(k)
    if v and not os.path.isabs(v):
        cfg['datasets']['emission_factors'][k] = os.path.join(MEGAN_DIR, v)
_eco = cfg['datasets']['emission_factors']['ecotype']
_eco_f = os.path.join(os.path.dirname(_eco), 'sc_ecotype_filled.nc')            # FIX-5, as notebooks 10 and 17
if not os.path.exists(_eco_f):
    from megan_sc import fixes as mfixes
    print('FIX-5: ecotype no-data pixels filled:', mfixes.fill_ecotype_nearest(_eco, _eco_f))
cfg['datasets']['emission_factors']['ecotype'] = _eco_f
dmn = mdom.build_domain(cfg)
assert np.allclose(dmn['lat'], LAT) and np.allclose(dmn['lon'], LON), 'engine grid differs from the run grid'
GF0 = dict(np.load(os.path.join(ROOT, str(YEAR), 'inputs', f'growthform_{YEAR}.npz')))
one, zero = np.ones((NY, NX)), np.zeros((NY, NX))
EFG = {g: mrun.build_ef_grid_multiform(cfg, dmn, {k: (one if k == g else zero) for k in landcover.GROWTH_FORMS}) for g in landcover.GROWTH_FORMS}
ef_run = np.load(os.path.join(ROOT, str(YEAR), 'inputs', f'ef_{YEAR}.npz'))['ef'].astype('float64')
ef_re = sum(GF0[g][None] * EFG[g] for g in landcover.GROWTH_FORMS)
print(f'EF rebuilt from the growth-form grids vs ef_{YEAR}.npz of the run: max rel dev {np.nanmax(np.abs(ef_re - ef_run) / np.maximum(ef_run, 1e-12)):.1e}')
WA = A_O * F_SC_O
fshare = {g: float((GF0[g] * WA).sum() / WA.sum()) for g in landcover.GROWTH_FORMS}
def ef_eff(g, c):
    w = GF0[g] * WA
    return float((EFG[g][I_C[c]] * w).sum() / w.sum()) if w.sum() > 0 else 0.0
EFE = {(g, c): ef_eff(g, c) for g in landcover.GROWTH_FORMS for c in CLASSES}
def s_T(c, T=298.15):
    s = I_C[c]
    f = lambda x: np.log(mg.gamtld(np.array([x]), np.array([x]), s))[0]
    ld = (f(T + 0.5) - f(T - 0.5)); li = msp.BETA[s]
    return msp.LDF[s] * ld + (1 - msp.LDF[s]) * li
ann_cls = {c: state_t(OT[:, I_C[c]]) for c in CLASSES}
prop = []
for it, cl in ITEMS.items():
    if it == 'NO' or it not in RATIO.index:
        continue
    efg = {g: sum(EFE[(g, c)] for c in cl) for g in landcover.GROWTH_FORMS}
    P = sum(fshare[g] * efg[g] for g in efg); nt = [g for g in efg if g != 'tree']
    other = sum(fshare[g] * efg[g] for g in nt) / sum(fshare[g] for g in nt)
    w = np.array([ann_cls[c] for c in cl]); w = w / w.sum()
    prop.append(dict(item=it, EF_tree=efg['tree'], EF_shrub=efg['shrub'], EF_herb=efg['herb'], EF_crop=efg['crop'],
                     tree_share=fshare['tree'] * efg['tree'] / P, gsens=(efg['tree'] - other) / P,
                     sT=float(np.dot(w, [s_T(c) for c in cl])), LDF=float(np.dot(w, [msp.LDF[I_C[c]] for c in cl])),
                     lnR=math.log(float(RATIO[it])), lnR_after_A=math.log(float(TA_.set_index('item').loc[it, 'all_three']))))
FP = pd.DataFrame(prop); FP.to_csv(os.path.join(OUT23, 'B_fingerprint_properties.csv'), index=False)
print('growth-form shares of the state:', {k: round(v, 3) for k, v in fshare.items()})
print(FP.round(4).to_string(index=False))
fp_rows = []
for ycol in ('lnR', 'lnR_after_A'):
    for xcol in ('tree_share', 'sT', 'LDF'):
        rp, pp = perm_corr(FP[xcol].values, FP[ycol].values)
        rs, ps = perm_corr(FP[xcol].values, FP[ycol].values, method='spearman', n=20_000)
        fp_rows.append(dict(response=ycol, property=xcol, pearson=rp, p_perm=pp, spearman=rs, p_perm_spearman=ps, n=len(FP)))
FPS = pd.DataFrame(fp_rows); FPS.to_csv(os.path.join(OUT23, 'B_fingerprint_correlations.csv'), index=False)
print(FPS.round(4).to_string(index=False))
import statsmodels.api as sm
for ycol in ('lnR', 'lnR_after_A'):
    for xs in (['sT'], ['tree_share'], ['sT', 'tree_share'], ['sT', 'tree_share', 'LDF']):
        o = sm.OLS(FP[ycol].values, sm.add_constant(FP[xs].values)).fit()
        loo = [sm.OLS(np.delete(FP[ycol].values, i), sm.add_constant(np.delete(FP[xs].values, i, axis=0))).fit().params[1] for i in range(len(FP))]
        print(f'{ycol:12s} ~ {" + ".join(xs):18s} R2 {o.rsquared:.2f} adj {o.rsquared_adj:.2f} | coefficients {np.round(o.params[1:], 3)} '
              f'(p {np.round(o.pvalues[1:], 3)}) | leave-one-out range of the first slope {min(loo):.3f} to {max(loo):.3f}')
UNI = FP[(FP[['EF_tree', 'EF_herb', 'EF_crop']].std(axis=1) / FP[['EF_tree', 'EF_herb', 'EF_crop']].mean(axis=1)) < 1e-6]
if len(UNI) >= 3:
    sl = np.polyfit(UNI.LDF, np.exp(UNI.lnR), 1)
    print(f'items with the same EF for trees, herbs and crops ({", ".join(UNI.item)}): ratio = {sl[1]:.3f} {sl[0]:+.3f} x LDF '
          f'(n = {len(UNI)}; light-independent end {sl[1]:.3f}, light-dependent end {sl[1] + sl[0]:.3f}; descriptive only)')

# %% [markdown]
# ## Cell 6 · C. Where the ratios depart: cell-level regression of ln(ratio)
#
# Response: ln(flux here / flux BRAIN) of isoprene (and monoterpenes, soil NO) on the BRAIN cells that lie fully inside
# the state and have both fluxes positive (the cells of notebook 18, cell 13). Predictors, all on the BRAIN grid:
# growth-form fractions (tree, herb, crop) and, within trees, the silviculture and forest-formation shares (MapBiomas 2019,
# 120 m); annual means of temperature, PAR and soil moisture used by this inventory; annual mean LAI; elevation and its
# departure from the ERA5-Land cell (notebook 20, if present) and the ERA5-Land surface pressure (terrain of the ERA5-Land
# cells, always available); distance to the Atlantic; share of cells whose meteorology
# was filled from a neighbour (FIX-1). Standardised predictors; OLS with spatial block-bootstrap intervals; LMG relative
# importance by predictor and by group; variance inflation factors; Moran's I of the response and of the residuals.
# The vegetation fractions nearly sum to one, so their individual coefficients are not separately interpretable (high
# VIF); the share of R² of the vegetation group is. A residual Moran's I well above 0 means a spatially organised factor
# is missing from the predictors (for example the difference between the WRF and ERA5-Land meteorology, which cannot be
# observed here).

# %%
# 6 · C. Cell regression -------------------------------------------------------------------------------------------------------------------------
import geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import from_origin
TRo = from_origin(LON[0] - RES / 2, LAT[-1] + RES / 2, RES, RES)
land = gpd.read_file(REPO_DIR + '/data/external/ibge_america_do_sul_2021_simplified.geojson').to_crs('EPSG:4326')
LAND = np.flipud(rasterize([(g, 1) for g in land.geometry], out_shape=(NY, NX), transform=TRo, fill=0, dtype='uint8')).astype(bool)
dx_km = RES * 111.32 * np.cos(np.deg2rad(LAT.mean())); dy_km = RES * 110.57
DIST = ndimage.distance_transform_edt(LAND, sampling=(dy_km, dx_km))
Z120 = np.load(REPO_DIR + '/data/external/mapbiomas_c10_lulc2019_SC_120m.npz')
def class_shares(Z, lon, lat, res, groups):
    A, T6 = Z['a'], Z['transform']
    ix = np.floor((T6[2] + T6[0] * (np.arange(A.shape[1]) + 0.5) - (lon[0] - res / 2)) / res).astype(np.int64)
    iy = np.floor((T6[5] + T6[4] * (np.arange(A.shape[0]) + 0.5) - (lat[0] - res / 2)) / res).astype(np.int64)
    ok = (iy[:, None] >= 0) & (iy[:, None] < lat.size) & (ix[None, :] >= 0) & (ix[None, :] < lon.size) & (A > 0)
    flat_ = (iy[:, None] * lon.size + ix[None, :])[ok]; a_ok = A[ok]; n = lat.size * lon.size
    n_all = np.bincount(flat_, minlength=n).reshape(lat.size, lon.size)
    return {g: np.bincount(flat_[np.isin(a_ok, codes)], minlength=n).reshape(lat.size, lon.size) / np.maximum(n_all, 1)
            for g, codes in groups.items()}, n_all
SH, N120 = class_shares(Z120, LON, LAT, RES, {'silv': [9], 'forest': [3]})
FILLED = (sum(f['filled_met_steps'] for f in FIXL) > 0).astype(float)
STAT20 = os.path.join(N20_DIR, f'static_{YEAR}.npz')
PRED = {'tree': GF0['tree'], 'herb': GF0['herb'], 'crop': GF0['crop'], 'silviculture': SH['silv'], 'forest_formation': SH['forest'],
        'T_mean': np.nanmean(MET['TEMP2'], 0), 'PAR_mean': np.nanmean(MET['PAR'], 0), 'soilm_mean': np.nanmean(MET['SOILM'], 0),
        'LAI_mean': np.nanmean(np.where(LAI_M > 0, LAI_M, np.nan), 0), 'coast_km': DIST, 'met_filled': FILLED,
        'pres_mean': np.nanmean(MET['PRES'], 0)}                  # ERA5-Land surface pressure: terrain of the ERA5-Land cells
if os.path.exists(STAT20):
    S20 = np.load(STAT20); PRED['elevation'] = S20['z_fine']; PRED['dz_vs_era5'] = S20['dz']
else:
    print('notebook 20 static file not found: elevation predictors skipped')
PB = {k: to_b(v) for k, v in PRED.items()}
GROUPS_C = {'vegetation': ['tree', 'herb', 'crop', 'silviculture', 'forest_formation'], 'LAI': ['LAI_mean'],
            'meteorology': ['T_mean', 'PAR_mean', 'soilm_mean'], 'terrain': [k for k in ('pres_mean', 'elevation', 'dz_vs_era5') if k in PB],
            'coast': ['coast_km', 'met_filled']}
GROUPS_C = {g: v for g, v in GROUPS_C.items() if v}
rr_, cc_ = np.nonzero(CELL)
def cell_table(label):
    fo_, fb_ = FLUX[label]; y = np.log(fo_[CELL] / fb_[CELL])
    X = np.column_stack([PB[k][CELL] for k in PB]); ok = np.isfinite(y) & np.isfinite(X).all(1)
    return y[ok], X[ok], rr_[ok], cc_[ok]
C_ROWS, LMG_ROWS, VIF_ROWS, MORAN_ROWS = [], [], [], []
for label in ('ISOP', 'MT (as TERP)', 'NO'):
    y, X, rw, cl = cell_table(label); names = list(PB)
    Xs = (X - X.mean(0)) / X.std(0)
    o = sm.OLS(y, sm.add_constant(Xs)).fit()
    bid = block_ids(rw, cl)
    reps = np.array(block_boot(lambda s: np.linalg.lstsq(np.column_stack([np.ones(s.size), Xs[s]]), y[s], rcond=None)[0][1:], len(y), bid))
    lo, hi = np.percentile(reps, [2.5, 97.5], axis=0)
    for j, nm in enumerate(names):
        C_ROWS.append(dict(item=label, predictor=nm, coef_per_sd=o.params[j + 1], boot_lo=lo[j], boot_hi=hi[j], n=len(y), r2=o.rsquared))
    sh, r2 = lmg(Xs, y, {nm: [j] for j, nm in enumerate(names)})
    gidx = {g: [names.index(v) for v in vs] for g, vs in GROUPS_C.items()}
    gsh, _ = lmg(Xs, y, gidx)
    greps = np.array([[d[g] for g in gidx] for d in (lmg(Xs[s], y[s], gidx)[0] for s in
                     (np.concatenate([np.nonzero(bid == b)[0] for b in bs]) for bs in
                      (np.random.default_rng(SEED + k).choice(np.unique(bid), np.unique(bid).size) for k in range(200))))])
    glo, ghi = np.percentile(greps, [2.5, 97.5], axis=0)
    for j, nm in enumerate(names):
        LMG_ROWS.append(dict(item=label, level='predictor', name=nm, r2_share=sh[nm], r2_total=r2))
    for k, g in enumerate(gidx):
        LMG_ROWS.append(dict(item=label, level='group', name=g, r2_share=gsh[g], boot_lo=glo[k], boot_hi=ghi[k], r2_total=r2))
    for j, nm in enumerate(names):
        others = np.delete(Xs, j, axis=1); VIF_ROWS.append(dict(item=label, predictor=nm, vif=1.0 / max(1.0 - r2_of(others, Xs[:, j]), 1e-12)))
    mi_y, p_y = moran(y, rw, cl); mi_e, p_e = moran(o.resid, rw, cl)
    print(f'{label}: n = {len(y)} cells, R2 = {r2:.3f}; Moran I of ln(ratio) {mi_y:.3f} (p {p_y:.3f}), of the residuals {mi_e:.3f} (p {p_e:.3f})')
    print('   LMG by group:', {g: f'{gsh[g]:.3f} [{glo[k]:.3f}, {ghi[k]:.3f}]' for k, g in enumerate(gidx)})
    MORAN_ROWS.append(dict(item=label, n=len(y), r2=r2, moran_ln_ratio=mi_y, p_ln_ratio=p_y, moran_residuals=mi_e, p_residuals=p_e))
CR = pd.DataFrame(C_ROWS); CR.to_csv(os.path.join(OUT23, 'C_cell_regression.csv'), index=False)
LM = pd.DataFrame(LMG_ROWS); LM.to_csv(os.path.join(OUT23, 'C_lmg.csv'), index=False)
VF = pd.DataFrame(VIF_ROWS); VF.to_csv(os.path.join(OUT23, 'C_vif.csv'), index=False)
MO = pd.DataFrame(MORAN_ROWS); MO.to_csv(os.path.join(OUT23, 'C_moran.csv'), index=False)
print(CR[CR.item == 'ISOP'].round(3).to_string(index=False))
print(LM[(LM.item == 'ISOP') & (LM.level == 'predictor')].sort_values('r2_share', ascending=False).round(3).to_string(index=False))
print(VF[VF.item == 'ISOP'].round(1).to_string(index=False))

# %% [markdown]
# ## Cell 7 · D. Fixed in time or seasonal? Two-way decomposition of the monthly ratios
#
# For the cells of cell 6 with positive isoprene in every compared month: L(c, m) = ln(flux here / flux BRAIN) =
# μ + a(c) + b(m) + e(c, m). The share of the sum of squares in a(c) (time-invariant spatial pattern: land cover, emission
# factors, and any lasting climate difference), in b(m) (statewide seasonal difference: meteorology, LAI phenology) and
# in e (cell-specific seasonal differences).

# %%
# 7 · D. Seasonal decomposition ----------------------------------------------------------------------------------------------------------------
def monthly_flux(s):
    mw = MECH_MW[s]
    fo = OB[:, I_M[s]] * mw * 1e-3 / (G['a_true'][None, :, None] * 1e-6)
    fb = BR[:, I_M[s]] * mw * 1e-3 / (G['a_nom'] * 1e-6)
    return fo, fb
D_ROWS = []
for label, s in (('ISOP', 'ISOP'), ('MT (as TERP)', 'TERP'), ('MEOH', 'MEOH')):
    fo, fb = monthly_flux(s); ms = np.nonzero(MSEL)[0]
    L = np.log(fo[ms][:, CELL] / fb[ms][:, CELL]).T                       # cells x months
    L = L[np.isfinite(L).all(1)]
    mu = L.mean(); a = L.mean(1) - mu; b = L.mean(0) - mu; e = L - mu - a[:, None] - b[None, :]
    ss = ((L - mu) ** 2).sum(); sa = L.shape[1] * (a ** 2).sum(); sb = L.shape[0] * (b ** 2).sum(); se = (e ** 2).sum()
    D_ROWS.append(dict(item=label, cells=L.shape[0], months=L.shape[1], share_spatial_fixed=sa / ss, share_seasonal_statewide=sb / ss,
                       share_cell_seasonal=se / ss, **{f'month_effect_{m + 1:02d}': float(np.exp(b[i])) for i, m in enumerate(ms)}))
DD = pd.DataFrame(D_ROWS); DD.to_csv(os.path.join(OUT23, 'D_seasonal_decomposition.csv'), index=False)
print(DD.round(3).to_string(index=False))

# %% [markdown]
# ## Cell 8 · E. Unresolved terrain
#
# If BRAIN's 4 km temperature followed the terrain while ERA5-Land holds the mean of its 9 km cell, a BRAIN cell standing
# Δz above its ERA5-Land cell would be colder by Γ Δz, and the ratio here/BRAIN would rise by about s_T Γ Δz.
# Predicted slope of ln(ratio) on Δz: s_T Γ, with s_T from the `lapse` run of notebook 20 (cell response of isoprene to
# ΔT) or, without it, 0.158 K⁻¹ (uniform warming of the leaf-level response at 25 °C); Γ from notebook 21 or 6.5 K km⁻¹.
# Observed slope: OLS of ln(ratio) on Δz alone and with the vegetation predictors, block-bootstrap intervals.

# %%
# 8 · E. Terrain ----------------------------------------------------------------------------------------------------------------------------------
if 'dz_vs_era5' in PB:
    J21 = json.load(open(N21_JSON)) if os.path.exists(N21_JSON) else {}
    gam = J21.get('gamma_for_notebook20_K_per_km', 6.5)
    fit21 = next((f for f in J21.get('fits', []) if f.get('label') == J21.get('chosen')), {})
    g_lo, g_hi = fit21.get('boot_lo', np.nan), fit21.get('boot_hi', np.nan)     # bootstrap interval of the station lapse rate
    sj = os.path.join(N20_OUT, 'lapse_temperature_sensitivity.json')
    sT_ = json.load(open(sj))['isop_dlnE_dT_per_K'] if os.path.exists(sj) else 0.158
    pred = sT_ * gam / 1000.0
    y, X, rw, cl = cell_table('ISOP'); names = list(PB); jz = names.index('dz_vs_era5')
    veg_cols = [names.index(k) for k in ('tree', 'herb', 'crop', 'silviculture', 'forest_formation')]
    def slope(s, cols):
        Xc = np.column_stack([np.ones(s.size), X[s][:, [jz] + cols]]); return np.linalg.lstsq(Xc, y[s], rcond=None)[0][1]
    bid = block_ids(rw, cl); idx = np.arange(len(y))
    E_ROWS = []
    for lab, cols in (('dz alone', []), ('dz + vegetation', veg_cols)):
        b_ = slope(idx, cols); reps = block_boot(lambda s: slope(s, cols), len(y), bid)
        lo, hi = np.percentile(reps, [2.5, 97.5])
        E_ROWS.append(dict(model=lab, observed_per_km=1000 * b_, lo=1000 * lo, hi=1000 * hi, predicted_per_km=1000 * pred,
                           predicted_lo=sT_ * g_lo, predicted_hi=sT_ * g_hi, sT=sT_, gamma=gam,
                           gamma_implied=1000 * b_ / sT_, gamma_implied_lo=1000 * lo / sT_, gamma_implied_hi=1000 * hi / sT_,
                           dz_sd_m=float(np.std(X[:, jz]))))
    ED = pd.DataFrame(E_ROWS); ED.to_csv(os.path.join(OUT23, 'E_terrain.csv'), index=False)
    print(ED.round(3).to_string(index=False))
    print('The terrain explanation is supported if the observed interval contains the predicted slope and excludes 0. '
          'predicted_lo/hi carry the bootstrap interval of the station lapse rate (notebook 21); gamma_implied is the lapse rate '
          'that the observed slope corresponds to (observed / s_T).')
else:
    print('E skipped: run notebook 20 first (elevation fields)')

# %% [markdown]
# ## Cell 9 · F. The sensitivity runs of notebook 20 against BRAIN
#
# Each variant changes this inventory by a factor v/base, computed per cell (on the BRAIN grid), per mesoregion and for
# the state from the runs of notebook 20; applied to the ratios of notebook 18. With `N20_MODE = 'design'` the factors
# come from 28 days and are indicative only. "Gap closed" is 1 − ln(ratio after) / ln(ratio before) for the state, and
# 1 − SD after / SD before for the spread of the mesoregion ratios.

# %%
# 9 · F. Sensitivity runs ---------------------------------------------------------------------------------------------------------------------------
if HAVE_N20:
    variants = sorted(d for d in os.listdir(N20_OUT) if os.path.isdir(os.path.join(N20_OUT, d)) and d != 'met_cache')
    base = n20_annual('base'); F_ROWS = []; CF20 = {}
    meso_r0 = TR_[TR_.item == 'ISOP'].set_index('mesoregion').ratio
    fo0, fb0 = FLUX['ISOP']
    for v in variants:
        ev = n20_annual(v)
        if ev is None:
            continue
        for it, cl in ITEMS.items():
            if it not in RATIO.index:
                continue
            fct = sum((ev[I_C[c]] * F_SC_O).sum() for c in cl) / sum((base[I_C[c]] * F_SC_O).sum() for c in cl)
            r1 = float(RATIO[it]) * fct
            F_ROWS.append(dict(variant=v, level='state', name=it, factor=fct, ratio=r1,
                               gap_closed=(1 - math.log(r1) / math.log(float(RATIO[it]))) if abs(math.log(float(RATIO[it]))) > 1e-3 else np.nan))
        ev_b = remap(ev[I_C['ISOP']]); e0_b = remap(base[I_C['ISOP']])           # tonnes per BRAIN cell (mass-conserving)
        cf = np.where(e0_b > 0, ev_b / np.where(e0_b > 0, e0_b, 1), np.nan)
        mr = {}
        for k, nm in enumerate(MESO):
            w = F_MESO_B[k]; fct = float((ev_b * w).sum() / (e0_b * w).sum()); mr[nm] = meso_r0[nm] * fct
            F_ROWS.append(dict(variant=v, level='mesoregion', name=nm, factor=fct, ratio=mr[nm]))
        sd0 = np.std(np.log(meso_r0.values)); sd1 = np.std(np.log(np.array([mr[n] for n in meso_r0.index])))
        F_ROWS.append(dict(variant=v, level='mesoregion spread', name='SD of ln(ratio)', factor=sd1 / sd0, ratio=sd1, gap_closed=1 - sd1 / sd0))
        CF20[v] = cf
        fo1 = fo0 * cf; a, b = fo1[CELL], fb0[CELL]; pos = (a > 0) & (b > 0); rt = a[pos] / b[pos]
        F_ROWS.append(dict(variant=v, level='cells', name='median ratio', ratio=float(np.median(rt))))
        F_ROWS.append(dict(variant=v, level='cells', name='within factor 2 (%)', ratio=float(100 * np.mean((rt >= 0.5) & (rt <= 2)))))
        F_ROWS.append(dict(variant=v, level='cells', name='r of log fluxes', ratio=float(np.corrcoef(np.log(a[pos]), np.log(b[pos]))[0, 1])))
    FV = pd.DataFrame(F_ROWS); FV.to_csv(os.path.join(OUT23, f'F_sensitivity_vs_brain_{N20_MODE}.csv'), index=False)
    if 'lapse' in CF20 and 'dz_vs_era5' in PB:
        # block E against the lapse run: slope of each series on dz (BRAIN cells), 95 % block-bootstrap intervals
        zb = PB['dz_vs_era5']; okl = CELL & (fo0 > 0) & (fb0 > 0) & np.isfinite(CF20['lapse']) & (CF20['lapse'] > 0) & np.isfinite(zb)
        y0 = np.log(fo0[okl] / fb0[okl]); yc = np.log(CF20['lapse'][okl]); zz = zb[okl] / 1000.0
        rl, cl_ = np.nonzero(okl); bidl = block_ids(rl, cl_); L_ROWS = []
        for lab, yv in (('ln(ratio) before the lapse run', y0), ('ln(lapse / base)', yc), ('ln(ratio) after the lapse run', y0 + yc)):
            sl_ = lambda s, yv=yv: np.polyfit(zz[s], yv[s], 1)[0]
            lo, hi = np.percentile(block_boot(sl_, zz.size, bidl), [2.5, 97.5])
            L_ROWS.append(dict(series=lab, cells=int(zz.size), slope_per_km=sl_(np.arange(zz.size)), lo=lo, hi=hi,
                               r_with_dz=float(np.corrcoef(yv, zz)[0, 1]), sd=float(yv.std())))
        FL = pd.DataFrame(L_ROWS); FL.to_csv(os.path.join(OUT23, f'F_lapse_vs_terrain_{N20_MODE}.csv'), index=False)
        print('isoprene, BRAIN cells: dependence on dz (BRAIN cell minus ERA5-Land cell elevation) before and after the lapse run')
        print(FL.round(3).to_string(index=False))
    print(FV[FV.level == 'state'].pivot(index='name', columns='variant', values='ratio').round(3).to_string())
    print(FV[FV.level == 'mesoregion'].pivot(index='name', columns='variant', values='ratio').round(3).to_string())
    print(FV[FV.level.isin(['mesoregion spread', 'cells'])].pivot(index='name', columns='variant', values='ratio').round(3).to_string())
else:
    print('F skipped: notebook 20 output not found in', N20_OUT)

# %% [markdown]
# ## Cell 10 · G. TROPOMI formaldehyde of 2019 as an independent check of the spatial pattern
#
# Formaldehyde columns follow isoprene emission with a displacement of some tens of kilometres and carry a background
# from methane oxidation, so the comparison uses the summer months of 2019 (January, February, December), the column
# above the monthly background of notebook 14, emissions smoothed with a 20 km Gaussian, and blocks of 5 × 5 BRAIN cells
# (20 km) with at least half of their area in the state. Both inventories are compared with the same TROPOMI field:
# the difference of their correlations gets a paired block-bootstrap interval (blocks of 3 × 3 aggregated cells). The
# partial correlation removes the part of each variable explained linearly by surface pressure (terrain: a column over
# high ground holds less air). The monthly cycle of the state mean column is correlated with the monthly isoprene of each
# inventory (12 months). This is a check of the pattern, not of the magnitude.

# %%
# 10 · G. TROPOMI ---------------------------------------------------------------------------------------------------------------------------------
if os.path.isdir(TROPOMI):
    SUMMER = [m for m in (1, 2, 12) if MSEL[m - 1]]
    def hcho(m):
        z = np.load(os.path.join(TROPOMI, f'hcho_{YEAR}_{m:02d}.npz')); return z['H'].astype('float64') - float(z['bg'])
    Hs = np.nanmean(np.stack([to_b(hcho(m)) for m in SUMMER]), 0)
    fo, fb = monthly_flux('ISOP')
    sig = 20.0 / 4.0
    Eo = ndimage.gaussian_filter(np.nan_to_num(fo[[m - 1 for m in SUMMER]].mean(0)), sig)
    Eb = ndimage.gaussian_filter(np.nan_to_num(fb[[m - 1 for m in SUMMER]].mean(0)), sig)
    Pb = to_b(np.nanmean(MET['PRES'][[m - 1 for m in SUMMER]], 0))
    def blk(a, k=5):
        ny, nx = (a.shape[0] // k) * k, (a.shape[1] // k) * k
        return np.nanmean(a[:ny, :nx].reshape(ny // k, k, nx // k, k), axis=(1, 3))
    fr = blk(F_SC_B); xo, xb, h, pr = blk(Eo), blk(Eb), blk(Hs), blk(Pb)
    ok = (fr >= 0.5) & np.isfinite(h) & np.isfinite(pr)
    r_, c_ = np.nonzero(ok); xo, xb, h, pr = xo[ok], xb[ok], h[ok], pr[ok]
    resid = lambda v: v - np.polyval(np.polyfit(pr, v, 1), pr)
    def part(a, s):
        return a[s] - np.polyval(np.polyfit(pr[s], a[s], 1), pr[s])
    STATS = {'summer 2019, Pearson r with TROPOMI': lambda s: (np.corrcoef(xo[s], h[s])[0, 1], np.corrcoef(xb[s], h[s])[0, 1]),
             'summer 2019, Spearman with TROPOMI': lambda s: (stats.spearmanr(xo[s], h[s])[0], stats.spearmanr(xb[s], h[s])[0]),
             'summer 2019, partial r (surface pressure removed)': lambda s: (np.corrcoef(part(xo, s), part(h, s))[0, 1],
                                                                             np.corrcoef(part(xb, s), part(h, s))[0, 1])}
    bid = block_ids(r_, c_, 3); allx = np.arange(len(h)); G_ROWS = []
    for lab, f in STATS.items():
        a0, b0 = f(allx); reps = [a - b for a, b in block_boot(f, len(h), bid, n=2000)]
        G_ROWS.append(dict(test=lab, this_inventory=a0, brain=b0, difference=a0 - b0, lo=np.percentile(reps, 2.5),
                           hi=np.percentile(reps, 97.5), n_blocks=len(h)))
    G_ROWS.append(dict(test='surface pressure vs TROPOMI and vs each inventory (r)', this_inventory=np.corrcoef(pr, xo)[0, 1],
                       brain=np.corrcoef(pr, xb)[0, 1], difference=np.corrcoef(pr, h)[0, 1], n_blocks=len(h)))
    SER = pd.DataFrame({'month': range(1, 13), 'hcho_state': [float(np.nansum(to_b(hcho(m)) * F_SC_B) / F_SC_B[np.isfinite(to_b(hcho(m)))].sum()) for m in range(1, 13)],
                        'isop_here_t': [float((OB[m, I_M['ISOP']] * F_SC_B).sum()) for m in range(12)],
                        'isop_brain_t': [float((BR[m, I_M['ISOP']] * F_SC_B * AREA_RATIO[:, None]).sum()) for m in range(12)]})
    SER = SER[MSEL]
    G_ROWS.append(dict(test='monthly cycle 2019, r with the state mean column', this_inventory=np.corrcoef(SER.isop_here_t, SER.hcho_state)[0, 1],
                       brain=np.corrcoef(SER.isop_brain_t, SER.hcho_state)[0, 1], n_blocks=len(SER)))
    GT = pd.DataFrame(G_ROWS); GT.to_csv(os.path.join(OUT23, 'G_tropomi_2019.csv'), index=False)
    print(GT.round(3).to_string(index=False))
else:
    print('G skipped: TROPOMI cache not found:', TROPOMI)

# %% [markdown]
# ## Cell 11 · H. Soil NO: where each inventory puts it
#
# BRAIN computes soil NO with BDSNP (fertiliser and deposition reservoirs, soil moisture and temperature, pulsing; Hudman
# et al. 2012, https://doi.org/10.5194/acp-12-7779-2012); this inventory with the MEGANv3.21 default (Yienger and Levy
# 1995 type, https://doi.org/10.1029/95JD00370). The table gives the share of the state total in crop-dominated cells and
# in each mesoregion, and the rank correlation of ln(ratio) with the crop fraction.

# %%
# 11 · H. Soil NO -----------------------------------------------------------------------------------------------------------------------------------
fo_no, fb_no = FLUX['NO']
crop_b = np.nan_to_num(GF['crop']) > DOMINANT
sh_o = float((fo_no * G['a_true'][:, None] * F_SC_B)[crop_b].sum() / (fo_no * G['a_true'][:, None] * F_SC_B).sum())
sh_b = float((fb_no * G['a_true'][:, None] * F_SC_B)[crop_b].sum() / (fb_no * G['a_true'][:, None] * F_SC_B).sum())
y, X, rw, cl = cell_table('NO'); rho_c = stats.spearmanr(y, X[:, list(PB).index('crop')])[0]
H_ = TR_[TR_.item == 'NO'][['mesoregion', 'brain_share_pct', 'this_share_pct', 'ratio']]
H_.to_csv(os.path.join(OUT23, 'H_soil_no_mesoregions.csv'), index=False)
print(f'soil NO in crop-dominated cells: {100 * sh_o:.1f} % of the state here, {100 * sh_b:.1f} % in BRAIN; '
      f'Spearman of ln(ratio) with the crop fraction {rho_c:.2f} (n = {len(y)})')
print(H_.round(2).to_string(index=False))

# %% [markdown]
# ## Cell 12 · Summary
#
# One line per explanation with its number; the wording for the texts is in the note that accompanies this notebook,
# conditional on these results.

# %%
# 12 · Summary ------------------------------------------------------------------------------------------------------------------------------------------
S_ROWS = [dict(block='A', item='ETOH', statistic='ratio without GAMBD', value=float(TA_.set_index('item').loc['ETOH', 'x_no_gambd'])),
          dict(block='A', item='ISOP', statistic='ratio without GAMSM', value=float(TA_.set_index('item').loc['ISOP', 'x_no_gamsm'])),
          dict(block='A', item='ISOP', statistic='ratio with PAR 0.5', value=float(TA_.set_index('item').loc['ISOP', 'x_par050'])),
          dict(block='A', item='ISOP', statistic='ratio with all three', value=float(TA_.set_index('item').loc['ISOP', 'all_three']))]
for _, r in FPS.iterrows():
    S_ROWS.append(dict(block='B', item=r.response, statistic=f'Pearson with {r.property} (permutation p)', value=r.pearson, p=r.p_perm))
for _, r in LM[LM.level == 'group'].iterrows():
    S_ROWS.append(dict(block='C', item=r['item'], statistic=f'LMG share of R2: {r["name"]}', value=r.r2_share, lo=r.get('boot_lo'), hi=r.get('boot_hi')))
for _, r in DD.iterrows():
    S_ROWS.append(dict(block='D', item=r['item'], statistic='share fixed in time (spatial)', value=r.share_spatial_fixed))
    S_ROWS.append(dict(block='D', item=r['item'], statistic='share statewide seasonal', value=r.share_seasonal_statewide))
if 'ED' in globals():
    for _, r in ED.iterrows():
        S_ROWS.append(dict(block='E', item='ISOP', statistic=f'slope per km ({r.model}); predicted {r.predicted_per_km:.3f}', value=r.observed_per_km, lo=r.lo, hi=r.hi))
if 'FV' in globals():
    for _, r in FV[(FV.level == 'state') & (FV.name == 'ISOP')].iterrows():
        S_ROWS.append(dict(block='F', item='ISOP', statistic=f'state ratio, variant {r.variant}', value=r.ratio, p=r.get('gap_closed')))
if 'GT' in globals():
    for _, r in GT.iterrows():
        S_ROWS.append(dict(block='G', item='ISOP', statistic=f'{r.test}: here minus BRAIN', value=r.this_inventory - r.brain, lo=r.get('lo'), hi=r.get('hi')))
S_ROWS.append(dict(block='H', item='NO', statistic='share in crop-dominated cells, here minus BRAIN', value=sh_o - sh_b))
SUMM = pd.DataFrame(S_ROWS); SUMM.to_csv(os.path.join(OUT23, 'summary.csv'), index=False)
print(SUMM.round(4).to_string(index=False))
json.dump(dict(notebook='23_brain_differences_2019', run_utc=datetime.datetime.utcnow().isoformat(timespec='seconds') + 'Z', nb18=NB18,
               n20=N20_OUT if HAVE_N20 else None, n21=N21_JSON if os.path.exists(N21_JSON) else None, tropomi=TROPOMI if os.path.isdir(TROPOMI) else None,
               n_boot=N_BOOT, n_perm=N_PERM, block_cells=BLOCK, seed=SEED), open(os.path.join(OUT23, 'provenance_notebook23.json'), 'w'), indent=1)
print('outputs in', OUT23)
