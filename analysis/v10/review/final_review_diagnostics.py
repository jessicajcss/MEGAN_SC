# -*- coding: utf-8 -*-
"""final_review_diagnostics.py - additional diagnostics of the Santa Catarina example (trends, LAI, TROPOMI magnitude, ...).

Every value written here is computed from files of this repository; nothing is typed in by hand except the literature
constants named below. Diagnostics defined after the first results were seen are marked POST HOC in the log.

Inputs (relative to the repository root, or to --results / --inputs):
  results/v10/notebook13_meso2024/  table2_annual_totals_SC_groups.csv, annual_totals_SC_by_class.csv,
                                    drivers_annual_SC.csv, table3_mesoregions_soa.csv, table4_municipalities_2025.csv
  results/v10/notebook14/           monthly_series_cf0.3_sza70_fire120.csv, tier1_spatial_cf0.3_sza70_fire120.csv
  results/v10/diagnostics/nb14/     nb14_diagnostics.csv (D2 state means)
  results/v10/diagnostics/c3/       c3_filled_cells_2003_2025.csv
  data/run_v10/<year>/inputs/       lai_<year>_<mm>.tif and fixlog_<year>_<mm>.npz (optional: LAI and filled-cell sections)
  data/boundaries/                  sc_state_ibge.gpkg (optional, with the inputs above)

Outputs: reproduced/final_review/final_review_numbers.csv (key, value, note) and final_review.log (published copy in
results/v10/final_review)

Literature constants:
  Palmer et al. (2003, JGR 108(D6) 4180, Table 3): mean HCHO column lifetime at 10-12 LT, North America, July 1996:
      1.48 to 1.76 h; final HCHO yields 0.29 (0.1 ppb NOx) and 0.43 per C (1 ppb NOx) from their Figure 2.
  Vigouroux et al. (2020, AMT 13, 3751-3767): TROPOMI = 1.10e15 + 0.64 FTIR (inverted for the bias-corrected case).
  MEGANv3.21 (Guenther, 2024), megvea.f: gamma_CO2 = ISmax - ISmax Ci^h / (C*^h + Ci^h), Ci = 0.7 CO2, ISmax = 1.344,
      h = 1.4614, C* = 585 (applied only if annual CO2 values are supplied with --co2, see below).

Usage:  python analysis/v10/review/final_review_diagnostics.py [--repo .] [--co2 path/to/co2_annmean_gl.txt]
The CO2 file is the NOAA GML global annual mean (https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_gl.txt); the
copy shipped in data/external/co2_annmean_gl.txt (see its header) is used by default, and without a CO2 file the
CO2-inhibition section is skipped.
Section 9 interprets the trends: the exact decomposition of the OLS trend of ln(annual emission) into the parts carried
by temperature, ln(LAI) and the residual year term of the regression of section 2, and the change of the land cover as
the model sees it (state area of each MEGAN growth form, data/run_v10/<year>/inputs/growthform_*.npz, and the state
emission potential sum(EF x cell area), inputs/ef_<year>.npz, which does not depend on meteorology or LAI).
"""
import argparse, glob, itertools, os, sys
import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument('--repo', default=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')))
ap.add_argument('--results', default=None, help='default: <repo>/results/v10')
ap.add_argument('--inputs', default=None, help='default: <repo>/data/run_v10 (per-year model inputs)')
ap.add_argument('--co2', default=None, help='NOAA GML co2_annmean_gl.txt (default: data/external/co2_annmean_gl.txt)')
ap.add_argument('--out', default=None, help='default: reproduced/final_review')
a = ap.parse_args()
if a.co2 is None:
    _c = os.path.join(a.repo, 'data', 'external', 'co2_annmean_gl.txt')
    a.co2 = _c if os.path.exists(_c) else None
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
RES = a.results or P.RESULTS
INP = a.inputs or P.RUN
OUT = a.out or P.out('final_review')
def pick(*rel):
    """The reproduced file if reproduce.py has written it, else the published one."""
    f = os.path.join(P.OUT, *rel); return f if os.path.exists(f) else os.path.join(RES, *rel)
os.makedirs(OUT, exist_ok=True)
TAG = 'cf0.3_sza70_fire120'
N13 = P.NB13; N14 = os.path.dirname(pick('notebook14', f'monthly_series_{TAG}.csv'))
ROWS, LOG = [], []


def put(key, value, note=''):
    ROWS.append(dict(key=key, value=value, note=note)); LOG.append(f'{key}: {value}' + (f'   [{note}]' if note else ''))


def section(t):
    LOG.append(''); LOG.append('== ' + t)


# ------------------------------------------------------------------ 1. trends: Hamed-Rao, false-discovery rate, start year
section('1. Trends (Mann-Kendall), Hamed-Rao, Benjamini-Hochberg, start year, 2025 excluded')
import pymannkendall as pmk
T = pd.read_csv(os.path.join(N13, 'table2_annual_totals_SC_groups.csv'))
GROUPS = ['BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'MBO', 'NO', 'CO']
rows = []
for g in GROUPS:
    y = T[g].values.astype(float); n = y.size
    o = pmk.original_test(y); hr = pmk.hamed_rao_modification_test(y)
    res = y - o.slope * np.arange(n); rk = stats.rankdata(res); x = rk - rk.mean()
    acf = np.correlate(x, x, 'full')[n - 1:] / np.sum(x * x)
    o05 = pmk.original_test(y[T.year.values >= 2005]); o24 = pmk.original_test(y[T.year.values <= 2024])
    rows.append(dict(group=g, MK_p=o.p, HR_p=hr.p, max_abs_rank_acf=float(np.max(np.abs(acf[1:]))), lag1_rank_acf=float(acf[1]),
                     p_from_2005=o05.p, pct_yr_2003_2024=100 * o24.slope / y[T.year.values <= 2024].mean(), p_2003_2024=o24.p))
TR = pd.DataFrame(rows)
# Benjamini-Hochberg step-up adjustment over the nine groups
p = TR.MK_p.values; order = np.argsort(p); m = p.size
q = np.empty(m); q[order] = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
TR['BH_q'] = np.minimum(q, 1)
TR.to_csv(os.path.join(OUT, 'trends_fdr_hr.csv'), index=False)
LOG.append(TR.round(4).to_string(index=False))
put('hr_threshold', round(1.96 / np.sqrt(len(T)), 3), 'pymannkendall keeps lags with |rank acf| > 1.96/sqrt(n)')
put('hr_max_abs_acf', round(TR.max_abs_rank_acf.max(), 3), 'largest |rank acf| of the detrended series, any lag, any group')
put('hr_changes_p', bool((TR.HR_p != TR.MK_p).any()))
# soil NO with its lag-1 autocorrelation forced into the Hamed-Rao variance (to show the size of the effect)
y = T['NO'].values.astype(float); n = y.size; o = pmk.original_test(y); rk = stats.rankdata(y - o.slope * np.arange(n))
x = rk - rk.mean(); r1 = float((np.correlate(x, x, 'full')[n - 1:] / np.sum(x * x))[1])
ratio = 1 + 2 / (n * (n - 1) * (n - 2)) * (n - 1) * (n - 2) * (n - 3) * r1
z = (o.s - np.sign(o.s)) / np.sqrt(n * (n - 1) * (2 * n + 5) / 18 * ratio)
put('NO_p_with_lag1', round(2 * stats.norm.sf(abs(z)), 4), f'lag-1 rank acf {r1:.3f}; n/n* = {ratio:.3f} (POST HOC)')
for g in ('ISOP', 'NO', 'BVOC', 'MT', 'SQT'):
    r = TR.set_index('group').loc[g]
    put(f'trend_{g}', f'MK p {r.MK_p:.4f}; BH q {r.BH_q:.4f}; p from 2005 {r.p_from_2005:.4f}; 2003-2024 {r.pct_yr_2003_2024:+.2f} %/yr (p {r.p_2003_2024:.4f})')

# ------------------------------------------------------------------ 2. drivers: LAI 2024-2025, partial correlations, residual year term
section('2. Drivers: LAI in 2024-2025, partial correlations (detrended), residual year term (POST HOC)')
D = pd.read_csv(os.path.join(N13, 'drivers_annual_SC.csv')).set_index('year')
base = D.loc[2003:2023, 'LAI']
for yy in (2024, 2025):
    v = D.loc[yy, 'LAI']
    put(f'lai_{yy}', f'{v:.3f} ({100 * (v / base.mean() - 1):+.1f} %; {(v - base.mean()) / base.std(ddof=1):+.1f} SD of 2003-2023)',
        f'2003-2023 mean {base.mean():.3f} +- {base.std(ddof=1):.3f}')
_dt = lambda v: v - np.polyval(np.polyfit(np.arange(len(v)), v, 1), np.arange(len(v)))
put('driver_interannual', f"PAR CV {100 * D.PAR.std(ddof=1) / D.PAR.mean():.1f} %; T SD {D.TEMP2_C.std(ddof=1):.2f} K; "
    f"detrended r(PAR, SOILM) {np.corrcoef(_dt(D.PAR.values), _dt(D.SOILM.values))[0, 1]:.2f}; "
    f"detrended r(T, SOILM) {np.corrcoef(_dt(D.TEMP2_C.values), _dt(D.SOILM.values))[0, 1]:.2f}",
    'annual state means 2003-2025 of drivers_annual_SC.csv (SD with ddof = 1)')
A = pd.read_csv(os.path.join(N13, 'annual_totals_SC_by_class.csv')).merge(D.reset_index(), on='year')
def det(v):
    t = np.arange(len(v)); return v - np.polyval(np.polyfit(t, v, 1), t)
def pcor(a_, b_, c_):
    a_, b_, c_ = det(a_), det(b_), det(c_); rab = np.corrcoef(a_, b_)[0, 1]; rac = np.corrcoef(a_, c_)[0, 1]; rbc = np.corrcoef(b_, c_)[0, 1]
    return (rab - rac * rbc) / np.sqrt((1 - rac ** 2) * (1 - rbc ** 2))
BVOC_CLASSES = [c for c in A.columns if c not in ('year', 'NO', 'TEMP2_C', 'PAR', 'SOILM', 'LAI')]
PC = pd.DataFrame([dict(cls=c, r_LAI_given_T=pcor(A[c].values, A.LAI.values, A.TEMP2_C.values),
                        r_T_given_LAI=pcor(A[c].values, A.TEMP2_C.values, A.LAI.values)) for c in BVOC_CLASSES + ['NO']])
PC.to_csv(os.path.join(OUT, 'partial_correlations_detrended.csv'), index=False)
LOG.append(PC.round(3).to_string(index=False))
b = PC[PC.cls != 'NO']
put('partial_LAI_given_T_BVOC', f'{b.r_LAI_given_T.min():.2f} to {b.r_LAI_given_T.max():.2f}',
    'top: ' + ', '.join(f'{r.cls} {r.r_LAI_given_T:.2f}' for r in b.sort_values('r_LAI_given_T', ascending=False).head(7).itertuples()))
put('r_detrended_T_LAI', round(np.corrcoef(det(A.TEMP2_C.values), det(A.LAI.values))[0, 1], 2))
import statsmodels.api as sm
for g, col in (('ISOP', 'ISOP'), ('NO', 'NO')):
    X = sm.add_constant(pd.DataFrame({'T': A.TEMP2_C, 'lnLAI': np.log(A.LAI), 'yr': A.year - 2003}))
    f = sm.OLS(np.log(A[col]), X).fit()
    put(f'residual_year_term_{g}', f'{100 * f.params.yr:+.2f} +- {100 * f.bse.yr:.2f} %/yr',
        f'OLS of ln(annual emission) on T, ln(LAI), year (n = 23); coefficient of T {f.params["T"]:.3f} K-1')
put('T_trend', f"{pmk.original_test(A.TEMP2_C.values).slope:.4f} K/yr (MK p {pmk.original_test(A.TEMP2_C.values).p:.3f})")

# ------------------------------------------------------------------ 3. TROPOMI: background control, seasons, annual, magnitude
section('3. TROPOMI: Atlantic-background control, seasons, annual means, magnitude (POST HOC except the annual mean)')
S = pd.read_csv(os.path.join(N14, f'monthly_series_{TAG}.csv'))
an = S.copy()
for c in ('ISOP_t', 'dHCHO_SC', 'T2_C', 'BG', 'HCHO_SC'):
    an[c] = S[c] - S.groupby('month')[c].transform('mean')
def neff_r(a_, b_, k=0):
    a_, b_ = np.asarray(a_, float), np.asarray(b_, float)
    r = np.corrcoef(a_, b_)[0, 1]; r1 = np.corrcoef(a_[:-1], a_[1:])[0, 1]; r2 = np.corrcoef(b_[:-1], b_[1:])[0, 1]
    ne = len(a_) * (1 - r1 * r2) / (1 + r1 * r2) - k; t = r * np.sqrt((ne - 2) / (1 - r ** 2))
    return r, ne, 2 * stats.t.sf(abs(t), ne - 2)
def resid(yv, xv):
    return yv - np.polyval(np.polyfit(xv, yv, 1), xv)
r, ne, pp = neff_r(an.ISOP_t, an.BG); put('anom_r_isop_background', f'{r:.3f} (N_eff {ne:.0f}, p {pp:.1e})')
r, ne, pp = neff_r(resid(an.ISOP_t.values, an.BG.values), resid(an.dHCHO_SC.values, an.BG.values), 1)
put('anom_partial_r_isop_dHCHO_given_background', f'{r:.3f} (N_eff {ne:.0f}, p {pp:.1e})')
r, ne, pp = neff_r(resid(an.ISOP_t.values, an.T2_C.values), resid(an.dHCHO_SC.values, an.T2_C.values), 1)
put('anom_partial_r_isop_dHCHO_given_T2', f'{r:.3f} (N_eff {ne:.0f}, p {pp:.1e})')
put('anom_r_dHCHO_background', round(np.corrcoef(an.dHCHO_SC, an.BG)[0, 1], 3))
r1 = np.corrcoef(an.ISOP_t[:-1], an.ISOP_t[1:])[0, 1]; r2 = np.corrcoef(an.dHCHO_SC[:-1], an.dHCHO_SC[1:])[0, 1]
r = np.corrcoef(an.ISOP_t, an.dHCHO_SC)[0, 1]; ne12 = (len(an) - 12) * (1 - r1 * r2) / (1 + r1 * r2)
put('anom_r_isop_dHCHO_neff_minus12', f'N_eff {ne12:.1f}, p {2 * stats.t.sf(r * np.sqrt((ne12 - 2) / (1 - r ** 2)), ne12 - 2):.1e}',
    'the 12 monthly means removed by deseasonalising taken out of N')
SEAS = {'DJF': (12, 1, 2), 'MAM': (3, 4, 5), 'JJA': (6, 7, 8), 'SON': (9, 10, 11)}
for k, mm in SEAS.items():
    s = an[an.month.isin(mm)]
    put(f'anom_{k}', f'n {len(s)}; r(isop, dHCHO) {np.corrcoef(s.ISOP_t, s.dHCHO_SC)[0, 1]:.2f}; r(isop, background) '
        f'{np.corrcoef(s.ISOP_t, s.BG)[0, 1]:.2f}')
ann = S.groupby('year')[['ISOP_t', 'dHCHO_SC', 'BG']].mean()
put('annual_r_isop_dHCHO', round(np.corrcoef(ann.ISOP_t, ann.dHCHO_SC)[0, 1], 3), 'n = 7; written by notebook 14 as descriptive')
put('annual_r_isop_background', round(np.corrcoef(ann.ISOP_t, ann.BG)[0, 1], 3), 'n = 7')
clim = S.groupby('month')[['ISOP_t', 'dHCHO_SC']].mean()
put('clim_min_month_isop_2019_2025', int(clim.ISOP_t.idxmin())); put('clim_min_month_dHCHO_2019_2025', int(clim.dHCHO_SC.idxmin()))
T1 = pd.read_csv(os.path.join(N14, f'tier1_spatial_{TAG}.csv'))
sm_ = T1[T1.smoothing_km > 0]
put('tier1_r_smoothed_10_40km', f'{sm_.pearson_r.min():.2f} to {sm_.pearson_r.max():.2f}')
DG = pd.read_csv(pick('diagnostics', 'nb14', 'nb14_diagnostics.csv'), index_col=0)['value']
E = float([v for k, v in DG.items() if k.startswith('D2 SC-mean model isoprene')][0])
dH = float([v for k, v in DG.items() if k.startswith('D2 SC-mean TROPOMI dHCHO')][0])
_bc = [v for k, v in DG.items() if k.startswith('D2 with BIAS_CORRECT')]
dHbc = float(str(_bc[0]).split(';')[0]) if _bc else dH / 0.64   # Vigouroux et al. (2020) inverted (offset cancels in a difference)
GRID = []
for tau in (1.0, 1.48, 1.5, 1.76, 2.0, 3.0, 4.0):
    for Y in (0.29, 0.43):
        GRID.append(dict(tau_h=tau, yield_perC=Y, ratio_obs_expected_raw=dH / (Y * E * tau * 3600), ratio_obs_expected_biascorr=dHbc / (Y * E * tau * 3600)))
G = pd.DataFrame(GRID); G.to_csv(os.path.join(OUT, 'tropomi_magnitude_grid.csv'), index=False); LOG.append(G.round(3).to_string(index=False))
for Y in (0.29, 0.43):
    put(f'implied_tau_h_Y{Y}', f'raw {dH / (Y * E * 3600):.3f}; bias-corrected {dHbc / (Y * E * 3600):.3f}')
g = G[G.tau_h.isin((1.48, 1.76))]
put('ratio_obs_expected_palmer_tau', f'raw {g.ratio_obs_expected_raw.min():.3f} to {g.ratio_obs_expected_raw.max():.3f}; '
    f'bias-corrected {g.ratio_obs_expected_biascorr.min():.3f} to {g.ratio_obs_expected_biascorr.max():.3f}',
    'tau 1.48-1.76 h (Palmer et al., 2003, Table 3), Y 0.29-0.43 per C')
lon_w, lon_e = -53.84, -48.33
put('overpass_local_solar_time', f'16:30 UTC = {16.5 + lon_w / 15:.2f} h to {16.5 + lon_e / 15:.2f} h local mean solar time across SC')

# ------------------------------------------------------------------ 4. mesoregions and municipalities
section('4. Mesoregions (Table 3) and municipalities (Table 4)')
M3 = pd.read_csv(os.path.join(N13, 'table3_mesoregions_soa.csv'))
xf, yi = M3.forest_pct.values, M3.intensity_kg_ha.values
r, pr = stats.pearsonr(xf, yi); zf = np.arctanh(r); se = 1 / np.sqrt(len(xf) - 3)
put('meso_r_fisher', f'r {r:.4f} (95% CI {np.tanh(zf - 1.96 * se):.4f} to {np.tanh(zf + 1.96 * se):.4f}; p {pr:.4f})')
rho = stats.spearmanr(xf, yi)[0]; rx, ry = stats.rankdata(xf), stats.rankdata(yi)
cnt = sum(abs(np.corrcoef(rx, pm)[0, 1]) >= abs(rho) - 1e-12 for pm in itertools.permutations(ry))
put('meso_spearman_exact', f'rho {rho:.4f}; exact two-sided p {cnt / 720:.4f} (asymptotic {stats.spearmanr(xf, yi)[1]:.4f})')
M4 = pd.read_csv(os.path.join(N13, 'table4_municipalities_2025.csv')); M4['ratio'] = M4.area_grid_km2 / M4.area_ibge_km2
put('mun_area_ratio', f'{M4.ratio.min():.4f} to {M4.ratio.max():.4f} (median {M4.ratio.median():.4f})',
    'area attributed by area weighting (sum of the cell fractions times the cell areas on the authalic sphere) / IBGE area')
M3a = pd.read_csv(os.path.join(N13, 'table3_mesoregions_soa.csv'))
put('meso_area_ratio', f'{(M3a.area_grid_ha / M3a.area_ibge_ha).min():.4f} to {(M3a.area_grid_ha / M3a.area_ibge_ha).max():.4f}',
    'mesoregions: area attributed by area weighting / sum of the IBGE areas of the municipalities')

# ------------------------------------------------------------------ 5. check C3 pooled over the 23 years
section('5. Check C3 pooled over 2003-2025 (descriptive: the same cells recur every year)')
C = pd.read_csv(pick('diagnostics', 'c3', 'c3_filled_cells_2003_2025.csv'))
for k, gg in [('all', C)] + list(C.groupby('cls')):
    nn = len(gg); kk = int(gg.exceed.sum())
    put(f'c3_pooled_{k}', f'{kk} of {nn} above the valid-cell p99.5 (expected {0.005 * nn:.1f}; binomial p {stats.binom.sf(kk - 1, nn, 0.005):.1e})')
put('c3_distinct_cells', int(C.groupby(['j', 'i']).ngroups))

# ------------------------------------------------------------------ 6. cell areas: authalic sphere vs WGS 84 ellipsoid
section('6. Cell areas: authalic sphere (megan_sc.domain) vs WGS 84 ellipsoid')
try:
    from pyproj import Geod
    geod = Geod(ellps='WGS84'); RA = 6371007.1810; d = 0.025; out = []
    for lat0 in (-29.5, -27.6, -25.9):
        ae = abs(geod.polygon_area_perimeter([0, d, d, 0], [lat0, lat0, lat0 + d, lat0 + d])[0])
        asph = RA ** 2 * np.deg2rad(d) * abs(np.sin(np.deg2rad(lat0 + d)) - np.sin(np.deg2rad(lat0)))
        out.append(100 * (asph / ae - 1))
    put('cell_area_sphere_minus_ellipsoid_pct', f'{min(out):.3f} to {max(out):.3f} % at 29.5-25.9 S')
except ImportError:
    LOG.append('pyproj not installed: section 6 skipped')

# ------------------------------------------------------------------ 7. model inputs: monthly LAI and filled coastal cells
section('7. Model inputs: monthly state-mean LAI (MODIS, as used) and filled coastal cells')
if os.path.isdir(INP):
    import geopandas as gpd, rasterio
    from rasterio.features import geometry_mask, rasterize
    from rasterio.transform import from_origin
    sc = gpd.read_file(os.path.join(P.BND, 'sc_state_ibge.gpkg')).to_crs(4326)
    geom = sc.geometry.union_all() if hasattr(sc.geometry, 'union_all') else sc.geometry.unary_union
    rows, mask = [], None
    for yy in range(2003, 2026):
        for mm in range(1, 13):
            f = os.path.join(INP, str(yy), 'inputs', f'lai_{yy}_{mm:02d}.tif')
            if not os.path.exists(f):
                continue
            with rasterio.open(f) as s:
                arr = s.read(1).astype(float)
                if mask is None:
                    mask = ~geometry_mask([geom], arr.shape, s.transform, all_touched=False)
            v = arr[mask]; rows.append(dict(year=yy, month=mm, lai_mean_vegetated=v[v > 0].mean(), n_cells=int((v > 0).sum())))
    if rows:
        L = pd.DataFrame(rows); L.to_csv(os.path.join(OUT, 'lai_monthly_SC.csv'), index=False)
        PV = L.pivot(index='year', columns='month', values='lai_mean_vegetated'); cl = PV.loc[2003:2023].mean(); sd = PV.loc[2003:2023].std()
        Z = (PV - cl) / sd; Z.round(2).to_csv(os.path.join(OUT, 'lai_monthly_SC_zscores.csv'))
        if 2025 in Z.index and 2024 in Z.index:
            seq = list(Z.loc[2024, [11, 12]]) + list(Z.loc[2025])
            put('lai_months_below_1sd_nov2024_dec2025', f'{sum(v < -1 for v in seq)} of {len(seq)}',
                'state-mean LAI of vegetated cells vs the same month in 2003-2023')
            put('lai_months_below_2sd_by_year', '; '.join(f'{yy} {int((Z.loc[yy] < -2).sum())}' for yy in Z.index))
        put('lai_max_cell_value', f'{max(rasterio.open(f).read(1).max() for f in sorted(glob.glob(os.path.join(INP, "*", "inputs", "lai_*.tif")))[::12]):.2f}',
            'MODIS fill values (24.9-25.5 after scaling) would show here')
    # filled cells (FIX-1): cells with any filled meteorology step in the year, inside SC
    fl = sorted(glob.glob(os.path.join(INP, '2019', 'inputs', 'fixlog_2019_*.npz')))
    if fl:
        from netCDF4 import Dataset
        with Dataset(glob.glob(os.path.join(INP, '2019', 'emissions', 'annual_emis_2019.nc'))[0]) as ds:
            vn = list(ds.variables); lat = np.array(ds.variables[[v for v in vn if v.lower().startswith('lat')][0]][:])
            lon = np.array(ds.variables[[v for v in vn if v.lower().startswith('lon')][0]][:])
        NY, NX = lat.size, lon.size; dd = float(lon[1] - lon[0]); TRF = from_origin(lon[0] - dd / 2, lat[-1] + dd / 2, dd, dd)
        SCm = np.flipud(rasterize([(geom, 1)], out_shape=(NY, NX), transform=TRF, fill=0, dtype='int32')) == 1
        AREA = (6371008.8 ** 2 * np.deg2rad(dd) * np.abs(np.sin(np.deg2rad(lat + dd / 2)) - np.sin(np.deg2rad(lat - dd / 2))))[:, None] * np.ones((1, NX))
        out = []
        for yy in range(2003, 2026):
            ff = sorted(glob.glob(os.path.join(INP, str(yy), 'inputs', f'fixlog_{yy}_*.npz')))
            if not ff:
                continue
            s_ = sum(np.load(f)['filled_met_steps'] for f in ff); m_ = (s_ > 0) & SCm
            out.append(dict(year=yy, filled_cells=int(m_.sum()), area_km2=float((AREA * m_).sum()) / 1e6))
        FC = pd.DataFrame(out); FC.to_csv(os.path.join(OUT, 'filled_cells_area.csv'), index=False)
        put('filled_cells_area', '; '.join(f'{int(n)} cells = {ar:.1f} km2' for n, ar in FC.groupby('filled_cells').area_km2.first().items()))
else:
    LOG.append(f'{INP} not found: sections on monthly LAI and filled cells skipped')

# ------------------------------------------------------------------ 8. CO2 inhibition of isoprene (needs NOAA annual means)
section('8. CO2 inhibition of isoprene applied to the annual totals (gamma_CO2 is spatially uniform, so this is exact)')
if a.co2 and os.path.exists(a.co2):
    co2 = pd.read_csv(a.co2, comment='#', sep=r'\s+', header=None, names=['year', 'co2', 'unc']).set_index('year').co2
    def gco2(c):
        ci = 0.7 * c; h = 1.4614; return 1.344 - 1.344 * ci ** h / (585 ** h + ci ** h)
    yrs = T.year.values; g_ = np.array([gco2(co2[yy]) for yy in yrs]); iso = T.ISOP.values * g_
    o = pmk.original_test(iso); o0 = pmk.original_test(T.ISOP.values)
    put('isop_trend_with_co2_inhibition', f'{100 * o.slope / iso.mean():+.2f} %/yr (MK p {o.p:.4f}); without: {100 * o0.slope / T.ISOP.mean():+.2f} %/yr',
        f'CO2 {co2[yrs[0]]:.1f} ppm ({yrs[0]}) to {co2[yrs[-1]]:.1f} ppm ({yrs[-1]}); gamma_CO2 {g_[0]:.3f} to {g_[-1]:.3f}')
    bv = T.BVOC.values - T.ISOP.values + iso; ob = pmk.original_test(bv)
    put('bvoc_trend_with_co2_inhibition', f'{100 * ob.slope / bv.mean():+.2f} %/yr (MK p {ob.p:.4f}); without: '
        f'{100 * pmk.original_test(T.BVOC.values).slope / T.BVOC.mean():+.2f} %/yr', 'isoprene of every year multiplied by gamma_CO2')
    put('isop_mean_with_co2_inhibition', f'{100 * (iso.mean() / T.ISOP.mean() - 1):+.2f} % of the 2003-2025 mean isoprene')
    # Benjamini-Hochberg over the nine groups with the isoprene and BVOC series carrying the inhibition
    p2 = TR.set_index('group').MK_p.copy(); p2['ISOP'] = o.p; p2['BVOC'] = ob.p
    pv = p2.values; od = np.argsort(pv); mm = pv.size
    q2 = np.empty(mm); q2[od] = np.minimum.accumulate((pv[od] * mm / np.arange(1, mm + 1))[::-1])[::-1]
    q2 = pd.Series(np.minimum(q2, 1), index=p2.index)
    put('fdr_with_co2_inhibition', f"BH q: isoprene {q2['ISOP']:.4f}, BVOC {q2['BVOC']:.4f}, soil NO {q2['NO']:.4f}",
        'MK p of the nine groups, with the isoprene and BVOC series multiplied by gamma_CO2')
else:
    LOG.append('no CO2 file given (--co2): section skipped')

# ------------------------------------------------------------------ 9. Interpretation of the trends (POST HOC)
section('9. Trends: decomposition of the OLS log trend (temperature, LAI, year term) and land cover as seen by the model')
for g in ('ISOP', 'BVOC', 'MT', 'SQT', 'NO'):
    X = sm.add_constant(pd.DataFrame({'T': A.TEMP2_C, 'lnLAI': np.log(A.LAI), 'yr': A.year - 2003}))
    f = sm.OLS(np.log(A[g] if g in A else T.set_index('year').loc[A.year, g].values), X).fit()
    sT = np.polyfit(A.year, A.TEMP2_C, 1)[0]; sL = np.polyfit(A.year, np.log(A.LAI), 1)[0]
    tot = np.polyfit(A.year, f.model.endog, 1)[0]              # OLS slope of ln(emission) on year = sum of the three parts
    put(f'ols_trend_decomposition_{g}', f'{100 * tot:+.2f} %/yr = T {100 * f.params["T"] * sT:+.2f} + LAI {100 * f.params["lnLAI"] * sL:+.2f} '
        f'+ year term {100 * f.params["yr"]:+.2f}', f'OLS slope of ln(annual emission) on year; coefficients of the regression of section 2 '
        f'(T {f.params["T"]:.3f} K-1, ln LAI {f.params["lnLAI"]:.3f}); OLS trends of T {sT:.4f} K/yr and ln LAI {100 * sL:+.3f} %/yr')
if os.path.isdir(INP):
    sys.path.insert(0, os.path.join(a.repo, 'analysis', 'v10')); import common as cm_
    la_, lo_ = cm_.grid(); SCm = cm_.sc_mask(la_, lo_); Ac = cm_.cell_area_m2(la_, lo_)
    rows9 = []
    for yy in range(2003, 2026):
        fi = os.path.join(INP, str(yy), 'inputs'); gf = sorted(x for x in os.listdir(fi) if x.startswith('growthform_') and x.endswith('.npz'))[0]
        zg = np.load(os.path.join(fi, gf)); ze = np.load(os.path.join(fi, f'ef_{yy}.npz')); ef = ze['ef']; cl = list(ze['classes'])
        r9 = dict(year=yy, growthform_file=gf, **{f'{k}_km2': float((zg[k] * Ac * SCm).sum() / 1e6) for k in ('tree', 'shrub', 'herb', 'crop', 'veg')})
        for nm, ks in (('ISOP', ['ISOP']), ('MT', [k for k in cl if k.startswith('MT_')]), ('SQT', ['SQT_HR', 'SQT_LR']), ('NO', ['NO']),
                       ('BVOC', [k for k in cl if k not in ('NO', 'CO')])):
            r9[f'EFxA_{nm}_mol_s'] = float(sum((np.nan_to_num(ef[cl.index(k)]) * Ac * SCm).sum() for k in ks) * 1e-9)
        rows9.append(r9)
    G9 = pd.DataFrame(rows9); G9.to_csv(os.path.join(OUT, 'landcover_ef_potential_SC.csv'), index=False)
    # LAI inputs: cell-months filled from the neighbouring months (FIX-4, fix logs) and the largest cell value of all months
    gaps = [(yy, mm, int(np.load(os.path.join(INP, str(yy), 'inputs', f'fixlog_{yy}_{mm:02d}.npz'))['lai_gap_cells']))
            for yy in range(2003, 2026) for mm in range(1, 13)]
    put('lai_gap_cell_months', f"{sum(g for *_, g in gaps)} in 2003-2025: " + '; '.join(f'{yy}-{mm:02d} {g}' for yy, mm, g in gaps if g),
        'vegetated grid cells without a valid MODIS pixel in the month, filled from the adjacent months (fixlog lai_gap_cells)')
    import rasterio
    put('lai_max_cell_all_months', f"{max(float(rasterio.open(f).read(1).max()) for f in glob.glob(os.path.join(INP, '*', 'inputs', 'lai_*.tif'))):.2f}",
        'largest monthly cell LAI of the 276 months (MODIS fill values, 24.9-25.5 after scaling, would show here)')
    y24 = G9[G9.year <= 2024]
    for k in ('tree', 'herb', 'crop', 'veg'):
        v = y24[f'{k}_km2'].values; o = pmk.original_test(v)
        put(f'growthform_area_{k}', f'{v[0]:,.0f} km2 (2003) to {v[-1]:,.0f} km2 (2024), {100 * (v[-1] / v[0] - 1):+.1f} %; '
            f'Theil-Sen {100 * o.slope / v.mean():+.2f} %/yr (MK p {o.p:.1e})', 'state area of the MEGAN growth form (MapBiomas 2003-2024; 2025 uses 2024)')
    for nm in ('ISOP', 'MT', 'SQT', 'BVOC', 'NO'):
        v = G9[f'EFxA_{nm}_mol_s'].values
        put(f'ef_potential_{nm}', f'{100 * (v[G9.year == 2024][0] / v[0] - 1):+.1f} % from 2003 to 2024; OLS trend of ln over 2003-2025 '
            f'{100 * np.polyfit(G9.year, np.log(v), 1)[0]:+.2f} %/yr', 'state sum of EF x cell area: change of the emission potential carried by land cover alone')
else:
    LOG.append(f'{INP} not found: section 9 (land cover) skipped')

pd.DataFrame(ROWS).to_csv(os.path.join(OUT, 'final_review_numbers.csv'), index=False)
open(os.path.join(OUT, 'final_review.log'), 'w', encoding='utf-8').write('\n'.join(LOG) + '\n')
print('\n'.join(LOG))
