#!/usr/bin/env python3
"""10_report_numbers_v10.py
The numbers of the Santa Catarina example (statistics by class, mesoregion and municipality) computed from the production run (data/run_v10, fixes.py sha
a933b29602929a82). Estimators:
  - standardised regression: annual series in z-scores (np.std, ddof = 0), no intercept, 4 drivers (TEMP2_C, PAR, LAI, SOILM);
  - temperature sensitivity: slope b of ln(E) on annual mean temperature, 100 (e^b - 1) % per K (raw and detrended);
  - correlations: Pearson r on linearly detrended annual series (from notebook 13, correlations_emissions_drivers.csv);
  - mesoregion and municipal sums: area weighting, i.e. the fraction of each cell inside each polygon (analysis/v10/common.py
    area_weights, as in notebook 13); intensities per IBGE area;
  - isoprene: times the canopy loss and production factor rho = 0.96 (Guenther et al., 2006; common.RHO), as in notebook 13.
Inputs: the outputs of notebook 13 (MEGAN_SC_NB13: tables and annual_grids_cache.npz) and the production run (data/run_v10:
        <year>/stats/stats_<year>_<mm>.csv, <year>/inputs/*). Outputs: reproduced/report/report_numbers_v10.csv and the CSVs below.
Usage: python3 10_report_numbers_v10.py
"""
import os, sys, glob, numpy as np, pandas as pd, rasterio
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
# paths from ../_paths.py (MEGAN_SC_* environment variables); common.py (grid, state mask, cell areas, polygon
# rasterisation) of analysis/v10
sys.path.insert(0, os.path.join(HERE, '..')); import common as c
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
V = P.RUN; A10 = P.NB13; P.grid_cache()
YEARS = list(range(2003, 2026)); OUT = P.out('report')
LOG = []
def rec(key, value, src):
    LOG.append(dict(key=key, value=value, source=src))
def pt(x, d=1):
    s = f"{x:,.{d}f}"
    if s.startswith('-') and not s.strip('-0.,'):     # a value that rounds to zero is printed without a sign
        s = s[1:]
    return s.replace(',', '§').replace('.', ',').replace('§', '.')
lat, lon = c.grid(); SC = c.sc_mask(lat, lon); AREA = c.cell_area_m2(lat, lon)
assert SC.sum() == 13970

# ---------------------------------------------------------------- 3.2 means, 2025, CV, extremes (Table 2 groups)
G = pd.read_csv(f"{A10}/table2_annual_totals_SC_groups.csv").set_index('year').loc[YEARS]
G['OTHERS'] = G.LVOC + G.STRESS_OTHER + G.MBO                     # group 'HC leves e demais compostos' (light HC and other compounds)
rec('mean_BVOC_CO_kt', pt((G.BVOC + G.CO).mean() / 1e3), 'table2')
rec('mean_BVOC_kt', pt(G.BVOC.mean() / 1e3), 'table2'); rec('mean_CO_kt', pt(G.CO.mean() / 1e3), 'table2'); rec('mean_NO_kt', pt(G.NO.mean() / 1e3, 2), 'table2')
for k in ('BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'OTHERS', 'NO', 'CO'):
    rec(f'y2025_{k}_kt', pt(G.loc[2025, k] / 1e3), 'table2 (2025)')
rec('check_2025_sum', pt((G.loc[2025, ['ISOP', 'MT', 'SQT', 'OVOC', 'OTHERS']].sum() - G.loc[2025, 'BVOC']) / 1e3, 3), 'groups add up to BVOC (kt)')
cv = G.BVOC.std(ddof=1) / G.BVOC.mean() * 100
rec('cv_BVOC_pct', pt(cv), 'table2, ddof = 1'); rec('peak', f"{G.BVOC.idxmax()} ({pt(G.BVOC.max() / 1e3, 0)} kt)", 'table2'); rec('trough', f"{G.BVOC.idxmin()} ({pt(G.BVOC.min() / 1e3, 0)} kt)", 'table2')
rec('isop_peak_year', str(G.ISOP.idxmax()), 'table2')
# ---------------------------------------------------------------- trends (notebook 13: Sen slope %/yr, 95% CI, Hamed-Rao p)
T = pd.read_csv(f"{A10}/trends.csv").set_index('series')
for k in ('ISOP', 'BVOC', 'MT', 'SQT', 'OVOC', 'LVOC', 'MBO', 'NO', 'CO'):
    r = T.loc[k]; rec(f'trend_{k}', f"{pt(r.pct_yr, 2)}% ao ano ({pt(r.pct_lo95, 2)} a {pt(r.pct_hi95, 2)}; p = {pt(r.HR_p, 3)})", 'trends.csv')
# ---------------------------------------------------------------- composition (2003-2025 means)
for k in ('ISOP', 'MT', 'OVOC', 'LVOC', 'SQT', 'STRESS_OTHER', 'MBO'):
    rec(f'share_{k}', pt(G[k].mean() / G.BVOC.mean() * 100), 'table2 means')
rec('isop_mt_ratio', pt(G.ISOP.mean() / G.MT.mean(), 2), 'table2 means')
# ---------------------------------------------------------------- 3.3 isoprene monthly climatology (stats files, SC totals)
rows = []
for y in YEARS:
    for m in range(1, 13):
        s = pd.read_csv(f"{V}/{y}/stats/stats_{y}_{m:02d}.csv"); s = s[s.species == 'ISOP']
        assert len(s) == 1 and int(s.cells.iloc[0]) == 13970
        rows.append(dict(year=y, month=m, ISOP_t=float(s.total_tonnes.iloc[0]) * c.rho('ISOP')))
MI = pd.DataFrame(rows); MI.to_csv(f"{OUT}/isop_monthly_v10.csv", index=False)
ann = MI.groupby('year').ISOP_t.sum(); assert np.allclose(ann.values, G.ISOP.values, rtol=1e-4), 'monthly stats do not add up to the annual totals'
for lab, yy in (('2015_2025', range(2015, 2026)), ('2003_2025', YEARS)):
    cl = MI[MI.year.isin(yy)].groupby('month').ISOP_t.mean()
    rec(f'isop_clim_{lab}', f"jan {pt(cl[1] / 1e3)} kt; jun {pt(cl[6] / 1e3)} kt; jan/jun {pt(cl[1] / cl[6])}; max mês {cl.idxmax()} ({pt(cl.max() / 1e3)}); min mês {cl.idxmin()} ({pt(cl.min() / 1e3)})", 'stats_<year>_<mm>.csv')
# ---------------------------------------------------------------- 3.4 / 3.5 mesoregions (2003-2025 means, Table 3 of notebook 13)
NB13 = P.NB13   # notebook 13 run with sc_mesoregions_2024.gpkg (islands included)
M3 = pd.read_csv(os.path.join(NB13, 'table3_mesoregions_soa.csv')).set_index('mesoregion')
rec('meso_top_isop_mass', ' > '.join(f"{n} {pt(v / 1e3)}" for n, v in M3.ISOP.sort_values(ascending=False).head(3).items()), 'table3 (kt)')
rec('meso_top_sqt_mass', ' > '.join(f"{n} {pt(v / 1e3)}" for n, v in M3.SQT.sort_values(ascending=False).head(3).items()), 'table3 (kt)')
rec('meso_top_no_mass', ' > '.join(f"{n} {pt(v / 1e3)}" for n, v in M3.NO.sort_values(ascending=False).head(3).items()), 'table3 (kt)')
rec('soa_total_kt', pt(M3.SOA.sum() / 1e3, 0) + ' (soma das mesorregiões); estado ' + pt((G.MT + G.SQT).mean() / 1e3, 0), 'table3 / table2')
rec('soa_share_top3', '; '.join(f"{n} {pt(v)}%" for n, v in M3.pct_SC.sort_values(ascending=False).head(3).items()), 'table3')
rec('soa_intensity_top3', '; '.join(f"{n} {pt(v)}" for n, v in M3.intensity_kg_ha.sort_values(ascending=False).head(3).items()), 'table3 (kg/ha/yr)')
rec('soa_intensity_min', '; '.join(f"{n} {pt(v)}" for n, v in M3.intensity_kg_ha.sort_values().head(2).items()), 'table3 (kg/ha/yr)')
rec('sqt_share_range', f"{pt(M3.SQT_share_pct.min())}% a {pt(M3.SQT_share_pct.max())}%", 'table3')
r_, p_ = stats.pearsonr(M3.forest_pct, M3.intensity_kg_ha); z = np.arctanh(r_); se = 1 / np.sqrt(len(M3) - 3)
rho, prho = stats.spearmanr(M3.forest_pct, M3.intensity_kg_ha)
rec('soa_forest_r', f"r = {pt(r_, 2)} (IC 95% {pt(np.tanh(z - 1.96 * se), 2)} a {pt(np.tanh(z + 1.96 * se), 2)}; p = {pt(p_, 3)}); rho = {pt(rho, 2)} (p = {pt(prho, 3)}); n = {len(M3)}", 'table3')
# ---------------------------------------------------------------- 4.1 mesoregions in 2025 (mass and intensity)
Z = np.load(f"{A10}/annual_grids_cache.npz")
g25 = {k: np.nan_to_num(Z[f'{k}_2025'].astype(float)) for k in c.CLASSES}
bv25 = sum(g25[k] for k in c.BVOC); assert abs((bv25 * SC).sum() - G.loc[2025, 'BVOC']) / G.loc[2025, 'BVOC'] < 1e-4
meso = c.read_gpkg(os.path.join(P.BND, 'sc_mesoregions_2024.gpkg')); names = [a['name_meso'] for a, _ in meso]   # mesoregions = unions of the 2024 municipalities (islands included)
W = c.area_weights([g for _, g in meso], lat, lon); HA = {a['name_meso']: a['area_ibge_km2'] * 100 for a, _ in meso}   # IBGE area (ha)
out = []
for i, n in enumerate(names):
    t = {k: float((W[i] @ f.ravel())[0]) for k, f in (('BVOC', bv25), ('ISOP', g25['ISOP']), ('NO', g25['NO']), ('CO', g25['CO']))}; ha = HA[n]
    out.append(dict(meso=n, **{f'{k}_t': v for k, v in t.items()}, **{f'{k}_kg_ha': v * 1e3 / ha for k, v in t.items()}))
MC = pd.DataFrame(out).set_index('meso'); MC.to_csv(f"{OUT}/meso_2025_v10.csv")
rec('meso2025_BVOC_intensity', '; '.join(f"{n} {pt(v, 0)}" for n, v in MC.BVOC_kg_ha.sort_values(ascending=False).items()), 'meso_2025_v10.csv (kg/ha/yr)')
rec('meso2025_NO_intensity', '; '.join(f"{n} {pt(v)}" for n, v in MC.NO_kg_ha.sort_values(ascending=False).items()), 'meso_2025_v10.csv')
rec('meso2025_CO_intensity', '; '.join(f"{n} {pt(v)}" for n, v in MC.CO_kg_ha.sort_values(ascending=False).items()), 'meso_2025_v10.csv')
rec('meso2025_ISOP_intensity_max', MC.ISOP_kg_ha.idxmax(), 'meso_2025_v10.csv')
rec('meso2025_mass_max', f"CO {MC.CO_t.idxmax()}; NO {MC.NO_t.idxmax()}; ISOP {MC.ISOP_t.idxmax()}; BVOC {MC.BVOC_t.idxmax()}", 'meso_2025_v10.csv')
# ---------------------------------------------------------------- 3.4 detrended correlations (notebook 13)
C = pd.read_csv(f"{A10}/correlations_emissions_drivers.csv")
CB = C[C.series.isin(c.BVOC)]
for d in ('TEMP2_C', 'LAI', 'PAR', 'SOILM'):
    s = CB[CB.driver == d]; rec(f'rdet_BVOC_{d}', f"{pt(s.r_detrended.min(), 2)} a {pt(s.r_detrended.max(), 2)}", 'correlations_emissions_drivers.csv')
    s2 = CB[CB.driver == d]; rec(f'rraw_BVOC_{d}', f"{pt(s2.r.min(), 2)} a {pt(s2.r.max(), 2)}", 'correlations_emissions_drivers.csv')
s = C[(C.series == 'NO') & (C.driver == 'TEMP2_C')].iloc[0]; rec('r_NO_T', f"bruto {pt(s.r, 2)}; sem tendência {pt(s.r_detrended, 2)}", 'correlations_emissions_drivers.csv')
for sp in ('SQT_HR', 'SQT_LR'):
    for d in ('TEMP2_C', 'LAI'):
        s = C[(C.series == sp) & (C.driver == d)].iloc[0]; rec(f'rdet_{sp}_{d}', pt(s.r_detrended, 2), 'correlations_emissions_drivers.csv')
# ---------------------------------------------------------------- 4.2 standardised regression and temperature sensitivity
drv = pd.read_csv(f"{A10}/drivers_annual_SC.csv").set_index('year').loc[YEARS]
CL = pd.read_csv(f"{A10}/annual_totals_SC_by_class.csv").set_index('year').loc[YEARS]
CL['SQT'] = CL.SQT_HR + CL.SQT_LR; CL['MT'] = CL[[k for k in CL.columns if k.startswith('MT_')]].sum(1)
zs = lambda a: (np.asarray(a, float) - np.mean(a)) / np.std(a)
keys = ['TEMP2_C', 'PAR', 'LAI', 'SOILM']; X = np.column_stack([zs(drv[k]) for k in keys])
VIF = {}
for i, k in enumerate(keys):
    o = np.delete(X, i, axis=1); XX = np.column_stack([np.ones(len(o)), o]); b = np.linalg.lstsq(XX, X[:, i], rcond=None)[0]
    r2 = 1 - np.sum((X[:, i] - XX @ b) ** 2) / np.sum((X[:, i] - X[:, i].mean()) ** 2); VIF[k] = 1 / (1 - r2)
rec('VIF', '; '.join(f"{k} {pt(v, 2)}" for k, v in VIF.items()), 'drivers_annual_SC.csv')
corr_drv = []
for a_ in keys:
    for b_ in keys:
        if a_ < b_:
            r, p = stats.pearsonr(drv[a_], drv[b_]); corr_drv.append(dict(a=a_, b=b_, r=r, p=p))
CD = pd.DataFrame(corr_drv); CD.to_csv(f"{OUT}/driver_correlations_v10.csv", index=False)
rec('driver_corr', '; '.join(f"{r.a}-{r.b} {pt(r.r, 2)} (p = {pt(r.p, 3)})" for r in CD.itertuples()), 'drivers_annual_SC.csv')
t = np.arange(23); det = lambda x: x - np.polyval(np.polyfit(t, x, 1), t)
r, p = stats.pearsonr(drv.TEMP2_C, drv.SOILM); rd, pd_ = stats.pearsonr(det(drv.TEMP2_C.values), det(drv.SOILM.values))
rec('r_T_SOILM', f"bruto {pt(r, 2)} (p = {pt(p, 2)}); sem tendência {pt(rd, 2)} (p = {pt(pd_, 2)})", 'drivers_annual_SC.csv')
REG, TS = [], []
for sp in list(c.CLASSES) + ['SQT', 'MT']:
    y = zs(CL[sp].values); b = np.linalg.lstsq(X, y, rcond=None)[0]; r2 = 1 - np.sum((y - X @ b) ** 2) / np.sum((y - y.mean()) ** 2)
    REG.append(dict(species=sp, **{f'beta_{k}': bb for k, bb in zip(keys, b)}, R2=r2))
    ly = np.log(CL[sp].values); tt = drv.TEMP2_C.values; sl = stats.linregress(tt, ly); sld = stats.linregress(det(tt), det(ly))
    TS.append(dict(species=sp, pct_per_K=100 * (np.exp(sl.slope) - 1), pct_per_K_se=100 * sl.stderr, pct_per_K_detrended=100 * (np.exp(sld.slope) - 1)))
REG = pd.DataFrame(REG).set_index('species'); REG.to_csv(f"{OUT}/std_regression_v10.csv")
TS = pd.DataFrame(TS).set_index('species'); TS.to_csv(f"{OUT}/temp_sensitivity_v10.csv")
RB = REG.loc[c.BVOC]
for k in keys:
    rec(f'beta_{k}_BVOC', f"{pt(RB[f'beta_{k}'].min(), 2)} a {pt(RB[f'beta_{k}'].max(), 2)} (máx. {RB[f'beta_{k}'].idxmax()}; mín. {RB[f'beta_{k}'].idxmin()})", 'std_regression_v10.csv')
rec('R2_BVOC', f"{pt(RB.R2.min(), 2)} a {pt(RB.R2.max(), 2)}", 'std_regression_v10.csv')
rec('beta_T_is_max_in_all_BVOC', str(bool((RB[[f'beta_{k}' for k in keys]].abs().idxmax(axis=1) == 'beta_TEMP2_C').all())), 'std_regression_v10.csv')
rec('top3_beta_LAI', ', '.join(RB.beta_LAI.sort_values(ascending=False).head(3).index), 'std_regression_v10.csv')
rec('SQT', f"beta_T {pt(REG.loc['SQT', 'beta_TEMP2_C'], 2)}; beta_LAI {pt(REG.loc['SQT', 'beta_LAI'], 2)}; R2 {pt(REG.loc['SQT', 'R2'], 2)}", 'std_regression_v10.csv')
rec('MT', f"beta_T {pt(REG.loc['MT', 'beta_TEMP2_C'], 2)}; beta_LAI {pt(REG.loc['MT', 'beta_LAI'], 2)}; R2 {pt(REG.loc['MT', 'R2'], 2)}", 'std_regression_v10.csv')
rec('NO_R2', f"{pt(REG.loc['NO', 'R2'], 2)}; betas T {pt(REG.loc['NO', 'beta_TEMP2_C'], 2)} PAR {pt(REG.loc['NO', 'beta_PAR'], 2)} LAI {pt(REG.loc['NO', 'beta_LAI'], 2)} SOILM {pt(REG.loc['NO', 'beta_SOILM'], 2)}", 'std_regression_v10.csv')
rec('CO', f"beta_T {pt(REG.loc['CO', 'beta_TEMP2_C'], 2)}; R2 {pt(REG.loc['CO', 'R2'], 2)}", 'std_regression_v10.csv')
rec('isop_pct_per_K', f"{pt(TS.loc['ISOP', 'pct_per_K'])} ± {pt(TS.loc['ISOP', 'pct_per_K_se'])}% por K; sem tendência {pt(TS.loc['ISOP', 'pct_per_K_detrended'])}% por K", 'temp_sensitivity_v10.csv')
# ---------------------------------------------------------------- soil-moisture factor: monthly SOILM below 0.16 and fix-log GAMSM
veg = None; n_below = 0; n_tot = 0; g_lt1 = 0; g_tot = 0; gmin = []
for y in YEARS:
    for m in range(1, 13):
        with rasterio.open(f"{V}/{y}/inputs/met_{y}_{m:02d}.tif") as d:
            sm = d.read(d.descriptions.index('SOILM') + 1)[::-1]
        z = np.load(f"{V}/{y}/inputs/fixlog_{y}_{m:02d}.npz"); gm = z['gamsm_mean']
        ok = SC & np.isfinite(sm); n_tot += ok.sum(); n_below += (sm[ok] < 0.16).sum()
        okg = SC & np.isfinite(gm) & (gm > 0); g_tot += okg.sum(); g_lt1 += (gm[okg] < 0.9999).sum(); gmin.append(np.nanmin(gm[okg]))
rec('soilm_below_016', f"{n_below} de {n_tot} combinações de célula e mês ({pt(100 * n_below / n_tot, 2)}%)", 'met_<year>_<mm>.tif, banda SOILM, células de SC com dado')
rec('gamsm_lt1', f"{g_lt1} de {g_tot} ({pt(100 * g_lt1 / g_tot, 2)}%); mínimo mensal por célula {pt(min(gmin), 3)}", 'fixlog_<year>_<mm>.npz, gamsm_mean')
# ---------------------------------------------------------------- 3.6 municipalities, 2025
MU = pd.read_csv(f"{A10}/table4_municipalities_2025.csv")
mg = c.read_gpkg(f"{c.BND}/SC_Municipios_2024.gpkg"); reg = {a['NM_MUN']: a['NM_RGINT'] for a, _ in mg}
MU['rgint'] = MU.mun.map(reg); MU = MU.sort_values('BVOC', ascending=False); MU.to_csv(f"{OUT}/muni_2025_ranked_v10.csv", index=False)
rec('n_mun', str(len(MU)), 'table4'); rec('mun_top5_BVOC', ', '.join(MU.mun.head(5)), 'table4')
rec('mun_top5_ISOP', ', '.join(MU.sort_values('ISOP', ascending=False).mun.head(5)), 'table4')
rec('mun_top5_same_set', str(set(MU.mun.head(5)) == set(MU.sort_values('ISOP', ascending=False).mun.head(5))), 'table4')
I = MU.intensity_kg_ha
rec('mun_intensity', f"{pt(I.min(), 0)} a {pt(I.max(), 0)} (mediana {pt(I.median(), 0)}; máx. {MU.loc[I.idxmax(), 'mun']}; mín. {MU.loc[I.idxmin(), 'mun']})", 'table4 (área IBGE)')
q80 = I.quantile(0.8); top = MU[I > q80].sort_values('intensity_kg_ha', ascending=False)
rec('mun_top_quintile', f"> {pt(q80, 0)}: " + ', '.join(f"{r.mun} ({pt(r.intensity_kg_ha, 0)}; {r.rgint})" for r in top.head(8).itertuples()), 'table4')
rec('mun_top_quintile_regions', '; '.join(f"{k} {v}" for k, v in top.rgint.value_counts().sort_index().sort_values(ascending=False, kind="stable").items()), 'table4 + NM_RGINT')
q20 = I.quantile(0.2); low = MU[I < q20]
rec('mun_bottom_quintile_regions', f"< {pt(q20, 0)}: " + '; '.join(f"{k} {v}" for k, v in low.rgint.value_counts().sort_index().sort_values(ascending=False, kind="stable").items()), 'table4 + NM_RGINT')
rec('mun_bottom5', ', '.join(f"{r.mun} ({pt(r.intensity_kg_ha, 0)}; {r.rgint})" for r in MU.sort_values('intensity_kg_ha').head(6).itertuples()), 'table4')
low_tot = MU.sort_values('BVOC').head(int(len(MU) * 0.2))
rec('mun_bottom_quintile_total_regions', '; '.join(f"{k} {v}" for k, v in low_tot.rgint.value_counts().sort_index().sort_values(ascending=False, kind="stable").items()), 'table4 + NM_RGINT (20% menores totais)')
top_tot = MU.head(int(len(MU) * 0.2))
rec('mun_top_quintile_total_regions', '; '.join(f"{k} {v}" for k, v in top_tot.rgint.value_counts().sort_index().sort_values(ascending=False, kind="stable").items()), 'table4 + NM_RGINT (20% maiores totais)')
rec('mun_cells_min', str(int(MU.n_cells.min())), 'table4 (células da grade que tocam cada município)')
# ---------------------------------------------------------------- flux, ratio and the 2019 comparison (notebook 13)
MN = pd.read_csv(f"{A10}/key_numbers.csv").set_index('key')
for k in ('isop_flux', 'cmp2019_ISOP', 'cmp2019_CO', 'cmp2019_NO', 'cmp2019_range', 'isop_mean_ratio', 'isop_mt_ratio', 'met_filled', 'lai_gaps',
          'mun_cells', 'mun_area_check', 'meso_sum'):
    rec(k, MN.loc[k, 'pt'], 'key_numbers.csv')
L = pd.DataFrame(LOG); L.to_csv(f"{OUT}/report_numbers_v10.csv", index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 200)
print(L[['key', 'value']].to_string(index=False))
