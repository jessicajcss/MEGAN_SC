"""Post hoc diagnostics of the notebook-14 results (TROPOMI HCHO vs model isoprene, 2019-2025), computed from the cache
files and the model outputs written by the Colab run. Defined after seeing the Tier 1-4 results; reported as such.
D1 reproducibility of the DJF dHCHO spatial pattern (split-half: DJF 2019-2021 vs 2022-2025) at 0.1 and 0.2 deg;
D2 magnitude: SC-mean DJF dHCHO against S*E for Y = 0.29-0.43 and lifetimes 1.5-4 h (mean-based, no regression dilution);
D3 seasonal: raw SC HCHO and Atlantic background separately; dHCHO vs fraction screened; wet-season-only climatology r;
D4 interannual: anomaly correlation without Aug-Oct (main fire season) and in Aug-Oct only (the linear bias correction of
D2 is applied to every column and to the background, so it rescales dHCHO and leaves correlations unchanged)."""
import os, glob, json, calendar, numpy as np, pandas as pd
from netCDF4 import Dataset
from scipy import stats
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import _paths as PATHS
R = PATHS.RUN; TAG = PATHS.TAG; DOUT = PATHS.out('diagnostics', 'nb14')
src = PATHS.notebook_source('14_tropomi_hcho_evaluation.py')
import types, sys; from unittest import mock
g = types.ModuleType('google'); col = types.ModuleType('google.colab'); col.drive = mock.MagicMock(); g.colab = col
sys.modules['google'] = g; sys.modules['google.colab'] = col
cells = [c for c in src.split('# %%') if c.strip()]
code = lambda c: '\n'.join(l for l in c.split('\n') if not l.lstrip().startswith(('%', '!')))
G = {'__name__': '__main__'}
exec(code(cells[1]), G)                                            # settings
G.update(MODEL_ROOT=R, BND=PATHS.BND, OUT_DIR=PATHS.out('_work', 'nb14_diag'))
os.makedirs(os.path.join(G['OUT_DIR'], 'cache_' + TAG), exist_ok=True)
exec(code(cells[2]), G)                                            # model side (MOD, SER, SC, ...)
for f in glob.glob(os.path.join(PATHS.TROPOMI, '*.npz')):             # the cache written by notebook 14 (data/tropomi_hcho_v10)
    import shutil; shutil.copyfile(f, os.path.join(G['OUT_DIR'], 'cache_' + TAG, os.path.basename(f)))
sys.modules['ee'] = mock.MagicMock(); exec(code(cells[3]), G)      # reads the cache only (84 cached)
SAT, MOD, SC, MONTHS, SER = G['SAT'], G['MOD'], G['SC'], G['MONTHS'], G['SER']
block = None; exec(code(cells[4]).split('T1 = []')[0], G); block = G['block']
out = {}
# D1 split-half
dj = lambda yrs: [(y, m) for (y, m) in MONTHS if m in (12, 1, 2) and y in yrs]
H = lambda keys: np.nanmean(np.stack([SAT[k]['H'] - float(SAT[k]['bg']) for k in keys]), axis=0)
a, b = H(dj(range(2019, 2022))), H(dj(range(2022, 2026)))
scf = SC.astype(float)
for k in (4, 8):
    frac = block(scf, k); x = np.where(frac >= 0.5, block(np.where(SC, a, np.nan), k), np.nan); y = np.where(frac >= 0.5, block(np.where(SC, b, np.nan), k), np.nan)
    ok = np.isfinite(x) & np.isfinite(y); r = np.corrcoef(x[ok], y[ok])[0, 1]
    out[f'D1 split-half r of DJF dHCHO pattern, {k * 0.025:.1f} deg blocks (n={ok.sum()})'] = round(r, 3)
    out[f'D1 Spearman-Brown reliability of the 7-year mean, {k * 0.025:.1f} deg'] = round(2 * r / (1 + r), 3)
# D2 magnitude
DJF = dj(range(2019, 2026))
E = np.nanmean(np.stack([MOD[k] for k in DJF]), axis=0); Hd = np.nanmean(np.stack([SAT[k]['H'] - float(SAT[k]['bg']) for k in DJF]), axis=0)
Em, Hm = float(np.nanmean(np.where(SC, E, np.nan))), float(np.nanmean(np.where(SC, Hd, np.nan)))
out['D2 SC-mean model isoprene at overpass, DJF (atoms C cm-2 s-1)'] = f'{Em:.3g}'
out['D2 SC-mean TROPOMI dHCHO, DJF (molec cm-2)'] = f'{Hm:.3g}'
out['D2 mean-based S = dHCHO/E (s); implied lifetime h for Y 0.29 / 0.43'] = f'{Hm / Em:.0f}; {Hm / Em / 0.29 / 3600:.2f} / {Hm / Em / 0.43 / 3600:.2f}'
for Y in (0.29, 0.43):
    for tau in (1.5, 3, 4):
        out[f'D2 expected dHCHO for Y={Y}, tau={tau} h (molec cm-2)'] = f'{Em * Y * tau * 3600:.3g}'
Hbc = np.nanmean(np.stack([(SAT[k]['H'] - 1.10e15) / 0.64 - (float(SAT[k]['bg']) - 1.10e15) / 0.64 for k in DJF]), axis=0)
Hbm = float(np.nanmean(np.where(SC, Hbc, np.nan)))
out['D2 with BIAS_CORRECT (Vigouroux et al., 2020): dHCHO; implied lifetime h, Y 0.29 / 0.43'] = f'{Hbm:.3g}; {Hbm / Em / 0.29 / 3600:.2f} / {Hbm / Em / 0.43 / 3600:.2f}'
Nm = np.nanmean(np.stack([SAT[k]['N'] for k in DJF]), axis=0)
out['D2 kept observations per native pixel per DJF month, SC median'] = round(float(np.nanmedian(np.where(SC, Nm, np.nan))), 1)
# D3 seasonal
S = SER.copy()
c = S.groupby('month')[['ISOP_t', 'HCHO_SC', 'BG', 'dHCHO_SC', 'frac_removed_screening']].mean()
out['D3 climatology r: model isoprene vs raw SC HCHO (no background)'] = round(stats.pearsonr(c.ISOP_t, c.HCHO_SC)[0], 3)
out['D3 climatology r: model isoprene vs Atlantic background'] = round(stats.pearsonr(c.ISOP_t, c.BG)[0], 3)
out['D3 climatology r: dHCHO vs fraction screened (fire + cloud + SZA)'] = round(stats.pearsonr(c.dHCHO_SC, c.frac_removed_screening)[0], 3)
wet = c.loc[[11, 12, 1, 2, 3, 4]]
out['D3 climatology r, Nov-Apr only (n=6): model isoprene vs dHCHO'] = round(stats.pearsonr(wet.ISOP_t, wet.dHCHO_SC)[0], 3)
# D4 anomalies without Aug-Oct and in Aug-Oct only
an = S.copy()
for col_ in ('ISOP_t', 'dHCHO_SC'): an[col_] = S[col_] - S.groupby('month')[col_].transform('mean')
k = ~an.month.isin([8, 9, 10]); rr = G['r_neff'] if 'r_neff' in G else None
exec(code(cells[5]).split('T3 = ')[0], G)
out['D4 anomaly r without Aug-Oct (n=63), r; N_eff; p'] = '; '.join(f'{v:.3f}' if isinstance(v, float) else str(v) for v in [G['r_neff'](an.ISOP_t[k], an.dHCHO_SC[k])[q] for q in ('r', 'n_eff', 'p_neff')])
out['D4 anomaly r, Aug-Oct only (n=21)'] = round(np.corrcoef(an.ISOP_t[~k], an.dHCHO_SC[~k])[0, 1], 3)
for kk, v in out.items(): print(f'{kk}: {v}')
pd.Series(out).to_csv(os.path.join(DOUT, 'nb14_diagnostics.csv'), header=['value'])
# D5 terrain: TROPOMI tropospheric columns scale with the air below the retrieval, so a stable spatial pattern can follow
# surface pressure (elevation) rather than emissions. DJF-mean surface pressure from the model input maps.
import rasterio
def band(path, name):
    with rasterio.open(path) as s_:
        return np.flipud(s_.read(list(s_.descriptions).index(name) + 1).astype('float64'))
P = np.nanmean(np.stack([band(f'{R}/{y}/inputs/met_{y}_{m:02d}.tif', 'PRES') for (y, m) in DJF]), axis=0)
rows5 = {}
for kb in (4, 8):
    frac = block(scf, kb)
    xs = {n: np.where(frac >= 0.5, block(np.where(SC, f_, np.nan), kb), np.nan).ravel() for n, f_ in (('H', Hd), ('E', E), ('P', P))}
    ok = np.isfinite(xs['H']) & np.isfinite(xs['E']) & np.isfinite(xs['P'])
    h, e, p = xs['H'][ok], xs['E'][ok], xs['P'][ok]
    r_hp = np.corrcoef(h, p)[0, 1]; r_he = np.corrcoef(h, e)[0, 1]; r_ep = np.corrcoef(e, p)[0, 1]
    pr = (r_he - r_hp * r_ep) / np.sqrt((1 - r_hp ** 2) * (1 - r_ep ** 2))
    X = np.column_stack([np.ones_like(p), (p - p.mean()) / p.std(), (e - e.mean()) / e.std()]); beta = np.linalg.lstsq(X, h, rcond=None)[0]
    rows5[f'{kb * 0.025:.1f} deg'] = dict(n=int(ok.sum()), r_dHCHO_pressure=round(r_hp, 3), r_dHCHO_isoprene=round(r_he, 3), r_isoprene_pressure=round(r_ep, 3),
                                          partial_r_dHCHO_isoprene_given_pressure=round(pr, 3), beta_pressure_per_sd=f'{beta[1]:.3g}', beta_isoprene_per_sd=f'{beta[2]:.3g}')
D5 = pd.DataFrame(rows5).T; print(D5.to_string()); D5.to_csv(os.path.join(DOUT, 'nb14_D5_terrain.csv'))
# D6 spatial magnitude at 0.2 deg: 1-SD spatial variation of S*E (tau 3 h) vs the observed spatial SD of dHCHO, and the
# observed dHCHO change per 1 SD of E after accounting for surface pressure (beta_isoprene_per_sd of D5).
frac = block(scf, 8)
eb = np.where(frac >= 0.5, block(np.where(SC, E, np.nan), 8), np.nan); hb = np.where(frac >= 0.5, block(np.where(SC, Hd, np.nan), 8), np.nan)
ok = np.isfinite(eb) & np.isfinite(hb)
for Y in (0.29, 0.43):
    print(f'D6 0.2 deg: SD of E {np.std(eb[ok]):.3g} atoms C cm-2 s-1 -> expected SD of dHCHO with Y={Y}, tau=3 h: {np.std(eb[ok]) * Y * 3 * 3600:.3g}')
print(f'D6 0.2 deg: observed SD of dHCHO {np.std(hb[ok]):.3g}; observed dHCHO per 1 SD of E given pressure: {rows5["0.2 deg"]["beta_isoprene_per_sd"]}')
# D7 smearing beyond Tier 1: partial r (dHCHO, smoothed E | surface pressure) at 0.2 deg for kernels up to 100 km.
from scipy import ndimage
frac = block(scf, 8); pb = np.where(frac >= 0.5, block(np.where(SC, P, np.nan), 8), np.nan)
for sig in (0, 20, 40, 60, 100):
    Es = ndimage.gaussian_filter(np.nan_to_num(E), sig / (G['RES'] * 111.0)) if sig else E
    eb = np.where(frac >= 0.5, block(np.where(SC, Es, np.nan), 8), np.nan); ok = np.isfinite(eb) & np.isfinite(hb) & np.isfinite(pb)
    h, e, p = hb[ok], eb[ok], pb[ok]; r_hp, r_he, r_ep = np.corrcoef(h, p)[0, 1], np.corrcoef(h, e)[0, 1], np.corrcoef(e, p)[0, 1]
    print(f'D7 smoothing {sig:3d} km: r(dHCHO,E) {r_he:+.3f} | partial r given pressure {(r_he - r_hp * r_ep) / np.sqrt((1 - r_hp ** 2) * (1 - r_ep ** 2)):+.3f}')
