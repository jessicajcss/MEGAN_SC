"""TROPOMI figure (PT and EN; report Figure A.10; its panels are article Figures 10 and 11, drawn in the journal format
by ijest/ijest_figures.py): TROPOMI HCHO against modelled isoprene, 2019-2025.
Panels a, b, d, e are the pre-specified notebook-14 views (DJF maps, seasonal cycle, deseasonalised anomalies);
panels c and f are diagnostics defined after the pre-specified tests (surface pressure as terrain proxy).
Inputs: the production run (data/run_v10), the TROPOMI cache of notebook 14 (data/tropomi_hcho_v10) and the boundaries.
Usage: python3 fig9_tropomi.py [out_dir]   (default reproduced/figures; writes <out_dir>/pt/fig9.png and <out_dir>/en/fig9.png)"""
import os, sys, glob, types, shutil
import numpy as np, pandas as pd, rasterio
from unittest import mock
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import geopandas as gpd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'maps')); import wf_style as W   # map frame of the biomass-burning article
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import _paths as PATHS
R = PATHS.RUN; TAG = PATHS.TAG; OUT = sys.argv[1] if len(sys.argv) > 1 else PATHS.out('figures')
src = open(os.path.join(PATHS.NB_SRC, '14_tropomi_hcho_evaluation.py'), encoding='utf-8').read()
g = types.ModuleType('google'); col = types.ModuleType('google.colab'); col.drive = mock.MagicMock(); g.colab = col
sys.modules['google'] = g; sys.modules['google.colab'] = col; sys.modules['ee'] = mock.MagicMock()
cells = [c for c in src.split('# %%') if c.strip()]
code = lambda c: '\n'.join(l for l in c.split('\n') if not l.lstrip().startswith(('%', '!')))
G = {'__name__': '__main__'}
exec(code(cells[1]), G)
work = PATHS.out('_work', 'fig9'); os.makedirs(os.path.join(work, 'cache_' + TAG), exist_ok=True)
G.update(MODEL_ROOT=R, BND=PATHS.BND, OUT_DIR=work)
exec(code(cells[2]), G)
for f in glob.glob(os.path.join(PATHS.TROPOMI, '*.npz')): shutil.copyfile(f, os.path.join(work, 'cache_' + TAG, os.path.basename(f)))
exec(code(cells[3]), G); exec(code(cells[4]).split('T1 = []')[0], G); exec(code(cells[5]).split('T3 = ')[0], G)
SAT, MOD, SC, MONTHS, SER, lat, lon, RES, block = (G[k] for k in ('SAT', 'MOD', 'SC', 'MONTHS', 'SER', 'lat', 'lon', 'RES', 'block'))
DJF = [(y, m) for (y, m) in MONTHS if m in (12, 1, 2)]
E = np.nanmean(np.stack([MOD[k] for k in DJF]), axis=0)
H = np.nanmean(np.stack([SAT[k]['H'] - float(SAT[k]['bg']) for k in DJF]), axis=0)
def band(path, name):
    with rasterio.open(path) as s_: return np.flipud(s_.read(list(s_.descriptions).index(name) + 1).astype('float64'))
P = np.nanmean(np.stack([band(f'{R}/{y}/inputs/met_{y}_{m:02d}.tif', 'PRES') for (y, m) in DJF]), axis=0) / 100.0   # hPa
S = SER.copy(); clim = S.groupby('month')[['ISOP_t', 'dHCHO_SC', 'BG']].mean()
an = S.copy()
for c in ('ISOP_t', 'dHCHO_SC'): an[c] = S[c] - S.groupby('month')[c].transform('mean')
rr = G['r_neff'](an.ISOP_t, an.dHCHO_SC)
r2 = np.corrcoef(clim.ISOP_t, clim.dHCHO_SC)[0, 1]
scf = SC.astype(float); frac = block(scf, 8)
hb = np.where(frac >= 0.5, block(np.where(SC, H, np.nan), 8), np.nan).ravel(); pb = np.where(frac >= 0.5, block(np.where(SC, P, np.nan), 8), np.nan).ravel()
eb = np.where(frac >= 0.5, block(np.where(SC, E, np.nan), 8), np.nan).ravel(); ok = np.isfinite(hb) & np.isfinite(pb) & np.isfinite(eb)
r_hp = np.corrcoef(hb[ok], pb[ok])[0, 1]
L = {'pt': dict(a='(a) Isopreno modelado na passagem do satélite,\nDJF 2019–2025 (10¹² átomos C cm⁻² s⁻¹)',
                b='(b) ΔHCHO do TROPOMI, DJF 2019–2025\n(10¹⁵ moléculas cm⁻²)',
                c='(c) Pressão à superfície, DJF (hPa)\n(diagnóstico definido após os testes)',
                d='(d) Ciclo sazonal, médias 2019–2025 (normalizadas pelo máximo)', e='(e) Anomalias mensais dessazonalizadas (escore z)',
                f='(f) ΔHCHO × pressão à superfície, blocos de 0,2°\n(diagnóstico definido após os testes)',
                iso='isopreno modelado', hcho='ΔHCHO do TROPOMI', bg='fundo atlântico (HCHO)', month='mês',
                cb='isopreno modelado (10¹² átomos C cm⁻² s⁻¹)', xp='pressão à superfície (hPa)', yh='ΔHCHO (10¹⁵ moléculas cm⁻²)',
                dec=','),
     'en': dict(a='(a) Modelled isoprene at overpass,\nDJF 2019–2025 (10¹² atoms C cm⁻² s⁻¹)',
                b='(b) TROPOMI ΔHCHO, DJF 2019–2025\n(10¹⁵ molecules cm⁻²)',
                c='(c) Surface pressure, DJF (hPa)\n(diagnostic defined after the tests)',
                d='(d) Seasonal cycle, 2019–2025 means (normalised by the maximum)', e='(e) Deseasonalised monthly anomalies (z-score)',
                f='(f) ΔHCHO vs surface pressure, 0.2° blocks\n(diagnostic defined after the tests)',
                iso='modelled isoprene', hcho='TROPOMI ΔHCHO', bg='Atlantic background (HCHO)', month='month',
                cb='modelled isoprene (10¹² atoms C cm⁻² s⁻¹)', xp='surface pressure (hPa)', yh='ΔHCHO (10¹⁵ molecules cm⁻²)',
                dec='.')}
def num(x, d, lang): return f'{x:.{d}f}'.replace('.', L[lang]['dec']).replace('-', '\u2212')
ext = [lon[0] - RES / 2, lon[-1] + RES / 2, lat[0] - RES / 2, lat[-1] + RES / 2]
state = W.state_from(gpd.read_file(os.path.join(PATHS.BND, 'SC_Municipios_2024.gpkg')))
for lang in ('pt', 'en'):
    l = L[lang]; plt.rcParams.update({'font.size': 9})
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.6))
    for ax, fld, ttl, cm in ((axs[0, 0], np.where(SC, E, np.nan) / 1e12, l['a'], 'Greens'), (axs[0, 1], np.where(SC, H, np.nan) / 1e15, l['b'], 'YlOrBr'),
                             (axs[0, 2], np.where(SC, P, np.nan), l['c'], 'cividis')):
        im = ax.imshow(fld, origin='lower', extent=ext, cmap=cm, interpolation='nearest'); fig.colorbar(im, ax=ax, shrink=0.72)
        state.boundary.plot(ax=ax, color='0.2', lw=0.6); W.extent(ax); ax.tick_params(labelsize=8)
        ax.xaxis.set_major_formatter(W.fmt_deg(lang)); ax.yaxis.set_major_formatter(W.fmt_deg(lang))
        W.north_arrow(ax, fs=7); W.scale_bar(ax, small=True); ax.set_title(ttl, loc='left', fontsize=9.5)
    ax = axs[1, 0]; mm = np.arange(1, 13)
    ax.plot(mm, clim.ISOP_t / clim.ISOP_t.max(), 'o-', color='#2e8b57', label=l['iso'])
    ax.plot(mm, clim.dHCHO_SC / clim.dHCHO_SC.max(), 's-', color='#b8860b', label=l['hcho'])
    ax.plot(mm, clim.BG / clim.BG.max(), '^--', color='0.45', label=l['bg'])
    ax.set_xticks(mm); ax.set_xlabel(l['month']); ax.set_title(l['d'] + f'\nr (isop., ΔHCHO) = {num(r2, 2, lang)}', loc='left', fontsize=9.5); ax.legend(frameon=False, fontsize=8)
    ax = axs[1, 1]; t = pd.to_datetime(dict(year=S.year, month=S.month, day=15)); z = lambda s: (s - s.mean()) / s.std()
    ax.plot(t, z(an.ISOP_t), color='#2e8b57', label=l['iso']); ax.plot(t, z(an.dHCHO_SC), color='#b8860b', label=l['hcho']); ax.axhline(0, color='0.6', lw=0.8)
    ptxt = 'p < 0' + l['dec'] + '001' if rr['p_neff'] < 0.001 else f"p = {num(rr['p_neff'], 3, lang)}"
    ax.set_title(l['e'] + f"\nr = {num(rr['r'], 2, lang)} (N_eff = {rr['n_eff']:.0f}; {ptxt})", loc='left', fontsize=9.5); ax.legend(frameon=False, fontsize=8)
    ax = axs[1, 2]; sc_ = ax.scatter(pb[ok], hb[ok] / 1e15, c=eb[ok] / 1e12, cmap='Greens', s=14, edgecolor='0.3', linewidth=0.3)
    fig.colorbar(sc_, ax=ax, shrink=0.72, label=l['cb']); ax.set_xlabel(l['xp']); ax.set_ylabel(l['yh'])
    ax.set_title(l['f'] + f'\nr = {num(r_hp, 2, lang)} (n = {int(ok.sum())})', loc='left', fontsize=9.5)
    if lang == 'pt':                                   # decimal comma on every numeric axis (not on the date axis of panel e)
        comma = matplotlib.ticker.FuncFormatter(lambda v, p: f'{v:g}'.replace('.', ',').replace('-', '\u2212'))
        for a_ in fig.axes:
            a_.yaxis.set_major_formatter(comma)
            if a_ is not axs[1, 1] and a_ not in (axs[0, 0], axs[0, 1], axs[0, 2]): a_.xaxis.set_major_formatter(comma)
        for a_ in (axs[0, 0], axs[0, 1], axs[0, 2]): a_.yaxis.set_major_formatter(W.fmt_deg(lang))
    fig.tight_layout(); os.makedirs(os.path.join(OUT, lang), exist_ok=True); fig.savefig(os.path.join(OUT, lang, 'fig9.png'), dpi=250); plt.close(fig)
    print('wrote', os.path.join(OUT, lang, 'fig9.png'))
print(f'check: r_clim {r2:.3f} | anomalies r {rr["r"]:.3f} N_eff {rr["n_eff"]:.1f} p {rr["p_neff"]:.5f} | r(dHCHO, p_s) 0.2 deg {r_hp:.3f} n {int(ok.sum())}')
