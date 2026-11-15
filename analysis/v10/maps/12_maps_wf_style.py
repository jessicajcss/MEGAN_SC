#!/usr/bin/env python3
"""12_maps_wf_style.py - every map of the biogenic manuscript and report in the layout of the biomass-burning article
(wf_style.py), with the colour ramps of the biogenic figures, the scale bar inside the frame at the bottom right of every
map, mesoregion labels at the interior point farthest from the edges, Jenks classes in the municipal maps, one source note
per figure and one legend per row when the panels of a row share their classes, and larger lettering in Figure 5; in fig8
the legend of each row is in its left map, the only place where the compound is named.

Manuscript (Portuguese and English; article Figures 7, 8 and 9, redrawn in the journal format by ijest/ijest_figures.py;
also report Figures 7, 9 and 10):
  fig5  SOA-precursor intensity by mesoregion, 2003-2025 mean (article Table 2, report Table 3; mesoregions with the coastal islands)
  fig7  BVOC by municipality, 2025: a) total, b) intensity (Jenks natural breaks, 5 classes)
  fig8  isoprene, CO and soil NO by municipality, 2019 and 2025 (Jenks natural breaks of both years together)
Report only (Portuguese):
  fig9   2025 intensity by group of compounds (quintile classes of the state cells, per group; report Figure 11)
  fig12  a) MapBiomas tree fraction, b) SOA-precursor intensity in 2025 (report Figure 14)
  figA2..figA8  monthly input fields of 2025 (equal-interval classes shared by the 12 months)
  figA9  BVOC intensity by year (quintile classes of the state cells of all years)
Inputs: the outputs of notebook 13 with sc_mesoregions_2024.gpkg (MEGAN_SC_NB13: tables and annual_grids_cache.npz),
data/boundaries, and the 2025 inputs of the production run (data/run_v10/2025/inputs).
Outputs: reproduced/figures/{pt,en}/fig5.png, fig7.png, fig8.png; reproduced/figures/pt_report/fig9.png, fig12.png,
figA2-figA9.png; class limits in reproduced/figures/maps_class_limits.csv (published copies in results/v10/figures)."""
import os, sys, textwrap
import numpy as np, pandas as pd, geopandas as gpd, rasterio, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import wf_style as W
# paths from analysis/v10/_paths.py (environment variables MEGAN_SC_*); common.py (grid, state mask, cell areas) of
# analysis/v10/
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'analysis', 'v10')); import common as c
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
V = P.RUN; NB13 = P.NB13; BND = P.BND
FIG = os.environ.get('MEGAN_SC_FIGS') or P.out('figures')
LIM = []                                                       # class limits written to maps_class_limits.csv

mun = gpd.read_file(os.path.join(BND, 'SC_Municipios_2024.gpkg')).to_crs(4326)
mesog = gpd.read_file(os.path.join(BND, 'sc_mesoregions_2024.gpkg')).to_crs(4326)
mesog['geometry'] = [W.fill_slivers(g) for g in mesog.geometry]
state = W.state_from(mun)
T3 = pd.read_csv(os.path.join(NB13, 'table3_mesoregions_soa.csv')).set_index('mesoregion')
MU = pd.read_csv(os.path.join(NB13, 'table4_municipalities_2025.csv')).set_index('mun')
MU19 = pd.read_csv(os.path.join(NB13, 'municipalities_2019.csv')).set_index('mun')
for d_ in (MU, MU19): assert set(d_.index) == set(mun.NM_MUN), 'municipality names'

def qbreaks(v, n=5):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; return np.quantile(v, np.linspace(0, 1, n + 1))
def dfor(x):                                                   # decimals for a class limit
    ax_ = abs(x); return 0 if ax_ >= 20 else (1 if ax_ >= 2 else 2)
def labels(b, lang, lo_open=False, hi_open=False, d=None):
    """One number of decimals per legend: from the first interior limit (or given, for equal-interval classes)."""
    a = ' a ' if lang == 'pt' else ' to '; out = []
    d = dfor(b[1]) if d is None else d
    for i in range(len(b) - 1):
        lo, hi = W.num(b[i], d, lang), W.num(b[i + 1], d, lang)
        if i == 0 and lo_open: out.append(f'< {hi}')
        elif i == len(b) - 2 and hi_open: out.append(f'≥ {lo}')
        elif i == len(b) - 2: out.append(f'{lo}{a}{hi}')
        else: out.append(f'{lo}{a}< {hi}')
    return out
def nice_step(lo, hi, n=7):
    raw = (hi - lo) / n; k = 10 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * k >= raw: return m * k
def ebreaks(lo, hi, n=7):
    s = nice_step(lo, hi, n); a = np.floor(lo / s) * s; b = np.ceil(hi / s) * s
    return np.round(np.arange(a, b + s / 2, s), 10), (0 if float(s).is_integer() else len(f'{s:g}'.split('.')[1]))

TX = {'pt': dict(yr='ano⁻¹', t='t ano⁻¹', kg='kg ha⁻¹ ano⁻¹', nmun='n.º de municípios',
                 f5cb='Intensidade de precursores de AOS (kg ha⁻¹ ano⁻¹)',
                 n5='MEGANv3.21; média de 2003 a 2025; intensidade pela área territorial (IBGE)\n'
                    'Mesorregiões: agregação dos municípios (malha municipal IBGE 2024), com as ilhas costeiras; coordenadas geográficas, WGS 84',
                 f7a='a) Emissão total de COVB (t ano⁻¹)', f7b='b) Intensidade de COVB (kg ha⁻¹ ano⁻¹)',
                 n7='MEGANv3.21, 2025; emissão de cada célula da grade de 0,025° repartida entre os municípios pela área contida em cada um\n'
                    'Classes: quebras naturais de Jenks; malha municipal IBGE 2024; coordenadas geográficas, WGS 84',
                 iso='Isopreno', co='CO', no='NO do solo',
                 y19='2019', y25='2025',
                 n8='MEGANv3.21. Classes: quebras naturais de Jenks, comuns a 2019 e 2025 em cada composto (legenda à esquerda de cada linha). Malha municipal IBGE 2024; coordenadas geográficas, WGS 84'),
      'en': dict(yr='yr⁻¹', t='t yr⁻¹', kg='kg ha⁻¹ yr⁻¹', nmun='no. of municipalities',
                 f5cb='SOA-precursor intensity (kg ha⁻¹ yr⁻¹)',
                 n5='MEGANv3.21; 2003-2025 mean; intensity over the territorial area (IBGE)\n'
                    'Mesoregions: union of the municipalities (IBGE 2024 municipal mesh), coastal islands included; geographic coordinates, WGS 84',
                 f7a='a) Total BVOC emission (t yr⁻¹)', f7b='b) BVOC intensity (kg ha⁻¹ yr⁻¹)',
                 n7='MEGANv3.21, 2025; emission of each 0.025° grid cell split among the municipalities by the area inside each one\n'
                    'Classes: Jenks natural breaks; IBGE 2024 municipal mesh; geographic coordinates, WGS 84',
                 iso='Isoprene', co='CO', no='Soil NO',
                 y19='2019', y25='2025',
                 n8='MEGANv3.21. Classes: Jenks natural breaks, shared by 2019 and 2025 for each compound (legend at the left of each row). IBGE 2024 municipal mesh; geographic coordinates, WGS 84')}
# label of each mesoregion at its interior point farthest from the edges (W.label_point); the two narrow coastal
# mesoregions are labelled over the Atlantic with a leader line (text centre, offset from the label point)
OCEAN_LABEL = {'Grande Florianópolis': (-47.72, -27.22), 'Sul Catarinense': (-48.05, -28.95)}
F5_XLIM = (-54.0, -47.0)                                       # Figure 5 only: room for the two ocean labels
for lang, L in TX.items():
    od = os.path.join(FIG, lang); os.makedirs(od, exist_ok=True)
    # ---- Figure 5: mesoregions (layout of Figure 3 of the biomass-burning article)
    g = mesog.assign(v=mesog.name_meso.map(T3.intensity_kg_ha), t=mesog.name_meso.map(T3.SOA)); assert g.v.notna().all()
    bnds = np.arange(3 * np.floor(g.v.min() / 3), 3 * np.ceil(g.v.max() / 3) + 1, 3)
    cmap = ListedColormap(W.ramp('YlOrRd', len(bnds) - 1)); norm = BoundaryNorm(bnds, cmap.N)
    # lettering sized for a figure printed at about 16 cm wide (readable at 100 % to 130 % zoom); each label holds the name and
    # the intensity (the totals are in Figure 4 and Table 3)
    fig, ax = plt.subplots(figsize=(8.6, 5.6), dpi=300)
    g.plot(ax=ax, column='v', cmap=cmap, norm=norm, edgecolor='0.25', lw=0.6); state.boundary.plot(ax=ax, color='0.2', lw=0.7)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax, shrink=0.85, format=FuncFormatter(lambda v, _: W.num(v, 0, lang)))
    cb.set_label(L['f5cb'], fontsize=11); cb.ax.tick_params(labelsize=10)
    for r in g.itertuples():
        x, y = W.label_point(r.geometry); dark = norm(r.v) >= len(bnds) - 2
        txt = f"{r.name_meso.replace('Grande ', 'Grande' + chr(10))}\n{W.num(r.v, 1, lang)}"
        if r.name_meso in OCEAN_LABEL:
            tx, ty = OCEAN_LABEL[r.name_meso]
            ax.annotate(txt, xy=(x, y), xytext=(tx, ty), ha='center', va='center', fontsize=10, color='black', linespacing=1.15,
                        arrowprops=dict(arrowstyle='-', color='0.25', lw=0.7, shrinkA=0, shrinkB=0))
            ax.plot([x], [y], 'o', ms=2.4, color='0.15', zorder=5)
        else:
            ax.text(x, y, txt, ha='center', va='center', fontsize=10, color='white' if dark else 'black', linespacing=1.15,
                    path_effects=None if dark else W.halo(lw=2.0))
    W.base(ax, L['n5'], lang, xlim=F5_XLIM, fs=9, note_fs=7.5, note_y=-0.16, small=True); ax.tick_params(labelsize=10)
    ax.xaxis.label.set_size(11); ax.yaxis.label.set_size(11)
    fig.tight_layout(); W.save(fig, os.path.join(od, 'fig5.png'))
    LIM.append(dict(figure='fig5', panel='SOA intensity', limits=list(bnds)))
    # ---- Figure 7: municipalities, 2025 (two stacked panels, as Figure 6 of the biomass-burning article)
    fig, axs = plt.subplots(2, 1, figsize=(8.6, 11.2), dpi=300)
    for ax, col, ttl, unit, cm in ((axs[0], 'BVOC', L['f7a'], L['t'], 'Greens'), (axs[1], 'intensity_kg_ha', L['f7b'], L['kg'], 'YlGn')):
        v = mun.NM_MUN.map(MU[col]).values; q = W.fisher_jenks(v, 5)
        W.class_choropleth(ax, mun, v, q, W.ramp(cm, 5), labels(q, lang), unit)
        W.outlines(ax, state); ax.set_title(ttl, loc='left', fontsize=9); W.base(ax, None, lang)
        LIM.append(dict(figure='fig7', panel=col, limits=[round(float(x), 3) for x in q]))
    axs[1].text(0.0, -0.105, L['n7'], transform=axs[1].transAxes, fontsize=6.3, va='top', ha='left')   # one note for both panels
    fig.tight_layout(h_pad=1.5); W.save(fig, os.path.join(od, 'fig7.png'))
    # ---- Figure 8: isoprene, CO and soil NO, 2019 and 2025 (grid as Figure 6 of the biomass-burning article; one compound per row,
    #      2019 and 2025 side by side; cartographic elements and source in the bottom row)
    fig = plt.figure(figsize=(8.6, 9.0), dpi=300); gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 0.16], hspace=0.03, wspace=0.03)
    for j, (k, cm, nm) in enumerate((('ISOP', 'Greens', L['iso']), ('CO', 'Blues', L['co']), ('NO', 'Reds', L['no']))):
        q = W.fisher_jenks(np.concatenate([mun.NM_MUN.map(MU19[k]).values, mun.NM_MUN.map(MU[k]).values]), 5)
        LIM.append(dict(figure='fig8', panel=k, limits=[round(float(x), 3) for x in q]))
        for i, (yy, T_) in enumerate(((2019, MU19), (2025, MU))):
            ax = fig.add_subplot(gs[j, i]); W.extent(ax, axis_off=True)
            v = mun.NM_MUN.map(T_[k]).values; cols = W.ramp(cm, 5); kk = np.clip(np.searchsorted(q, v, side='right') - 1, 0, 4)
            for ci in range(5):
                sub = mun[kk == ci]
                if len(sub): sub.plot(ax=ax, color=cols[ci], edgecolor=cols[ci], lw=0.15)
            state.boundary.plot(ax=ax, color='0.1', lw=0.8)
            if j == 0: ax.set_title(L['y19'] if yy == 2019 else L['y25'], fontsize=10, fontweight='bold')
            if i == 0:                                         # one legend per row, in the left map: both years share the classes
                hs = [Patch(facecolor=cols[ci], edgecolor='0.45', lw=0.4, label=lb) for ci, lb in enumerate(labels(q, lang))]
                ax.legend(handles=hs, title=f"{nm} ({L['t']})", loc='lower left', bbox_to_anchor=(-0.02, -0.04), frameon=False,
                          fontsize=6.6, title_fontsize=8.0, handlelength=1.6, handleheight=0.9, labelspacing=0.28, alignment='left')
            if j == 2 and i == 1:                              # north arrow and scale bar inside the last map, over the Atlantic
                W.north_arrow(ax, x=-48.35, y0=-28.95, y1=-28.45, fs=8); W.scale_bar(ax, right=-48.25, lat_sb=-29.4, small=True, fs=6.0)
    ax = fig.add_subplot(gs[3, :]); ax.axis('off'); ax.text(0.0, 0.9, L['n8'], fontsize=6.4, va='top', ha='left', transform=ax.transAxes, wrap=True)
    W.save(fig, os.path.join(od, 'fig8.png'))
    print(lang, 'figures 5, 7, 8 ->', od)

# ======================= report-only maps (Portuguese) on the model grid
lat, lon = c.grid(); SC = c.sc_mask(lat, lon); A = c.cell_area_m2(lat, lon)
EXT = [lon[0] - 0.0125, lon[-1] + 0.0125, lat[0] - 0.0125, lat[-1] + 0.0125]
Z = np.load(P.grid_cache())
def grid_year(y): return {k: np.nan_to_num(Z[f'{k}_{y}'].astype(float)) for k in c.CLASSES}
def groups(g):
    return {'COVB total': sum(g[k] for k in c.BVOC), 'Isopreno': g['ISOP'], 'Monoterpenos': sum(g[k] for k in c.CLASSES if k.startswith('MT_')),
            'Sesquiterpenos': g['SQT_HR'] + g['SQT_LR'], 'COV oxigenados': sum(g[k] for k in ('MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD')),
            'HC leves e outros': g['LVOC'] + g['STRESS'] + g['OTHER'] + g['MBO'], 'NO do solo': g['NO'], 'CO': g['CO']}
def gridmap(ax, fld, b, cmap_name, lw=0.6):
    cmap = ListedColormap(W.ramp(cmap_name, len(b) - 1)); cmap.set_bad((1, 1, 1, 0)); norm = BoundaryNorm(b, cmap.N, extend='neither')
    f = np.ma.masked_invalid(np.clip(fld, b[0], b[-1]))
    ax.imshow(f, origin='lower', extent=EXT, cmap=cmap, norm=norm, interpolation='nearest')
    state.boundary.plot(ax=ax, color='0.1', lw=lw); return cmap
def legend_grid(ax, b, cmap_name, title, lang='pt', lo_open=False, hi_open=False, fs=6.0, tfs=7.2, anchor=(-0.02, -0.04), loc='lower left'):
    cols = W.ramp(cmap_name, len(b) - 1)
    hs = [Patch(facecolor=cols[i], edgecolor='0.45', lw=0.4, label=lb) for i, lb in enumerate(labels(b, lang, lo_open, hi_open))]
    ax.legend(handles=hs, title=title, loc=loc, bbox_to_anchor=anchor, frameon=False, fontsize=fs, title_fontsize=tfs,
              handlelength=1.4, handleheight=0.8, labelspacing=0.22, alignment='left')
def note_below(fig, text, y=0.0, width=170):
    fig.text(0.02, y, '\n'.join(textwrap.wrap(text, width)), fontsize=6.3, va='top', ha='left')
OUTR = os.path.join(FIG, 'pt_report'); os.makedirs(OUTR, exist_ok=True)
NOTE_GRID = 'Contorno: malha municipal IBGE 2024; coordenadas geográficas, WGS 84'

# ---- Figure 9: 2025 intensity by group (3 x 3 grid; cartographic elements and source in the free panel)
g25 = groups(grid_year(2025)); cms = ['YlGn', 'Blues', 'Oranges', 'YlOrBr', 'Greens', 'RdPu', 'Reds', 'BuGn']
fig, axs = plt.subplots(3, 3, figsize=(8.6, 6.6), dpi=300)
for ax in axs.ravel(): W.extent(ax, axis_off=True)
for ax, (k, v), cm in zip(axs.ravel(), g25.items(), cms):
    inten = np.where(SC, v * 1e3 / (A / 1e4), np.nan); q = qbreaks(inten[SC])
    gridmap(ax, inten, q, cm); legend_grid(ax, q, cm, f"{k}\n{W.num(np.nansum(v * SC) / 1e3, 1)} Gg ano⁻¹ no Estado")
    LIM.append(dict(figure='report fig9', panel=k, limits=[round(float(x), 4) for x in q]))
ax = axs[2, 2]; W.north_arrow(ax, x=-48.45, y0=-28.6, y1=-27.9, fs=9); W.scale_bar(ax, right=-48.2, lat_sb=-29.3, small=True)
ax.text(-54.0, -26.0, 'MEGANv3.21, 2025\nIntensidade por célula da\ngrade de 0,025° (kg ha⁻¹ ano⁻¹)\nClasses: quintis das células\nde SC em cada grupo\nMalha municipal IBGE 2024;\nWGS 84', fontsize=6.4, va='top', ha='left')
fig.tight_layout(w_pad=0.2, h_pad=0.2); W.save(fig, os.path.join(OUTR, 'fig9.png'))

# ---- Figure 12: tree fraction and SOA-precursor intensity, 2025 (two stacked panels with geographic frame)
tree = np.load(os.path.join(V, '2025', 'inputs', 'growthform_2024.npz'))['tree']; assert tree.shape == SC.shape
soa = np.where(SC, (g25['Monoterpenos'] + g25['Sesquiterpenos']) * 1e3 / (A / 1e4), np.nan)
fig, axs = plt.subplots(2, 1, figsize=(8.6, 11.2), dpi=300)
for ax, fld, b, cm, ttl, lt, note in (
        (axs[0], np.where(SC, tree, np.nan), np.linspace(0, 1, 6), 'Greens', 'a) Fração arbórea (MapBiomas 30 m, 2024)', 'Fração de árvore na célula',
         'a) MapBiomas Coleção 10 (30 m), 2024, agregado à grade de 0,025° (fração de árvore usada pelo modelo). '),
        (axs[1], soa, qbreaks(soa[SC]), 'YlOrBr', 'b) Precursores de AOS (monoterpenos + sesquiterpenos), 2025', 'kg ha⁻¹ ano⁻¹ (quintis das células)',
         'b) MEGANv3.21, 2025; intensidade por célula da grade de 0,025°.')):
    gridmap(ax, fld, b, cm); mesog.boundary.plot(ax=ax, color='0.35', lw=0.45); state.boundary.plot(ax=ax, color='0.2', lw=0.6)
    legend_grid(ax, b, cm, lt, fs=6.3, tfs=6.8, anchor=(0.0, 0.02)); ax.set_title(ttl, loc='left', fontsize=9); W.base(ax, None)
    LIM.append(dict(figure='report fig12', panel=ttl[:2], limits=[round(float(x), 4) for x in b])); NOTE12 = (NOTE12 if ax is not axs[0] else '') + note
axs[1].text(0.0, -0.105, NOTE12 + '\n' + NOTE_GRID, transform=axs[1].transAxes, fontsize=6.3, va='top', ha='left')   # one note for both panels
fig.tight_layout(h_pad=1.5); W.save(fig, os.path.join(OUTR, 'fig12.png'))

# ---- Figures A.2 to A.8: monthly input fields of 2025 (3 x 4 maps, legend and cartographic elements below)
MES = ['jan', 'fev', 'mar', 'abr', 'mai', 'jun', 'jul', 'ago', 'set', 'out', 'nov', 'dez']
def monthly(var):
    out = []
    for m in range(1, 13):
        f = os.path.join(V, '2025', 'inputs', f"{'lai' if var == 'LAI' else 'met'}_2025_{m:02d}.tif")
        with rasterio.open(f) as d:
            a = d.read((d.descriptions.index(var) + 1) if var != 'LAI' else 1).astype(float)[::-1]
        out.append(np.where(SC, a, np.nan))
    return out
NOTE_MET = 'Células em branco: sem dado do ERA5-Land (faixa costeira); no modelo, recebem a meteorologia da célula válida mais próxima (Seção 2.4).'
PANELS = [('figA2', 'LAI', lambda a: a, 'IAF (m² m⁻²)', 'YlGn',
           'IAF mensal usado na execução. Em 2025, IAF do MODIS corrigido com o VIIRS: média de 2013–2023 do MODIS na célula para o mês, multiplicada pela razão de anomalia do VIIRS (Suomi NPP) da zona da célula (Seção 2.3.3 do relatório). IAF nulo nas células sem vegetação.'),
          ('figA3', 'TEMP2', lambda a: a - 273.15, 'Temperatura (°C)', 'RdYlBu_r', NOTE_MET),
          ('figA4', 'PAR', lambda a: a, 'PAR (W m⁻²)', 'magma', 'PAR = 0,45 × radiação de onda curta descendente, com o passo das 00 UTC desacumulado. ' + NOTE_MET),
          ('figA5', 'SOILM', lambda a: a, 'Umidade do solo (m³ m⁻³)', 'BrBG', NOTE_MET),
          ('figA6', 'WINDSPD', lambda a: a, 'Vento (m s⁻¹)', 'viridis', NOTE_MET),
          ('figA7', 'QV', lambda a: a * 1e3, 'Razão de mistura (g kg⁻¹)', 'GnBu', 'Derivada do ponto de orvalho a 2 m e da pressão à superfície. ' + NOTE_MET),
          ('figA8', 'PRES', lambda a: a / 100, 'Pressão (hPa)', 'cividis', NOTE_MET)]
ND = {'LAI': 2, 'TEMP2': 1, 'PAR': 1, 'SOILM': 2, 'WINDSPD': 2, 'QV': 1, 'PRES': 1}
for name, var, conv, lbl, cmap_name, note in PANELS:
    fields = [conv(a) for a in monthly(var)]
    allv = np.concatenate([f[SC] for f in fields]); allv = allv[np.isfinite(allv)]
    lo, hi = np.percentile(allv, 0.5), np.percentile(allv, 99.5)
    if var == 'LAI': lo = 0
    b, dd = ebreaks(lo, hi); b[0] = min(b[0], np.nanmin(allv)); b[-1] = max(b[-1], np.nanmax(allv))   # open first and last classes
    fig = plt.figure(figsize=(8.6, 5.9), dpi=300); gs = fig.add_gridspec(4, 4, height_ratios=[1, 1, 1, 0.6], hspace=0.02, wspace=0.03)
    for i, (f, mes) in enumerate(zip(fields, MES)):
        ax = fig.add_subplot(gs[i // 4, i % 4]); W.extent(ax, axis_off=True); gridmap(ax, f, b, cmap_name, lw=0.5)
        ax.text(-53.9, -29.4, f"{mes}\nmédia {W.num(np.nanmean(f), ND[var])}", fontsize=6.4, fontweight='bold', va='bottom', ha='left')
    ax = fig.add_subplot(gs[3, 0:2]); ax.axis('off')
    cols = W.ramp(cmap_name, len(b) - 1)
    hs = [Patch(facecolor=cols[i], edgecolor='0.45', lw=0.4, label=lb) for i, lb in enumerate(labels(b, 'pt', lo_open=True, hi_open=True, d=dd))]
    ax.legend(handles=hs, title=lbl, loc='upper left', ncol=3, frameon=False, fontsize=6.0, title_fontsize=6.8, handlelength=1.4, alignment='left')
    ax = fig.add_subplot(gs[3, 3]); W.extent(ax, axis_off=True); ax.set_ylim(-29.5, -27.0)   # bottom right: arrow and bar
    W.north_arrow(ax, x=-48.4, y0=-28.6, y1=-27.4, fs=8); W.scale_bar(ax, right=-48.2, lat_sb=-29.2, small=True)
    note_below(fig, note + (' Médias de todos os passos de 3 h do mês, inclusive os noturnos.' if var != 'LAI' else '') +
               ' Valores nos painéis: média das células de SC com dado. ' + NOTE_GRID + '.', y=0.075, width=175)
    W.save(fig, os.path.join(OUTR, f'{name}.png'))
    LIM.append(dict(figure=name, panel=var, limits=[round(float(x), 4) for x in b]))
    print(name, var, 'classes', [round(float(x), 3) for x in b], 'monthly SC means', [round(float(np.nanmean(f)), 2) for f in fields])

# ---- Figure A.9: BVOC intensity by year (4 x 6 grid; the 24th panel holds the legend and the cartographic elements)
Y = list(range(2003, 2026)); maps = {}
for y in Y:
    g_ = grid_year(y); bv = sum(g_[k] for k in c.BVOC); maps[y] = (np.where(SC, bv * 1e3 / (A / 1e4), np.nan), np.nansum(bv * SC) / 1e3)
q = qbreaks(np.concatenate([m[0][SC] for m in maps.values()]))
fig, axs = plt.subplots(4, 6, figsize=(8.6, 4.9), dpi=300)
for ax in axs.ravel(): W.extent(ax, axis_off=True)
for ax, y in zip(axs.ravel(), Y):
    gridmap(ax, maps[y][0], q, 'YlGnBu', lw=0.4)
    ax.text(-53.9, -29.4, f"{y}\n{W.num(maps[y][1], 0)} Gg", fontsize=5.6, fontweight='bold', va='bottom', ha='left')
ax = axs[3, 5]; legend_grid(ax, q, 'YlGnBu', 'COVB total\n(kg ha⁻¹ ano⁻¹)', fs=5.2, tfs=5.8, anchor=(-0.05, 0.25))
W.north_arrow(ax, x=-48.5, y0=-28.6, y1=-27.7, fs=6); W.scale_bar(ax, right=-48.1, lat_sb=-29.35, small=True)
fig.tight_layout(w_pad=0.1, h_pad=0.1)
note_below(fig, 'MEGANv3.21; intensidade por célula da grade de 0,025°; classes por quintis das células de SC, '
                'comuns a todos os anos; total estadual de cada ano em cada painel. ' + NOTE_GRID + '.', y=0.0, width=175)
W.save(fig, os.path.join(OUTR, 'figA9.png'))
LIM.append(dict(figure='figA9', panel='BVOC', limits=[round(float(x), 4) for x in q]))
pd.DataFrame(LIM).to_csv(os.path.join(FIG, 'maps_class_limits.csv'), index=False)
print('report maps ->', OUTR)
