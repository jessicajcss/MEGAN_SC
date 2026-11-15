#!/usr/bin/env python3
"""study_area_map.py - study-area figure of the article (Figure 1): position of Santa Catarina in South America, land use
and land cover of MapBiomas Collection 10 in 2019 (the year compared with the previous state inventory and the year of the
native-forest shares of article Table 2 and report Table 3), and the MEGAN growth forms that the model derives from it on the 0.025° grid.

Panel a follows tools/study_area_map.py of the biomass-burning pipeline (same data, classes and colours of the MapBiomas
Collection 10 legend); panel b shows the tree fraction of each grid cell in 2019 (data/run_v10/2019/inputs/growthform_2019.npz)
and the state shares of the four growth forms.

Inputs: data/external/mapbiomas_c10_lulc2019_SC_120m.npz and mapbiomas_c10_lulc2019_SC_class_areas.csv (MapBiomas Collection
10, 2019, Santa Catarina window at about 120 m and class areas inside the state), data/external/ibge_br_uf_2024_simplified.geojson
and ibge_america_do_sul_2021_simplified.geojson (IBGE), data/boundaries (municipalities and mesoregions of 2024), and the
2019 growth forms of the production run.
Outputs: <MEGAN_SC_OUT>/figures/{pt,en}/fig_study_area.png and <MEGAN_SC_OUT>/figures/study_area_shares.csv (published
copies in results/v10/figures)."""
import math, os, sys
import numpy as np, pandas as pd, geopandas as gpd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); import _paths as P
sys.path.insert(0, HERE); import wf_style as W
REPO = P.REPO; EXT = os.path.join(REPO, 'data', 'external'); BND = P.BND
YEAR = 2019; LAT0 = -27.6
GROUPS = [  # key, MapBiomas class ids, colour (Collection 10 legend), PT label, EN label
    ('forest', [3, 5, 49], '#1f8d49', 'Formação florestal', 'Forest formation'),
    ('grass', [12], '#d6bc74', 'Formação campestre', 'Grassland'),
    ('wet', [11, 29, 32, 50], '#519799', 'Outras formações não florestais', 'Other non-forest formations'),
    ('silv', [9], '#7a5900', 'Silvicultura', 'Forest plantation'),
    ('past', [15], '#edde8e', 'Pastagem', 'Pasture'),
    ('crop', [39, 40, 41], '#f54ca9', 'Lavouras temporárias', 'Temporary crops'),
    ('mosaic', [21], '#ffefc3', 'Mosaico de usos', 'Mosaic of uses'),
    ('urban', [24], '#d4271e', 'Área urbanizada', 'Urban area'),
    ('nonveg', [23, 25, 30], '#ffa07a', 'Outras áreas não vegetadas', 'Other non-vegetated areas'),
    ('water', [33, 31], '#2532e4', 'Água', 'Water'),
]
GF = [('tree', 'Árvore', 'Tree', '#1f8d49'), ('herb', 'Herbácea', 'Herbaceous', '#d6bc74'), ('crop', 'Cultivo', 'Crop', '#f54ca9'),
      ('shrub', 'Arbusto', 'Shrub', '#8fbc8f'), ('none', 'Sem vegetação', 'No vegetation', '#bdbdbd')]
TXT = {
    'pt': dict(a=f'a) Uso e cobertura da terra (MapBiomas), {YEAR}', legend='Classes (% da área classificada\ndo Estado)', sa='América do Sul', ocean='Oceano\nAtlântico',
               capital='Capital', city='Município (sede aproximada)', meso='Limite de mesorregião', BR='BRASIL', AR='ARGENTINA',
               b='b) Fração de árvore na grade do modelo (0,025°), 2019', cb='Fração da célula ocupada\npela forma árvore',
               shares='Formas de crescimento do MEGAN\nno Estado, 2019 (% da área)',
               note=f'MapBiomas, Coleção 10 ({YEAR}); formas de crescimento do MEGAN obtidas do MapBiomas a 30 m e agregadas à grade de 0,025°. '
                    'Limites: IBGE (municípios e UF de 2024; América do Sul de 2021). Coordenadas geográficas, SIRGAS 2000/WGS 84.'),
    'en': dict(a=f'a) Land use and land cover (MapBiomas), {YEAR}', legend='Classes (% of the classified\nstate area)', sa='South America', ocean='Atlantic\nOcean',
               capital='State capital', city='Municipality (approximate seat)', meso='Mesoregion boundary', BR='BRAZIL', AR='ARGENTINA',
               b='b) Tree fraction on the model grid (0.025°), 2019', cb='Fraction of the cell covered\nby the tree growth form',
               shares='MEGAN growth forms in the\nstate, 2019 (% of the area)',
               note=f'MapBiomas, Collection 10 ({YEAR}); MEGAN growth forms derived from MapBiomas at 30 m and aggregated to the 0.025° grid. '
                    'Boundaries: IBGE (municipalities and states, 2024; South America, 2021). Geographic coordinates, SIRGAS 2000/WGS 84.'),
}
CITIES = ['Florianópolis', 'Joinville', 'Blumenau', 'Chapecó', 'Lages', 'Criciúma', 'Caçador', 'São Joaquim']
OCEAN, LAND, EDGE = '#dfe9f3', '#ebebeb', '#ffffff'
halo = [pe.withStroke(linewidth=1.6, foreground='white')]


def lulc():
    z = np.load(os.path.join(EXT, 'mapbiomas_c10_lulc2019_SC_120m.npz'))
    a, inside, t = z['a'], z['inside'], z['transform']
    lut = np.full(256, -1, dtype=int)
    for k, (_, ids, *_r) in enumerate(GROUPS):
        lut[ids] = k
    g = lut[a].astype(float); g[(~inside) | (g < 0)] = np.nan
    extent = (t[2], t[2] + t[0] * a.shape[1], t[5] + t[4] * a.shape[0], t[5])
    areas = pd.read_csv(os.path.join(EXT, 'mapbiomas_c10_lulc2019_SC_class_areas.csv'))
    pct = {key: areas[areas.class_id.isin(ids)].pct_classified.sum() for key, ids, *_r in GROUPS}
    assert abs(sum(pct.values()) - 100) < 0.01, sum(pct.values())
    return g, extent, pct


def growth_forms():
    sys.path.insert(0, os.path.join(REPO, 'analysis', 'v10')); import common as c
    lat, lon = c.grid(); SC = c.sc_mask(lat, lon); A = c.cell_area_m2(lat, lon)
    z = np.load(os.path.join(P.RUN, str(YEAR), 'inputs', f'growthform_{YEAR}.npz'))
    tot = (A * SC).sum(); sh = {k: 100 * float((z[k] * A * SC).sum() / tot) for k in ('tree', 'herb', 'crop', 'shrub')}
    sh['none'] = 100 - sum(sh.values())
    ext = [lon[0] - 0.0125, lon[-1] + 0.0125, lat[0] - 0.0125, lat[-1] + 0.0125]
    return np.where(SC, z['tree'], np.nan), ext, sh


def deg(lang):
    return FuncFormatter(lambda v, _: (f'{v:g}°').replace('.', ',' if lang == 'pt' else '.').replace('-', '−'))


def frame(ax, lang):
    ax.set_xlim(-54.0, -48.2); ax.set_ylim(-29.45, -25.8); ax.set_aspect(1 / math.cos(math.radians(LAT0)))
    ax.xaxis.set_major_formatter(deg(lang)); ax.yaxis.set_major_formatter(deg(lang))
    ax.set_xticks(np.arange(-54, -48.5, 1)); ax.set_yticks(np.arange(-29, -25.9, 1))
    ax.tick_params(labelsize=6.5, direction='in', top=True, right=True, length=3)
    ax.grid(color='0.55', lw=0.25, ls=':', zorder=1)


def main(lang, g, extent, pct, tree, gext, sh, mun, state, mesog, uf, sa):
    L = TXT[lang]; num = (lambda v: f'{v:.1f}'.replace('.', ',')) if lang == 'pt' else (lambda v: f'{v:.1f}')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7})
    fig = plt.figure(figsize=(7.2, 8.6), dpi=400)
    # ---- a) land use and land cover (top half)
    ax = fig.add_axes([0.305, 0.555, 0.685, 0.43]); ax.set_facecolor(OCEAN)
    sa.plot(ax=ax, color=LAND, edgecolor=EDGE, lw=0.6); state.plot(ax=ax, color=OCEAN, lw=0, zorder=1.5)
    uf[uf.SIGLA_UF != 'SC'].plot(ax=ax, color=LAND, edgecolor='#bdbdbd', lw=0.5)
    ax.imshow(g, extent=extent, cmap=ListedColormap([c for _, _, c, *_r in GROUPS]), vmin=-0.5, vmax=len(GROUPS) - 0.5,
              interpolation='nearest', zorder=2)
    mesog.boundary.plot(ax=ax, color='0.15', lw=0.45, zorder=3); state.boundary.plot(ax=ax, color='0.05', lw=0.8, zorder=4)
    pts = mun.set_index('NM_MUN').loc[CITIES].geometry.representative_point()
    for name, p in pts.items():
        cap = name == 'Florianópolis'
        ax.plot(p.x, p.y, marker='*' if cap else 'o', ms=7 if cap else 3.2, mfc='#ffd400' if cap else 'k', mec='k', mew=0.5, zorder=7)
        dx, dy, ha = (0.07, 0.02, 'left') if name not in ('Florianópolis', 'Joinville') else (-0.07, 0.02, 'right')
        if cap: dx, dy, ha = (-0.08, -0.13, 'right')
        ax.text(p.x + dx, p.y + dy, name, fontsize=6.3, ha=ha, va='bottom', zorder=7, path_effects=halo)
    for lab, (x, y) in (('PR', (-51.3, -25.95)), ('RS', (-52.3, -29.25)), (L['AR'], (-53.93, -27.65))):
        ax.text(x, y, lab, fontsize=6.5, color='0.45', ha='left', va='center', zorder=5)
    ax.text(-48.62, -28.95, L['ocean'], fontsize=6.5, style='italic', color='#3a6a96', ha='center', va='center', zorder=5)
    frame(ax, lang); ax.set_ylabel('Latitude (°)', fontsize=7, labelpad=2)
    ax.text(0.0, 1.015, L['a'], transform=ax.transAxes, fontsize=7.2, fontweight='bold', ha='left', va='bottom')
    W.north_arrow(ax, x=-48.45, y0=-26.55, y1=-26.15, fs=8); W.scale_bar(ax, right=-48.35, lat_sb=-29.36, fs=5.8)
    axi = fig.add_axes([0.015, 0.775, 0.255, 0.21]); axi.set_facecolor(OCEAN)          # South America
    sa.plot(ax=axi, color='#d9d9d9', edgecolor='white', lw=0.3); uf.plot(ax=axi, color='#bdbdbd', edgecolor='white', lw=0.25)
    uf[uf.SIGLA_UF == 'SC'].plot(ax=axi, color='#c0392b', edgecolor='#c0392b', lw=0.4)
    axi.add_patch(Rectangle((-54.0, -29.45), 5.8, 3.65, fill=False, edgecolor='k', lw=0.6))
    axi.set_xlim(-82, -34); axi.set_ylim(-56, 13.5); axi.set_aspect(1 / math.cos(math.radians(-20))); axi.set_xticks([]); axi.set_yticks([])
    axi.text(0.04, 0.04, L['sa'], transform=axi.transAxes, fontsize=6.3, va='bottom', ha='left'); axi.text(-52.5, -9.0, L['BR'], fontsize=6, color='0.3', ha='center')
    axl = fig.add_axes([0.015, 0.545, 0.27, 0.215]); axl.axis('off')
    order = sorted(range(len(GROUPS)), key=lambda k: -pct[GROUPS[k][0]])
    hs = [Patch(facecolor=GROUPS[k][2], edgecolor='0.4', lw=0.3, label=f'{GROUPS[k][3] if lang == "pt" else GROUPS[k][4]} ({num(pct[GROUPS[k][0]])}%)') for k in order]
    hs += [Line2D([], [], color='0.15', lw=0.6, label=L['meso']),
           Line2D([], [], marker='*', ls='', ms=7, mfc='#ffd400', mec='k', mew=0.5, label=L['capital']),
           Line2D([], [], marker='o', ls='', ms=3.2, mfc='k', mec='k', label=L['city'])]
    axl.legend(handles=hs, loc='upper left', frameon=False, fontsize=6.0, title=L['legend'], title_fontsize=6.6, alignment='left',
               handlelength=1.6, handleheight=1.0, labelspacing=0.36, borderaxespad=0)
    # ---- b) tree fraction on the model grid and shares of the growth forms (bottom half)
    ax = fig.add_axes([0.305, 0.085, 0.685, 0.43]); ax.set_facecolor(OCEAN)
    sa.plot(ax=ax, color=LAND, edgecolor=EDGE, lw=0.6); uf[uf.SIGLA_UF != 'SC'].plot(ax=ax, color=LAND, edgecolor='#bdbdbd', lw=0.5)
    state.plot(ax=ax, color='white', lw=0, zorder=1.5)
    bnd = np.linspace(0, 1, 11); cmap = ListedColormap(W.ramp('Greens', 10, 0.08, 0.95)); norm = BoundaryNorm(bnd, cmap.N)
    ax.imshow(np.ma.masked_invalid(tree), origin='lower', extent=gext, cmap=cmap, norm=norm, interpolation='nearest', zorder=2)
    mesog.boundary.plot(ax=ax, color='0.15', lw=0.45, zorder=3); state.boundary.plot(ax=ax, color='0.05', lw=0.8, zorder=4)
    frame(ax, lang); ax.set_xlabel('Longitude (°)', fontsize=7, labelpad=2); ax.set_ylabel('Latitude (°)', fontsize=7, labelpad=2)
    ax.text(0.0, 1.015, L['b'], transform=ax.transAxes, fontsize=7.2, fontweight='bold', ha='left', va='bottom')
    W.north_arrow(ax, x=-48.45, y0=-26.55, y1=-26.15, fs=8); W.scale_bar(ax, right=-48.35, lat_sb=-29.36, fs=5.8)
    cax = fig.add_axes([0.03, 0.40, 0.2, 0.014])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation='horizontal', ticks=[0, 0.25, 0.5, 0.75, 1])
    cb.ax.tick_params(labelsize=6); cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: (f'{v:g}').replace('.', ',' if lang == 'pt' else '.')))
    cb.set_label(L['cb'], fontsize=6.3)
    axs = fig.add_axes([0.105, 0.115, 0.125, 0.17]); y = np.arange(len(GF))[::-1]
    for yi, (k, pt_, en_, col) in zip(y, GF):
        axs.barh(yi, sh[k], color=col, edgecolor='0.4', lw=0.3, height=0.7)
        axs.text(sh[k] + 1.2, yi, f'{num(sh[k])}%', va='center', fontsize=6)
    axs.set_yticks(y); axs.set_yticklabels([pt_ if lang == 'pt' else en_ for _, pt_, en_, _ in GF], fontsize=6.2)
    axs.set_xlim(0, 75); axs.set_xticks([]); [axs.spines[s].set_visible(False) for s in ('top', 'right', 'bottom')]
    fig.text(0.02, 0.30, L['shares'], fontsize=6.3, ha='left', va='bottom')
    fig.text(0.305, 0.012, L['note'], fontsize=5.6, ha='left', va='bottom', wrap=True)
    out = P.out('figures', lang, 'fig_study_area.png'); fig.savefig(out, dpi=400, facecolor='white'); plt.close(fig)
    print('written', out)


if __name__ == '__main__':
    mun = gpd.read_file(os.path.join(BND, 'SC_Municipios_2024.gpkg')).to_crs(4674)
    mesog = gpd.read_file(os.path.join(BND, 'sc_mesoregions_2024.gpkg')).to_crs(4674)
    state = W.state_from(mun).to_crs(4674)
    uf = gpd.read_file(os.path.join(EXT, 'ibge_br_uf_2024_simplified.geojson')).to_crs(4674)
    sa = gpd.read_file(os.path.join(EXT, 'ibge_america_do_sul_2021_simplified.geojson')).to_crs(4674)
    g, extent, pct = lulc(); tree, gext, sh = growth_forms()
    for lang in ('pt', 'en'):
        main(lang, g, extent, pct, tree, gext, sh, mun, state, mesog, uf, sa)
    rows = [dict(source='MapBiomas C10 2019, % of the classified state area', group=k, pct=round(pct[k], 3)) for k, *_ in GROUPS]
    rows += [dict(source='MEGAN growth forms 2019, % of the state grid area', group=k, pct=round(sh[k], 3)) for k, *_ in GF]
    pd.DataFrame(rows).to_csv(P.out('figures', 'study_area_shares.csv'), index=False)
    print(pd.DataFrame(rows).to_string(index=False))
