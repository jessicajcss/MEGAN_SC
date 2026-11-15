#!/usr/bin/env python3
"""ijest_figures.py - Figures 1-11 and S4-S5 of the English manuscript prepared for the International Journal of Environmental
Science and Technology (Springer), drawn from the result files of this repository.

Same data and the same plotting choices as the published figures of the article (notebook 13 for the charts,
maps/12_maps_wf_style.py and maps/study_area_map.py for the maps, schematics/verification_flowchart.py for the
verification figure and fig9_tropomi.py for the TROPOMI comparison), redrawn for the journal's artwork rules:
  * no figure titles, no source notes and no second heading line (period, units, statistics) inside the figures; the
    captions in the manuscript hold them;
  * final width 174 mm (Fig. 5: 84 mm), height at most 234 mm; sans-serif lettering (Liberation Sans, metrically equal to
    Arial) of about 6.5 to 9 pt at the final size; PNG at 600 dpi;
  * American English and the mesoregion names used in the manuscript (Table 2);
  * the TROPOMI figure of the article split into a spatial (Fig. 10) and a temporal (Fig. 11) figure;
  * emissions in Gg (1 Gg = 1 kt); each PNG is also written as an LZW-compressed TIFF.
No number is typed in: every value is read from the outputs of notebooks 13 and 14, the benchmark and the run, and the
correlations drawn in Fig. S4 are asserted against correlations_emissions_drivers.csv.

Figure numbers (manuscript; FigS = Online Resource 1), the number of the function that draws each one, and sources:
  Fig1   fig1   study area (MapBiomas 2019, growth forms)          maps/study_area_map.py
  FigS4  fig2   detrended correlation example                       notebook 13, fig_detrend
  Fig2   fig3   verification and evaluation                         schematics/verification_flowchart.py
  Fig3   fig4   annual series                                       notebook 13, fig1
  Fig4   fig5   composition and seasonal cycles                     notebook 13, fig2
  Fig5   fig6   detrended correlations with drivers                 notebook 13, fig3
  Fig6   fig7   SOA precursors by mesoregion                        notebook 13, fig4
  Fig7   fig8   SOA precursor intensity by mesoregion (map)         maps/12_maps_wf_style.py, fig5
  FigS5  fig9   intensity versus native forest                      notebook 13, fig6
  Fig8   fig10  BVOC by municipality, 2025                          maps/12_maps_wf_style.py, fig7
  Fig9   fig11  isoprene, CO and soil NO by municipality            maps/12_maps_wf_style.py, fig8
  Fig10  fig12  TROPOMI dHCHO and modeled isoprene: maps            fig9_tropomi.py, panels a, b, c, f
  Fig11  fig13  TROPOMI dHCHO and modeled isoprene: time series     fig9_tropomi.py, panels d, e
Usage: python analysis/v10/ijest/ijest_figures.py [1 2 ... 12]   (function numbers; default: all; 12 draws fig12 and fig13)
Output: <MEGAN_SC_OUT>/figures/ijest/FigN.png and FigN.tif (published copy in results/v10/figures/ijest)
"""
import glob, math, os, shutil, sys, textwrap, types
import numpy as np, pandas as pd, geopandas as gpd, rasterio, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, FancyBboxPatch, Circle, Polygon
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, TextArea, VPacker
from matplotlib.ticker import FuncFormatter, MultipleLocator
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); V10 = os.path.dirname(HERE)
sys.path.insert(0, V10); import _paths as P
sys.path.insert(0, os.path.join(V10, 'maps')); import wf_style as W
sys.path.insert(0, os.path.join(V10, 'schematics'))

W2, W1, HMAX = 174 / 25.4, 84 / 25.4, 234 / 25.4                 # Springer widths and maximum height (inches)
DPI = 600
FONT = ['Liberation Sans', 'DejaVu Sans']                       # DejaVu Sans only for glyphs missing in Liberation Sans


def pin_fonts():
    """Fix the font files. With fonts-liberation and fonts-liberation2 both installed, matplotlib scores the two
    Liberation Sans files of each style equally and takes the first one in its font list, whose order follows Python's
    string hashing whenever the font cache is rebuilt. Sorting the list (matplotlib's own fonts first, then Liberation
    Sans 2, then by path) makes the choice independent of the cache. Same result where only one version is installed."""
    from matplotlib import font_manager as fm
    rank = lambda f: (0 if 'mpl-data' in f.fname else 1 if 'liberation2' in f.fname else 2, f.fname)
    fm.fontManager.ttflist.sort(key=rank)
    fm.FontManager._findfont_cached.cache_clear()


def style():
    plt.rcParams.update({'font.family': FONT, 'font.size': 7.5, 'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'legend.fontsize': 7, 'axes.edgecolor': '#4d4d4d', 'axes.linewidth': 0.7, 'text.color': '#1a1a1a',
                         'axes.labelcolor': '#1a1a1a', 'xtick.color': '#4d4d4d', 'ytick.color': '#4d4d4d',
                         'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'lines.linewidth': 1.2, 'axes.unicode_minus': True})


pin_fonts()
style()
OUT = P.out('figures', 'ijest'); os.makedirs(OUT, exist_ok=True)
NB13, BND, RUN = P.NB13, P.BND, P.RUN
MESO_EN = {'Oeste Catarinense': 'Western SC', 'Norte Catarinense': 'Northern SC', 'Vale do Itajaí': 'Vale do Itajaí',
           'Serrana': 'Mountain region', 'Grande Florianópolis': 'Greater Florianópolis', 'Sul Catarinense': 'Southern SC'}
PALC = dict(tot="#3d4a57", iso="#0072B2", mt="#D55E00", sqt="#E69F00", ov="#009E73", lv="#CC79A7", oth="#999999", mbo="#56B4E9", no="#C0392B")
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
YEARS = list(range(2003, 2026))
MINUS = lambda s: s.replace('-', '−')


def clean(ax):
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
def grid_y(ax): ax.yaxis.grid(True, color='#e8e8e6', lw=0.6, zorder=0); ax.set_axisbelow(True)
def grid_x(ax): ax.xaxis.grid(True, color='#e8e8e6', lw=0.6, zorder=0); ax.set_axisbelow(True)
def nice(v, step): return step * np.ceil(v * 1.08 / step)
def panel(ax, txt, x=0.0, y=1.02, fs=8):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=fs, fontweight='bold', ha='left', va='bottom')


def ptitle(ax, title, sub=None, fs=8.6, fs2=7.0, pad=3.0, x=0.0, xref='axes fraction'):
    """Panel heading in bold, left-aligned above the axes. x is in axes fraction, or in figure fraction with
    xref='figure fraction'. The second line `sub` (units, period, statistics) is drawn only when SUBTITLES is True, which it
    is not for the journal figures; its content is in the figure caption."""
    off = pad
    if sub and SUBTITLES:
        ax.annotate(sub, xy=(x, 1), xycoords=(xref, 'axes fraction'), xytext=(0, pad), textcoords='offset points',
                    ha='left', va='bottom', fontsize=fs2, color='0.25', annotation_clip=False)
        off = pad + fs2 * 1.38
    ax.annotate(title, xy=(x, 1), xycoords=(xref, 'axes fraction'), xytext=(0, off), textcoords='offset points',
                ha='left', va='bottom', fontsize=fs, fontweight='bold', annotation_clip=False)


def north(ax, x, y, w=7.5, h=11.0, fs=8.5):
    """North arrow of fixed size in points: bold N above an arrowhead split in a black and a white half, with the top of the
    N at (x, y) in data coordinates."""
    da = DrawingArea(w, h, 0, 0)
    da.add_artist(Polygon([(w / 2, h), (0, 0), (w / 2, h * 0.32)], closed=True, fc='black', ec='black', lw=0.6))
    da.add_artist(Polygon([(w / 2, h), (w, 0), (w / 2, h * 0.32)], closed=True, fc='white', ec='black', lw=0.6))
    box = VPacker(children=[TextArea('N', textprops=dict(fontsize=fs, fontweight='bold')), da], align='center', pad=0, sep=0.8)
    ab = AnchoredOffsetbox(loc='upper center', child=box, pad=0, borderpad=0, frameon=False, bbox_to_anchor=(x, y),
                           bbox_transform=ax.transData)
    ab.set_zorder(10); ax.add_artist(ab)


SUBTITLES = False
NEW = {1: '1', 2: 'S4', 3: '2', 4: '3', 5: '4', 6: '5', 7: '6', 8: '7', 9: 'S5', 10: '8', 11: '9', 12: '10', 13: '11'}


def save(fig, n, tight=True):
    n = NEW[n]
    f = os.path.join(OUT, f'Fig{n}.png')
    if tight: fig.savefig(f, dpi=DPI, facecolor='white', bbox_inches='tight', pad_inches=0.04)
    else: fig.savefig(f, dpi=DPI, facecolor='white')
    plt.close(fig)
    from PIL import Image
    with Image.open(f) as im: wpx, hpx = im.size
    w_in, h_in = wpx / DPI, hpx / DPI
    print(f'Fig{n}: {wpx} x {hpx} px = {w_in * 25.4:.0f} x {h_in * 25.4:.0f} mm at {DPI} dpi -> {f}')
    with Image.open(f) as im:
        im.convert('RGB').save(f[:-4] + '.tif', compression='tiff_lzw', dpi=(DPI, DPI))
    assert h_in <= HMAX + 0.05 and w_in <= W2 + 0.12, ('size', n, w_in, h_in)


# ---------------------------------------------------------------------------- tables written by notebook 13
G = pd.read_csv(os.path.join(NB13, 'table2_annual_totals_SC_groups.csv')).set_index('year')
COR = pd.read_csv(os.path.join(NB13, 'correlations_emissions_drivers.csv'))
DRV = pd.read_csv(os.path.join(NB13, 'drivers_annual_SC.csv')).set_index('year')
SEAS = pd.read_csv(os.path.join(NB13, 'seasonal_climatology.csv'), header=[0, 1], index_col=0)
CLIM, CMIN, CMAX = SEAS['mean'], SEAS['min'], SEAS['max']
T3 = pd.read_csv(os.path.join(NB13, 'table3_mesoregions_soa.csv')).set_index('mesoregion')
CLASSES = ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM', 'NO', 'SQT_HR', 'SQT_LR',
           'MEOH', 'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO']


def fig2():
    """Detrended correlation, soil NO and air temperature (notebook 13, fig_detrend)."""
    yr = np.array(YEARS); x0 = np.arange(len(yr)); tt = DRV['TEMP2_C'].loc[yr].values.astype(float); nn = (G.NO.loc[yr] / 1e3).values.astype(float)
    fit_t = np.polyval(np.polyfit(x0, tt, 1), x0); fit_n = np.polyval(np.polyfit(x0, nn, 1), x0); rt, rn = tt - fit_t, nn - fit_n
    r_raw, r_det = stats.pearsonr(tt, nn)[0], stats.pearsonr(rt, rn)[0]
    c0 = COR[(COR.series == 'NO') & (COR.driver == 'TEMP2_C')].iloc[0]
    assert abs(r_raw - c0.r) < 1e-9 and abs(r_det - c0.r_detrended) < 1e-9
    cm0 = plt.get_cmap('viridis'); nr0 = mpl.colors.Normalize(yr.min(), yr.max()); col0 = cm0(nr0(yr))
    fig, axs = plt.subplots(2, 2, figsize=(W2, 5.7)); fig.subplots_adjust(left=0.08, right=0.875, top=0.905, bottom=0.085, hspace=0.62, wspace=0.32)
    for ax, s, f, lab, ttl, sub in ((axs[0, 0], tt, fit_t, 'Temperature (°C)', '(a) Mean 2-m air temperature', 'Statewide annual mean, ERA5-Land'),
                                    (axs[0, 1], nn, fit_n, 'Soil NO (Gg yr⁻¹)', '(b) Soil NO', 'Statewide annual total')):
        ax.vlines(yr, f, s, color='#9a9a9a', lw=0.7, zorder=2); ax.plot(yr, f, color='#3d4a57', lw=1.1, ls='--', zorder=3)
        ax.scatter(yr, s, c=col0, s=18, edgecolor='white', linewidth=0.4, zorder=4)
        ax.set_ylabel(lab); ax.set_xticks(np.arange(2003, 2026, 4)); clean(ax); grid_y(ax); ptitle(ax, ttl, sub)
    hh = [Line2D([], [], color='#3d4a57', lw=1.1, ls='--'), Line2D([], [], color='#9a9a9a', lw=0.7),
          Line2D([], [], marker='o', ls='', color=cm0(0.5), ms=4, markeredgecolor='white')]
    axs[0, 0].legend(hh, ['least-squares line', 'deviation from the line', 'annual value (color: year)'], fontsize=6.5, loc='upper left', frameon=False)
    for ax, xv, yv_, xl, yl, ttl, sub, r in ((axs[1, 0], tt, nn, 'Temperature (°C)', 'Soil NO (Gg yr⁻¹)', '(c) Original series', 'Soil NO versus temperature', r_raw),
                                              (axs[1, 1], rt, rn, 'Temperature deviation (°C)', 'Soil NO deviation (Gg yr⁻¹)', '(d) Detrended residuals',
                                               'Deviations from the least-squares lines', r_det)):
        b1, b0 = np.polyfit(xv, yv_, 1); xx = np.linspace(xv.min(), xv.max(), 50)
        ax.plot(xx, b1 * xx + b0, color='#C0392B', lw=1.0, zorder=3); ax.scatter(xv, yv_, c=col0, s=18, edgecolor='white', linewidth=0.4, zorder=4)
        ax.set_xlabel(xl); ax.set_ylabel(yl); clean(ax); ptitle(ax, ttl, sub)
        ax.text(0.03, 0.97, f'r = {r:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=8, fontweight='bold')
    axs[1, 1].axhline(0, color='0.78', lw=0.6, zorder=1); axs[1, 1].axvline(0, color='0.78', lw=0.6, zorder=1)
    cax = fig.add_axes([0.905, 0.12, 0.014, 0.76]); cb0 = fig.colorbar(mpl.cm.ScalarMappable(norm=nr0, cmap=cm0), cax=cax)
    cb0.set_ticks(np.arange(2004, 2026, 4)); cb0.set_label('Year'); cb0.ax.tick_params(labelsize=7)
    save(fig, 2, tight=False)


def fig4():
    """Annual series (notebook 13, fig1); only panels b and c share the vertical scale."""
    yr = np.array(YEARS); common_ax = nice(max(G.ISOP.max(), G.MT.max()) / 1e3, 50)
    panels = [('a', 'Total BVOCs', G.BVOC / 1e3, PALC['tot'], nice(G.BVOC.max() / 1e3, 100), 200),
              ('b', 'Isoprene', G.ISOP / 1e3, PALC['iso'], common_ax, 50), ('c', 'Monoterpenes', G.MT / 1e3, PALC['mt'], common_ax, 50),
              ('d', 'Sesquiterpenes', G.SQT / 1e3, PALC['sqt'], nice(G.SQT.max() / 1e3, 5), 5), ('e', 'Soil NO', G.NO / 1e3, PALC['no'], nice(G.NO.max() / 1e3, 1), 3)]
    fig, axes = plt.subplots(5, 1, figsize=(W2, 8.0), sharex=True); fig.subplots_adjust(left=0.09, right=0.975, top=0.965, bottom=0.045, hspace=0.42)
    for ax, (lab, name, s, col, ymax, ymaj) in zip(axes, panels):
        ax.bar(yr, s.loc[yr].values, width=0.72, color=col, edgecolor='white', linewidth=0.3, zorder=3)
        ax.set_ylim(0, ymax); ax.yaxis.set_major_locator(MultipleLocator(ymaj)); ax.set_ylabel('Gg yr⁻¹')
        clean(ax); grid_y(ax); ptitle(ax, f'({lab}) {name}')
    axes[-1].set_xlim(yr.min() - 0.7, yr.max() + 0.7); axes[-1].set_xticks(np.arange(YEARS[0], YEARS[-1] + 1, 2))
    save(fig, 4, tight=False)


def fig5():
    """Composition (a) and seasonal cycles (b, c) (notebook 13, fig2)."""
    fam = [('Isoprene', G.ISOP.mean(), PALC['iso']), ('Monoterpenes', G.MT.mean(), PALC['mt']), ('Oxygenated VOCs', G.OVOC.mean(), PALC['ov']),
           ('Light hydrocarbons (LVOC)', G.LVOC.mean(), PALC['lv']), ('Sesquiterpenes', G.SQT.mean(), PALC['sqt']),
           ('Other/stress', G.STRESS_OTHER.mean(), PALC['oth']), ('MBO', G.MBO.mean(), PALC['mbo'])]
    fam = sorted([(n, v / 1e3, c) for n, v, c in fam], key=lambda t: t[1], reverse=True); total = sum(v for _, v, _ in fam)
    assert abs(total - G.BVOC.mean() / 1e3) < 1e-6
    fig = plt.figure(figsize=(W2, 7.2))
    gsA = fig.add_gridspec(1, 1, left=0.25, right=0.97, top=0.925, bottom=0.615)
    gsB = fig.add_gridspec(1, 2, left=0.08, right=0.985, top=0.475, bottom=0.165, wspace=0.27)
    ax = fig.add_subplot(gsA[0, 0]); y_ = np.arange(len(fam))[::-1]
    for yi, (n, v, c) in zip(y_, fam):
        ax.barh(yi, v, color=c, edgecolor='white', height=0.72, zorder=3)
        ax.text(v + total * 0.006, yi, (f'{v:.0f} Gg ({100 * v / total:.1f}%)' if v >= 1 else f'{v:.1f} Gg ({100 * v / total:.1f}%)'), va='center', ha='left', fontsize=7)
    ax.set_yticks(y_); ax.set_yticklabels([n for n, _, _ in fam]); ax.set_xlim(0, nice(fam[0][1] * 1.25, 50)); ax.set_xlabel('Mean emission 2003–2025 (Gg yr⁻¹)')
    grid_x(ax); clean(ax)
    ptitle(ax, '(a) Mean annual composition, 2003–2025', f'Total BVOCs {total:.0f} Gg yr⁻¹; share of each family in parentheses', x=0.08, xref='figure fraction')
    mm_ = np.arange(1, 13); HL = {}
    ax = fig.add_subplot(gsB[0, 0])
    ax.fill_between(mm_, CMIN.BVOC.values / 1e3, CMAX.BVOC.values / 1e3, color=PALC['tot'], alpha=0.15, lw=0, label='Total BVOCs: minimum–maximum')
    for k, col, nm in (('BVOC', PALC['tot'], 'Total BVOCs'), ('ISOP', PALC['iso'], 'Isoprene'), ('MT', PALC['mt'], 'Monoterpenes'),
                       ('OVOC', PALC['ov'], 'Oxygenated VOCs'), ('LVOC', PALC['lv'], 'Light hydrocarbons (LVOC)')):
        ax.plot(mm_, CLIM[k].values / 1e3, 'o-', ms=2.6, lw=1.2, color=col, label=nm)
    ax.set_xticks(mm_); ax.set_xticklabels(MON, fontsize=6.5); ax.set_ylabel('Monthly emission (Gg month⁻¹)'); ax.set_ylim(0, None); clean(ax); grid_y(ax)
    ptitle(ax, '(b) Seasonal cycle', 'Monthly emissions, 2003–2025 mean; shading: range of total BVOCs'); HL.update({l_: h_ for h_, l_ in zip(*ax.get_legend_handles_labels())})
    ax = fig.add_subplot(gsB[0, 1])
    for k, col, nm in (('ISOP', PALC['iso'], 'Isoprene'), ('SQT', PALC['sqt'], 'Sesquiterpenes'), ('CO', '#8c6d31', 'CO'), ('MT', PALC['mt'], 'Monoterpenes'),
                       ('LVOC', PALC['lv'], 'Light hydrocarbons (LVOC)'), ('NO', PALC['no'], 'Soil NO')):
        ax.plot(mm_, (CLIM[k] / CLIM[k].mean()).values, marker='s' if k == 'CO' else 'o', ms=2.4, lw=1.1, color=col, label=nm, ls='--' if k == 'NO' else '-')
    ax.axhline(1, color='0.6', lw=0.7, label='monthly mean (= 1)'); ax.set_xticks(mm_); ax.set_xticklabels(MON, fontsize=6.5)
    ax.set_ylabel('Emission / monthly mean'); clean(ax); grid_y(ax); ax.yaxis.set_major_locator(MultipleLocator(0.5))
    for h_, l_ in zip(*ax.get_legend_handles_labels()): HL.setdefault(l_, h_)
    ptitle(ax, '(c) Normalized seasonal cycles', 'Each series divided by its mean monthly emission')
    order = ['Total BVOCs', 'Total BVOCs: minimum–maximum', 'Isoprene', 'Monoterpenes', 'Sesquiterpenes', 'Oxygenated VOCs',
             'Light hydrocarbons (LVOC)', 'CO', 'Soil NO', 'monthly mean (= 1)']
    fig.legend([HL[k] for k in order], order, loc='upper center', bbox_to_anchor=(0.53, 0.095), ncol=4, fontsize=6.8, frameon=False,
               columnspacing=1.2, handlelength=2.0, handletextpad=0.5)
    save(fig, 5, tight=False)


def fig6():
    """Detrended correlations between annual emissions and drivers (notebook 13, fig3); 84 mm wide."""
    Rm = np.array([[COR[(COR.series == k) & (COR.driver == d)].r_detrended.iloc[0] for d in ('TEMP2_C', 'PAR', 'LAI', 'SOILM')] for k in CLASSES])
    fig = plt.figure(figsize=(W1, 5.2)); ax = fig.add_axes([0.20, 0.135, 0.56, 0.855]); cax = fig.add_axes([0.80, 0.30, 0.04, 0.50])
    im = ax.imshow(Rm, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i in range(Rm.shape[0]):
        for j in range(Rm.shape[1]):
            ax.text(j, i, MINUS(f'{Rm[i, j]:.2f}'), ha='center', va='center', fontsize=6.5, color='white' if abs(Rm[i, j]) > 0.6 else 'black')
    ax.set_xticks(range(4)); ax.set_xticklabels(['Temperature (2 m)', 'PAR', 'LAI', 'Soil moisture'], rotation=35, ha='right', fontsize=6.8)
    ax.set_yticks(range(len(CLASSES))); ax.set_yticklabels(CLASSES, fontsize=6.6); ax.tick_params(length=2)
    cb = fig.colorbar(im, cax=cax); cb.set_label('Pearson r', fontsize=7); cb.ax.tick_params(labelsize=6.5)
    cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f'{v:g}')))
    save(fig, 6)


def fig7():
    """SOA precursors by mesoregion (notebook 13, fig4)."""
    S = T3.sort_values('SOA', ascending=False)
    fig, ax = plt.subplots(figsize=(W2, 3.0)); fig.subplots_adjust(left=0.19, right=0.975, top=0.97, bottom=0.16)
    yy = np.arange(len(S))[::-1]; mtv = S.MT.values / 1e3; sqv = S.SQT.values / 1e3
    ax.barh(yy, mtv, color=PALC['mt'], edgecolor='white', height=0.68, zorder=3, label='Monoterpenes')
    ax.barh(yy, sqv, left=mtv, color=PALC['sqt'], edgecolor='white', height=0.68, zorder=3, label='Sesquiterpenes')
    for yi, (idx, r_) in zip(yy, S.iterrows()):
        ax.text(r_.SOA / 1e3 + S.SOA.max() / 1e3 * 0.015, yi, f'{r_.pct_SC:.0f}%  |  {r_.intensity_kg_ha:.0f} kg ha⁻¹ yr⁻¹', va='center', ha='left', fontsize=7)
    ax.set_yticks(yy); ax.set_yticklabels([MESO_EN[k] for k in S.index]); ax.set_xlim(0, nice(S.SOA.max() / 1e3 * 1.45, 10))
    ax.set_xlabel('SOA precursor emissions (Gg yr⁻¹)'); grid_x(ax); clean(ax); ax.legend(loc='lower right', fontsize=7, frameon=False)
    save(fig, 7, tight=False)


def fig9():
    """SOA precursor intensity versus native forest fraction (notebook 13, fig6)."""
    x, yv = T3.forest_pct.values.astype(float), T3.intensity_kg_ha.values
    r, p = stats.pearsonr(x, yv); rho, _ = stats.spearmanr(x, yv); zz = np.arctanh(r); ci = np.tanh([zz - 1.96 / np.sqrt(len(x) - 3), zz + 1.96 / np.sqrt(len(x) - 3)])
    sl, icp = np.polyfit(x, yv, 1)
    fig, ax = plt.subplots(figsize=(W2, 3.2)); fig.subplots_adjust(left=0.09, right=0.975, top=0.97, bottom=0.15)
    xs = np.array([x.min() - 2, x.max() + 2]); ax.plot(xs, icp + sl * xs, ls='--', color='#9a9a9a', lw=1.1, zorder=2)
    ax.scatter(x, yv, s=55, color=PALC['mt'], edgecolor='white', lw=0.6, zorder=3)
    for xi, yi, idx in zip(x, yv, T3.index): ax.annotate(MESO_EN[idx], (xi, yi), xytext=(6, 3), textcoords='offset points', fontsize=7)
    ax.set_xlim(x.min() - 4, x.max() + 14); ax.set_xlabel('Native forest fraction, MapBiomas 2019 (%)'); ax.set_ylabel('SOA precursor intensity\n(kg ha⁻¹ yr⁻¹)')
    clean(ax); grid_y(ax)
    ax.text(0.02, 0.97, f'Pearson r = {r:.2f} (95% CI {ci[0]:.2f}–{ci[1]:.3f}); Spearman ρ = {rho:.2f}; n = {len(x)}', transform=ax.transAxes,
            ha='left', va='top', fontsize=7)
    save(fig, 9, tight=False)


# ---------------------------------------------------------------------------- maps
_B = {}


def bounds():
    if not _B:
        mun = gpd.read_file(os.path.join(BND, 'SC_Municipios_2024.gpkg')).to_crs(4326)
        mesog = gpd.read_file(os.path.join(BND, 'sc_mesoregions_2024.gpkg')).to_crs(4326)
        mesog['geometry'] = [W.fill_slivers(g) for g in mesog.geometry]
        _B.update(mun=mun, mesog=mesog, state=W.state_from(mun),
                  MU=pd.read_csv(os.path.join(NB13, 'table4_municipalities_2025.csv')).set_index('mun'),
                  MU19=pd.read_csv(os.path.join(NB13, 'municipalities_2019.csv')).set_index('mun'))
        for d_ in (_B['MU'], _B['MU19']): assert set(d_.index) == set(mun.NM_MUN), 'municipality names'
    return _B


def dfor(x): ax_ = abs(x); return 0 if ax_ >= 20 else (1 if ax_ >= 2 else 2)


def labels(b, d=None):
    d = dfor(b[1]) if d is None else d; out = []
    for i in range(len(b) - 1):
        lo, hi = W.num(b[i], d, 'en'), W.num(b[i + 1], d, 'en')
        out.append(f'{lo} to {hi}' if i == len(b) - 2 else f'{lo} to < {hi}')
    return out


def geo_axes(ax, xlim=None, fs=6.8, small=False):
    W.extent(ax)
    if xlim is not None: ax.set_xlim(*xlim)
    ax.xaxis.set_major_formatter(W.fmt_deg('en')); ax.yaxis.set_major_formatter(W.fmt_deg('en'))
    ax.set_xlabel('Longitude (°)', fontsize=fs + 0.4); ax.set_ylabel('Latitude (°)', fontsize=fs + 0.4); ax.tick_params(labelsize=fs)
    right = (xlim or W.XLIM)[1]
    north(ax, right - 0.35 if xlim is not None else right - 0.2, -25.86); W.scale_bar(ax, right=right - 0.12, fs=fs - 0.4, small=small)


def fig8():
    """SOA precursor intensity by mesoregion, 2003-2025 mean (maps/12_maps_wf_style.py, fig5)."""
    B = bounds(); mesog, state = B['mesog'], B['state']
    g = mesog.assign(v=mesog.name_meso.map(T3.intensity_kg_ha)); assert g.v.notna().all()
    bnds = np.arange(3 * np.floor(g.v.min() / 3), 3 * np.ceil(g.v.max() / 3) + 1, 3)
    cmap = ListedColormap(W.ramp('YlOrRd', len(bnds) - 1)); norm = BoundaryNorm(bnds, cmap.N)
    OCEAN_LABEL = {'Grande Florianópolis': (-47.62, -27.22), 'Sul Catarinense': (-47.95, -28.95)}
    fig, ax = plt.subplots(figsize=(W2, 5.1))
    g.plot(ax=ax, column='v', cmap=cmap, norm=norm, edgecolor='0.25', lw=0.5); state.boundary.plot(ax=ax, color='0.2', lw=0.6)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax, shrink=0.85, pad=0.02, format=FuncFormatter(lambda v, _: W.num(v, 0, 'en')))
    cb.set_label('SOA precursor intensity (kg ha⁻¹ yr⁻¹)', fontsize=7.5); cb.ax.tick_params(labelsize=7)
    for r in g.itertuples():
        x, y = W.label_point(r.geometry); dark = norm(r.v) >= len(bnds) - 2
        nm = MESO_EN[r.name_meso].replace('Greater ', 'Greater\n').replace(' region', '\nregion')
        txt = f'{nm}\n{W.num(r.v, 1, "en")}'
        if r.name_meso in OCEAN_LABEL:
            tx, ty = OCEAN_LABEL[r.name_meso]
            ax.annotate(txt, xy=(x, y), xytext=(tx, ty), ha='center', va='center', fontsize=7.2, color='black', linespacing=1.1,
                        arrowprops=dict(arrowstyle='-', color='0.25', lw=0.6, shrinkA=0, shrinkB=0))
            ax.plot([x], [y], 'o', ms=2.0, color='0.15', zorder=5)
        else:
            ax.text(x, y, txt, ha='center', va='center', fontsize=7.2, color='white' if dark else 'black', linespacing=1.1,
                    path_effects=None if dark else W.halo(lw=1.6))
    geo_axes(ax, xlim=(-54.0, -46.9), small=True)
    fig.canvas.draw(); pa = ax.get_position(); pc = cb.ax.get_position()       # colour bar as tall as the map
    cb.ax.set_position([pc.x0, pa.y0, pc.width, pa.height])
    save(fig, 8)


def fig10():
    """BVOC by municipality, 2025: total and intensity, Jenks natural breaks (maps/12_maps_wf_style.py, fig7)."""
    B = bounds(); mun, state, MU = B['mun'], B['state'], B['MU']
    fig, axs = plt.subplots(2, 1, figsize=(W2, 8.4))
    for ax, col, ttl, sb, unit, cm in ((axs[0], 'BVOC', '(a) Total BVOC emissions, 2025', 'Municipal totals (t yr⁻¹); Jenks natural breaks', 't yr⁻¹', 'Greens'),
                                        (axs[1], 'intensity_kg_ha', '(b) BVOC intensity, 2025', 'Emission per unit municipal area (kg ha⁻¹ yr⁻¹); Jenks natural breaks',
                                         'kg ha⁻¹ yr⁻¹', 'YlGn')):
        v = mun.NM_MUN.map(MU[col]).values; q = W.fisher_jenks(v, 5)
        cols = W.ramp(cm, 5); k = np.clip(np.searchsorted(q, v, side='right') - 1, 0, 4); gg = mun.assign(_c=k)
        for i in range(5):
            sub = gg[gg._c == i]
            if len(sub): sub.plot(ax=ax, color=cols[i], edgecolor=cols[i], lw=0.15)
        hs = [Patch(facecolor=cols[i], edgecolor='0.5', lw=0.3, label=lb) for i, lb in enumerate(labels(q))]
        ax.legend(handles=hs, title=unit, loc='lower left', fontsize=6.8, title_fontsize=7.2, frameon=False, bbox_to_anchor=(0.0, 0.02),
                  handlelength=1.4, labelspacing=0.35, alignment='left')
        state.boundary.plot(ax=ax, color='0.2', lw=0.6); geo_axes(ax); ptitle(ax, ttl, sb)
    fig.tight_layout(h_pad=1.2); save(fig, 10)


def fig11():
    """Isoprene, CO and soil NO by municipality, 2019 and 2025, Jenks breaks of both years (maps/12_maps_wf_style.py, fig8)."""
    B = bounds(); mun, state, MU, MU19 = B['mun'], B['state'], B['MU'], B['MU19']
    fig = plt.figure(figsize=(W2, 7.4)); gs = fig.add_gridspec(3, 2, hspace=0.04, wspace=0.03, left=0.0, right=1.0, top=0.96, bottom=0.0)
    for j, (k, cm, nm) in enumerate((('ISOP', 'Greens', 'Isoprene'), ('CO', 'Blues', 'CO'), ('NO', 'Reds', 'Soil NO'))):
        q = W.fisher_jenks(np.concatenate([mun.NM_MUN.map(MU19[k]).values, mun.NM_MUN.map(MU[k]).values]), 5)
        for i, (yy, T_) in enumerate(((2019, MU19), (2025, MU))):
            ax = fig.add_subplot(gs[j, i]); W.extent(ax, axis_off=True)
            v = mun.NM_MUN.map(T_[k]).values; cols = W.ramp(cm, 5); kk = np.clip(np.searchsorted(q, v, side='right') - 1, 0, 4)
            for ci in range(5):
                sub = mun[kk == ci]
                if len(sub): sub.plot(ax=ax, color=cols[ci], edgecolor=cols[ci], lw=0.15)
            state.boundary.plot(ax=ax, color='0.1', lw=0.6)
            if j == 0: ax.set_title(str(yy), fontsize=8.5, fontweight='bold')
            if i == 0:
                hs = [Patch(facecolor=cols[ci], edgecolor='0.45', lw=0.4, label=lb) for ci, lb in enumerate(labels(q))]
                ax.legend(handles=hs, title=f'{nm} (t yr⁻¹)', loc='lower left', bbox_to_anchor=(-0.02, -0.03), frameon=False,
                          fontsize=6.6, title_fontsize=7.2, handlelength=1.5, handleheight=0.9, labelspacing=0.25, alignment='left')
            if j == 2 and i == 1:
                north(ax, -48.35, -28.25); W.scale_bar(ax, right=-48.25, lat_sb=-29.4, small=True, fs=6.3)
    save(fig, 11)


def fig1():
    """Study area (maps/study_area_map.py) without the source note, at the journal width."""
    import study_area_map as SA
    mun = gpd.read_file(os.path.join(BND, 'SC_Municipios_2024.gpkg')).to_crs(4674)
    mesog = gpd.read_file(os.path.join(BND, 'sc_mesoregions_2024.gpkg')).to_crs(4674)
    state = W.state_from(mun).to_crs(4674)
    uf = gpd.read_file(os.path.join(SA.EXT, 'ibge_br_uf_2024_simplified.geojson')).to_crs(4674)
    sa = gpd.read_file(os.path.join(SA.EXT, 'ibge_america_do_sul_2021_simplified.geojson')).to_crs(4674)
    g, extent, pct = SA.lulc(); tree, gext, sh = SA.growth_forms()
    GROUPS, GF, OCEAN, LAND, EDGE, halo = SA.GROUPS, SA.GF, SA.OCEAN, SA.LAND, SA.EDGE, SA.halo
    style(); num = lambda v: f'{v:.1f}'; fs0 = 6.8
    fig = plt.figure(figsize=(W2, 8.6))
    ax = fig.add_axes([0.305, 0.535, 0.685, 0.425]); ax.set_facecolor(OCEAN)
    sa.plot(ax=ax, color=LAND, edgecolor=EDGE, lw=0.6); state.plot(ax=ax, color=OCEAN, lw=0, zorder=1.5)
    uf[uf.SIGLA_UF != 'SC'].plot(ax=ax, color=LAND, edgecolor='#bdbdbd', lw=0.5)
    ax.imshow(g, extent=extent, cmap=ListedColormap([c for _, _, c, *_r in GROUPS]), vmin=-0.5, vmax=len(GROUPS) - 0.5, interpolation='nearest', zorder=2)
    mesog.boundary.plot(ax=ax, color='0.15', lw=0.45, zorder=3); state.boundary.plot(ax=ax, color='0.05', lw=0.8, zorder=4)
    pts = mun.set_index('NM_MUN').loc[SA.CITIES].geometry.representative_point()
    for name, p in pts.items():
        cap = name == 'Florianópolis'
        ax.plot(p.x, p.y, marker='*' if cap else 'o', ms=7 if cap else 3.0, mfc='#ffd400' if cap else 'k', mec='k', mew=0.5, zorder=7)
        dx, dy, ha = (0.07, 0.02, 'left') if name not in ('Florianópolis', 'Joinville') else (-0.07, 0.02, 'right')
        if cap: dx, dy, ha = (-0.08, -0.13, 'right')
        ax.text(p.x + dx, p.y + dy, name, fontsize=fs0, ha=ha, va='bottom', zorder=7, path_effects=halo)
    for lab, (x, y) in (('PR', (-51.3, -25.95)), ('RS', (-52.3, -29.25)), ('ARGENTINA', (-53.93, -27.65))):
        ax.text(x, y, lab, fontsize=fs0, color='0.45', ha='left', va='center', zorder=5)
    ax.text(-48.62, -28.95, 'Atlantic\nOcean', fontsize=fs0, style='italic', color='#3a6a96', ha='center', va='center', zorder=5)
    SA.frame(ax, 'en'); ax.tick_params(labelsize=fs0); ax.set_ylabel('Latitude (°)', fontsize=fs0 + 0.3, labelpad=2)
    ptitle(ax, '(a) Land use and land cover, 2019', 'MapBiomas Collection 10 (30 m), mesoregion boundaries, and reference municipalities')
    north(ax, -48.31, -25.88); W.scale_bar(ax, right=-48.35, lat_sb=-29.36, fs=6.4)
    axi = fig.add_axes([0.012, 0.765, 0.255, 0.215]); axi.set_facecolor(OCEAN)
    sa.plot(ax=axi, color='#d9d9d9', edgecolor='white', lw=0.3); uf.plot(ax=axi, color='#bdbdbd', edgecolor='white', lw=0.25)
    uf[uf.SIGLA_UF == 'SC'].plot(ax=axi, color='#c0392b', edgecolor='#c0392b', lw=0.4)
    axi.add_patch(Rectangle((-54.0, -29.45), 5.8, 3.65, fill=False, edgecolor='k', lw=0.6))
    axi.set_xlim(-82, -34); axi.set_ylim(-56, 13.5); axi.set_aspect(1 / math.cos(math.radians(-20))); axi.set_xticks([]); axi.set_yticks([])
    axi.text(0.04, 0.04, 'South America', transform=axi.transAxes, fontsize=fs0, va='bottom', ha='left'); axi.text(-52.5, -9.0, 'BRAZIL', fontsize=6.5, color='0.3', ha='center')
    axl = fig.add_axes([0.012, 0.525, 0.28, 0.225]); axl.axis('off')
    order = sorted(range(len(GROUPS)), key=lambda k: -pct[GROUPS[k][0]])
    hs = [Patch(facecolor=GROUPS[k][2], edgecolor='0.4', lw=0.3, label=f'{GROUPS[k][4]} ({num(pct[GROUPS[k][0]])}%)') for k in order]
    hs += [Line2D([], [], color='0.15', lw=0.6, label='Mesoregion boundary'),
           Line2D([], [], marker='*', ls='', ms=7, mfc='#ffd400', mec='k', mew=0.5, label='State capital'),
           Line2D([], [], marker='o', ls='', ms=3.0, mfc='k', mec='k', label='Municipality (approximate seat)')]
    axl.legend(handles=hs, loc='upper left', frameon=False, fontsize=6.5, title='Classes (% of the classified\nstate area)', title_fontsize=6.9,
               alignment='left', handlelength=1.5, handleheight=1.0, labelspacing=0.34, borderaxespad=0)
    ax = fig.add_axes([0.305, 0.05, 0.685, 0.425]); ax.set_facecolor(OCEAN)
    sa.plot(ax=ax, color=LAND, edgecolor=EDGE, lw=0.6); uf[uf.SIGLA_UF != 'SC'].plot(ax=ax, color=LAND, edgecolor='#bdbdbd', lw=0.5)
    state.plot(ax=ax, color='white', lw=0, zorder=1.5)
    bnd = np.linspace(0, 1, 11); cmap = ListedColormap(W.ramp('Greens', 10, 0.08, 0.95)); norm = BoundaryNorm(bnd, cmap.N)
    ax.imshow(np.ma.masked_invalid(tree), origin='lower', extent=gext, cmap=cmap, norm=norm, interpolation='nearest', zorder=2)
    mesog.boundary.plot(ax=ax, color='0.15', lw=0.45, zorder=3); state.boundary.plot(ax=ax, color='0.05', lw=0.8, zorder=4)
    SA.frame(ax, 'en'); ax.tick_params(labelsize=fs0)
    ax.set_xlabel('Longitude (°)', fontsize=fs0 + 0.3, labelpad=2); ax.set_ylabel('Latitude (°)', fontsize=fs0 + 0.3, labelpad=2)
    ptitle(ax, '(b) Tree fraction on the model grid, 2019', 'Fraction of each 0.025° cell covered by the MEGAN tree growth form')
    north(ax, -48.31, -25.88); W.scale_bar(ax, right=-48.35, lat_sb=-29.36, fs=6.4)
    cax = fig.add_axes([0.03, 0.385, 0.2, 0.013])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation='horizontal', ticks=[0, 0.25, 0.5, 0.75, 1])
    cb.ax.tick_params(labelsize=6.5); cb.set_label('Fraction of the cell covered\nby the tree growth form', fontsize=6.8)
    axs = fig.add_axes([0.115, 0.08, 0.12, 0.18]); y = np.arange(len(GF))[::-1]
    for yi, (k, pt_, en_, col) in zip(y, GF):
        axs.barh(yi, sh[k], color=col, edgecolor='0.4', lw=0.3, height=0.7); axs.text(sh[k] + 1.2, yi, f'{num(sh[k])}%', va='center', fontsize=6.5)
    axs.set_yticks(y); axs.set_yticklabels([en_ for _, pt_, en_, _ in GF], fontsize=6.6)
    axs.set_xlim(0, 80); axs.set_xticks([]); [axs.spines[s].set_visible(False) for s in ('top', 'right', 'bottom')]
    fig.text(0.012, 0.275, 'MEGAN growth forms in the\nstate, 2019 (% of the area)', fontsize=6.8, ha='left', va='bottom')
    save(fig, 1, tight=False)


def fig3():
    """Verification and observational evaluation (schematics/verification_flowchart.py), American English, journal width."""
    import verification_flowchart as VF
    N = VF.numbers(); H, S = VF.texts(N, 'en')
    rep = [('south-eastern USA', 'southeastern United States'), ('modelled', 'modeled'), ('Automatic checks', 'Automated tests'),
           (' checks', ' tests'), ('3-hourly against hourly steps', '3-h against 1-h steps'), ('soil-moisture', 'soil moisture'),
           ('column-emission', 'column–emission')]
    def am(s):
        for a, b in rep: s = s.replace(a, b)
        return s
    S = [{k: (am(v) if isinstance(v, str) else v) for k, v in s.items()} for s in S]
    for s in S:
        if 'tests' in s: s['tests'] = [(am(a), b, c) for a, b, c in s['tests']]
    H = {k: (am(v) if isinstance(v, str) else tuple(am(x) for x in v)) for k, v in H.items()}
    style(); fs = 6.9; fsq = 7.1
    W_IN, U = W2, 10.0
    xs = dict(card=6.4, how=7.8, res=28.2, mean=48.6, end=W_IN * U - 0.4, badge=3.2)
    cw = dict(how=38, res=35, mean=38)
    lh = fs * 1.28 / 72 * U
    fig = plt.figure(figsize=(W_IN, 12.0)); ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_IN * U); ax.set_ylim(120.0, 0); ax.axis('off')
    R = fig.canvas.get_renderer(); inv = ax.transData.inverted(); ext = lambda t: t.get_window_extent(renderer=R).transformed(inv)
    for k in ('how', 'res', 'mean'):
        ax.text(xs[k], 2.4, H[k], fontsize=6.8, color=VF.GREY, fontweight='bold', va='center', ha='left')
    y = 4.2; centres = []
    for i, s in enumerate(S):
        ph = 0 if i < len(S) - 1 else 1; col = (VF.C_VER, VF.C_EVA)[ph]; fill = (VF.F_VER, VF.F_EVA)[ph]
        if i in (0, len(S) - 1):
            if ph: y += 0.6
            t0 = ax.text(xs['card'] + 0.4, y + 2.0, H['ph'][ph], fontsize=7.6, fontweight='bold', color=col, va='center')
            ax.text(ext(t0).x1 + 0.8, y + 2.0, '· ' + H['phs'][ph], fontsize=7.2, color=VF.GREY, va='center', style='italic')
            y += 4.0
        arts = [ax.text(xs['how'], y + 2.1, s['t'], fontsize=8.4, fontweight='bold', color=VF.INK, va='center')]
        ax.text(xs['end'] - 1.0, y + 2.1, ('✓ ' if s['st'] == 'ok' else '◐ ') + s['chip'], fontsize=7.0,
                color=VF.OKC if s['st'] == 'ok' else col, ha='right', va='center', fontweight='bold')
        q = ax.text(xs['how'], y + 3.6, s['q'], fontsize=fsq, color=VF.GREY, va='top', style='italic'); arts.append(q)
        yb = ext(q).y0 + 1.2 if ext(q).y0 > ext(q).y1 else ext(q).y1 + 1.2
        ax.plot([xs['how'], xs['end'] - 1.0], [yb - 0.6, yb - 0.6], color=VF.LINE, lw=0.5, zorder=2)
        arts.append(ax.text(xs['how'], yb, textwrap.fill(s['how'], cw['how']), fontsize=fs, color=VF.INK, va='top', linespacing=1.28))
        arts.append(ax.text(xs['mean'], yb, textwrap.fill(s['mean'], cw['mean']), fontsize=fs, color=VF.INK, va='top', linespacing=1.28))
        if 'tests' in s:
            yy = yb
            for nm, val, ok in s['tests']:
                t0 = ax.text(xs['res'], yy, ('● ' if ok else '○ ') + nm + ':', fontsize=fs, color=VF.OKC if ok else VF.NCC, va='top', fontweight='bold')
                x1 = ext(t0).x1 + 0.5; t1 = ax.text(x1, yy, val, fontsize=fs, color=VF.INK if ok else VF.GREY, va='top'); arts += [t0, t1]
                if ext(t1).x1 > xs['mean'] - 1.0 and ' (' in val:
                    a_, b_ = val.split(' (', 1); t1.set_text(a_); yy += 1.15 * lh
                    arts.append(ax.text(x1, yy, '(' + b_, fontsize=fs, color=VF.INK if ok else VF.GREY, va='top'))
                yy += 1.15 * lh
        else:
            k_ = ax.text(xs['res'], yb - 0.1, s['key'], fontsize=9.4, fontweight='bold', color=col, va='top'); arts.append(k_)
            arts.append(ax.text(xs['res'], max(ext(k_).y0, ext(k_).y1) + 0.6, textwrap.fill(s['det'], cw['res']), fontsize=fs,
                                color=VF.GREY, va='top', linespacing=1.28))
        bottom = max(max(ext(a).y0, ext(a).y1) for a in arts) + 1.1
        ax.add_patch(FancyBboxPatch((xs['card'], y), xs['end'] - xs['card'], bottom - y, boxstyle='round,pad=0,rounding_size=1.0', fc=fill, ec='none', zorder=0))
        cy = y + 2.1; centres.append((cy, col))
        ax.add_patch(Circle((xs['badge'], cy), 1.85, fc=col, ec='white', lw=1.2, zorder=4))
        ax.text(xs['badge'], cy + 0.05, str(i + 1), color='white', fontsize=8, fontweight='bold', ha='center', va='center', zorder=5)
        y = bottom + 1.2
    for (y0, c0), (y1, c1) in zip(centres[:-1], centres[1:]):
        ax.plot([xs['badge']] * 2, [y0 + 1.9, y1 - 1.9], color=c0 if c0 == c1 else '#b9b4d9', lw=1.4, zorder=3, solid_capstyle='butt')
    ax.set_ylim(y + 0.2, 0); fig.set_size_inches(W_IN, (y + 0.2) / U)
    save(fig, 3, tight=False)


def tropomi_data():
    """Data of the TROPOMI figure, prepared exactly as in fig9_tropomi.py (cells 1-5 of notebook 14)."""
    from unittest import mock
    TAG = P.TAG; src = open(os.path.join(P.NB_SRC, '14_tropomi_hcho_evaluation.py'), encoding='utf-8').read()
    g = types.ModuleType('google'); col = types.ModuleType('google.colab'); col.drive = mock.MagicMock(); g.colab = col
    sys.modules['google'] = g; sys.modules['google.colab'] = col; sys.modules['ee'] = mock.MagicMock()
    cells = [c for c in src.split('# %%') if c.strip()]
    code = lambda c: '\n'.join(l for l in c.split('\n') if not l.lstrip().startswith(('%', '!')))
    Gd = {'__name__': '__main__'}; exec(code(cells[1]), Gd)
    work = P.out('_work', 'fig9'); os.makedirs(os.path.join(work, 'cache_' + TAG), exist_ok=True)
    Gd.update(MODEL_ROOT=RUN, BND=BND, OUT_DIR=work); exec(code(cells[2]), Gd)
    for f in glob.glob(os.path.join(P.TROPOMI, '*.npz')): shutil.copyfile(f, os.path.join(work, 'cache_' + TAG, os.path.basename(f)))
    exec(code(cells[3]), Gd); exec(code(cells[4]).split('T1 = []')[0], Gd); exec(code(cells[5]).split('T3 = ')[0], Gd)
    SAT, MOD, SC, MONTHS, SER, lat, lon, RES, block = (Gd[k] for k in ('SAT', 'MOD', 'SC', 'MONTHS', 'SER', 'lat', 'lon', 'RES', 'block'))
    DJF = [(y, m) for (y, m) in MONTHS if m in (12, 1, 2)]
    E = np.nanmean(np.stack([MOD[k] for k in DJF]), axis=0)
    Hh = np.nanmean(np.stack([SAT[k]['H'] - float(SAT[k]['bg']) for k in DJF]), axis=0)
    def band(path, name):
        with rasterio.open(path) as s_: return np.flipud(s_.read(list(s_.descriptions).index(name) + 1).astype('float64'))
    Pp = np.nanmean(np.stack([band(f'{RUN}/{y}/inputs/met_{y}_{m:02d}.tif', 'PRES') for (y, m) in DJF]), axis=0) / 100.0
    S = SER.copy(); clim = S.groupby('month')[['ISOP_t', 'dHCHO_SC', 'BG']].mean(); an = S.copy()
    for c in ('ISOP_t', 'dHCHO_SC'): an[c] = S[c] - S.groupby('month')[c].transform('mean')
    rr = Gd['r_neff'](an.ISOP_t, an.dHCHO_SC); r2 = np.corrcoef(clim.ISOP_t, clim.dHCHO_SC)[0, 1]
    frac = block(SC.astype(float), 8)
    hb = np.where(frac >= 0.5, block(np.where(SC, Hh, np.nan), 8), np.nan).ravel(); pb = np.where(frac >= 0.5, block(np.where(SC, Pp, np.nan), 8), np.nan).ravel()
    eb = np.where(frac >= 0.5, block(np.where(SC, E, np.nan), 8), np.nan).ravel(); ok = np.isfinite(hb) & np.isfinite(pb) & np.isfinite(eb)
    r_hp = np.corrcoef(hb[ok], pb[ok])[0, 1]
    print(f'check: r_clim {r2:.3f} | anomalies r {rr["r"]:.3f} N_eff {rr["n_eff"]:.1f} p {rr["p_neff"]:.5f} | r(dHCHO, p_s) 0.2 deg {r_hp:.3f} n {int(ok.sum())}')
    return dict(E=E, H=Hh, P=Pp, SC=SC, S=S, clim=clim, an=an, rr=rr, r2=r2, hb=hb, pb=pb, eb=eb, ok=ok, r_hp=r_hp, lat=lat, lon=lon, RES=RES)


def fig12_13():
    D = tropomi_data(); lat, lon, RES, SC = D['lat'], D['lon'], D['RES'], D['SC']; style()
    ext = [lon[0] - RES / 2, lon[-1] + RES / 2, lat[0] - RES / 2, lat[-1] + RES / 2]
    state = W.state_from(gpd.read_file(os.path.join(BND, 'SC_Municipios_2024.gpkg')))
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.75)); fig.subplots_adjust(left=0.075, right=0.95, top=0.935, bottom=0.075, hspace=0.30, wspace=0.30)
    for ax, fld, ttl, sub, cm, cbl in ((axs[0, 0], np.where(SC, D['E'], np.nan) / 1e12, '(a) Modeled isoprene at satellite overpass',
                                        'December–February mean, 2019–2025', 'Greens', '10¹² C atoms cm⁻² s⁻¹'),
                                       (axs[0, 1], np.where(SC, D['H'], np.nan) / 1e15, '(b) TROPOMI ΔHCHO', 'December–February mean, 2019–2025',
                                        'YlOrBr', '10¹⁵ molecules cm⁻²'),
                                       (axs[1, 0], np.where(SC, D['P'], np.nan), '(c) ERA5-Land surface pressure',
                                        'December–February mean, 2019–2025; post hoc diagnostic', 'cividis', 'hPa')):
        im = ax.imshow(fld, origin='lower', extent=ext, cmap=cm, interpolation='nearest'); cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(cbl, fontsize=6.8); cb.ax.tick_params(labelsize=6.5)
        state.boundary.plot(ax=ax, color='0.2', lw=0.5); W.extent(ax); ax.tick_params(labelsize=6.5)
        ax.xaxis.set_major_formatter(W.fmt_deg('en')); ax.yaxis.set_major_formatter(W.fmt_deg('en'))
        north(ax, W.NORTH_X, -25.86, w=6.5, h=9.5, fs=7.8); W.scale_bar(ax, small=True, fs=6.0); ptitle(ax, ttl, sub, fs=8.2, fs2=6.8)
    ax = axs[1, 1]; hb, pb, eb, ok = D['hb'], D['pb'], D['eb'], D['ok']; ax.set_box_aspect(axs[1, 0].get_box_aspect() or 0.66)
    sc_ = ax.scatter(pb[ok], hb[ok] / 1e15, c=eb[ok] / 1e12, cmap='Greens', s=8, edgecolor='0.3', linewidth=0.25)
    cb = fig.colorbar(sc_, ax=ax, shrink=0.85, pad=0.02); cb.set_label('Modeled isoprene\n(10¹² C atoms cm⁻² s⁻¹)', fontsize=6.8); cb.ax.tick_params(labelsize=6.5)
    ax.set_xlabel('Surface pressure (hPa)'); ax.set_ylabel('ΔHCHO (10¹⁵ molecules cm⁻²)'); clean(ax)
    ptitle(ax, '(d) ΔHCHO versus surface pressure', f'0.2° blocks; r = {D["r_hp"]:.2f}, n = {int(ok.sum())}; post hoc diagnostic', fs=8.2, fs2=6.8)
    save(fig, 12, tight=False)
    clim, S, an, rr, r2 = D['clim'], D['S'], D['an'], D['rr'], D['r2']
    fig, axs = plt.subplots(1, 2, figsize=(W2, 3.05)); fig.subplots_adjust(left=0.065, right=0.985, top=0.89, bottom=0.25, wspace=0.22)
    ax = axs[0]; mm = np.arange(1, 13)
    ax.plot(mm, clim.ISOP_t / clim.ISOP_t.max(), 'o-', ms=3, color='#2e8b57', label='Modeled isoprene')
    ax.plot(mm, clim.dHCHO_SC / clim.dHCHO_SC.max(), 's-', ms=3, color='#b8860b', label='TROPOMI ΔHCHO')
    ax.plot(mm, clim.BG / clim.BG.max(), '^--', ms=3, color='0.45', label='Atlantic HCHO background')
    ax.set_xticks(mm); ax.set_xticklabels(MON, fontsize=6.3); ax.set_ylabel('Normalized by the maximum'); clean(ax); grid_y(ax)
    ptitle(ax, '(a) Seasonal climatology, 2019–2025', MINUS(f'Normalized by each series maximum; r (isoprene, ΔHCHO) = {r2:.2f}'))
    ax.legend(frameon=False, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=3, columnspacing=0.9)
    ax = axs[1]; t = pd.to_datetime(dict(year=S.year, month=S.month, day=15)); z = lambda s: (s - s.mean()) / s.std()
    ax.plot(t, z(an.ISOP_t), color='#2e8b57', lw=0.9, label='Modeled isoprene'); ax.plot(t, z(an.dHCHO_SC), color='#b8860b', lw=0.9, label='TROPOMI ΔHCHO')
    ax.axhline(0, color='0.6', lw=0.6); ax.set_ylabel('z-score'); clean(ax)
    ptxt = 'p < 0.001' if rr['p_neff'] < 0.001 else f"p = {rr['p_neff']:.3f}"
    ptitle(ax, '(b) Deseasonalized monthly anomalies', f"z-scores; r = {rr['r']:.2f}, effective N = {rr['n_eff']:.0f}, {ptxt}")
    ax.legend(frameon=False, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    save(fig, 13, tight=False)


FIGS = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6, 7: fig7, 8: fig8, 9: fig9, 10: fig10, 11: fig11, 12: fig12_13}
if __name__ == '__main__':
    want = [int(a) for a in sys.argv[1:]] or sorted(FIGS)
    for n in want:
        style(); FIGS[12 if n == 13 else n]()
