"""wf_style.py - map style of the biomass-burning article (wildfire/FINN_SC_pipeline_v2/python/analysis/a10_figures.py and
tools/study_area_map.py), shared by the map figures of the biogenic manuscript and report.

Frame: geographic axes (−54 to −48 °E, −29,5 to −25,8 °N; aspect 1/cos 27,6°), degree labels with decimal comma and true
minus in Portuguese; north arrow over the Atlantic; 0-50-100 km scale bar inside the frame at the bottom right of every
map; one source note per figure, below the frames; state outline (0.2 grey, 0.6 pt) and, on request, mesoregion lines (0.35 grey, 0.45 pt); municipalities
filled without outlines; class legends at lower left, without frame; no figure title (the caption is in the document);
DejaVu Sans 8 pt. Multi-panel grids: no axes, one north arrow and one scale bar in the free panel at the bottom right.
Municipal classes: Jenks natural breaks computed with Fisher's exact algorithm (fisher_jenks); labels of polygons at
the interior point farthest from the edges (label_point)."""
import math
import numpy as np, matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter

LAT0 = -27.6; NORTH_X = -48.28
XLIM, YLIM = (-54, -48), (-29.5, -25.8)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8})

def num(v, d=0, lang='pt'):
    s = f'{v:,.{d}f}'
    s = s.replace(',', 'X').replace('.', ',').replace('X', '.') if lang == 'pt' else s
    return s.replace('-', '−')

def fmt_deg(lang):
    return FuncFormatter(lambda v, _: (f'{v:g}').replace('.', ',' if lang == 'pt' else '.').replace('-', '−'))

def north_arrow(ax, x=NORTH_X, y0=-26.50, y1=-26.08, fs=8):
    ax.annotate('', xy=(x, y1), xytext=(x, y0), arrowprops=dict(arrowstyle='-|>', color='k', lw=1.0))
    ax.text(x, y1 + 0.04, 'N', ha='center', va='bottom', fontsize=fs, fontweight='bold')

def scale_bar(ax, x0=None, lat_sb=-29.33, small=False, right=None, fs=None):
    """0-50-100 km bar; by default at the bottom right of the frame (right edge at XLIM[1] - 0.12, or at `right`). The first
    label is left-aligned and the last right-aligned on the ends of the bar, so that bar and labels stay inside the frame."""
    dlon = 100 / (111.32 * math.cos(math.radians(lat_sb)))
    if x0 is None: x0 = (XLIM[1] - 0.12 if right is None else right) - dlon
    ax.add_patch(Rectangle((x0, lat_sb), dlon / 2, 0.035, facecolor='k', edgecolor='k', lw=0.4))
    ax.add_patch(Rectangle((x0 + dlon / 2, lat_sb), dlon / 2, 0.035, facecolor='w', edgecolor='k', lw=0.4))
    fs = fs or (5.5 if small else 6.5)
    for k, (lab, ha) in enumerate((('0', 'left'), ('50', 'center'), ('100 km', 'right'))):
        if small and k == 1: continue
        ax.text(x0 + k * dlon / 2, lat_sb + 0.07, lab, ha=ha, va='bottom', fontsize=fs)

def extent(ax, axis_off=False):
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); ax.set_aspect(1 / math.cos(math.radians(LAT0)))
    if axis_off: ax.axis('off')

def base(ax, note, lang='pt', small=False, xlim=None, fs=None, note_fs=6.3, note_y=None):
    """Geographic frame, north arrow, scale bar (bottom right, inside the frame) and the source note below the frame
    (`note` = None for panels whose note is given once for the whole figure); `xlim` widens the frame."""
    extent(ax)
    if xlim is not None: ax.set_xlim(*xlim)
    right = (xlim or XLIM)[1]
    ax.xaxis.set_major_formatter(fmt_deg(lang)); ax.yaxis.set_major_formatter(fmt_deg(lang))
    ax.set_xlabel('Longitude (°)'); ax.set_ylabel('Latitude (°)')
    north_arrow(ax, x=right - 0.35 if xlim is not None else NORTH_X, fs=(fs + 1.5) if fs else 8); scale_bar(ax, small=small, right=right - 0.12, fs=fs)
    if note:
        ax.text(0.0, note_y if note_y is not None else (-0.105 if not small else -0.09), note, transform=ax.transAxes, fontsize=note_fs, va='top', ha='left')

def outlines(ax, state, mesog=None, lw_state=0.6, color_state='0.2'):
    if mesog is not None: mesog.boundary.plot(ax=ax, color='0.35', lw=0.45)
    state.boundary.plot(ax=ax, color=color_state, lw=lw_state)

def magma_palette(n, zero=True):
    cm = plt.get_cmap('magma_r'); cols = [cm(x) for x in np.linspace(0.10, 0.92, n)]
    return (['#eeeeee'] + cols) if zero else cols

def ramp(cmap, n, lo=0.12, hi=0.95):
    cm = plt.get_cmap(cmap); return [cm(x) for x in np.linspace(lo, hi, n)]

def class_choropleth(ax, gdf, values, bins, colors, labels, legend_title, small=False, legend_loc=(0.0, 0.02)):
    """bins: class limits b0 < b1 < ... < bn; class i = [b_i, b_(i+1)) and the last class includes bn.
    Municipalities filled without outlines; legend with the class limits only (every map shows all 295 municipalities)."""
    v = np.asarray(values, float); k = np.clip(np.searchsorted(bins, v, side='right') - 1, 0, len(bins) - 2)
    g = gdf.assign(_c=k)
    for i in range(len(bins) - 1):
        sub = g[g._c == i]
        if len(sub): sub.plot(ax=ax, color=colors[i], edgecolor=colors[i], lw=0.15)
    counts = [int((k == i).sum()) for i in range(len(bins) - 1)]
    hs = [Patch(facecolor=colors[i], edgecolor='0.5', lw=0.3, label=labels[i]) for i in range(len(bins) - 1)]
    ax.legend(handles=hs, title=legend_title, loc='lower left', fontsize=5.0 if small else 6.3, title_fontsize=5.4 if small else 6.8,
              frameon=False, bbox_to_anchor=legend_loc, handlelength=1.4, labelspacing=0.25 if small else 0.4, alignment='left')
    return counts

def save(fig, path, dpi=300):
    fig.savefig(path, dpi=dpi, bbox_inches='tight', pad_inches=0.08, facecolor='white'); plt.close(fig)

def fill_slivers(geom, max_km2=1.0):
    """Remove interior rings smaller than max_km2: gaps left by the union of the simplified municipal polygons of
    SC_Municipios_2024.gpkg (each municipality was simplified on its own)."""
    import geopandas as gpd
    from shapely.geometry import Polygon, MultiPolygon
    eq = gpd.GeoSeries([geom], crs='EPSG:4326').to_crs('ESRI:102033').iloc[0]
    parts = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
    parts_eq = list(eq.geoms) if isinstance(eq, MultiPolygon) else [eq]
    out = [Polygon(p.exterior, [i for i, ie in zip(p.interiors, pe.interiors) if Polygon(ie).area / 1e6 >= max_km2])
           for p, pe in zip(parts, parts_eq)]
    return MultiPolygon(out) if len(out) > 1 else out[0]

def state_from(mun):
    import geopandas as gpd
    u = mun.to_crs(4326).geometry.union_all() if hasattr(mun.geometry, 'union_all') else mun.to_crs(4326).geometry.unary_union
    return gpd.GeoDataFrame(geometry=[fill_slivers(u)], crs='EPSG:4326')

def fisher_jenks(values, k=5):
    """Jenks natural breaks by Fisher's exact dynamic programme: the partition of the sorted values into k contiguous classes
    with the smallest total within-class sum of squared deviations. Returns the k + 1 limits (minimum, k - 1 upper limits of
    the first classes, maximum); with class_choropleth, class i holds b_i <= v < b_(i+1) and the last class includes the
    maximum. Same limits as mapclassify.FisherJenks (checked on the values of Figures 7 and 8)."""
    x = np.sort(np.asarray(values, float)); x = x[np.isfinite(x)]; n = x.size
    if n <= k: return np.unique(x)
    c1 = np.concatenate([[0.0], np.cumsum(x)]); c2 = np.concatenate([[0.0], np.cumsum(x * x)])
    D = np.full((k + 1, n + 1), np.inf); B = np.zeros((k + 1, n + 1), int); D[0, 0] = 0.0
    for m in range(1, k + 1):
        for j in range(m, n + 1):
            i = np.arange(m - 1, j); s1 = c1[j] - c1[i]; cost = D[m - 1, i] + (c2[j] - c2[i]) - s1 * s1 / (j - i)
            a = int(np.argmin(cost)); D[m, j], B[m, j] = cost[a], i[a]
    ends, j = [], n
    for m in range(k, 0, -1):
        ends.append(j); j = B[m, j]
    ends = sorted(ends)
    # upper limit of class i = first value of class i + 1, so that searchsorted(..., side='right') - 1 keeps each value in its class
    return np.array([x[0]] + [x[e] for e in ends[:-1]] + [x[-1]])

def label_point(geom, lat0=LAT0):
    """Interior point of the largest part of a polygon that is farthest from its edges (shapely polylabel), computed in
    coordinates scaled by cos(lat0) so that distances are isotropic on the map."""
    from shapely.geometry import MultiPolygon
    from shapely import affinity
    from shapely.ops import polylabel
    part = max(geom.geoms, key=lambda g: g.area) if isinstance(geom, MultiPolygon) else geom
    c = math.cos(math.radians(lat0)); p = polylabel(affinity.scale(part, xfact=c, yfact=1.0, origin=(0, 0)), tolerance=0.002)
    return p.x / c, p.y

def halo(color='white', lw=2.2):
    return [patheffects.withStroke(linewidth=lw, foreground=color)]
