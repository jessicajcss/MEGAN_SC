# -*- coding: utf-8 -*-
"""landuse_context.py - land use of the mesoregions and municipalities of Santa Catarina (POST HOC, 16 November 2026).

Describes the land use behind the municipal and mesoregion emissions; it does not change any emission.
1. MapBiomas Collection 10, 2024: area of each land-cover group of maps/study_area_map.py (native forest, grassland,
   wetland, plantation forestry, pasture, temporary crops, mixed use, urban, other nonvegetated, water) in each
   municipality (columns mb_*), from the 120-m map of data/external (pixels assigned to the municipality that holds the pixel centre,
   pixel areas on the authalic sphere, as in maps/landcover_mesoregion.py). Shares are given against the classified area
   (class > 0). Only 2024 is used: the 120-m map of 2019 in data/external is a different export and is not compared with
   it.
2. Growth forms of the production run (inputs/growthform_<year>.npz, built every year from the MapBiomas map of that
   year): area of tree, herbaceous and crop growth forms and of the vegetated area (their sum) in each municipality in
   2003 and 2024 (columns gf_*), with each cell allocated to the municipalities in proportion to the overlap area, as for the municipal
   emissions. The nonvegetated classes (urban areas, water, beaches and rocky outcrops) are the rest of the area.

Usage: python analysis/v10/review/landuse_context.py   (writes <MEGAN_SC_OUT>/posthoc_review/landuse_context_*.csv;
published copy in results/v10/posthoc_review; also run by `python reproduce.py posthoc`)
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
from affine import Affine
from rasterio.features import rasterize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(HERE, '..', 'maps'))
import _paths as P
import common as C
import study_area_map as S                       # GROUPS (no side effects on import)

OUT = os.path.join(P.OUT, 'posthoc_review'); os.makedirs(OUT, exist_ok=True)
YEAR_LC, Y0, Y1 = 2024, 2003, 2024
FORMS = ('tree', 'herb', 'crop', 'veg')

mun = gpd.read_file(os.path.join(P.BND, 'SC_Municipios_2024.gpkg')).to_crs(4326)
attr = pd.read_csv(os.path.join(P.BND, 'sc_municipios_2024_mesorregiao.csv'), dtype={'CD_MUN': str})
mun['CD_MUN'] = mun.CD_MUN.astype(str)
mun['meso'] = mun.CD_MUN.map(attr.set_index('CD_MUN').mesorregiao)
assert len(mun) == 295 and mun.meso.notna().all()
tab = mun[['CD_MUN', 'NM_MUN', 'meso']].copy()

# 1. MapBiomas 2024 groups by municipality
z = np.load(os.path.join(P.REPO, 'data', 'external', f'mapbiomas_c10_lulc{YEAR_LC}_SC_120m.npz'))
a, inside, t = z['a'], z['inside'], Affine(*z['transform'][:6])
lab = rasterize([(g, k + 1) for k, g in enumerate(mun.geometry)], out_shape=a.shape, transform=t, fill=0)
R = 6371007.2
lat_p = t.f + (np.arange(a.shape[0]) + 0.5) * t.e
dl, dp = np.radians(t.a), np.radians(-t.e)
pa = np.broadcast_to(((R ** 2) * dl * (np.sin(np.radians(lat_p) + dp / 2) - np.sin(np.radians(lat_p) - dp / 2)))[:, None], a.shape)
m = inside & (lab > 0) & (a > 0)
L, A, PA = lab[m], a[m], pa[m]
n = len(mun) + 1
tab[f'mb_classified_km2_{YEAR_LC}'] = np.bincount(L, weights=PA, minlength=n)[1:] / 1e6
for key, ids, *_ in S.GROUPS:
    sel = np.isin(A, ids)
    tab[f'mb_{key}_km2_{YEAR_LC}'] = np.bincount(L[sel], weights=PA[sel], minlength=n)[1:] / 1e6

# 2. growth forms of the run, 2003 and 2024
lat, lon = C.grid(); area = C.cell_area_m2(lat, lon).ravel() / 1e6
Wm = C.area_weights(list(mun.geometry), lat, lon)
tab['gf_grid_km2'] = Wm @ area
for year in (Y0, Y1):
    d = np.load(os.path.join(P.RUN, str(year), 'inputs', f'growthform_{year}.npz'))
    for f in FORMS:
        tab[f'gf_{f}_km2_{year}'] = Wm @ (np.nan_to_num(d[f]).ravel() * area)

tab = tab.round(3)
tab.to_csv(os.path.join(OUT, 'landuse_context_municipality.csv'), index=False)
num = [c for c in tab.columns if c not in ('CD_MUN', 'NM_MUN', 'meso')]
meso = tab.groupby('meso')[num].sum()
meso.loc['Santa Catarina (sum of municipalities)'] = tab[num].sum()
meso.round(3).to_csv(os.path.join(OUT, 'landuse_context_mesoregion.csv'))

# printed summary (shares in % of the classified area; changes of the growth forms in km2 and %)
pct = lambda d, k: 100 * d[f'mb_{k}_km2_{YEAR_LC}'] / d[f'mb_classified_km2_{YEAR_LC}']
show = pd.DataFrame({k + '%': pct(meso, k) for k in ('forest', 'silv', 'urban', 'crop', 'past', 'mosaic')})
for f in FORMS:
    show['d' + f + '_km2'] = meso[f'gf_{f}_km2_{Y1}'] - meso[f'gf_{f}_km2_{Y0}']
show['dveg_%'] = 100 * show['dveg_km2'] / meso[f'gf_veg_km2_{Y0}']
print(show.round(2).to_string())
mt = tab.set_index('NM_MUN')
dveg = 100 * (mt[f'gf_veg_km2_{Y1}'] - mt[f'gf_veg_km2_{Y0}']) / mt[f'gf_veg_km2_{Y0}']
print('largest relative losses of vegetated area 2003-2024 (%):', dveg.sort_values().head(8).round(2).to_dict())
print('largest urban areas 2024 (km2):', mt[f'mb_urban_km2_{YEAR_LC}'].sort_values(ascending=False).head(6).round(1).to_dict())
print('written', os.path.join(OUT, 'landuse_context_municipality.csv'), 'and landuse_context_mesoregion.csv')
