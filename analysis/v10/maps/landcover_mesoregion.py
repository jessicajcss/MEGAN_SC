#!/usr/bin/env python3
"""landcover_mesoregion.py - share of each MapBiomas land-cover group by mesoregion of Santa Catarina, year YEAR.

The 2019 table used until 16 November 2026 (data/boundaries/landcover2019_mesoregion.csv) was made with the same method
before this script was added: the ~120 m MapBiomas window, pixel areas by
latitude on the authalic sphere, pixels assigned to the mesoregion polygon (union of its 2024 municipalities) that holds
the pixel centre, and
    pct_of_mesoregion = area of the group in the mesoregion / classified area (class > 0) of the mesoregion
    pct_of_group      = area of the group in the mesoregion / area of the group in the six mesoregions.
Groups and year are those of study_area_map.py (FIG1_YEAR; default 2024, the last year of MapBiomas Collection 10).

Inputs : data/external/mapbiomas_c10_lulc<YEAR>_SC_120m.npz (extract_mapbiomas_sc.py), data/boundaries/SC_Municipios_2024.gpkg
         and data/boundaries/sc_municipios_2024_mesorregiao.csv (make_sc_mesoregions_2024.py). MUN_SHP may point to another
         copy of the IBGE 2024 municipal mesh (the 2019 table was made with the full-resolution shapefile; the gpkg
         of this repository is a simplified copy); the municipalities are matched by CD_MUN.
Output : <MEGAN_SC_BND_OUT or MEGAN_SC_OUT/boundaries>/landcover<YEAR>_mesoregion.csv; the copy shipped in data/boundaries
         was made by this script (MEGAN_SC_BND_OUT=data/boundaries).
"""
import os, sys
import numpy as np, pandas as pd, geopandas as gpd
from affine import Affine
from rasterio.features import rasterize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import study_area_map as S                       # GROUPS, YEAR, EXT, BND (no side effects on import)
sys.path.insert(0, os.path.join(HERE, '..')); import _paths as P

BOUT = os.environ.get('MEGAN_SC_BND_OUT') or P.out('boundaries')
os.makedirs(BOUT, exist_ok=True)


def table(year=S.YEAR, mun_file=None):
    z = np.load(os.path.join(S.EXT, f'mapbiomas_c10_lulc{year}_SC_120m.npz'))
    a, inside, t = z['a'], z['inside'], Affine(*z['transform'][:6])
    mun = gpd.read_file(mun_file or os.path.join(S.BND, 'SC_Municipios_2024.gpkg')).to_crs(4326)
    attr = pd.read_csv(os.path.join(S.BND, 'sc_municipios_2024_mesorregiao.csv'), dtype={'CD_MUN': str})
    mun['meso'] = mun.CD_MUN.astype(str).map(attr.set_index('CD_MUN').mesorregiao)
    assert len(mun) == 295 and mun.meso.notna().all()
    mg = mun[['meso', 'geometry']].dissolve(by='meso').reset_index()
    lab = rasterize([(g, k + 1) for k, g in enumerate(mg.geometry)], out_shape=a.shape, transform=t, fill=0)
    R = 6371007.2
    lat = t.f + (np.arange(a.shape[0]) + 0.5) * t.e
    dl, dp = np.radians(t.a), np.radians(-t.e)
    pa = np.broadcast_to(((R ** 2) * dl * (np.sin(np.radians(lat) + dp / 2) - np.sin(np.radians(lat) - dp / 2)))[:, None], a.shape)
    rows = []
    for key, ids, *_r in S.GROUPS:
        sel = np.isin(a, ids) & inside & (lab > 0)
        tot = pa[sel].sum()
        for k, name in enumerate(mg.meso):
            cl = inside & (lab == k + 1) & (a > 0)
            g = pa[sel & (lab == k + 1)].sum()
            rows.append(dict(group=key, mesorregiao=name, km2=round(g / 1e6, 1), pct_of_group=round(100 * g / tot, 2),
                             pct_of_mesoregion=round(100 * g / pa[cl].sum(), 2)))
    return pd.DataFrame(rows)


if __name__ == '__main__':
    d = table(mun_file=os.environ.get('MUN_SHP'))
    out = os.path.join(BOUT, f'landcover{S.YEAR}_mesoregion.csv')
    d.to_csv(out, index=False)
    print('written', out)
    print(d.pivot(index='mesorregiao', columns='group', values='pct_of_mesoregion').round(2).to_string())
