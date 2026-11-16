#!/usr/bin/env python3
"""extract_mapbiomas_sc.py - Santa Catarina window of a MapBiomas Collection 10 land use and land cover GeoTIFF, for
the study-area figure (study_area_map.py) and the land-cover shares by mesoregion (landcover_mesoregion.py).

Writes, for one year:
  data/external/mapbiomas_c10_lulc<year>_SC_120m.npz        arrays `a` (class id), `inside` (pixel centre inside the
      state) and `transform` (affine, EPSG:4326), read at 1/4 of the native 30 m resolution (about 120 m) with
      nearest-neighbour resampling;
  data/external/mapbiomas_c10_lulc<year>_SC_class_areas.csv class areas inside the state (km2) and their share of the
      classified area (class 0 excluded), read at 1/2 of the native resolution (about 60 m), pixel areas by latitude on
      the WGS 84 authalic sphere.
The state is the union of the IBGE 2024 municipalities (data/boundaries/SC_Municipios_2024.gpkg).
The 2019 files of data/external were made with the same method before this script was added. Input: the GeoTIFF exported by mapbiomas_lulc_SC_export.ipynb (this folder; Collection 10,
native grid) or the Brazil coverage file of the collection.
Usage: python analysis/v10/maps/extract_mapbiomas_sc.py <geotiff> <year>
"""
import os, sys
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.windows import from_bounds
from rasterio.features import geometry_mask
from rasterio.enums import Resampling

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); import _paths as P
EXT = os.path.join(P.REPO, 'data', 'external')
NAMES = {0: 'Sem classificação (águas costeiras dentro do limite do IBGE)', 3: 'Formação florestal', 5: 'Mangue',
         9: 'Silvicultura', 11: 'Campo alagado e área pantanosa', 12: 'Formação campestre', 15: 'Pastagem',
         21: 'Mosaico de usos', 23: 'Praia, duna e areal', 24: 'Área urbanizada', 25: 'Outras áreas não vegetadas',
         29: 'Afloramento rochoso', 30: 'Mineração', 31: 'Aquicultura', 32: 'Apicum', 33: 'Rio, lago e oceano',
         39: 'Soja', 40: 'Arroz', 41: 'Outras lavouras temporárias', 49: 'Restinga arbórea', 50: 'Restinga herbácea',
         75: 'Usina fotovoltaica (beta)'}


def authalic_row_areas(t, h):
    """Area (m2) of one pixel in each row of a north-up EPSG:4326 grid, on the WGS 84 authalic sphere."""
    a_, f_ = 6378137.0, 1 / 298.257223563
    e2 = f_ * (2 - f_); e_ = np.sqrt(e2)
    q = lambda phi: (1 - e2) * (np.sin(phi) / (1 - e2 * np.sin(phi) ** 2)
                                - np.log((1 - e_ * np.sin(phi)) / (1 + e_ * np.sin(phi))) / (2 * e_))
    qp = q(np.pi / 2); ra = a_ * np.sqrt(qp / 2)
    lat_top = t.f + t.e * np.arange(h); lat_bot = lat_top + t.e
    beta = lambda lat: np.arcsin(q(np.radians(lat)) / qp)
    return ra ** 2 * np.radians(abs(t.a)) * np.abs(np.sin(beta(lat_top)) - np.sin(beta(lat_bot)))


def main(src, year):
    mun = gpd.read_file(os.path.join(P.BND, 'SC_Municipios_2024.gpkg')).to_crs(4326)
    sc = mun.dissolve().geometry.iloc[0]
    minx, miny, maxx, maxy = sc.bounds
    with rasterio.open(src) as r:
        win = from_bounds(minx - 0.05, miny - 0.05, maxx + 0.05, maxy + 0.05, r.transform).round_offsets().round_lengths()
        win = win.intersection(rasterio.windows.Window(0, 0, r.width, r.height))
        out = {}
        for fact in (4, 2):
            h, w = int(win.height // fact), int(win.width // fact)
            a = r.read(1, window=win, out_shape=(h, w), resampling=Resampling.nearest)
            t = r.window_transform(win) * rasterio.Affine.scale(win.width / w, win.height / h)
            out[fact] = (a, t, geometry_mask([sc], out_shape=a.shape, transform=t, invert=True))
    a, t, inside = out[4]
    npz = os.path.join(EXT, f'mapbiomas_c10_lulc{year}_SC_120m.npz')
    np.savez_compressed(npz, a=a, inside=inside, transform=np.array(t)[:6])
    a, t, inside = out[2]
    area = np.repeat(authalic_row_areas(t, a.shape[0])[:, None], a.shape[1], axis=1)
    cls = np.unique(a[inside])
    unknown = sorted(set(cls.tolist()) - set(NAMES))
    assert not unknown, f'classes without a name: {unknown}'
    tab = pd.DataFrame({'class_id': cls, 'class_name': [NAMES[int(c)] for c in cls],
                        'area_km2': [area[inside & (a == c)].sum() / 1e6 for c in cls]})
    classified = tab[tab.class_id != 0].area_km2.sum()
    tab['pct_classified'] = np.where(tab.class_id != 0, 100 * tab.area_km2 / classified, np.nan)
    csv = os.path.join(EXT, f'mapbiomas_c10_lulc{year}_SC_class_areas.csv')
    tab.round({'area_km2': 2, 'pct_classified': 3}).to_csv(csv, index=False)
    print('wrote', npz, 'and', csv); print(tab.round(3).to_string(index=False))
    print('state area (km2):', round(tab.area_km2.sum(), 1), '| classified:', round(classified, 1))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
