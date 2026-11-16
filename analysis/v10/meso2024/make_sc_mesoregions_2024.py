#!/usr/bin/env python3
"""make_sc_mesoregions_2024.py - mesoregions of Santa Catarina as unions of the IBGE 2024 municipalities.

Why: the IBGE mesoregion polygons of data/boundaries/sc_mesoregions.gpkg leave out the coastal islands (Santa Catarina
Island, in Florianópolis, and São Francisco do Sul Island), whose 83 grid cells would then belong to no mesoregion.
Building each mesoregion from the 2024 municipal mesh puts the islands in their mesoregions, makes the mesoregions
coincide with the municipalities of Table 4, and uses the same polygons as the MapBiomas forest shares
(analysis/v10/maps/landcover_mesoregion.py).

Attribution: each municipality goes to the IBGE mesoregion polygon that holds the largest part of its area
(South America Albers equal-area, ESRI:102033). If the environment variable WILDFIRE_ATTR names a second attribution
table (CSV with CD_MUN and mesorregiao), the two attributions are compared and must agree; the shipped files were compared
with an attribution by the majority of 400 random points against another IBGE layer, with no disagreement (the .log of
this folder).

Inputs : data/boundaries/SC_Municipios_2024.gpkg, data/boundaries/sc_mesoregions.gpkg (IBGE mesoregions, islands left out),
         data/boundaries/landcover2024_mesoregion.csv (MapBiomas Collection 10, 2024, written by analysis/v10/maps/landcover_mesoregion.py;
         until 16 November 2026 forest_pct came from landcover2019_mesoregion.csv, 2019)
Outputs: reproduced/boundaries/sc_mesoregions_2024.gpkg (EPSG:4326; name_meso, n_mun, area_ibge_km2, forest_pct) and
         reproduced/boundaries/sc_municipios_2024_mesorregiao.csv (CD_MUN, NM_MUN, mesorregiao, majority_fraction);
         the copies shipped in data/boundaries were made by this script (set MEGAN_SC_BND_OUT=data/boundaries to overwrite them)
"""
import os, sys
import pandas as pd, geopandas as gpd

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
BND = P.BND; BOUT = os.environ.get('MEGAN_SC_BND_OUT') or P.out('boundaries')
mun = gpd.read_file(os.path.join(BND, 'SC_Municipios_2024.gpkg'))
old = gpd.read_file(os.path.join(BND, 'sc_mesoregions.gpkg'))
assert len(mun) == 295 and len(old) == 6
m = mun.to_crs('ESRI:102033'); o = old.to_crs('ESRI:102033')
inter = gpd.overlay(m[['CD_MUN', 'NM_MUN', 'geometry']], o[['name_meso', 'geometry']], how='intersection', keep_geom_type=True)
inter['a'] = inter.area
best = inter.sort_values('a').groupby('CD_MUN').tail(1).set_index('CD_MUN')
tab = mun[['CD_MUN', 'NM_MUN']].copy(); tab['CD_MUN'] = tab.CD_MUN.astype(str); best.index = best.index.astype(str)
tab['mesorregiao'] = tab.CD_MUN.map(best.name_meso)
tab['majority_fraction'] = tab.CD_MUN.map(best.a / m.set_index(m.CD_MUN.astype(str)).area).round(4)
assert tab.mesorregiao.notna().all()
wf = os.environ.get('WILDFIRE_ATTR', '')
if wf and os.path.exists(wf):
    w = pd.read_csv(wf, dtype={'CD_MUN': str}).set_index('CD_MUN').mesorregiao
    dis = tab[tab.mesorregiao != tab.CD_MUN.map(w)]
    print(f'compared with {wf}: {len(dis)} disagreements'); assert dis.empty, dis
tab.sort_values('CD_MUN').to_csv(os.path.join(BOUT, 'sc_municipios_2024_mesorregiao.csv'), index=False)
print(tab.mesorregiao.value_counts().to_string()); print('lowest majority fractions:'); print(tab.nsmallest(6, 'majority_fraction').to_string(index=False))

g = mun.assign(mesorregiao=mun.CD_MUN.astype(str).map(tab.set_index('CD_MUN').mesorregiao))
d = g.dissolve(by='mesorregiao', aggfunc={'CD_MUN': 'count', 'AREA_KM2': 'sum'}).reset_index()
d = d.rename(columns={'mesorregiao': 'name_meso', 'CD_MUN': 'n_mun', 'AREA_KM2': 'area_ibge_km2'})
# SC_Municipios_2024.gpkg is a simplified copy of the IBGE mesh (about 82 thousand vertices, against 2.5 million in the
# original shapefile), and each municipality was simplified on its own, so the union leaves thousands of small gaps between
# neighbours (inside the state: 6,850 holes, 22.9 km2 in all, the largest 0.44 km2). Interior rings smaller than 1 km2
# are those gaps and are filled; none of the IBGE municipalities encloses a larger area outside itself.
from shapely.geometry import Polygon, MultiPolygon
def fill_slivers(geom, max_km2=1.0):
    eq = gpd.GeoSeries([geom], crs='EPSG:4326').to_crs('ESRI:102033').iloc[0]
    parts_eq = list(eq.geoms) if isinstance(eq, MultiPolygon) else [eq]
    parts = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
    out, n = [], 0
    for p, pe in zip(parts, parts_eq):
        keep = [i for i, ie in zip(p.interiors, pe.interiors) if Polygon(ie).area / 1e6 >= max_km2]
        n += len(p.interiors) - len(keep); out.append(Polygon(p.exterior, keep))
    return (MultiPolygon(out) if len(out) > 1 else out[0]), n
d = d.to_crs('EPSG:4326')
filled = [fill_slivers(g) for g in d.geometry]
d['geometry'] = [g for g, _ in filled]
print('interior rings filled (< 1 km2) per mesoregion:', dict(zip(d.name_meso, [n for _, n in filled])))
LC_YEAR = 2024   # MapBiomas year of forest_pct: the last year of Collection 10 (2019 until 16 November 2026)
F = pd.read_csv(os.path.join(BND, f'landcover{LC_YEAR}_mesoregion.csv'))
d['forest_pct'] = d.name_meso.map(F[F.group == 'forest'].set_index('mesorregiao').pct_of_mesoregion)
assert d.forest_pct.notna().all()
d = d[['name_meso', 'n_mun', 'area_ibge_km2', 'forest_pct', 'geometry']].to_crs('EPSG:4326')
out = os.path.join(BOUT, 'sc_mesoregions_2024.gpkg')
import tempfile, shutil                     # SQLite (GPKG) is written locally and copied: some synced/mounted folders
tmp = os.path.join(tempfile.mkdtemp(), 'sc_mesoregions_2024.gpkg')   # do not support its locking
d.to_file(tmp, driver='GPKG'); shutil.copyfile(tmp, out)
d2 = d.to_crs('ESRI:102033')
print(pd.DataFrame({'name_meso': d.name_meso, 'n_mun': d.n_mun, 'area_ibge_km2': d.area_ibge_km2.round(1),
                    'polygon_km2': (d2.area / 1e6).round(1), 'forest_pct': d.forest_pct}).to_string(index=False))
