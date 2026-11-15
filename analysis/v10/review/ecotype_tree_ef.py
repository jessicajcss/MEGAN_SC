# -*- coding: utf-8 -*-
"""ecotype_tree_ef.py - the MEGAN v3 ecotypes of Santa Catarina and the tree emission factors they carry (Online Resource 1,
Table S8).

Reads only files distributed with the repository:
  data/sc_ecotype.nc                                    ecotype map clipped from EVT3d.nc (1/120 degree)
  data/megan_inputs/EFP240120/SpeciationTree230921.csv  tree composition of each ecotype (VegID, treeSpecfrac)
  data/megan_inputs/EFP240120/EFv240120.csv             emission factors EF1-EF19 of each VegID (nmol m-2 s-1)
  data/boundaries/sc_state_ibge.gpkg                    IBGE state polygon
The ecotype map and the two tables come from the MEGAN 3.2 input files (Guenther 2024, Zenodo,
https://doi.org/10.5281/zenodo.10939297).
Output: reproduced/ijest_review/ecotype_tree_ef.csv (published copy in results/v10/ijest_review): for each ecotype present
in the state, its share of the state area (ecotype cells with centre inside the polygon), the tree entries and their
fractions, and the composed tree emission factor (sum of fraction x EF) of isoprene (EF1), MT_PINE (EF3) and MBO (EF2), with
the same factors of the two entries. Code 65535 (no ecotype; a coastal strip that the run fills with the nearest valid
ecotype, Section 2.4 of the article) is listed with its share only.
"""
import os, sys
import numpy as np, pandas as pd, xarray as xr, geopandas as gpd, shapely

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
sys.path.insert(0, os.path.join(P.REPO, 'src'))
from megan_sc import efproc

D = os.path.join(P.REPO, 'data', 'megan_inputs', 'EFP240120')
ds = xr.open_dataset(os.path.join(P.REPO, 'data', 'sc_ecotype.nc')); evt = ds.EVT.values
LO, LA = np.meshgrid(ds.lon.values, ds.lat.values)
state = gpd.read_file(os.path.join(P.BND, 'sc_state_ibge.gpkg')).to_crs(4326).geometry.union_all()
inside = shapely.contains_xy(state, LO.ravel(), LA.ravel()).reshape(evt.shape)
ids, counts = np.unique(evt[inside], return_counts=True)
tab = efproc.load_ef_table(os.path.join(D, 'EFv240120.csv'))
tree = pd.read_csv(os.path.join(D, 'SpeciationTree230921.csv'), encoding='latin-1')
rows = []
for e, c in zip(ids, counts):
    r = dict(ecotype=int(e), state_area_pct=100 * c / counts.sum())
    sub = tree[tree.EcoTypeID == e]
    if len(sub):
        r['tree_entries'] = '; '.join(f'{v} {f:.3f}' for v, f in zip(sub.VegID, sub.treeSpecfrac))
        r['tree_entry_types'] = '; '.join(str(tab.loc[v, 'CommonName']) for v in sub.VegID)
        for k, nm in (('EF1', 'isop'), ('EF3', 'mt_pine'), ('EF2', 'mbo')):
            ef = np.array([float(tab.loc[v, k]) for v in sub.VegID]); fr = sub.treeSpecfrac.to_numpy(float)
            r[f'{nm}_entries'] = '; '.join(f'{x:g}' for x in ef); r[f'{nm}_tree'] = float((fr * ef).sum() / fr.sum())
        refs = [str(tab.loc[v, 'References']) for v in sub.VegID]
        r['references_field'] = '; '.join('(empty)' if x in ('nan', '') else x for x in refs)
    rows.append(r)
out = P.out('ijest_review', 'ecotype_tree_ef.csv')
pd.DataFrame(rows).to_csv(out, index=False, float_format='%.6g')
print(pd.DataFrame(rows).to_string()); print('written', out)
