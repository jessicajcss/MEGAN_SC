"""Shared helpers of the analysis scripts of analysis/v10: GeoPackage reader, model grid, state mask, cell areas and
polygon rasterisation on the south-first model grid. The grid is read from the production run (data/run_v10, or
$MEGAN_SC_RUN; see _paths.py)."""
import os, sqlite3, struct
import numpy as np
from shapely import wkb
import _paths as P

BND = P.BND
CLASSES = ['ISOP','MBO','MT_PINE','MT_ACYC','MT_CAMP','MT_SABI','MT_AROM','NO','SQT_HR','SQT_LR',
           'MEOH','ACTO','ETOH','ACID','LVOC','OXPROD','STRESS','OTHER','CO']
BVOC = [c for c in CLASSES if c not in ('NO', 'CO')]
LDF1 = ['ISOP', 'MBO', 'CO']          # classes with LDF = 1 (purely light-dependent) in MEGVEA.EXT
# Canopy loss and production factor rho of MEGAN (E = epsilon * gamma * rho; Guenther et al., 2006, Sect. 3.3): 0.96 for
# isoprene; the MEGANv3.21 code, and so the output of notebook 10, applies none (rho = 1). Applied when the run is read.
RHO = {'ISOP': 0.96}

def rho(cls):
    return RHO.get(cls, 1.0)

def read_gpkg(path):
    """Minimal GeoPackage reader -> list of (attribute dict, shapely geometry)."""
    con = sqlite3.connect(path)
    tname = con.execute("select table_name from gpkg_contents where data_type='features'").fetchone()[0]
    gcol = con.execute("select column_name from gpkg_geometry_columns where table_name=?", (tname,)).fetchone()[0]
    cols = [r[1] for r in con.execute(f"pragma table_info('{tname}')")]
    out = []
    for row in con.execute(f"select * from '{tname}'"):
        rec = dict(zip(cols, row)); blob = rec.pop(gcol)
        flags = blob[3]; env = (flags >> 1) & 0b111
        hdr = 8 + {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]
        out.append((rec, wkb.loads(bytes(blob[hdr:]))))
    return out

def grid():
    """Cell-centre latitudes (south first) and longitudes of the model grid, from an annual file of the production run."""
    import xarray as xr
    with xr.open_dataset(os.path.join(P.RUN, '2025', 'emissions', 'annual_emis_2025.nc')) as ds:
        return ds.lat.values, ds.lon.values

def transform_for(lat, lon):
    from rasterio.transform import Affine
    rx = float(lon[1] - lon[0]); ry = float(lat[1] - lat[0])
    return Affine(rx, 0.0, float(lon[0]) - rx / 2, 0.0, -ry, float(lat[-1]) + ry / 2)

def rasterize_ids(geoms_ids, lat, lon, all_touched=False):
    """Rasterize (geom, id) on the south-first model grid; cell-centre rule by default."""
    from rasterio.features import rasterize
    r = rasterize(geoms_ids, out_shape=(lat.size, lon.size), transform=transform_for(lat, lon),
                  fill=0, dtype='int32', all_touched=all_touched)
    return np.flipud(r)

def cell_area_m2(lat, lon):
    R = 6371008.8; d = np.deg2rad(lon[1] - lon[0]); h = (lat[1] - lat[0]) / 2
    a = R**2 * d * np.abs(np.sin(np.deg2rad(lat + h)) - np.sin(np.deg2rad(lat - h)))
    return a[:, None] * np.ones((1, lon.size))

def area_weights(geoms, lat, lon):
    """Sparse matrix (polygons x cells, cells south-first row by row) of the fraction of each grid cell inside each polygon
    (area weighting, as in notebook 13). The fraction is taken in longitude-latitude; inside a 0.025° cell it differs from
    the equal-area fraction by less than 0.02 %."""
    import shapely
    from scipy import sparse
    res = float(lon[1] - lon[0]); ny, nx = lat.size, lon.size
    x0 = np.broadcast_to(lon[None, :] - res / 2, (ny, nx)); y0 = np.broadcast_to(lat[:, None] - res / 2, (ny, nx))
    cells = shapely.box(x0.ravel(), y0.ravel(), (x0 + res).ravel(), (y0 + res).ravel()); tree = shapely.STRtree(cells)
    r_, c_, v_ = [], [], []
    for k, g in enumerate(geoms):
        idx = tree.query(g, predicate='intersects'); f = shapely.area(shapely.intersection(cells[idx], g)) / res ** 2
        keep = f > 0; r_ += [k] * int(keep.sum()); c_ += list(idx[keep]); v_ += list(f[keep])
    return sparse.csr_matrix((v_, (r_, c_)), shape=(len(geoms), ny * nx))

def sc_mask(lat, lon):
    g = read_gpkg(f"{BND}/sc_state_ibge.gpkg")
    return rasterize_ids([(gg, 1) for _, gg in g], lat, lon) == 1
