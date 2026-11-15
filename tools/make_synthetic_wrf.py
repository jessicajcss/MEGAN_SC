"""Write a minimal, standards-correct ``wrfout``-like NetCDF for testing io/wrf.py.

This is NOT scientific data — it is a structural stand-in so the WRF reader can be
exercised (and regression-tested) before a real multi-GB ``wrfout`` is downloaded. It
reproduces exactly the pieces the reader touches:

* the ``Times`` char array (WRF's ``YYYY-MM-DD_HH:MM:SS`` format),
* a genuinely 2-D *curvilinear* coordinate pair (XLAT/XLONG) with a shear term, so the
  reader's triangulated regridding is actually tested (not a trivial rectilinear case),
* the met variables with WRF-ARW ranks: (Time, sn, we) for 2-D fields and
  (Time, soil_layers_stag, sn, we) for SMOIS.

Every met field is a **known linear function of (lon, lat)** so a test can assert the
interpolation is exact (linear fields are reproduced by barycentric interpolation) and
that no lat/lon axis got flipped.
"""
from __future__ import annotations

import numpy as np


# Linear field definitions f(lon, lat) — chosen so tests can predict interpolated values.
def temp_field(lon, lat):   return 295.0 + 0.5 * (lon + 52.0) - 0.8 * (lat + 27.0)   # K
def swdown_field(lon, lat): return 600.0 + 0.0 * lon + 0.0 * lat                     # W/m2 (const)
def u10_field(lon, lat):    return 3.0 + 0.1 * (lon + 52.0)                          # m/s
def v10_field(lon, lat):    return -4.0 + 0.05 * (lat + 27.0)                        # m/s
def psfc_field(lon, lat):   return 99000.0 + 20.0 * (lat + 27.0)                     # Pa
def q2_field(lon, lat):     return 0.012 + 0.0005 * (lon + 52.0)                     # kg/kg
def smois_field(lon, lat):  return 0.28 + 0.01 * (lat + 27.0)                        # m3/m3


# Vegetation fields, banded by longitude so tests can predict them:
#   west of -52  -> evergreen broadleaf forest (IVGTYP 2, tree), LAI 5
#   -52 .. -50.5 -> cropland (IVGTYP 12, crop),                  LAI 1.5
#   east of -50.5-> water (IVGTYP 17, masked),                   LAI 0
def ivgtyp_field(lon, lat):
    out = np.full(lon.shape, 12, dtype="i4")   # cropland default
    out[lon < -52.0] = 2                        # forest
    out[lon > -50.5] = 17                       # water
    return out


def lai_wrf_field(lon, lat):
    lai = np.full(lon.shape, 1.5)
    lai[lon < -52.0] = 5.0
    lai[lon > -50.5] = 0.0
    return lai


def xland_field(lon, lat):
    x = np.ones(lon.shape)
    x[lon > -50.5] = 2.0                         # 2 = water
    return x


def build_curvilinear(nx=20, ny=16, lon0=-54.5, lat0=-30.0, dlon=0.35, dlat=0.30,
                      shear=0.08):
    """A sheared (curvilinear) lon/lat mesh spanning well beyond the SC test box."""
    j, i = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    lon = lon0 + i * dlon + shear * j          # shear makes columns non-vertical -> 2-D XLONG
    lat = lat0 + j * dlat + 0.5 * shear * i     # and rows non-horizontal -> 2-D XLAT
    return lon.astype("f4"), lat.astype("f4")


def write(path, times=("2019-01-15_15:00:00", "2019-01-15_18:00:00"),
          nx=20, ny=16):
    from netCDF4 import Dataset

    lon2d, lat2d = build_curvilinear(nx=nx, ny=ny)
    nt = len(times)

    ds = Dataset(path, "w", format="NETCDF4")
    ds.TITLE = "SYNTHETIC WRF (test fixture — not real data)"
    ds.createDimension("Time", nt)
    ds.createDimension("DateStrLen", 19)
    ds.createDimension("south_north", ny)
    ds.createDimension("west_east", nx)
    ds.createDimension("soil_layers_stag", 4)

    vt = ds.createVariable("Times", "S1", ("Time", "DateStrLen"))
    for k, s in enumerate(times):
        vt[k] = np.array(list(s.ljust(19)[:19]), dtype="S1")

    def add2d(name, arr):
        v = ds.createVariable(name, "f4", ("Time", "south_north", "west_east"))
        v[:] = np.broadcast_to(arr, (nt, arr.shape[0], arr.shape[1]))

    add2d("XLAT", lat2d)
    add2d("XLONG", lon2d)
    add2d("T2", temp_field(lon2d, lat2d))
    add2d("SWDOWN", swdown_field(lon2d, lat2d))
    add2d("U10", u10_field(lon2d, lat2d))
    add2d("V10", v10_field(lon2d, lat2d))
    add2d("PSFC", psfc_field(lon2d, lat2d))
    add2d("Q2", q2_field(lon2d, lat2d))

    vs = ds.createVariable("SMOIS", "f4",
                           ("Time", "soil_layers_stag", "south_north", "west_east"))
    top = smois_field(lon2d, lat2d)
    layers = np.stack([top, top * 0.9, top * 0.8, top * 0.7])   # top layer is index 0
    vs[:] = np.broadcast_to(layers, (nt, 4, ny, nx))

    # vegetation fields (for read_wrf_vegetation — the no-Earth-Engine path)
    add2d("LAI", lai_wrf_field(lon2d, lat2d))
    add2d("VEGFRA", np.where(xland_field(lon2d, lat2d) < 1.5, 70.0, 0.0))
    add2d("XLAND", xland_field(lon2d, lat2d))
    vi = ds.createVariable("IVGTYP", "i4", ("Time", "south_north", "west_east"))
    vi[:] = np.broadcast_to(ivgtyp_field(lon2d, lat2d), (nt, ny, nx))
    ds.close()
    return path


if __name__ == "__main__":
    import sys
    out = sys.argv[1] if len(sys.argv) > 1 else "synthetic_wrfout.nc"
    write(out)
    print("wrote", out)
