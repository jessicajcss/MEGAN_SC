"""Write a minimal, standards-correct IOAPI (Lambert Conformal) NetCDF for testing.

Structural stand-in for Hoinaski's MEGAN emission files so io/ioapi.py can be
regression-tested before the real file is downloaded. Grid attributes follow the
Models-3 IOAPI convention; two species (ISOP, TERP) are written with a known spatial
pattern and explicit ``units`` so tests can verify both the grid projection math and the
header parsing.
"""
from __future__ import annotations

import numpy as np


def write(path, ncols=25, nrows=22, xcell=4000.0, ycell=4000.0,
          xcent=-51.0, ycent=-27.0, p_alp=-10.0, p_bet=-40.0, ntime=2):
    from netCDF4 import Dataset
    from pyproj import Transformer

    # Put the grid's SW corner near SC: convert a lon/lat anchor to LCC metres.
    fwd = Transformer.from_crs(
        f"+proj=lcc +lat_1={p_alp} +lat_2={p_bet} +lat_0={ycent} +lon_0={xcent} "
        f"+x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m", "EPSG:4326", always_xy=True)
    inv = Transformer.from_crs(
        "EPSG:4326",
        f"+proj=lcc +lat_1={p_alp} +lat_2={p_bet} +lat_0={ycent} +lon_0={xcent} "
        f"+x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m", always_xy=True)
    x_sw, y_sw = inv.transform(-53.8, -29.4)      # SW-ish corner of SC

    ds = Dataset(path, "w", format="NETCDF4")
    ds.createDimension("TSTEP", None)
    ds.createDimension("DATE-TIME", 2)
    ds.createDimension("LAY", 1)
    ds.createDimension("VAR", 2)
    ds.createDimension("ROW", nrows)
    ds.createDimension("COL", ncols)

    for k, v in dict(GDTYP=2, NCOLS=ncols, NROWS=nrows, NLAYS=1, NVARS=2,
                     XORIG=float(x_sw), YORIG=float(y_sw), XCELL=xcell, YCELL=ycell,
                     XCENT=xcent, YCENT=ycent, P_ALP=p_alp, P_BET=p_bet, P_GAM=xcent,
                     SDATE=2019015, STIME=0, TSTEP=10000).items():
        setattr(ds, k, v)
    ds.setncattr("VAR-LIST", "ISOP            TERP            ")

    tflag = ds.createVariable("TFLAG", "i4", ("TSTEP", "VAR", "DATE-TIME"))
    for t in range(ntime):
        tflag[t, :, 0] = 2019015
        tflag[t, :, 1] = t * 10000                 # 00:00, 01:00, ...

    # cell-centre lon/lat to build a recognisable pattern
    x = float(x_sw) + (np.arange(ncols) + 0.5) * xcell
    y = float(y_sw) + (np.arange(nrows) + 0.5) * ycell
    X, Y = np.meshgrid(x, y)
    lon, lat = fwd.transform(X, Y)

    for name, patt in (("ISOP", 5.0 + 0.0 * lon), ("TERP", 1.0 + 0.0 * lon)):
        v = ds.createVariable(name, "f4", ("TSTEP", "LAY", "ROW", "COL"))
        v.units = "moles/s"
        v.long_name = name.ljust(16)
        v.var_desc = (name + " emission").ljust(80)
        v[:] = np.broadcast_to(patt, (ntime, 1, nrows, ncols))
    ds.close()
    return path


if __name__ == "__main__":
    import sys
    write(sys.argv[1] if len(sys.argv) > 1 else "synthetic_ioapi.nc")
    print("wrote it")
