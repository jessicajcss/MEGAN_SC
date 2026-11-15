"""Read Models-3 **IOAPI** NetCDF (the format of Hoinaski's MEGAN / CMAQ emission files).

Used only for *validation* — to line our biogenic fields up against the reference
gridded MEGAN output on its own grid. IOAPI is a stable, fully-documented convention
(https://www.cmascenter.org/ioapi/), so this reader is built to that spec and unit-tested
with a synthetic Lambert-grid file. The one thing that varies between products and can NOT
be guessed — the exact species-variable names and their ``units`` — is read straight from
the file header (see :func:`describe`), never hard-coded. Run ``scripts/inspect_reference.py``
on the real file first and let the numbers speak.

Grid geometry
-------------
IOAPI stores the grid in global attributes, not as coordinate variables:
``GDTYP`` (2 = Lambert Conformal), ``XORIG/YORIG`` (SW-corner, projection metres),
``XCELL/YCELL`` (cell size), ``XCENT/YCENT`` (projection origin lon/lat),
``P_ALP/P_BET`` (the two standard parallels). Cell (col, row) centre is
``x = XORIG + (col+0.5)*XCELL``, ``y = YORIG + (row+0.5)*YCELL``; we project those to
lon/lat with pyproj on IOAPI's default 6 370 000 m sphere.
"""
from __future__ import annotations

import numpy as np

# IOAPI reserved variable names that are not chemical species.
_NON_SPECIES = {"TFLAG", "ETFLAG", "LATD", "LOND"}


def _proj4(attrs):
    """Build the pyproj CRS string for the file's grid type."""
    g = int(attrs["GDTYP"])
    a = 6370000.0                              # IOAPI default sphere radius (m)
    if g == 1:                                 # geographic lat/lon (degrees)
        return "+proj=longlat +a=%f +b=%f" % (a, a)
    if g == 2:                                 # Lambert Conformal Conic
        return ("+proj=lcc +lat_1=%(P_ALP)f +lat_2=%(P_BET)f +lat_0=%(YCENT)f "
                "+lon_0=%(XCENT)f +x_0=0 +y_0=0 +a=%(a)f +b=%(a)f +units=m"
                % {**attrs, "a": a})
    if g == 7:                                 # (equatorial) Mercator
        return ("+proj=merc +lat_ts=%(P_ALP)f +lon_0=%(XCENT)f +x_0=0 +y_0=0 "
                "+a=%(a)f +b=%(a)f +units=m" % {**attrs, "a": a})
    if g == 6:                                 # polar stereographic
        return ("+proj=stere +lat_0=%(P_BET)f +lat_ts=%(P_BET)f +lon_0=%(P_GAM)f "
                "+x_0=0 +y_0=0 +a=%(a)f +b=%(a)f +units=m" % {**attrs, "a": a})
    raise NotImplementedError(f"IOAPI GDTYP={g} not supported")


def _global_attrs(ds):
    out = {}
    for k in ("GDTYP", "NCOLS", "NROWS", "NLAYS", "NVARS",
              "XORIG", "YORIG", "XCELL", "YCELL", "XCENT", "YCENT",
              "P_ALP", "P_BET", "P_GAM", "SDATE", "STIME", "TSTEP"):
        if hasattr(ds, k):
            out[k] = getattr(ds, k)
    return out


def cell_centers_lonlat(attrs):
    """2-D (NROWS, NCOLS) longitude/latitude of every grid-cell centre."""
    from pyproj import Transformer

    ncols, nrows = int(attrs["NCOLS"]), int(attrs["NROWS"])
    x = attrs["XORIG"] + (np.arange(ncols) + 0.5) * attrs["XCELL"]
    y = attrs["YORIG"] + (np.arange(nrows) + 0.5) * attrs["YCELL"]
    X, Y = np.meshgrid(x, y)                                    # (row, col)
    tr = Transformer.from_crs(_proj4(attrs), "EPSG:4326", always_xy=True)
    lon, lat = tr.transform(X, Y)
    return np.asarray(lon), np.asarray(lat)


def species_vars(ds):
    """List of chemical-species variable names present in the file."""
    return [v for v in ds.variables if v not in _NON_SPECIES
            and ds.variables[v].ndim == 4]     # (TSTEP, LAY, ROW, COL)


def describe(path):
    """Return a dict summarising the file: grid attrs, species (+units), time range.

    This is the honest first step — call it (or ``scripts/inspect_reference.py``) on the
    real reference file so the species names and units are *read*, never assumed.
    """
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        attrs = _global_attrs(ds)
        sp = species_vars(ds)
        units = {s: getattr(ds.variables[s], "units", "").strip() for s in sp}
        times = read_tflag(ds)
        return {"grid": attrs, "species": sp, "units": units,
                "ntimes": len(times), "times": times}


def read_tflag(ds):
    """Decode TFLAG into a list of (YYYYDDD, HHMMSS) per timestep (first variable)."""
    if "TFLAG" not in ds.variables:
        return []
    tf = np.asarray(ds.variables["TFLAG"][:])          # (TSTEP, VAR, 2)
    return [(int(tf[t, 0, 0]), int(tf[t, 0, 1])) for t in range(tf.shape[0])]


def read_species(path, names=None, tstep=None, layer=0):
    """Read species arrays from an IOAPI file.

    Parameters
    ----------
    names : list[str] | None
        Species to read (default: all). Case-insensitive; unknown names are skipped
        with the available list returned under key ``_available``.
    tstep : int | None
        Time index (default: sum over all timesteps → a daily total in the file's units
        times the number of steps; use an explicit index for an instantaneous field).
    layer : int
        Vertical layer (default 0, the surface).

    Returns
    -------
    dict: species_name -> 2-D (NROWS, NCOLS) array, plus ``_lon``/``_lat`` (2-D centres),
    ``_units`` (dict), ``_available`` (all species in the file).
    """
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        attrs = _global_attrs(ds)
        avail = species_vars(ds)
        upper = {s.upper(): s for s in avail}
        want = avail if names is None else [upper[n.upper()] for n in names if n.upper() in upper]

        lon, lat = cell_centers_lonlat(attrs)
        out = {"_lon": lon, "_lat": lat, "_available": avail,
               "_units": {s: getattr(ds.variables[s], "units", "").strip() for s in want}}
        for s in want:
            a = np.asarray(ds.variables[s][:], dtype="float64")   # (TSTEP, LAY, ROW, COL)
            out[s] = a[tstep, layer] if tstep is not None else a[:, layer].sum(axis=0)
        return out
