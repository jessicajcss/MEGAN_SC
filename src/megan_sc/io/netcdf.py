"""CF-compliant NetCDF output for MEGAN-SC emissions."""
from __future__ import annotations

import numpy as np


def write_emissions(path, arr, species_names, lon, lat, *, units,
                    long_name_prefix="Biogenic emission of", global_attrs=None):
    """Write a (n_species, n_lat, n_lon) emission stack as CF-NetCDF.

    One variable per species. lon/lat are 1-D cell-centre coordinates.
    """
    from netCDF4 import Dataset

    arr = np.asarray(arr, dtype=np.float32)
    ns, nlat, nlon = arr.shape
    assert nlat == len(lat) and nlon == len(lon), "grid mismatch"

    ds = Dataset(path, "w", format="NETCDF4")
    ds.Conventions = "CF-1.8"
    ds.title = "MEGAN-SC biogenic emissions"
    for k, v in (global_attrs or {}).items():
        setattr(ds, k, str(v))

    ds.createDimension("lat", nlat)
    ds.createDimension("lon", nlon)
    vlat = ds.createVariable("lat", "f8", ("lat",))
    vlat.units = "degrees_north"; vlat.standard_name = "latitude"; vlat[:] = np.asarray(lat)
    vlon = ds.createVariable("lon", "f8", ("lon",))
    vlon.units = "degrees_east"; vlon.standard_name = "longitude"; vlon[:] = np.asarray(lon)

    for i, name in enumerate(species_names):
        safe = str(name).replace(" ", "_")
        v = ds.createVariable(safe, "f4", ("lat", "lon"), zlib=True, complevel=4,
                              fill_value=np.float32(-9999.0))
        v.units = units
        v.long_name = f"{long_name_prefix} {name}"
        v.coordinates = "lon lat"
        v[:] = arr[i]
    ds.close()
    return path
