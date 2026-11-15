"""Write analysis-grid arrays as georeferenced **GeoTIFF** (EPSG:4326) for QGIS.

Every input and output layer can be dumped to a GeoTIFF that loads directly in QGIS and lines up
with a municipality shapefile for zonal statistics. The domain grid is south-first
(lat ascending); GeoTIFF is north-up (row 0 = north), so arrays are flipped on write.
"""
from __future__ import annotations

import numpy as np


def _transform(lon, lat):
    """North-up affine for the analysis grid (top-left corner + pixel size)."""
    from rasterio.transform import Affine
    lon = np.asarray(lon, dtype="float64"); lat = np.asarray(lat, dtype="float64")
    rx = float(lon[1] - lon[0]) if lon.size > 1 else 0.025
    ry = float(lat[1] - lat[0]) if lat.size > 1 else 0.025
    west = float(lon[0]) - rx / 2.0
    north = float(lat[-1]) + ry / 2.0
    return Affine(rx, 0.0, west, 0.0, -ry, north)


def write_geotiff(path, arr, lon, lat, band_names=None, nodata=None):
    """Write a 2-D (ny,nx) or 3-D (nbands,ny,nx) array as a GeoTIFF on the analysis grid.

    ``band_names`` (optional) are written as band descriptions so QGIS shows them. Returns
    ``path``. CRS is EPSG:4326; the array is flipped to north-up so it overlays a shapefile.
    """
    import rasterio

    a = np.asarray(arr, dtype="float32")
    if a.ndim == 2:
        a = a[None, :, :]
    nbands, ny, nx = a.shape
    assert ny == len(lat) and nx == len(lon), "array/grid mismatch"
    a = a[:, ::-1, :]                                   # south-first -> north-up
    profile = dict(driver="GTiff", height=ny, width=nx, count=nbands, dtype="float32",
                   crs="EPSG:4326", transform=_transform(lon, lat), compress="deflate")
    if nodata is not None:
        profile["nodata"] = nodata
    with rasterio.open(path, "w", **profile) as dst:
        for i in range(nbands):
            dst.write(a[i], i + 1)
            if band_names is not None and i < len(band_names):
                dst.set_band_description(i + 1, str(band_names[i]))
    return path


def write_species_stack(path, arr, lon, lat, species_names, nodata=-9999.0):
    """Convenience: write a (n_species, ny, nx) emission stack with species as band names."""
    return write_geotiff(path, arr, lon, lat, band_names=list(species_names), nodata=nodata)
