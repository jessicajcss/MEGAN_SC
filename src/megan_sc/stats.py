"""Descriptive statistics, per-municipality zonal stats, and run checkpointing.

Pure-numpy where possible (testable); geopandas/rasterio only inside the zonal helpers.
Used by the production inventory notebook to (a) summarise each species/month, (b) aggregate
to the IBGE municipality shapefile for QGIS-style descriptive stats, and (c) checkpoint the
month accumulator so a Colab interruption loses nothing.
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# Descriptive statistics
# --------------------------------------------------------------------------- #
def field_descriptive(field, mask=None):
    """min / mean / max / std / sum / n over finite (and optionally masked) cells."""
    a = np.asarray(field, dtype="float64")
    m = np.isfinite(a)
    if mask is not None:
        m &= np.asarray(mask, dtype=bool)
    v = a[m]
    if v.size == 0:
        return dict(n=0, sum=0.0, mean=np.nan, min=np.nan, max=np.nan, std=np.nan)
    return dict(n=int(v.size), sum=float(v.sum()), mean=float(v.mean()),
                min=float(v.min()), max=float(v.max()), std=float(v.std()))


def species_stats_rows(stack_tonnes, names, year, month, lon=None, lat=None, mask=None):
    """One descriptive-stats row per species (monthly tonnes stack -> list of dicts)."""
    rows = []
    for i, nm in enumerate(names):
        s = field_descriptive(stack_tonnes[i], mask=mask)
        rows.append({"year": year, "month": month, "species": str(nm),
                     "total_tonnes": s["sum"], "cells": s["n"],
                     "max_tonnes_cell": s["max"], "mean_tonnes_cell": s["mean"]})
    return rows


# --------------------------------------------------------------------------- #
# Per-municipality zonal statistics
# --------------------------------------------------------------------------- #
def rasterize_municipalities(gdf, lon, lat, id_field):
    """Rasterise a municipality GeoDataFrame onto the analysis grid.

    Returns ``(zone, ids)`` where ``zone`` is a (ny,nx) int array of zone indices (0..K-1, or
    -1 outside any municipality) aligned to the south-first grid, and ``ids`` is the list of
    ``id_field`` values in zone-index order. Fast: one rasterisation, then group with bincount.
    """
    import geopandas as gpd  # noqa
    from rasterio.features import rasterize
    from rasterio.transform import Affine

    lon = np.asarray(lon, dtype="float64"); lat = np.asarray(lat, dtype="float64")
    rx = float(lon[1] - lon[0]); ry = float(lat[1] - lat[0])
    transform = Affine(rx, 0.0, float(lon[0]) - rx / 2, 0.0, -ry, float(lat[-1]) + ry / 2)
    ids = list(gdf[id_field].values)
    shapes = [(geom, i) for i, geom in enumerate(gdf.geometry.values)]
    north_up = rasterize(shapes, out_shape=(len(lat), len(lon)), transform=transform,
                         fill=-1, dtype="int32", all_touched=False)
    return north_up[::-1, :], ids                         # flip north-up -> south-first grid


def zonal_sums(field, zone, n_zones):
    """Per-zone sum of a field (zone from :func:`rasterize_municipalities`)."""
    a = np.asarray(field, dtype="float64").ravel()
    z = np.asarray(zone, dtype="int64").ravel()
    good = (z >= 0) & np.isfinite(a)
    return np.bincount(z[good], weights=a[good], minlength=n_zones)


def municipality_table(stack_tonnes, names, zone, ids, year, month):
    """Long-format per-municipality per-species monthly tonnes (list of dicts, for a DataFrame)."""
    n = len(ids)
    rows = []
    for i, nm in enumerate(names):
        sums = zonal_sums(stack_tonnes[i], zone, n)
        for k, mid in enumerate(ids):
            rows.append({"year": year, "month": month, "municipality": mid,
                         "species": str(nm), "tonnes": float(sums[k])})
    return rows


# --------------------------------------------------------------------------- #
# Checkpointing (survive Colab interruptions)
# --------------------------------------------------------------------------- #
def save_checkpoint(path, arrays, done_steps, meta=None):
    """Save named accumulator arrays + completed-timestep set to a .npz (atomic replace).

    ``arrays`` is a dict name→ndarray (arrays may have different shapes — e.g. tonnes cube,
    met cube, diurnal matrix). Written to a temp file then ``os.replace`` so an interrupt during
    the write can't corrupt the checkpoint.
    """
    import os
    tmp = str(path) + ".tmp.npz"
    np.savez_compressed(tmp, done=np.asarray(list(done_steps)),
                        meta=np.asarray(str(meta or "")),
                        **{k: np.asarray(v, dtype="float64") for k, v in arrays.items()})
    os.replace(tmp, path)
    return path


def load_checkpoint(path):
    """Load a checkpoint -> (arrays_dict, set_of_done_steps) or None if absent/unreadable."""
    import os
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path, allow_pickle=False)
        arrays = {k: z[k] for k in z.files if k not in ("done", "meta")}
        return arrays, set(int(x) for x in z["done"].tolist())
    except Exception:
        return None
