"""Compare MEGAN-SC output against a reference gridded inventory (Hoinaski MEGAN).

Pure-numpy comparison helpers — no Earth Engine, no I/O — so they are unit-testable and
reusable from both the notebook and ``scripts/compare_reference_day.py``. The philosophy
is *diagnosis, not fitting*: these functions quantify the gap
between our fields and the reference; they never rescale ours to match.

Two things they do:

1. **Unit-convert the reference to tonnes** for a magnitude check. MEGAN IOAPI biogenic
   files are per-cell emission *rates* (typically ``moles/s``); tonnes over a period =
   rate x molar-weight x seconds. The molar weights below are the standard ones; the
   ``units`` string read from the file selects the path (moles vs grams) so we never guess.
2. **Put both fields on one grid and score the spatial pattern** (bias, RMSE, correlation,
   ratio of domain totals). The reference (curvilinear Lambert) is interpolated onto our
   regular analysis grid for the map comparison; totals are computed on each native grid so
   the ratio is independent of the regridding.
"""
from __future__ import annotations

import numpy as np

# Molar weights (g/mol) for the species we compare. Monoterpenes lumped as C10H16.
MW = {
    "ISOP": 68.12, "ISOPRENE": 68.12,
    "TERP": 136.23, "MONOTERP": 136.23, "MT": 136.23,
    "SESQ": 204.35, "SQT": 204.35,
    "CO": 28.01, "NO": 30.01,
    "MEOH": 32.04, "ACET": 58.08, "CH3OH": 32.04,
}


def reference_tonnes(rate_per_cell, units, mw, seconds):
    """Convert a per-cell reference rate field to **tonnes** emitted over ``seconds``.

    ``units`` is the file's own units string. Supported: ``moles/s`` (or ``mol/s``,
    ``moles s-1``) → x mw → grams; ``g/s`` / ``g s-1`` → already grams. Anything else
    raises, rather than silently applying the wrong factor.
    """
    u = units.lower().replace(" ", "").replace("^", "")
    r = np.asarray(rate_per_cell, dtype="float64")
    if u in ("moles/s", "mol/s", "moles-1", "moless-1", "moles.s-1", "mol.s-1",
             "moless1", "mols-1"):
        grams = r * mw * seconds
    elif u in ("g/s", "gs-1", "g.s-1", "grams/s"):
        grams = r * seconds
    else:
        raise ValueError(f"unrecognised reference units {units!r}; add the conversion explicitly")
    return grams * 1e-6            # g -> tonnes


def regrid_to_regular(field, src_lon2d, src_lat2d, tgt_lon1d, tgt_lat1d, method="linear"):
    """Interpolate a curvilinear field onto a regular (tgt_lat x tgt_lon) grid.

    Linear where a Delaunay triangulation exists, nearest-neighbour to fill any holes
    outside the convex hull; if the source is degenerate (too few or collinear points)
    the whole field falls back to nearest rather than raising.
    """
    from scipy.interpolate import griddata
    from scipy.spatial import QhullError

    pts = np.column_stack([np.asarray(src_lon2d).ravel(), np.asarray(src_lat2d).ravel()])
    vals = np.asarray(field, dtype="float64").ravel()
    good = np.isfinite(vals)
    pts, vals = pts[good], vals[good]
    TLon, TLat = np.meshgrid(np.asarray(tgt_lon1d), np.asarray(tgt_lat1d))
    tgt = np.column_stack([TLon.ravel(), TLat.ravel()])
    try:
        out = griddata(pts, vals, tgt, method=method)
    except QhullError:
        out = griddata(pts, vals, tgt, method="nearest")
    holes = ~np.isfinite(out)
    if holes.any():
        out[holes] = griddata(pts, vals, tgt[holes], method="nearest")
    return out.reshape(TLat.shape)


def polygon_mask(lon2d, lat2d, geom):
    """Boolean mask (shape of lon2d) of cell centres inside a shapely polygon.

    Used to clip a field to the actual **Santa Catarina state** boundary — essential
    when comparing against the SC-only state inventory, since the reference d02 grid
    spans far more than SC (a lon/lat bounding box overcounts by ~2-3x). Uses the
    vectorised point-in-polygon test when available, else a prepared-geometry loop.
    """
    lon = np.asarray(lon2d, dtype="float64")
    lat = np.asarray(lat2d, dtype="float64")
    try:
        from shapely import vectorized
        return vectorized.contains(geom, lon, lat)
    except Exception:
        from shapely.geometry import Point
        from shapely.prepared import prep
        pg = prep(geom)
        out = np.zeros(lon.shape, dtype=bool)
        for idx in np.ndindex(lon.shape):
            out[idx] = pg.contains(Point(lon[idx], lat[idx]))
        return out


def bbox_mask(lon2d, lat2d, lon_min, lon_max, lat_min, lat_max):
    """Boolean mask of cell centres inside a lon/lat bounding box.

    Use this to clip the reference to *our* analysis footprint before summing tonnes:
    tonnes is an **extensive** quantity, so totals must be summed on each native grid over
    the same geographic area — never interpolated across resolutions and then summed.
    """
    lon = np.asarray(lon2d, dtype="float64")
    lat = np.asarray(lat2d, dtype="float64")
    return (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)


def masked_total(field, mask=None):
    """Sum of a field over True cells (all finite cells if mask is None)."""
    a = np.asarray(field, dtype="float64")
    m = np.isfinite(a)
    if mask is not None:
        m &= np.asarray(mask, dtype=bool)
    return float(a[m].sum())


def field_stats(ours, ref, mask=None):
    """Spatial-agreement stats between two fields on the *same* grid.

    Returns bias (ours-ref mean), rmse, Pearson correlation of the spatial pattern,
    and the ratio of sums (ours/ref). NaNs (and any ``mask==False`` cells) are excluded.
    """
    a = np.asarray(ours, dtype="float64")
    b = np.asarray(ref, dtype="float64")
    m = np.isfinite(a) & np.isfinite(b)
    if mask is not None:
        m &= np.asarray(mask, dtype=bool)
    a, b = a[m], b[m]
    if a.size == 0:
        return {"n": 0, "bias": np.nan, "rmse": np.nan, "corr": np.nan, "ratio": np.nan}
    bias = float(np.mean(a - b))
    rmse = float(np.sqrt(np.mean((a - b) ** 2)))
    corr = float(np.corrcoef(a, b)[0, 1]) if a.size > 1 and a.std() > 0 and b.std() > 0 else np.nan
    ratio = float(a.sum() / b.sum()) if b.sum() != 0 else np.nan
    return {"n": int(a.size), "bias": bias, "rmse": rmse, "corr": corr, "ratio": ratio}


def compare_species(our_map, our_lon, our_lat, our_area, ref_rate, ref_lon2d, ref_lat2d,
                    ref_units, ref_area_m2, mw, seconds, *, our_units="nanomol/m2/s"):
    """Full magnitude + pattern comparison for one species.

    ``our_map``   : our per-cell flux on our regular grid, in ``our_units``
                    (default nanomol m-2 s-1) — converted to tonnes with ``our_area``.
    ``ref_rate``  : reference per-cell rate on its native grid, in ``ref_units``.
    ``ref_area_m2``: reference cell area (scalar or array), used only if ref is a flux.
    Returns totals (tonnes) for both, the ratio, and spatial stats on our grid.
    """
    # our tonnes: nanomol/m2/s -> mol/m2/s (1e-9) -> mol/s (xarea) -> g (xmw) -> t (1e-6)
    ou = our_units.lower().replace(" ", "")
    if ou in ("nanomol/m2/s", "nanomol/m^2/s", "nanomolm-2s-1"):
        our_mol_s = np.asarray(our_map) * 1e-9 * np.asarray(our_area)
    elif ou in ("mol/m2/s", "molm-2s-1"):
        our_mol_s = np.asarray(our_map) * np.asarray(our_area)
    else:
        raise ValueError(f"unrecognised our_units {our_units!r}")
    our_tonnes = our_mol_s * mw * seconds * 1e-6

    ref_tonnes = reference_tonnes(ref_rate, ref_units, mw, seconds)

    ref_on_ours = regrid_to_regular(ref_tonnes, ref_lon2d, ref_lat2d, our_lon, our_lat)
    stats = field_stats(our_tonnes, ref_on_ours)
    return {
        "our_total_t": float(np.nansum(our_tonnes)),
        "ref_total_t": float(np.nansum(ref_tonnes)),
        "total_ratio": float(np.nansum(our_tonnes) / np.nansum(ref_tonnes))
        if np.nansum(ref_tonnes) else np.nan,
        "spatial": stats,
        "our_tonnes": our_tonnes,
        "ref_tonnes_on_our_grid": ref_on_ours,
    }
