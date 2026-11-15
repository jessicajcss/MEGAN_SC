"""Analysis domain & grid — user-defined by bounding box or shapefile.

Builds the regular lat/lon grid the whole pipeline runs on and computes
**geodesic** cell areas (m^2) — the quantity that turns per-area fluxes into
per-cell masses and tonnes. Areas vary with latitude and are NOT approximated by
Δlat·Δlon. A spherical (authalic-radius) formula is the fast default; an exact
ellipsoidal area via pyproj is available for validation.
"""
from __future__ import annotations

import numpy as np

# WGS84 authalic (equal-area) sphere radius (m) — matches ellipsoid surface area.
R_AUTHALIC = 6371007.1810


def make_axes(lon_min, lon_max, lat_min, lat_max, resolution_deg):
    """Cell-centred 1-D lon/lat axes for the bounding box at ``resolution_deg``."""
    res = float(resolution_deg)
    nlon = int(round((lon_max - lon_min) / res))
    nlat = int(round((lat_max - lat_min) / res))
    if nlon < 1 or nlat < 1:
        raise ValueError("resolution too coarse for the bounding box")
    lon = lon_min + (np.arange(nlon) + 0.5) * res
    lat = lat_min + (np.arange(nlat) + 0.5) * res
    return lon, lat


def cell_area_spherical(lat_centers, resolution_deg, radius=R_AUTHALIC):
    """Per-cell area (m^2) for each latitude band (1-D over lat).

    Exact on a sphere: A = R^2 · Δλ · (sin φ_north − sin φ_south).
    """
    res = float(resolution_deg)
    dlon = np.radians(res)
    lat_n = np.radians(np.asarray(lat_centers, float) + res / 2.0)
    lat_s = np.radians(np.asarray(lat_centers, float) - res / 2.0)
    return radius ** 2 * dlon * (np.sin(lat_n) - np.sin(lat_s))


def cell_area_grid(lat, resolution_deg, nlon, radius=R_AUTHALIC):
    """2-D (nlat, nlon) geodesic cell-area field (m^2)."""
    a1d = cell_area_spherical(lat, resolution_deg, radius)
    return np.repeat(a1d[:, None], nlon, axis=1)


def cell_area_ellipsoidal(lat_centers, resolution_deg):
    """Exact ellipsoidal per-cell area (m^2) per latitude band, via pyproj.Geod.

    Slower; use to validate the spherical areas (agreement within ~0.1-0.3%).
    """
    from pyproj import Geod

    geod = Geod(ellps="WGS84")
    res = float(resolution_deg)
    out = np.empty(len(lat_centers), dtype=float)
    for i, la in enumerate(np.asarray(lat_centers, float)):
        lons = [0.0, res, res, 0.0]
        lats = [la - res / 2, la - res / 2, la + res / 2, la + res / 2]
        area, _ = geod.polygon_area_perimeter(lons, lats)
        out[i] = abs(area)
    return out


def regrid_block_mean(field, flon, flat, tlon, tlat, res):
    """Area-agnostic block-mean regrid of a fine regular grid onto a coarser one.

    ``field`` is (..., ny_fine, nx_fine) on fine axes ``flon``/``flat``; returns
    (..., ny_target, nx_target) on target axes ``tlon``/``tlat`` (cell spacing
    ``res``). Each target cell = mean of the fine cells whose centres fall in it;
    empty target cells are NaN.
    """
    field = np.asarray(field, dtype=np.float64)
    lead = field.shape[:-2]
    ny_f, nx_f = field.shape[-2:]
    ix = np.clip(np.floor((np.asarray(flon) - (tlon[0] - res / 2)) / res).astype(int), 0, tlon.size - 1)
    iy = np.clip(np.floor((np.asarray(flat) - (tlat[0] - res / 2)) / res).astype(int), 0, tlat.size - 1)
    IX, IY = np.meshgrid(ix, iy)
    idx = (IY * tlon.size + IX).ravel()
    ncell = tlat.size * tlon.size
    cnt = np.bincount(idx, minlength=ncell)
    denom = np.where(cnt > 0, cnt, np.nan)
    stack = field.reshape(-1, ny_f, nx_f)
    out = np.empty((stack.shape[0], ncell), dtype=np.float64)
    for k in range(stack.shape[0]):
        s = np.bincount(idx, weights=stack[k].ravel(), minlength=ncell)
        with np.errstate(invalid="ignore"):
            out[k] = s / denom
    return out.reshape(lead + (tlat.size, tlon.size))


def build_domain(config):
    """Build the analysis grid from a config dict's ``domain`` block.

    Expects keys: ``mode`` ('bbox'|'shapefile'), ``resolution_deg``, and either
    ``bbox`` (lon_min/lon_max/lat_min/lat_max) or ``shapefile`` (path). When a
    shapefile is given, the grid spans its bounds and a boolean ``mask`` marks
    cells whose centre falls inside the polygon.

    Returns a dict: lon, lat (1-D), LON, LAT (2-D), area (2-D m^2), mask (2-D or None),
    plus resolution and bounds.
    """
    d = config["domain"]
    res = d["resolution_deg"]
    mask = None

    if d.get("mode") == "shapefile" and d.get("shapefile"):
        import geopandas as gpd
        from shapely.geometry import Point
        from shapely.prepared import prep

        gdf = gpd.read_file(d["shapefile"]).to_crs("EPSG:4326")
        geom = gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.geometry.unary_union
        lon_min, lat_min, lon_max, lat_max = geom.bounds
        lon, lat = make_axes(lon_min, lon_max, lat_min, lat_max, res)
        LON, LAT = np.meshgrid(lon, lat)
        pgeom = prep(geom)
        mask = np.array([[pgeom.contains(Point(x, y)) for x in lon] for y in lat], dtype=bool)
    else:
        b = d["bbox"]
        lon_min, lon_max = b["lon_min"], b["lon_max"]
        lat_min, lat_max = b["lat_min"], b["lat_max"]
        lon, lat = make_axes(lon_min, lon_max, lat_min, lat_max, res)
        LON, LAT = np.meshgrid(lon, lat)

    area = cell_area_grid(lat, res, lon.size)
    return {
        "lon": lon, "lat": lat, "LON": LON, "LAT": LAT,
        "area": area, "mask": mask, "resolution_deg": res,
        "bounds": (float(lon.min() - res / 2), float(lat.min() - res / 2),
                   float(lon.max() + res / 2), float(lat.max() + res / 2)),
    }
