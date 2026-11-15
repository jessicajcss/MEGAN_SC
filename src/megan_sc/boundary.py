"""Load (or one-time fetch) the Santa Catarina state boundary polygon.

The SC state polygon is needed to clip fields to the *actual* state for the state-inventory
comparison — the WRF d02 / MEGAN reference grid covers far more than SC (much of southern
Brazil), so a lon/lat bounding box overcounts SC's ~95,000 km² by 2-3x.

Offline-first: if the GeoPackage already exists it is loaded and no network is used. If it is
missing, :func:`fetch_sc_boundary` pulls the official IBGE boundary via ``geobr`` (needs
internet — run once on a networked machine; the sandbox proxy blocks it) and caches it.
"""
from __future__ import annotations

import os

DEFAULT_PATH = "data/boundaries/sc_state_ibge.gpkg"
SC_CODE = 42            # IBGE UF code for Santa Catarina


def load_sc_geometry(path=DEFAULT_PATH):
    """Return the SC boundary as a single (multi)polygon in EPSG:4326.

    Loads ``path`` if present; otherwise tries a one-time geobr fetch and caches it.
    """
    import geopandas as gpd

    if not os.path.exists(path):
        fetch_sc_boundary(path)
    gdf = gpd.read_file(path).to_crs("EPSG:4326")
    return gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") \
        else gdf.geometry.unary_union


def fetch_sc_boundary(path=DEFAULT_PATH, year=2020):
    """One-time download of the IBGE SC polygon via geobr; save as GeoPackage.

    Raises a clear error (not a silent fallback) if there is no network, so the caller
    knows to run it where internet is available.
    """
    try:
        import geobr
    except ImportError as e:
        raise RuntimeError("pip install geobr to fetch the SC boundary") from e
    try:
        sc = geobr.read_state(code_state="SC", year=year).to_crs("EPSG:4326")
    except Exception as e:
        raise RuntimeError(
            "Could not download the SC boundary (no internet / blocked proxy). "
            "Run this once on a networked machine, or drop an SC state polygon at "
            f"{path!r}.") from e
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    sc.to_file(path, driver="GPKG")
    return path
