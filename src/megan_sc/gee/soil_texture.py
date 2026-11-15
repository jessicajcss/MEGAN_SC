"""MapBiomas **Solo** (Soil) granulometry -> a texture-aware wilting point for GAMSM.

OPTIONAL and OFF by default (`datasets.soil.wilting_point_source: mapbiomas-texture`). The
soil-moisture activity factor GAMSM = f(theta, theta_wilt) uses a wilting point that depends on
soil texture; by default we use MEGAN's constant/soil-type value. This module lets that wilting
point instead come from **MapBiomas Soil Collection 3** sand/clay fractions (30 m), mapped
through the USDA texture triangle to MEGAN's own ``WWLT`` table
(:func:`megan_sc.engine.soil.wilting_point_from_texture`) — so it stays faithful to MEGAN (no
fitted pedotransfer coefficients).

Impact is modest in humid SC (soil rarely below wilting point, so GAMSM ≈ 1 most of the year);
it matters mainly in dry spells / the drier west. Kept behind a flag for that reason.

**Asset note (confirm on first use).** MapBiomas Solo granulometry is published on Earth Engine
under ``projects/mapbiomas-public/assets/brazil/soil/collection3/…`` (see
https://brasil.mapbiomas.org/en/dados-mapbiomas-solo/). The exact asset IDs, band names and
depth handling can change between collections, so they are **config parameters**
(``datasets.soil.clay_asset`` / ``sand_asset`` / ``texture_band``) rather than hard-coded, and
should be checked against the current MapBiomas Soil GEE toolkit before a production run.
"""
from __future__ import annotations

# Defaults from the MapBiomas Soil download page (Collection 3). Override in the config.
CLAY_ASSET_DEFAULT = ("projects/mapbiomas-public/assets/brazil/soil/collection3/"
                      "mapbiomas_brazil_collection3_soil_clay_fraction_v1")
SAND_ASSET_DEFAULT = ("projects/mapbiomas-public/assets/brazil/soil/collection3/"
                      "mapbiomas_brazil_collection3_soil_sand_fraction_v1")


def _first_band(ee, asset, band=None):
    """Load an asset as an Image (or the first image of a collection) and pick a band."""
    obj = ee.data.getAsset(asset) if hasattr(ee.data, "getAsset") else None
    typ = (obj or {}).get("type", "")
    img = ee.ImageCollection(asset).first() if typ == "IMAGE_COLLECTION" else ee.Image(asset)
    return img.select([band]) if band else img.select(0)


def sand_clay_image(ee, region, clay_asset=CLAY_ASSET_DEFAULT, sand_asset=SAND_ASSET_DEFAULT,
                    band=None):
    """ee.Image with 'sand' and 'clay' bands (surface-layer %) clipped to ``region``.

    ``band`` selects the depth/attribute band if the asset is multi-band (e.g. the 0-10 cm
    layer); leave None to take the first band. Values are expected as percentages (0-100).
    """
    clay = _first_band(ee, clay_asset, band).rename("clay")
    sand = _first_band(ee, sand_asset, band).rename("sand")
    return ee.Image.cat([sand, clay]).clip(region).toFloat()


def export_sand_clay_to_domain(ee, region, domain, folder,
                               clay_asset=CLAY_ASSET_DEFAULT, sand_asset=SAND_ASSET_DEFAULT,
                               band=None, max_pixels=65536):
    """Export sand/clay aggregated to the analysis grid, to Drive (server-side, no hang).

    Mirrors ``mapbiomas.export_growth_forms_to_domain``: 30 m soil → analysis grid via
    reduceResolution, onto the exact grid (so it aligns with met/LAI). Caller starts/polls the
    task, then reads with :func:`read_sand_clay_geotiff`.
    """
    from . import mapbiomas
    res_m = float(domain["resolution_deg"]) * 111000.0
    img = (sand_clay_image(ee, region, clay_asset, sand_asset, band)
           .reduceResolution(ee.Reducer.mean(), maxPixels=max_pixels)
           .reproject(crs="EPSG:4326", scale=res_m)
           .toFloat())
    return ee.batch.Export.image.toDrive(
        image=img, description="megan_soiltex", folder=folder,
        crs="EPSG:4326", crsTransform=mapbiomas.domain_crs_transform(domain),
        dimensions=f"{int(domain['lon'].size)}x{int(domain['lat'].size)}",
        maxPixels=int(1e13), fileNamePrefix="megan_soiltex",
    )


def read_sand_clay_geotiff(path, domain):
    """Read the exported sand/clay GeoTIFF -> dict {'sand','clay'} on the analysis grid."""
    import numpy as np
    import rasterio

    with rasterio.open(path) as src:
        arr = src.read().astype("float64")            # bands: sand, clay
    ny, nx = domain["lat"].size, domain["lon"].size
    assert arr.shape[1:] == (ny, nx), f"GeoTIFF {arr.shape[1:]} != grid {(ny, nx)}"
    return {"sand": np.flipud(arr[0]), "clay": np.flipud(arr[1])}


def wilting_point_grid(sand_pct, clay_pct):
    """Wilting-point grid (m^3/m^3) from sand/clay % — USDA class → MEGAN WWLT table."""
    from ..engine.soil import wilting_point_from_texture
    return wilting_point_from_texture(sand_pct, clay_pct)
