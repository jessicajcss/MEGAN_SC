"""MapBiomas land cover -> MEGAN growth-form / canopy fractions, via Earth Engine.

Remaps the 30 m MapBiomas classification to per-pixel growth-form weights using
the documented crosswalk (megan_sc.landcover), then aggregates to the analysis
grid with reduceResolution. Produces tree/shrub/herb/crop fraction bands plus the
vegetated fraction — the high-resolution vegetation structure that weights the
ecotype emission factors and drives the canopy model.

Asset: ``projects/mapbiomas-public/assets/brazil/lulc/v1`` (band ``classification``, 30 m, one image per year and
collection; Collection 10 covers 1985-2024). Runs under the user's EE project.
"""
from __future__ import annotations

from .. import landcover

ASSET = "projects/mapbiomas-public/assets/brazil/lulc/v1"


def _remap(gf):
    codes = sorted(landcover.CROSSWALK)
    i = landcover._GF_IDX[gf]
    return codes, [float(landcover.CROSSWALK[c][i]) for c in codes]


def classification_image(ee, year, collection_id=None):
    """The MapBiomas ``classification`` image for a year (deterministic pick).

    The asset is an ImageCollection with one image per (year, collection_id, version).
    We filter by ``year`` and take the **newest collection** (sort collection_id desc) so
    a run is reproducible; pass ``collection_id`` to pin an exact collection.
    """
    col = ee.ImageCollection(ASSET).filter(ee.Filter.eq("year", int(year)))
    if collection_id is not None:
        col = col.filter(ee.Filter.eq("collection_id", collection_id))
    return ee.Image(col.sort("collection_id", False).first()).select("classification")


def resolve_collection_id(ee, year, collection=None):
    """The ``collection_id`` property value of the image to use for ``year``.

    ``collection=None`` returns None, so :func:`classification_image` takes the newest collection in the asset (the
    behaviour of the Santa Catarina production run, when Collection 10 was the newest). ``collection=10`` returns the
    newest release of Collection 10 (property values 10, 10.0, 10.1, ...), read from the asset itself so that the value
    keeps the type stored in Earth Engine. A run pinned this way does not change when a new collection is added.
    """
    if collection is None:
        return None
    col = ee.ImageCollection(ASSET).filter(ee.Filter.eq("year", int(year)))
    ids = [i for i in col.aggregate_array("collection_id").getInfo() if int(float(i)) == int(collection)]
    if not ids:
        raise ValueError(f"MapBiomas asset {ASSET} has no image of Collection {collection} for {year}")
    return max(ids, key=float)


def growth_form_fractions(ee, region, year, out_scale_m, max_pixels=65536,
                          collection_id=None):
    """ee.Image with tree/shrub/herb/crop fraction bands on the analysis grid.

    ``out_scale_m`` is the analysis-cell size in metres (e.g. ~2780 for 0.025°).
    Each band is the mean 30 m growth-form weight aggregated to that cell.

    ``max_pixels`` caps the number of native (30 m) pixels aggregated per output
    cell in ``reduceResolution``; it must be >= (out_scale_m / 30)^2. The default
    65536 (the EE maximum) covers analysis grids down to ~0.07° (~30x downscale on
    each axis → up to ~65k input pixels). For much coarser grids, or if interactive
    ``computePixels`` hits a memory limit, use :func:`export_year` (Drive) instead.
    """
    img = classification_image(ee, year, collection_id).clip(region)
    bands = []
    for gf in landcover.GROWTH_FORMS:
        codes, weights = _remap(gf)
        w = img.remap(codes, weights, 0).rename(gf).toFloat()
        wr = (w.reduceResolution(ee.Reducer.mean(), maxPixels=max_pixels)
              .reproject(crs="EPSG:4326", scale=out_scale_m))
        bands.append(wr)
    out = ee.Image.cat(bands)
    return out.addBands(out.reduce(ee.Reducer.sum()).rename("veg"))


def export_year(ee, region, year, out_scale_m, folder, max_pixels=65536):
    """Export growth-form fraction bands for a year to Drive (caller starts the task).

    Robust for large domains / coarse grids where interactive computePixels would hit
    memory limits — Earth Engine does the aggregation server-side.
    """
    img = growth_form_fractions(ee, region, year, out_scale_m, max_pixels).toFloat()
    return ee.batch.Export.image.toDrive(
        image=img, description=f"megan_gf_{year}", folder=folder,
        region=region, crs="EPSG:4326", scale=out_scale_m, maxPixels=int(1e13),
        fileNamePrefix=f"megan_gf_{year}",
    )


def _growth_form_native(ee, year, region, collection_id=None):
    """4 growth-form weight bands (tree/shrub/herb/crop) + veg, at native 30 m."""
    cls = classification_image(ee, year, collection_id).clip(region)
    bands = []
    for gf in landcover.GROWTH_FORMS:
        codes, weights = _remap(gf)
        bands.append(cls.remap(codes, weights, 0).rename(gf).toFloat())
    out = ee.Image.cat(bands)
    return out.addBands(out.reduce(ee.Reducer.sum()).rename("veg"))


def domain_crs_transform(domain):
    """The EE crsTransform [xScale,xShear,xTrans,yShear,yScale,yTrans] for a domain grid.

    Matches the fixed affine used by the notebooks' interactive ``sample()`` so an exported
    GeoTIFF lines up **pixel-for-pixel** with the interactively sampled met and LAI.
    """
    r = float(domain["resolution_deg"])
    return [r, 0.0, float(domain["lon"][0] - r / 2),
            0.0, -r, float(domain["lat"][-1] + r / 2)]


def export_growth_forms_to_domain(ee, region, year, domain, folder,
                                  collection_id=None, max_pixels=65536):
    """Export growth-form fractions aggregated to the EXACT analysis grid, to Drive.

    This is the **production** path for 30 m MapBiomas: the 30 m→analysis-grid aggregation is
    done server-side (no interactive ``computePixels`` hang), onto the same grid as the met/LAI
    (via :func:`domain_crs_transform`). Caller does ``task = …; task.start()`` and polls
    ``task.status()``; when done, read the GeoTIFF with :func:`read_growth_form_geotiff`.
    """
    res_m = float(domain["resolution_deg"]) * 111000.0
    img = (_growth_form_native(ee, year, region, collection_id)
           .reduceResolution(ee.Reducer.mean(), maxPixels=max_pixels)
           .reproject(crs="EPSG:4326", scale=res_m)
           .toFloat())
    return ee.batch.Export.image.toDrive(
        image=img, description=f"megan_gf_{int(year)}", folder=folder,
        crs="EPSG:4326", crsTransform=domain_crs_transform(domain),
        dimensions=f"{int(domain['lon'].size)}x{int(domain['lat'].size)}",
        maxPixels=int(1e13), fileNamePrefix=f"megan_gf_{int(year)}",
    )


def _tile_edges(n, ntiles):
    import numpy as np
    e = np.linspace(0, n, ntiles + 1).astype(int)
    return list(zip(e[:-1], e[1:]))


def sample_image_tiled(ee, img, domain, bands, ntiles=8, retries=3):
    """Sample an EE image onto the analysis grid **in tiles**, in memory (no Drive/Export).

    Aggregating 30 m MapBiomas to the ~2.8 km grid over the whole state in one interactive
    ``computePixels`` hangs; exporting to Drive then reading back is fragile (duplicate Drive
    folders, sync lag). This splits the grid into ``ntiles``×``ntiles`` blocks and pulls each
    with its own fixed affine, then stitches — so the result comes straight back as numpy with
    no Drive dependency. Each block carries far fewer native pixels, so it fits the compute
    budget. Returns ``{band: (ny,nx) array}`` on the south-first domain grid.
    """
    import numpy as np

    res = float(domain["resolution_deg"])
    lon = np.asarray(domain["lon"], dtype="float64"); lat = np.asarray(domain["lat"], dtype="float64")
    ny, nx = lat.size, lon.size
    top_x = float(lon[0] - res / 2); top_y = float(lat[-1] + res / 2)   # north-up top-left
    out = {b: np.full((ny, nx), np.nan) for b in bands}                 # north-up buffer
    sel = img.select(bands)
    for (y0, y1) in _tile_edges(ny, ntiles):
        for (x0, x1) in _tile_edges(nx, ntiles):
            w, h = int(x1 - x0), int(y1 - y0)
            if w == 0 or h == 0:
                continue
            grid = {"dimensions": {"width": w, "height": h},
                    "affineTransform": {"scaleX": res, "shearX": 0.0, "translateX": top_x + x0 * res,
                                        "shearY": 0.0, "scaleY": -res, "translateY": top_y - y0 * res},
                    "crsCode": "EPSG:4326"}
            for attempt in range(retries):
                try:
                    st = ee.data.computePixels({"expression": sel, "grid": grid,
                                                "fileFormat": "NUMPY_NDARRAY"})
                    break
                except Exception:
                    if attempt == retries - 1:
                        raise
            for b in bands:
                out[b][y0:y1, x0:x1] = np.asarray(st[b], dtype="float64")   # north-up tile -> north-up buffer
    return {b: np.flipud(out[b]) for b in bands}                           # -> south-first


def growth_form_fractions_gridded(ee, region, domain, year, ntiles=8, collection_id=None,
                                  max_pixels=65536):
    """Growth-form fractions on the analysis grid, **in memory** (tiled; no Drive/Export).

    The robust production path: 30 m MapBiomas → growth-form weights → reduceResolution to the
    grid, pulled tile-by-tile via :func:`sample_image_tiled`. Ocean / no-data cells come back
    NaN and are set to 0 (non-vegetated). Returns ``{tree,shrub,herb,crop,veg}``.
    """
    import numpy as np

    res_m = float(domain["resolution_deg"]) * 111000.0
    img = (_growth_form_native(ee, year, region, collection_id)
           .reduceResolution(ee.Reducer.mean(), maxPixels=max_pixels)
           .reproject(crs="EPSG:4326", scale=res_m))
    bands = list(landcover.GROWTH_FORMS) + ["veg"]
    g = sample_image_tiled(ee, img, domain, bands, ntiles=ntiles)
    return {b: np.nan_to_num(g[b], nan=0.0) for b in bands}


def read_growth_form_geotiff(path, domain):
    """Read an exported growth-form GeoTIFF into ``{tree,shrub,herb,crop,veg}`` on the grid.

    Bands are written in GROWTH_FORMS order (+ veg). The GeoTIFF's row 0 is north (negative
    yScale), so we flip to the domain's south-first (lat-ascending) convention — matching the
    ``np.flipud`` the interactive ``sample()`` applies.
    """
    import numpy as np
    import rasterio

    names = list(landcover.GROWTH_FORMS) + ["veg"]
    with rasterio.open(path) as src:
        arr = src.read().astype("float64")            # (bands, H, W)
    ny, nx = domain["lat"].size, domain["lon"].size
    assert arr.shape[1:] == (ny, nx), f"GeoTIFF {arr.shape[1:]} != domain grid {(ny, nx)}"
    # ocean / no-MapBiomas-coverage cells come back as NaN → non-vegetated (0), not NaN,
    # so growth-form means are real and no NaN propagates into the EF grid.
    out = {names[i]: np.nan_to_num(np.flipud(arr[i]), nan=0.0)
           for i in range(min(len(names), arr.shape[0]))}
    return out
