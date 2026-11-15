"""MODIS land cover (MCD12Q1 PFT) -> MEGAN growth-form fractions, via Earth Engine.

An alternative to MapBiomas for the vegetation structure. MODIS ``MCD12Q1``'s
``LC_Type5`` band is a Plant-Functional-Type classification (500 m, annual) — a
direct fit for MEGAN's growth forms. At 500 m it aggregates to the analysis grid
with only ~30 native pixels per cell, so interactive sampling is cheap (unlike the
30 m MapBiomas aggregation). Global and consistent; coarser than MapBiomas but at a
~2.8 km analysis grid the difference is small.

Asset: ``MODIS/061/MCD12Q1`` (band ``LC_Type5``). Crosswalk in ``landcover.MODIS_PFT_CROSSWALK``.
"""
from __future__ import annotations

from .. import landcover

ASSET = "MODIS/061/MCD12Q1"
BAND = "LC_Type5"


def _remap(gf):
    codes = sorted(landcover.MODIS_PFT_CROSSWALK)
    i = landcover._GF_IDX[gf]
    return codes, [float(landcover.MODIS_PFT_CROSSWALK[c][i]) for c in codes]


def growth_form_fractions(ee, region, year, out_scale_m, max_pixels=1024):
    """ee.Image with tree/shrub/herb/crop fraction bands on the analysis grid."""
    img = (ee.ImageCollection(ASSET)
           .filterDate(f"{int(year)}-01-01", f"{int(year) + 1}-01-01")
           .first().select(BAND).clip(region))
    bands = []
    for gf in landcover.GROWTH_FORMS:
        codes, weights = _remap(gf)
        w = img.remap(codes, weights, 0).rename(gf).toFloat()
        wr = (w.reduceResolution(ee.Reducer.mean(), maxPixels=max_pixels)
              .reproject(crs="EPSG:4326", scale=out_scale_m))
        bands.append(wr)
    out = ee.Image.cat(bands)
    return out.addBands(out.reduce(ee.Reducer.sum()).rename("veg"))


def export_year(ee, region, year, out_scale_m, folder):
    """Export growth-form fraction bands for a year to Drive (caller starts the task)."""
    img = growth_form_fractions(ee, region, year, out_scale_m).toFloat()
    return ee.batch.Export.image.toDrive(
        image=img, description=f"megan_gf_modis_{year}", folder=folder,
        region=region, crs="EPSG:4326", scale=out_scale_m, maxPixels=int(1e13),
        fileNamePrefix=f"megan_gf_modis_{year}",
    )
