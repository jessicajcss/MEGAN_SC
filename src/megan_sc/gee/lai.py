"""MODIS leaf-area index -> MEGAN LAI, via Earth Engine.

MEGAN uses current and previous-period LAI (megvea GAMLA, megcan). We use
MODIS ``MCD15A3H`` (Collection 6.1, 500 m, 4-daily). The ``Lai`` band has scale
factor 0.1. For a monthly inventory we take the monthly-mean LAI per cell; the
engine's ``nlai`` and time-step logic treats consecutive months as current/previous.

Dataset: ``MODIS/061/MCD15A3H``.
"""
from __future__ import annotations

COLLECTION = "MODIS/061/MCD15A3H"
LAI_SCALE = 0.1


def monthly_lai(ee, region, year, month):
    """Monthly-mean LAI image (unitless m2/m2) for the region, quality-masked."""
    start = ee.Date.fromYMD(year, month, 1)
    end = start.advance(1, "month")

    def _prep(img):
        # FparLai_QC bit 0 == 0 means good quality; keep those pixels.
        qc = img.select("FparLai_QC")
        good = qc.bitwiseAnd(1).eq(0)
        return img.select("Lai").multiply(LAI_SCALE).updateMask(good).rename("LAI")

    col = (ee.ImageCollection(COLLECTION)
           .filterBounds(region).filterDate(start, end).map(_prep))
    return col.mean().rename("LAI").set("year", year, "month", month)


def export_month(ee, region, year, month, folder, resolution_deg=0.05):
    """Export monthly-mean LAI to Drive at the analysis resolution.

    Reprojection to the analysis grid is done once here (avoid repeated reproject).
    Returns the ee.batch.Task (caller starts/polls it).
    """
    img = monthly_lai(ee, region, year, month).clip(region).toFloat()
    return ee.batch.Export.image.toDrive(
        image=img, description=f"megan_lai_{year}{month:02d}", folder=folder,
        region=region, crs="EPSG:4326", scale=None, maxPixels=int(1e13),
        fileNamePrefix=f"megan_lai_{year}{month:02d}",
    )
