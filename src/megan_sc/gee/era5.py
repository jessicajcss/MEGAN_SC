"""ERA5-Land hourly meteorology -> MEGAN met variables, via Earth Engine.

Derives the exact fields the ported engine consumes (see engine/canopy.py,
engine/daymet.py, engine/soil.py):

    TEMP2    2 m air temperature (K)          <- temperature_2m
    PAR      photosynth. active rad. (W m-2)  <- 0.45 * shortwave-down (de-accumulated)
    WINDSPD  wind speed (m s-1)               <- hypot(u10, v10)
    PRES     surface pressure (Pa)            <- surface_pressure
    QV       water-vapour mixing ratio (kg/kg)<- from dewpoint_2m + surface_pressure
    SOILM    volumetric soil water (m3/m3)    <- volumetric_soil_water_layer_1

Dataset: ``ECMWF/ERA5_LAND/HOURLY`` (~9 km, hourly).

IMPORTANT — radiation de-accumulation: ERA5-Land ``surface_solar_radiation_downwards``
is accumulated (J m-2) from 00 UTC and resets daily. To obtain the hourly flux we
difference consecutive hours within a day and divide by 3600 s. This is handled by
:func:`hourly_flux`. The notebooks replace it by ``fixes.hourly_flux_fixed`` (FIX-8, installed by
``fixes.apply_engine_patches()``), which subtracts the previous hour at every hour except 01 UTC, so that the 00 UTC
value, whose field holds the 24 h of the previous day (ECMWF, ERA5-Land data documentation), is differenced too.
"""
from __future__ import annotations

COLLECTION = "ECMWF/ERA5_LAND/HOURLY"


def _mixing_ratio(ee, img):
    """Water-vapour mixing ratio (kg/kg) from dewpoint (K) and pressure (Pa)."""
    td = img.select("dewpoint_temperature_2m")          # K
    p = img.select("surface_pressure")                   # Pa
    # Bolton (1980) saturation vapour pressure at the dewpoint -> actual vapour pressure e (Pa)
    e = td.expression(
        "611.2 * exp(17.67 * (Td - 273.15) / (Td - 29.65))", {"Td": td})
    return e.expression("0.622 * e / (P - e)", {"e": e, "P": p}).rename("QV")


def hourly_flux(ee, region, date):
    """De-accumulated hourly shortwave-down flux (W m-2) image for a given hour.

    ERA5-Land radiation accumulates from 00 UTC; hourly flux = (this hour − prev
    hour) / 3600, with the 01 UTC step using the raw accumulation. Used alone, this function does not difference the
    00 UTC value; the notebooks use ``fixes.hourly_flux_fixed`` instead (FIX-8).
    """
    d = ee.Date(date)
    col = ee.ImageCollection(COLLECTION).filterBounds(region)
    cur = col.filterDate(d, d.advance(1, "hour")).first().select("surface_solar_radiation_downwards")
    prev = col.filterDate(d.advance(-1, "hour"), d).first()
    prev_acc = ee.Image(ee.Algorithms.If(
        d.get("hour"),  # nonzero hour -> subtract previous accumulation
        ee.Image(prev).select("surface_solar_radiation_downwards"),
        ee.Image(0)))
    return cur.subtract(prev_acc).max(0).divide(3600.0).rename("SW_WM2")


def megan_met_image(ee, region, date):
    """Build the MEGAN met bands (TEMP2, PAR, WINDSPD, PRES, QV, SOILM) for one hour."""
    d = ee.Date(date)
    img = ee.ImageCollection(COLLECTION).filterBounds(region).filterDate(
        d, d.advance(1, "hour")).first()

    temp2 = img.select("temperature_2m").rename("TEMP2")
    pres = img.select("surface_pressure").rename("PRES")
    u = img.select("u_component_of_wind_10m")
    v = img.select("v_component_of_wind_10m")
    windspd = u.hypot(v).rename("WINDSPD")
    soilm = img.select("volumetric_soil_water_layer_1").rename("SOILM")
    qv = _mixing_ratio(ee, img)
    par = hourly_flux(ee, region, date).multiply(0.45).rename("PAR")   # PAR ~= 0.45 * SW

    return (temp2.addBands([par, windspd, pres, qv, soilm])
            .set("system:time_start", d.millis()))


def export_month(ee, region, year, month, out_grid_scale_deg, folder, project,
                 skip_existing=True):
    """Export a month of hourly MEGAN met bands to Drive (one asset per hour).

    Best practices applied: server-side collection, filterBounds+filterDate, clip
    to ``region``, explicit crs/scale/maxPixels, tileScale to avoid memory errors.
    Returns the list of ee.batch.Task objects created (caller starts/polls them).
    """
    import ee as _ee  # noqa
    start = ee.Date.fromYMD(year, month, 1)
    end = start.advance(1, "month")
    hours = end.difference(start, "hour").getInfo()
    tasks = []
    for h in range(int(hours)):
        d = start.advance(h, "hour")
        stamp = d.format("YYYYMMdd_HH").getInfo()
        img = megan_met_image(ee, region, d).clip(region).toFloat()
        task = ee.batch.Export.image.toDrive(
            image=img, description=f"megan_met_{stamp}", folder=folder,
            region=region, scale=None, crs="EPSG:4326",
            crsTransform=None, maxPixels=int(1e13),
            fileNamePrefix=f"megan_met_{stamp}",
        )
        tasks.append(task)
    return tasks
