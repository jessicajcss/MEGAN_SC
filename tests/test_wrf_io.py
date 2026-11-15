"""Tests for the WRF-ARW met reader (io/wrf.py) using a synthetic wrfout fixture.

The fixture's met fields are exact linear functions of (lon, lat), so barycentric
interpolation should reproduce them at interior target points to ~float precision —
which simultaneously verifies the variable mapping, the unit passthroughs, and that no
latitude/longitude axis was flipped during regridding.
"""
import datetime as dt
import sys, os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import make_synthetic_wrf as mkw           # noqa: E402
from megan_sc.io import wrf                 # noqa: E402


@pytest.fixture()
def wrf_file(tmp_path):
    return mkw.write(str(tmp_path / "wrfout_d01_2019-01-15_18_00_00.nc"))


def _domain(res=0.1):
    # a small regular grid well INSIDE the synthetic WRF footprint (avoid edge extrapolation)
    lon = np.arange(-52.5, -50.5 + 1e-9, res)
    lat = np.arange(-28.5, -27.0 + 1e-9, res)
    LON, LAT = np.meshgrid(lon, lat)
    return {"lon": lon, "lat": lat, "LON": LON, "LAT": LAT, "resolution_deg": res}


def test_list_times_parses_wrf_format(wrf_file):
    times = wrf.list_times(wrf_file)
    assert times[0] == dt.datetime(2019, 1, 15, 15, 0, 0)
    assert times[1] == dt.datetime(2019, 1, 15, 18, 0, 0)


def test_read_returns_expected_keys_and_shape(wrf_file):
    dom = _domain()
    met = wrf.read_wrf_met(wrf_file, dom, index=1)
    for k in ("temp", "par", "wind", "pres", "qv", "soilm", "time"):
        assert k in met
    assert met["temp"].shape == (dom["lat"].size, dom["lon"].size)
    assert met["time"] == dt.datetime(2019, 1, 15, 18)


def test_time_selection_by_when_picks_nearest(wrf_file):
    dom = _domain()
    met = wrf.read_wrf_met(wrf_file, dom, when="2019-01-15T17:10:00")
    assert met["time"] == dt.datetime(2019, 1, 15, 18)   # 17:10 is closer to 18:00 than 15:00


def test_linear_fields_interpolate_exactly(wrf_file):
    dom = _domain()
    met = wrf.read_wrf_met(wrf_file, dom, index=0)
    LON, LAT = dom["LON"], dom["LAT"]

    # TEMP2 passes through T2 unchanged, then is regridded; compare to the analytic field
    assert np.allclose(met["temp"], mkw.temp_field(LON, LAT), atol=1e-3)
    # PAR = 0.45 * SWDOWN (SWDOWN constant 600 -> 270)
    assert np.allclose(met["par"], 0.45 * mkw.swdown_field(LON, LAT), atol=1e-3)
    # WINDSPD = hypot(U10, V10)
    assert np.allclose(met["wind"],
                       np.hypot(mkw.u10_field(LON, LAT), mkw.v10_field(LON, LAT)), atol=1e-3)
    # PRES passthrough, QV passthrough (mixing ratio), SOILM top layer passthrough
    assert np.allclose(met["pres"], mkw.psfc_field(LON, LAT), atol=1e-2)
    assert np.allclose(met["qv"], mkw.q2_field(LON, LAT), atol=1e-6)
    assert np.allclose(met["soilm"], mkw.smois_field(LON, LAT), atol=1e-4)


def test_no_axis_flip(wrf_file):
    # temp decreases with latitude (coef -0.8) and increases with longitude (coef +0.5);
    # confirm the regridded field has those gradients in the right direction.
    dom = _domain()
    met = wrf.read_wrf_met(wrf_file, dom, index=0)
    assert met["temp"][0, -1] > met["temp"][0, 0]     # east (higher lon) warmer
    assert met["temp"][0, 0] > met["temp"][-1, 0]     # south (lower lat) warmer


def test_regrid_false_returns_native(wrf_file):
    dom = _domain()
    nat = wrf.read_wrf_met(wrf_file, dom, index=0, regrid=False)
    assert "lat2d" in nat and "lon2d" in nat
    assert nat["temp"].shape == nat["lat2d"].shape       # WRF native grid, not target


def test_par_fraction_override(wrf_file):
    dom = _domain()
    met = wrf.read_wrf_met(wrf_file, dom, index=0, par_fraction=0.5)
    assert np.allclose(met["par"], 0.5 * 600.0, atol=1e-3)


def test_wrf_vegetation_no_gee_path(wrf_file):
    """IVGTYP + LAI read straight from the WRF file -> growth forms + LAI on the grid."""
    import numpy as np
    # forest band (lon < -52): expect tree-dominant, LAI ~5
    dom_forest = {"lon": np.array([-53.2, -53.0, -52.8]),
                  "lat": np.array([-27.5, -27.3])}
    dom_forest["LON"], dom_forest["LAT"] = np.meshgrid(dom_forest["lon"], dom_forest["lat"])
    veg = wrf.read_wrf_vegetation(wrf_file, dom_forest, index=0)
    for k in ("lai", "tree", "shrub", "herb", "crop", "veg"):
        assert k in veg and veg[k].shape == dom_forest["LON"].shape
    assert veg["tree"].mean() > 0.8            # forest -> tree
    assert veg["crop"].mean() < 0.2
    assert abs(veg["lai"].mean() - 5.0) < 0.6  # WRF LAI ~5 over forest


def test_wrf_vegetation_crop_and_water(wrf_file):
    import numpy as np
    # crop band (-52..-50.5): crop-dominant, LAI ~1.5
    lon = np.array([-51.4, -51.2, -51.0]); lat = np.array([-27.4, -27.2])
    LON, LAT = np.meshgrid(lon, lat)
    veg = wrf.read_wrf_vegetation(wrf_file, {"lon": lon, "lat": lat, "LON": LON, "LAT": LAT})
    assert veg["crop"].mean() > 0.7
    assert abs(veg["lai"].mean() - 1.5) < 0.5
    # water band (lon > -50.5): masked -> veg ~0, LAI 0
    lon2 = np.array([-49.5, -49.2]); lat2 = np.array([-27.4, -27.2])
    LON2, LAT2 = np.meshgrid(lon2, lat2)
    vw = wrf.read_wrf_vegetation(wrf_file, {"lon": lon2, "lat": lat2, "LON": LON2, "LAT": LAT2})
    assert vw["veg"].mean() < 0.2
    assert vw["lai"].mean() < 0.2               # LAI masked over water
