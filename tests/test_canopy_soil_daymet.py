"""Tests for the soil-moisture, daily-met, and canopy ports.

Soil and daymet assertions are exact (simple, fully specified formulas).
Canopy assertions are physical-sanity + shape/finiteness checks; exact
cell-by-cell validation happens later against the bundled SE-US CANMET output.
"""
import numpy as np
import pytest

from megan_sc.engine import soil, daymet, canopy


# ---- soil moisture (GAMSM) ---------------------------------------------------
def test_gamma_sm_three_regimes():
    wilt = 0.10
    assert soil.gamma_sm(0.05, wilt) == 0.0                       # below wilting
    assert soil.gamma_sm(0.12, wilt) == pytest.approx(0.5)        # ramp: (0.12-0.10)/0.04
    assert soil.gamma_sm(0.30, wilt) == 1.0                       # saturated ceiling


def test_gamma_sm_at_wilting_point_is_zero_edge():
    assert soil.gamma_sm(0.10, 0.10) == pytest.approx(0.0)


def test_wilting_point_lookup():
    assert soil.wilting_point(1) == pytest.approx(0.068)
    assert soil.wilting_point(3) == pytest.approx(0.114)
    assert soil.wilting_point(8) == pytest.approx(0.250)


# ---- daily met ---------------------------------------------------------------
def test_daily_aggregate_values():
    hours = 24
    temp = np.linspace(290.0, 300.0, hours)[:, None, None] * np.ones((hours, 2, 2))
    par = np.full((hours, 2, 2), 100.0)
    wind = np.linspace(0.0, 5.0, hours)[:, None, None] * np.ones((hours, 2, 2))
    out = daymet.daily_aggregate(temp, par, wind)
    assert out["D_TEMP"] == pytest.approx(np.full((2, 2), 295.0))
    assert out["D_PPFD"] == pytest.approx(np.full((2, 2), 100.0 * 4.5))  # mean PAR * 4.5
    assert out["MaxT"] == pytest.approx(np.full((2, 2), 300.0))
    assert out["MinT"] == pytest.approx(np.full((2, 2), 290.0))
    assert out["MaxWS"] == pytest.approx(np.full((2, 2), 5.0))


# ---- canopy: solar geometry --------------------------------------------------
def test_calc_beta_overhead_sun_near_zenith():
    # Equator, ~equinox (day 80), local solar noon -> sun near zenith (beta ~ 90)
    beta = canopy.calc_beta(80.0, 0.0, 12.0)
    assert 80.0 < float(beta) <= 90.0


def test_calc_beta_in_valid_range():
    beta = canopy.calc_beta(np.array([1, 100, 200, 300]), -27.0, 12.0)
    assert np.all(beta >= -90.0) and np.all(beta <= 90.0)


# ---- canopy: radiation night branch -----------------------------------------
def test_canopy_rad_night_defaults():
    lai = np.array([4.0])
    sinbeta = np.array([-0.3])   # sun below horizon -> night branch
    z = np.array([0.0])
    sf, snp, shp, *_ = canopy.canopy_rad(lai, sinbeta, z, z, z, z, ctype=2)
    assert np.allclose(sf, 0.2)          # night sun fraction
    assert np.allclose(snp, 0.0) and np.allclose(shp, 0.0)


# ---- canopy: full driver -----------------------------------------------------
def _grid(val, shape=(2, 2)):
    return np.full(shape, val, dtype=float)


def test_megcan_daytime_physical_ranges():
    shape = (2, 2)
    ctf = np.zeros((6, *shape))
    ctf[2] = 100.0  # 100% temperate broadleaf
    out = canopy.megcan(
        temp=_grid(298.0), par_wm2=_grid(400.0), wind=_grid(2.0),
        pres=_grid(1.0e5), qv=_grid(0.010), lat=_grid(-27.0), lon=_grid(-50.0),
        lai=_grid(4.0), ctf=ctf, day=15, hour_utc=15.0,
    )
    for k in ("SunT", "ShaT", "SunP", "ShaP", "SunF"):
        assert out[k].shape == (2, 2, 5)
        assert np.all(np.isfinite(out[k]))
    # leaf temps within a physically plausible band of air temperature
    assert np.all(np.abs(out["SunT"] - 298.0) <= 15.0)
    assert np.all(np.abs(out["ShaT"] - 298.0) <= 15.0)
    # sun leaves receive at least as much light as shade leaves
    assert np.all(out["SunP"] >= out["ShaP"] - 1e-6)
    # sun fraction is a fraction
    assert np.all((out["SunF"] >= 0.0) & (out["SunF"] <= 1.0))


def test_megcan_no_vegetation_falls_back_to_air():
    shape = (2, 2)
    ctf = np.zeros((6, *shape))  # no canopy
    out = canopy.megcan(
        temp=_grid(300.0), par_wm2=_grid(500.0), wind=_grid(1.0),
        pres=_grid(1.0e5), qv=_grid(0.01), lat=_grid(-27.0), lon=_grid(-50.0),
        lai=_grid(0.0), ctf=ctf, day=15, hour_utc=15.0,
    )
    assert np.allclose(out["SunT"], 300.0)
    assert np.allclose(out["SunP"], 500.0 * 4.5)   # PAR -> PPFD
    assert np.allclose(out["SunF"], 1.0)


def test_megcan_nighttime_finite():
    shape = (2, 2)
    ctf = np.zeros((6, *shape))
    ctf[2] = 100.0
    out = canopy.megcan(
        temp=_grid(290.0), par_wm2=_grid(0.0), wind=_grid(1.5),
        pres=_grid(1.0e5), qv=_grid(0.008), lat=_grid(-27.0), lon=_grid(-50.0),
        lai=_grid(3.0), ctf=ctf, day=200, hour_utc=3.0,
    )
    for k in out:
        assert np.all(np.isfinite(out[k]))
    assert np.allclose(out["SunP"], 0.0)     # no light at night
    assert np.allclose(out["SunF"], 0.2)     # night sun fraction
