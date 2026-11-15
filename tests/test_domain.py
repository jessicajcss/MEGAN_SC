"""Tests for the domain / gridding module (geodesic cell areas)."""
import numpy as np
import pytest

from megan_sc import domain


def test_make_axes_shapes():
    lon, lat = domain.make_axes(-54.0, -48.0, -30.0, -26.0, 0.5)
    assert lon.size == 12 and lat.size == 8
    # cell-centred: first centre is half a cell in from the edge
    assert lon[0] == pytest.approx(-54.0 + 0.25)


def test_cell_area_sums_to_analytic_bbox_area():
    lon_min, lon_max, lat_min, lat_max, res = -54.0, -48.0, -30.0, -26.0, 0.25
    lon, lat = domain.make_axes(lon_min, lon_max, lat_min, lat_max, res)
    area = domain.cell_area_grid(lat, res, lon.size)
    R = domain.R_AUTHALIC
    analytic = (R ** 2 * np.radians(lon_max - lon_min)
                * (np.sin(np.radians(lat_max)) - np.sin(np.radians(lat_min))))
    assert area.sum() == pytest.approx(analytic, rel=1e-9)


def test_area_positive_and_shrinks_toward_pole():
    lon, lat = domain.make_axes(-54.0, -48.0, -30.0, -26.0, 0.5)
    a = domain.cell_area_spherical(lat, 0.5)
    assert np.all(a > 0)
    # southern hemisphere: more-negative latitude is closer to the pole -> smaller area
    assert a[0] < a[-1]


def test_spherical_matches_ellipsoidal_within_tolerance():
    lat = np.array([-29.0, -27.0, -25.0])
    sph = domain.cell_area_spherical(lat, 0.1)
    ell = domain.cell_area_ellipsoidal(lat, 0.1)
    assert np.allclose(sph, ell, rtol=3e-3)   # agree within 0.3%


def test_build_domain_bbox():
    cfg = {"domain": {"mode": "bbox", "resolution_deg": 0.1,
                      "bbox": {"lon_min": -53.9, "lon_max": -48.3,
                               "lat_min": -29.5, "lat_max": -25.9}}}
    dom = domain.build_domain(cfg)
    assert dom["LON"].shape == dom["area"].shape == (dom["lat"].size, dom["lon"].size)
    assert dom["area"].min() > 0 and np.all(np.isfinite(dom["area"]))
    # SC bbox is ~56x36 cells at 0.1 deg
    assert dom["lon"].size == 56 and dom["lat"].size == 36
