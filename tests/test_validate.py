"""Tests for the reference-comparison helpers (validate.py) — unit conversions & stats."""
import numpy as np
import pytest

from megan_sc import validate as val


def test_reference_tonnes_moles_path():
    # 1 mol/s of isoprene for 1 hour = 68.12 g/s * 3600 s = 245232 g = 0.245232 t
    t = val.reference_tonnes(np.array([[1.0]]), "moles/s", val.MW["ISOP"], 3600.0)
    assert t[0, 0] == pytest.approx(68.12 * 3600 * 1e-6)


def test_reference_tonnes_grams_path():
    t = val.reference_tonnes(np.array([[2.0]]), "g/s", 68.12, 3600.0)
    assert t[0, 0] == pytest.approx(2.0 * 3600 * 1e-6)


def test_reference_tonnes_rejects_unknown_units():
    with pytest.raises(ValueError):
        val.reference_tonnes(np.array([1.0]), "kg/hr", 68.12, 3600.0)


def test_field_stats_perfect_match():
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    s = val.field_stats(a, a.copy())
    assert s["bias"] == pytest.approx(0.0)
    assert s["rmse"] == pytest.approx(0.0)
    assert s["corr"] == pytest.approx(1.0)
    assert s["ratio"] == pytest.approx(1.0)


def test_field_stats_ratio_and_bias():
    a = np.array([2.0, 4.0, 6.0])
    b = np.array([1.0, 2.0, 3.0])
    s = val.field_stats(a, b)
    assert s["ratio"] == pytest.approx(2.0)
    assert s["bias"] == pytest.approx(2.0)         # mean(a-b) = mean(1,2,3)=2
    assert s["corr"] == pytest.approx(1.0)         # perfectly correlated


def test_our_flux_to_tonnes_conversion():
    # 1 nanomol/m2/s over area 1e6 m2 for 3600 s, MW 68.12:
    #   1e-9 mol/m2/s * 1e6 m2 = 1e-3 mol/s; *68.12 g/mol *3600 s = 245.2 g = 2.452e-4 t
    our_map = np.array([[1.0]])
    area = np.array([[1.0e6]])
    ref_rate = np.array([[0.0]])
    ref_lon = np.array([[-51.0]]); ref_lat = np.array([[-27.0]])
    out = val.compare_species(our_map, np.array([-51.0]), np.array([-27.0]), area,
                              ref_rate, ref_lon, ref_lat, "moles/s", 1.6e7,
                              val.MW["ISOP"], 3600.0)
    assert out["our_total_t"] == pytest.approx(1e-9 * 1e6 * 68.12 * 3600 * 1e-6)


def test_regrid_constant_field_is_preserved():
    ref = np.full((5, 5), 7.0)
    rlon, rlat = np.meshgrid(np.linspace(-53, -49, 5), np.linspace(-29, -26, 5))
    out = val.regrid_to_regular(ref, rlon, rlat, np.linspace(-52, -50, 4),
                                np.linspace(-28, -27, 3))
    assert np.allclose(out, 7.0)


def test_polygon_mask_square():
    from shapely.geometry import box
    from megan_sc import validate as val
    import numpy as np
    lon, lat = np.meshgrid(np.linspace(-54, -48, 7), np.linspace(-30, -25, 6))
    geom = box(-52, -29, -50, -27)                 # a 2x2 deg square
    m = val.polygon_mask(lon, lat, geom)
    assert m.shape == lon.shape
    # every masked cell must be inside the square; at least some cells qualify
    assert m.sum() > 0
    assert np.all((lon[m] >= -52) & (lon[m] <= -50) & (lat[m] >= -29) & (lat[m] <= -27))


def test_masked_total_matches_manual():
    from megan_sc import validate as val
    import numpy as np
    f = np.array([[1.0, 2.0], [3.0, 4.0]])
    m = np.array([[True, False], [True, True]])
    assert val.masked_total(f, m) == 8.0            # 1+3+4
    assert val.masked_total(f) == 10.0


def test_bbox_mask_selects_interior():
    from megan_sc import validate as val
    import numpy as np
    lon, lat = np.meshgrid(np.linspace(-56, -46, 11), np.linspace(-31, -21, 11))
    m = val.bbox_mask(lon, lat, -54, -48, -30, -26)
    assert m.shape == lon.shape
    assert np.all((lon[m] >= -54) & (lon[m] <= -48) & (lat[m] >= -30) & (lat[m] <= -26))
    # a point clearly outside is excluded
    assert not m[0, 0]        # lon -56, lat -31 -> outside
