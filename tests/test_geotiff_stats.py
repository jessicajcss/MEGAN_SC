"""Tests for GeoTIFF export, descriptive/zonal stats, and checkpointing."""
import numpy as np
import pytest

from megan_sc.io import geotiff
from megan_sc import stats


def test_geotiff_roundtrip_georef(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    lon = np.array([-53.0, -52.0, -51.0]); lat = np.array([-29.0, -28.0])
    arr = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])       # (ny=2, nx=3), south-first
    p = geotiff.write_geotiff(str(tmp_path / "t.tif"), arr, lon, lat, band_names=["x"])
    with rasterio.open(p) as src:
        assert src.crs.to_epsg() == 4326
        assert src.width == 3 and src.height == 2
        data = src.read(1)
        # north-up: row 0 is the northern lat (-28) row -> [4,5,6]
        assert np.allclose(data[0], [4, 5, 6])
        assert np.allclose(data[-1], [1, 2, 3])
        # top-left corner at (lon0-res/2, lat_max+res/2)
        assert src.transform.c == pytest.approx(-53.5)
        assert src.transform.f == pytest.approx(-27.5)


def test_field_descriptive():
    a = np.array([[1.0, 2.0], [3.0, np.nan]])
    s = stats.field_descriptive(a)
    assert s["n"] == 3 and s["sum"] == 6.0 and s["max"] == 3.0


def test_zonal_sums():
    field = np.array([[1.0, 2.0], [3.0, 4.0]])
    zone = np.array([[0, 0], [1, -1]])       # zone0: cells 1,2 ; zone1: cell 3 ; -1 excluded
    z = stats.zonal_sums(field, zone, 2)
    assert z[0] == 3.0 and z[1] == 3.0


def test_checkpoint_roundtrip(tmp_path):
    arrays = {"t": np.arange(12, dtype=float).reshape(3, 2, 2),
              "di": np.ones((3, 4))}                      # different shapes
    p = stats.save_checkpoint(str(tmp_path / "ck.npz"), arrays, {0, 1, 2})
    loaded = stats.load_checkpoint(p)
    assert loaded is not None
    a2, done = loaded
    assert np.allclose(a2["t"], arrays["t"]) and a2["di"].shape == (3, 4)
    assert done == {0, 1, 2}
    assert stats.load_checkpoint(str(tmp_path / "nope.npz")) is None
