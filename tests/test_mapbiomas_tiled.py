"""Tests for the tiled in-memory MapBiomas sampler (no Earth Engine needed).

``sample_image_tiled`` splits the analysis grid into blocks, pulls each with its own fixed affine
via ``ee.data.computePixels``, and stitches them back — replacing the fragile Drive/Export route.
The correctness that matters: every tile must land at the right grid position (right affine), and
the north-up EE buffer must be flipped to the domain's south-first (lat-ascending) convention.

We mock ``computePixels`` to return, for each requested tile, an array that *encodes the geographic
coordinate of every pixel* (``lat*1000 + lon``) computed from the tile's own affineTransform. So if
any tile were placed at the wrong offset, requested with the wrong affine, or flipped the wrong way,
the stitched field would not equal the analytic ground truth ``lat[row]*1000 + lon[col]``.
"""
import numpy as np
import pytest

from megan_sc.gee import mapbiomas as mb


class _MockImg:
    """Stand-in for an ee.Image: remembers the encoding fn and which bands were selected."""
    def __init__(self, fn, bands=None):
        self.fn = fn
        self.bands = bands

    def select(self, bands):
        return _MockImg(self.fn, list(bands))


class _MockData:
    """Stand-in for ee.data: computePixels honours the tile's affineTransform exactly."""
    def computePixels(self, params):
        img = params["expression"]
        grid = params["grid"]
        w = grid["dimensions"]["width"]
        h = grid["dimensions"]["height"]
        aff = grid["affineTransform"]
        tx, ty, res = aff["translateX"], aff["translateY"], aff["scaleX"]
        assert grid["crsCode"] == "EPSG:4326"
        assert aff["scaleY"] == -res            # north-up
        out = {}
        for b in img.bands:
            arr = np.empty((h, w), dtype="float64")
            for k in range(h):                  # row 0 = north (top of tile)
                for j in range(w):
                    la = ty - (k + 0.5) * res    # pixel-centre latitude
                    lo = tx + (j + 0.5) * res    # pixel-centre longitude
                    arr[k, j] = img.fn(la, lo, b)
            out[b] = arr
        return out


class _MockEE:
    def __init__(self):
        self.data = _MockData()


def _domain(res, nx, ny, lon0=-53.0, lat0=-29.0):
    lon = lon0 + res * np.arange(nx)            # ascending
    lat = lat0 + res * np.arange(ny)            # ascending (south-first)
    return {"resolution_deg": res, "lon": lon, "lat": lat}


def _encode(la, lo, b):
    return la * 1000.0 + lo + (0.5 if b == "y" else 0.0)


@pytest.mark.parametrize("ntiles", [1, 2, 3, 4])
def test_tiled_stitch_matches_ground_truth(ntiles):
    res = 0.1
    dom = _domain(res, nx=10, ny=8)             # 10x8, tile edges won't divide evenly for 3/4
    lon, lat = dom["lon"], dom["lat"]
    got = mb.sample_image_tiled(_MockEE(), _MockImg(_encode), dom, ["x", "y"], ntiles=ntiles)
    ny, nx = lat.size, lon.size
    exp = np.array([[lat[r] * 1000.0 + lon[c] for c in range(nx)] for r in range(ny)])
    assert got["x"].shape == (ny, nx)
    assert np.allclose(got["x"], exp, atol=1e-9)          # right position + right affine everywhere
    assert np.allclose(got["y"], exp + 0.5, atol=1e-9)    # per-band selection carried through


def test_south_first_orientation():
    # row 0 of the returned array must be the SOUTHERNMOST latitude (domain convention).
    res = 0.1
    dom = _domain(res, nx=4, ny=5)
    got = mb.sample_image_tiled(_MockEE(), _MockImg(_encode), dom, ["x"], ntiles=2)
    lat, lon = dom["lat"], dom["lon"]
    # bottom row encodes lat[0] (south), top row encodes lat[-1] (north)
    assert got["x"][0, 0] == pytest.approx(lat[0] * 1000.0 + lon[0], abs=1e-9)
    assert got["x"][-1, 0] == pytest.approx(lat[-1] * 1000.0 + lon[0], abs=1e-9)
    assert got["x"][0, 0] < got["x"][-1, 0]               # south latitude is more negative


def test_tile_edges_partition_exactly():
    # blocks must tile [0,n) with no gap/overlap for any ntiles
    for n in (1, 7, 8, 10, 144, 224):
        for ntiles in (1, 3, 8):
            edges = mb._tile_edges(n, ntiles)
            assert edges[0][0] == 0 and edges[-1][1] == n
            for (a, b), (c, d) in zip(edges, edges[1:]):
                assert b == c                               # contiguous, no gap/overlap
            assert sum(b - a for a, b in edges) == n        # covers everything


def test_retry_then_succeed(monkeypatch):
    # a transient computePixels failure is retried, not fatal
    res = 0.1
    dom = _domain(res, nx=4, ny=4)
    ee = _MockEE()
    calls = {"n": 0}
    real = ee.data.computePixels

    def flaky(params):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("transient EE hiccup")
        return real(params)

    ee.data.computePixels = flaky
    got = mb.sample_image_tiled(ee, _MockImg(_encode), dom, ["x"], ntiles=1, retries=3)
    assert calls["n"] == 2                                   # failed once, retried, succeeded
    assert np.isfinite(got["x"]).all()
