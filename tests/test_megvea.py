"""Tests for the megvea ER assembly."""
import numpy as np
import pytest

from megan_sc.engine import species as sp
from megan_sc.engine import gamma as g
from megan_sc.engine import megvea


def _canopy(grid=(2, 2), layers=5, tleaf=300.0, ppfd=800.0, sunf=0.6):
    shp = grid + (layers,)
    return {
        "SunT": np.full(shp, tleaf), "ShaT": np.full(shp, tleaf - 1.0),
        "SunP": np.full(shp, ppfd), "ShaP": np.full(shp, ppfd * 0.3),
        "SunF": np.full(shp, sunf),
    }


def test_er_shape_and_nonnegative():
    grid = (3, 4)
    er = megvea.emission_activity(_canopy(grid), laic=np.full(grid, 4.0),
                                  laip=np.full(grid, 4.0), d_temp=np.full(grid, 298.0),
                                  d_ppfd=np.full(grid, 1000.0))
    assert er.shape == (19,) + grid
    assert np.all(er >= 0.0) and np.all(np.isfinite(er))


def test_er_zero_where_no_leaves():
    grid = (2, 2)
    er = megvea.emission_activity(_canopy(grid), laic=np.zeros(grid),
                                  laip=np.zeros(grid), d_temp=np.full(grid, 298.0),
                                  d_ppfd=np.full(grid, 1000.0))
    assert np.all(er == 0.0)   # LAIc == 0 => ER == 0


def test_er_matches_manual_formula_for_plain_class():
    # MT_AROM (class 7, idx 6) is not special-cased: ER = LAIc*GAMTP*GAMLA*GAMSM
    grid = (2, 2)
    can = _canopy(grid)
    laic = np.full(grid, 3.5); laip = np.full(grid, 3.0)
    dt = np.full(grid, 297.0); dp = np.full(grid, 900.0)
    gamsm = np.full(grid, 0.8)
    s = sp.index_of("MT_AROM")
    er = megvea.emission_activity(can, laic, laip, dt, dp, gamsm=gamsm)
    ldf = np.full(grid, sp.LDF[s])
    expect = (laic
              * g.gamtp(laic, can["SunT"], can["ShaT"], can["SunP"], can["ShaP"],
                        can["SunF"], dt, dp, ldf, s)
              * g.gamma_a(laip, laic, dt, s, 46)
              * gamsm)
    assert np.allclose(er[s], np.where(expect > 0, expect, 0.0))


def test_isoprene_co2_off_is_unity_multiplier():
    grid = (2, 2)
    can = _canopy(grid)
    laic = np.full(grid, 4.0); laip = np.full(grid, 4.0)
    dt = np.full(grid, 298.0); dp = np.full(grid, 1000.0)
    er_off = megvea.emission_activity(can, laic, laip, dt, dp, apply_co2=False)
    er_on = megvea.emission_activity(can, laic, laip, dt, dp, apply_co2=True)
    # GAMCO2 == 1 at 400 ppm, so toggling it makes no difference
    assert np.allclose(er_off[sp.IDX_ISOP], er_on[sp.IDX_ISOP])
