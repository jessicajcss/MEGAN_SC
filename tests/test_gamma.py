"""Unit tests for the MEGANv3.21 activity-factor port.

Every assertion is grounded in the Fortran source (megvea.f / MEGVEA.EXT), not
in an external reference. These lock the port to the model's exact behaviour.
"""
import numpy as np
import pytest

from megan_sc.engine import gamma as g
from megan_sc.engine import species as sp


# ---- GAMTLI: light-independent temperature response --------------------------
def test_gamtli_is_unity_at_reference_temperature():
    # exp(beta*(303.15-303.15)) == 1 for every class
    for s in range(sp.NCLASS):
        assert g.gamtli(sp.TS_LI, s) == pytest.approx(1.0, abs=1e-12)


def test_gamtli_monotonic_increasing_in_temperature():
    temps = np.array([290.0, 300.0, 310.0, 320.0])
    vals = g.gamtli(temps, sp.index_of("ISOP"))
    assert np.all(np.diff(vals) > 0)


# ---- GAMP: light response ----------------------------------------------------
def test_gamp_zero_when_no_daily_light():
    assert g.gamp(1000.0, 0.0) == 0.0          # PPFD24 < 0.01 gate


def test_gamp_zero_at_zero_leaf_ppfd():
    assert g.gamp(0.0, 1000.0) == 0.0


def test_gamp_saturates_to_c1():
    # As PPFD1 -> inf, GAMP -> C1 = 1.03
    assert g.gamp(1e9, 1000.0) == pytest.approx(sp.C1_LIGHT, rel=1e-6)


# ---- GAMCO2: CO2 inhibition (isoprene) --------------------------------------
def test_gamma_co2_unity_at_400ppm():
    assert g.gamma_co2(400.0) == pytest.approx(1.0)


def test_gamma_co2_below_one_above_400ppm():
    assert g.gamma_co2(600.0) < 1.0


# ---- GAMLA: leaf age ---------------------------------------------------------
def test_gamma_a_unity_for_flat_age_class():
    # ACTO (class 12, idx 11) has Anew=Agro=Amat=Aold=1 -> GAMLA==1 for any LAI
    idx = sp.index_of("ACTO")
    laip = np.array([1.0, 3.0, 5.0])
    laic = np.array([2.0, 3.0, 4.0])   # growth, steady, senescence
    temp = np.array([295.0, 305.0, 300.0])
    out = g.gamma_a(laip, laic, temp, idx, nlai=46)
    assert np.allclose(out, 1.0)


def test_gamma_a_isoprene_new_leaves_suppressed():
    # Strong growth from low LAIp -> large Fnew -> GAMLA pulled toward Anew(=0.05)<1
    idx = sp.index_of("ISOP")
    out = g.gamma_a(np.array([0.1]), np.array([5.0]), np.array([300.0]), idx, 46)
    assert out[0] < 1.0


# ---- GAMMA_CD: canopy depth --------------------------------------------------
def test_gamma_cd_matches_formula_and_caps_at_three():
    cdea = g.gamma_cd(np.array([5.0]), layers=5)[0]
    # layer 3 node = 0.5 -> LAIdepth 2.5 -> -0.2*2.5+1.3 = 0.8
    assert cdea[2] == pytest.approx(0.8, abs=1e-6)
    # layer 5 node -> LAIdepth 4.76 capped to 3 -> -0.2*3+1.3 = 0.7
    assert cdea[4] == pytest.approx(0.7, abs=1e-6)


# ---- GAMTLD: light-dependent temperature ------------------------------------
def test_gamtld_zero_below_260K():
    assert g.gamtld(250.0, 300.0, sp.index_of("ISOP")) == 0.0


def test_gamtld_positive_in_normal_range():
    assert g.gamtld(303.0, 300.0, sp.index_of("ISOP")) > 0.0


# ---- GAMTP: canopy integration ----------------------------------------------
def test_canopy_weights_sum_to_one():
    _, w = sp.canopy_nodes_weights(5)
    assert w.sum() == pytest.approx(1.0, abs=1e-6)


def test_gamtp_light_independent_unity_at_reference_temp():
    # LDF=0 (pure light-independent) with every leaf temp at 303.15 K:
    # Ea2L == 1 in all layers, so GAMTP == sum(VPGWT) == 1.
    layers = 5
    shape = (3,)  # three grid cells
    T = np.full(shape + (layers,), sp.TS_LI)
    P = np.full(shape + (layers,), 500.0)
    F = np.full(shape + (layers,), 0.5)
    lai = np.full(shape, 4.0)
    out = g.gamtp(lai, T, T, P, P, F, np.full(shape, 300.0),
                  np.full(shape, 400.0), ldf=np.zeros(shape), s=sp.index_of("MT_AROM"),
                  layers=layers)
    assert np.allclose(out, 1.0, atol=1e-4)


def test_gamtp_preserves_grid_shape():
    layers, shape = 5, (4, 6)
    T = np.full(shape + (layers,), 300.0)
    P = np.full(shape + (layers,), 800.0)
    F = np.full(shape + (layers,), 0.6)
    lai = np.full(shape, 3.0)
    out = g.gamtp(lai, T, T, P, P, F, np.full(shape, 298.0),
                  np.full(shape, 450.0), ldf=np.full(shape, 1.0),
                  s=sp.index_of("ISOP"), layers=layers)
    assert out.shape == shape
    assert np.all(np.isfinite(out))
