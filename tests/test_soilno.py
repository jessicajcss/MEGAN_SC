"""Tests for the soil-NO (Yienger-Levy 1995) port."""
import numpy as np
import pytest

from megan_sc.engine import soil


def test_veg_adj_unity_at_zero_lai_and_decreasing():
    assert soil.veg_adj(0.0) == pytest.approx(1.0)      # (e^0 + e^0)/2
    lais = np.array([0.0, 1.0, 3.0, 6.0])
    assert np.all(np.diff(soil.veg_adj(lais)) < 0)      # canopy reduces soil-NO with LAI


def test_growseason_sc_faithful_quirk():
    lat = np.array([-27.5])                             # Santa Catarina (SH temperate)
    # FAITHFUL source behaviour: day-of-year 15 -> non-growing; 200 -> GDAY = jday + 61
    assert soil.growseason(2020015, lat, corrected=False)[0][0] == 0
    assert soil.growseason(2020200, lat, corrected=False)[0][0] == pytest.approx(200 + 61)


def test_growseason_corrected_differs_in_january():
    lat = np.array([-27.5])
    faithful, _ = soil.growseason(2020015, lat, corrected=False)
    corrected, _ = soil.growseason(2020015, lat, corrected=True)
    assert faithful[0] == 0            # faithful: Jan is non-growing
    assert corrected[0] > 0            # corrected: Jan IS in the SH growing season (Nov-May)


def test_fertilizer_adj_seasonal_sh():
    lat = np.array([-27.5])
    # corrected (default) SH season: July (~doy 200) is winter -> 0; January is summer -> >0
    assert soil.fertilizer_adj(2020200, lat)[0] == 0.0
    assert soil.fertilizer_adj(2020015, lat)[0] > 0.0
    assert soil.fertilizer_adj(2020015, lat)[0] <= 1.0     # clamped, physical


def test_soilnox_grass_factor_warm_vs_cold():
    warm = soil.soilnox(np.array([300.0]), np.array([2.0]), np.array([-27.5]), 2020200)[1]
    cold = soil.soilnox(np.array([260.0]), np.array([2.0]), np.array([-27.5]), 2020200)[1]
    assert warm[0] > 0.0
    assert cold[0] == 0.0            # CFNOG = 0 below 268.869 K


def test_soilnox_soil_path_finite():
    cfno, cfnog = soil.soilnox(
        np.array([298.0]), np.array([3.0]), np.array([-27.5]), 2020200,
        soilt=np.array([295.0]), soilm=np.array([0.3]), soiltype=np.array([6]))
    assert np.all(np.isfinite(cfno)) and np.all(np.isfinite(cfnog))
    assert cfno[0] >= 0.0


def test_gamma_no_nongrowing_uses_grass():
    ctf = np.zeros((6, 1, 1)); ctf[5] = 100.0            # 100% crop
    cfno = np.array([[5.0]]); cfnog = np.array([[2.0]])
    # non-growing -> CFNOG everywhere regardless of crop fraction
    g0 = soil.gamma_no(cfno, cfnog, ctf, gday=np.array([[0.0]]))
    assert g0[0, 0] == pytest.approx(2.0)
    # growing, 100% crop -> CFNO
    g1 = soil.gamma_no(cfno, cfnog, ctf, gday=np.array([[100.0]]))
    assert g1[0, 0] == pytest.approx(5.0)
