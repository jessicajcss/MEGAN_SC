"""Tests for the MapBiomas -> MEGAN vegetation crosswalk."""
import numpy as np
import pytest

from megan_sc import landcover as lc


def test_growth_form_weights_key_classes():
    # forest -> tree; grassland/pasture -> herb; soybean -> crop; urban/water -> none
    assert tuple(lc.growth_form_weights(3)) == (1, 0, 0, 0)     # Forest Formation
    assert tuple(lc.growth_form_weights(12)) == (0, 0, 1, 0)    # Grassland
    assert tuple(lc.growth_form_weights(15)) == (0, 0, 1, 0)    # Pasture
    assert tuple(lc.growth_form_weights(39)) == (0, 0, 0, 1)    # Soybean
    assert tuple(lc.growth_form_weights(24)) == (0, 0, 0, 0)    # Urban -> masked
    assert tuple(lc.growth_form_weights(33)) == (0, 0, 0, 0)    # Ocean -> masked


def test_savanna_and_mosaic_are_mixed():
    assert tuple(lc.growth_form_weights(4)) == (0.5, 0.5, 0, 0)   # Savanna tree/shrub
    assert tuple(lc.growth_form_weights(21)) == (0, 0, 0.5, 0.5)  # Mosaic herb/crop


def test_unknown_code_is_nonveg():
    assert tuple(lc.growth_form_weights(999)) == (0, 0, 0, 0)


def test_fractions_on_grid_half_forest_half_water():
    # fine grid: left half forest (3), right half ocean (33); aggregate to 1 cell
    fine = np.full((4, 4), 3, dtype=int)
    fine[:, 2:] = 33
    flon = np.array([-53.0, -52.99, -52.98, -52.97])
    flat = np.array([-28.0, -27.99, -27.98, -27.97])
    tlon = np.array([-52.985]); tlat = np.array([-27.985])
    fr = lc.fractions_on_grid(fine, flon, flat, tlon, tlat, 0.04)
    assert fr["tree"][0, 0] == pytest.approx(0.5)
    assert fr["herb"][0, 0] == pytest.approx(0.0)
    assert fr["veg"][0, 0] == pytest.approx(0.5)   # half the cell is water -> unvegetated


def test_modis_pft_crosswalk():
    cw = lc.MODIS_PFT_CROSSWALK
    assert tuple(lc.growth_form_weights(2, cw)) == (1, 0, 0, 0)   # evergreen broadleaf -> tree
    assert tuple(lc.growth_form_weights(1, cw)) == (1, 0, 0, 0)   # evergreen needleleaf -> tree
    assert tuple(lc.growth_form_weights(5, cw)) == (0, 1, 0, 0)   # shrub
    assert tuple(lc.growth_form_weights(6, cw)) == (0, 0, 1, 0)   # grass -> herb
    assert tuple(lc.growth_form_weights(7, cw)) == (0, 0, 0, 1)   # cereal crop
    assert tuple(lc.growth_form_weights(9, cw)) == (0, 0, 0, 0)   # urban -> masked
    assert tuple(lc.growth_form_weights(0, cw)) == (0, 0, 0, 0)   # water -> masked


def test_canopy_type_fractions_map_growth_forms():
    gf = {"tree": np.array([[0.6]]), "shrub": np.array([[0.1]]),
          "herb": np.array([[0.2]]), "crop": np.array([[0.1]])}
    ctf = lc.canopy_type_fractions(gf, tree_canopy_type=3)
    assert ctf.shape == (6, 1, 1)
    assert ctf[2, 0, 0] == pytest.approx(60.0)   # tree -> type 3 (temperate broadleaf)
    assert ctf[3, 0, 0] == pytest.approx(10.0)   # shrub -> type 4
    assert ctf[4, 0, 0] == pytest.approx(20.0)   # herb -> type 5
    assert ctf[5, 0, 0] == pytest.approx(10.0)   # crop -> type 6
    assert ctf.sum() == pytest.approx(100.0)     # fully vegetated cell
