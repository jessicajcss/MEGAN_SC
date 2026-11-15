"""Tests for the optional soil-texture -> wilting-point path (USDA triangle -> MEGAN WWLT)."""
import numpy as np
from megan_sc.engine import soil


def test_usda_texture_canonical_points():
    # (sand%, clay%) -> expected USDA/Noah class index
    cases = {
        (92, 3): 1,    # sand
        (80, 8): 2,    # loamy sand
        (65, 10): 3,   # sandy loam
        (20, 15): 4,   # silt loam
        (6, 6): 5,     # silt
        (40, 18): 6,   # loam
        (60, 25): 7,   # sandy clay loam
        (10, 33): 8,   # silty clay loam
        (32, 33): 9,   # clay loam
        (52, 42): 10,  # sandy clay
        (6, 48): 11,   # silty clay
        (20, 60): 12,  # clay
    }
    for (s, c), expect in cases.items():
        got = int(soil.usda_texture_class(s, c))
        assert got == expect, f"sand={s} clay={c}: got {got}, expected {expect}"


def test_wilting_point_from_texture_monotonic():
    # sandy soils have a lower wilting point than clayey soils
    wp_sand = float(soil.wilting_point_from_texture(92, 3))
    wp_loam = float(soil.wilting_point_from_texture(40, 18))
    wp_clay = float(soil.wilting_point_from_texture(20, 60))
    assert wp_sand < wp_loam < wp_clay
    assert wp_sand == soil.WWLT[0]      # class 1
    assert wp_clay == soil.WWLT[11]     # class 12


def test_texture_class_vectorised():
    sand = np.array([[92.0, 20.0], [40.0, 6.0]])
    clay = np.array([[3.0, 60.0], [18.0, 48.0]])
    cls = soil.usda_texture_class(sand, clay)
    assert cls.shape == (2, 2)
    assert cls[0, 0] == 1 and cls[0, 1] == 12 and cls[1, 0] == 6 and cls[1, 1] == 11
    wp = soil.wilting_point_from_texture(sand, clay)
    assert wp.shape == (2, 2) and np.all(wp > 0)
