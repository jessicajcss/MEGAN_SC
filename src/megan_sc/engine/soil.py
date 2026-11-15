"""Soil-moisture activity factor GAMSM — port of megsea.f (soil-moisture block).

Source: ``reference/src/MEGSEA/megsea.f`` (L602-615) and ``MEGSEA.EXT``.
GAMSM scales isoprene emission down under soil-moisture stress. Soil NO (GAMNO,
the Yienger-Levy 1995 algorithm) is a separate, larger routine handled in a
dedicated module; this file covers the soil-moisture factor used by megvea.
"""
from __future__ import annotations

import numpy as np

# Wilting point (m^3/m^3) by USDA/MCIP soil type index 1..16 (MEGSEA.EXT, JN90).
WWLT = np.array(
    [0.068, 0.075, 0.114, 0.179, 0.155, 0.175, 0.218, 0.250,
     0.219, 0.283, 0.286, 0.286, 0.286, 0.286, 0.286, 0.286], dtype=np.float64
)
D1 = 0.04   # soil-moisture ramp width above wilting point (MEGSEA.EXT)


def wilting_point(soil_type) -> np.ndarray:
    """Map integer soil type (1..16) to wilting point via the WWLT table."""
    idx = np.asarray(soil_type, dtype=int) - 1
    idx = np.clip(idx, 0, WWLT.size - 1)
    return WWLT[idx]


def usda_texture_class(sand_pct, clay_pct):
    """USDA soil-texture class (1..12) from sand & clay PERCENTAGES (silt = 100−sand−clay).

    Uses the **definitional** USDA texture-triangle boundaries (not a fitted pedotransfer),
    so there are no uncertain regression coefficients. Class order matches the USDA/Noah
    16-type indexing used by :data:`WWLT` (1 sand … 12 clay), so the result feeds
    :func:`wilting_point` directly. Vectorised.

    Classes: 1 sand, 2 loamy sand, 3 sandy loam, 4 silt loam, 5 silt, 6 loam,
    7 sandy clay loam, 8 silty clay loam, 9 clay loam, 10 sandy clay, 11 silty clay, 12 clay.
    """
    S = np.asarray(sand_pct, dtype=np.float64)
    C = np.asarray(clay_pct, dtype=np.float64)
    Si = 100.0 - S - C
    conds = [
        (C >= 40) & (Si >= 40),                                   # silty clay (11)
        (C >= 40) & (S > 45),                                     # -> sandy clay region handled next
        (C >= 40),                                                # clay (12)
        (C >= 35) & (S > 45),                                     # sandy clay (10)
        (C >= 27) & (S <= 20),                                    # silty clay loam (8)
        (C >= 27) & (S <= 45),                                    # clay loam (9)
        (C >= 20) & (C < 27) & (S > 45) & (Si < 28),              # sandy clay loam (7)
        (C >= 7) & (C < 27) & (Si >= 28) & (Si < 50) & (S <= 52),  # loam (6)
        (Si >= 80) & (C < 12),                                    # silt (5)
        (Si >= 50) & (C >= 12) & (C < 27),                        # silt loam (4)
        (Si >= 50) & (C < 12),                                    # silt loam (4)  (50<=Si<80, C<12)
        (Si + 1.5 * C < 15),                                      # sand (1)
        (Si + 2.0 * C < 30),                                      # loamy sand (2)
    ]
    choices = [11, 10, 12, 10, 8, 9, 7, 6, 5, 4, 4, 1, 2]
    # default: sandy clay loam if clayey-sandy, else sandy loam (3)
    default = np.where((C >= 20) & (C < 35) & (S > 45), 7, 3)
    return np.select(conds, choices, default=default).astype(int)


def wilting_point_from_texture(sand_pct, clay_pct):
    """Wilting point (m^3/m^3) from sand & clay % via USDA class → MEGAN ``WWLT`` table.

    The MEGAN-faithful, no-fitted-coefficient route for a texture-aware wilting point (used
    by the optional MapBiomas-Solo soil path; off by default).
    """
    return wilting_point(usda_texture_class(sand_pct, clay_pct))


def gamma_sm(soil_moisture, wilt):
    """Soil-moisture activity factor (megsea.f).

        GAMSM = 0                        if theta <  wilt
              = (theta - wilt) / D1      if wilt <= theta < wilt + D1
              = 1                        if theta >= wilt + D1

    ``soil_moisture`` (theta) and ``wilt`` are volumetric water content
    (m^3/m^3); ``wilt`` may be a scalar, a per-cell wilting-point map, or the
    result of :func:`wilting_point`.
    """
    theta = np.asarray(soil_moisture, dtype=np.float64)
    wilt = np.asarray(wilt, dtype=np.float64)
    t1 = wilt + D1
    ramp = (theta - wilt) / D1
    return np.where(theta < wilt, 0.0, np.where(theta < t1, ramp, 1.0))


# ===================================================================== #
# Soil NO (Yienger & Levy 1995) — port of soilnox.F / soilnox_fx.F / megsea.f
# ===================================================================== #
import datetime as _dt

# Saturation soil water by soil type 1..12 (soilnox.F, Jacquemin & Noilhan 1990).
SATURATION = np.array(
    [0.395, 0.410, 0.435, 0.485, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.482, 0.482],
    dtype=np.float64,
)


def _g2j(year, mmdd):
    """MMDD (e.g. 1101) -> day-of-year for that year (soilnox_fx.F G2J)."""
    mm, dd = mmdd // 100, mmdd % 100
    return _dt.date(int(year), mm, dd).timetuple().tm_yday


def growseason(date_yyyyddd, lat, corrected=True):
    """Growing-season day (GDAY) and length (GLEN) per cell (soilnox_fx.F GROWSEASON).

    ``date_yyyyddd`` is an I/O-API Julian date (YYYYDDD). Returns (gday, glen) arrays.

    SOURCE QUIRK & DEFAULT CHOICE: the Fortran compares the day-of-year (DDD, 1-366)
    against month-day literals (e.g. 1101), so in the Southern Hemisphere the Nov-Dec
    branch never triggers and the "growing season" degenerates to day-of-year >= 101,
    giving GDAY = DDD + 61 (which exceeds GLEN and drives the fertilizer factor
    negative for much of the year). That is not physically defensible for an SC
    inventory, so we DEFAULT to ``corrected=True`` — the intended SH season (Nov-May,
    summer growing / winter dormant). Pass ``corrected=False`` to reproduce the exact
    (broken) MEGANv3.21 behaviour for bit-comparison. See docs/ALGORITHM_SPEC.md.
    """
    year = int(date_yyyyddd // 1000)
    jday = int(date_yyyyddd - year * 1000)
    lat = np.asarray(lat, dtype=np.float64)
    gday = np.zeros_like(lat)
    glen = np.zeros_like(lat)

    j0101 = _g2j(year, 101)
    j1231 = _g2j(year, 1231)
    may31 = _g2j(year, 531)
    nov1 = _g2j(year, 1101)
    glen_sh = 30 + 31 + may31 - j0101 + 1  # 212 (213 leap)

    # tropical (|lat| <= 23): year-round
    trop = (lat <= 23.0) & (lat >= -23.0)
    gday = np.where(trop, jday - j0101 + 1, gday)
    glen = np.where(trop, j1231 - j0101 + 1, glen)

    # southern hemisphere temperate (-60 <= lat < -23)
    sh = (lat < -23.0) & (lat >= -60.0)
    if corrected:
        gday_sh = np.where(jday >= nov1, jday - nov1 + 1,
                           np.where(jday <= may31, jday - j0101 + 1 + 61, 0.0))
    else:  # faithful: 1101 branch never fires; jday in [101,531] -> jday+61
        gday_sh = np.where((jday >= 101) & (jday <= 531), jday - j0101 + 1 + 61, 0.0)
    gday = np.where(sh, gday_sh, gday)
    glen = np.where(sh, glen_sh, glen)

    # northern hemisphere temperate (23 < lat <= 65)
    nh = (lat > 23.0) & (lat <= 65.0)
    start = np.floor((lat - 23.0) * 4.5)
    end = j1231 - np.floor((lat - 23.0) * 3.3)
    gday_nh = np.where((jday >= start) & (jday <= end), jday - start + 1, 0.0)
    gday = np.where(nh, gday_nh, gday)
    glen = np.where(nh, end - start + 1, glen)
    # polar (|lat| > 65 N / < -60 S) stay 0
    return gday, glen


def fertilizer_adj(date_yyyyddd, lat, corrected=True):
    """Fertilizer adjustment factor (soilnox_fx.F FERTLZ_ADJ)."""
    gday, glen = growseason(date_yyyyddd, lat, corrected)
    with np.errstate(invalid="ignore", divide="ignore"):
        late = 1.0 + 30.0 / glen - gday / glen
    fert = np.where(gday == 0, 0.0, np.where(gday < 30, 1.0, late))
    return np.maximum(fert, 0.0)


def veg_adj(lai):
    """Vegetation (canopy reduction) adjustment factor (soilnox_fx.F VEG_ADJ)."""
    lai = np.asarray(lai, dtype=np.float64)
    return (np.exp(-0.24 * lai) + np.exp(-0.0525 * lai)) * 0.5


def soilnox(tair, lai, lat, date_yyyyddd, *, soilt=None, soilm=None, soiltype=None,
            precadj=1.0, corrected=True):
    """Soil-NO crop (CFNO) and grass (CFNOG) factors (soilnox.F SOILNOX).

    Provide ``soilt``/``soilm``/``soiltype`` for the PX (soil) path; omit them to use
    the non-soil path (CFNODRY/CFNOWET from air temperature only). Returns (cfno, cfnog).
    """
    tair = np.minimum(np.asarray(tair, dtype=np.float64), 315.0)
    # CFNOG (grass): cap at 303 K
    tair_g = np.minimum(tair, 303.0)
    cfnog = np.where(tair_g > 268.8690, np.exp(0.04686 * tair_g - 14.30579), 0.0)

    use_soil = soilt is not None and soilm is not None and soiltype is not None
    tsoi = np.asarray(soilt, np.float64) if use_soil else 0.72 * tair + 82.28
    tsoi = np.clip(tsoi, 273.16, 303.16)
    cfnodry = (1.0 / 3.0) * (1.0 / 30.0) * (tsoi - 273.16)
    cfnowet = np.where(tsoi <= 283.16,
                       (tsoi - 273.16) * np.exp(-0.103 * 30.0) * 0.28,
                       np.exp(0.103 * (tsoi - 273.16)) * np.exp(-0.103 * 30.0))
    if use_soil:
        st = np.asarray(soiltype, dtype=int)
        sat = SATURATION[np.clip(st - 1, 0, SATURATION.size - 1)]
        ratio = np.asarray(soilm, np.float64) / sat
        cf = ratio * cfnowet + (1.0 - ratio) * cfnodry
        cf = np.where((st > 0) & (st <= SATURATION.size), cf, 0.0)
    else:
        cf = 0.5 * cfnowet + 0.5 * cfnodry

    cfno = cf * fertilizer_adj(date_yyyyddd, lat, corrected) * veg_adj(lai) * precadj
    return cfno, cfnog


def gamma_no(cfno, cfnog, ctf, gday):
    """Soil-NO activity factor GAMNO (megsea.f assembly).

    Non-growing (gday==0): GAMNO = CFNOG everywhere. Growing: canopy-type-weighted,
    with crops (MEGAN canopy type 6) using CFNO and all others using CFNOG.
    ``ctf`` is the 6-type canopy fraction array (any consistent units; ratios cancel).
    """
    ctf = np.asarray(ctf, dtype=np.float64)
    tmo1 = ctf[:5].sum(axis=0) + ctf[5]
    tmo2 = ctf[:5].sum(axis=0) * cfnog + ctf[5] * cfno
    with np.errstate(invalid="ignore", divide="ignore"):
        grow = np.where(tmo1 > 0.0, tmo2 / tmo1, 0.0)
    return np.where(np.asarray(gday) > 0.0, grow, cfnog)
