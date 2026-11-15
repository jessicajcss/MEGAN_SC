"""Emission-activity assembly ER — port of the megvea.f main loop.

Combines the canopy light/temperature integration (GAMTP), leaf age (GAMLA),
soil moisture (GAMSM), CO2, bidirectional-LAI and stress factors into the
per-class, non-dimensional emission activity ER used by mgn2mech/speciation.

Source: ``reference/src/MEGVEA/megvea.f`` (L646-773).
"""
from __future__ import annotations

import numpy as np

from . import species as sp
from . import gamma as g

# Classes (0-based) whose LDF is read from a gridded map, not the constant.
LDF_MAP_CLASSES = sp.LDF_FROM_MAP   # (2, 3, 4, 5)


def emission_activity(canopy, laic, laip, d_temp, d_ppfd,
                      gamsm=None, ldf_map=None, nlai=46, layers=5,
                      apply_co2=False, apply_bidir=True,
                      maxt=None, mint=None, maxws=None, aqi=None):
    """Compute the 19-class emission activity ER (megvea.f).

    Parameters
    ----------
    canopy : dict with per-layer arrays SunT, ShaT, SunP, ShaP, SunF (..., layers)
    laic, laip : current / previous LAI (per cell)
    d_temp, d_ppfd : daily-average temperature (K) and PPFD (umol m-2 s-1)
    gamsm  : soil-moisture activity (per cell); None => 1 (no stress)
    ldf_map: optional (19, ...) LDF override for the map classes; None => constants
    apply_co2  : apply GAMCO2 to isoprene (default off; ==1 at 400 ppm anyway)
    apply_bidir: apply GAMBD to ethanol/acetaldehyde (class 13) (default on)
    maxt/mint/maxws/aqi : enable the corresponding stress factor if provided

    Returns ER array of shape (19, ...), clipped at 0.
    """
    SunT, ShaT = canopy["SunT"], canopy["ShaT"]
    SunP, ShaP, SunF = canopy["SunP"], canopy["ShaP"], canopy["SunF"]
    laic = np.asarray(laic, dtype=np.float64)
    laip = np.asarray(laip, dtype=np.float64)
    d_temp = np.asarray(d_temp, dtype=np.float64)
    d_ppfd = np.asarray(d_ppfd, dtype=np.float64)
    grid = laic.shape
    one = np.ones(grid, dtype=np.float64)
    gamsm = one if gamsm is None else np.asarray(gamsm, dtype=np.float64)

    er = np.zeros((sp.NCLASS,) + grid, dtype=np.float64)
    for s in range(sp.NCLASS):
        # LDF: constant, or gridded map for the four MT classes
        if ldf_map is not None and s in LDF_MAP_CLASSES:
            ldf = np.asarray(ldf_map[s], dtype=np.float64)
        else:
            ldf = np.full(grid, sp.LDF[s], dtype=np.float64)

        gamla = g.gamma_a(laip, laic, d_temp, s, nlai)
        gamtp = g.gamtp(laic, SunT, ShaT, SunP, ShaP, SunF, d_temp, d_ppfd, ldf, s, layers)

        gamht = g.gamma_ht(maxt, s) if maxt is not None else one
        gamlt = g.gamma_lt(mint, s) if mint is not None else one
        gamhw = g.gamma_hw(maxws, s) if maxws is not None else one
        gamaq = g.gamma_aq(aqi, s) if aqi is not None else one

        base = laic * gamtp * gamla * gamhw * gamaq * gamht * gamlt * gamsm
        if s == sp.IDX_ISOP:
            gamco2 = g.gamma_co2() if apply_co2 else one
            val = base * gamco2
        elif s == sp.IDX_ETOH and apply_bidir:
            val = base * g.gamma_laibidir(laic)
        else:
            val = base
        er[s] = np.where(val > 0.0, val, 0.0)
    return er
