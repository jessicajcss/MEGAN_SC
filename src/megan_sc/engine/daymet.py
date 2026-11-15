"""Daily meteorology aggregation — port of daymet.f.

Source: ``reference/src/DAYMET/daymet.f`` (L286-320).
Produces the daily fields megvea needs: daily-mean temperature and PPFD, and
daily max/min temperature and max wind speed. PPFD is derived from incoming PAR
(W m-2) by the MEGAN factor 4.5 (umol m-2 s-1 per W m-2), matching megcan.
"""
from __future__ import annotations

import numpy as np

PAR_WM2_TO_PPFD = 4.5   # daymet.f L320 / megcan.f L437


def daily_aggregate(temp_k, par_wm2, wind_ms, axis=0):
    """Aggregate hourly met to MEGAN daily fields along ``axis`` (default 0).

    Parameters
    ----------
    temp_k   : hourly 2 m temperature (K), time along ``axis``
    par_wm2  : hourly incoming PAR (W m-2)
    wind_ms  : hourly wind speed (m s-1)

    Returns a dict with keys: ``D_TEMP`` (K), ``D_PPFD`` (umol m-2 s-1),
    ``MaxT`` (K), ``MinT`` (K), ``MaxWS`` (m s-1). Each has the time axis reduced.
    """
    temp_k = np.asarray(temp_k, dtype=np.float64)
    par_wm2 = np.asarray(par_wm2, dtype=np.float64)
    wind_ms = np.asarray(wind_ms, dtype=np.float64)
    return {
        "D_TEMP": temp_k.mean(axis=axis),
        "D_PPFD": par_wm2.mean(axis=axis) * PAR_WM2_TO_PPFD,
        "MaxT": temp_k.max(axis=axis),
        "MinT": temp_k.min(axis=axis),
        "MaxWS": wind_ms.max(axis=axis),
    }
