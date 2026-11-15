"""MEGANv3.21 activity-factor functions — pure-Python port of megvea.f.

Each function reproduces the corresponding Fortran routine exactly (same
constants, same branching). All maths is done in float64 and vectorised with
numpy so a function accepts a scalar or an ndarray of any shape.

Reference: ``reference/src/MEGVEA/megvea.f`` and ``docs/ALGORITHM_SPEC.md``.
Class index ``s`` is 0-based (Fortran ``S`` == ``s + 1``).
"""
from __future__ import annotations

import numpy as np

from . import species as sp

_F = np.float64


def _arr(x) -> np.ndarray:
    return np.asarray(x, dtype=_F)


# ---------------------------------------------------------------------------
# Temperature responses
# ---------------------------------------------------------------------------
def gamtld(t_leaf, t24, s: int):
    """Light-dependent temperature response (megvea.f FUNCTION GAMTLD).

    T240 is set equal to T24 in this version (no separate 240-h history).
    Returns 0 where leaf temperature < 260 K.
    """
    t1 = _arr(t_leaf)
    t24 = _arr(t24)
    t240 = t24
    topt = 312.5 + 0.6 * (t240 - 297.0)
    x = ((1.0 / topt) - (1.0 / t1)) / 0.00831
    eopt = sp.CLEO[s] * np.exp(0.05 * (t24 - 297.0)) * np.exp(0.05 * (t240 - 297.0))
    ct1 = sp.CT1[s]
    ct2 = sp.CT2
    g = eopt * ct2 * np.exp(ct1 * x) / (ct2 - ct1 * (1.0 - np.exp(ct2 * x)))
    return np.where(t1 < 260.0, 0.0, g)


def gamtli(temp, s: int):
    """Light-independent temperature response (megvea.f FUNCTION GAMTLI).

    exp(beta_s * (T - 303.15)).  Equals 1.0 exactly at T = 303.15 K.
    """
    return np.exp(sp.BETA[s] * (_arr(temp) - sp.TS_LI))


# ---------------------------------------------------------------------------
# Light response
# ---------------------------------------------------------------------------
def gamp(ppfd_leaf, ppfd24):
    """Light response (megvea.f FUNCTION GAMP).

    Returns 0 where the daily-average PPFD < 0.01 (dark/no-data gate).
    C1 is the fixed constant 1.03 in this version (the PPFD24-dependent form is
    commented out in the source). Saturates to C1 (=1.03) as PPFD -> inf.
    """
    p1 = _arr(ppfd_leaf)
    p24 = _arr(ppfd24)
    a = sp.ALPHA_LIGHT
    c1 = sp.C1_LIGHT
    g = (a * c1 * p1) / np.sqrt(1.0 + (a ** 2) * (p1 ** 2))
    return np.where(p24 < 0.01, 0.0, g)


# ---------------------------------------------------------------------------
# Canopy-depth response
# ---------------------------------------------------------------------------
def gamma_cd(lai, layers: int = 5):
    """Canopy-depth emission response per layer (megvea.f SUBROUTINE GAMMA_CD).

    Returns array with a trailing axis of length ``layers``:
        CDEA_k = CCD1 * min(LAI * cdepth_k, 3) + CCD2
    """
    lai = _arr(lai)
    cdepth, _ = sp.canopy_nodes_weights(layers)
    laidepth = np.minimum(lai[..., None] * cdepth, 3.0)
    return sp.CCD1 * laidepth + sp.CCD2


# ---------------------------------------------------------------------------
# Canopy integration of light + temperature (GAMTP)
# ---------------------------------------------------------------------------
def gamtp(lai, sun_t, sha_t, sun_p, sha_p, sun_f, d_temp, d_ppfd, ldf, s: int,
          layers: int = 5):
    """Canopy-average light+temperature activity factor (megvea.f, GAMTP block).

    Computes CDEA internally from ``lai`` (per cell) via :func:`gamma_cd`, then
    integrates over the canopy exactly as megvea.f does::

        Ea1L_k = CDEA_k * [GAMTLD(SunT)*GAMP(SunP)*SunF
                           + GAMTLD(ShaT)*GAMP(ShaP)*(1-SunF)]   # light-dependent
        Ea2L_k = GAMTLI(SunT)*SunF + GAMTLI(ShaT)*(1-SunF)       # light-independent
        GAMTP  = sum_k [ (Ea1L_k*LDF + Ea2L_k*(1-LDF)) * VPGWT_k ]

    Per-layer inputs (sun_t, sha_t, sun_p, sha_p, sun_f) carry a trailing axis of
    length ``layers``; ``lai``/``d_temp``/``d_ppfd``/``ldf`` are per cell.
    """
    _, vpgwt = sp.canopy_nodes_weights(layers)
    cdea = gamma_cd(lai, layers)                     # (..., layers)
    d_temp_b = _arr(d_temp)[..., None]
    d_ppfd_b = _arr(d_ppfd)[..., None]
    ldf_b = _arr(ldf)[..., None]

    ea1 = cdea * (
        gamtld(sun_t, d_temp_b, s) * gamp(sun_p, d_ppfd_b) * sun_f
        + gamtld(sha_t, d_temp_b, s) * gamp(sha_p, d_ppfd_b) * (1.0 - sun_f)
    )
    ea2 = gamtli(sun_t, s) * sun_f + gamtli(sha_t, s) * (1.0 - sun_f)
    integrand = (ea1 * ldf_b + ea2 * (1.0 - ldf_b)) * vpgwt
    return integrand.sum(axis=-1)


# ---------------------------------------------------------------------------
# Leaf-age response
# ---------------------------------------------------------------------------
def gamma_a(laip, laic, d_temp, s: int, nlai: int = 46):
    """Leaf-age activity factor (megvea.f SUBROUTINE GAMMA_A).

    Handles the three regimes (growth / steady / senescence). Division by LAIc is
    guarded so LAIc == 0 yields GAMLA = 1 (moot because emission is scaled by LAIc).
    """
    laip = _arr(laip)
    laic = _arr(laic)
    tt = _arr(d_temp)
    t = _F(sp.TSTLEN[nlai])

    ratio = np.divide(laip, laic, out=np.ones_like(laic), where=laic > 0)

    # growth branch (LAIp < LAIc)
    ti = np.where(tt <= 303.0, 5.0 + 0.7 * (300.0 - tt), 2.9)
    tm = 2.3 * ti
    fnew_g = np.where(ti >= t, 1.0 - ratio, (ti / t) * (1.0 - ratio))
    fmat_g = np.where(tm >= t, ratio, ratio + ((t - tm) / t) * (1.0 - ratio))
    fgro_g = 1.0 - fnew_g - fmat_g
    fold_g = np.zeros_like(ratio)

    # senescence branch (LAIp > LAIc)
    fold_s = np.divide(laip - laic, laip, out=np.zeros_like(laip), where=laip > 0)
    fmat_s = 1.0 - fold_s
    fnew_s = np.zeros_like(ratio)
    fgro_s = np.zeros_like(ratio)

    grow = laip < laic
    sen = laip > laic
    # steady default (LAIp == LAIc)
    fnew = np.where(grow, fnew_g, np.where(sen, fnew_s, 0.0))
    fgro = np.where(grow, fgro_g, np.where(sen, fgro_s, 0.1))
    fmat = np.where(grow, fmat_g, np.where(sen, fmat_s, 0.8))
    fold = np.where(grow, fold_g, np.where(sen, fold_s, 0.1))

    return fnew * sp.ANEW[s] + fgro * sp.AGRO[s] + fmat * sp.AMAT[s] + fold * sp.AOLD[s]


# ---------------------------------------------------------------------------
# CO2 inhibition (isoprene only)
# ---------------------------------------------------------------------------
def gamma_co2(co2_ppm=sp.CO2_STD):
    """CO2 inhibition factor (megvea.f SUBROUTINE GAMMA_CO2). ==1 at CO2=400."""
    co2 = _arr(co2_ppm)
    ci = 0.7 * co2
    g = sp.ISMAX - (sp.ISMAX * ci ** sp.CO2H) / (sp.CSTAR ** sp.CO2H + ci ** sp.CO2H)
    return np.where(co2 == 400.0, 1.0, g)


# ---------------------------------------------------------------------------
# Bidirectional LAI response (class 13 only)
# ---------------------------------------------------------------------------
def gamma_laibidir(lai):
    """Bidirectional-exchange LAI factor (megvea.f SUBROUTINE GAMMA_LAIbidir)."""
    lai = _arr(lai)
    return np.where(lai < 2.0, 0.5 * lai,
                    np.where(lai <= 6.0, 1.0 - 0.0625 * (lai - 2.0), 0.75))


# ---------------------------------------------------------------------------
# Stress responses (default OFF; piecewise-linear ramps)
# ---------------------------------------------------------------------------
def gamma_ht(maxt, s: int):
    """High-temperature stress (megvea.f SUBROUTINE GAMMA_HT). maxt in K."""
    maxt = _arr(maxt)
    thtk = 273.15 + sp.THT[s]
    t1 = thtk + sp.DTHT[s]
    ramp = 1.0 + (sp.CHT[s] - 1.0) * (maxt - thtk) / sp.DTHT[s]
    return np.where(maxt <= thtk, 1.0, np.where(maxt < t1, ramp, sp.CHT[s]))


def gamma_lt(mint, s: int):
    """Low-temperature stress (megvea.f SUBROUTINE GAMMA_LT). mint in K."""
    mint = _arr(mint)
    tltk = 273.15 + sp.TLT[s]
    t1 = tltk - sp.DTLT[s]
    ramp = 1.0 + (sp.CLT[s] - 1.0) * (tltk - mint) / sp.DTLT[s]
    return np.where(mint >= tltk, 1.0, np.where(mint > t1, ramp, sp.CLT[s]))


def gamma_hw(maxws, s: int):
    """High-wind stress (megvea.f SUBROUTINE GAMMA_HW). maxws in m/s."""
    maxws = _arr(maxws)
    t1 = sp.THW[s] + sp.DTHW[s]
    ramp = 1.0 + (sp.CHW[s] - 1.0) * (maxws - sp.THW[s]) / sp.DTHW[s]
    return np.where(maxws <= sp.THW[s], 1.0, np.where(maxws < t1, ramp, sp.CHW[s]))


def gamma_aq(aqi, s: int):
    """Air-quality (W126) stress (megvea.f SUBROUTINE GAMMA_AQ)."""
    aqi = _arr(aqi)
    t1 = sp.TAQ[s] + sp.DTAQ[s]
    ramp = 1.0 + (sp.CAQ[s] - 1.0) * (aqi - sp.TAQ[s]) / sp.DTAQ[s]
    return np.where(aqi <= sp.TAQ[s], 1.0, np.where(aqi < t1, ramp, sp.CAQ[s]))
