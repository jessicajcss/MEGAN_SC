"""MEGCAN canopy environment model — pure-Python port of megcan.f.

Converts above-canopy meteorology into the per-layer sunlit/shaded leaf
temperatures and PPFD that megvea needs (SunT, ShaT, SunP, ShaP, SunF).

Faithful translation of ``reference/src/MEGCAN/megcan.f`` and ``MEGCAN.EXT``:
solar geometry (Calcbeta/CalcEccentricity), the Spitters/Goudriaan solar-fraction
split, the two-stream canopy radiation model (CanopyRad/CalcExtCoeff/
CalcRadComponents), and the iterative leaf energy balance (CanopyEB/LeafEB and
its helper functions). Everything is float64 and vectorised over the grid; the
6 canopy types are looped (as in the Fortran) and area-weighted by the canopy
type fractions.
"""
from __future__ import annotations

import numpy as np

from .species import CDEPTH5

_F = np.float64

# ---- constants (MEGCAN.EXT / megcan.f) --------------------------------------
LAYERS = 5
SOLAR_CONSTANT = 1367.0                 # W m-2
WATER_AIR_RATIO = 18.016 / 28.97
SB = 5.67e-8                            # Stefan-Boltzmann, W m-2 K-4
CONVERT_SHADE_PPFD = 4.6
CONVERT_SUN_PPFD = 4.0
PAR_WM2_TO_PPFD = 4.5                   # megcan.f L437 (PPFD = PAR * 4.5)
PPFD_TO_SOLAR = 2.25                    # megcan.f L513 (Solar = PPFD / 2.25)
RPI180 = 57.29578                       # degrees per radian (as used in source)

# ---- canopy characteristics (17 chars x 6 types), MEGCAN.EXT ----------------
# Types: 1 needleleaf, 2 tropical, 3 temperate broadleaf, 4 shrub, 5 herb, 6 crop
CANOPYCHAR = np.array([
    [16.0, 16.0, 16.0, 1.0, 0.5, 1.0],       # 1  canopy depth (m)
    [0.005, 0.05, 0.05, 0.015, 0.01, 0.02],  # 2  leaf width (m)
    [0.1, 0.1, 0.1, 0.1, 0.15, 0.15],        # 3  leaf length (m)
    [24.0, 24.0, 24.0, 2.0, 0.5, 1.0],       # 4  canopy height (m)
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],          # 5  scattering coeff PPFD
    [0.8, 0.8, 0.8, 0.8, 0.8, 0.8],          # 6  scattering coeff near-IR
    [0.057, 0.057, 0.057, 0.057, 0.057, 0.057],  # 7  refl diffuse PPFD
    [0.389, 0.389, 0.389, 0.389, 0.389, 0.389],  # 8  refl diffuse near-IR
    [0.85, 1.1, 0.9, 0.85, 0.7, 0.65],       # 9  clustering coefficient
    [0.95, 0.95, 0.95, 0.95, 0.95, 0.95],    # 10 leaf IR emissivity
    [1.25, 1.25, 1.25, 1.0, 1.25, 1.25],     # 11 stomata/cuticle factor
    [0.06, 0.06, 0.06, 0.06, 0.06, 0.06],    # 12 daytime lapse rate (K m-1)
    [-0.06, -0.06, -0.06, -0.06, -0.06, -0.06],  # 13 nighttime lapse rate
    [700.0, 700.0, 700.0, 700.0, 700.0, 700.0],  # 14 warm humidity change (Pa)
    [150.0, 150.0, 150.0, 150.0, 150.0, 150.0],  # 15 cool humidity change (Pa)
    [0.7, 0.7, 0.7, 0.7, 0.7, 0.7],          # 16 normalized depth wind negligible
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],          # 17 canopy transparency
], dtype=_F)
NRTYP = 6
# 0-based characteristic indices
C_CDEPTH, C_LWIDTH, C_LLENGTH, C_CHEIGHT, C_SCATV, C_SCATN = 0, 1, 2, 3, 4, 5
C_REFLDV, C_REFLDN, C_CLUSTER, C_EPS, C_TRANSPIRE = 6, 7, 8, 9, 10
C_LAPSE_DAY, C_LAPSE_NIGHT, C_HUM_WARM, C_HUM_COOL, C_WIND_NEGL, C_CANTRAN = 11, 12, 13, 14, 15, 16


# ---- solar geometry ---------------------------------------------------------
def calc_beta(day, lat, hour):
    """Solar elevation angle in degrees (megcan.f FUNCTION Calcbeta)."""
    day = np.asarray(day, dtype=_F)
    lat = np.asarray(lat, dtype=_F)
    hour = np.asarray(hour, dtype=_F)
    sindelta = -np.sin(0.40907) * np.cos(6.28 * (day + 10.0) / 365.0)
    cosdelta = np.sqrt(1.0 - sindelta ** 2)
    a = np.sin(lat / RPI180) * sindelta
    b = np.cos(lat / RPI180) * cosdelta
    sinbeta = a + b * np.cos(2.0 * np.pi * (hour - 12.0) / 24.0)
    return np.arcsin(np.clip(sinbeta, -1.0, 1.0)) * RPI180


def calc_eccentricity(day):
    day = np.asarray(day, dtype=_F)
    return 1.0 + 0.033 * np.cos(2.0 * 3.14 * (day - 10.0) / 365.0)


# ---- humidity / thermodynamic helpers --------------------------------------
def water_vap_pres(dens, pres):
    """Water vapour pressure (Pa) from mixing ratio (megcan.f WaterVapPres)."""
    dens = np.asarray(dens, dtype=_F)
    return (dens / (dens + WATER_AIR_RATIO)) * np.asarray(pres, dtype=_F)


def convert_humidity_pa2kgm3(pa, tk):
    return 0.002165 * np.asarray(pa, dtype=_F) / np.asarray(tk, dtype=_F)


def svd_tk(tk):
    """Saturation vapour density (kg m-3) (megcan.f SvdTk)."""
    tk = np.asarray(tk, dtype=_F)
    svp = 10.0 ** ((-2937.4 / tk) - 4.9283 * np.log10(tk) + 23.5518)  # mbar
    return 0.2165 * svp / tk


def lhv(tk):
    """Latent heat of vaporization (J kg-1) (megcan.f LHV)."""
    return 2501000.0 - 2370.0 * (np.asarray(tk, dtype=_F) - 273.15)


def leaf_ir(tk, eps):
    """IR emission (W m-2), both leaf sides (megcan.f LeafIR)."""
    return np.asarray(eps, dtype=_F) * SB * (2.0 * np.asarray(tk, dtype=_F) ** 4)


def res_sc(par):
    """Stomatal resistance s m-1 (megcan.f ResSC). par = PPFD umol m-2 s-1."""
    par = np.asarray(par, dtype=_F)
    scadj = (0.0027 * 1.066 * par) / np.sqrt(1.0 + 0.0027 ** 2 * par ** 2)
    return np.where(scadj < 0.1, 2000.0, 200.0 / np.where(scadj < 0.1, 1.0, scadj))


def leaf_le(tleaf, ambvap, lathv, gh, stomres, transpire_type):
    """Latent-heat flux (megcan.f LeafLE)."""
    leafres = (1.0 / (1.075 * (gh / 1231.0))) + stomres
    vapdeficit = svd_tk(tleaf) - ambvap
    le = transpire_type * (1.0 / leafres) * lathv * vapdeficit
    return np.where(le < 0.0, 0.0, le)


def leaf_blc(ghforced, tdelta, llength):
    """Boundary-layer conductance (megcan.f LeafBLC)."""
    tdelta = np.asarray(tdelta, dtype=_F)
    ghfree = np.where(
        tdelta >= 0.0,
        0.5 * 0.00253 * ((160000000.0 * np.abs(tdelta) / llength ** 3) ** 0.25) / llength,
        0.0,
    )
    return ghforced + ghfree


def leaf_h(tdelta, gh):
    """Sensible-heat flux, both sides (megcan.f LeafH)."""
    return 2.0 * gh * tdelta


def stability(solar, ctype):
    """Canopy temperature lapse rate (megcan.f FUNCTION Stability)."""
    solar = np.asarray(solar, dtype=_F)
    day_lapse = CANOPYCHAR[C_LAPSE_DAY, ctype]
    night_lapse = CANOPYCHAR[C_LAPSE_NIGHT, ctype]
    tb = 500.0
    mid = day_lapse - ((tb - solar) / tb) * (day_lapse - night_lapse)
    return np.where(solar > tb, day_lapse, np.where(solar > 0.0, mid, night_lapse))


# ---- solar fraction split ---------------------------------------------------
def solar_fractions(solar, maxsolar):
    """Split solar into diffuse/beam visible & near-IR (megcan.f SolarFractions)."""
    solar = np.asarray(solar, dtype=_F)
    maxsolar = np.asarray(maxsolar, dtype=_F)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(maxsolar > 0.0, solar / np.where(maxsolar > 0.0, maxsolar, 1.0), 0.5)
    transmis = np.where(maxsolar <= 0.0, 0.5, np.where(maxsolar < solar, 1.0, ratio))
    fracdiff = 0.156 + 0.86 / (1.0 + np.exp(11.1 * (transmis - 0.53)))
    ppfdfrac = 0.55 - transmis * 0.12
    ppfddiffrac = np.minimum(fracdiff * (1.06 + transmis * 0.4), 1.0)
    qv = ppfdfrac * solar
    qdiffv = qv * ppfddiffrac
    qbeamv = qv - qdiffv
    qn = solar - qv
    qdiffn = qn * fracdiff
    qbeamn = qn - qdiffn
    return qdiffv, qbeamv, qdiffn, qbeamn


# ---- two-stream canopy radiation -------------------------------------------
def _calc_ext_coeff(qbeam, scat, kb, kd):
    p = (1.0 - scat) ** 0.5
    reflb = 1.0 - np.exp((-2.0 * ((1.0 - p) / (1.0 + p)) * kb) / (1.0 + kb))
    kbp = kb * p
    kdp = kd * p
    qbeam_absorb = kb * qbeam * (1.0 - scat)
    return reflb, kbp, kdp, qbeam_absorb


def _calc_rad_components(qdiff, qbeam, kdp, kbp, kb, scat, refld, reflb, laidepth):
    qdabs = qdiff * kdp * (1.0 - refld) * np.exp(-kdp * laidepth)
    qsabs = qbeam * ((kbp * (1.0 - reflb) * np.exp(-kbp * laidepth))
                     - (kb * (1.0 - scat) * np.exp(-kb * laidepth)))
    return qdabs, qsabs


def canopy_rad(lai, sinbeta, qbeamv, qdiffv, qbeamn, qdiffn, ctype, layers=LAYERS):
    """Per-layer sunlit fraction and sun/shade PPFD & absorbed radiation.

    Faithful to megcan.f SUBROUTINE CanopyRad. Scalar canopy-char values for
    ``ctype``; grid arrays for lai/sinbeta/q*. Returns arrays shaped (..., layers):
    sunfrac, sun_ppfd, shade_ppfd, sunqv, shadeqv, sunqn, shadeqn.
    """
    nodes = CDEPTH5 if layers == 5 else (np.arange(1, layers + 1) - 0.5) / layers
    scatv = CANOPYCHAR[C_SCATV, ctype]
    scatn = CANOPYCHAR[C_SCATN, ctype]
    refldv = CANOPYCHAR[C_REFLDV, ctype]
    refldn = CANOPYCHAR[C_REFLDN, ctype]
    cluster = CANOPYCHAR[C_CLUSTER, ctype]
    cantran = CANOPYCHAR[C_CANTRAN, ctype]

    lai = np.asarray(lai, dtype=_F)
    sinbeta = np.asarray(sinbeta, dtype=_F)
    laiadj = lai / (1.0 - cantran)

    day = ((qbeamv + qdiffv) > 0.001) & (sinbeta > 0.002) & (laiadj > 0.001)
    sinbeta_safe = np.where(sinbeta > 0.002, sinbeta, 1.0)
    kb = cluster * 0.5 / sinbeta_safe
    kd = 0.8 * cluster
    reflbv, kbpv, kdpv, qbabsv = _calc_ext_coeff(qbeamv, scatv, kb, kd)
    reflbn, kbpn, kdpn, qbabsn = _calc_ext_coeff(qbeamn, scatn, kb, kd)

    day3 = day[..., None]
    shape = lai.shape + (layers,)
    sunfrac = np.empty(shape, dtype=_F)
    sun_ppfd = np.empty(shape, dtype=_F)
    shade_ppfd = np.empty(shape, dtype=_F)
    sunqv = np.empty(shape, dtype=_F)
    shadeqv = np.empty(shape, dtype=_F)
    sunqn = np.empty(shape, dtype=_F)
    shadeqn = np.empty(shape, dtype=_F)

    for i, node in enumerate(nodes):
        laidepth = laiadj * node
        sf = np.exp(-kb * laidepth)
        qdabsvl, qsabsvl = _calc_rad_components(qdiffv, qbeamv, kdpv, kbpv, kb,
                                                scatv, refldv, reflbv, laidepth)
        qdabsnl, qsabsnl = _calc_rad_components(qdiffn, qbeamn, kdpn, kbpn, kb,
                                                scatn, refldn, reflbn, laidepth)
        shp = (qdabsvl + qsabsvl) * CONVERT_SHADE_PPFD / (1.0 - scatv)
        snp = shp + qbabsv * CONVERT_SUN_PPFD / (1.0 - scatv)
        shqv = qdabsvl + qsabsvl
        snqv = shqv + qbabsv
        shqn = qdabsnl + qsabsnl
        snqn = shqn + qbabsn

        d = day3[..., 0]
        sunfrac[..., i] = np.where(d, sf, 0.2)
        sun_ppfd[..., i] = np.where(d, snp, 0.0)
        shade_ppfd[..., i] = np.where(d, shp, 0.0)
        sunqv[..., i] = np.where(d, snqv, 0.0)
        shadeqv[..., i] = np.where(d, shqv, 0.0)
        sunqn[..., i] = np.where(d, snqn, 0.0)
        shadeqn[..., i] = np.where(d, shqn, 0.0)

    return sunfrac, sun_ppfd, shade_ppfd, sunqv, shadeqv, sunqn, shadeqn


# ---- leaf energy balance ----------------------------------------------------
def _leaf_eb(ppfd, q, irin, eps, transpire_type, lwidth, llength,
             tairk, humidairpa, ws):
    """Iterative leaf temperature (megcan.f SUBROUTINE LeafEB). Returns Tleaf (K)."""
    ws_used = np.where(ws <= 0.0, 0.001, ws)
    humidair_kgm3 = convert_humidity_pa2kgm3(humidairpa, tairk)
    ghforced = 0.0259 / (0.004 * np.sqrt(llength / ws_used))
    stomres = res_sc(ppfd)
    irout_airt = leaf_ir(tairk, eps)
    lathv = lhv(tairk)
    lh_airt = leaf_le(tairk, humidair_kgm3, lathv, ghforced, stomres, transpire_type)

    e1 = q + irin - irout_airt - lh_airt
    e1 = np.where(e1 == 0.0, -1.0, e1)

    tdelt = np.ones_like(tairk)
    balance = np.full_like(tairk, 10.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        for _ in range(10):
            active = np.abs(balance) > 2.0
            gh1 = leaf_blc(ghforced, tdelt, llength)
            sh1 = leaf_h(tdelt, gh1)
            lathv_i = lhv(tairk + tdelt)
            lh = leaf_le(tairk + tdelt, humidair_kgm3, lathv_i, gh1, stomres, transpire_type)
            lh1 = lh - lh_airt
            irout = leaf_ir(tairk + tdelt, eps)
            irout1 = irout - irout_airt
            denom = (sh1 + lh1 + irout1) / tdelt
            new_tdelt = np.where(denom != 0.0, e1 / np.where(denom != 0.0, denom, 1.0), tdelt)
            new_balance = q + irin - irout - sh1 - lh
            tdelt = np.where(active, new_tdelt, tdelt)
            balance = np.where(active, new_balance, balance)

    tdelt = np.clip(tdelt, -10.0, 10.0)
    tleaf = tairk + tdelt
    return np.where(np.isfinite(tleaf), tleaf, tairk)


def canopy_eb(trate, tair0, humidairpa0, ws0, sun_ppfd, shade_ppfd,
              sunqv, shadeqv, sunqn, shadeqn, ctype, layers=LAYERS):
    """Per-layer sun/shade leaf temperatures (megcan.f SUBROUTINE CanopyEB)."""
    nodes = CDEPTH5 if layers == 5 else (np.arange(1, layers + 1) - 0.5) / layers
    cdepth = CANOPYCHAR[C_CDEPTH, ctype]
    lwidth = CANOPYCHAR[C_LWIDTH, ctype]
    llength = CANOPYCHAR[C_LLENGTH, ctype]
    cheight = CANOPYCHAR[C_CHEIGHT, ctype]
    eps = CANOPYCHAR[C_EPS, ctype]
    transpire = CANOPYCHAR[C_TRANSPIRE, ctype]
    hum_warm = CANOPYCHAR[C_HUM_WARM, ctype]
    hum_cool = CANOPYCHAR[C_HUM_COOL, ctype]
    wind_negl = CANOPYCHAR[C_WIND_NEGL, ctype]

    tair0 = np.asarray(tair0, dtype=_F)
    deltah = np.where(
        tair0 > 288.0, hum_warm / cheight,
        np.where(tair0 > 278.0,
                 (hum_warm - ((288.0 - tair0) / 10.0) * (hum_warm - hum_cool)) / cheight,
                 hum_cool / cheight),
    )

    shape = tair0.shape + (layers,)
    sun_tk = np.empty(shape, dtype=_F)
    shade_tk = np.empty(shape, dtype=_F)
    log_ref = np.log(cheight - wind_negl * cheight)

    for i, node in enumerate(nodes):
        ldepth = cdepth * node
        tairk = tair0 + trate * ldepth
        humidairpa = humidairpa0 + deltah * ldepth
        wsh = (cheight - ldepth) - (wind_negl * cheight)
        wsh_safe = np.where(wsh < 0.001, 1.0, wsh)
        ws = ws0 * np.log(wsh_safe) / log_ref
        ws = np.where(wsh < 0.001, 0.05, ws)

        emissatm = 0.642 * (humidairpa / tairk) ** (1.0 / 7.0)
        irin = leaf_ir(tairk, emissatm)

        sun_tk[..., i] = _leaf_eb(sun_ppfd[..., i], sunqv[..., i] + sunqn[..., i],
                                  irin, eps, transpire, lwidth, llength,
                                  tairk, humidairpa, ws)
        shade_tk[..., i] = _leaf_eb(shade_ppfd[..., i], shadeqv[..., i] + shadeqn[..., i],
                                    irin, eps, transpire, lwidth, llength,
                                    tairk, humidairpa, ws)
    return sun_tk, shade_tk


# ---- top-level driver -------------------------------------------------------
def megcan(temp, par_wm2, wind, pres, qv, lat, lon, lai, ctf, day, hour_utc,
           layers=LAYERS):
    """Run the canopy model for one timestep (megcan.f main loop).

    Parameters (all grids unless noted)
    -----------------------------------
    temp    : 2 m air temperature (K)
    par_wm2 : incoming PAR (W m-2)
    wind    : wind speed (m s-1)
    pres    : surface pressure (Pa)
    qv      : water-vapour mixing ratio (kg/kg)
    lat,lon : latitude/longitude (deg)
    lai     : leaf area index (m2/m2)
    ctf     : canopy-type fractions, shape (6, ...) in PERCENT (0-100)
    day     : Julian day (scalar int)
    hour_utc: UTC hour of day (scalar float)

    Returns dict of arrays shaped (..., layers): SunT, ShaT, SunP, ShaP, SunF.
    """
    temp = np.asarray(temp, dtype=_F)
    ppfd = np.asarray(par_wm2, dtype=_F) * PAR_WM2_TO_PPFD
    lat = np.asarray(lat, dtype=_F)
    lon = np.asarray(lon, dtype=_F)
    lai = np.asarray(lai, dtype=_F)
    ctf = np.asarray(ctf, dtype=_F)

    grid = temp.shape
    totalct = np.sum(ctf * 0.01, axis=0)

    # solar geometry (solar hour + day rollover, as in megcan.f L500-503)
    hour = hour_utc + lon / 15.0
    day_arr = np.where(hour < 0.0, day - 1, day).astype(_F)
    hour = np.where(hour < 0.0, hour + 24.0, hour)
    beta = calc_beta(day_arr, lat, hour)
    sinbeta = np.sin(beta / RPI180)
    solar = ppfd / PPFD_TO_SOLAR
    maxsolar = sinbeta * SOLAR_CONSTANT * calc_eccentricity(day_arr)
    qdiffv, qbeamv, qdiffn, qbeamn = solar_fractions(solar, maxsolar)
    humidairpa0 = water_vap_pres(qv, pres)

    out = grid + (layers,)
    sun_tk_tot = np.zeros(out, dtype=_F)
    shade_tk_tot = np.zeros(out, dtype=_F)
    sun_ppfd_tot = np.zeros(out, dtype=_F)
    shade_ppfd_tot = np.zeros(out, dtype=_F)
    sun_frac_tot = np.zeros(out, dtype=_F)

    for ict in range(min(NRTYP, ctf.shape[0])):
        w = 0.01 * ctf[ict]
        if not np.any(w > 0.0):
            continue
        sf, snp, shp, snqv, shqv, snqn, shqn = canopy_rad(
            lai, sinbeta, qbeamv, qdiffv, qbeamn, qdiffn, ict, layers)
        trate = stability(solar, ict)
        stk, shtk = canopy_eb(trate, temp, humidairpa0, wind, snp, shp,
                              snqv, shqv, snqn, shqn, ict, layers)
        w3 = w[..., None]
        sun_tk_tot += w3 * stk
        shade_tk_tot += w3 * shtk
        sun_ppfd_tot += w3 * snp
        shade_ppfd_tot += w3 * shp
        sun_frac_tot += w3 * sf

    valid = ((totalct > 0.0) & (lai > 0.0))[..., None]
    tct = np.where(totalct > 0.0, totalct, 1.0)[..., None]
    temp3 = temp[..., None]
    ppfd3 = ppfd[..., None]
    return {
        "SunT": np.where(valid, sun_tk_tot / tct, temp3),
        "ShaT": np.where(valid, shade_tk_tot / tct, temp3),
        "SunP": np.where(valid, sun_ppfd_tot / tct, ppfd3),
        "ShaP": np.where(valid, shade_ppfd_tot / tct, ppfd3),
        "SunF": np.where(valid, sun_frac_tot / tct, 1.0),
    }
