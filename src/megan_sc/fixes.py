"""megan_sc.fixes -- corrections to the MEGAN-SC production run (adversarial review, Oct 2026).

Copy this file to  megan-sc/src/megan_sc/fixes.py .  It does not modify any existing module; the
corrected production cells (colab_cells_v10.py) call `run_timestep_fixed` instead of
`run.run_timestep`.  Every change is traceable to a MEGANv3.21 source line or to a diagnosed
artefact in the delivered 2003-2025 output.

FIX-1  Coastal / no-data meteorology (CRITICAL).
       ERA5-Land is land-only.  At cells it does not cover, Earth Engine returns a non-finite or a
       very large fill value.  `run.run_timestep` sanitised these with nan_to_num (PAR -> 0) and a
       clip of temperature to [150, 340] K, assuming such cells carry LAI = 0.  In the delivered
       run, 120 coastal cells inside Santa Catarina (821 km2) carry MODIS LAI and MapBiomas
       vegetation, so they were run every hour with T = 340 K and PAR = 0: ISOP, MBO and CO (LDF = 1)
       were exactly 0 and the light-independent part of every other class was multiplied by
       exp(beta*(340-303.15)) (x40 for monoterpenes, x525 for sesquiterpenes).  This produced
       16-22 % of the state BVOC total and ~75 % of the sesquiterpenes in every year.
       Fix: replace every invalid met value by the value of the nearest valid cell (all bands taken
       from the same donor cell) BEFORE the engine runs, and refuse to run (`strict=True`) if a
       vegetated cell is still invalid.
FIX-2  Daily-mean meteorology (MAJOR).
       `run.run_timestep` defaults to D_TEMP = T(hour) - 2 K and D_PPFD = 0.5*PPFD(hour) when the
       caller passes none; the production notebook passed none.  In MEGANv3.21 D_TEMP and D_PPFD
       are the 24-h means of the simulation day (daymet.f L287-320; read by megvea.f L515-523) and
       D_TEMP enters GAMTLD through Eopt = Cleo*exp(0.1*(T24-297)) and Topt (megvea.f GAMTLD).
       Fix: `run_timestep_fixed` REQUIRES d_temp and d_ppfd; the Colab cells compute D_TEMP as the
       mean of the 24 hourly ERA5-Land 2 m temperatures of the UTC day and D_PPFD from the day's PAR.
FIX-3  Leaf age (MODERATE).
       `run.run_timestep` called megvea with LAIp = LAIc and NLAI = 46, which pins GAMMA_A to its
       steady-state branch.  MEGANv3.21 reads the previous LAI period (megvea.f L585-600, FINDLAI).
       Fix: pass the previous month's LAI and NLAI = 12 (TSTLEN = 30 d).
FIX-5  Ecotype no-data on coastal land (MINOR).
       The MEGAN v3 ecotype map (EVT) carries the no-data code 65535 on a thin coastal strip.  Vegetated
       analysis cells covered only by no-data get EF = 0 and emit nothing (28 SC cells, 192 km2,
       ~0.15 % of state BVOC; e.g. the whole municipality of Governador Celso Ramos).
       Fix: `fill_ecotype_nearest` writes a copy of the ecotype file in which no-data pixels take the
       nearest valid ecotype; point the config's emission_factors.ecotype to that copy.
FIX-4  MODIS LAI monthly gaps (MINOR).
       Vegetated cells with no good-quality MCD15A3H retrieval in a month received LAI = 0 and
       therefore zero emission (0-1.4 % of annual isoprene, concentrated in cloudy summer months).
       Fix: `fill_lai_gaps` fills such cell-months from the adjacent valid months.
FIX-6  Canopy-depth factor in GAMTP (MAJOR for fidelity; found by the Fortran benchmark, Oct 2026).
       megvea.f (L725-731) computes
           Ea1L(K) = CDEA(K)*GAMTLD(SunT)*GAMP(SunP)*SunF + GAMTLD(ShaT)*GAMP(ShaP)*(1-SunF)
       i.e. by Fortran operator precedence CDEA multiplies ONLY the sunlit term.  The Python port
       (gamma.gamtp) multiplied both terms by CDEA.  On the official SE_US test case this changed the
       light-dependent ER by up to 25 % per cell-hour (+2 % domain total).  Fix: `gamtp_megan321`
       reproduces the reference code exactly (benchmark: max relative difference 4e-6).
FIX-7  Soil-NO air-temperature path (MODERATE; found by the Fortran benchmark).
       soilnox.F caps TAIR at 303 K BEFORE both the grass factor and the soil-temperature proxy
       TSOI = 0.72*TAIR + 82.28 (no-soil path).  The port capped only the grass branch, so crop soil NO
       was too high whenever T2 > 30 C.  Fix: `soilnox_megan321` applies the cap as soilnox.F does.
FIX-8  ERA5-Land radiation de-accumulation at 00 and 01 UTC (MODERATE; found by the notebook test, Oct 2026).
       ERA5-Land accumulations run from 00 UTC, and the field valid at 00 UTC holds the 24 h ending at
       00 UTC (ECMWF, ERA5-Land data documentation).  gee/era5.hourly_flux did not subtract anything at
       00 UTC (so PAR at 00 UTC was the previous day's 24-h total / 3600 x 0.45, ~2700 W m-2 in summer) and
       subtracted the 00 UTC field at 01 UTC (PAR clipped to 0).  In SC 00 UTC is night (about 21 h local
       solar time), so leaf light stayed 0, but the canopy stability/lapse term saw the spurious radiation
       (+4 to +8 % light-independent ER in that sample) and the PAR maps/driver series were inflated.
       Fix: `hourly_flux_fixed` subtracts the previous hour always except at 01 UTC.
NOTE   Southern-Hemisphere growing season (deliberate, documented deviation, kept).
       growseason.f compares the day of year (1-366) with MMDD constants (1101, 0101, 0531) in the
       lat < -23 branch, so south of 23 S the Fortran growing season becomes DOY 101-366.  The port's
       soil.growseason(corrected=True) applies the intended 1 Nov - 31 May season.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .engine import canopy, megvea, soil, speciation as spc, gamma as _gamma, species as _sp

MET_KEYS = ("temp", "par", "wind", "pres", "qv", "soilm")
# physically admissible ranges for the ERA5-Land-derived MEGAN fields (SI units used by the engine)
VALID_RANGE = {
    "temp": (150.0, 340.0),       # K (2 m air temperature)
    "par": (0.0, 1000.0),         # W m-2 (0.45 x shortwave down)
    "wind": (0.0, 75.0),          # m s-1
    "pres": (4.0e4, 1.1e5),       # Pa
    "qv": (0.0, 0.06),            # kg kg-1
    "soilm": (0.0, 1.0),          # m3 m-3
}


# --------------------------------------------------------------------------- FIX-1
def invalid_met_mask(met: dict) -> np.ndarray:
    """True where ANY engine met band is non-finite or outside its physical range."""
    bad = None
    for k in MET_KEYS:
        if k not in met or met[k] is None:
            continue
        x = np.asarray(met[k], dtype=np.float64)
        lo, hi = VALID_RANGE[k]
        b = ~np.isfinite(x) | (x < lo) | (x > hi)
        bad = b if bad is None else (bad | b)
    if bad is None:
        raise ValueError("met dict has none of the expected keys " + str(MET_KEYS))
    return bad


def nearest_valid_index(bad: np.ndarray):
    """(ii, jj) index arrays pointing every cell to its nearest valid cell (itself if valid)."""
    if bad.all():
        raise RuntimeError("FIX-1: no valid meteorology anywhere on the grid")
    _, (ii, jj) = ndimage.distance_transform_edt(bad, return_indices=True)
    return ii, jj


def fill_met_nearest(met: dict, bad: np.ndarray | None = None, keys=MET_KEYS) -> tuple[dict, np.ndarray]:
    """Return (filled_met, bad_mask).  All bands of an invalid cell come from ONE donor cell."""
    bad = invalid_met_mask(met) if bad is None else bad
    out = dict(met)
    if not bad.any():
        return out, bad
    ii, jj = nearest_valid_index(bad)
    for k in keys:
        if k not in met or met[k] is None:
            continue
        x = np.array(met[k], dtype=np.float64, copy=True)
        x[bad] = x[ii[bad], jj[bad]]
        out[k] = x
    return out, bad


def fill_field_nearest(field, lo=-np.inf, hi=np.inf) -> np.ndarray:
    """Nearest-valid fill for a single 2-D field (e.g. daily D_TEMP)."""
    x = np.array(field, dtype=np.float64, copy=True)
    bad = ~np.isfinite(x) | (x < lo) | (x > hi)
    if bad.any():
        ii, jj = nearest_valid_index(bad)
        x[bad] = x[ii[bad], jj[bad]]
    return x


def assert_valid_on_vegetation(met: dict, lai) -> None:
    """Raise if any cell with LAI > 0 still has invalid meteorology (guards against FIX-1 regressions)."""
    lai = np.nan_to_num(np.asarray(lai, dtype=np.float64), nan=0.0)
    bad = invalid_met_mask(met) & (lai > 0)
    if bad.any():
        raise RuntimeError(f"FIX-1: {int(bad.sum())} vegetated cells still have invalid meteorology")


# --------------------------------------------------------------------------- FIX-4
def fill_lai_gaps(lai_m, lai_prev, lai_next, veg_frac, veg_min=0.05) -> tuple[np.ndarray, np.ndarray]:
    """Fill vegetated cells whose monthly LAI is missing/zero from the adjacent valid months.

    Uses the mean of the previous and next month where both are valid, otherwise whichever is
    valid.  Cells without vegetation (veg_frac <= veg_min) are left untouched.  Returns
    (filled_lai, filled_mask)."""
    cur = np.nan_to_num(np.asarray(lai_m, dtype=np.float64), nan=0.0)
    p = np.nan_to_num(np.asarray(lai_prev, dtype=np.float64), nan=0.0)
    n = np.nan_to_num(np.asarray(lai_next, dtype=np.float64), nan=0.0)
    gap = (cur <= 0) & (np.asarray(veg_frac) > veg_min) & ((p > 0) | (n > 0))
    both = (p > 0) & (n > 0)
    repl = np.where(both, 0.5 * (p + n), np.where(p > 0, p, n))
    out = np.where(gap, repl, cur)
    return out, gap


# --------------------------------------------------------------------------- FIX-2 helpers
PAR_WM2_TO_PPFD = 4.5      # daymet.f L320 / megcan.f


def daily_means_from_samples(temp_samples, par_samples):
    """D_TEMP (K) and D_PPFD (umol m-2 s-1) from a stack of same-day samples (time on axis 0).

    Prefer a 24-hourly stack for D_TEMP (MEGAN daymet uses 24 records); with 3-hourly samples the
    mean of the 8 values is an unbiased estimator of the daily mean for a smooth diurnal cycle.
    D_PPFD only acts as a night gate in MEGANv3.21 (C1 is constant in GAMP)."""
    t = np.asarray(temp_samples, dtype=np.float64)
    p = np.asarray(par_samples, dtype=np.float64)
    return t.mean(axis=0), p.mean(axis=0) * PAR_WM2_TO_PPFD


# --------------------------------------------------------------------------- corrected timestep
def run_timestep_fixed(domain, ef, met, *, day, hour_utc, canopy_type_fractions,
                       d_temp, d_ppfd, laip, nlai=12, wilt=0.12, date_yyyyddd=None,
                       soiltype=None, strict=True):
    """Drop-in replacement for run.run_timestep with FIX-1, FIX-2 and FIX-3 applied.

    `met` holds temp, par, wind, pres, qv, soilm, lai (current-month LAI, already gap-filled).
    `d_temp`, `d_ppfd` are the daily means of the simulation day (required).
    `laip` is the previous month's LAI (required); `nlai` = 12 for monthly LAI.
    Returns the same dict as run.run_timestep plus 'filled' (cells whose met was replaced).
    """
    if d_temp is None or d_ppfd is None:
        raise ValueError("FIX-2: d_temp and d_ppfd (daily means of the simulation day) are required")
    if laip is None:
        raise ValueError("FIX-3: laip (previous-period LAI) is required")
    lat, lon = domain["LAT"], domain["LON"]
    lai = np.nan_to_num(np.asarray(met["lai"], dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    laip = np.nan_to_num(np.asarray(laip, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)

    met_f, bad = fill_met_nearest(met)                                    # FIX-1
    if strict:
        assert_valid_on_vegetation(met_f, lai)
    d_temp = fill_field_nearest(d_temp, *VALID_RANGE["temp"])             # same rule for daily T
    d_ppfd = np.nan_to_num(np.asarray(d_ppfd, dtype=np.float64), nan=0.0)

    can = canopy.megcan(met_f["temp"], met_f["par"], met_f["wind"], met_f["pres"],
                        met_f["qv"], lat, lon, lai, canopy_type_fractions, day, hour_utc)
    gamsm = soil.gamma_sm(met_f["soilm"], wilt)
    er = megvea.emission_activity(can, laic=lai, laip=laip, d_temp=d_temp, d_ppfd=d_ppfd,
                                  gamsm=gamsm, nlai=nlai)                 # FIX-2, FIX-3

    gamno = None
    if date_yyyyddd is not None:
        kw = {}
        if soiltype is not None and "soilt" in met_f and met_f.get("soilt") is not None:
            kw = dict(soilt=met_f["soilt"], soilm=met_f["soilm"], soiltype=soiltype)
        cfno, cfnog = soil.soilnox(met_f["temp"], lai, lat, date_yyyyddd, **kw)
        gday, _ = soil.growseason(date_yyyyddd, lat)
        gamno = soil.gamma_no(cfno, cfnog, canopy_type_fractions, gday)

    compounds = spc.speciate_to_compounds(er, ef, gamno=gamno)
    classes = spc.compounds_to_classes(compounds)
    return {"compounds": compounds, "classes": classes, "er": er, "gamno": gamno,
            "filled": bad & (lai > 0)}


# --------------------------------------------------------------------------- FIX-5
EVT_NODATA = 65535


def fill_ecotype_nearest(src_nc: str, dst_nc: str, nodata: int = EVT_NODATA, var: str = "EVT") -> int:
    """Copy the ecotype NetCDF replacing no-data pixels by the nearest valid ecotype.

    Returns the number of pixels filled.  Ocean pixels are filled too, which is harmless because
    ocean cells carry no MapBiomas vegetation (zero growth-form fractions -> zero EF weight)."""
    from netCDF4 import Dataset
    import shutil
    shutil.copyfile(src_nc, dst_nc)
    with Dataset(dst_nc, "r+") as nc:
        evt = np.asarray(nc.variables[var][:])
        bad = evt == nodata
        if bad.any():
            ii, jj = nearest_valid_index(bad)
            evt = evt.copy(); evt[bad] = evt[ii[bad], jj[bad]]
            nc.variables[var][:] = evt
        return int(bad.sum())


# --------------------------------------------------------------------------- FIX-6
def gamtp_megan321(lai, sun_t, sha_t, sun_p, sha_p, sun_f, d_temp, d_ppfd, ldf, s, layers=5):
    """GAMTP exactly as megvea.f L725-736 (CDEA multiplies the sunlit light-dependent term only)."""
    g = _gamma
    _, vpgwt = _sp.canopy_nodes_weights(layers)
    cdea = g.gamma_cd(lai, layers)
    dt = np.asarray(d_temp, dtype=np.float64)[..., None]
    dp = np.asarray(d_ppfd, dtype=np.float64)[..., None]
    ldf_b = np.asarray(ldf, dtype=np.float64)[..., None]
    ea1 = (cdea * g.gamtld(sun_t, dt, s) * g.gamp(sun_p, dp) * sun_f
           + g.gamtld(sha_t, dt, s) * g.gamp(sha_p, dp) * (1.0 - sun_f))
    ea2 = g.gamtli(sun_t, s) * sun_f + g.gamtli(sha_t, s) * (1.0 - sun_f)
    return ((ea1 * ldf_b + ea2 * (1.0 - ldf_b)) * vpgwt).sum(axis=-1)


# --------------------------------------------------------------------------- FIX-7
def soilnox_megan321(tair, lai, lat, date_yyyyddd, *, soilt=None, soilm=None, soiltype=None,
                     precadj=1.0, corrected=True):
    """soil.soilnox with the TAIR <= 303 K cap applied before both branches, as in soilnox.F."""
    t = np.minimum(np.asarray(tair, dtype=np.float64), 315.0)
    t = np.minimum(t, 303.0)
    return _ORIG["soilnox"](t, lai, lat, date_yyyyddd, soilt=soilt, soilm=soilm, soiltype=soiltype,
                            precadj=precadj, corrected=corrected)


# --------------------------------------------------------------------------- FIX-8
def deaccumulate(acc_cur, acc_prev, hour_utc):
    """Hourly mean flux (W m-2) over (h-1, h] from ERA5-Land accumulations (J m-2) - numpy reference.

    acc_cur = field valid at hour h, acc_prev = field valid at h-1.  At 01 UTC the field already is the
    1-h accumulation; at 00 UTC it is the 24-h total ending at 00 UTC, so the previous hour (23 UTC of the
    previous day) must be subtracted, like at every other hour."""
    prev = 0.0 if int(hour_utc) == 1 else np.asarray(acc_prev, dtype=np.float64)
    return np.maximum(np.asarray(acc_cur, dtype=np.float64) - prev, 0.0) / 3600.0


def hourly_flux_fixed(ee, region, date):
    """Earth Engine version of `deaccumulate` for surface_solar_radiation_downwards (replaces era5.hourly_flux)."""
    from .gee import era5 as _era5
    d = ee.Date(date)
    col = ee.ImageCollection(_era5.COLLECTION).filterBounds(region)
    band = "surface_solar_radiation_downwards"
    cur = col.filterDate(d, d.advance(1, "hour")).first().select(band)
    prev = ee.Image(col.filterDate(d.advance(-1, "hour"), d).first()).select(band)
    prev_acc = ee.Image(ee.Algorithms.If(ee.Number(d.get("hour")).eq(1), ee.Image(0), prev))
    return cur.subtract(prev_acc).max(0).divide(3600.0).rename("SW_WM2")


_ORIG = {}
def apply_engine_patches():
    """Install FIX-6, FIX-7 and FIX-8 in the running engine (idempotent).  Returns the list of patches."""
    from .gee import era5 as _era5
    if "gamtp" not in _ORIG:
        _ORIG["gamtp"] = _gamma.gamtp
        _ORIG["soilnox"] = soil.soilnox
        _ORIG["hourly_flux"] = _era5.hourly_flux
    _gamma.gamtp = gamtp_megan321          # megvea.emission_activity calls gamma.gamtp at run time
    soil.soilnox = soilnox_megan321        # run_timestep_fixed calls soil.soilnox at run time
    _era5.hourly_flux = hourly_flux_fixed  # era5.megan_met_image calls hourly_flux at run time
    return ["FIX-6 gamtp_megan321", "FIX-7 soilnox_megan321", "FIX-8 hourly_flux_fixed"]


def patches_active() -> bool:
    from .gee import era5 as _era5
    return (_gamma.gamtp is gamtp_megan321 and soil.soilnox is soilnox_megan321
            and _era5.hourly_flux is hourly_flux_fixed)
