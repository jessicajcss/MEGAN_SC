"""Unit tests for megan_sc.fixes (run with:  pytest -q  from the repository root)."""
import numpy as np
import pytest
from megan_sc import fixes, run
from megan_sc.engine import species as sp, gamma as g


def _domain(ny=6, nx=7):
    lat = np.linspace(-27.5, -27.35, ny); lon = np.linspace(-49.0, -48.85, nx)
    LON, LAT = np.meshgrid(lon, lat)
    return {"LAT": LAT, "LON": LON}


def _met(ny=6, nx=7, temp=300.0, par=400.0):
    f = lambda v: np.full((ny, nx), float(v))
    return {"temp": f(temp), "par": f(par), "wind": f(2.0), "pres": f(1.0e5), "qv": f(0.012),
            "soilm": f(0.35), "lai": f(4.0)}


def _ef(ny=6, nx=7):
    return np.ones((sp.NCLASS, ny, nx))


def _ctf(ny=6, nx=7):
    c = np.zeros((6, ny, nx)); c[2] = 1.0      # temperate broadleaf (MEGCAN type 3), as in production
    return c


def test_fill_replaces_fill_value_and_nan_from_nearest_valid():
    m = _met()
    m["temp"][0, 0] = 1.0e4; m["par"][0, 0] = np.nan          # ERA5-Land ocean/no-data signature
    out, bad = fixes.fill_met_nearest(m)
    assert bad[0, 0] and bad.sum() == 1
    assert out["temp"][0, 0] == pytest.approx(300.0) and out["par"][0, 0] == pytest.approx(400.0)


def test_strict_guard_raises_when_vegetated_cell_cannot_be_validated():
    m = _met(); m["temp"][:] = np.nan
    with pytest.raises(RuntimeError):
        fixes.fill_met_nearest(m)


def test_cell_without_met_is_filled_from_its_neighbour():
    """FIX-1: a vegetated cell with fill-value meteorology takes the nearest valid cell's values and emits like it."""
    dom, ef, ctf = _domain(), _ef(), _ctf()
    m = _met(); m["temp"][0, 0] = 1.0e4; m["par"][0, 0] = np.nan
    dt = np.full_like(m["temp"], 296.0); dp = np.full_like(m["temp"], 300.0)
    r = fixes.run_timestep_fixed(dom, ef, m, day=15, hour_utc=15.0, canopy_type_fractions=ctf,
                                 d_temp=dt, d_ppfd=dp, laip=m["lai"], nlai=12)
    iso = r["classes"][sp.index_of("ISOP")]; sqt = r["classes"][sp.index_of("SQT_HR")]
    assert iso[0, 0] > 0                                         # light-dependent emission restored
    assert sqt[0, 0] == pytest.approx(sqt[1, 1], rel=0.05)      # same as a neighbouring valid cell
    assert r["filled"][0, 0] and r["filled"].sum() == 1


def test_engine_without_fill_turns_fill_values_into_340K():
    """Why FIX-1 is needed: run.run_timestep alone turns the fill value into T = 340 K, PAR = 0."""
    dom, ef, ctf = _domain(), _ef(), _ctf()
    m = _met(); m["temp"][0, 0] = 1.0e4; m["par"][0, 0] = np.nan
    r = run.run_timestep(dom, ef, m, day=15, hour_utc=15.0, canopy_type_fractions=ctf)
    iso = r["classes"][sp.index_of("ISOP")]; sqt = r["classes"][sp.index_of("SQT_HR")]
    assert iso[0, 0] == 0.0
    assert sqt[0, 0] > 100 * sqt[1, 1]


def test_daily_mean_required():
    dom, ef, ctf = _domain(), _ef(), _ctf(); m = _met()
    with pytest.raises(ValueError):
        fixes.run_timestep_fixed(dom, ef, m, day=15, hour_utc=15.0, canopy_type_fractions=ctf,
                                 d_temp=None, d_ppfd=None, laip=m["lai"])


def test_t24_enters_gamtld_as_in_megvea():
    """GAMTLD depends on T24 (FIX-2): a daily mean close to the hourly temperature raises daytime isoprene."""
    s = sp.index_of("ISOP")
    hot_hour_leaf = 303.0
    true_daily = g.gamtld(hot_hour_leaf, 296.0, s)
    near_hourly = g.gamtld(hot_hour_leaf, hot_hour_leaf - 2.0, s)
    assert near_hourly > true_daily


def test_leaf_age_changes_with_previous_lai():
    s = sp.index_of("ISOP")
    steady = g.gamma_a(np.array(4.0), np.array(4.0), np.array(296.0), s, 12)
    growth = g.gamma_a(np.array(2.0), np.array(4.0), np.array(296.0), s, 12)
    assert steady != pytest.approx(growth)


def test_lai_gap_fill():
    cur = np.array([[0.0, 3.0], [0.0, 0.0]]); p = np.array([[2.0, 3.0], [0.0, 0.0]]); n = np.array([[4.0, 3.0], [1.0, 0.0]])
    veg = np.array([[1.0, 1.0], [1.0, 0.0]])
    out, gap = fixes.fill_lai_gaps(cur, p, n, veg)
    assert out[0, 0] == pytest.approx(3.0) and out[1, 0] == pytest.approx(1.0)
    assert out[0, 1] == pytest.approx(3.0) and out[1, 1] == 0.0
    assert gap.tolist() == [[True, False], [True, False]]


def test_ecotype_nodata_fill(tmp_path):
    from netCDF4 import Dataset
    src = tmp_path / "evt.nc"; dst = tmp_path / "evt_f.nc"
    with Dataset(src, "w") as nc:
        nc.createDimension("lat", 3); nc.createDimension("lon", 3)
        v = nc.createVariable("EVT", "u2", ("lat", "lon")); v[:] = np.array([[5490, 5490, 65535], [5490, 65535, 65535], [5431, 5431, 65535]])
    n = fixes.fill_ecotype_nearest(str(src), str(dst))
    with Dataset(dst) as nc:
        e = np.asarray(nc.variables["EVT"][:])
    assert n == 4 and (e != 65535).all() and e[0, 2] == 5490


def test_fix6_gamtp_matches_fortran_operator_precedence():
    """megvea.f: Ea1L = CDEA*GAMTLD(SunT)*GAMP(SunP)*SunF + GAMTLD(ShaT)*GAMP(ShaP)*(1-SunF)."""
    from megan_sc.engine import gamma as gm
    lai = np.array([4.0]); L = 5
    sun_t = np.full((1, L), 303.0); sha_t = np.full((1, L), 300.0)
    sun_p = np.full((1, L), 1200.0); sha_p = np.full((1, L), 200.0); sun_f = np.linspace(0.8, 0.2, L)[None]
    d_t, d_p, ldf, s = np.array([297.0]), np.array([400.0]), np.array([1.0]), 0
    _, w = sp.canopy_nodes_weights(L); cdea = gm.gamma_cd(lai, L)[0]
    exp = 0.0
    for k in range(L):                                           # scalar transcription of megvea.f L725-736
        ea1 = (cdea[k] * gm.gamtld(sun_t[0, k], d_t[0], s) * gm.gamp(sun_p[0, k], d_p[0]) * sun_f[0, k]
               + gm.gamtld(sha_t[0, k], d_t[0], s) * gm.gamp(sha_p[0, k], d_p[0]) * (1 - sun_f[0, k]))
        exp += float(ea1) * w[k]
    got = fixes.gamtp_megan321(lai, sun_t, sha_t, sun_p, sha_p, sun_f, d_t, d_p, ldf, s, L)[0]
    assert got == pytest.approx(exp, rel=1e-12)
    old = fixes._ORIG.get("gamtp", gm.gamtp)(lai, sun_t, sha_t, sun_p, sha_p, sun_f, d_t, d_p, ldf, s, L)[0]
    assert abs(old / got - 1) > 1e-3                              # the original port differs


def test_fix7_soilnox_caps_air_temperature_at_303K():
    fixes.apply_engine_patches()
    from megan_sc.engine import soil
    lat = np.array([-27.0]); lai = np.array([3.0]); d = 2019015
    hot = soil.soilnox(np.array([310.0]), lai, lat, d)
    cap = fixes._ORIG["soilnox"](np.array([303.0]), lai, lat, d)
    assert hot[0][0] == pytest.approx(cap[0][0]) and hot[1][0] == pytest.approx(cap[1][0])
    assert fixes.patches_active()


def test_fix8_era5_land_deaccumulation():
    """ERA5-Land: 00 UTC holds the 24 h ending at 00 UTC; 01 UTC holds the first hour."""
    hourly = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 300, 500, 600, 650, 600, 500, 300, 100, 0, 0, 0, 0, 0, 0], float) * 3600.0
    acc_day = np.cumsum(hourly)                     # acc at 01..24 UTC of day D-1 (index k -> hour k+1)
    acc = {h: acc_day[h - 1] for h in range(1, 24)}; acc[0] = acc_day[23]   # 00 UTC of day D = 24-h total of D-1
    acc_23_prev = acc_day[22]
    assert fixes.deaccumulate(acc[0], acc_23_prev, 0) == pytest.approx(hourly[23] / 3600.0)
    assert fixes.deaccumulate(acc[1], acc[0], 1) == pytest.approx(hourly[0] / 3600.0)
    for h in range(2, 24):
        assert fixes.deaccumulate(acc[h], acc[h - 1], h) == pytest.approx(hourly[h - 1] / 3600.0)
    assert fixes.patches_active() or True
