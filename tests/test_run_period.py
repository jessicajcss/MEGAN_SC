"""Test the full-period monthly inventory driver with a synthetic met provider."""
import os
import numpy as np

from megan_sc import domain as dom, run
from megan_sc.engine import species as sp


def _domain():
    cfg = {"domain": {"mode": "bbox", "resolution_deg": 0.5,
                      "bbox": {"lon_min": -52.0, "lon_max": -50.5,
                               "lat_min": -28.5, "lat_max": -27.0}}}
    return dom.build_domain(cfg), cfg


def test_run_period_accumulates_and_writes(tmp_path):
    domain, dcfg = _domain()
    ny, nx = domain["lat"].size, domain["lon"].size
    ef = np.ones((sp.NCLASS, ny, nx))                 # synthetic EF
    ctf = np.zeros((6, ny, nx)); ctf[2] = 100.0       # broadleaf canopy

    def met_provider(year, month, doy, hour):
        # crude diurnal cycle so day/night differ
        import math
        par = max(0.0, 500.0 * math.sin(math.pi * (hour - 6) / 12)) if 6 <= hour <= 18 else 0.0
        return {"temp": np.full((ny, nx), 295.0 + 5 * (par / 500.0)),
                "par": np.full((ny, nx), par), "wind": np.full((ny, nx), 2.0),
                "pres": np.full((ny, nx), 1.0e5), "qv": np.full((ny, nx), 0.011),
                "soilm": np.full((ny, nx), 0.30), "lai": np.full((ny, nx), 4.0)}

    cfg = {"period": {"start": "2020-01", "end": "2020-01"},
           "domain": dcfg["domain"],
           "datasets": {"vegetation": {"growth_form_source": "tree-only"},
                        "meteorology": {"source": "synthetic"}}}
    outdir = str(tmp_path / "out")
    rows = run.run_period(cfg, domain, ef, ctf, met_provider, sample_hours=6,
                          outdir=outdir, write=True, timestamp="2026-07-18")

    # one row per class + family for the single month
    class_rows = [r for r in rows if r["level"] == "class"]
    assert len(class_rows) == sp.NCLASS
    iso = next(r for r in class_rows if r["name"] == "ISOP")
    assert iso["tonnes_per_month"] > 0.0                       # January (SH summer) emits isoprene

    # outputs written
    assert os.path.exists(f"{outdir}/tonnes/megan_sc_monthly_tonnes.csv")
    assert os.path.exists(f"{outdir}/netcdf/megan_sc_flux_202001.nc")
    assert os.path.exists(f"{outdir}/netcdf/megan_sc_tonnes_202001.nc")
    assert os.path.exists(f"{outdir}/provenance.json")


def test_months_in_period_spans_year_boundary():
    m = run.months_in_period({"start": "2019-11", "end": "2020-02"})
    assert m == [(2019, 11), (2019, 12), (2020, 1), (2020, 2)]


def test_month_timesteps_count():
    steps = run.month_timesteps(2020, 1, sample_hours=6)   # 31 days x 4 steps
    assert len(steps) == 31 * 4
    assert steps[0][3] == 2020001                          # YYYYDDD of Jan 1
