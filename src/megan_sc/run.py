"""Config-driven runner: domain -> EF -> met -> LAI -> engine -> outputs.

Ties the whole pipeline together. Real meteorology and LAI come from the Earth
Engine layer (gee/); until those exports are available on a given machine, the
``synthetic_met`` helper lets the full chain run for a pipeline demonstration
(clearly labelled). Emission factors are the real ecotype composition.
"""
from __future__ import annotations

import os
import csv
import json
import calendar
import datetime as _dt
import numpy as np
import yaml

from . import domain as dom
from . import efproc
from .engine import canopy, megvea, soil, speciation as spc, species as sp


# --------------------------------------------------------------------------- #
def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def _load_ecotype(path):
    from netCDF4 import Dataset
    nc = Dataset(path)
    evt = np.asarray(nc.variables["EVT"][:])
    flon = np.asarray(nc.variables["lon"][:])
    flat = np.asarray(nc.variables["lat"][:])
    nc.close()
    return evt, flon, flat


def build_ef_grid(config, domain, growth_form="tree"):
    """Gridded per-class EF (19, nlat, nlon) on the analysis grid.

    Composes ecotype -> species -> EF at native resolution, then block-mean
    regrids to the analysis grid. NOTE: uses a single growth form (default tree).
    The production run (notebook 10) uses build_ef_grid_multiform with the
    MapBiomas growth-form fractions. Emission factors are those of the MEGANv3.21
    tables and are not calibrated against any other inventory.
    """
    ef_cfg = config["datasets"]["emission_factors"]
    evt, flon, flat = _load_ecotype(ef_cfg["ecotype"])
    ef_tab = efproc.load_ef_table(ef_cfg["ef_table"])
    by_eco = efproc.ef_by_ecotype(ef_tab, ef_cfg["speciation_dir"], growth_form)
    ef_fine = efproc.map_ecotype_ef(evt, by_eco)                     # (19, ny_f, nx_f)
    ef_grid = dom.regrid_block_mean(ef_fine, flon, flat,
                                    domain["lon"], domain["lat"], domain["resolution_deg"])
    return np.nan_to_num(ef_grid, nan=0.0)


def build_ef_grid_multiform(config, domain, gf_fractions):
    """Growth-form-weighted gridded EF (19, nlat, nlon) — the coherent MapBiomas path.

    ``gf_fractions``: dict tree/shrub/herb/crop -> (nlat, nlon) fractions on the
    analysis grid (from MapBiomas via Earth Engine, or GFv7nc). EF for each cell =
    Σ_growthform  frac · (ecotype -> species -> EF for that growth form). Keeps the
    MEGAN ecotype/species EF science; MapBiomas only supplies the growth-form mix.
    """
    from .landcover import GROWTH_FORMS
    ef_cfg = config["datasets"]["emission_factors"]
    evt, flon, flat = _load_ecotype(ef_cfg["ecotype"])
    ef_tab = efproc.load_ef_table(ef_cfg["ef_table"])
    total = np.zeros((19, domain["lat"].size, domain["lon"].size), dtype=np.float64)
    for gf in GROWTH_FORMS:
        by_eco = efproc.ef_by_ecotype(ef_tab, ef_cfg["speciation_dir"], gf)
        ef_fine = efproc.map_ecotype_ef(evt, by_eco)
        ef_g = dom.regrid_block_mean(ef_fine, flon, flat, domain["lon"],
                                     domain["lat"], domain["resolution_deg"])
        total += np.nan_to_num(gf_fractions[gf], nan=0.0)[None, :, :] * np.nan_to_num(ef_g, nan=0.0)
    return total


def synthetic_met(domain, *, temp=298.0, par=450.0, wind=2.0, pres=1.0e5,
                  qv=0.012, soilm=0.30, lai=4.0):
    """SYNTHETIC met/LAI for a pipeline demo (NOT real weather). Uniform fields."""
    shp = (domain["lat"].size, domain["lon"].size)
    g = lambda v: np.full(shp, float(v))
    return {
        "temp": g(temp), "par": g(par), "wind": g(wind), "pres": g(pres),
        "qv": g(qv), "soilm": g(soilm), "lai": g(lai),
    }


def run_timestep(domain, ef, met, *, day, hour_utc, canopy_type_fractions,
                 d_temp=None, d_ppfd=None, wilt=0.12, gamno=None, date_yyyyddd=None,
                 soiltype=None):
    """One timestep: canopy -> ER -> speciation -> per-compound flux (nanomol m-2 s-1).

    Soil NO (class 8): if ``gamno`` is not supplied and ``date_yyyyddd`` is given,
    it is computed from the Yienger-Levy soil-NO routine (soil path when ``soiltype``
    + soil moisture are available, else the air-temperature path).

    Returns dict with 'compounds' (201, ny, nx), 'classes' (19, ny, nx),
    'er' (19, ny, nx), and 'gamno'.
    """
    lat, lon = domain["LAT"], domain["LON"]

    # Sanitise the met inputs. ERA5-Land is a LAND-only reanalysis, so ocean / no-data cells in
    # the SC bounding box come back masked or as a large fill value; those overflow the
    # activity-factor exponentials (exp of a huge temperature) and then contaminate domain-total
    # aggregates (0*inf -> NaN). We replace non-finite values with 0 and clip temperature to a
    # physical range. Land cells (valid met) are unaffected — their temperatures never approach
    # these bounds — so results over land are identical; only ocean/no-data cells are made a
    # clean zero-emission (they carry LAI=0 anyway). This also silences the RuntimeWarnings.
    met = {k: (np.nan_to_num(np.asarray(v, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
               if k in ("temp", "par", "wind", "pres", "qv", "soilm", "lai") and v is not None else v)
           for k, v in met.items()}
    met["temp"] = np.clip(met["temp"], 150.0, 340.0)          # K; physical bounds

    lai = met["lai"]
    can = canopy.megcan(met["temp"], met["par"], met["wind"], met["pres"],
                        met["qv"], lat, lon, lai, canopy_type_fractions, day, hour_utc)
    d_temp = met["temp"] - 2.0 if d_temp is None else d_temp
    d_ppfd = met["par"] * 4.5 * 0.5 if d_ppfd is None else d_ppfd   # crude daily-mean PPFD
    gamsm = soil.gamma_sm(met["soilm"], wilt)
    er = megvea.emission_activity(can, laic=lai, laip=lai, d_temp=d_temp,
                                  d_ppfd=d_ppfd, gamsm=gamsm)

    if gamno is None and date_yyyyddd is not None:
        kw = {}
        if soiltype is not None and "soilm" in met:
            kw = dict(soilt=met.get("soilt"), soilm=met["soilm"], soiltype=soiltype)
        cfno, cfnog = soil.soilnox(met["temp"], lai, lat, date_yyyyddd, **kw)
        gday, _ = soil.growseason(date_yyyyddd, lat)
        gamno = soil.gamma_no(cfno, cfnog, canopy_type_fractions, gday)

    compounds = spc.speciate_to_compounds(er, ef, gamno=gamno)
    classes = spc.compounds_to_classes(compounds)
    return {"compounds": compounds, "classes": classes, "er": er, "gamno": gamno}


def state_totals(compounds, area, seconds=3600.0):
    """Per-compound domain-total emission (mol s-1 and tonnes over `seconds`)."""
    mol_s = spc.flux_to_mol_s(compounds, area[None, :, :])
    tonnes = spc.compound_tonnes(mol_s, seconds)
    return mol_s.reshape(spc.N_SPCA, -1).sum(1), tonnes.reshape(spc.N_SPCA, -1).sum(1)


# --------------------------------------------------------------------------- #
# Full-period (monthly) inventory driver
# --------------------------------------------------------------------------- #
def months_in_period(period):
    """List of (year, month) tuples from a config period block (inclusive)."""
    y0, m0 = map(int, str(period["start"]).split("-"))
    y1, m1 = map(int, str(period["end"]).split("-"))
    out, y, m = [], y0, m0
    while (y, m) <= (y1, m1):
        out.append((y, m)); m += 1
        if m > 12: m, y = 1, y + 1
    return out


def month_timesteps(year, month, sample_hours):
    """(day_of_year, hour_utc, seconds_represented, date_YYYYDDD) for a month.

    ``sample_hours`` sub-samples the diurnal cycle (e.g. 3 → 8 steps/day); each
    step represents ``sample_hours`` hours for time integration.
    """
    ndays = calendar.monthrange(year, month)[1]
    seconds = sample_hours * 3600.0
    steps = []
    for d in range(1, ndays + 1):
        doy = _dt.date(year, month, d).timetuple().tm_yday
        date = year * 1000 + doy
        for h in range(0, 24, sample_hours):
            steps.append((doy, float(h), seconds, date))
    return steps


def provenance(config, timestamp=None):
    """Reproducibility record written alongside outputs."""
    veg = config["datasets"]["vegetation"]
    return {
        "code_version": "0.1.0",
        "timestamp": timestamp,
        "growth_form_source": veg.get("growth_form_source"),
        "resolution_deg": config["domain"]["resolution_deg"],
        "meteorology": config["datasets"]["meteorology"].get("source"),
        "notes": ("Soil-NO SH growing season = corrected (1 Nov - 31 May) by default. "
                  "Emission factors from the MEGANv3.21 tables, not calibrated against "
                  "other inventories; ERA5-Land radiation de-accumulated as in fixes.py (FIX-8)."),
        "config": config,
    }


def run_period(config, domain, ef, ctf, met_provider, *, wilt=0.12, soiltype=None,
               sample_hours=3, outdir="outputs", write=True, timestamp=None):
    """Accumulate the monthly inventory over the config period.

    ``met_provider(year, month, doy, hour_utc)`` returns a met dict (temp, par, wind,
    pres, qv, soilm, lai) for that timestep — GEE-backed for real runs, synthetic for
    demos. Writes per-month CF-NetCDF (monthly-mean class flux + monthly tonnes) and a
    combined tonnes CSV + provenance. Returns the per-species-per-month tonnes rows.
    """
    from .io import netcdf
    area, lon, lat = domain["area"], domain["lon"], domain["lat"]
    ny, nx = lat.size, lon.size
    rows = []
    for (year, month) in months_in_period(config["period"]):
        steps = month_timesteps(year, month, sample_hours)
        comp_tonnes = np.zeros((spc.N_SPCA, ny, nx), dtype=np.float64)
        flux_classes = np.zeros((sp.NCLASS, ny, nx), dtype=np.float64)
        for (doy, hour, seconds, date) in steps:
            met = met_provider(year, month, doy, hour)
            res = run_timestep(domain, ef, met, day=doy, hour_utc=hour,
                               canopy_type_fractions=ctf, wilt=wilt,
                               date_yyyyddd=date, soiltype=soiltype)
            mol_s = spc.flux_to_mol_s(res["compounds"], area[None, :, :])
            comp_tonnes += spc.compound_tonnes(mol_s, seconds)
            flux_classes += res["classes"]
        mean_flux = flux_classes / max(len(steps), 1)
        class_tonnes = spc.compounds_to_classes(comp_tonnes)
        fam = spc.classes_to_families(class_tonnes)
        for i, nm in enumerate(sp.NAMES):
            rows.append({"year": year, "month": month, "level": "class",
                         "name": nm, "tonnes_per_month": float(class_tonnes[i].sum())})
        for k, v in fam.items():
            rows.append({"year": year, "month": month, "level": "family",
                         "name": k, "tonnes_per_month": float(v.sum())})
        if write:
            tag = f"{year}{month:02d}"
            os.makedirs(f"{outdir}/netcdf", exist_ok=True)
            netcdf.write_emissions(f"{outdir}/netcdf/megan_sc_flux_{tag}.nc", mean_flux,
                                   sp.NAMES, lon, lat, units="nanomol m-2 s-1",
                                   global_attrs={"year": year, "month": month})
            netcdf.write_emissions(f"{outdir}/netcdf/megan_sc_tonnes_{tag}.nc", class_tonnes,
                                   sp.NAMES, lon, lat, units="tonnes/month",
                                   long_name_prefix="Monthly emission of")
    if write:
        os.makedirs(f"{outdir}/tonnes", exist_ok=True)
        with open(f"{outdir}/tonnes/megan_sc_monthly_tonnes.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["year", "month", "level", "name", "tonnes_per_month"])
            w.writeheader(); w.writerows(rows)
        with open(f"{outdir}/provenance.json", "w") as f:
            json.dump(provenance(config, timestamp), f, indent=2, default=str)
    return rows
