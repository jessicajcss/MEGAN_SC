"""Read WRF-ARW output (``wrfout``) into MEGAN met variables on the analysis grid.

This is the ``meteorology.source == "wrf"`` path — an alternative to ERA5-Land
(``gee/era5.py``). It lets the ported engine run on the *same* meteorology used by a
reference inventory (e.g. Hoinaski et al. 2024, whose daily ``wrfout_d01`` files drive
their MEGAN v3.1 run), so a model-vs-model comparison isolates differences in the
emission model from differences in the driving met.

Variable mapping (WRF-ARW standard names → the fields the engine consumes)
--------------------------------------------------------------------------
    TEMP2    2 m air temperature (K)            <- T2
    PAR      photosynth. active radiation (W m-2) <- par_fraction * SWDOWN
    WINDSPD  wind speed (m s-1)                 <- hypot(U10, V10)
    PRES     surface pressure (Pa)              <- PSFC
    QV       water-vapour mixing ratio (kg/kg)  <- Q2   (already a mixing ratio)
    SOILM    volumetric soil water (m3/m3)      <- SMOIS, top layer

Deliberate choices, kept consistent with the ERA5 path and documented so nothing is
silently assumed (see docs/ALGORITHM_SPEC.md):

* **PAR = 0.45 x SWDOWN.** Same broadband-shortwave → PAR fraction used for ERA5.
  Unlike ERA5-Land's accumulated radiation, WRF ``SWDOWN`` is an *instantaneous* surface
  downward shortwave flux (W m-2), so there is **no de-accumulation** — a cleaner input.
* **QV = Q2 directly.** WRF ``Q2`` is the 2 m water-vapour mixing ratio (kg kg-1), exactly
  what MEGCAN wants; no dewpoint→mixing-ratio conversion is needed (ERA5 needed one).
* **SOILM = SMOIS[top layer].** Noah/Noah-MP top layer is 0–10 cm; ERA5-Land layer 1 is
  0–7 cm. The small depth mismatch is noted but not corrected — both are the shallow layer
  the soil-moisture activity factor (GAMSM) expects.
* **Regridding.** WRF runs on a projected, *curvilinear* grid (2-D XLAT/XLONG). Each field
  is interpolated to the regular lat/lon analysis grid with :func:`scipy.interpolate.griddata`
  (linear, nearest-neighbour fallback outside the convex hull). Because the reference d01 is
  ~20 km and the SC analysis grid is finer, this is a downscaling interpolation — it adds no
  information beyond WRF's resolution, which is stated plainly rather than implied.

The reader tolerates the common WRF variable-name spellings (aliases below) but does **not**
invent data: a required field that is absent raises, rather than being filled with a guess.
"""
from __future__ import annotations

import datetime as _dt

import numpy as np

# WRF-ARW variable names, with the occasional alias seen in post-processed files.
# First present name wins. Values are required unless listed in _OPTIONAL.
_ALIASES = {
    "T2":     ("T2", "T2M"),
    # SWDOWN and SWDNB are both *downward* surface shortwave; GSW is *net* (down×(1-albedo))
    # and is deliberately excluded so PAR is never silently biased low by surface albedo.
    "SWDOWN": ("SWDOWN", "SWDNB"),
    "U10":    ("U10", "U10M"),
    "V10":    ("V10", "V10M"),
    "PSFC":   ("PSFC", "PSF"),
    "Q2":     ("Q2", "Q2M", "QV2M"),
    "SMOIS":  ("SMOIS", "SOILM", "SH2O"),
    "XLAT":   ("XLAT", "XLAT_M", "CLAT"),
    "XLONG":  ("XLONG", "XLONG_M", "CLONG"),
    # vegetation fields (used by read_wrf_vegetation — the no-Earth-Engine path)
    "IVGTYP": ("IVGTYP", "LU_INDEX"),
    "LAI":    ("LAI",),
    "VEGFRA": ("VEGFRA",),
    "XLAND":  ("XLAND",),
}
_OPTIONAL = ()   # all of the above are required for a MEGAN run

PAR_FRACTION = 0.45   # PAR (W m-2) ~= 0.45 * broadband shortwave-down; matches the ERA5 path


# --------------------------------------------------------------------------- #
# Time handling
# --------------------------------------------------------------------------- #
def _decode_times(ds):
    """Parse the WRF ``Times`` char array into a list of naive UTC datetimes."""
    raw = np.ma.getdata(ds.variables["Times"][:])   # (Time, DateStrLen) of |S1 (or char codes)
    out = []
    for row in np.asarray(raw):
        if row.dtype.kind == "S":
            s = row.tobytes().decode("ascii", errors="ignore")
        elif row.dtype.kind == "U":
            s = "".join(row.tolist())
        else:                                        # numeric character codes
            s = bytes(int(c) for c in np.asarray(row).ravel()).decode("ascii", errors="ignore")
        s = s.strip().replace("\x00", "")
        out.append(_dt.datetime.strptime(s, "%Y-%m-%d_%H:%M:%S"))
    return out


def list_times(path):
    """Return the list of timestamps (naive UTC ``datetime``) in a wrfout file."""
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        return _decode_times(ds)


def _select_time_index(times, when=None, index=None):
    """Resolve a requested time to a Times index (nearest match if not exact)."""
    if index is not None:
        if not (-len(times) <= index < len(times)):
            raise IndexError(f"time index {index} out of range (n={len(times)})")
        return index % len(times)
    if when is None:
        return 0
    if isinstance(when, str):
        when = _dt.datetime.fromisoformat(when)
    deltas = [abs((t - when).total_seconds()) for t in times]
    return int(np.argmin(deltas))


# --------------------------------------------------------------------------- #
# Variable access
# --------------------------------------------------------------------------- #
def _resolve(ds, key):
    for name in _ALIASES[key]:
        if name in ds.variables:
            return name
    raise KeyError(
        f"WRF file has none of {_ALIASES[key]} for '{key}'. "
        f"Available: {sorted(ds.variables)[:40]}...")


def _field(ds, key, ti):
    """A 2-D (south_north, west_east) slice of a WRF variable at time index ``ti``.

    Handles the usual rank layouts:
      (Time, sn, we)                  -> [ti]
      (Time, soil_layers, sn, we)     -> [ti, 0]   (top soil layer)
      (sn, we)                        -> as-is (static field, e.g. coords)
    """
    v = ds.variables[_resolve(ds, key)]
    a = np.asarray(v[:], dtype="float64")
    if a.ndim == 2:
        return a
    if a.ndim == 3:                      # (Time, sn, we)
        return a[ti]
    if a.ndim == 4:                      # (Time, layer, sn, we) -> top layer
        return a[ti, 0]
    raise ValueError(f"unexpected rank {a.ndim} for {key}")


def _coords(ds, ti):
    lat = _field(ds, "XLAT", ti)
    lon = _field(ds, "XLONG", ti)
    return lat, lon


# --------------------------------------------------------------------------- #
# Regridding: curvilinear WRF -> regular analysis grid
# --------------------------------------------------------------------------- #
def interp_to_points(field, src_lon, src_lat, tgt_lon2d, tgt_lat2d):
    """Interpolate a 2-D curvilinear field onto arbitrary target points (2-D arrays).

    Linear (barycentric) interpolation via Delaunay triangulation of the WRF cell
    centres, with a nearest-neighbour fill for any target point outside the source
    convex hull. Target arrays may be a regular mesh OR a curvilinear (e.g. Lambert)
    reference grid — this is what lets us run on Hoinaski's SC grid unchanged.
    Returns an array shaped like ``tgt_lon2d``.
    """
    from scipy.interpolate import griddata

    pts = np.column_stack([np.asarray(src_lon).ravel(), np.asarray(src_lat).ravel()])
    vals = np.asarray(field, dtype="float64").ravel()
    good = np.isfinite(vals)
    pts, vals = pts[good], vals[good]

    TLon = np.asarray(tgt_lon2d, dtype="float64")
    TLat = np.asarray(tgt_lat2d, dtype="float64")
    tgt = np.column_stack([TLon.ravel(), TLat.ravel()])

    out = griddata(pts, vals, tgt, method="linear")
    holes = ~np.isfinite(out)
    if holes.any():                      # edges outside the WRF convex hull
        out[holes] = griddata(pts, vals, tgt[holes], method="nearest")
    return out.reshape(TLon.shape)


def interp_to_grid(field, src_lon, src_lat, tgt_lon, tgt_lat):
    """Interpolate onto a regular (tgt_lat x tgt_lon) grid — 1-D axes convenience wrapper."""
    TLon, TLat = np.meshgrid(np.asarray(tgt_lon), np.asarray(tgt_lat))
    return interp_to_points(field, src_lon, src_lat, TLon, TLat)


def _target_mesh(domain):
    """(TLon, TLat) 2-D target arrays from a domain dict.

    Uses the 2-D ``LON``/``LAT`` if present (supports curvilinear reference grids);
    otherwise meshes the 1-D ``lon``/``lat`` (the usual regular analysis grid).
    """
    if domain.get("LON") is not None and np.ndim(domain["LON"]) == 2:
        return np.asarray(domain["LON"]), np.asarray(domain["LAT"])
    return np.meshgrid(np.asarray(domain["lon"]), np.asarray(domain["lat"]))


# --------------------------------------------------------------------------- #
# Public reader
# --------------------------------------------------------------------------- #
def read_wrf_met(path, domain, *, when=None, index=None, par_fraction=PAR_FRACTION,
                 regrid=True):
    """Read one WRF timestep into the MEGAN met dict on the analysis grid.

    Parameters
    ----------
    path : str
        wrfout NetCDF file (unzipped).
    domain : dict
        As built by :func:`megan_sc.domain.build_domain`; uses ``lon``/``lat`` (1-D
        cell centres). If ``regrid`` is False the WRF native grid is returned instead.
    when : str | datetime, optional
        Target time; the nearest ``Times`` entry is used. Mutually informative with
        ``index`` (``index`` wins if given).
    index : int, optional
        Time index into the file (default 0 / nearest to ``when``).
    par_fraction : float
        Shortwave→PAR fraction (default 0.45).
    regrid : bool
        Interpolate onto the analysis grid (default). False returns fields on the WRF
        grid plus ``lat2d``/``lon2d`` — useful for diagnostics.

    Returns
    -------
    dict with keys ``temp, par, wind, pres, qv, soilm`` (2-D float64) and ``time``
    (the selected ``datetime``). ``lai`` is intentionally NOT set here — LAI comes from
    MODIS (``gee/lai.py``), same as the ERA5 path; the caller merges it in.
    """
    from netCDF4 import Dataset

    with Dataset(path) as ds:
        times = _decode_times(ds)
        ti = _select_time_index(times, when=when, index=index)

        t2 = _field(ds, "T2", ti)
        sw = np.maximum(_field(ds, "SWDOWN", ti), 0.0)
        u10 = _field(ds, "U10", ti)
        v10 = _field(ds, "V10", ti)
        psfc = _field(ds, "PSFC", ti)
        q2 = np.maximum(_field(ds, "Q2", ti), 0.0)
        smois = np.clip(_field(ds, "SMOIS", ti), 0.0, 1.0)
        slat, slon = _coords(ds, ti)

    native = {
        "temp": t2,
        "par": par_fraction * sw,
        "wind": np.hypot(u10, v10),
        "pres": psfc,
        "qv": q2,
        "soilm": smois,
        "time": times[ti],
    }
    if not regrid:
        native["lat2d"], native["lon2d"] = slat, slon
        return native

    TLon, TLat = _target_mesh(domain)
    out = {"time": times[ti]}
    for k in ("temp", "par", "wind", "pres", "qv", "soilm"):
        out[k] = interp_to_points(native[k], slon, slat, TLon, TLat)
    return out


def read_wrf_vegetation(path, domain, *, when=None, index=0, crosswalk=None,
                        mask_lai=True):
    """Growth-form fractions + LAI from the WRF file itself — the **no-Earth-Engine** path.

    Everything MEGAN needs about the land surface is already in the wrfout file: ``IVGTYP``
    (MODIFIED_IGBP_MODIS_NOAH land-use class) gives the vegetation type, and ``LAI`` gives
    the leaf area. Because this is the *same* land surface the reference MEGAN run was driven
    with (via MCIP), it makes the comparison a cleaner apples-to-apples test than pulling
    MODIS from Earth Engine — and it removes the GEE authentication step entirely.

    Returns a dict on the analysis grid: ``lai`` and the growth-form fractions
    ``tree/shrub/herb/crop`` (+ ``veg`` = their sum), ready for
    ``run.build_ef_grid_multiform`` and ``landcover.canopy_type_fractions``.

    ``crosswalk`` defaults to :data:`landcover.IGBP_MODIS_CROSSWALK`. ``mask_lai`` zeroes LAI
    over non-vegetated cells (water/urban/bare), protecting data quality.
    """
    from netCDF4 import Dataset

    from .. import landcover

    cw = landcover.IGBP_MODIS_CROSSWALK if crosswalk is None else crosswalk
    with Dataset(path) as ds:
        times = _decode_times(ds)
        ti = _select_time_index(times, when=when, index=index)
        ivg = np.rint(_field(ds, "IVGTYP", ti)).astype(int)     # dominant land-use class
        lai = np.maximum(_field(ds, "LAI", ti), 0.0)
        slat, slon = _coords(ds, ti)

    w = landcover.growth_form_weights(ivg, cw)                  # (ny, nx, 4) on WRF grid
    TLon, TLat = _target_mesh(domain)
    out = {}
    for g, i in landcover._GF_IDX.items():
        out[g] = np.clip(interp_to_points(w[..., i], slon, slat, TLon, TLat), 0.0, 1.0)
    out["veg"] = sum(out[g] for g in landcover.GROWTH_FORMS)

    lai_g = np.clip(interp_to_points(lai, slon, slat, TLon, TLat), 0.0, None)
    out["lai"] = lai_g * (out["veg"] > 0) if mask_lai else lai_g
    return out


def met_provider(path_for, domain, *, par_fraction=PAR_FRACTION):
    """Build a ``run_period``-compatible met provider backed by WRF files + MODIS LAI.

    ``path_for(year, month, doy, hour) -> wrfout_path`` locates the file for a step;
    the returned callable reads that file at the requested hour and returns the met
    dict (LAI must be injected by the caller, e.g. via a MODIS closure, exactly as in
    the ERA5 notebook).
    """
    def provider(year, month, doy, hour):
        date = _dt.date(year, 1, 1) + _dt.timedelta(days=int(doy) - 1)
        when = _dt.datetime(date.year, date.month, date.day, int(hour))
        return read_wrf_met(path_for(year, month, doy, hour), domain,
                            when=when, par_fraction=par_fraction)
    return provider
