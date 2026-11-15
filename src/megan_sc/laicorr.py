"""Correction of the monthly MODIS LAI with the VIIRS anomaly (used for 2024 and 2025).

Why. The monthly LAI of the production run (notebook 10: MODIS MCD15A3H Collection 6.1, main-algorithm retrievals only,
area mean of the 500 m pixels in each 0.025-degree cell, gaps filled from the adjacent months by FIX-4) fell in 2024
and 2025 by more than in any earlier year. The LAI of the VIIRS instrument on Suomi NPP (VNP15A2H Collection 2), whose
anomalies followed those of MODIS in 2013-2023, did not show the decrease (notebook 16;
results/v10/lai_viirs_check). The decrease is therefore treated as a defect of the recent MODIS record.

What. For a corrected year Y, cell c and calendar month m,

    LAI*(c, Y, m) = M(c, m) * r(c, Y, m),        r(c, Y, m) = V(c, Y, m) / V(c, m)

M and V are the reference means (2013-2023) of the MODIS and VIIRS LAI of that cell and calendar month, each computed
over the years in which the cell has a valid (positive) value and kept only where at least MIN_FRAC of the reference
years are valid. The corrected field keeps the level and the seasonal cycle of the MODIS record used in 2003-2023 and
takes the departure of year Y from the reference from VIIRS. Because r is a ratio of two VIIRS values, it does not
depend on the scale factor applied to the VIIRS band.

How r is formed (`method`):
  'cell'    r of each cell;
  'smooth'  sum of V(c, Y, m) over a (2k+1) x (2k+1) window divided by the sum of V(c, m) over the same cells (only
            cells valid in both enter either sum);
  'state'   one ratio per month: sum of V(c, Y, m) over the sum of V(c, m), over the cells of `ref_mask` valid in both;
  'clim'    r = 1 (MODIS reference mean only; the baseline of the hindcast).
Cells where r cannot be formed take the 'state' ratio of the month (fallback). `clip = (lo, hi)` bounds r.
Cells without a MODIS reference mean keep their own MODIS value multiplied by the state factor of the month (state mean
of LAI* over the state mean of the MODIS value, over the cells of `ref_mask` valid in both).

The choice of `method`, `k` and `clip` is made by the leave-one-year-out hindcast of
analysis/v10/review/lai_viirs_correction.py (each year of 2013-2023 predicted from the other ten years and compared
with the MODIS LAI actually used in that year). The settings chosen there are DEFAULT below.
"""
import numpy as np

REF_YEARS = tuple(range(2013, 2024))
MIN_FRAC = 0.5
DOMINANT = 0.6
# Selected by the hindcast (results/v10/lai_viirs_check/lai_correction_hindcast_2013_2023.csv): one ratio per month
# for each zone of growth_form_zones(): tree-, herbaceous- and crop-dominated cells and the remaining cells of the state.
DEFAULT = dict(method='zone', k=0, clip=None)


def growth_form_zones(gf, sc_mask, dominant=DOMINANT):
    """Zones of the selected correction: 1 tree-, 2 herbaceous-, 3 crop-dominated cells (growth-form fraction above
    `dominant`, MapBiomas growth forms of the run), 4 the other cells of sc_mask, 0 outside sc_mask."""
    z = np.where(sc_mask, 4, 0).astype('int32')
    z[sc_mask & (gf['tree'] > dominant)] = 1
    z[sc_mask & (gf['herb'] > dominant)] = 2
    z[sc_mask & (gf['crop'] > dominant)] = 3
    return z


def valid(a):
    """Positive finite values; NaN elsewhere (0 is a masked or missing cell in both products as sampled here)."""
    a = np.asarray(a, dtype='float64')
    return np.where(np.isfinite(a) & (a > 0), a, np.nan)


def monthly_reference(fields, min_frac=MIN_FRAC):
    """Mean over years of one calendar month. fields: 2-D arrays, one per reference year. Returns (mean, n_valid);
    the mean is NaN where fewer than ceil(min_frac * len(fields)) years are valid."""
    A = np.stack([valid(f) for f in fields])
    n = np.isfinite(A).sum(0)
    need = int(np.ceil(min_frac * len(fields)))
    mu = np.where(n > 0, np.nansum(A, 0) / np.maximum(n, 1), np.nan)
    mu[n < need] = np.nan
    return mu, n


def box_sum(a, k):
    """Sum over the (2k+1) x (2k+1) window centred on each cell; NaN counts as 0 and the edges are padded with 0."""
    a0 = np.nan_to_num(np.asarray(a, dtype='float64'), nan=0.0)
    p = np.pad(a0, k)
    out = np.zeros_like(a0)
    for di in range(2 * k + 1):
        for dj in range(2 * k + 1):
            out += p[di:di + a0.shape[0], dj:dj + a0.shape[1]]
    return out


def viirs_ratio(v, vref, method='cell', k=1, ref_mask=None, clip=None, zones=None, min_cells=25):
    """Ratio r of the VIIRS LAI of the month to its reference mean (see the module docstring). With method 'zone', r is
    formed as for 'state' inside each zone of the integer map `zones` (0 = no zone); zones with fewer than `min_cells`
    cells valid in both take the state ratio. Returns (r, fallback, clipped, state_ratio)."""
    v = valid(v); vr = valid(vref)
    both = np.isfinite(v) & np.isfinite(vr)
    sel = both if ref_mask is None else both & ref_mask
    state = float(v[sel].sum() / vr[sel].sum()) if sel.any() else 1.0
    fb = np.zeros(vr.shape, bool)
    if method == 'clim':
        r = np.ones(vr.shape)
    elif method == 'state':
        r = np.full(vr.shape, state)
    elif method == 'zone':
        r = np.full(vr.shape, state); fb = np.ones(vr.shape, bool)
        for z in np.unique(zones[zones > 0]):
            zs = sel & (zones == z)
            if zs.sum() >= min_cells:
                r[zones == z] = v[zs].sum() / vr[zs].sum(); fb[zones == z] = False
    elif method in ('cell', 'smooth'):
        if method == 'cell':
            r = np.where(both, v / np.where(both, vr, 1.0), np.nan)
        else:
            num = box_sum(np.where(both, v, 0.0), k); den = box_sum(np.where(both, vr, 0.0), k)
            r = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
        fb = ~np.isfinite(r)
        r = np.where(fb, state, r)
    else:
        raise ValueError(method)
    cl = np.zeros(vr.shape, bool)
    if clip is not None and method not in ('clim', 'state', 'zone'):
        lo, hi = clip
        cl = (r < lo) | (r > hi)
        r = np.clip(r, lo, hi)
    return r, fb, cl, state


def corrected_lai(mref, r, modis, ref_mask=None):
    """LAI* = mref * r, with the fallback for cells without a MODIS reference mean (module docstring).
    Returns (lai, fallback); lai is 0 where neither the reference mean nor the MODIS value exists."""
    mr = valid(mref); mo = valid(modis)
    out = mr * r
    ok = np.isfinite(out) & np.isfinite(mo)
    sel = ok if ref_mask is None else ok & ref_mask
    f = float(out[sel].sum() / mo[sel].sum()) if sel.any() else 1.0
    fb = ~np.isfinite(out) & np.isfinite(mo)
    out = np.where(fb, mo * f, out)
    return np.where(np.isfinite(out), out, 0.0), fb


def correct_month(modis_ref_fields, viirs_ref_fields, modis_month, viirs_month, ref_mask=None, method=None, k=None,
                  clip='default', min_frac=MIN_FRAC, zones=None):
    """Corrected LAI of one month from the reference fields of the same calendar month (one per reference year) and
    the MODIS and VIIRS fields of the month to correct. Returns (lai, info)."""
    method = DEFAULT['method'] if method is None else method
    k = DEFAULT['k'] if k is None else k
    clip = DEFAULT['clip'] if clip == 'default' else clip
    if method == 'zone' and zones is None:
        raise ValueError("method 'zone' needs the zone map (growth_form_zones)")
    mref, _ = monthly_reference(modis_ref_fields, min_frac)
    vref, _ = monthly_reference(viirs_ref_fields, min_frac)
    r, fb_r, cl, state = viirs_ratio(viirs_month, vref, method, k, ref_mask, clip, zones)
    lai, fb_m = corrected_lai(mref, r, modis_month, ref_mask)
    m = np.ones(lai.shape, bool) if ref_mask is None else ref_mask
    info = dict(method=method, k=k, clip=clip, state_ratio=state, ratio_fallback_cells=int((fb_r & m).sum()),
                clipped_cells=int((cl & m).sum()), modis_fallback_cells=int((fb_m & m).sum()))
    return lai, info


def read_south_first(path):
    """First band of a GeoTIFF written by the run (north-up), returned south first like the run arrays."""
    import rasterio
    with rasterio.open(path) as s:
        return np.flipud(s.read(1).astype('float64'))


def corrected_run_lai(root, cache, sc_mask, years=(2024, 2025), ref_years=REF_YEARS, dominant=DOMINANT):
    """Corrected monthly LAI of `years` with the DEFAULT settings.
    root   folder of the production run (notebook 10): <year>/inputs/lai_<year>_<mm>.tif (the LAI used by the run) and
           <year>/inputs/growthform_*.npz (zones);
    cache  lai_viirs_modis_cache.npz of notebook 16 (keys VIIRS_<year>_<mm>_LAI, same grid, south first);
    sc_mask  Santa Catarina cells (south first).
    Returns ({(year, month): corrected field}, [info of each month])."""
    import os, glob
    C = np.load(cache)
    memo = {}
    def M(y, m):
        if (y, m) not in memo:
            memo[(y, m)] = read_south_first(os.path.join(root, str(y), 'inputs', f'lai_{y}_{m:02d}.tif'))
        return memo[(y, m)]
    def V(y, m):
        return C[f'VIIRS_{y}_{m:02d}_LAI'].astype('float64')
    out, infos = {}, []
    for y in years:
        gf = dict(np.load(sorted(glob.glob(os.path.join(root, str(y), 'inputs', 'growthform_*.npz')))[0]))
        zones = growth_form_zones(gf, sc_mask, dominant)
        for m in range(1, 13):
            lai, info = correct_month([M(r, m) for r in ref_years], [V(r, m) for r in ref_years], M(y, m), V(y, m),
                                      ref_mask=sc_mask, zones=zones)
            info.update(year=y, month=m)
            out[(y, m)] = lai; infos.append(info)
    return out, infos
