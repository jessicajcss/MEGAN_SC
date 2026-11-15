"""Minimal IOAPI-netCDF helpers for the MEGANv3.21 benchmark (read; write a lat-lon IOAPI file)."""
import numpy as np
from netCDF4 import Dataset

def rd(path, var, t=None, squeeze_lay=True):
    with Dataset(path) as d:
        a = np.asarray(d.variables[var][:], dtype=np.float64)          # (TSTEP, LAY, ROW, COL)
    if t is not None: a = a[t]
    if squeeze_lay and a.ndim >= 3 and a.shape[-3] == 1: a = np.squeeze(a, axis=-3)
    return a

def attrs(path):
    with Dataset(path) as d:
        return {k: d.getncattr(k) for k in d.ncattrs()}

def compare(py, ft, mask=None, tiny=1e-12):
    """Agreement statistics between a Python array and the Fortran reference (same shape)."""
    py = np.asarray(py, np.float64); ft = np.asarray(ft, np.float64)
    if mask is not None:
        mask = np.broadcast_to(mask, py.shape); py = py[mask]; ft = ft[mask]
    d = py - ft
    big = np.abs(ft) > max(tiny, 1e-6 * np.nanmax(np.abs(ft)) if ft.size else tiny)
    rel = np.abs(d[big]) / np.abs(ft[big]) if big.any() else np.array([0.0])
    s_ft = np.nansum(ft); s_py = np.nansum(py)
    return dict(n=int(py.size), max_abs=float(np.nanmax(np.abs(d))) if d.size else 0.0,
                max_rel=float(np.nanmax(rel)), p99_rel=float(np.nanpercentile(rel, 99)),
                frac_rel_gt_1e4=float(np.mean(rel > 1e-4)), sum_ratio=float(s_py / s_ft) if s_ft != 0 else np.nan,
                r=float(np.corrcoef(py.ravel(), ft.ravel())[0, 1]) if py.size > 2 and np.std(ft) > 0 else np.nan)

IOAPI_FIXED = ['IOAPI_VERSION', 'EXEC_ID', 'FTYPE', 'CDATE', 'CTIME', 'WDATE', 'WTIME', 'SDATE', 'STIME', 'TSTEP', 'NTHIK',
               'NCOLS', 'NROWS', 'NLAYS', 'NVARS', 'GDTYP', 'P_ALP', 'P_BET', 'P_GAM', 'XCENT', 'YCENT', 'XORIG', 'YORIG',
               'XCELL', 'YCELL', 'VGTYP', 'VGTOP', 'VGLVLS', 'GDNAM', 'UPNAM', 'VAR-LIST', 'FILEDESC', 'HISTORY']

def write_ioapi(path, template, data, *, sdate, stime, tstep, units=None, lon0=None, lat0=None, dlon=None, dlat=None,
                gdnam='SC_BENCH', nlays=1):
    """Write an IOAPI-3 lat-lon file. data: {var: array (T, LAY, ROW, COL) or (T, ROW, COL)}; template supplies attrs."""
    t = attrs(template)
    names = list(data)
    arrs = {k: (np.asarray(v, np.float32)[:, None] if np.asarray(v).ndim == 3 else np.asarray(v, np.float32)) for k, v in data.items()}
    T, L, R, C = arrs[names[0]].shape
    with Dataset(path, 'w', format='NETCDF3_64BIT_OFFSET') as d:
        d.createDimension('TSTEP', None); d.createDimension('DATE-TIME', 2); d.createDimension('LAY', L)
        d.createDimension('VAR', len(names)); d.createDimension('ROW', R); d.createDimension('COL', C)
        g = dict(t)
        g.update(SDATE=np.int32(sdate), STIME=np.int32(stime), TSTEP=np.int32(tstep), NCOLS=np.int32(C), NROWS=np.int32(R),
                 NLAYS=np.int32(L), NVARS=np.int32(len(names)), GDTYP=np.int32(1), P_ALP=0.0, P_BET=0.0, P_GAM=0.0,
                 XCENT=0.0, YCENT=0.0, XORIG=float(lon0), YORIG=float(lat0), XCELL=float(dlon), YCELL=float(dlat),
                 GDNAM=f'{gdnam:<16s}', UPNAM=f'{"BENCH":<16s}', FILEDESC='MEGAN-SC benchmark input',
                 HISTORY='written by megan_sc benchmark',
                 **{'VAR-LIST': ''.join(f'{n:<16s}' for n in names)})
        g['VGLVLS'] = np.array(t['VGLVLS'], dtype=np.float32)[:2] if np.size(t['VGLVLS']) >= 2 else np.array([0, 0], np.float32)
        for k in IOAPI_FIXED:
            if k in g: d.setncattr(k, g[k])
        tf = d.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        tf.units = '<YYYYDDD,HHMMSS>'; tf.long_name = 'TFLAG'.ljust(16); tf.var_desc = 'Timestep-valid flags'.ljust(80)
        for k in names:
            v = d.createVariable(k, 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'))
            v.long_name = k.ljust(16); v.units = ((units or {}).get(k, 'none')).ljust(16); v.var_desc = k.ljust(80)
        for ti in range(T):
            if tstep == 0:
                date, time = sdate, stime
            else:
                hh = (stime // 10000) + ti * (tstep // 10000); date, time = sdate + hh // 24, (hh % 24) * 10000
            for vi, k in enumerate(names):
                tf[ti, vi, :] = [date, time]
                d.variables[k][ti] = arrs[k][ti]


def drop_vars(src, dst, drop):
    """Copy an IOAPI file without the variables in ``drop`` (header NVARS, VAR-LIST and TFLAG adjusted).
    Used to make the no-soil-data MGNMET for MEGSEA (SOIM1, SOIT1, SLTYP removed)."""
    with Dataset(src) as s, Dataset(dst, 'w', format=s.data_model) as d:
        keep = [v for v in s.variables if v != 'TFLAG' and v not in drop]
        for n, dim in s.dimensions.items():
            d.createDimension(n, None if dim.isunlimited() else (len(keep) if n == 'VAR' else len(dim)))
        for k in s.ncattrs():
            v = s.getncattr(k)
            if k == 'NVARS': v = np.int32(len(keep))
            if k == 'VAR-LIST': v = ''.join(f'{n:<16s}' for n in keep)
            d.setncattr(k, v)
        idx = [list(s.variables)[1:].index(n) for n in keep]
        for n in ['TFLAG'] + keep:
            sv = s.variables[n]; dv = d.createVariable(n, sv.dtype, sv.dimensions)
            dv.setncatts({a: sv.getncattr(a) for a in sv.ncattrs()})
            dv[:] = sv[:][:, idx, :] if n == 'TFLAG' else sv[:]
