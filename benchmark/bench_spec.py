"""Speciation + EF application + unit chain: Python (speciation.py) vs official mgn2mech (RUN_CONVERSION=N, SPCTONHR=Y).
usage: bench_spec.py <out csv>; env ENGINE_SRC, MGN, FD as in bench_soil.py"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.environ.get('ENGINE_SRC', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from megan_sc.engine import speciation as spc, species as sp
from ioapi import rd
from netCDF4 import Dataset
M = os.environ.get('MGN', os.path.expanduser('~/scratch/megan321/MEGANv3.21')); FD = os.environ.get('FD', 'se_us/fort')
EF = np.stack([rd(f'{M}/Input/MAP/EFMAP.2019b.SE_US.ncf', f'EF_{n}')[0] for n in sp.NAMES])
AREA = 9000.0 * 9000.0
out = Dataset(f'{FD}/MGN201_tonhr.ncf'); names = spc.SPCA_NAMES; fvars = [k for k in out.variables][1:]; assert all(a[:16] == b for a, b in zip(names, fvars))
ERS = {n: rd(f'{FD}/MGNERS.ncf', n) for n in sp.NAMES}; GNO = rd(f'{FD}/MGNSEA_YL_nosoil.ncf', 'GAMNO')
sum_py = np.zeros(len(names)); sum_ft = np.zeros(len(names)); maxrel = np.zeros(len(names))
for t in range(24):
    er = np.stack([ERS[n][t] for n in sp.NAMES])
    fl = spc.speciate_to_compounds(er, EF, gamno=GNO[t])
    ton = spc.compound_tonnes(spc.flux_to_mol_s(fl, AREA), 3600.0)
    for i, n in enumerate(names):
        f = np.asarray(out[fvars[i]][t, 0], np.float64); p = ton[i]
        sum_py[i] += p.sum(); sum_ft[i] += f.sum()
        big = np.abs(f) > 1e-6 * max(np.abs(f).max(), 1e-30)
        if big.any(): maxrel[i] = max(maxrel[i], float(np.max(np.abs(p[big] - f[big]) / np.abs(f[big]))))
d = pd.DataFrame(dict(species=names, parent=[sp.NAMES[c] if c >= 0 else '-' for c in spc.SPEC_CLASS], tonnes_py=sum_py, tonnes_fortran=sum_ft, max_rel_cellhour=maxrel))
d['ratio'] = d.tonnes_py / d.tonnes_fortran.replace(0, np.nan)
d.to_csv(sys.argv[1] if len(sys.argv) > 1 else 'se_us/bench_speciation.csv', index=False)
cl = d.groupby('parent')[['tonnes_py', 'tonnes_fortran']].sum(); cl['ratio'] = cl.tonnes_py / cl.tonnes_fortran
print(cl.to_string(float_format=lambda v: f'{v:.6g}'))
print('species with |ratio-1|>1e-5:', d[(d.ratio - 1).abs() > 1e-5][['species', 'parent', 'tonnes_py', 'tonnes_fortran', 'ratio']].to_string())
print('max cell-hour rel diff over all species: %.3g' % d.max_rel_cellhour.max(), '| species zero in both:', int(((d.tonnes_py == 0) & (d.tonnes_fortran == 0)).sum()))
