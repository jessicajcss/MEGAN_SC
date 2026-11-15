"""Soil NO (Yienger-Levy path) and GAMSM: Python port vs official MEGSEA (BDSNP_YN=N), SE_US 2018-05-02.
usage: bench_soil.py <out csv>; env ENGINE_SRC (megan_sc src), MGN (MEGANv3.21 dir), FD (Fortran output dir), MEGAN_PATCH=0|1"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.environ.get('ENGINE_SRC', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from megan_sc.engine import soil
if os.environ.get('MEGAN_PATCH', '0') == '1':
    from megan_sc import fixes; print('engine patches:', fixes.apply_engine_patches())
from ioapi import rd, compare
M = os.environ.get('MGN', os.path.expanduser('~/scratch/megan321/MEGANv3.21')); FD = os.environ.get('FD', 'se_us/fort'); MET = f'{M}/Input/MGNMET/MET.MEGAN.SE_US.rad45.2018122.ncf'
L = f'{M}/Input/MAP/LAI3_SE_US.ncf'; lat = rd(L, 'LAT')[0]; lai = rd(L, 'LAI16')[0]; cts = rd(f'{M}/work/CANTYPFILE', 'CTS')
D = 2018122; rows = []
gday, _ = soil.growseason(D, lat, corrected=False)
for path in ('soil', 'nosoil'):
    ref = rd(f'{FD}/MGNSEA_YL_{path}.ncf', 'GAMNO'); py = [dict() for _ in range(24)]
    for t in range(24):
        T = rd(MET, 'TEMP2', t); pa = rd(MET, 'PREC_ADJ', t)
        kw = dict(soilt=rd(MET, 'SOIT1', t), soilm=rd(MET, 'SOIM1', t), soiltype=np.rint(rd(MET, 'SLTYP', t)).astype(int)) if path == 'soil' else {}
        for corr in (False, True):
            cfno, cfnog = soil.soilnox(T, lai, lat, D, precadj=pa, corrected=corr, **kw)
            g = soil.gamma_no(cfno, cfnog, cts, soil.growseason(D, lat, corrected=corr)[0])
            py[t][corr] = g
    for corr in (False, True):
        rows.append(dict(test=f'GAMNO YL {path}-path, corrected={corr}', **compare(np.stack([p[corr] for p in py]), ref)))
gsm_ref = rd(f'{FD}/MGNSEA_YL_soil.ncf', 'GAMSM')
gsm_py = np.stack([soil.gamma_sm(rd(MET, 'SOIM1', t), soil.wilting_point(np.rint(rd(MET, 'SLTYP', t)).astype(int))) for t in range(24)])
rows.append(dict(test='GAMSM (WWLT by soil type)', **compare(gsm_py, gsm_ref)))
d = pd.DataFrame(rows); print(d.to_string(float_format=lambda v: f'{v:.3g}')); d.to_csv(sys.argv[1] if len(sys.argv) > 1 else 'se_us/bench_soil.csv', index=False)
