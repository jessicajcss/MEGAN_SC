"""Stage-by-stage benchmark of the MEGAN-SC Python port against the official MEGANv3.21 Fortran
(Guenther 2024, Zenodo 10.5281/zenodo.10526206), on the official SE_US test case (2018-05-02).
usage: python bench_se_us.py <megan_sc src dir> <MEGANv3.21 dir> <fortran out dir> <out csv> [GAMSM_YN] [GAMBD_YN]"""
import sys, os, numpy as np, pandas as pd
SRC, MGN, FD, OUTCSV = sys.argv[1:5]; SM = (sys.argv[5] if len(sys.argv) > 5 else 'N') == 'Y'; BD = (sys.argv[6] if len(sys.argv) > 6 else 'N') == 'Y'
sys.path.insert(0, SRC); sys.path.insert(0, os.path.dirname(__file__))
from megan_sc.engine import canopy, megvea, daymet, species as sp
PATCH = os.environ.get('MEGAN_PATCH', '0') == '1'
if PATCH:
    from megan_sc import fixes; print('engine patches:', fixes.apply_engine_patches())
from ioapi import rd, compare
INP = f'{MGN}/Input'; MET = f'{INP}/MGNMET/MET.MEGAN.SE_US.rad45.2018122.ncf'; LAIF = f'{INP}/MAP/LAI3_SE_US.ncf'
JJJ = 122; NT = 24
lat = rd(LAIF, 'LAT')[0]; lon = rd(LAIF, 'LONG')[0]
lai_c = rd(LAIF, f'LAI{int(np.ceil(JJJ / 8.0)):02d}')[0]; lai_p = rd(LAIF, f'LAI{int(np.ceil(JJJ / 8.0)) - 1:02d}')[0]   # findlai.f, NLAI=46
cts = rd(f'{MGN}/work/CANTYPFILE', 'CTS')                       # (6, ROW, COL), percent
met = {k: rd(MET, k)[:NT] for k in ('TEMP2', 'PAR', 'WINDSPD', 'PRES', 'QV')}
rows = []
def add(stage, var, st): rows.append(dict(stage=stage, var=var, **st))
# ---- S1 DAYMET
# daymet.f averages AVEBY = MIN(24, MXREC) records; the official run script (RLENG = 230000) gives 23 (00-22 UTC).
# MinT is not compared: daymet.f allocates TEMP(24,...) and takes MINVAL over all 24 rows, so with 23 records read
# MinT = 0 (unused unless GAMLT_YN = Y).
NAVG = int(os.environ.get('DAYMET_NREC', '23'))
dm = daymet.daily_aggregate(met['TEMP2'][:NAVG], met['PAR'][:NAVG], met['WINDSPD'][:NAVG], axis=0)
for k in ('D_TEMP', 'D_PPFD', 'MaxT', 'MaxWS'):
    add('S1 daymet', k, compare(dm[k], rd(f'{FD}/DAYMET.ncf', k)[0]))
# ---- S2 MEGCAN (Python canopy model vs CANMET), hour by hour
CANV = [('SunleafTK', 'SunT'), ('ShadeleafTK', 'ShaT'), ('SunPPFD', 'SunP'), ('ShadePPFD', 'ShaP'), ('SunFrac', 'SunF')]
can_py = {}; acc = {f: ([], []) for f, _ in CANV}
CANF = {f: rd(f'{FD}/CANMET.ncf', f, squeeze_lay=False).astype(np.float32) for f, _ in CANV}   # (T, LAY, ROW, COL)
ERSF = {n: rd(f'{FD}/MGNERS.ncf', n).astype(np.float32) for n in sp.NAMES}
veg = (lai_c > 0) & (cts.sum(0) > 0)
for t in range(NT):
    c = canopy.megcan(met['TEMP2'][t], met['PAR'][t], met['WINDSPD'][t], met['PRES'][t], met['QV'][t],
                      lat, lon, lai_c, cts, JJJ, float(t))
    can_py[t] = c
    for f, p in CANV:
        acc[f][0].append(np.moveaxis(c[p], -1, 0).astype(np.float32)); acc[f][1].append(CANF[f][t])
for f, _ in CANV:
    add('S2 megcan', f, compare(np.stack(acc[f][0]), np.stack(acc[f][1]), mask=np.broadcast_to(veg, np.stack(acc[f][1]).shape)))
# ---- S3 MEGVEA isolated (Fortran CANMET + DAYMET fed to Python megvea) and S4 end-to-end (Python megcan -> Python megvea)
ldf = np.zeros((19,) + lat.shape)
for s in sp.LDF_FROM_MAP: ldf[s] = rd(f'{INP}/MAP/LDF_SE_US.2019b.ncf', f'LDF{s + 1:02d}')[0]
dT = rd(f'{FD}/DAYMET.ncf', 'D_TEMP')[0]; dP = rd(f'{FD}/DAYMET.ncf', 'D_PPFD')[0]
smf = os.environ.get('SMFILE', f'{MGN}/Output/INT/MGNSEA.SE_US.2018122.ncf'); GSM = rd(smf, 'GAMSM') if SM else None
er3 = []; er4 = []; erf = []
for t in range(NT):
    canf = {p: np.moveaxis(CANF[f][t].astype(np.float64), 0, -1) for f, p in CANV}
    gsm = GSM[t] if SM else None
    kw = dict(laic=lai_c, laip=lai_p, d_temp=dT, d_ppfd=dP, gamsm=gsm, ldf_map=ldf, nlai=46, apply_co2=False, apply_bidir=BD)
    er3.append(megvea.emission_activity(canf, **kw).astype(np.float32))
    er4.append(megvea.emission_activity(can_py[t], **kw).astype(np.float32))
    erf.append(np.stack([ERSF[n][t] for n in sp.NAMES]))
    print('hour', t, flush=True)
er3, er4, erf = np.stack(er3), np.stack(er4), np.stack(erf)        # (T, 19, ROW, COL)
for s, n in enumerate(sp.NAMES):
    if n == 'NO':
        continue                                                   # soil NO is handled by MEGSEA (BDSNP in the reference)
    add('S3 megvea (Fortran canopy in)', n, compare(er3[:, s], erf[:, s]))
    add('S4 end-to-end ER', n, compare(er4[:, s], erf[:, s]))
d = pd.DataFrame(rows); d.to_csv(OUTCSV, index=False)
pd.set_option('display.width', 220)
print(d.groupby('stage').agg(max_rel=('max_rel', 'max'), p99_rel=('p99_rel', 'max'), worst_sum_ratio_dev=('sum_ratio', lambda x: float(np.nanmax(np.abs(x - 1))))).to_string())
print(d.to_string(float_format=lambda v: f'{v:.3g}'))
