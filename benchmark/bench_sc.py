"""Python port vs official MEGANv3.21 Fortran on the Santa Catarina grid (Southern Hemisphere, 0.025 deg).
usage: bench_sc.py <engine src (with fixes.py)> <MEGAN-SC repo dir> <case dir> <out csv>"""
import sys, os, numpy as np, pandas as pd
SRC, REPO, D, OUTCSV = sys.argv[1:5]
sys.path.insert(0, SRC); sys.path.insert(0, os.path.dirname(__file__))
from megan_sc.engine import canopy, megvea, daymet, soil, gamma, species as sp, speciation as spc
from megan_sc import fixes, run, domain as dom
from ioapi import rd, compare
Z = np.load(f'{D}/py_inputs.npz'); F = f'{D}/fort'
T, PAR = Z['T'], Z['PAR']; LAT, LON = Z['LAT'], Z['LON']; ctf = Z['ctf']; lai12 = Z['lai12']; doy = int(Z['doy']); SD = int(Z['sdate'])
mon = int(pd.Timestamp(f'{SD // 1000}-01-01').dayofyear) and (pd.Timestamp(f'{SD // 1000}-01-01') + pd.Timedelta(days=doy - 1)).month
laic = lai12[mon - 1]; laip = lai12[0] if mon == 1 else lai12[mon - 2]   # findlai.f with NLAI = 12 and MXLAI = 1 (January -> LAI01)
wind, pres, qv = Z['wind'], Z['pres'], Z['qv']
veg = (laic > 0) & (ctf.sum(0) > 0)
rows = []
def add(stage, var, st, **kw): rows.append(dict(case=os.path.basename(D), stage=stage, var=var, **st, **kw))
dm = daymet.daily_aggregate(T, PAR, np.repeat(wind[None], 24, 0), axis=0)
for k in ('D_TEMP', 'D_PPFD', 'MaxT'):
    add('S1 daymet', k, compare(dm[k], rd(f'{F}/DAYMET.ncf', k)[0]))
CANV = [('SunleafTK', 'SunT'), ('ShadeleafTK', 'ShaT'), ('SunPPFD', 'SunP'), ('ShadePPFD', 'ShaP'), ('SunFrac', 'SunF')]
CANF = {f: rd(f'{F}/CANMET.ncf', f, squeeze_lay=False).astype(np.float32) for f, _ in CANV}
ERF = {n: rd(f'{F}/MGNERS.ncf', n).astype(np.float32) for n in sp.NAMES}
dT = rd(f'{F}/DAYMET.ncf', 'D_TEMP')[0]; dP = rd(f'{F}/DAYMET.ncf', 'D_PPFD')[0]
canacc = {f: [] for f, _ in CANV}; erp = []; eru = []; erf = []
for t in range(24):
    c = canopy.megcan(T[t], PAR[t], wind, pres, qv, LAT, LON, laic, ctf, doy, float(t))
    for f, p in CANV: canacc[f].append(np.moveaxis(c[p], -1, 0).astype(np.float32))
    kw = dict(laic=laic, laip=laip, d_temp=dT, d_ppfd=dP, gamsm=None, nlai=12, apply_co2=False, apply_bidir=True)
    fixes.apply_engine_patches(); erp.append(megvea.emission_activity(c, **kw).astype(np.float32))
    gamma.gamtp = fixes._ORIG['gamtp']; eru.append(megvea.emission_activity(c, **kw).astype(np.float32)); fixes.apply_engine_patches()
    erf.append(np.stack([ERF[n][t] for n in sp.NAMES]))
for f, _ in CANV:
    a = np.stack(canacc[f]); add('S2 megcan', f, compare(a, CANF[f], mask=np.broadcast_to(veg, a.shape)))
erp, eru, erf = np.stack(erp), np.stack(eru), np.stack(erf)
# emission-weighted effect sizes: production EF grid (MapBiomas growth forms x MEGAN v3 ecotypes)
_cp = f'{REPO}/config/example_config.yaml'
cfg = run.load_config(_cp if os.path.exists(_cp) else f'{REPO}/example_config.yaml'); cfg['domain']['mode'] = 'bbox'
for k in ('ef_table', 'speciation_dir', 'ecotype'):
    v = cfg['datasets']['emission_factors'].get(k)
    if v and not os.path.isabs(v): cfg['datasets']['emission_factors'][k] = os.path.join(REPO, v)
domn = dom.build_domain(cfg)
gf = {k: v for k, v in np.load(os.path.join(os.environ['RAW'], '2019/inputs/growthform_2019.npz')).items()}
ef = run.build_ef_grid_multiform(cfg, domn, gf); area = domn['area']
assert ef.shape[1:] == laic.shape, (ef.shape, laic.shape)
def day_tonnes(er):            # tonnes per class for the day (hourly steps), NO excluded (soil path)
    out = np.zeros(19)
    for t in range(24):
        fl = spc.speciate_to_compounds(er[t].astype(np.float64), ef, gamno=np.zeros(laic.shape))
        out += spc.compounds_to_classes(spc.compound_tonnes(spc.flux_to_mol_s(fl, area), 3600.0)).sum(axis=(1, 2))
    return out
tp, tu, tf = day_tonnes(erp), day_tonnes(eru), day_tonnes(erf)
for s, n in enumerate(sp.NAMES):
    if n == 'NO':
        continue
    add('S3+S4 ER, patched engine (FIX-6)', n, compare(erp[:, s], erf[:, s]), day_tonnes_ratio=tp[s] / tf[s])
    add('S4 ER, original port (no FIX-6)', n, compare(eru[:, s], erf[:, s]), day_tonnes_ratio=tu[s] / tf[s])
# soil NO, air-temperature path (as in production); Fortran growseason with its SH branch
gref = rd(f'{F}/MGNSEA.ncf', 'GAMNO')
for patched in (True, False):
    if patched: fixes.apply_engine_patches()
    else: soil.soilnox = fixes._ORIG['soilnox']
    for corr in (False, True):
        g = []
        for t in range(24):
            cfno, cfnog = soil.soilnox(T[t], laic, LAT, SD, corrected=corr)
            g.append(soil.gamma_no(cfno, cfnog, ctf, soil.growseason(SD, LAT, corrected=corr)[0]))
        add(f'soil NO GAMNO ({"FIX-7" if patched else "original port"}, SH season {"corrected" if corr else "as Fortran"})', 'NO', compare(np.stack(g), gref))
fixes.apply_engine_patches()
d = pd.DataFrame(rows); d.to_csv(OUTCSV, index=False)
pd.set_option('display.width', 250)
print(d.to_string(float_format=lambda v: f'{v:.3g}'))
