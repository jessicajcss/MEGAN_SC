"""Integrity of the production run (data/run_v10, or MEGAN_SC_RUN):
(a) every month accumulated exactly ndays x 8 time steps (counter of cells filled at every step, fixlog_<y>_<m>.npz);
(b) annual_emis_<y>.nc equals the sum of the 12 monthly files;
(c) the checkpoint left in 2025/checkpoints is a prefix of the finished month (no cell above the final value)."""
import os, sys, glob, calendar, numpy as np
from netCDF4 import Dataset
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import _paths as P
R = P.RUN
CL = ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM', 'NO', 'SQT_HR', 'SQT_LR', 'MEOH', 'ACTO', 'ETOH', 'ACID',
      'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO']
def nc(f):
    with Dataset(f) as ds: return {c: np.array(ds.variables[c][:], dtype='float64') for c in CL}
bad, rel = [], 0.0
for y in range(2003, 2026):
    S = {c: 0.0 for c in CL}
    for m in range(1, 13):
        mx = int(np.load(f'{R}/{y}/inputs/fixlog_{y}_{m:02d}.npz')['filled_met_steps'].max())
        if mx != calendar.monthrange(y, m)[1] * 8: bad.append((y, m, mx))
        e = nc(f'{R}/{y}/emissions/emis_{y}_{m:02d}.nc'); S = {c: S[c] + e[c] for c in CL}
    a = nc(f'{R}/{y}/emissions/annual_emis_{y}.nc')
    rel = max(rel, max(abs(a[c].sum() - S[c].sum()) / max(abs(S[c].sum()), 1e-30) for c in CL))
print(f'(a) months with steps != ndays*8: {bad or "none"} (276 months checked)')
print(f'(b) max relative |annual - sum of monthly| over years and classes: {rel:.1e}')
for ck in sorted(glob.glob(f'{R}/*/checkpoints/*.npz')):
    z = np.load(ck, allow_pickle=True); tag = os.path.basename(ck)[:-4]; y = tag[:4]
    fin = nc(f'{R}/{y}/emissions/emis_{tag}.nc'); done = z['done']
    over = sum(int((z['t'][k] > fin[c] * (1 + 1e-5) + 1e-9).sum()) for k, c in enumerate(CL))
    print(f'(c) {tag}: {done.size} steps saved (steps {done.min()}-{done.max()}, {done.size / 8:.0f} days); cells above final month: {over}; '
          f'saved/final ISOP {z["t"][0].sum() / fin["ISOP"].sum():.3f}, MEOH {z["t"][CL.index("MEOH")].sum() / fin["MEOH"].sum():.3f}')
