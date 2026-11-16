# Post-hoc diagnostic (not the C3 criterion): remove the LAI-contrast effect from dlog gamma with a line fitted on the
# valid cells, then count filled-cell exceedances of the valid-cell 99.5th percentile of the residuals.
import os, sys, glob, json, re
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
src = P.notebook_source('11_verify_run.py')
code = src[src.index('import json, yaml\nimport numpy as np, pandas as pd, rasterio'):src.index("years = sorted(int(os.path.basename(p))")]
code = re.sub(r'^# %%.*$', '', code, flags=re.M)
NEW_ROOT = P.RUN
import common as C                                    # RHO_ISOP of notebook 11 (cell 1) = common.RHO; the test uses ratios
g = dict(REPO_DIR=P.REPO, RUN_CONFIG=os.path.join('config', 'santa_catarina_2003_2025.yaml'), NEW_ROOT=NEW_ROOT,
         BOUNDARY=os.path.join(P.BND, 'sc_state_ibge.gpkg'), REPORT_DIR=None, os=os, glob=glob, json=json, RHO_ISOP=C.rho('ISOP'))
exec(code, g)
rows = []; share = []
for y in range(2003, 2026):
    D = os.path.join(NEW_ROOT, str(y)); E, (lat, lon), _ = g['read_annual'](NEW_ROOT, y)
    SC = g['region_mask'](lat, lon); A = g['cell_area'](lat, lon)
    filled = np.zeros_like(SC)
    for f in glob.glob(os.path.join(D, 'inputs', f'fixlog_{y}_??.npz')): filled |= np.load(f)['filled_met_steps'] > 0
    Z = np.load(os.path.join(D, 'inputs', f'ef_{y}.npz')); EF = Z['ef'].astype('float64'); ecl = [str(c) for c in Z['classes']]
    lai = np.mean([g['south_first_tif'](f) for f in sorted(glob.glob(os.path.join(D, 'inputs', f'lai_{y}_[01][0-9].tif')))], axis=0)
    laim = g['nbr_median'](lai, SC & ~filled & (lai > 0))
    llai = np.log(np.where((lai > 0) & (laim > 0), lai, 1) / np.where((lai > 0) & (laim > 0), laim, 1))
    match = (lai > 0) & (laim > 0) & (np.abs(llai) < np.log(1.25))
    bv = sum(E[c] for c in E if c not in ('NO', 'CO'))
    share.append(dict(year=y, filled_cells=int((SC & filled).sum()), **{f'share_{c}': (E[c] * SC * filled).sum() / (E[c] * SC).sum() for c in ('ISOP', 'MT_PINE', 'MEOH', 'SQT_HR', 'NO')},
                      share_BVOC=(bv * SC * filled).sum() / (bv * SC).sum()))
    for c in ('ISOP', 'MT_PINE', 'MEOH'):
        ef_c = EF[ecl.index(c)]; okc = SC & (ef_c > 0) & (E[c] > 0)
        lg = np.log(np.where(okc, E[c], 1.0) / np.where(okc, ef_c * A, 1.0)); lg[~okc] = np.nan
        d = lg - g['nbr_median'](lg, okc & ~filled)
        v = okc & ~filled & match & np.isfinite(d); fs = okc & filled & match & np.isfinite(d)
        b, a = np.polyfit(llai[v], d[v], 1); res = d - (a + b * llai)
        thr = np.percentile(np.abs(res[v]), 99.5); n = int(fs.sum()); k = int((np.abs(res[fs]) > thr).sum())
        rows.append(dict(year=y, cls=c, slope=b, r_valid=np.corrcoef(llai[v], d[v])[0, 1], n=n, k_resid=k, lim=int(stats.binom.ppf(0.99, n, 0.005)),
                         max_abs_resid_filled=np.abs(res[fs]).max(), thr_resid=thr, sd_llai_valid=llai[v].std(), sd_llai_filled=llai[fs].std()))
R = pd.DataFrame(rows); S = pd.DataFrame(share)
R.to_csv(P.out('diagnostics', 'c3', 'c3_lai_adjusted.csv'), index=False); S.to_csv(P.out('diagnostics', 'c3', 'filled_cells_share.csv'), index=False)
pd.set_option('display.width', 250)
print('tests', len(R), '| fails after LAI adjustment', int((R.k_resid > R.lim).sum()), '| total k', int(R.k_resid.sum()), 'of', int(R.n.sum()), 'cell-tests')
print(R.groupby('cls')[['slope', 'r_valid', 'sd_llai_valid', 'sd_llai_filled', 'max_abs_resid_filled']].agg(['mean', 'max']).round(3))
print(R[R.k_resid > 0][['year', 'cls', 'n', 'k_resid', 'lim', 'max_abs_resid_filled', 'thr_resid']].round(3).to_string())
print(S.describe().loc[['mean', 'min', 'max']].round(4).to_string())
