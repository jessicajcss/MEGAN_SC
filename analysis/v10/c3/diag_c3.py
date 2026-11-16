import os, sys, glob, json, re
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
src = P.notebook_source('11_verify_run.py')
# take the configuration and helper definitions only (between the settings cell and the main loop)
start = src.index('import json, yaml\nimport numpy as np, pandas as pd, rasterio'); end = src.index("years = sorted(int(os.path.basename(p))")
code = src[start:end]
code = re.sub(r'^# %%.*$', '', code, flags=re.M)
NEW_ROOT = P.RUN
SC_GPKG = os.path.join(P.BND, 'sc_state_ibge.gpkg')
import common as C                                    # RHO_ISOP of notebook 11 (cell 1) = common.RHO; C3 is a ratio test
g = dict(REPO_DIR=P.REPO, RUN_CONFIG=os.path.join('config', 'santa_catarina_2003_2025.yaml'), NEW_ROOT=NEW_ROOT,
         BOUNDARY=SC_GPKG, REPORT_DIR=None, os=os, glob=glob, json=json, RHO_ISOP=C.rho('ISOP'))
exec(code, g)
read_annual, south_first_tif, sc_mask, nbr_median, cell_area = (g[k] for k in ('read_annual','south_first_tif','region_mask','nbr_median','cell_area'))
from scipy import stats
years = [int(a) for a in sys.argv[1:]] or list(range(2003, 2026))
rows = []
for y in years:
    D = os.path.join(NEW_ROOT, str(y))
    E, (lat, lon), _ = read_annual(NEW_ROOT, y)
    SC = sc_mask(lat, lon); A = cell_area(lat, lon)
    filled = np.zeros_like(SC); logs = sorted(glob.glob(os.path.join(D, 'inputs', f'fixlog_{y}_??.npz')))
    nfill = np.zeros(SC.shape)
    for f in logs:
        fm = np.load(f)['filled_met_steps']; filled |= fm > 0; nfill += fm
    Z = np.load(os.path.join(D, 'inputs', f'ef_{y}.npz')); EF = Z['ef'].astype('float64'); ecl = [str(c) for c in Z['classes']]
    laif = sorted(glob.glob(os.path.join(D, 'inputs', f'lai_{y}_[01][0-9].tif')))
    lai = np.mean([south_first_tif(f) for f in laif], axis=0)
    laim = nbr_median(lai, SC & ~filled & (lai > 0))
    match = (lai > 0) & (laim > 0) & (np.abs(np.log(np.where(lai > 0, lai, 1) / np.where(laim > 0, laim, 1))) < np.log(1.25))
    for c in ('ISOP', 'MT_PINE', 'MEOH'):
        ef_c = EF[ecl.index(c)]; okc = SC & (ef_c > 0) & (E[c] > 0)
        lg = np.log(np.where(okc, E[c], 1.0) / np.where(okc, ef_c * A, 1.0)); lg[~okc] = np.nan
        d = lg - nbr_median(lg, okc & ~filled)
        v = okc & ~filled & match & np.isfinite(d); fsel = okc & filled & match & np.isfinite(d)
        thr = np.percentile(np.abs(d[v]), 99.5); n = int(fsel.sum()); lim = int(stats.binom.ppf(0.99, n, 0.005)) if n else 0
        for j, i in zip(*np.where(fsel)):
            rows.append(dict(year=y, cls=c, j=int(j), i=int(i), lat=round(float(lat[j]), 4), lon=round(float(lon[i]), 4), dlog=float(d[j, i]),
                             thr=float(thr), exceed=bool(abs(d[j, i]) > thr), n=n, lim=lim, filled_steps=int(nfill[j, i]),
                             lai=float(lai[j, i]), lai_nbr=float(laim[j, i]),
                             valid_q50=float(np.median(np.abs(d[v]))), valid_q99=float(np.percentile(np.abs(d[v]), 99)), valid_max=float(np.abs(d[v]).max())))
    print(y, 'done', flush=True)
out = P.out('diagnostics', 'c3', 'c3_filled_cells_2003_2025.csv' if years == list(range(2003, 2026)) else
            'c3_cells_' + '_'.join(map(str, years)) + '.csv')
pd.DataFrame(rows).to_csv(out, index=False); print('wrote', out)
