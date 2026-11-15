#!/usr/bin/env python3
"""s1_compare.py - test S1 of notebook 15 on the test runs of notebook 10 (TEST_DAYS = 7, months 1, 4, 7, 10 of 2019):
3-hourly sampling (VARIANT 'base', production setting) against hourly sampling (VARIANT 'hourly').
Both runs: fixes.py sha a933b29602929a82, same inputs, D_TEMP = 24-h mean in both; D_PPFD acts only as a night gate.
Outputs (reproduced/notebook15/s1/; published copies in results/v10/notebook15/s1/): s1_class_totals.csv (state totals by class and month, ratio 3-hourly/hourly), s1_groups.csv,
s1_cell_ratios.csv (cell-level ratio percentiles), s1_diurnal_check.csv, S1_isoprene_ratio_2019_28d.png, s1_compare.log.
Usage: python3 s1_compare.py
"""
import os, sys, numpy as np, pandas as pd
from netCDF4 import Dataset
HERE = os.path.dirname(os.path.abspath(__file__))
# MEGAN_SC_S1 (default data/s1_test) = folder holding the two test runs (megan_sc_outputs_v10/ and megan_sc_outputs_v10_hourly/);
# common.py (grid, state mask, cell areas) of analysis/v10
sys.path.insert(0, HERE); import common as c
sys.path.insert(0, os.path.join(HERE)); import _paths as P
T = P.S1; HERE = P.out('notebook15', 's1')   # outputs
RUNS = {'3h': os.path.join(T, 'megan_sc_outputs_v10', '2019'), '1h': os.path.join(T, 'megan_sc_outputs_v10_hourly', '2019')}
MONTHS = [1, 4, 7, 10]
lat, lon = c.grid(); SC = c.sc_mask(lat, lon)
def grid(run, m):
    with Dataset(os.path.join(RUNS[run], 'emissions', f'emis_2019_{m:02d}_test7d.nc')) as d:
        assert np.allclose(d.variables['lat'][:], lat) and np.allclose(d.variables['lon'][:], lon)
        return {k: np.nan_to_num(np.array(d.variables[k][:], dtype='float64')) * c.rho(k) for k in c.CLASSES}   # isoprene times rho
G = {(r, m): grid(r, m) for r in RUNS for m in MONTHS}
GROUPS = {'BVOC': c.BVOC, 'ISOP': ['ISOP'], 'MT': [k for k in c.CLASSES if k.startswith('MT_')], 'SQT': ['SQT_HR', 'SQT_LR'],
          'OVOC': ['MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD'], 'LVOC': ['LVOC'], 'NO': ['NO'], 'CO': ['CO']}
rows = []
for m in MONTHS + ['28d']:
    ms = MONTHS if m == '28d' else [m]
    for k in c.CLASSES:
        a = sum((G[('3h', mm)][k] * SC).sum() for mm in ms); b = sum((G[('1h', mm)][k] * SC).sum() for mm in ms)
        rows.append(dict(month=m, cls=k, t_3h=a, t_1h=b, ratio_3h_1h=a / b if b > 0 else np.nan))
CT = pd.DataFrame(rows); CT.to_csv(os.path.join(HERE, 's1_class_totals.csv'), index=False)
gr = []
for m in MONTHS + ['28d']:
    s = CT[CT.month == m].set_index('cls')
    for g, ks in GROUPS.items():
        a, b = s.loc[ks, 't_3h'].sum(), s.loc[ks, 't_1h'].sum(); gr.append(dict(month=m, group=g, t_3h=a, t_1h=b, ratio_3h_1h=a / b))
GR = pd.DataFrame(gr); GR.to_csv(os.path.join(HERE, 's1_groups.csv'), index=False)
cr = []
for k in ('ISOP', 'MT_PINE', 'SQT_HR', 'MEOH', 'NO', 'CO'):
    a = sum(G[('3h', m)][k] for m in MONTHS); b = sum(G[('1h', m)][k] for m in MONTHS)
    ok = SC & (b > 0) & (b >= np.percentile(b[SC & (b > 0)], 5))       # drop the 5 % smallest emitters (ratio of tiny numbers)
    r = a[ok] / b[ok]
    cr.append(dict(cls=k, n_cells=int(ok.sum()), p01=np.percentile(r, 1), p05=np.percentile(r, 5), median=np.median(r),
                   p95=np.percentile(r, 95), p99=np.percentile(r, 99), min=r.min(), max=r.max()))
    if k == 'ISOP': RATIO = np.where(ok, a / np.where(b > 0, b, 1), np.nan)
CR = pd.DataFrame(cr); CR.to_csv(os.path.join(HERE, 's1_cell_ratios.csv'), index=False)
dc = []
for m in MONTHS:
    d3 = pd.read_csv(os.path.join(RUNS['3h'], 'hourly', f'diurnal_2019_{m:02d}_test7d.csv'), index_col=0)
    d1 = pd.read_csv(os.path.join(RUNS['1h'], 'hourly', f'diurnal_2019_{m:02d}_test7d.csv'), index_col=0)
    common = d3.columns
    rel = (d3[common] - d1[common]).abs() / d1[common].abs().where(d1[common].abs() > 0)
    for k in ('ISOP', 'MT_PINE', 'SQT_HR', 'NO', 'CO'):
        dc.append(dict(month=m, cls=k, max_rel_diff_common_hours=float(np.nanmax(rel.loc[k].values)) if np.isfinite(rel.loc[k].values).any() else 0.0,
                       left_point_3h_over_hourly_sum=3 * d1.loc[k, common].sum() / d1.loc[k].sum()))
DC = pd.DataFrame(dc); DC.to_csv(os.path.join(HERE, 's1_diurnal_check.csv'), index=False)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
fig, ax = plt.subplots(figsize=(8, 5.5)); EXT = [lon[0] - 0.0125, lon[-1] + 0.0125, lat[0] - 0.0125, lat[-1] + 0.0125]
im = ax.imshow(RATIO, origin='lower', cmap='RdBu_r', vmin=0.95, vmax=1.05, extent=EXT); ax.set_xlim(-54.0, -48.2); ax.set_ylim(-29.5, -25.8)
ax.set_aspect(1 / np.cos(np.deg2rad(-27.3))); ax.axis('off')
fig.colorbar(im, ax=ax, shrink=0.8, label='passo de 3 h / passo horário', format=FuncFormatter(lambda x, p: f'{x:.2f}'.replace('.', ',')))
ax.set_title('Isopreno, 28 dias de 2019 (1 a 7 de jan., abr., jul. e out.): passo de 3 h / passo horário', fontsize=10)
fig.savefig(os.path.join(HERE, 'S1_isoprene_ratio_2019_28d.png'), dpi=200, bbox_inches='tight', facecolor='white'); plt.close(fig)
pd.set_option('display.width', 200)
with open(os.path.join(HERE, 's1_compare.log'), 'w') as f:
    for title, df in (('groups (state totals, t over the 7 days of each month)', GR), ('cell-level ratios 3h/1h, 28 days', CR), ('diurnal check', DC)):
        s = f'--- {title}\n{df.round(5).to_string(index=False)}\n'; print(s); f.write(s)
