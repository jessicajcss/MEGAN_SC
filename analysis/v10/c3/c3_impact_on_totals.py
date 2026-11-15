"""How much the C3 deviations weigh in the state totals: emission of the exceeding cells, and the change in the state
total if each LAI-matched filled cell had the neighbour-median emission activity (E * exp(-dlog))."""
import os, sys, numpy as np, pandas as pd
from netCDF4 import Dataset
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
R = P.RUN
_f = P.out('diagnostics', 'c3', 'c3_filled_cells_2003_2025.csv')
r = pd.read_csv(_f if os.path.exists(_f) else os.path.join(P.RESULTS, 'diagnostics', 'c3', 'c3_filled_cells_2003_2025.csv'))
T = pd.read_csv(f'{R}/annual_totals_SC_v10_all_classes.csv').set_index('year')
for (y, c) in [(2015, 'MEOH'), (2015, 'MT_PINE')]:
    with Dataset(f'{R}/{y}/emissions/annual_emis_{y}.nc') as ds: E = np.array(ds.variables[c][:], dtype='float64')
    x = r[(r.year == y) & (r.cls == c) & r.exceed]; e = np.array([E[j, i] for j, i in zip(x.j, x.i)])
    print(y, c, len(x), 'cells; emission t', e.round(1), f'share {e.sum() / T.loc[y, c]:.5%}',
          f'change if neighbour-median activity {(e * np.exp(-x.dlog.values) - e).sum() / T.loc[y, c]:+.5%}')
out = []
for y in range(2003, 2026):
    with Dataset(f'{R}/{y}/emissions/annual_emis_{y}.nc') as ds:
        for c in ('ISOP', 'MT_PINE', 'MEOH'):
            E = np.array(ds.variables[c][:], dtype='float64'); x = r[(r.year == y) & (r.cls == c)]
            e = np.array([E[j, i] for j, i in zip(x.j, x.i)]); out.append((e * np.exp(-x.dlog.values) - e).sum() / T.loc[y, c])
print(f'all LAI-matched filled cells at neighbour-median activity: change in state totals {min(out):+.3%} to {max(out):+.3%}')
