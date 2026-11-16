#!/usr/bin/env python3
"""11_report_figures_v10.py
Additional Portuguese figure of the Santa Catarina example that is not a map, from the production run (data/run_v10):
  fig10  standardised regression by class (std_regression_v10.csv, from 10_report_numbers_v10.py).
The additional Portuguese maps (pt_report/fig9, fig12 and figA2 to figA9) are drawn by ../maps/12_maps_wf_style.py.
Output: reproduced/figures/pt_report/fig10.png (published copy in results/v10/figures/pt_report/).
Usage: python3 11_report_figures_v10.py (after 10_report_numbers_v10.py, whose std_regression_v10.csv it reads)
"""
import os, sys, numpy as np, pandas as pd, rasterio, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
HERE = os.path.dirname(os.path.abspath(__file__))
# paths from ../_paths.py (MEGAN_SC_* environment variables); common.py (grid, state mask, cell areas) of analysis/v10
sys.path.insert(0, os.path.join(HERE, '..')); import common as c
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P
OUTD = P.out('figures', 'pt_report')
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.formatter.use_locale": False})
lat, lon = c.grid(); SC = c.sc_mask(lat, lon); A = c.cell_area_m2(lat, lon)
EXT = [lon[0] - 0.0125, lon[-1] + 0.0125, lat[0] - 0.0125, lat[-1] + 0.0125]; ASPECT = 1 / np.cos(np.deg2rad(-27.3))
sc_geom = [g for _, g in c.read_gpkg(f"{c.BND}/sc_state_ibge.gpkg")]
dec = lambda s: s.replace('.', ',').replace('-', '−')
DECF = FuncFormatter(lambda x, p: dec(f"{x:g}"))
def outline(ax):
    for g in sc_geom:
        for p in getattr(g, 'geoms', [g]):
            x, y = p.exterior.xy; ax.plot(x, y, color='0.3', lw=0.5)
def frame(ax):
    ax.set_xlim(-54.0, -48.2); ax.set_ylim(-29.5, -25.8); ax.set_aspect(ASPECT); ax.axis('off')
Z = np.load(P.grid_cache())
def grid_year(y): return {k: np.nan_to_num(Z[f'{k}_{y}'].astype(float)) for k in c.CLASSES}
def groups(g):
    return {'COVB total': sum(g[k] for k in c.BVOC), 'Isopreno': g['ISOP'], 'Monoterpenos': sum(g[k] for k in c.CLASSES if k.startswith('MT_')),
            'Sesquiterpenos': g['SQT_HR'] + g['SQT_LR'], 'COV oxigenados': sum(g[k] for k in ('MEOH', 'ACTO', 'ETOH', 'ACID', 'OXPROD')),
            'HC leves e outros': g['LVOC'] + g['STRESS'] + g['OTHER'] + g['MBO'], 'NO do solo': g['NO'], 'CO': g['CO']}

g25 = None   # the maps (Figures 9, 12 and A.2 to A.9) are drawn by ../maps/12_maps_wf_style.py

# ---- Figure 10: standardised regression
R = pd.read_csv(P.out("report", "std_regression_v10.csv") if os.path.exists(P.out("report", "std_regression_v10.csv")) else os.path.join(P.RESULTS, "report", "std_regression_v10.csv")).set_index('species').loc[c.CLASSES]
lab = {'ISOP': 'Isopreno', 'MBO': 'MBO', 'MT_PINE': 'MT pineno', 'MT_ACYC': 'MT acíclico', 'MT_CAMP': 'MT canfeno', 'MT_SABI': 'MT sabineno',
       'MT_AROM': 'MT aromático', 'NO': 'NO do solo', 'SQT_HR': 'SQT alta reat.', 'SQT_LR': 'SQT baixa reat.', 'MEOH': 'Metanol', 'ACTO': 'Acetona',
       'ETOH': 'Acetald./etanol', 'ACID': 'Ácidos', 'LVOC': 'HC leves', 'OXPROD': 'Prod. oxidação', 'STRESS': 'Estresse', 'OTHER': 'Outros', 'CO': 'CO'}
Bm = R[['beta_TEMP2_C', 'beta_PAR', 'beta_LAI', 'beta_SOILM']].values
fig, ax = plt.subplots(figsize=(8.6, 9.6)); im = ax.imshow(Bm, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
for i in range(Bm.shape[0]):
    for j in range(Bm.shape[1]):
        ax.text(j, i, dec(f"{Bm[i, j]:+.2f}"), ha='center', va='center', fontsize=8.5, color='white' if abs(Bm[i, j]) > 0.6 else 'black')
    ax.text(4.1, i, dec(f"{R.R2.iloc[i]:.2f}"), ha='left', va='center', fontsize=8.5)
ax.text(4.1, -0.6, 'R²', fontweight='bold', va='bottom'); ax.set_xlim(-0.5, 4.8)
ax.set_xticks(range(4)); ax.set_xticklabels(['Temperatura', 'PAR', 'IAF', 'Umid. solo']); ax.set_yticks(range(len(c.CLASSES))); ax.set_yticklabels([lab[k] for k in c.CLASSES])
cb = fig.colorbar(im, ax=ax, shrink=0.5, format=DECF); cb.set_label('Coeficiente padronizado (β)')
ax.set_title("Variáveis de controle das emissões: regressão padronizada por classe\n(2003–2025)", fontsize=11.5, fontweight='bold', pad=22)
fig.text(0.02, 0.01, "Séries anuais em escores z, sem intercepto; PAR e umidade do solo são colineares (r = −0,90).", fontsize=7.5, color='0.4')
fig.tight_layout(rect=[0, 0.03, 1, 1]); fig.savefig(f"{OUTD}/fig10.png", dpi=220, facecolor='white'); plt.close(fig)

print('figures ->', OUTD, sorted(os.listdir(OUTD)))
