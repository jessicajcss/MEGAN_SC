#!/usr/bin/env python3
"""run_notebook.py - run one of the notebooks of notebooks/ outside Google Colab.

    python tools/run_notebook.py notebooks/13_santa_catarina_analysis.ipynb "ROOT='data/run_v10';OUT_DIR='reproduced/x'"

The code cells are executed in order in one namespace, as Colab would. Four things differ from Colab:
  * google.colab is replaced by a stub (the Drive mount becomes a no-op) and `%pip` / `!` lines are skipped, so the
    packages of environment.yml must already be installed;
  * the overrides `NAME=python_expression;...` are evaluated and assigned right after the first code cell, which is the
    settings cell of every notebook, so they replace the Colab paths of that cell;
  * matplotlib uses the Agg backend (figures are written to files, plt.show() does nothing);
  * the environment variables MEGAN_SC_RUNNER=1 and MEGAN_SC_REPO=<this repository> tell the setup block of the
    notebooks (notebooks/src/_setup_repo.py) that it runs outside Colab and where the repository is.
Earth Engine is not stubbed: a notebook that needs it (10, or 14 without its cache) stops with the Earth Engine error.
Each cell's run time is printed; the exit code is 1 if a cell fails.
"""
import json, os, sys, time, traceback, types
from unittest import mock

if len(sys.argv) < 2:
    sys.exit(__doc__)
NB = sys.argv[1]
OV = [kv.split('=', 1) for kv in sys.argv[2].split(';') if kv.strip()] if len(sys.argv) > 2 else []
cells = [''.join(c['source']) for c in json.load(open(NB, encoding='utf-8'))['cells'] if c['cell_type'] == 'code']
g = types.ModuleType('google'); col = types.ModuleType('google.colab'); col.drive = mock.MagicMock(); g.colab = col
sys.modules['google'] = g; sys.modules['google.colab'] = col
os.environ['MEGAN_SC_RUNNER'] = '1'                  # the setup block of the notebooks then skips the Drive mount
os.environ.setdefault('MEGAN_SC_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use('Agg')
G = {'__name__': '__main__'}
t_all = time.time()
for i, c in enumerate(cells):
    code = '\n'.join(l for l in c.split('\n') if not l.lstrip().startswith(('%', '!')))
    t0 = time.time()
    try:
        exec(compile(code, f'{os.path.basename(NB)}:cell{i + 1}', 'exec'), G)
    except Exception:
        traceback.print_exc(); print(f'### cell {i + 1} FAILED after {time.time() - t0:.1f} s', flush=True); sys.exit(1)
    print(f'### cell {i + 1} ok {time.time() - t0:.1f} s', flush=True)
    if i == 0:
        for k, v in OV:
            G[k.strip()] = eval(v, G); print(f'### override {k.strip()} = {G[k.strip()]!r}', flush=True)
print(f'### {os.path.basename(NB)}: all {len(cells)} code cells ok in {time.time() - t_all:.0f} s')
