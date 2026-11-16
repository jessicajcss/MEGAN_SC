#!/usr/bin/env python3
"""offline_run.py - run notebook 10 (or the region wizard's engine) without Earth Engine, on synthetic inputs.

    python tools/offline_run.py <notebook.ipynb> <work_dir> "NAME=python_expression;..."

Earth Engine is replaced by a stand-in that returns deterministic synthetic fields: hourly meteorology with a daily
cycle (and a strip of cells without data on the east edge, so that FIX-1 acts), monthly LAI, and growth-form fractions.
The engine, the procedures FIX-1 to FIX-8, the emission factors and every file written are the real ones. The numbers
mean nothing physically; the run checks that the notebook goes through, writes its files, and gives the same results
for the same inputs. Each override NAME=expression is applied after the first code cell, as in tools/run_notebook.py;
NAME@k=expression applies it after code cell k (for a setting defined in a later cell). Set MEGAN_SC_REPO to run the
notebook on another copy of the repository (notebook 10 writes fixes.py and the filled ecotype map into it).
"""
import json, os, sys, time, traceback, types
from unittest import mock
import numpy as np

NB, WORK = sys.argv[1], os.path.abspath(sys.argv[2])
OV = []
for kv in (sys.argv[3].split(';') if len(sys.argv) > 3 else []):
    if kv.strip():
        k, v = kv.split('=', 1); k, _, at = k.strip().partition('@')
        OV.append((k, v, int(at) if at else 1))
os.makedirs(WORK, exist_ok=True)
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ['MEGAN_SC_RUNNER'] = '1'
os.environ.setdefault('MEGAN_SC_REPO', REPO)
sys.path.insert(0, os.path.join(os.environ['MEGAN_SC_REPO'], 'src'))
g = types.ModuleType('google'); col = types.ModuleType('google.colab'); col.drive = mock.MagicMock(); g.colab = col
sys.modules['google'] = g; sys.modules['google.colab'] = col
import matplotlib
matplotlib.use('Agg')

# ---- the Earth Engine stand-in -------------------------------------------------------------------------------------
fake_ee = mock.MagicMock(name='ee')
fake_ee.Date = lambda s, *a, **k: s                       # dates travel as strings
sys.modules['ee'] = fake_ee


def met_fields(date_iso, lon, lat):
    """Synthetic ERA5-Land-like fields at one UTC hour on the grid (south-first)."""
    day, hh = date_iso.split('T')
    h = int(hh[:2]); doy = int(np.datetime64(day, 'D').astype(int) % 365)
    LON, LAT = np.meshgrid(lon, lat)
    solar = h + LON / 15.0                                # local mean solar hour
    t = 293.0 + 6.0 * np.sin(2 * np.pi * (solar - 9) / 24) - 0.4 * (LAT - LAT.mean()) + 2.0 * np.cos(2 * np.pi * doy / 365)
    par = np.clip(450.0 * np.sin(np.pi * (solar - 6) / 12), 0, None)
    f = {'TEMP2': t, 'PAR': par, 'WINDSPD': 2.0 + 0.1 * (LON - LON.min()), 'PRES': 95000.0 + 50 * (LAT - LAT.min()),
         'QV': np.full_like(t, 0.012), 'SOILM': 0.30 + 0.002 * (LON - LON.min())}
    edge = LON > lon[-1] - 3 * (lon[1] - lon[0])          # no reanalysis data on the 3 easternmost columns (FIX-1)
    for k in f:
        f[k] = np.where(edge, np.nan, f[k])
    return f


def growth_forms(domain, year):
    rng = np.random.default_rng(int(year))
    ny, nx = domain['lat'].size, domain['lon'].size
    tree = rng.uniform(0.2, 0.7, (ny, nx)); herb = rng.uniform(0.0, 0.2, (ny, nx)); crop = rng.uniform(0.0, 0.1, (ny, nx))
    shrub = np.zeros((ny, nx)); vegf = tree + herb + crop
    return {'tree': tree, 'shrub': shrub, 'herb': herb, 'crop': crop, 'veg': vegf}


def install(G):
    """Replace, in the namespace of the notebook, everything that would call Earth Engine (run after every cell)."""
    lon, lat = G.get('lon'), G.get('lat')
    if 'era5' in G:
        G['era5'].megan_met_image = lambda ee, region, date: ('met', date)
    if 'daily_tmean_image' in G:
        G['daily_tmean_image'] = lambda day_iso: ('dtemp', day_iso)
    if 'lai_image' in G:
        G['lai_image'] = lambda y, m: ('lai', y, m)
    if 'laimod' in G:
        G['laimod'].monthly_lai = lambda ee, region, y, m: ('lai', y, m)
    if 'mapbiomas' in G:
        G['mapbiomas'].growth_form_fractions_gridded = lambda ee, region, domain, year, **k: growth_forms(domain, year)
        G['mapbiomas'].resolve_collection_id = lambda ee, year, collection=None: collection
        G['mapbiomas'].sample_image_tiled = lambda ee, img, domain, bands, **k: {b: v for b, v in img.items() if b in bands}
    if 'modis_lc' in G:
        G['modis_lc'].growth_form_fractions = lambda ee, region, year, out_scale_m, **k: growth_forms(G['domain'], year)
    if 'gee' in G and isinstance(G['gee'], types.ModuleType):
        G['gee'].initialize = lambda project=None, **k: fake_ee
    if lon is None or lat is None:
        return

    def sample_multi(img, bands):
        kind = img[0]
        if kind == 'met':
            f = met_fields(img[1], lon, lat)
        elif kind == 'dtemp':
            f = {'D_TEMP': np.nanmean([met_fields(f'{img[1]}T{h:02d}:00', lon, lat)['TEMP2'] for h in range(24)], axis=0)}
        elif kind == 'lai':
            LON, LAT = np.meshgrid(lon, lat)
            f = {'LAI': 3.0 + 1.0 * np.sin(2 * np.pi * (img[2] - 1) / 12) + 0.05 * (LAT - LAT.min())}
        else:
            raise ValueError(kind)
        return {b: np.array(f[b], dtype='float64') for b in bands}
    G['sample_multi'] = sample_multi


cells = [''.join(c['source']) for c in json.load(open(NB, encoding='utf-8'))['cells'] if c['cell_type'] == 'code']
G = {'__name__': '__main__'}
os.chdir(WORK)
t_all = time.time()
for i, c in enumerate(cells):
    code = '\n'.join(l for l in c.split('\n') if not l.lstrip().startswith(('%', '!')))
    t0 = time.time()
    try:
        exec(compile(code, f'{os.path.basename(NB)}:cell{i + 1}', 'exec'), G)
    except SystemExit as e:
        print(f'### cell {i + 1} stopped: {e}', flush=True); sys.exit(1)
    except Exception:
        traceback.print_exc(); print(f'### cell {i + 1} FAILED after {time.time() - t0:.1f} s', flush=True); sys.exit(1)
    for k, v, at in OV:
        if at == i + 1:
            G[k] = eval(v, G); print(f'### override after cell {at}: {k} = {G[k]!r}', flush=True)
    install(G)
    print(f'### cell {i + 1} ok {time.time() - t0:.1f} s', flush=True)
print(f'### {os.path.basename(NB)}: all {len(cells)} code cells ok in {time.time() - t_all:.0f} s')
