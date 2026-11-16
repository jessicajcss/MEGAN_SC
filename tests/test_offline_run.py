"""Notebooks 10, 11 and 25 for a region other than the Santa Catarina example, without Earth Engine.

tools/offline_run.py replaces Earth Engine by synthetic fields; the engine, the procedures, the emission factors and the
files written are the real ones. The run is made on a copy of the repository, because notebook 10 writes fixes.py and
the filled ecotype map into the repository it runs from. The Santa Catarina configuration is checked against the grid of
the shipped run."""
import os, shutil, subprocess, sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
gpd = pytest.importorskip('geopandas')
yaml = pytest.importorskip('yaml')


def test_santa_catarina_configuration_gives_the_grid_of_the_shipped_run():
    sys.path.insert(0, os.path.join(REPO, 'src'))
    from megan_sc import domain as dom
    from megan_sc.io import geotiff
    from rasterio.features import rasterize
    cfg = yaml.safe_load(open(os.path.join(REPO, 'config', 'santa_catarina_2003_2025.yaml'), encoding='utf-8'))
    d = dom.build_domain(cfg)
    assert (d['lat'].size, d['lon'].size) == (144, 224)
    assert np.isclose(d['lat'][0], -29.4875) and np.isclose(d['lon'][0], -53.8875)
    g = gpd.read_file(os.path.join(REPO, cfg['region']['boundary'])).to_crs('EPSG:4326')
    geom = g.geometry.union_all() if hasattr(g.geometry, 'union_all') else g.geometry.unary_union
    mask = np.flipud(rasterize([(geom, 1)], out_shape=(144, 224), transform=geotiff._transform(d['lon'], d['lat']),
                               fill=0, dtype='uint8')).astype(bool)
    assert int(mask.sum()) == 13970


def _copy_repo(dst):
    for sub in ('src', 'notebooks', 'tools', 'config'):
        shutil.copytree(os.path.join(REPO, sub), os.path.join(dst, sub),
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    os.makedirs(os.path.join(dst, 'data'))
    shutil.copytree(os.path.join(REPO, 'data', 'megan_inputs'), os.path.join(dst, 'data', 'megan_inputs'))
    shutil.copy(os.path.join(REPO, 'data', 'sc_ecotype.nc'), os.path.join(dst, 'data'))


def test_notebooks_10_11_25_offline_for_another_region(tmp_path):
    repo = tmp_path / 'MEGAN_SC'; _copy_repo(str(repo))
    mun = gpd.read_file(os.path.join(REPO, 'data', 'boundaries', 'SC_Municipios_2024.gpkg'))
    sel = mun[mun.NM_MUN.isin(['Joinville', 'Araquari'])][['CD_MUN', 'NM_MUN', 'geometry']]
    sel.to_file(tmp_path / 'region.gpkg', driver='GPKG')
    cfg = {'run': {'name': 'test_region'}, 'period': {'start': '2019-01', 'end': '2019-01', 'sample_hours': 6},
           'domain': {'mode': 'shapefile', 'shapefile': str(tmp_path / 'region.gpkg'), 'resolution_deg': 0.1},
           'region': {'name': 'Test region', 'abbr': 'TR', 'zones': str(tmp_path / 'region.gpkg'), 'zones_id': 'CD_MUN'},
           'datasets': {'meteorology': {'source': 'era5-land'}, 'lai': {'source': 'modis'},
                        'vegetation': {'growth_form_source': 'modis', 'tree_canopy_type': 3, 'lai_mask': True},
                        'soil': {'constant_wilt': 0.12},
                        'emission_factors': {'ef_table': 'data/megan_inputs/EFP240120/EFv240120.csv',
                                             'speciation_dir': 'data/megan_inputs/EFP240120/', 'ecotype': 'data/sc_ecotype.nc'}},
           'earth_engine': {'project': 'offline-test'}, 'output': {'dir': 'outputs_test_region', 'tag': 'test_region'}}
    yaml.safe_dump(cfg, open(tmp_path / 'cfg.yaml', 'w'))
    env = dict(os.environ, MEGAN_SC_REPO=str(repo))
    out = tmp_path / 'out'
    r = subprocess.run([sys.executable, str(repo / 'tools' / 'offline_run.py'), str(repo / 'notebooks' / '10_production_run.ipynb'),
                        str(tmp_path / 'work'), f"RUN_CONFIG@3='{tmp_path / 'cfg.yaml'}';OUT_BASE@3='{out}'"],
                       env=env, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    run = out / 'outputs_test_region' / '2019'
    for f in ('emissions/emis_2019_01.nc', 'stats/stats_2019_01.csv', 'stats/municipality_2019_01.csv',
              'inputs/fixlog_2019_01.npz', 'inputs/growthform_2019.npz', 'provenance_v10.json'):
        assert (run / f).exists(), f
    for nb in ('11_verify_run.ipynb', '25_results_any_region.ipynb'):
        r = subprocess.run([sys.executable, str(repo / 'tools' / 'run_notebook.py'), str(repo / 'notebooks' / nb),
                            f"RUN_CONFIG='{tmp_path / 'cfg.yaml'}';OUT_BASE='{out}'"], env=env, capture_output=True, text=True,
                           timeout=600)
        assert r.returncode == 0, (nb, r.stdout[-3000:] + r.stderr[-3000:])
    assert (out / 'outputs_test_region' / 'results' / 'monthly_totals_by_class.csv').exists()
