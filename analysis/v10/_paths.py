"""Paths shared by the analysis scripts of analysis/v10 (repository layout of biogenic_repo).

Every path has a default inside the repository and can be overridden with an environment variable, so the same scripts
run on a laptop, in Google Colab or in the Docker image. `python reproduce.py paths` prints the values in use.

  MEGAN_SC_RUN      production run, one folder per year (<year>/emissions, inputs, ...); 2024 and 2025  data/run_v10
                    with the VIIRS-corrected LAI (notebook 17)
  MEGAN_SC_RUN_MODIS_ARCHIVE  2024 and 2025 computed with the MODIS LAI                 data/run_v10_modis_lai_2024_2025
  MEGAN_SC_RUN_MODIS  the MODIS-LAI record of every year (view built by tools/lai_variant_view.py; read by the
                    step laicorr and by the MODIS control of notebook 19)                    data/run_v10_modis_view
  MEGAN_SC_BND      IBGE boundaries (state, mesoregions with the islands, 2024 municipalities)          data/boundaries
  MEGAN_SC_TROPOMI  monthly TROPOMI HCHO fields on the model grid, cached by notebook 14                 data/tropomi_hcho_v10
  MEGAN_SC_S1       the two 28-day test runs of test S1 (3-hourly and hourly sampling)                 data/s1_test
  MEGAN_SC_LAIVIIRS monthly MODIS and VIIRS LAI on the model grid, 2013-2025, cached by notebook 16      data/lai_viirs_v10
  MEGAN_SC_RESULTS  published results, the reference for `reproduce.py compare`                        results/v10
  MEGAN_SC_OUT      where reproduce.py writes the recomputed results                                     reproduced
  MEGAN_SC_NB13     outputs of notebook 13 read by the later scripts (tables and the annual grid cache)  <OUT>/notebook13_meso2024
                    if it exists, else results/v10/notebook13_meso2024
"""
import os

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))


def _p(env, *default):
    return os.path.abspath(os.environ.get(env) or os.path.join(REPO, *default))


RUN = _p('MEGAN_SC_RUN', 'data', 'run_v10')
RUN_MODIS_ARCHIVE = _p('MEGAN_SC_RUN_MODIS_ARCHIVE', 'data', 'run_v10_modis_lai_2024_2025')
RUN_MODIS = _p('MEGAN_SC_RUN_MODIS', 'data', 'run_v10_modis_view')
BND = _p('MEGAN_SC_BND', 'data', 'boundaries')
TROPOMI = _p('MEGAN_SC_TROPOMI', 'data', 'tropomi_hcho_v10')
S1 = _p('MEGAN_SC_S1', 'data', 's1_test')
LAIVIIRS = _p('MEGAN_SC_LAIVIIRS', 'data', 'lai_viirs_v10')
RESULTS = _p('MEGAN_SC_RESULTS', 'results', 'v10')
OUT = _p('MEGAN_SC_OUT', 'reproduced')
_nb13_out = os.path.join(OUT, 'notebook13_meso2024')
NB13 = os.path.abspath(os.environ.get('MEGAN_SC_NB13') or
                       (_nb13_out if os.path.exists(os.path.join(_nb13_out, 'key_numbers.csv'))
                        else os.path.join(RESULTS, 'notebook13_meso2024')))
NB_SRC = os.path.join(REPO, 'notebooks', 'src')
TAG = 'cf0.3_sza70_fire120'                       # notebook 14 screening: cloud fraction < 0.3, SZA < 70°, 20 km fire mask


FILE_EXT = ('.csv', '.png', '.log', '.npz', '.nc', '.gpkg', '.json', '.txt', '.tif', '.md')


def out(*parts, mkdir=True):
    """Path under MEGAN_SC_OUT; the folder (or, for a file name, its parent folder) is created."""
    p = os.path.join(OUT, *parts)
    if mkdir:
        os.makedirs(os.path.dirname(p) if p.lower().endswith(FILE_EXT) else p, exist_ok=True)
    return p


def notebook_source(name):
    """Source of notebooks/src/<name>, for the scripts that run cells of a notebook outside it. The block that opens every
    notebook (notebooks/src/_setup_repo.py: it finds the repository and, in Colab, mounts Google Drive) is replaced by the
    repository folder of this copy."""
    src = open(os.path.join(NB_SRC, name), encoding='utf-8').read()
    return src.replace('__FILE:notebooks/src/_setup_repo.py__', f'REPO_DIR = {REPO!r}; IN_COLAB = False')


def grid_cache():
    f = os.path.join(NB13, 'annual_grids_cache.npz')
    if not os.path.exists(f):
        raise FileNotFoundError(f'{f} not found: run notebook 13 first (python reproduce.py nb13), or set MEGAN_SC_NB13 to its '
                                'output folder.')
    return f


if __name__ == '__main__':
    for k in ('REPO', 'RUN', 'RUN_MODIS_ARCHIVE', 'RUN_MODIS', 'BND', 'TROPOMI', 'S1', 'RESULTS', 'OUT', 'NB13'):
        v = globals()[k]
        print(f'{k:8s} {v}  {"(exists)" if os.path.exists(v) else "(missing)"}')
