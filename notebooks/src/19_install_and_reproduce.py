# %% [markdown]
# # 19 · Install the corrected 2024-2025, update the CO₂ file, reproduce everything twice and prepare the release
#
# **What this notebook does.** It closes the LAI correction of 2024 and 2025 (notebooks 16 and 17) in six checked steps,
# without Earth Engine:
#
# 1. **Inputs.** It copies the data the reproduction reads (production run, TROPOMI fields, the two S1 test runs, the LAI
#    cache of notebook 16) from your Drive to the Colab disk and checks every file against `data/DATA_SHA256SUMS.txt`.
#    A match proves that the files in your Drive are the files of the repository.
# 2. **Checks of the rerun** (`megan_sc_outputs_v10_lai_viirs`, notebook 17): complete years, same `fixes.py` and
#    `laicorr.py` as the repository, identical meteorology, no 340 K signature, twelve fix logs; the corrected LAI
#    recomputed with `laicorr.corrected_run_lai` must equal the published fields (cell 7a of notebook 17, which was not
#    run) and the LAI the rerun actually used.
# 3. **Installation** (`INSTALL`, on by default, only when every check passes; safe to run again): the MODIS-LAI years
#    2024 and 2025 move to `megan_sc_outputs_v10_modis_lai_2024_2025`, and the corrected years are copied into the
#    production folder, exactly as cell 12 of notebook 17.
# 4. **CO₂.** It downloads the official NOAA GML global annual means (`co2_annmean_gl.txt`), records the file
#    creation line, and compares it with the copy shipped in `data/external` (taken from a mirror).
# 5. **Reproduction, three passes** with the notebooks and scripts of this repository (`reproduce.py`), never with older
#    copies in Drive:
#    * **A, MODIS control**: run with the MODIS LAI and the shipped CO₂ file. Every CSV must equal the published copy
#      in `results/v10`; this proves that the Colab environment reproduces the published results.
#    * **A2, CO₂ file only**: step `final` with the MODIS LAI and the official CO₂ file, to isolate the effect of the
#      new file on the CO₂-inhibition numbers.
#    * **B, corrected**: everything with the corrected 2024-2025 and the official CO₂ file (notebooks 11, 13 and 14
#      included). Differences from `results/v10` are expected wherever 2024-2025 or CO₂ enter.
#
#    The step `laicorr` always reads the MODIS record of every year (`tools/lai_variant_view.py`), because it documents
#    the MODIS decrease and builds the correction.
# 6. **Summary and release.** A side-by-side table of every number of notebook 13, of the additional diagnostics and of
#    the trends with the MODIS and with the corrected LAI (`lai_sensitivity_summary.csv` and `.md`). With
#    `PROMOTE = True` the corrected results become `results/v10` of your repository copy (the MODIS set is kept as
#    `results/v10_modis_lai_2024_2025`), the official CO₂ file replaces the mirror copy, and the checksums and the
#    manifest are rewritten. With `PACK = True` the five release archives are written (the archive of the MODIS-LAI
#    2024-2025 is the new fifth one). Creating the GitHub release and the Zenodo deposit is done by hand (last cell).
#
# **Before running.** Notebook 17 finished (both years). Put this version of the repository, unzipped, in Drive
# (`REPO_DRIVE`); it must contain `tools/lai_variant_view.py`. Check the paths of cell 1. The whole notebook takes about
# 30-60 min in Colab (most of it copying files from Drive); the three passes take about 5 min each.
#
# **Running it again** in the same Colab runtime (for example with `PROMOTE` and `PACK` switched on after a first run)
is safe: cell 2 restores the MODIS-LAI 2024-2025 on the Colab disk, pass A always uses the shipped CO₂ file, and cell 7
replaces the local archive of an earlier run. After a promotion, the repository copy is final and the notebook stops in
cell 2 if run on it again.

**What it changes in Drive.** Cell 4 moves and copies year folders (only with `INSTALL = True`); cell 9 writes into
# `REPO_DRIVE` and `RELEASE_DIR` (only with `PROMOTE` / `PACK`). Everything else is written to the Colab disk (`WORK`)
# and to `OUT_DRIVE`.

# %%
# 1 · Settings ---------------------------------------------------------------------------------------------------------
import os, sys, json, glob, shutil, hashlib, subprocess, datetime, time, re, urllib.request
try:
    from google.colab import drive; drive.mount('/content/drive'); IN_COLAB = True
except Exception:
    IN_COLAB = False
REPO_DRIVE    = '/content/drive/MyDrive/MEGAN_SC'                                # <-- this repository, unzipped
ROOT_V10      = '/content/drive/MyDrive/megan_sc_outputs_v10'                    # <-- production run (notebook 10)
LAIV_ROOT     = '/content/drive/MyDrive/megan_sc_outputs_v10_lai_viirs'          # <-- rerun of notebook 17
ARCHIVE       = '/content/drive/MyDrive/megan_sc_outputs_v10_modis_lai_2024_2025'  # MODIS-LAI 2024-2025 after install
VIIRS_DIR     = ROOT_V10 + '/analysis_v10/lai_viirs_check'                       # <-- outputs of notebook 16 (cache)
TROPOMI_CACHE = ROOT_V10 + '/analysis_v10/tropomi_hcho/cache_cf0.3_sza70_fire120'  # <-- monthly HCHO of notebook 14
DATA_ZIPS     = '/content/drive/MyDrive/biogenic_repo_release_assets'   # <-- the three release zips of October 2026
                                                  #     (megan_sc_run_v10_2003-2025.zip, tropomi_hcho_v10_2019-2025.zip,
                                                  #      megan_sc_s1_test_2019.zip)
S1_PARENT     = '/content/drive/MyDrive'          # only without the zips: holds megan_sc_outputs_v10/2019/*_test7d* and megan_sc_outputs_v10_hourly/2019
WORK          = '/content/megan_sc_r21'           # Colab disk: repository copy, data copy, reproduction outputs
OUT_DRIVE     = ROOT_V10 + '/analysis_v10/notebook19'   # summaries and logs kept in Drive
INSTALL       = True                              # install the corrected years (only if every check of cell 3 passes)
PROMOTE       = False                             # write the corrected results, CO2 file, sums and manifest into REPO_DRIVE
PACK          = False                             # with PROMOTE: write the release archives into RELEASE_DIR
RELEASE_DIR   = '/content/drive/MyDrive/MEGAN_SC_release_v1.0.0'
NOAA_URL      = 'https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_gl.txt'
YEARS_CORR    = (2024, 2025)
FIXES_SHA16   = 'a933b29602929a82'                # fixes.py of every year of the run
REPO = os.path.join(WORK, 'MEGAN_SC'); DATA = os.path.join(REPO, 'data')
os.makedirs(WORK, exist_ok=True); os.makedirs(OUT_DRIVE, exist_ok=True)
LARGE = {'run_v10': None, 'tropomi_hcho_v10': TROPOMI_CACHE, 's1_test': S1_PARENT, 'lai_viirs_v10': VIIRS_DIR}
def sha256(p, n=None):
    d = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            d.update(b)
    return d.hexdigest()[:n] if n else d.hexdigest()
def variant(root, y):
    p = os.path.join(root, str(y), 'provenance_v10.json')
    return json.load(open(p)).get('variant') if os.path.exists(p) else None
print('Colab' if IN_COLAB else 'local run', '| repository:', REPO_DRIVE)
for k in ('REPO_DRIVE', 'ROOT_V10', 'LAIV_ROOT', 'VIIRS_DIR', 'TROPOMI_CACHE'):
    print(f'  {k:13s} {globals()[k]}  {"ok" if os.path.exists(globals()[k]) else "NOT FOUND"}')

# %% [markdown]
# ## Cell 2 · Local copy of the repository and of the data
#
# The repository is copied to the Colab disk without its large data folders. The data come from the three release zips
# of October 2026 in `DATA_ZIPS` (each zip is checked against `data/RELEASE_ASSETS_SHA256SUMS.txt` before it is
# unpacked) and, for `lai_viirs_v10/`, from the outputs of notebook 16 in Drive. A folder whose zip is not in `DATA_ZIPS`
# is copied file by file from the production folders in Drive (for `run_v10/2024` and `/2025` the MODIS-LAI years:
# from `ARCHIVE` once the corrected years are installed, from `ROOT_V10` before). Then every file listed in
# `data/DATA_SHA256SUMS.txt` under the four folders is checked. A different checksum stops the notebook. The S1 test
# runs share their folder with the production run of 2019 in Drive, so their `provenance_v10.json` there is the
# production one: use the zip for `s1_test`.

# %%
# 2 · Local copy of the repository and of the data ----------------------------------------------------------------------
%pip -q install -r {REPO_DRIVE}/requirements-reproduce.txt
assert os.path.exists(os.path.join(REPO_DRIVE, 'tools', 'lai_variant_view.py')), \
    'REPO_DRIVE is an older version of the repository (tools/lai_variant_view.py missing): copy the current one'
if not os.path.exists(REPO):
    shutil.copytree(REPO_DRIVE, REPO, ignore=shutil.ignore_patterns(*LARGE, 'run_v10_modis_*', '_download', 'reproduced',
                                                                     '__pycache__', '.ipynb_checkpoints'))
sys.path.insert(0, os.path.join(REPO, 'src')); sys.path.insert(0, os.path.join(REPO, 'tools'))
SUMS = {}
for line in open(os.path.join(DATA, 'DATA_SHA256SUMS.txt'), encoding='utf-8'):
    h, rel = line.rstrip('\n').split('  ', 1); SUMS[rel] = h
assert not any(r.startswith('run_v10_modis_lai_2024_2025/') for r in SUMS), \
    ('this copy of the repository was already promoted by cell 9 (its checksum list holds the corrected run_v10): the '
     'reproduction is now `python reproduce.py all` in the repository; to repeat notebook 19, start from the repository '
     'copy that was in Drive before the promotion')
def modis_year_source(y):
    """Drive folder holding the MODIS-LAI run of year y."""
    if variant(ARCHIVE, y) == 'base':
        return os.path.join(ARCHIVE, str(y))
    assert variant(ROOT_V10, y) == 'base', f'{y}: MODIS-LAI year found neither in {ARCHIVE} nor in {ROOT_V10}'
    return os.path.join(ROOT_V10, str(y))
def drive_source(rel):
    top, rest = rel.split('/', 1)
    if top == 'run_v10':
        y = rest.split('/', 1)[0]
        if y.isdigit() and int(y) in YEARS_CORR:
            return os.path.join(modis_year_source(int(y)), rest.split('/', 1)[1])
        return os.path.join(ROOT_V10, rest)
    return os.path.join(LARGE[top], rest)
import zipfile
ZIPS = {'run_v10': 'megan_sc_run_v10_2003-2025.zip', 'tropomi_hcho_v10': 'tropomi_hcho_v10_2019-2025.zip',
        's1_test': 'megan_sc_s1_test_2019.zip'}
ZSUM = dict(reversed(l.split()) for l in open(os.path.join(DATA, 'RELEASE_ASSETS_SHA256SUMS.txt')) if l.strip())
def folder_ok(top):
    return all(os.path.exists(os.path.join(DATA, r)) for r in SUMS if r.split('/', 1)[0] == top)
print('release zips in', DATA_ZIPS, ':', sorted(f for f in os.listdir(DATA_ZIPS) if f.endswith('.zip')) if os.path.isdir(DATA_ZIPS) else 'FOLDER NOT FOUND')
t0 = time.time()
for top, zname in ZIPS.items():
    zp = os.path.join(DATA_ZIPS, zname)
    if folder_ok(top):
        print(f'{top}: already on the Colab disk (checked below)'); continue
    if not os.path.exists(zp):
        if top == 's1_test':                       # its files cannot be taken from Drive (see the text above this cell)
            raise FileNotFoundError(
                f'{zname} not found in {DATA_ZIPS}. The two S1 test runs are not in the production folders of Drive '
                f'(they share the 2019 folder with the production run). Upload {zname} (39 MB; on your computer in '
                f'the folder of the release zips) to {DATA_ZIPS} and run this cell again.')
        print(f'{top}: zip not found, copied file by file from Drive'); continue
    zl = os.path.join(WORK, zname); shutil.copyfile(zp, zl)
    assert sha256(zl) == ZSUM[zname], f'{zname}: checksum differs from data/RELEASE_ASSETS_SHA256SUMS.txt'
    with zipfile.ZipFile(zl) as zf:
        zf.extractall(DATA)
    os.remove(zl); print(f'{top}: unpacked from {zname} (checksum ok) after {time.time() - t0:.0f} s')
bad, n_copied = [], 0
for rel, h in SUMS.items():
    top = rel.split('/', 1)[0]
    if top not in LARGE:
        continue
    dst = os.path.join(DATA, rel)
    if not (os.path.exists(dst) and sha256(dst) == h):
        src = drive_source(rel)
        if not os.path.exists(src):
            bad.append((rel, 'missing in Drive: ' + src)); continue
        os.makedirs(os.path.dirname(dst), exist_ok=True); shutil.copyfile(src, dst); n_copied += 1
        if sha256(dst) != h:
            bad.append((rel, 'checksum differs: ' + src))
print(f'{sum(r.split("/", 1)[0] in LARGE for r in SUMS)} data files checked, {n_copied} copied in {time.time() - t0:.0f} s, '
      f'{len(bad)} problems')
for b in bad[:30]:
    print('  ', *b)
assert not bad, 'the data in Drive differ from data/DATA_SHA256SUMS.txt: see the list above'
print('the production run (MODIS LAI), TROPOMI fields, S1 runs and LAI cache in Drive equal the repository data')

# %% [markdown]
# ## Cell 3 · Checks of the rerun of notebook 17
#
# | Check | Must be |
# |---|---|
# | twelve monthly files and the annual file of each year | present |
# | `provenance_v10.json`: variant, `fixes.py` hash, `laicorr.py` hash | `lai_viirs`, a933b29602929a82, the hash of `src/megan_sc/laicorr.py` of this repository |
# | monthly means of the six meteorological inputs, state cells | identical to the MODIS-LAI run (difference 0) |
# | cells with trees, BVOC > 0 and isoprene = 0 (the 340 K signature) | 0 |
# | fix logs per year | 12 |
# | corrected LAI recomputed now vs the published fields (`results/v10/lai_viirs_check/lai_viirs_corrected_2024_2025.npz`) | max abs difference < 1e-4 (cell 7a of notebook 17) |
# | corrected LAI vs the LAI the rerun used (`inputs/lai_<year>_<mm>.tif`, vegetated cells) | max abs difference < 1e-4 (float32 storage) |

# %%
# 3 · Checks of the rerun ------------------------------------------------------------------------------------------------
import numpy as np, pandas as pd, rasterio, geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import from_origin
from netCDF4 import Dataset
from megan_sc import laicorr
LAICORR_SHA16 = sha256(os.path.join(REPO, 'src', 'megan_sc', 'laicorr.py'), 16)
RUN_M = os.path.join(DATA, 'run_v10')                     # local copy, MODIS LAI in every year
with Dataset(os.path.join(RUN_M, '2019', 'emissions', 'annual_emis_2019.nc')) as ds:
    LAT = np.asarray(ds.variables['lat'][:], 'float64'); LON = np.asarray(ds.variables['lon'][:], 'float64')
RES = float(LON[1] - LON[0])
SC = gpd.read_file(os.path.join(DATA, 'boundaries', 'sc_state_ibge.gpkg')).to_crs('EPSG:4326')
SC_POLY = SC.geometry.union_all() if hasattr(SC.geometry, 'union_all') else SC.geometry.unary_union
SC_MASK = np.flipud(rasterize([(SC_POLY, 1)], out_shape=(LAT.size, LON.size),
                              transform=from_origin(LON[0] - RES / 2, LAT[-1] + RES / 2, RES, RES), fill=0, dtype='uint8')).astype(bool)
print('state cells:', int(SC_MASK.sum()), '(13,970 expected)'); assert int(SC_MASK.sum()) == 13970
def bands(p):                                               # all bands, south first
    with rasterio.open(p) as s:
        return np.flip(s.read().astype('float64'), axis=1)
CLS = ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM', 'NO', 'SQT_HR', 'SQT_LR', 'MEOH', 'ACTO',
       'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO']
CHECKS = []
def check(name, ok, detail=''):
    CHECKS.append(dict(check=name, ok=bool(ok), detail=str(detail))); print(('ok    ' if ok else 'FAILED'), name, detail)
for y in YEARS_CORR:
    d = os.path.join(LAIV_ROOT, str(y))
    nc = [os.path.join(d, 'emissions', f'emis_{y}_{m:02d}.nc') for m in range(1, 13)] + [os.path.join(d, 'emissions', f'annual_emis_{y}.nc')]
    check(f'{y}: 12 monthly + annual files', all(map(os.path.exists, nc)), sum(map(os.path.exists, nc)))
    P = json.load(open(os.path.join(d, 'provenance_v10.json')))
    lc = P.get('lai_correction', {})
    check(f'{y}: provenance', P.get('variant') == 'lai_viirs' and P.get('fixes_sha256') == FIXES_SHA16 and
          lc.get('laicorr_sha256') == LAICORR_SHA16,
          f"variant {P.get('variant')}, fixes {P.get('fixes_sha256')}, laicorr {lc.get('laicorr_sha256')} (repository {LAICORR_SHA16})")
    mx = 0.0
    for m in range(1, 13):
        a = bands(os.path.join(RUN_M, str(y), 'inputs', f'met_{y}_{m:02d}.tif')); b = bands(os.path.join(d, 'inputs', f'met_{y}_{m:02d}.tif'))
        ok_ = SC_MASK[None] & np.isfinite(a) & np.isfinite(b); mx = max(mx, float(np.abs(a[ok_] - b[ok_]).max()))
    check(f'{y}: meteorology identical to the MODIS-LAI run', mx == 0.0, f'max abs difference {mx:.3g}')
    with Dataset(nc[-1]) as ds:
        E = {c: np.nan_to_num(np.asarray(ds.variables[c][:], 'float64')) for c in CLS}
    bv = sum(E[c] for c in CLS if c not in ('NO', 'CO'))
    gf = dict(np.load(os.path.join(d, 'inputs', f'growthform_{min(y, 2024)}.npz')))
    sig = SC_MASK & (gf['tree'] > 0.1) & (bv > 0) & (E['ISOP'] == 0)
    check(f'{y}: 340 K signature cells', int(sig.sum()) == 0, int(sig.sum()))
    nlog = len(glob.glob(os.path.join(d, 'inputs', f'fixlog_{y}_*.npz')))
    check(f'{y}: fix logs', nlog == 12, nlog)
LAI_CORR, LAI_LOG = laicorr.corrected_run_lai(RUN_M, os.path.join(DATA, 'lai_viirs_v10', 'lai_viirs_modis_cache.npz'), SC_MASK, years=YEARS_CORR)
PUB = np.load(os.path.join(REPO, 'results', 'v10', 'lai_viirs_check', 'lai_viirs_corrected_2024_2025.npz'))
d_pub = max(float(np.abs(PUB[f'LAI_{y}_{m:02d}'] - LAI_CORR[(y, m)]).max()) for (y, m) in LAI_CORR)
check('corrected LAI recomputed vs published fields (cell 7a of notebook 17)', d_pub < 1e-4, f'max abs difference {d_pub:.2e}')
d_use = 0.0
for (y, m), a in LAI_CORR.items():
    used = np.flipud(rasterio.open(os.path.join(LAIV_ROOT, str(y), 'inputs', f'lai_{y}_{m:02d}.tif')).read(1).astype('float64'))
    gf = dict(np.load(os.path.join(LAIV_ROOT, str(y), 'inputs', f'growthform_{min(y, 2024)}.npz')))
    veg = sum(gf[k] for k in ('tree', 'shrub', 'herb', 'crop')) > 0
    d_use = max(d_use, float(np.abs(used[veg] - a[veg]).max()))
check('corrected LAI vs the LAI used by the rerun (vegetated cells)', d_use < 1e-4, f'max abs difference {d_use:.2e}')
rows = []
for y in YEARS_CORR:
    t0_ = pd.read_csv(os.path.join(RUN_M, str(y), 'stats', f'annual_by_family_{y}.csv')).set_index('species')['annual_tonnes_SC']
    t1_ = pd.read_csv(os.path.join(LAIV_ROOT, str(y), 'stats', f'annual_by_family_{y}.csv')).set_index('species')['annual_tonnes_SC']
    for g, cl in {'BVOC': [c for c in CLS if c not in ('NO', 'CO')], 'ISOP': ['ISOP'], 'NO': ['NO'], 'CO': ['CO']}.items():
        rows.append(dict(year=y, group=g, modis_t=round(t0_[cl].sum()), corrected_t=round(t1_[cl].sum()),
                         change_pct=round(100 * (t1_[cl].sum() / t0_[cl].sum() - 1), 2)))
print(pd.DataFrame(rows).to_string(index=False))
CHECKS_OK = all(c['ok'] for c in CHECKS)
pd.DataFrame(CHECKS).to_csv(os.path.join(OUT_DRIVE, 'rerun_checks.csv'), index=False)
print('\nALL CHECKS PASSED' if CHECKS_OK else '\nSOME CHECKS FAILED: the corrected years will not be installed')

# %% [markdown]
# ## Cell 4 · Installation of the corrected years (as cell 12 of notebook 17)
#
# Safe to run again: a year already installed (variant `lai_viirs` in `ROOT_V10`, variant `base` in `ARCHIVE`) is only
# checked (same files and checksums of the NetCDF outputs as the rerun). Any other state stops the cell.

# %%
# 4 · Installation ---------------------------------------------------------------------------------------------------------
def same_tree(a, b, pattern='emissions/*.nc'):
    fa = sorted(os.path.relpath(p, a) for p in glob.glob(os.path.join(a, pattern)))
    fb = sorted(os.path.relpath(p, b) for p in glob.glob(os.path.join(b, pattern)))
    return fa == fb and all(sha256(os.path.join(a, f)) == sha256(os.path.join(b, f)) for f in fa)
INSTALLED = {}
for y in YEARS_CORR:
    vr, va = variant(ROOT_V10, y), variant(ARCHIVE, y)
    src_new, dst, arc = os.path.join(LAIV_ROOT, str(y)), os.path.join(ROOT_V10, str(y)), os.path.join(ARCHIVE, str(y))
    if vr == 'lai_viirs' and va == 'base':
        assert same_tree(src_new, dst), f'{dst} is a corrected year but differs from {src_new}'
        INSTALLED[y] = 'already installed (checked)'
    elif vr == 'base' and va is None:
        if INSTALL and CHECKS_OK:
            os.makedirs(ARCHIVE, exist_ok=True)
            shutil.move(dst, arc); shutil.copytree(src_new, dst)
            assert variant(ROOT_V10, y) == 'lai_viirs' and variant(ARCHIVE, y) == 'base' and same_tree(src_new, dst)
            INSTALLED[y] = 'installed now'
        else:
            INSTALLED[y] = 'not installed (' + ('INSTALL = False' if not INSTALL else 'checks failed') + ')'
    else:
        raise RuntimeError(f'{y}: unexpected state, ROOT_V10 variant {vr!r}, ARCHIVE variant {va!r}; fix by hand before continuing')
print(INSTALLED)
CORRECTED_IN_DRIVE = all(variant(ROOT_V10, y) == 'lai_viirs' for y in YEARS_CORR)
print('production folder', ROOT_V10, 'now holds the corrected 2024-2025' if CORRECTED_IN_DRIVE else 'still holds the MODIS-LAI 2024-2025')

# %% [markdown]
# ## Cell 5 · Official NOAA CO₂ file
#
# The file is downloaded from NOAA GML (`NOAA_URL`). Its header carries the creation date of the file, recorded here,
# because NOAA revises the most recent years. The values are compared with the copy shipped in `data/external`,
# which came from a mirror (header of that file). The γCO₂ factor of MEGANv3.21 is computed for both, so the size of
# the change is visible before pass A2.

# %%
# 5 · Official NOAA CO2 file ------------------------------------------------------------------------------------------------
CO2_NEW = os.path.join(WORK, 'co2_annmean_gl_noaa.txt'); CO2_OLD = os.path.join(DATA, 'external', 'co2_annmean_gl.txt')
CO2_OLD_KEEP = os.path.join(WORK, 'co2_annmean_gl_shipped.txt')
if not os.path.exists(CO2_OLD_KEEP):                  # first run in this runtime: the copy shipped with the repository
    shutil.copyfile(CO2_OLD, CO2_OLD_KEEP)
assert 'datasets/co2-ppm' in open(CO2_OLD_KEEP, encoding='utf-8', errors='replace').read(), \
    f'{CO2_OLD_KEEP} is not the CO2 file shipped with the repository (mirror copy of 13 October 2026)'
req = urllib.request.Request(NOAA_URL, headers={'User-Agent': 'MEGAN_SC notebook 19'})
raw = urllib.request.urlopen(req, timeout=60).read()
open(CO2_NEW, 'wb').write(raw)
RETRIEVED = datetime.datetime.utcnow().isoformat(timespec='seconds') + 'Z'
txt = raw.decode('utf-8', errors='replace')
created = [l for l in txt.splitlines() if re.search(r'creat', l, re.I) and l.startswith('#')]
def read_co2(path):
    rows = [re.match(r'^\s*(\d{4})\s+([0-9.]+)\s+([0-9.]+)', l) for l in open(path, encoding='utf-8', errors='replace') if not l.startswith('#')]
    return pd.DataFrame([(int(m[1]), float(m[2]), float(m[3])) for m in rows if m], columns=['year', 'co2', 'unc']).set_index('year')
NEW, OLD = read_co2(CO2_NEW), read_co2(CO2_OLD_KEEP)
assert set(range(2003, 2026)) <= set(NEW.index), f'the NOAA file lacks years of 2003-2025: {sorted(set(range(2003, 2026)) - set(NEW.index))}'
gco2 = lambda c: 1.344 - 1.344 * (0.7 * c) ** 1.4614 / (585 ** 1.4614 + (0.7 * c) ** 1.4614)   # megvea.f, Ci = 0.7 CO2
Y = list(range(2003, 2026))
CMPC = pd.DataFrame({'noaa_ppm': NEW.co2.reindex(Y), 'shipped_ppm': OLD.co2.reindex(Y)})
CMPC['diff_ppm'] = CMPC.noaa_ppm - CMPC.shipped_ppm
CMPC['gamma_noaa'] = gco2(CMPC.noaa_ppm); CMPC['gamma_shipped'] = gco2(CMPC.shipped_ppm)
CMPC['gamma_diff'] = CMPC.gamma_noaa - CMPC.gamma_shipped
CMPC.to_csv(os.path.join(OUT_DRIVE, 'co2_noaa_vs_shipped.csv'))
print('NOAA file:', NOAA_URL, '| retrieved', RETRIEVED, '| sha256', sha256(CO2_NEW, 16))
print('creation line(s):', created if created else 'none found; header follows'); print('\n'.join(txt.splitlines()[:6]) if not created else '')
print(f'2003-2025: max |NOAA - shipped| = {CMPC.diff_ppm.abs().max():.3f} ppm, max |gamma difference| = {CMPC.gamma_diff.abs().max():.2e}')
print(CMPC.round(4).to_string())

# %% [markdown]
# ## Cell 6 · Pass A: MODIS control (must give back `results/v10`)

# %%
# 6 · Pass A: MODIS control -------------------------------------------------------------------------------------------------
def reproduce(out, steps=('all',), **env_extra):
    env = dict(os.environ, MEGAN_SC_OUT=out, MPLBACKEND='Agg', **env_extra)
    env.pop('MEGAN_SC_NB13', None)
    t0 = time.time()
    r = subprocess.run([sys.executable, 'reproduce.py', *steps], cwd=REPO, env=env, capture_output=True, text=True)
    print(r.stdout[-4000:], r.stderr[-2000:]); print(f'exit {r.returncode} in {time.time() - t0:.0f} s')
    return r.returncode
def compare_table(out):
    log = os.path.join(out, 'logs', 'compare.log'); rows = []
    for l in open(log, encoding='utf-8', errors='replace'):
        m = re.match(r'^(\S+)\s{2,}(IDENTICAL|EQUAL|DIFFERENT|NEW|PIXEL-IDENTICAL|PNG DIFF|PNG SIZE|PNG\?|PRESENT)\s*(.*)$', l.rstrip('\n'))
        if m:
            rows.append(dict(file=m[1], status=m[2], info=m[3].strip()))
    return pd.DataFrame(rows)
assert all(variant(RUN_M, y) == 'base' for y in YEARS_CORR), 'the local run_v10 must hold the MODIS-LAI years for pass A'
shutil.copyfile(CO2_OLD_KEEP, CO2_OLD)          # pass A uses the shipped CO2 file (an earlier run of cell 7 replaced it)
RES_A = os.path.join(REPO, 'results', 'v10_modis_lai_2024_2025')
RES_A = RES_A if os.path.isdir(RES_A) else os.path.join(REPO, 'results', 'v10')    # the published MODIS-LAI results
OUT_A = os.path.join(WORK, 'reproduced_modis')
rcA = reproduce(OUT_A, MEGAN_SC_RESULTS=RES_A)
CA = compare_table(OUT_A); CA.to_csv(os.path.join(OUT_DRIVE, 'passA_compare.csv'), index=False)
DIFF_A = CA[(CA.status == 'DIFFERENT') & CA.file.str.endswith('.csv')]
print(CA.status.value_counts().to_string())
print('PASS A REPRODUCES THE PUBLISHED RESULTS' if rcA == 0 and DIFF_A.empty else 'PASS A DIFFERS FROM results/v10:\n' + DIFF_A.to_string(index=False))

# %% [markdown]
# ## Cell 7 · Pass A2 (official CO₂ file, MODIS LAI, step `final`) and pass B (corrected LAI, official CO₂ file)
#
# Before pass B the local `data/run_v10/2024` and `/2025` move to `data/run_v10_modis_lai_2024_2025` and the corrected
# years take their place: the same layout the release will have.

# %%
# 7 · Passes A2 and B ------------------------------------------------------------------------------------------------------
shutil.copyfile(CO2_NEW, CO2_OLD)                                     # official NOAA file from here on
OUT_A2 = os.path.join(WORK, 'reproduced_modis_co2noaa')
rcA2 = reproduce(OUT_A2, ('final',))
ARCH_L = os.path.join(DATA, 'run_v10_modis_lai_2024_2025')
for y in YEARS_CORR:
    if variant(RUN_M, y) == 'base':                 # cell 2 restores the MODIS-LAI year here on every run
        a_y = os.path.join(ARCH_L, str(y))
        if os.path.isdir(a_y):
            shutil.rmtree(a_y)                      # left by an earlier run of this cell; replaced, never nested
        os.makedirs(ARCH_L, exist_ok=True); shutil.move(os.path.join(RUN_M, str(y)), a_y)
        shutil.copytree(os.path.join(LAIV_ROOT, str(y)), os.path.join(RUN_M, str(y)),
                        ignore=shutil.ignore_patterns('checkpoints', 'maps', '*.tif.aux.xml'))
    assert variant(RUN_M, y) == 'lai_viirs' and variant(ARCH_L, y) == 'base'
# keep exactly the files the release lists for every year (the rerun folders also hold GeoTIFF stacks of the emissions)
listed = {r for r in SUMS if r.startswith('run_v10/')}
for y in YEARS_CORR:
    for p in glob.glob(os.path.join(RUN_M, str(y), '**', '*'), recursive=True):
        rel = 'run_v10/' + os.path.relpath(p, RUN_M).replace(os.sep, '/')
        if os.path.isfile(p) and rel not in listed:
            os.remove(p)
    missing = [r for r in listed if r.startswith(f'run_v10/{y}/') and not os.path.exists(os.path.join(DATA, r))]
    assert not missing, f'{y}: files of the release layout missing in the rerun: {missing[:5]}'
OUT_B = os.path.join(WORK, 'reproduced_corrected')
STEPS = [l.split()[0] for l in subprocess.run([sys.executable, 'reproduce.py', 'list'], cwd=REPO, capture_output=True,
                                              text=True).stdout.splitlines() if l[:1].strip()]
assert STEPS and STEPS[-1] == 'compare', STEPS
rcB = reproduce(OUT_B, tuple(s for s in STEPS if s != 'compare'))   # every analysis step: must end with exit 0
# the comparison with the MODIS results runs on its own: pass B differs from results/v10 wherever 2024-2025 or LAI enter,
# so compare_results.py ends with exit 1 here by design; cell 9 only requires that the comparison finished
rcBcmp = reproduce(OUT_B, ('compare',))
CMP_B_DONE = 'CSV files:' in open(os.path.join(OUT_B, 'logs', 'compare.log'), encoding='utf-8', errors='replace').read()
print(f'pass B: analysis steps exit {rcB}; comparison with the MODIS results exit {rcBcmp} '
      f'({"finished" if CMP_B_DONE else "DID NOT FINISH"}; exit 1 only means that some files differ, as expected)')
CB = compare_table(OUT_B); CB.to_csv(os.path.join(OUT_DRIVE, 'passB_compare_vs_published_modis.csv'), index=False)
print(CB.status.value_counts().to_string())
print('CSV files that differ from the MODIS results (expected wherever 2024-2025 or CO2 enter):')
print(CB[(CB.status == 'DIFFERENT') & CB.file.str.endswith('.csv')].to_string(index=False))

# %% [markdown]
# ## Cell 8 · Side-by-side summary: MODIS LAI and corrected LAI
#
# Every number of `manuscript_numbers.csv` (notebook 13), `final_review_numbers.csv` and `ijest_review_numbers.csv`, the
# trends of `trends.csv` and `trends_fdr_hr.csv`, and the annual totals of 2024 and 2025, from pass A (MODIS LAI) and
# pass B (corrected LAI). The CO₂-inhibition lines are also given from pass A2, which isolates the new CO₂ file.

# %%
# 8 · Side-by-side summary ----------------------------------------------------------------------------------------------------
def kv(out, rel, key='key', val='value'):
    p = os.path.join(out, rel)
    if not os.path.exists(p):
        return pd.Series(dtype=str)
    t = pd.read_csv(p, dtype=str, keep_default_na=False).set_index(key)
    x = t[val].where(t[val] != '', t['en']) if 'en' in t.columns else t[val]     # text-only numbers of notebook 13
    return x[~x.index.duplicated()]
parts = []
for rel, val in (('notebook13_meso2024/manuscript_numbers.csv', 'value'), ('final_review/final_review_numbers.csv', 'value'),
                 ('ijest_review/ijest_review_numbers.csv', 'value')):
    a, b = kv(OUT_A, rel, val=val), kv(OUT_B, rel, val=val)
    df = pd.DataFrame({'modis_lai': a, 'corrected_lai': b}); df.index.name = 'key'; df['source'] = rel
    parts.append(df.reset_index())
A2 = kv(OUT_A2, 'final_review/final_review_numbers.csv')
S = pd.concat(parts, ignore_index=True)
S['modis_lai_noaa_co2'] = [A2.get(k, '') if (r == 'final_review/final_review_numbers.csv' and 'co2' in k.lower()) else ''
                           for k, r in zip(S.key, S.source)]
S['changed'] = S.modis_lai != S.corrected_lai
TR = []
for rel, cols in (('notebook13_meso2024/trends.csv', ['pct_yr', 'pct_lo95', 'pct_hi95', 'MK_p', 'HR_p']),
                  ('final_review/trends_fdr_hr.csv', ['MK_p', 'HR_p', 'BH_q', 'p_from_2005', 'pct_yr_2003_2024', 'p_2003_2024'])):
    a = pd.read_csv(os.path.join(OUT_A, rel)).set_index(pd.read_csv(os.path.join(OUT_A, rel)).columns[0])
    b = pd.read_csv(os.path.join(OUT_B, rel)).set_index(pd.read_csv(os.path.join(OUT_B, rel)).columns[0])
    for g in a.index:
        for c in cols:
            if c in a.columns and g in b.index:
                TR.append(dict(source=rel, series=g, stat=c, modis_lai=a.loc[g, c], corrected_lai=b.loc[g, c]))
TR = pd.DataFrame(TR)
T2a = pd.read_csv(os.path.join(OUT_A, 'notebook13_meso2024', 'table2_annual_totals_SC_groups.csv')).set_index('year')
T2b = pd.read_csv(os.path.join(OUT_B, 'notebook13_meso2024', 'table2_annual_totals_SC_groups.csv')).set_index('year')
T2 = pd.concat({'modis_lai': T2a.loc[list(YEARS_CORR)], 'corrected_lai': T2b.loc[list(YEARS_CORR)]}, axis=0)
S.to_csv(os.path.join(OUT_DRIVE, 'lai_sensitivity_summary.csv'), index=False)
TR.to_csv(os.path.join(OUT_DRIVE, 'lai_sensitivity_trends.csv'), index=False)
T2.to_csv(os.path.join(OUT_DRIVE, 'lai_sensitivity_table2_2024_2025.csv'))
pv = TR[(TR.source.str.startswith('notebook13')) & TR.stat.isin(['pct_yr', 'pct_lo95', 'pct_hi95', 'MK_p'])].pivot(
    index='series', columns='stat', values=['modis_lai', 'corrected_lai'])
md = ['# Trends 2003-2025 with the MODIS LAI and with the corrected LAI of 2024-2025 (notebook 19)', '',
      f'Pass A: MODIS LAI, shipped CO2 file. Pass B: corrected LAI, NOAA CO2 file retrieved {RETRIEVED}.', '',
      '| series | MODIS %/yr (95 % CI) | MK p | corrected %/yr (95 % CI) | MK p |', '|---|---|---|---|---|']
for g in pv.index:
    f = lambda k, s: float(pv.loc[g, (k, s)])
    md.append(f"| {g} | {f('modis_lai', 'pct_yr'):+.2f} ({f('modis_lai', 'pct_lo95'):+.2f} to {f('modis_lai', 'pct_hi95'):+.2f}) | "
              f"{f('modis_lai', 'MK_p'):.4f} | {f('corrected_lai', 'pct_yr'):+.2f} ({f('corrected_lai', 'pct_lo95'):+.2f} to "
              f"{f('corrected_lai', 'pct_hi95'):+.2f}) | {f('corrected_lai', 'MK_p'):.4f} |")
md += ['', f'Numbers that change: {int(S.changed.sum())} of {len(S)} (lai_sensitivity_summary.csv).']
open(os.path.join(OUT_DRIVE, 'lai_sensitivity_summary.md'), 'w', encoding='utf-8').write('\n'.join(md) + '\n')
print('\n'.join(md)); print(S[S.changed].to_string(index=False, max_colwidth=60))
for o in (OUT_A, OUT_A2, OUT_B):                                 # logs of the three passes kept in Drive
    shutil.copytree(os.path.join(o, 'logs'), os.path.join(OUT_DRIVE, os.path.basename(o) + '_logs'), dirs_exist_ok=True)
json.dump(dict(notebook='19_install_and_reproduce', run_utc=datetime.datetime.utcnow().isoformat(timespec='seconds') + 'Z',
               checks_ok=CHECKS_OK, installed=INSTALLED, noaa_url=NOAA_URL, noaa_retrieved=RETRIEVED, noaa_sha256=sha256(CO2_NEW),
               noaa_creation_lines=created, passA_exit=rcA, passA_csv_different=int(len(DIFF_A)), passA2_exit=rcA2, passB_exit_analysis_steps=rcB,
               passB_compare_exit=rcBcmp, passB_compare_finished=CMP_B_DONE,
               laicorr_sha256_16=LAICORR_SHA16, fixes_sha256_16=FIXES_SHA16),
          open(os.path.join(OUT_DRIVE, 'provenance_notebook19.json'), 'w'), indent=1, default=str)

# %% [markdown]
# ## Cell 9 · Optional: promote the corrected results into your repository copy and pack the release
#
# Runs only with `PROMOTE = True`, and only if pass A reproduced the published results, every analysis step of pass B
# ended without error and its comparison with `results/v10` finished (the differences there are the expected ones).
# In `REPO_DRIVE`:
# * `results/v10` → copied to `results/v10_modis_lai_2024_2025` (once; kept afterwards);
# * every file of pass B that has a published counterpart replaces it in `results/v10` (logs, caches and new files are
#   not copied; the list is printed);
# * `data/external/co2_annmean_gl.txt` ← the NOAA file, with `co2_annmean_gl.SOURCE.txt` (address, retrieval time,
#   creation line, SHA-256);
# * `data/DATA_SHA256SUMS.txt` and `results/MANIFEST.sha256` are rewritten.
#
# With `PACK = True` as well, the five zips (`tools/data_release.py pack`) and `SHA256SUMS_assets.txt` go to `RELEASE_DIR`,
# and the checksum list is copied to `data/RELEASE_ASSETS_SHA256SUMS.txt`.

# %%
# 9 · Promote and pack ---------------------------------------------------------------------------------------------------------
def run(cmd):
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True); print(' '.join(cmd[1:]), '->', r.returncode)
    print(r.stdout[-2500:], r.stderr[-1500:]); assert r.returncode == 0, cmd
if PROMOTE:
    why = ([f'pass A ended with exit code {rcA}'] if rcA else []) + \
          ([f'{len(DIFF_A)} CSV files of pass A differ from the published copies: ' + ', '.join(DIFF_A.file)] if len(DIFF_A) else []) + \
          ([f'an analysis step of pass B failed (exit code {rcB})'] if rcB else []) + \
          ([] if CMP_B_DONE else ['the comparison of pass B with results/v10 did not finish (see its compare.log)'])
    for o in (OUT_A, OUT_B):                          # the log of every step that failed
        for f in sorted(glob.glob(os.path.join(o, 'logs', '*.log'))):
            t = open(f, encoding='utf-8', errors='replace').read()
            if 'Traceback' in t or 'FAILED' in t:
                print('=====', os.path.relpath(f, WORK)); print(t[-2500:])
    assert not why, 'not promoted: ' + '; '.join(why)
    RES = os.path.join(REPO, 'results')
    if not os.path.exists(os.path.join(RES, 'v10_modis_lai_2024_2025')):
        shutil.copytree(os.path.join(RES, 'v10'), os.path.join(RES, 'v10_modis_lai_2024_2025'))
    promoted, skipped = [], []
    for root, _, files in os.walk(OUT_B):
        rr = os.path.relpath(root, OUT_B)
        if rr.split(os.sep)[0] in ('logs', '_work') or 'cache_' in rr:
            continue
        for f in files:
            rel = os.path.normpath(os.path.join(rr, f)); dst = os.path.join(RES, 'v10', rel)
            (promoted if os.path.exists(dst) else skipped).append(rel)
            if os.path.exists(dst):
                shutil.copyfile(os.path.join(OUT_B, rel), dst)
    print(f'{len(promoted)} files of results/v10 replaced by pass B; {len(skipped)} files of pass B without a published copy (not copied):')
    print('  ' + '\n  '.join(sorted(skipped)[:40]))
    open(os.path.join(DATA, 'external', 'co2_annmean_gl.SOURCE.txt'), 'w', encoding='utf-8').write(
        f'co2_annmean_gl.txt: NOAA GML global annual mean CO2, downloaded by notebook 19 from {NOAA_URL}\n'
        f'retrieved: {RETRIEVED}\nfile creation line(s): {" | ".join(created) if created else "(none in the header)"}\n'
        f'sha256: {sha256(CO2_NEW)}\nreplaces the copy taken on 13 October 2026 from the mirror https://github.com/datasets/co2-ppm\n')
    run([sys.executable, 'tools/data_release.py', 'sums'])
    run([sys.executable, 'reproduce.py', 'check-data'])
    run([sys.executable, 'tools/verify_manifest.py', '--write'])
    if PACK:
        rel_local = os.path.join(WORK, 'release'); run([sys.executable, 'tools/data_release.py', 'pack', rel_local])
        shutil.copyfile(os.path.join(rel_local, 'SHA256SUMS_assets.txt'), os.path.join(DATA, 'RELEASE_ASSETS_SHA256SUMS.txt'))
        os.makedirs(RELEASE_DIR, exist_ok=True)
        for f in sorted(os.listdir(rel_local)):
            shutil.copyfile(os.path.join(rel_local, f), os.path.join(RELEASE_DIR, f)); print('->', os.path.join(RELEASE_DIR, f))
    for sub in ('results', os.path.join('data', 'external'), os.path.join('data', 'DATA_SHA256SUMS.txt'),
                os.path.join('data', 'RELEASE_ASSETS_SHA256SUMS.txt')):
        s, d = os.path.join(REPO, sub), os.path.join(REPO_DRIVE, sub)
        if os.path.isdir(s):
            shutil.copytree(s, d, dirs_exist_ok=True)
        else:
            shutil.copyfile(s, d)
    print('REPO_DRIVE updated: results/, data/external/, data/DATA_SHA256SUMS.txt, data/RELEASE_ASSETS_SHA256SUMS.txt')
else:
    print('PROMOTE = False: nothing written to', REPO_DRIVE)

# %% [markdown]
# ## Cell 10 · Release (by hand)
#
# 1. Read `lai_sensitivity_summary.md` and update the texts (green passages of the report and the manuscript).
# 2. Commit `REPO_DRIVE` to GitHub (`results/`, `data/DATA_SHA256SUMS.txt`, `data/RELEASE_ASSETS_SHA256SUMS.txt`,
#    `data/external/`).
# 3. GitHub → Releases → *Draft a new release*, tag `v1.0.0`; attach the five zips and `SHA256SUMS_assets.txt` of
#    `RELEASE_DIR`. Check the download with `python tools/data_release.py fetch https://github.com/<owner>/MEGAN_SC/releases/download/v1.0.0`
#    and `python reproduce.py check-data` in a fresh clone.
# 4. Zenodo: either the GitHub integration (Zenodo → GitHub → switch the repository on before publishing the release)
#    or a manual upload of the same files. Then write the DOI into `CITATION.cff`, `README.md` and
#    `DATA_AVAILABILITY.md`, and the release date into `CITATION.cff`.
#
# To repeat the MODIS control later, from the repository root:
#
#     python tools/lai_variant_view.py
#     MEGAN_SC_RUN=data/run_v10_modis_view MEGAN_SC_RESULTS=results/v10_modis_lai_2024_2025 python reproduce.py all

# %%
# 10 · Status ----------------------------------------------------------------------------------------------------------------
print(pd.DataFrame(CHECKS).to_string(index=False))
print('installed:', INSTALLED, '| pass A exit', rcA, '| pass A2 exit', rcA2, '| pass B analysis steps exit', rcB,
      '| pass B comparison exit', rcBcmp, '(1 = differs from the MODIS results, expected)')
print('summaries in', OUT_DRIVE)
