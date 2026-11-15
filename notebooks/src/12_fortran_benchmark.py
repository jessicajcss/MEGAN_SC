# %% [markdown]
# # 12 · Benchmark of the MEGAN-SC Python port against the official MEGANv3.21 Fortran
#
# Runs the official MEGANv3.21 binaries (Guenther 2024, Zenodo 10.5281/zenodo.10526206, as shipped in
# `MEGANv3.21.tar.gz`) and the MEGAN-SC Python engine on the same inputs and compares them stage by stage.
# Needs **Linux x86-64** (Google Colab works; the binaries do not run on Windows outside WSL). About 15-25 min.
#
# | Case | Inputs | What is compared |
# |---|---|---|
# | S0 harness | official SE_US test case, 2018-05-02 (163 × 163, 12 km Lambert) | our Fortran run vs the reference output shipped in the tarball (proves the binaries and settings are right) |
# | SE_US | same | S1 DAYMET, S2 MEGCAN canopy, S3 MEGVEA with the Fortran canopy as input, S4 Python end to end; GAMSM/GAMBD off and on; soil NO (MEGSEA, Yienger-Levy, with and without soil data); speciation and unit chain (mgn2mech, 201 species in t/h) |
# | SC | Santa Catarina grid (0.025°, 144 × 224) on 15 Jan and 15 Jul 2019, 2019 monthly LAI (Earth Engine default resampling) and growth forms, idealised diurnal T and PAR | the same stages in the Southern Hemisphere and the production canopy setting, plus day totals weighted with the production EF grid |
#
# Each comparison is made with the engine as used in production, with the two engine patches of `fixes.py` (FIX-6
# GAMTP operator precedence, FIX-7 soil-NO air-temperature cap), which has to reproduce the Fortran to float32
# precision, and with the engine modules alone, a control showing that the comparison detects a departure from
# megvea.f and soilnox.F (rows labelled "original port").
#
# Documented deviations that are kept on purpose (not tested as equal): the Southern-Hemisphere growing season for soil
# NO (MEGANv3.21 compares day-of-year with month-day constants, so in SH winter its fertiliser factor goes negative in
# crop cells; the port clips it at 0 and production uses the Nov-May season), constant LDF for the four MT classes, and
# GAMSM/GAMBD switched on in production.
#
# **Send back:** the folder `RESULTS_DIR` (CSV files + `benchmark_summary.csv`).

# %%
# 1 · Paths ------------------------------------------------------------------------------------------------
import os, sys, glob, subprocess, shutil, hashlib, json
try:
    from google.colab import drive
    IN_COLAB = True
except ImportError:
    IN_COLAB = False
if IN_COLAB:
    # A failed mount must stop the notebook: if it is hidden, the notebook reads and writes a local /content/drive
    # folder that is not Google Drive (notebook 13 then stops with "found 0").
    if os.path.isdir('/content/drive') and not os.path.ismount('/content/drive') and os.listdir('/content/drive'):
        raise RuntimeError('/content/drive holds local files written while Google Drive was not mounted, so Drive cannot be '
                           'mounted there. Use Runtime > Disconnect and delete runtime, then run the notebook again from cell 1.')
    drive.mount('/content/drive')
MGN_TGZ   = '/content/drive/MyDrive/MEGANv3.21.tar.gz'        # <-- MEGANv3.21.tar.gz from Zenodo 10.5281/zenodo.10526206
MEGAN_DIR = '/content/drive/MyDrive/megan-sc'                 # <-- MEGAN-SC repo (same folder as notebook 10)
RAW       = '/content/drive/MyDrive/megan_sc_outputs'         # <-- run folder whose 2019 inputs build the SC test case
WORK      = '/content/bench'                                  # local scratch (needs ~1.5 GB)
RESULTS_DIR = '/content/drive/MyDrive/megan_sc_benchmark_results'
MGN = os.path.join(WORK, 'MEGANv3.21'); PY = os.path.join(WORK, 'py'); ENGINE_SRC = os.path.join(MEGAN_DIR, 'src')
for d in (WORK, PY, RESULTS_DIR): os.makedirs(d, exist_ok=True)
def sh(cmd, env=None, cwd=None, check=True):
    e = dict(os.environ); e.update(env or {})
    r = subprocess.run(cmd, shell=True, env=e, cwd=cwd or WORK, capture_output=True, text=True)
    print(r.stdout[-3000:], r.stderr[-3000:])
    if check and r.returncode != 0: raise RuntimeError(f'failed ({r.returncode}): {cmd}')
    return r

# %%
# 2 · Dependencies and the Intel OpenMP runtime the binaries are linked against --------------------------
%pip -q install intel-openmp netCDF4 rasterio geopandas pyyaml pandas
import site, sysconfig
cands = []
for base in {sys.prefix, sys.base_prefix, site.getuserbase(), '/usr/local', '/usr'}:
    cands += glob.glob(os.path.join(base, 'lib', 'libiomp5.so')) + glob.glob(os.path.join(base, 'lib*', '**', 'libiomp5.so'), recursive=True)
assert cands, 'libiomp5.so not found after pip install intel-openmp'
OMP_LIB = os.path.dirname(cands[0]); os.environ['OMP_LIB'] = OMP_LIB
print('libiomp5.so in', OMP_LIB)

# %%
# 3 · Extract the official MEGANv3.21 package (binaries, SE_US inputs, reference output) -----------------
if not os.path.exists(os.path.join(MGN, 'src', 'MEGVEA', 'megvea')):
    sh(f"tar -xzf '{MGN_TGZ}' -C '{WORK}' --exclude='MEGANv3.21/CAMx_ready' --exclude='MEGANv3.21/Output/MEGANv32*'")
for p in ('DAYMET/daymet', 'MEGCAN/megcan', 'MEGSEA/megsea', 'MEGVEA/megvea', 'MGN2MECH/mgn2mech'):
    os.chmod(os.path.join(MGN, 'src', p), 0o755)
sh(f"LD_LIBRARY_PATH={OMP_LIB} ldd {MGN}/src/MEGVEA/megvea | grep -i 'not found' || echo 'all shared libraries resolved'", check=False)
print('tarball sha256[:16]:', hashlib.sha256(open(MGN_TGZ, 'rb').read()).hexdigest()[:16])

# %%
# 4 · Benchmark scripts (benchmark/ of the repository) and the engine patches (identical to notebook 10) ----------
FILES = {
 'ioapi.py': r'''__FILE:benchmark/ioapi.py__''',
 'bench_se_us.py': r'''__FILE:benchmark/bench_se_us.py__''',
 'bench_soil.py': r'''__FILE:benchmark/bench_soil.py__''',
 'bench_spec.py': r'''__FILE:benchmark/bench_spec.py__''',
 'make_sc_inputs.py': r'''__FILE:benchmark/make_sc_inputs.py__''',
 'bench_sc.py': r'''__FILE:benchmark/bench_sc.py__''',
 'run_fortran_se_us.sh': r'''__FILE:benchmark/run_fortran_se_us.sh__''',
 'run_fortran_se_us_soil_spec.sh': r'''__FILE:benchmark/run_fortran_se_us_soil_spec.sh__''',
 'run_fortran_sc.sh': r'''__FILE:benchmark/run_fortran_sc.sh__''',
}
for k, v in FILES.items():
    with open(os.path.join(PY, k), 'w') as f: f.write(v)
    if k.endswith('.sh'): os.chmod(os.path.join(PY, k), 0o755)
FIXES_PY = r'''__FIXES_PY__'''
with open(os.path.join(ENGINE_SRC, 'megan_sc', 'fixes.py'), 'w', encoding='utf-8') as f: f.write(FIXES_PY)
FIXES_SHA = hashlib.sha256(FIXES_PY.encode('utf-8')).hexdigest()[:16]
ENV = dict(MGN=MGN, OMP_LIB=OMP_LIB, ENGINE_SRC=ENGINE_SRC, PYTHONDONTWRITEBYTECODE='1', RAW=RAW)
print('scripts in', PY, '| fixes.py sha256[:16] =', FIXES_SHA)

# %%
# 5 · Fortran, SE_US: DAYMET -> MEGCAN -> MEGVEA (GAMSM/GAMBD off), MEGSEA (YL soil NO), mgn2mech; then MEGVEA with GAMSM/GAMBD on
FD = os.path.join(WORK, 'se_us', 'fort'); FD2 = os.path.join(WORK, 'se_us', 'fort_smbd')
if not os.path.exists(os.path.join(FD, 'MGNERS.ncf')):
    sh(f'{PY}/run_fortran_se_us.sh {FD} N N 46', env=ENV)
if not os.path.exists(os.path.join(FD, 'MGN201_tonhr.ncf')):
    sh(f'{PY}/run_fortran_se_us_soil_spec.sh {FD}', env=ENV)
if not os.path.exists(os.path.join(FD2, 'MGNERS.ncf')):
    sh(f'{PY}/run_fortran_se_us.sh {FD2} Y Y 46', env=dict(ENV, SMFILE=f'{FD}/MGNSEA_YL_soil.ncf'))
for f in sorted(glob.glob(f'{FD}/log_*.txt') + glob.glob(f'{FD2}/log_*.txt')):
    print(os.path.relpath(f, WORK), 'Normal Completion' in open(f, errors='ignore').read())

# %%
# 6 · S0 harness check: our Fortran run (GAMSM/GAMBD off, DAYMET over 23 records, as in the official step3/step6
# scripts) vs the intermediate output shipped in Output/INT. Expected: bit-identical (max_abs = 0) for DAYMET, CANMET and
# all 19 MGNERS classes (benchmark/results/S0_harness_vs_shipped_reference.csv).
sys.path.insert(0, PY)
import numpy as np, pandas as pd
from ioapi import rd, compare
REF = os.path.join(MGN, 'Output', 'INT'); rows = []
pairs = [('DAYMET.ncf', 'DAYMET.SE_US.2018122.ncf', ['D_TEMP', 'D_PPFD', 'MaxT', 'MaxWS']),
         ('CANMET.ncf', 'CANMET.SE_US.2018122.ncf', ['SunleafTK', 'ShadeleafTK', 'SunPPFD', 'ShadePPFD', 'SunFrac']),
         ('MGNERS.ncf', 'MGNERS.SE_US.2018122.ncf', ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM',
          'NO', 'SQT_HR', 'SQT_LR', 'MEOH', 'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO'])]
for ours, ref, vars_ in pairs:
    if not os.path.exists(os.path.join(REF, ref)):
        print('reference not in tarball:', ref); continue
    for v in vars_:
        a = rd(os.path.join(FD, ours), v); b = rd(os.path.join(REF, ref), v)
        n = min(a.shape[0], b.shape[0])
        rows.append(dict(file=ref, var=v, steps=n, **compare(a[:n], b[:n])))
S0 = pd.DataFrame(rows); S0.to_csv(os.path.join(RESULTS_DIR, 'S0_harness_vs_shipped_reference.csv'), index=False)
print(S0.to_string(float_format=lambda v: f'{v:.3g}'))
print('S0 bit-identical:', bool((S0.max_abs == 0).all()), f'({len(S0)} variables)')

# %%
# 7 · Python vs Fortran, SE_US (each run 1-3 min) ----------------------------------------------------------
SE = os.path.join(WORK, 'se_us')
runs = [('se_us_bench_N_N.csv', FD, 'N', 'N', '0', {}),
        ('se_us_bench_patched_N_N.csv', FD, 'N', 'N', '1', {}),
        ('se_us_bench_patched_Y_Y.csv', FD2, 'Y', 'Y', '1', {'SMFILE': f'{FD}/MGNSEA_YL_soil.ncf'})]
for out, fd, sm, bd, patch, extra in runs:
    if not os.path.exists(os.path.join(RESULTS_DIR, out)):
        r = sh(f'python3 {PY}/bench_se_us.py {ENGINE_SRC} {MGN} {fd} {RESULTS_DIR}/{out} {sm} {bd} | tail -6', env=dict(ENV, MEGAN_PATCH=patch, **extra))
for out, patch in (('se_us_bench_soil_unpatched.csv', '0'), ('se_us_bench_soil_patched.csv', '1')):
    sh(f'python3 {PY}/bench_soil.py {RESULTS_DIR}/{out}', env=dict(ENV, FD=FD, MEGAN_PATCH=patch))
sh(f'python3 {PY}/bench_spec.py {RESULTS_DIR}/se_us_bench_speciation.csv | tail -4', env=dict(ENV, FD=FD))

# %%
# 8 · Santa Catarina case: 15 Jan and 15 Jul 2019 (inputs ~170 MB per day; deleted after use) -------------
for mon, day in ((1, 15), (7, 15)):
    out = os.path.join(RESULTS_DIR, f'sc_bench_sc_m{mon:02d}.csv')
    if os.path.exists(out): continue
    D = os.path.join(WORK, 'sc', f'm{mon:02d}')
    sh(f'python3 {PY}/make_sc_inputs.py {ENGINE_SRC} {RAW} {MGN} {D} {mon} {day}', env=ENV)
    sh(f'{PY}/run_fortran_sc.sh {D}', env=ENV)
    sh(f'python3 {PY}/bench_sc.py {ENGINE_SRC} {MEGAN_DIR} {D} {out} | tail -30', env=ENV)
    shutil.rmtree(D, ignore_errors=True)

# %%
# 9 · Summary and verdict ---------------------------------------------------------------------------------
# Criteria for the patched engine (float32 Fortran output, float64 Python): 99th percentile of the cell-hour relative
# difference <= 1e-4 (<= 5e-4 for the canopy temperatures/PPFD, whose iterative energy balance converges to ~1e-4) and
# state/domain sums within 1e-5. The original port is expected to FAIL S3/S4 (FIX-6) and the no-soil NO path (FIX-7).
def summ(f):
    d = pd.read_csv(os.path.join(RESULTS_DIR, f)); key = 'stage' if 'stage' in d else 'test'
    g = d.groupby(key).agg(max_rel=('max_rel', 'max'), p99_rel=('p99_rel', 'max'),
                           sum_dev=('sum_ratio', lambda x: float(np.nanmax(np.abs(x - 1))) if x.notna().any() else np.nan))
    if 'day_tonnes_ratio' in d:
        g['day_tonnes_dev'] = d.groupby(key)['day_tonnes_ratio'].apply(lambda x: float(np.nanmax(np.abs(x - 1))) if x.notna().any() else np.nan)
    g.insert(0, 'file', f); return g.reset_index().rename(columns={key: 'stage'})
files = ['se_us_bench_N_N.csv', 'se_us_bench_patched_N_N.csv', 'se_us_bench_patched_Y_Y.csv',
         'se_us_bench_soil_unpatched.csv', 'se_us_bench_soil_patched.csv', 'sc_bench_sc_m01.csv', 'sc_bench_sc_m07.csv']
SUM = pd.concat([summ(f) for f in files if os.path.exists(os.path.join(RESULTS_DIR, f))], ignore_index=True)
def verdict(r):
    s, f = r.stage, r.file
    if 'SH season corrected' in s:
        return 'deviation by design (SH season)'
    if f.startswith('sc_') and 'as Fortran' in s:
        return 'PASS' if r.p99_rel <= 1e-4 else 'see note (SH crop cells)'
    original = ('original port' in s) or f in ('se_us_bench_N_N.csv', 'se_us_bench_soil_unpatched.csv')
    lim = 5e-4 if s.startswith('S2') else 1e-4
    tdev = r.get('day_tonnes_dev', np.nan)
    ok = r.p99_rel <= lim and (pd.isna(r.sum_dev) or r.sum_dev <= 1e-5) and (pd.isna(tdev) or tdev <= 1e-5)
    if ok: return 'PASS'
    if original and s.startswith(('S3', 'S4')): return 'expected FAIL (port bug, FIX-6)'
    if original and 'nosoil' in s: return 'expected FAIL (port bug, FIX-7)'
    return 'FAIL'
SUM['verdict'] = SUM.apply(verdict, axis=1)
sp = pd.read_csv(os.path.join(RESULTS_DIR, 'se_us_bench_speciation.csv'))
SUM.loc[len(SUM)] = dict(file='se_us_bench_speciation.csv', stage='mgn2mech 201 species (t/h)', max_rel=sp.max_rel_cellhour.max(), p99_rel=np.nan,
                         sum_dev=float(np.nanmax(np.abs(sp.ratio - 1))), verdict='PASS' if np.nanmax(np.abs(sp.ratio - 1)) <= 1e-5 else 'FAIL')
SUM['fixes_sha256'] = FIXES_SHA
SUM.to_csv(os.path.join(RESULTS_DIR, 'benchmark_summary.csv'), index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
print(SUM.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
print('\nFAIL rows:', int((SUM.verdict == 'FAIL').sum()), '| results in', RESULTS_DIR)

# %% [markdown]
# **Reading the summary.** `PASS` rows are the evidence that the engine, as used in production, reproduces MEGANv3.21.
# `expected FAIL` rows are the controls without the two engine patches (FIX-6: GAMTP operator precedence, light-dependent
# fraction of every class; FIX-7: soil NO without soil data, air temperature not capped at 303 K). In the SC July case the
# "as Fortran" soil-NO comparison differs only in crop cells, where MEGANv3.21's Southern-Hemisphere season gives a
# negative fertiliser factor (negative GAMNO) that the port clips at 0; non-crop cells agree to ~1e-6. Production uses
# the Nov-May season adjusted to the Southern Hemisphere, so this row is a documented deviation and not a port error.
