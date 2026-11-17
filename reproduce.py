#!/usr/bin/env python3
"""reproduce.py - recompute every table, figure and number of the Santa Catarina example (2003-2025) from the model
run shipped in data/run_v10, and compare them with the published copies in results/v10.

    python reproduce.py list              # the steps, in order, with what each one writes
    python reproduce.py paths             # the folders in use (environment variables MEGAN_SC_* override them)
    python reproduce.py check-data        # SHA-256 of every data file against data/DATA_SHA256SUMS.txt
    python reproduce.py all               # every step below, then `compare`
    python reproduce.py nb13 maps compare # any subset, in the order given

Every step writes under reproduced/ (MEGAN_SC_OUT) and logs to reproduced/logs/<step>.log; nothing in data/ or
results/ is modified. The tutorial (docs/tutorial/05_reproduce_the_example.md) explains each step and its outputs.
Re-running the model itself (notebook 10, for Santa Catarina or any other region and period) needs Google Earth Engine
and is described in docs/tutorial/06_run_the_model.md.
"""
import glob, hashlib, os, shutil, subprocess, sys, time

REPO = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(REPO, 'analysis', 'v10'))
import _paths as P  # noqa: E402

PY = sys.executable
NBDIR = os.path.join(REPO, 'notebooks')
RUNNER = os.path.join(REPO, 'tools', 'run_notebook.py')


def q(p):
    """A path as a Python string literal for the notebook overrides."""
    return repr(os.path.abspath(p))


def prepare_tropomi_cache():
    """Notebook 14 skips Earth Engine when every month is in OUT_DIR/cache_<TAG>: the cache shipped in data/."""
    d = P.out('notebook14', 'cache_' + P.TAG)
    for f in sorted(glob.glob(os.path.join(P.TROPOMI, 'hcho_*.npz'))):
        t = os.path.join(d, os.path.basename(f))
        if not os.path.exists(t):
            shutil.copyfile(f, t)
    return os.path.dirname(d)


def prepare_modis_root():
    """Folder with the MODIS-LAI record of every year (tools/lai_variant_view.py): data/run_v10 itself while its 2024
    and 2025 are the MODIS years, else the view that takes those two years from data/run_v10_modis_lai_2024_2025."""
    sys.path.insert(0, os.path.join(REPO, 'tools'))
    import lai_variant_view as V
    return V.modis_root(P.RUN, P.RUN_MODIS_ARCHIVE, P.RUN_MODIS)


def nb(name, overrides):
    return [PY, RUNNER, os.path.join(NBDIR, name), ';'.join(f'{k}={v}' for k, v in overrides.items())]


def script(rel, *args):
    return [PY, os.path.join(REPO, rel), *args]


# name: (description, writes, command factory)
STEPS = {
    'integrity': ('steps per day of every month and annual = sum of monthly files (notebook-10 run)',
                  'logs/integrity.log', lambda: script('analysis/v10/check_run_integrity.py')),
    'nb11': ('notebook 11: the six automatic checks C0-C5 for every year',
             'notebook11/verify_v10_report.csv, annual_totals_SC_v10_all_classes.csv',
             lambda: nb('11_verify_run.ipynb', dict(NEW_ROOT=q(P.RUN),
                                                    BOUNDARY=q(os.path.join(P.BND, 'sc_state_ibge.gpkg')),
                                                    REPORT_DIR=q(P.out('notebook11'))))),
    'nb13': ('notebook 13: annual totals by class and group, trends, correlations, drivers, seasonal cycle, mesoregions, municipalities, comparison of 2019 with the previous state inventory, key_numbers.csv, files fig1-8 and fig_detrend (drafts of the maps), grid cache',
             'notebook13_meso2024/*',
             lambda: nb('13_santa_catarina_analysis.ipynb', dict(ROOT=q(P.RUN), BND=q(P.BND),
                                                                OUT_DIR=q(P.out('notebook13_meso2024'))))),
    'nb14': ('notebook 14: the four TROPOMI tests (from the cached monthly HCHO fields, no Earth Engine)',
             'notebook14/*',
             lambda: nb('14_tropomi_hcho_evaluation.ipynb', dict(MODEL_ROOT=q(P.RUN), BND=q(P.BND),
                                                                OUT_DIR=q(prepare_tropomi_cache())))),
    'nb15': ('notebook 15: sensitivity test S3 (soil-moisture factor) from the fix logs',
             'notebook15/S3_fixlog_*.csv',
             lambda: nb('15_sensitivity_2019.ipynb', dict(BASE_ROOT=q(P.RUN), HOURLY_ROOT=q(P.out('_work', 'none')),
                                                         LAIDEF_ROOT=q(P.out('_work', 'none')), BND=q(P.BND),
                                                         OUT_DIR=q(P.out('notebook15'))))),
    's1': ('test S1: 3-hourly against hourly sampling on 28 days of 2019 (two test runs of notebook 10)',
           'notebook15/s1/*', lambda: script('analysis/v10/s1_compare.py')),
    'c3': ('check C3 in detail: filled coastal cells, LAI-adjusted test, effect on the totals',
           'diagnostics/c3/*', None),
    'nb14diag': ('post hoc diagnostics of the TROPOMI comparison (D1-D7)', 'diagnostics/nb14/*',
                 lambda: script('analysis/v10/nb14_diagnostics.py')),
    'report': ('numbers by class, mesoregion and municipality, standardised regression and temperature sensitivity, and the figure of the regression (pt_report/fig10)', 'report/*, figures/pt_report/fig10.png', None),
    'maps': ('every map (fig5, fig7, fig8 in PT and EN; pt_report/fig9, fig12, figA2-figA9) and the TROPOMI figure (fig9)',
             'figures/{pt,en,pt_report}/*, figures/maps_class_limits.csv', None),
    'meso2024': ('MapBiomas 2024 land-cover shares by mesoregion, and the mesoregions as unions of the 2024 municipalities (islands included)', 'boundaries/*',
                 None),
    'final': ('additional diagnostics (FDR control, Hamed-Rao check, LAI of 2024-2025, parts of the trends, CO2 inhibition, robustness of the trends, UTC month boundaries, ...)',
              'final_review/*', lambda: script('analysis/v10/review/final_review_diagnostics.py')),
    'extra': ('study-area map (MapBiomas 2024 and the model growth forms) and the verification flowchart',
              'figures/{pt,en}/fig_study_area.png, fig_verification.png, figures/study_area_shares.csv', None),
    'posthoc': ('post hoc diagnostics (planted-forest emission factors, LAI of 2024-2025), the ecotype emission factors of the state and the land use of the municipalities',
                'posthoc_review/*', None),
    'laicorr': ('VIIRS check of the 2024-2025 MODIS LAI (notebook 16 cache): hindcast of the candidate corrections, corrected LAI of 2024-2025 used by notebook 17; reads the MODIS-LAI record of every year (tools/lai_variant_view.py)',
                'lai_viirs_check/lai_correction*, lai_viirs_corrected_2024_2025.npz',
                lambda: script('analysis/v10/review/lai_viirs_correction.py', prepare_modis_root())),
    'compare': ('compare every reproduced file with its published copy in results/v10', 'logs/compare.log', None),
}
MULTI = {
    'c3': ['analysis/v10/c3/diag_c3.py', 'analysis/v10/c3/diag_c3_lai.py', 'analysis/v10/c3/c3_impact_on_totals.py'],
    'report': ['analysis/v10/report/10_report_numbers_v10.py', 'analysis/v10/report/11_report_figures_v10.py'],
    'meso2024': ['analysis/v10/maps/landcover_mesoregion.py', 'analysis/v10/meso2024/make_sc_mesoregions_2024.py'],
    'maps': ['analysis/v10/maps/12_maps_wf_style.py', 'analysis/v10/fig9_tropomi.py'],
    'extra': ['analysis/v10/maps/study_area_map.py', 'analysis/v10/schematics/verification_flowchart.py'],
    'posthoc': ['analysis/v10/review/posthoc_review_diagnostics.py', 'analysis/v10/review/ecotype_tree_ef.py',
                'analysis/v10/review/landuse_context.py'],
    'final': ['analysis/v10/review/final_review_diagnostics.py', 'analysis/v10/review/trend_robustness.py',
              'analysis/v10/review/month_boundary.py'],
}
ORDER = ['integrity', 'nb11', 'nb13', 'nb14', 'nb15', 's1', 'c3', 'nb14diag', 'report', 'maps', 'meso2024', 'final',
         'extra', 'posthoc', 'laicorr', 'compare']


def env():
    e = dict(os.environ)
    for k, v in (('MEGAN_SC_RUN', P.RUN), ('MEGAN_SC_RUN_MODIS_ARCHIVE', P.RUN_MODIS_ARCHIVE), ('MEGAN_SC_RUN_MODIS', P.RUN_MODIS),
                 ('MEGAN_SC_BND', P.BND), ('MEGAN_SC_TROPOMI', P.TROPOMI), ('MEGAN_SC_S1', P.S1), ('MEGAN_SC_LAIVIIRS', P.LAIVIIRS),
                 ('MEGAN_SC_RESULTS', P.RESULTS), ('MEGAN_SC_OUT', P.OUT),
                 ('MPLBACKEND', 'Agg')):
        e.setdefault(k, v)
    e['PYTHONPATH'] = os.pathsep.join([os.path.join(REPO, 'src'), e.get('PYTHONPATH', '')]).strip(os.pathsep)
    return e


def run(step):
    log = P.out('logs', f'{step}.log')
    if step == 'compare':
        cmds = [[PY, os.path.join(REPO, 'tools', 'compare_results.py'), P.OUT, P.RESULTS, f'boundaries={P.BND}']]
    elif step in MULTI:
        cmds = [script(s) for s in MULTI[step]]
    else:
        try:
            cmds = [STEPS[step][2]()]
        except Exception as e:                       # preparing the inputs of the step failed
            print(f'{step:10s} FAILED while preparing its inputs: {e!r}'); return False
    t0 = time.time(); ok = True
    with open(log, 'w', encoding='utf-8') as fh:
        for c in cmds:
            fh.write('$ ' + ' '.join(c) + '\n'); fh.flush()
            r = subprocess.run(c, cwd=REPO, env=env(), stdout=fh, stderr=subprocess.STDOUT)
            if r.returncode != 0:
                ok = False; break
    print(f'{step:10s} {"ok    " if ok else "FAILED"} {time.time() - t0:7.1f} s   log: {os.path.relpath(log, REPO)}', flush=True)
    if not ok:
        print(open(log, encoding='utf-8', errors='replace').read()[-3000:])
    return ok


def check_data():
    sums = os.path.join(REPO, 'data', 'DATA_SHA256SUMS.txt')
    bad = n = 0
    for line in open(sums, encoding='utf-8'):
        h, rel = line.rstrip('\n').split('  ', 1); f = os.path.join(REPO, 'data', rel); n += 1
        if not os.path.exists(f):
            print('missing', rel); bad += 1; continue
        d = hashlib.sha256()
        with open(f, 'rb') as fh:
            for b in iter(lambda: fh.read(1 << 20), b''):
                d.update(b)
        if d.hexdigest() != h:
            print('checksum differs', rel); bad += 1
    print(f'{n} data files checked, {bad} missing or different')
    return bad == 0


def main(argv):
    if not argv or argv[0] in ('-h', '--help', 'help'):
        print(__doc__); return 0
    if argv == ['list']:
        for k in ORDER:
            print(f'{k:10s} {STEPS[k][0]}\n{"":10s} -> reproduced/{STEPS[k][1]}')
        return 0
    if argv == ['paths']:
        for k in ('REPO', 'RUN', 'RUN_MODIS_ARCHIVE', 'RUN_MODIS', 'BND', 'TROPOMI', 'S1', 'RESULTS', 'OUT'):
            print(f'{k:8s} {getattr(P, k)}')
        return 0
    if argv == ['check-data']:
        return 0 if check_data() else 1
    steps = ORDER if argv == ['all'] else argv
    unknown = [s for s in steps if s not in STEPS]
    if unknown:
        print('unknown step(s):', unknown, '| known:', ' '.join(ORDER)); return 2
    os.makedirs(P.OUT, exist_ok=True)
    failed = [s for s in steps if not run(s)]
    print('all steps ok' if not failed else f'failed: {failed}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
