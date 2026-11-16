"""The region wizard (run_your_region/region_setup.ipynb): the committed notebook is the output of its generator, every
code cell parses, and the period of STEP 1 (first and last month and year) is checked and expanded as documented."""
import ast, json, os, shutil, subprocess, sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NB = os.path.join(REPO, 'run_your_region', 'region_setup.ipynb')
GEN = os.path.join(REPO, 'scripts', 'make_region_setup_notebook.py')


def cells():
    return json.load(open(NB, encoding='utf-8'))['cells']


def settings_source():
    src = [''.join(c['source']) for c in cells() if c['cell_type'] == 'code']
    return next(s for s in src if 'START_MONTH' in s and 'PERIOD_MONTHS' in s)


def run_settings(**values):
    src = settings_source()
    for k, v in values.items():                      # replace the value typed in STEP 1
        lines = src.splitlines()
        i = next(j for j, l in enumerate(lines) if l.startswith(k + ' ') or l.startswith(k + '='))
        lines[i] = f'{k} = {v!r}'
        src = '\n'.join(lines)
    ns = {}
    exec(compile(src, 'step1', 'exec'), ns)
    return ns


def test_notebook_is_generator_output(tmp_path):
    os.makedirs(tmp_path / 'scripts'); os.makedirs(tmp_path / 'notebooks' / 'src')
    shutil.copy(GEN, tmp_path / 'scripts')
    shutil.copy(os.path.join(REPO, 'notebooks', 'src', '_setup_repo.py'), tmp_path / 'notebooks' / 'src')
    subprocess.run([sys.executable, str(tmp_path / 'scripts' / os.path.basename(GEN))], check=True, capture_output=True)
    assert open(tmp_path / 'run_your_region' / 'region_setup.ipynb', 'rb').read() == open(NB, 'rb').read()


def test_code_cells_parse():
    for c in cells():
        if c['cell_type'] == 'code':
            ast.parse(''.join(c['source']))


def test_period_across_a_year_boundary():
    ns = run_settings(START_YEAR=2019, START_MONTH=11, END_YEAR=2020, END_MONTH=2)
    assert ns['PERIOD_MONTHS'] == [(2019, 11), (2019, 12), (2020, 1), (2020, 2)]
    assert ns['PERIOD_DAYS'] == 30 + 31 + 31 + 29
    assert (ns['PERIOD_START'], ns['PERIOD_END']) == ('2019-11', '2020-02')


def test_single_month():
    ns = run_settings(START_YEAR=2021, START_MONTH=7, END_YEAR=2021, END_MONTH=7)
    assert ns['PERIOD_MONTHS'] == [(2021, 7)] and ns['PERIOD_DAYS'] == 31


@pytest.mark.parametrize('values', [dict(START_MONTH=13), dict(END_MONTH=0), dict(START_MONTH='3'),
                                    dict(START_YEAR=2021, END_YEAR=2020), dict(START_MONTH=6, END_MONTH=5)])
def test_invalid_period_stops_with_a_message(values):
    with pytest.raises(SystemExit) as e:
        run_settings(**values)
    assert str(e.value).startswith('❌')


@pytest.mark.parametrize('mean_lat,test_month_setting,period,expected', [
    (-27.0, 'auto', [(2020, m) for m in range(1, 13)], (2020, 1)),       # south: January
    (45.0, 'auto', [(2020, m) for m in range(1, 13)], (2020, 7)),        # north: July
    (-27.0, 'auto', [(2021, 6), (2021, 7), (2021, 8)], (2021, 6)),       # no January in the period: its first month
    (-27.0, 3, [(2019, 11), (2019, 12), (2020, 1), (2020, 2), (2020, 3)], (2020, 3)),
])
def test_test_month_is_inside_the_period(mean_lat, test_month_setting, period, expected):
    src = next(''.join(c['source']) for c in cells() if c['cell_type'] == 'code' and 'warm_month' in ''.join(c['source']))
    lines = src.splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.startswith('warm_month'))
    i1 = next(i for i, l in enumerate(lines) if l.startswith('test_year, test_month'))
    ns = dict(mean_lat=mean_lat, TEST_MONTH=test_month_setting, PERIOD_MONTHS=period)
    exec('\n'.join(lines[i0:i1 + 1]), ns)
    assert (ns['test_year'], ns['test_month']) == expected
