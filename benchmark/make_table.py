"""make_table.py - prints the Markdown table of VALIDATION.md from results_colab/benchmark_summary.csv (no number is typed by hand).
Usage: python benchmark/make_table.py > /tmp/table.md
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(HERE, 'results_colab', 'benchmark_summary.csv'))
case = {'se_us_bench_N_N.csv': 'SE_US, original port', 'se_us_bench_patched_N_N.csv': 'SE_US, patched',
        'se_us_bench_patched_Y_Y.csv': 'SE_US, patched, GAMSM and GAMBD on', 'se_us_bench_soil_unpatched.csv': 'SE_US soil, original port',
        'se_us_bench_soil_patched.csv': 'SE_US soil, patched', 'sc_bench_sc_m01.csv': 'SC grid, 15 Jan 2019',
        'sc_bench_sc_m07.csv': 'SC grid, 15 Jul 2019', 'se_us_bench_speciation.csv': 'SE_US speciation'}
def f(x):
    return '' if pd.isna(x) else f'{x:.1e}'
print('| Case | Stage | max rel. diff. | p99 rel. diff. | domain-sum rel. diff. | day-tonnes rel. diff. | Verdict |')
print('|---|---|---|---|---|---|---|')
for r in d.itertuples():
    print(f"| {case.get(r.file, r.file)} | {r.stage} | {f(r.max_rel)} | {f(r.p99_rel)} | {f(r.sum_dev)} | {f(r.day_tonnes_dev)} | {r.verdict} |")
print(f"\nfixes.py sha256[:16] = {d.fixes_sha256.iloc[0]}; rows: {len(d)}; "
      + '; '.join(f'{k}: {v}' for k, v in d.verdict.value_counts().items()))
