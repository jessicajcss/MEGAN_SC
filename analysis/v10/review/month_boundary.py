# -*- coding: utf-8 -*-
"""month_boundary.py - effect of counting months in UTC instead of local time (POST HOC).

The run integrates every month from 00 UTC of its first day to 00 UTC of the first day of the next month, as the
ERA5-Land data are in UTC. In local time at UTC-3, a month would run from 03 UTC to 03 UTC. Each 3-h step h (00, 03,
..., 21 UTC) stands for [h, h+3) UTC in the monthly integration, so the local month equals the UTC month minus the 00 UTC
step of its first day plus the 00 UTC step of the first day of the next month. The 00 UTC step falls at about
20:24-20:47 local mean solar time across the state (cell centres 53.89-48.31 deg W), at night.

Input  : <MEGAN_SC_RUN>/<year>/hourly/diurnal_<year>_<mm>.csv, the mean over the days of the month of each 3-h step
         (domain sum of the class flux, nmol m-2 s-1 summed over cells), written by notebook 10
Output : reproduced/final_review/month_boundary_effect.csv (one row per month and group) and month_boundary_effect.log
         (published copy in results/v10/final_review)

The first-day values are not stored, so the 00 UTC step of the first day of each month is approximated by the mean
00 UTC step of that month. The groups are summed in moles; all VOC classes means the 17 classes other than NO and CO.
"""
import calendar, glob, os, sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P

OUT = P.out('final_review')


def load(run_dir):
    d = {}
    for f in sorted(glob.glob(os.path.join(run_dir, '*', 'hourly', 'diurnal_*.csv'))):
        tag = os.path.basename(f)[len('diurnal_'):-len('.csv')]
        y, m = tag.split('_')[:2]
        if not (y.isdigit() and m.isdigit()) or len(tag.split('_')) != 2:   # skip test runs (_testNd suffix)
            continue
        d[(int(y), int(m))] = pd.read_csv(f, index_col=0)
    if not d:
        sys.exit(f'no diurnal_<year>_<mm>.csv files under {run_dir}')
    return d


def main():
    D = load(P.RUN)
    cls = list(next(iter(D.values())).index)
    groups = {'isoprene': ['ISOP'], 'monoterpenes': [c for c in cls if c.startswith('MT_')], 'soil NO': ['NO'],
              'all VOC classes': [c for c in cls if c not in ('NO', 'CO')]}
    rows = []
    for (y, m), d in sorted(D.items()):
        nxt = (y + (m == 12), m % 12 + 1)
        if nxt not in D:
            continue
        n = calendar.monthrange(y, m)[1]
        for g, cl in groups.items():
            day = d.loc[cl].values.sum()
            s00, s00n = d.loc[cl, '00h'].sum(), D[nxt].loc[cl, '00h'].sum()
            rows.append(dict(year=y, month=m, group=g, share_00utc_step_of_day=s00 / day if day else 0.0,
                             local_minus_utc_pct=100.0 * (s00n - s00) / (n * day) if day else 0.0))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, 'month_boundary_effect.csv'), index=False, float_format='%.6g')
    years = sorted({y for y, _ in D})
    A = []
    for y in years:
        if (y + 1, 1) not in D or not all((y, m) in D for m in range(1, 13)):
            continue
        for g, cl in groups.items():
            tot = sum(calendar.monthrange(y, m)[1] * D[(y, m)].loc[cl].values.sum() for m in range(1, 13))
            A.append(dict(year=y, group=g, local_minus_utc_pct=100.0 * (D[(y + 1, 1)].loc[cl, '00h'].sum()
                                                                          - D[(y, 1)].loc[cl, '00h'].sum()) / tot))
    A = pd.DataFrame(A)
    S = R.groupby('group').agg(months=('month', 'size'), share_00utc_mean=('share_00utc_step_of_day', 'mean'),
                               share_00utc_max=('share_00utc_step_of_day', 'max'),
                               abs_effect_p95_pct=('local_minus_utc_pct', lambda x: x.abs().quantile(0.95)),
                               abs_effect_max_pct=('local_minus_utc_pct', lambda x: x.abs().max()))
    S['annual_abs_effect_max_pct'] = A.groupby('group').local_minus_utc_pct.agg(lambda x: x.abs().max())
    with open(os.path.join(OUT, 'month_boundary_effect.log'), 'w', encoding='utf-8') as fh:
        fh.write(f'months counted from 03 UTC (local midnight at UTC-3) instead of 00 UTC; run folder {os.path.basename(P.RUN)}\n')
        fh.write(f'monthly effects: {len(R) // len(groups)} months with a following month; annual: {len(A) // len(groups)} years\n')
        fh.write(S.to_string(float_format=lambda v: f'{v:.4f}') + '\n')
    print(open(os.path.join(OUT, 'month_boundary_effect.log'), encoding='utf-8').read())


if __name__ == '__main__':
    main()
