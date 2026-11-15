# -*- coding: utf-8 -*-
"""trend_robustness.py - robustness of the annual trends of the nine compound groups (POST HOC).

Tests whether the significance of the trends depends on the two years whose LAI was corrected with VIIRS (2024-2025), on
the definition of the multiple-testing family and on the way serial correlation is handled. Nothing here changes the inventory or the published trends; the script
recomputes the same series four ways.

Input  : notebook13_meso2024/table2_annual_totals_SC_groups.csv (annual state totals by group, t/yr)
Output : reproduced/final_review/trend_robustness.csv and trend_robustness.log (published copy in results/v10/final_review)

For each group and each period (2003-2025, the published period; 2003-2023, before the LAI correction applies):
  sen_pct_yr, sen_lo95, sen_hi95  Theil-Sen slope in % of the period mean, with the 95 % interval of scipy (independence)
  mk_p                            Mann-Kendall p value (Mann 1945; Kendall 1975), as published for 2003-2025
  tfpw_p                          Mann-Kendall p value after trend-free pre-whitening (Yue et al. 2002), which removes the
                                  lag-1 autocorrelation of the detrended series without a significance pre-test
  bh_q, by_q                      false-discovery-rate adjusted p values over the nine groups: Benjamini and Hochberg (1995),
                                  valid under independence or positive dependence, and Benjamini and Yekutieli (2001), valid
                                  under any dependence
  bh_q_8                          Benjamini-Hochberg over the eight groups without the total BVOC (the total is a sum of the
                                  others)
A moving-block bootstrap of the Theil-Sen slope is not used: with n = 23 and negative lag-1
autocorrelation its intervals were narrower than the published Theil-Sen intervals and excluded zero for groups whose
Mann-Kendall p exceeds 0.05, so it understated the uncertainty and is not reported.
Packages: pymannkendall 1.4.3 (original_test, trend_free_pre_whitening_modification_test), scipy, statsmodels.
"""
import os, sys
import numpy as np
import pandas as pd
import pymannkendall as pmk
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import _paths as P

GROUPS = ['BVOC', 'ISOP', 'MT', 'SQT', 'OVOC', 'LVOC', 'MBO', 'NO', 'CO']
PERIODS = [(2003, 2025), (2003, 2023)]
OUT = P.out('final_review')


def main():
    src = os.path.join(P.NB13, 'table2_annual_totals_SC_groups.csv')
    t = pd.read_csv(src).set_index('year')
    rows, log = [], [f'input: {os.path.relpath(src, P.REPO)}']
    for y0, y1 in PERIODS:
        d = t.loc[y0:y1]; part = []
        for g in GROUPS:
            y = d[g].to_numpy(float); x = np.arange(len(y)); m = y.mean()
            sl, ic, lo, hi = stats.theilslopes(y, x, 0.95)
            o = pmk.original_test(y); tf = pmk.trend_free_pre_whitening_modification_test(y)
            part.append(dict(period=f'{y0}-{y1}', group=g, n=len(y), sen_pct_yr=100 * sl / m, sen_lo95=100 * lo / m,
                             sen_hi95=100 * hi / m,
                             mk_p=o.p, tfpw_p=tf.p))
        df = pd.DataFrame(part)
        df['bh_q'] = multipletests(df.mk_p, method='fdr_bh')[1]
        df['by_q'] = multipletests(df.mk_p, method='fdr_by')[1]
        sub = df.group != 'BVOC'
        df.loc[sub, 'bh_q_8'] = multipletests(df.loc[sub, 'mk_p'], method='fdr_bh')[1]
        rows.append(df)
        log.append(f'{y0}-{y1}: ' + '; '.join(f'{r.group} {r.sen_pct_yr:+.2f} %/yr MK p {r.mk_p:.4f} TFPW p {r.tfpw_p:.4f} '
                                              f'BH q {r.bh_q:.4f} BY q {r.by_q:.4f}' for r in df.itertuples()))
    res = pd.concat(rows, ignore_index=True)
    res.to_csv(os.path.join(OUT, 'trend_robustness.csv'), index=False, float_format='%.6g')
    with open(os.path.join(OUT, 'trend_robustness.log'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(log) + '\n')
    print('\n'.join(log)); print('written:', os.path.join(OUT, 'trend_robustness.csv'))


if __name__ == '__main__':
    main()
