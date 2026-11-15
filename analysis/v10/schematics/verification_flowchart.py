#!/usr/bin/env python3
"""verification_flowchart.py - the verification and observational evaluation of the inventory (article and report Figure 2): one
card per step with the question it answers, how it was tested, its result and what the result means, on a numbered
timeline (verification, steps 1-4; evaluation, step 5).

Every number in the figure is read from the result files (reproduced copy if reproduce.py has written it, else the
published one in results/v10):
  benchmark/results_colab/benchmark_summary.csv, sc_bench_sc_m01.csv, sc_bench_sc_m07.csv   code equivalence (notebook 12)
  diagnostics/check_run_integrity.log                                                        steps per day, annual = sum
  notebook11/verify_v10_report.csv, diagnostics/c3/c3_impact_on_totals.log                  checks C0-C5
  notebook15/s1/s1_groups.csv, notebook15/S3_fixlog_annual.csv                               configuration tests
  sensitivity_2019/full/variants_summary.csv (notebook 20), ijest_review/plantation_ef_sensitivity_2019.csv
                                                                                             sensitivity tests of 2019
  notebook14/tier{1,2,3,4}_spatial|seasonal|anomalies|slope_cf0.3_sza70_fire120.csv          TROPOMI tests
Outputs: <MEGAN_SC_OUT>/figures/{pt,en}/fig_verification.png (published copies in results/v10/figures).
Step 1 states that the emission-factor grid entered both codes and was not part of the equivalence test; step 4 reports
the one-at-a-time sensitivity runs of 2019 (notebook 20) and the planted-forest emission factors next to the time-step and
soil-moisture tests."""
import math, os, re, sys, textwrap
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..')); import _paths as P
TAG = P.TAG


def pick(*rel):
    f = os.path.join(P.OUT, *rel); return f if os.path.exists(f) else os.path.join(P.RESULTS, *rel)


def up1(x):
    """x rounded up to one significant figure, as mantissa and exponent."""
    e = math.floor(math.log10(x)); m = math.ceil(x / 10 ** e)
    if m == 10: m, e = 1, e + 1
    return m, e


def numbers():
    N = {}
    bench = os.path.join(P.REPO, 'benchmark', 'results_colab')
    s = pd.read_csv(os.path.join(bench, 'benchmark_summary.csv'))
    p99 = list(s[(s.verdict == 'PASS') & s.stage.str.contains('end-to-end')].p99_rel)
    for m in ('m01', 'm07'):
        b = pd.read_csv(os.path.join(bench, f'sc_bench_sc_{m}.csv')); p99 += list(b[b.stage.str.contains('patched engine')].p99_rel)
    N['p99'] = up1(max(p99))
    log = open(pick('diagnostics', 'check_run_integrity.log'), encoding='utf-8').read()
    N['months'] = int(re.search(r'\((\d+) months checked\)', log).group(1)); assert 'none' in log.split('(a)')[1].split('\n')[0]
    N['annual_dev'] = up1(float(re.search(r'sum of monthly\| over years and classes: ([0-9.eE+-]+)', log).group(1)))
    v = pd.read_csv(pick('notebook11', 'verify_v10_report.csv'))
    N['n_checks'], N['n_pass'] = len(v), int((v.result == 'PASS').sum())
    f = v[v.result != 'PASS'].iloc[0]; N['fail'] = (f.check.split()[0], int(f.year))
    c3 = open(pick('diagnostics', 'c3', 'c3_impact_on_totals.log'), encoding='utf-8').read()
    m = re.search(rf'{N["fail"][1]} MEOH (\d+) cells;.*?change if neighbour-median activity ([+-][0-9.]+)%', c3)
    N['fail_cells'], N['fail_effect'] = int(m.group(1)), float(m.group(2))
    g = pd.read_csv(pick('notebook15', 's1', 's1_groups.csv')); b = g[g.group == 'BVOC']
    N['s1'] = 100 * (b.t_3h.sum() / b.t_1h.sum() - 1)
    a = pd.read_csv(pick('notebook15', 'S3_fixlog_annual.csv'))
    N['gamsm'] = 100 * a.bvoc_reduction_by_gamsm_approx.mean()
    vs = pd.read_csv(pick('sensitivity_2019', 'full', 'variants_summary.csv'))
    vs = vs[(vs.region == 'Santa Catarina') & (vs.group == 'BVOC') & (vs.variant != 'base')]
    ch = list(vs.change_vs_base_pct) + [N['s1'], -N['gamsm']]                 # gamsm is a reduction (positive number)
    N['sens'] = (min(ch), max(ch))
    pl = pd.read_csv(pick('ijest_review', 'plantation_ef_sensitivity_2019.csv'))
    pl = pl[(pl.group == 'ISOP') & pl.case.str.startswith('plantations as')].change_pct
    N['plant'] = (pl.min(), pl.max())
    t1 = pd.read_csv(pick('notebook14', f'tier1_spatial_{TAG}.csv')); t1 = t1[(t1.block_deg.round(2) == 0.1) & (t1.smoothing_km == 0)].iloc[0]
    t2 = pd.read_csv(pick('notebook14', f'tier2_seasonal_{TAG}.csv'), index_col=0).iloc[:, 0]
    t3 = pd.read_csv(pick('notebook14', f'tier3_anomalies_{TAG}.csv'), index_col=0).iloc[0]
    t4 = pd.read_csv(pick('notebook14', f'tier4_slope_{TAG}.csv')); t4 = t4[t4.smoothing_km == 0].iloc[0]
    N['t1'], N['t2'], N['t3'], N['t4'] = (t1.pearson_r, t1.r_lo95, t1.r_hi95), (t2['r_climatology'], t2['p']), (t3.r, t3.p_neff), t4.r
    return N


def fmt(v, d, lang):
    s = f'{v:.{d}f}'; s = s.replace('.', ',') if lang == 'pt' else s
    return s.replace('-', '−')


def sci(me):
    m, e = me; return f'{m}\u00a0×\u00a010{str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))}'


def texts(N, lang):
    f = lambda v, d: fmt(v, d, lang)
    N = dict(N, bound=math.ceil(10 * max(abs(N['s1']), abs(N['gamsm'])) - 1e-9) / 10)      # 0.1 % steps
    p3 = (lambda sep: f'p < 0{sep}001' if N['t3'][1] < 0.001 else f'p = {f(N["t3"][1], 3)}')(',' if lang == 'pt' else '.')
    if lang == 'pt':
        H = dict(how='COMO FOI TESTADO', res='RESULTADO', mean='O QUE SIGNIFICA', ph=('VERIFICAÇÃO', 'AVALIAÇÃO'),
                 phs=('o inventário foi calculado como pretendido?', 'o resultado concorda com observações independentes?'))
        S = [
         dict(t='Equivalência do código', q='A implementação em Python reproduz o código Fortran do MEGANv3.21?',
              how='Mesmas entradas nos dois códigos: caso de teste oficial (sudeste dos EUA, 2 de maio de 2018) e dois dias na grade de SC.',
              key=f'< {sci(N["p99"])}', det='percentil 99 da diferença relativa por célula e hora',
              mean='Com as mesmas entradas, o motor reproduz as emissões do Fortran. A grade de fatores de emissão entrou nos dois códigos e não fez parte do teste.', st='ok', chip='aprovado'),
         dict(t='Integridade da execução', q='Todos os passos de tempo foram calculados e somados?',
              how='Cada mês deve ter 8 passos por dia, e o arquivo anual deve ser a soma dos 12 arquivos mensais.',
              key=f'{N["months"]} de {N["months"]} meses', det=f'completos; anual menos a soma dos mensais < {sci(N["annual_dev"])}',
              mean='Nenhum mês incompleto ou duplicado entra nos totais.', st='ok', chip='aprovado'),
         dict(t='Testes automáticos C0 a C5', q='As saídas de cada ano exibem alguma assinatura de falha conhecida?',
              how='C1: uma célula com árvores que emite COVB sem emitir isopreno indica cálculo sem luz, isto é, sem meteorologia válida.',
              key=f'{N["n_pass"]} de {N["n_checks"]} testes', det=f'a falha ({N["fail"][0]}, metanol, {N["fail"][1]}, {N["fail_cells"]} células) altera o metanol estadual em +{f(N["fail_effect"], 3)}%',
              mean='Os totais não contêm células sem meteorologia, valores extremos nem meses faltantes.', st='ok', chip='aprovado, 1 falha local'),
         dict(t='Testes de sensibilidade', q='Quanto as configurações e os fatores de emissão das florestas plantadas mudam os totais?',
              how='Uma escolha alterada por vez: passo de tempo (28 dias), fator de umidade do solo (2003–2025) e, em 2019, fração da PAR, temperatura no relevo, tipo de dossel, ponto de murcha e fatores de Pinus ou Eucalyptus nas florestas plantadas.',
              key=f'{f(N["sens"][0], 1)}% a +{f(N["sens"][1], 1)}%', det=f'no COVB, pelas configurações; {f(N["plant"][0], 0)}% a +{f(N["plant"][1], 0)}% no isopreno, pelos fatores das florestas plantadas',
              mean='As configurações alteram os totais em poucos por cento; os fatores de emissão pesam mais. Nenhuma observação restringe a magnitude.', st='part', chip=f'configurações < {f(math.ceil(max(abs(N["sens"][0]), abs(N["sens"][1]))), 0)}%'),
         dict(t='Avaliação com o TROPOMI (HCHO)', q='As colunas de formaldeído observadas por satélite acompanham o isopreno modelado?',
              how='Quatro testes definidos antes da execução: padrão espacial de verão, ciclo sazonal, anomalias mensais e inclinação coluna-emissão.',
              tests=[('Anomalias', f'r = {f(N["t3"][0], 2)} ({p3})', True),
                     ('Padrão espacial', f'r = {f(N["t1"][0], 2)} (IC 95% {f(N["t1"][1], 2)} a {f(N["t1"][2], 2)})', False),
                     ('Ciclo sazonal', f'r = {f(N["t2"][0], 2)}', False), ('Inclinação', f'r = {f(N["t4"], 2)}', False)],
              mean='As anomalias mensais concordam, o que testa a resposta às variações do tempo. O padrão espacial, o ciclo sazonal e a magnitude ficam sem confirmação observacional.', st='part', chip='consistente em 1 de 4 testes'),
        ]
    else:
        H = dict(how='HOW IT WAS TESTED', res='RESULT', mean='WHAT IT MEANS', ph=('VERIFICATION', 'EVALUATION'),
                 phs=('was the inventory computed as intended?', 'does the result agree with independent observations?'))
        S = [
         dict(t='Code equivalence', q='Does the Python implementation reproduce the Fortran MEGANv3.21 code?',
              how='Same inputs in both codes: the official test case (south-eastern USA, 2 May 2018) and two days on the SC grid.',
              key=f'< {sci(N["p99"])}', det='99th percentile of the relative difference per cell and hour',
              mean='With the same inputs the engine reproduces the Fortran emissions. The emission-factor grid entered both codes and was not part of the test.', st='ok', chip='passed'),
         dict(t='Run integrity', q='Were all time steps computed and summed?',
              how='Each month must hold 8 steps per day, and the annual file must be the sum of the 12 monthly files.',
              key=f'{N["months"]} of {N["months"]} months', det=f'complete; annual minus the sum of the monthly files < {sci(N["annual_dev"])}',
              mean='No incomplete or duplicated month enters the totals.', st='ok', chip='passed'),
         dict(t='Automatic checks C0 to C5', q='Do the outputs of each year show a known failure signature?',
              how='C1: a cell with trees that emits BVOC but no isoprene was computed without light, that is, without valid meteorology.',
              key=f'{N["n_pass"]} of {N["n_checks"]} checks', det=f'the failure ({N["fail"][0]}, methanol, {N["fail"][1]}, {N["fail_cells"]} cells) changes state methanol by +{f(N["fail_effect"], 3)}%',
              mean='The totals contain no cells without meteorology, no extreme values and no missing months.', st='ok', chip='passed, 1 local failure'),
         dict(t='Sensitivity tests', q='How much do the settings and the planted-forest emission factors change the totals?',
              how='One choice changed at a time: time step (28 days), soil-moisture factor (2003–2025) and, for 2019, PAR fraction, terrain temperature, canopy type, wilting point, and Pinus or Eucalyptus factors for planted forests.',
              key=f'{f(N["sens"][0], 1)}% to +{f(N["sens"][1], 1)}%', det=f'in BVOC, from the settings; {f(N["plant"][0], 0)}% to +{f(N["plant"][1], 0)}% in isoprene, from the planted-forest factors',
              mean='The settings change the totals by a few percent; the emission factors matter more. No observation constrains the magnitude.', st='part', chip=f'settings < {f(math.ceil(max(abs(N["sens"][0]), abs(N["sens"][1]))), 0)}%'),
         dict(t='Evaluation with TROPOMI (HCHO)', q='Do satellite formaldehyde columns follow the modelled isoprene?',
              how='Four tests defined before the evaluation was run: summer spatial pattern, seasonal cycle, monthly anomalies and column-emission slope.',
              tests=[('Anomalies', f'r = {f(N["t3"][0], 2)} ({p3})', True),
                     ('Spatial pattern', f'r = {f(N["t1"][0], 2)} (95% CI {f(N["t1"][1], 2)} to {f(N["t1"][2], 2)})', False),
                     ('Seasonal cycle', f'r = {f(N["t2"][0], 2)}', False), ('Slope', f'r = {f(N["t4"], 2)}', False)],
              mean='Monthly anomalies agree, which tests the response to weather. The spatial pattern, the seasonal cycle and the magnitude remain without observational confirmation.', st='part', chip='consistent in 1 of 4 tests'),
        ]
    return H, S


C_VER, C_EVA = '#0f6e67', '#5a46b0'          # teal (verification), violet (evaluation)
F_VER, F_EVA = '#eef6f5', '#f2f0fa'          # card fills
INK, GREY, LINE = '#1f2328', '#5b6470', '#d0d7de'
OKC, NCC = '#1a7f37', '#8c959f'              # passed / not confirmed


def draw(N, lang):
    H, S = texts(N, lang)
    plt.rcParams.update({'font.family': 'DejaVu Sans'})
    W_IN, U = 7.2, 10.0                                             # figure width (in); data units per inch on both axes
    xs = dict(card=6.6, how=8.0, res=29.4, mean=50.2, end=W_IN * U - 0.5, badge=3.3)
    cw = dict(how=40, res=37, mean=39)                              # characters per line of 6.3 pt text in each column
    lh = 6.3 * 1.28 / 72 * U                                        # height of one line of 6.3 pt text (data units)
    fig = plt.figure(figsize=(W_IN, 12.0), dpi=300); ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_IN * U); ax.set_ylim(120.0, 0); ax.axis('off')
    R = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    ext = lambda t: t.get_window_extent(renderer=R).transformed(inv)
    for k in ('how', 'res', 'mean'):
        ax.text(xs[k], 2.4, H[k], fontsize=6.2, color=GREY, fontweight='bold', va='center', ha='left')
    y = 4.2; centres = []
    for i, s in enumerate(S):
        ph = 0 if i < len(S) - 1 else 1; col = (C_VER, C_EVA)[ph]; fill = (F_VER, F_EVA)[ph]
        if i in (0, len(S) - 1):                                    # phase heading
            if ph: y += 0.6
            t0 = ax.text(xs['card'] + 0.4, y + 2.0, H['ph'][ph], fontsize=7.0, fontweight='bold', color=col, va='center')
            ax.text(ext(t0).x1 + 0.8, y + 2.0, '· ' + H['phs'][ph], fontsize=6.6, color=GREY, va='center', style='italic')
            y += 4.0
        arts = []
        arts.append(ax.text(xs['how'], y + 2.1, s['t'], fontsize=8.0, fontweight='bold', color=INK, va='center'))
        ax.text(xs['end'] - 1.0, y + 2.1, ('✓ ' if s['st'] == 'ok' else '◐ ') + s['chip'], fontsize=6.4,
                color=OKC if s['st'] == 'ok' else col, ha='right', va='center', fontweight='bold')
        q = ax.text(xs['how'], y + 3.6, s['q'], fontsize=6.6, color=GREY, va='top', style='italic'); arts.append(q)
        yb = ext(q).y0 + 1.2 if ext(q).y0 > ext(q).y1 else ext(q).y1 + 1.2
        ax.plot([xs['how'], xs['end'] - 1.0], [yb - 0.6, yb - 0.6], color=LINE, lw=0.5, zorder=2)
        arts.append(ax.text(xs['how'], yb, textwrap.fill(s['how'], cw['how']), fontsize=6.3, color=INK, va='top', linespacing=1.28))
        arts.append(ax.text(xs['mean'], yb, textwrap.fill(s['mean'], cw['mean']), fontsize=6.3, color=INK, va='top', linespacing=1.28))
        if 'tests' in s:                                            # the four TROPOMI tests, one per line
            yy = yb
            for nm, val, ok in s['tests']:
                t0 = ax.text(xs['res'], yy, ('● ' if ok else '○ ') + nm + ':', fontsize=6.3, color=OKC if ok else NCC, va='top', fontweight='bold')
                x1 = ext(t0).x1 + 0.5; t1 = ax.text(x1, yy, val, fontsize=6.3, color=INK if ok else GREY, va='top'); arts += [t0, t1]
                if ext(t1).x1 > xs['mean'] - 1.0 and ' (' in val:          # a confidence interval goes to the next line
                    a_, b_ = val.split(' (', 1); t1.set_text(a_); yy += 1.15 * lh
                    arts.append(ax.text(x1, yy, '(' + b_, fontsize=6.3, color=INK if ok else GREY, va='top'))
                yy += 1.15 * lh
        else:
            k_ = ax.text(xs['res'], yb - 0.1, s['key'], fontsize=9.2, fontweight='bold', color=col, va='top'); arts.append(k_)
            arts.append(ax.text(xs['res'], max(ext(k_).y0, ext(k_).y1) + 0.6, textwrap.fill(s['det'], cw['res']), fontsize=6.3,
                                color=GREY, va='top', linespacing=1.28))
        bottom = max(max(ext(a).y0, ext(a).y1) for a in arts) + 1.1
        ax.add_patch(FancyBboxPatch((xs['card'], y), xs['end'] - xs['card'], bottom - y, boxstyle='round,pad=0,rounding_size=1.0',
                                    fc=fill, ec='none', zorder=0))
        cy = y + 2.1; centres.append((cy, col))
        ax.add_patch(Circle((xs['badge'], cy), 1.85, fc=col, ec='white', lw=1.2, zorder=4))
        ax.text(xs['badge'], cy + 0.05, str(i + 1), color='white', fontsize=8, fontweight='bold', ha='center', va='center', zorder=5)
        y = bottom + 1.2
    for (y0, c0), (y1, c1) in zip(centres[:-1], centres[1:]):      # timeline behind the badges
        ax.plot([xs['badge']] * 2, [y0 + 1.9, y1 - 1.9], color=c0 if c0 == c1 else '#b9b4d9', lw=1.4, zorder=3, solid_capstyle='butt')
    ax.set_ylim(y + 0.2, 0); fig.set_size_inches(W_IN, (y + 0.2) / U)
    out = P.out('figures', lang, 'fig_verification.png'); fig.savefig(out, dpi=300, facecolor='white'); plt.close(fig)
    print('written', out)


if __name__ == '__main__':
    N = numbers(); print(N)
    for lang in ('pt', 'en'):
        draw(N, lang)
