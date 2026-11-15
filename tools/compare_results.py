#!/usr/bin/env python3
"""compare_results.py <reproduced_dir> <published_dir> [name=folder ...] - compare every reproduced file with its published copy.

For each file under <reproduced_dir> with the same relative path under <published_dir> (results/v10):
  * CSV: same columns and rows; numbers equal within a relative tolerance of 1e-6 (absolute 1e-9), text equal after
    normalising NumPy 2 reprs (np.float64(x) -> x). Status: IDENTICAL (bytes), EQUAL (within tolerance) or DIFFERENT.
  * PNG: pixel-identical, or the share of differing pixels (figures depend on fonts and library versions, so a PNG
    difference is reported but does not count as a failure).
  * GPKG: same attribute table, and every geometry equal within 1e-9 of its area (symmetric difference).
  * other text files (.log): reported as present; logs carry timings and paths and are not compared.
A `name=folder` argument maps the reproduced sub-folder `name` to another published folder (reproduce.py passes
boundaries=data/boundaries, whose files are inputs of the repository rather than results).
Files without a published copy are listed as NEW. The table goes to stdout; the exit code is 1 if any CSV differs.
"""
import os, re, sys
import numpy as np, pandas as pd

TOL_R, TOL_A = 1e-6, 1e-9
NP = re.compile(r'np\.(?:float|int)\d*\((.*?)\)')


def norm_text(x):
    return NP.sub(r'\1', str(x)).strip()


def cmp_csv(a, b):
    if open(a, 'rb').read() == open(b, 'rb').read():
        return 'IDENTICAL', 0.0
    A = pd.read_csv(a, dtype=str, keep_default_na=False); B = pd.read_csv(b, dtype=str, keep_default_na=False)
    if list(A.columns) != list(B.columns) or A.shape != B.shape:
        return 'DIFFERENT', f'shape/columns {A.shape} {list(A.columns)[:4]} vs {B.shape} {list(B.columns)[:4]}'
    worst = 0.0
    for c in A.columns:
        for x, y in zip(A[c], B[c]):
            x, y = norm_text(x), norm_text(y)
            if x == y:
                continue
            try:
                fx, fy = float(x), float(y)
            except ValueError:
                return 'DIFFERENT', f'column {c}: {x[:60]!r} vs {y[:60]!r}'
            if np.isnan(fx) and np.isnan(fy):
                continue
            d = abs(fx - fy); rel = d / max(abs(fy), 1e-300)
            if d > TOL_A and rel > TOL_R:
                return 'DIFFERENT', f'column {c}: {fx!r} vs {fy!r} (rel {rel:.1e})'
            worst = max(worst, rel if d > TOL_A else 0.0)
    return 'EQUAL', worst


def cmp_png(a, b):
    if open(a, 'rb').read() == open(b, 'rb').read():
        return 'IDENTICAL', 0.0
    try:
        from PIL import Image
        A = np.asarray(Image.open(a).convert('RGB')); B = np.asarray(Image.open(b).convert('RGB'))
    except Exception as e:
        return 'PNG?', str(e)
    if A.shape != B.shape:
        return 'PNG SIZE', f'{A.shape} vs {B.shape}'
    frac = float((np.abs(A.astype(int) - B.astype(int)).sum(axis=2) > 0).mean())
    return ('PIXEL-IDENTICAL' if frac == 0 else 'PNG DIFF'), frac


def cmp_gpkg(a, b):
    import geopandas as gpd
    A = gpd.read_file(a); B = gpd.read_file(b)
    ca = [c for c in A.columns if c != 'geometry']; cb = [c for c in B.columns if c != 'geometry']
    if ca != cb or len(A) != len(B):
        return 'DIFFERENT', f'attributes {ca} vs {cb}, rows {len(A)} vs {len(B)}'
    fa = A[ca].astype(str).reset_index(drop=True); fb = B[cb].astype(str).reset_index(drop=True)
    if not fa.equals(fb):
        return 'DIFFERENT', 'attribute values'
    worst = max(float(ga.symmetric_difference(gb).area / max(gb.area, 1e-300)) for ga, gb in zip(A.geometry, B.geometry))
    return ('EQUAL' if worst < 1e-9 else 'DIFFERENT'), worst


def main(rep, pub, *maps):
    MAP = dict(m.split('=', 1) for m in maps)
    rows = []
    for root, _, files in os.walk(rep):
        rel_root = os.path.relpath(root, rep)
        if rel_root.split(os.sep)[0] in ('_work', 'logs') or os.sep + 'cache_' in os.sep + rel_root or rel_root.startswith('notebook14' + os.sep + 'cache_'):
            continue
        for f in sorted(files):
            rel = os.path.normpath(os.path.join(rel_root, f)); a = os.path.join(rep, rel); b = os.path.join(pub, rel)
            top = rel.split(os.sep)[0]
            if top in MAP:
                b = os.path.join(MAP[top], os.path.relpath(rel, top))
            ext = os.path.splitext(f)[1].lower()
            if not os.path.exists(b):
                rows.append((rel, 'NEW', '')); continue
            if ext == '.csv':
                st, info = cmp_csv(a, b)
            elif ext == '.png':
                st, info = cmp_png(a, b)
            elif ext == '.gpkg':
                st, info = cmp_gpkg(a, b)
            else:
                st, info = 'PRESENT', ''
            rows.append((rel, st, info))
    w = max(len(r[0]) for r in rows) if rows else 10
    for rel, st, info in sorted(rows):
        print(f'{rel:{w}s}  {st:15s}  {info if isinstance(info, str) else f"{info:.2e}"}')
    csv = [r for r in rows if r[0].endswith('.csv')]
    bad = [r for r in csv if r[1] == 'DIFFERENT']
    print(f'\nCSV files: {len(csv)} compared, {sum(r[1] == "IDENTICAL" for r in csv)} identical, '
          f'{sum(r[1] == "EQUAL" for r in csv)} equal within tolerance, {len(bad)} different, '
          f'{sum(r[1] == "NEW" for r in rows if r[0].endswith(".csv"))} without a published copy')
    png = [r for r in rows if r[0].endswith('.png') and r[1] != 'NEW']
    print(f'PNG files: {len(png)} compared, {sum(r[1] in ("IDENTICAL", "PIXEL-IDENTICAL") for r in png)} pixel-identical')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
