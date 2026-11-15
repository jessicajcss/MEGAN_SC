#!/usr/bin/env python3
"""lai_variant_view.py - the production run with the MODIS LAI in every year (the "MODIS view").

`data/run_v10` holds 2024 and 2025 with the VIIRS-corrected LAI (notebook 17), and `data/run_v10_modis_lai_2024_2025`
holds the same two years computed with the MODIS LAI. Two uses need the MODIS record of every year:
  * the step `laicorr` of reproduce.py (analysis/v10/review/lai_viirs_correction.py), which compares the MODIS LAI of
    2024-2025 with VIIRS and builds the correction;
  * the MODIS control (REPRODUCE.md, section B), which recomputes results/v10_modis_lai_2024_2025.
The view is a folder with the layout of run_v10 in which every year is a symbolic link: 2024 and 2025 to the archive,
the other years to run_v10 (copies with --copy, or where the file system has no symbolic links).

    python tools/lai_variant_view.py                       # data/run_v10 + data/run_v10_modis_lai_2024_2025
    python tools/lai_variant_view.py --run R --archive A --out V [--copy]

When run_v10 still holds the MODIS years (provenance variant 'base' in 2024 and 2025), run_v10 itself is the MODIS
record and no view is needed: `modis_root()` then returns run_v10.
"""
import argparse, json, os, shutil, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
YEARS_REPLACED = (2024, 2025)


def variant(run, year):
    """Variant recorded in <run>/<year>/provenance_v10.json ('base', 'lai_viirs', ...), or None."""
    p = os.path.join(run, str(year), 'provenance_v10.json')
    if not os.path.exists(p):
        return None
    with open(p, encoding='utf-8') as fh:
        return json.load(fh).get('variant')


def _link(src, dst, copy):
    if os.path.lexists(dst):
        if os.path.islink(dst) and os.path.realpath(dst) == os.path.realpath(src):
            return 'kept'
        if os.path.islink(dst) or os.path.isfile(dst):
            os.remove(dst)
        else:
            shutil.rmtree(dst)
    if not copy:
        try:
            os.symlink(os.path.abspath(src), dst, target_is_directory=os.path.isdir(src))
            return 'link'
        except (OSError, NotImplementedError):
            pass
    if os.path.isdir(src):
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)
    return 'copy'


def build_view(run, archive, out, copy=False, years=YEARS_REPLACED):
    """Create (or refresh) the MODIS view `out`; returns {entry: 'link'|'copy'|'kept'}."""
    for y in years:
        a = os.path.join(archive, str(y))
        if not os.path.isdir(a):
            raise FileNotFoundError(f'{a} not found: the MODIS-LAI year {y} is needed for the MODIS view')
        v = variant(archive, y)
        if v not in (None, 'base'):
            raise ValueError(f'{a}: provenance variant {v!r}, expected the MODIS-LAI run (variant base)')
    os.makedirs(out, exist_ok=True)
    done = {}
    for e in sorted(os.listdir(run)):
        src = os.path.join(archive, e) if e in {str(y) for y in years} else os.path.join(run, e)
        done[e] = _link(src, os.path.join(out, e), copy)
    return done


def modis_root(run, archive, out, copy=False):
    """Folder holding the MODIS-LAI record of every year: run itself if its 2024-2025 are the MODIS years, else the view."""
    if all(variant(run, y) in (None, 'base') for y in YEARS_REPLACED):
        return run
    build_view(run, archive, out, copy)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--run', default=os.path.join(REPO, 'data', 'run_v10'))
    ap.add_argument('--archive', default=os.path.join(REPO, 'data', 'run_v10_modis_lai_2024_2025'))
    ap.add_argument('--out', default=os.path.join(REPO, 'data', 'run_v10_modis_view'))
    ap.add_argument('--copy', action='store_true', help='copy the year folders instead of linking them')
    a = ap.parse_args()
    r = modis_root(a.run, a.archive, a.out, a.copy)
    if r == a.run:
        print(f'{a.run}: 2024 and 2025 are the MODIS-LAI years; no view needed')
    else:
        print(f'MODIS view ready: {r}')
    sys.exit(0)
