#!/usr/bin/env python3
"""anonymize_release.py - copy of the release archives without the Earth Engine project of the run (double-blind review).

Every model year of the run writes provenance_v10.json, whose field "ee_project" names the Google Cloud project that ran
Earth Engine. That name can identify the authors, so for double-blind review the field is set to "withheld" in the three
archives that hold these files (megan_sc_run_v10_2003-2025, megan_sc_run_v10_modis_lai_2024_2025, megan_sc_s1_test_2019).
Nothing else changes: every other member is copied with the same bytes after decompression (checked through its CRC-32
and size), and the two archives without provenance files are copied as they are.

    python tools/anonymize_release.py <dir_with_the_five_zips> <out_dir>

<out_dir> receives the five zips and SHA256SUMS_assets.txt (upload both to the release). The SHA-256 of the rewritten
provenance files are those of data/DATA_SHA256SUMS.txt; the script stops if one differs.
"""
import hashlib, os, re, shutil, sys, zipfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'tools'))
from data_release import ASSETS, SUMS, sha256   # noqa: E402

FIELD = re.compile(rb'("ee_project"\s*:\s*)"[^"]*"')
WITHHELD = b'"withheld"'


def anonymize_bytes(data):
    """The provenance file with the value of "ee_project" replaced; the rest of the bytes unchanged."""
    out, n = FIELD.subn(lambda m: m.group(1) + WITHHELD, data)
    if n != 1:
        raise ValueError(f'expected one "ee_project" field, found {n}')
    return out


def expected_sums():
    rows = {}
    for line in open(os.path.join(REPO, 'data', 'DATA_SHA256SUMS.txt'), encoding='utf-8'):
        if line.strip():
            h, rel = line.split(None, 1); rows[rel.strip()] = h
    return rows


def rewrite(src, dst, expected):
    n_changed = n_same = 0
    with zipfile.ZipFile(src) as zi, zipfile.ZipFile(dst, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zo:
        for info in zi.infolist():
            data = zi.read(info)
            if info.filename.endswith('/provenance_v10.json'):
                data = anonymize_bytes(data)
                h = hashlib.sha256(data).hexdigest()
                if expected.get(info.filename) != h:
                    raise SystemExit(f'{info.filename}: SHA-256 {h[:16]} differs from data/DATA_SHA256SUMS.txt')
                n_changed += 1
            else:
                n_same += 1
            new = zipfile.ZipInfo(info.filename, date_time=info.date_time)
            new.compress_type = zipfile.ZIP_DEFLATED; new.external_attr = info.external_attr
            zo.writestr(new, data)
    with zipfile.ZipFile(src) as a, zipfile.ZipFile(dst) as b:          # every unchanged member: same CRC-32 and size
        ia = {i.filename: i for i in a.infolist()}; ib = {i.filename: i for i in b.infolist()}
        assert set(ia) == set(ib), 'member lists differ'
        for name, i in ia.items():
            if not name.endswith('/provenance_v10.json'):
                assert (i.CRC, i.file_size) == (ib[name].CRC, ib[name].file_size), name
        assert b.testzip() is None, 'corrupt member in ' + dst
    return n_changed, n_same


def main(src_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True); expected = expected_sums(); lines = []
    for name in ASSETS:
        src = os.path.join(src_dir, name + '.zip'); dst = os.path.join(out_dir, name + '.zip')
        if not os.path.exists(src):
            sys.exit(f'{src} not found')
        with zipfile.ZipFile(src) as z:
            has = any(n.endswith('/provenance_v10.json') for n in z.namelist())
        if has:
            c, s = rewrite(src, dst, expected); print(f'{name}.zip: {c} provenance files rewritten, {s} members unchanged')
        else:
            shutil.copyfile(src, dst); print(f'{name}.zip: copied (no provenance file)')
        lines.append(f'{sha256(dst)}  {name}.zip')
    open(os.path.join(out_dir, SUMS), 'w', newline='\n').write('\n'.join(lines) + '\n')
    print('->', os.path.join(out_dir, SUMS)); print('\n'.join(lines))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
