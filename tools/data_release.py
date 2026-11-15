#!/usr/bin/env python3
"""data_release.py - checksums, release archives and download of the data folders that are kept out of git.

The large data folders are distributed as release assets (GitHub Releases or Zenodo) and unpacked into data/: the
production run (2024 and 2025 with the VIIRS-corrected LAI once notebook 19 has installed them), the TROPOMI fields, the
two runs of test S1, the LAI cache of notebook 16 and the archive of 2024 and 2025 computed with the MODIS LAI. The small inputs (data/boundaries, data/megan_inputs, data/sc_ecotype.nc) are in git.

    python tools/data_release.py sums               # write data/DATA_SHA256SUMS.txt (every file under data/)
    python tools/data_release.py pack <out_dir>     # one zip per large folder + <out_dir>/SHA256SUMS_assets.txt
    python tools/data_release.py fetch <base_url>   # download the zips from <base_url>/<name>.zip, check, unpack into data/
    python tools/data_release.py fetch <zip_dir>    # the same from a local folder holding the zips

<base_url> is, for a GitHub release, https://github.com/<owner>/<repo>/releases/download/<tag> ; for a Zenodo record,
https://zenodo.org/records/<record>/files (the script adds ?download=1). Every zip is checked before it is unpacked,
against SHA256SUMS_assets.txt next to the zips or, failing that, data/RELEASE_ASSETS_SHA256SUMS.txt (kept in git).
After `fetch`, run `python reproduce.py check-data`,
which checks every unpacked file against data/DATA_SHA256SUMS.txt (kept in git).
"""
import hashlib, os, shutil, sys, urllib.request, zipfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, 'data')
ASSETS = {  # zip name -> folder under data/ (what each one holds is described in data/README.md)
    'megan_sc_run_v10_2003-2025': 'run_v10',
    'tropomi_hcho_v10_2019-2025': 'tropomi_hcho_v10',
    'megan_sc_s1_test_2019': 's1_test',
    'lai_viirs_v10_2013-2025': 'lai_viirs_v10',
    'megan_sc_run_v10_modis_lai_2024_2025': 'run_v10_modis_lai_2024_2025',   # MODIS-LAI 2024-2025 (notebook 19 archive)
}
VIEWS = ('run_v10_modis_view',)          # folders of links built by tools/lai_variant_view.py: never summed or packed
SUMS = 'SHA256SUMS_assets.txt'


def sha256(path):
    d = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            d.update(b)
    return d.hexdigest()


def sums():
    rows = []
    for root, dirs, files in os.walk(DATA):
        dirs.sort()
        dirs[:] = [x for x in dirs if not x.startswith(('.', '_download')) and not (root == DATA and x in VIEWS)]
        for f in sorted(files):
            p = os.path.join(root, f); rel = os.path.relpath(p, DATA).replace(os.sep, '/')
            if rel in ('DATA_SHA256SUMS.txt', 'README.md', 'RELEASE_ASSETS_SHA256SUMS.txt') or f.startswith('.'):
                continue
            rows.append(f'{sha256(p)}  {rel}')
    open(os.path.join(DATA, 'DATA_SHA256SUMS.txt'), 'w', encoding='utf-8', newline='\n').write('\n'.join(rows) + '\n')
    print(f'{len(rows)} files -> data/DATA_SHA256SUMS.txt')


def pack(out):
    os.makedirs(out, exist_ok=True); lines = []
    for name, folder in ASSETS.items():
        src = os.path.join(DATA, folder); z = os.path.join(out, name + '.zip')
        if not os.path.isdir(src):
            sys.exit(f'{src} not found')
        with zipfile.ZipFile(z, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            for root, dirs, files in os.walk(src):
                dirs.sort()
                for f in sorted(files):
                    p = os.path.join(root, f); zf.write(p, os.path.relpath(p, DATA).replace(os.sep, '/'))
        h = sha256(z); lines.append(f'{h}  {name}.zip'); print(f'{name}.zip  {os.path.getsize(z) / 1e6:,.0f} MB  {h[:16]}')
    open(os.path.join(out, SUMS), 'w', newline='\n').write('\n'.join(lines) + '\n')
    print('->', os.path.join(out, SUMS))


def _get(src, name, dst):
    """Copy (local folder) or download (URL) one release file; return False if it is not there."""
    if os.path.isdir(src):
        p = os.path.join(src, name)
        if not os.path.exists(p):
            return False
        shutil.copyfile(p, dst); return True
    url = src.rstrip('/') + '/' + name + ('?download=1' if 'zenodo.org' in src else '')
    try:
        print('downloading', url, flush=True)
        urllib.request.urlretrieve(url, dst); return True
    except Exception as e:                       # noqa: BLE001  (report and let the caller decide)
        print('   could not download:', e)
        return False


def fetch(src):
    tmp = os.path.join(DATA, '_download'); os.makedirs(tmp, exist_ok=True)
    expected = {}
    local = os.path.join(DATA, 'RELEASE_ASSETS_SHA256SUMS.txt')        # copy kept in git
    sums_file = os.path.join(tmp, SUMS) if _get(src, SUMS, os.path.join(tmp, SUMS)) else (local if os.path.exists(local) else None)
    if sums_file:
        for line in open(sums_file, encoding='utf-8'):
            if line.strip():
                h, n = line.split(None, 1); expected[n.strip()] = h
    else:
        print(f'{SUMS} not found: the zips are unpacked without checking their checksums')
    for name in ASSETS:
        z = os.path.join(tmp, name + '.zip')
        if not _get(src, name + '.zip', z):
            sys.exit(f'{name}.zip not found at {src}')
        if expected:
            h = sha256(z)
            if expected.get(name + '.zip') != h:
                sys.exit(f'{name}.zip: checksum {h[:16]} differs from {SUMS}; download again')
            print(f'{name}.zip checksum ok')
        with zipfile.ZipFile(z) as zf:
            zf.extractall(DATA)
        print('unpacked', name, flush=True); os.remove(z)
    shutil.rmtree(tmp, ignore_errors=True)
    print('done: now run  python reproduce.py check-data')


if __name__ == '__main__':
    a = sys.argv[1:]
    if a[:1] == ['sums']:
        sums()
    elif a[:1] == ['pack'] and len(a) == 2:
        pack(a[1])
    elif a[:1] == ['fetch'] and len(a) == 2:
        fetch(a[1])
    else:
        sys.exit(__doc__)
