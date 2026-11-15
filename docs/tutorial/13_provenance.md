# 13 Provenance: identifying the code, the data and the results

Every number of the article and of the state inventory report can be traced to the code that produced the model
output, to the data files that hold it, and to the script that computed the number. This chapter shows how.

## 13.1 The code

Notebook 10 writes `src/megan_sc/fixes.py` (the input and run procedures FIX-1 to FIX-8, chapter 1, section 1.9) into
the engine folder before importing it, and records the first 16 characters of its SHA-256, a933b29602929a82:

* in the global attribute `megan_sc_fixes_sha256` of every monthly and annual NetCDF file of `data/run_v10/`;
* in the key `fixes_sha256` of `data/run_v10/<year>/provenance_v10.json`, together with the sampling interval, the LAI
  aggregation, the wilting point, the engine patches, the ecotype file and the Earth Engine project of that year.

Check C0 of notebook 11 verifies that every year carries the same hash and the three engine patches. You can check it
yourself from the repository root:

```python
import hashlib, json
from netCDF4 import Dataset
print(hashlib.sha256(open('src/megan_sc/fixes.py', encoding='utf-8').read().encode('utf-8')).hexdigest()[:16])
with Dataset('data/run_v10/2025/emissions/annual_emis_2025.nc') as ds:
    print(ds.megan_sc_fixes_sha256)
print(json.load(open('data/run_v10/2025/provenance_v10.json'))['fixes_sha256'])
```

The three lines print the same value. A run made with another version of `fixes.py` carries another hash, so its
output can never be mistaken for the one shipped here.

## 13.2 The data

`data/DATA_SHA256SUMS.txt` holds the SHA-256 of every file under `data/` (2,435 files, the README and the two checksum lists excepted), and `python reproduce.py check-data` compares the files on your disk with it. The release
archives are checked before they are unpacked, against `SHA256SUMS_assets.txt` next to the zip files or, failing that,
`data/RELEASE_ASSETS_SHA256SUMS.txt` (chapter 3). The monthly inputs as sampled by notebook 10 (meteorology, LAI,
growth forms, emission factors and the fix logs, `data/run_v10/<year>/inputs/`) will be part of the release, so the run can
be inspected without Earth Engine (chapter 2).

## 13.3 The results

`results/v10/` holds the published copy of every table, number and figure, and chapter 8 maps each number of both
documents to its file and key. `python reproduce.py all` recomputes them from `data/` into `reproduced/` and compares
each file with its published copy (chapter 5); `docs/tutorial/reproduction_log.txt` keeps the record of the run of
15 November 2026.

## 13.4 Citing

`CITATION.cff` holds the metadata of the repository and of the article (GitHub shows it under "Cite this
repository"); the DOI of the archived version will be added there after deposit. Cite also the MEGANv3.21 code (Guenther, 2024b,
Zenodo 10.5281/zenodo.10526206), the MEGAN 3.2 input files (Guenther, 2024a, Zenodo 10.5281/zenodo.10939297) and the
data products listed in `DATA_AVAILABILITY.md`. `CHANGELOG.md` lists what each release contains.
