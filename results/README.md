# Results

* `v10/`: results of the Santa Catarina example with the VIIRS-corrected LAI of 2024 and 2025 (notebooks 11, 13, 14,
  15 and 18-23 and the scripts of `../analysis/v10/`), and the reference of `python reproduce.py compare`; see
  `v10/README.md`.
* `v10_modis_lai_2024_2025/`: the same analyses with the uncorrected MODIS LAI of 2024 and 2025, the sensitivity of the
  results to the LAI of those two years.
* `MANIFEST.sha256`: SHA-256 of every file of both folders (`python tools/verify_manifest.py`).

A run of another region writes its results next to its own output (`<run>/results/`, notebook 25), not here.
