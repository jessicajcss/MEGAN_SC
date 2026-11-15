# Changelog

## 1.0.0 (first public release)

* Engine: Python implementation of the MEGANv3.21 emission code (canopy, activity factors, soil NO, speciation) with
  the input and run procedures FIX-1 to FIX-8 of `src/megan_sc/fixes.py` (SHA-256 prefix a933b29602929a82, stamped in
  every output file), compared stage by stage with the official Fortran executables (`VALIDATION.md`, notebook 12).
* Production run of Santa Catarina, 2003-2025, 0.025°, every 3 h (notebook 10), with 2024 and 2025 computed with the
  VIIRS-corrected LAI (notebooks 16, 17 and 19); the run and the MODIS-LAI years are distributed as release archives
  (`DATA_AVAILABILITY.md`).
* Results of the article, its Online Resource 1 and the state inventory report in `results/v10`, and the sensitivity
  case with the MODIS LAI of 2024-2025 in `results/v10_modis_lai_2024_2025`; checksums in `results/MANIFEST.sha256`.
* Checks and evaluations: notebook 11 (checks C0-C5), TROPOMI formaldehyde (notebook 14), time step and soil moisture
  (notebook 15), 2019 sensitivity runs (notebook 20), ERA5-Land against INMET (notebook 21), VIIRS LAI of three
  satellites (notebooks 16 and 22), comparison with BRAIN 2019 (notebooks 18 and 23), and the time conventions of the
  inputs and of the model (`docs/TIME_CONVENTIONS.md`).
* `reproduce.py`: every number, table and figure recomputed from the release archives in 16 steps and compared with the
  published copies (`docs/tutorial/reproduction_log.txt`).
* Region wizard (`run_your_region/`) for another state, province or country, with the period given as first and last
  month; tutorials in Portuguese and English (`TUTORIAL_MEGAN_SC_pt.pdf`, `TUTORIAL_MEGAN_SC_en.pdf`) and 13 technical
  chapters (`docs/tutorial/`); 127 unit tests.
* Double-blind review: author names, e-mail, ORCID and affiliations are withheld in `CITATION.cff` and in the documents,
  and the field `ee_project` of the provenance files of the release archives reads `"withheld"`
  (`tools/anonymize_release.py`).
