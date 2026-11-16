# Results with the uncorrected MODIS LAI of 2024 and 2025

The same analyses as `../v10/`, from the production run with the MODIS LAI in every year: `data/run_v10` for 2003-2023
and `data/run_v10_modis_lai_2024_2025` for 2024 and 2025, read together through the view of `tools/lai_variant_view.py`.
This folder measures the sensitivity of the results of the example to the LAI of 2024-2025; every other result of the
example comes from `../v10/`. The folders, files and scripts are those of
`../v10/README.md`, without the outputs of the notebooks that do not depend on the LAI of 2024-2025 (BRAIN, INMET,
NOAA-20/21, the 2019 sensitivity runs, notebook 19) and without `posthoc_review/ecotype_tree_ef.csv`.
Isoprene includes the canopy loss and production factor ρ = 0.96 (Guenther et al., 2006), applied when the run is read.

`REPRODUCE.md` (section B) gives the command that recomputes this folder from the shipped data and compares the result
with these files.
