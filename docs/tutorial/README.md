# MEGAN-SC tutorial: from the MEGAN equations to every number of the article

This tutorial is written for a reader who wants to check, rerun or extend the biogenic-emissions inventory of Santa
Catarina. It assumes some Python and some atmospheric chemistry, and explains the rest. Each chapter can be read on
its own; the reading paths below pick the chapters for a given purpose.

## Illustrated tutorial (PDF)

`TUTORIAL_MEGAN_SC_pt.pdf` (Portuguese) and `TUTORIAL_MEGAN_SC_en.pdf` (English), at the repository root, are the
illustrated tutorial for scientists and lay readers: what the inventory is, MEGANv3.21 in plain words, the inputs,
the repository, reproduction, notebooks 10 to 23, Google Colab with Earth Engine, reading the results, your own
region, limitations, troubleshooting, glossary, references and a printable checklist. Their sources are in this
folder: `tutorial.html`, `tutorial_en.html`, `style.css` and `img/` (copies of figures of `results/v10/`).
`python docs/tutorial/render.py` renders both PDFs with Chromium through the Python package `playwright` (not part of
the requirements files), and `python docs/tutorial/check.py "$PWD/docs/tutorial/tutorial_en.html"` checks that every
piece of text lies inside its box.

## Chapters

| # | Chapter | What you will learn |
|---|---|---|
| 1 | [The science: MEGAN in one chapter](01_science.md) | what the model computes, every activity factor with its equation and constants, the canopy model, soil NO, speciation and units, the settings of the Santa Catarina run, and the eight input and run procedures (FIX-1 to FIX-8) |
| 2 | [The data](02_data.md) | the inputs (ERA5-Land, MODIS, MapBiomas, MEGAN tables, IBGE boundaries), the shipped model output, file formats, and how to open each file in Python or QGIS |
| 3 | [Setting up](03_setup.md) | Python environment, data download and checksum, unit tests, Docker, Google Colab and Earth Engine |
| 4 | [The pipeline](04_pipeline.md) | how the notebooks and scripts connect, from Earth Engine to the last figure, and what each one reads and writes |
| 5 | [Level A: reproduce every result](05_level_A_reproduce.md) | `reproduce.py` step by step, the expected output of each step, and how to read the comparison with the published files |
| 6 | [Level B: run the model again](06_level_B_rerun_the_model.md) | notebook 10 in Colab, the Fortran benchmark, the TROPOMI notebook, the sensitivity runs, and how to compare a rerun with the published run |
| 7 | [The statistics](07_statistics.md) | every estimator used in the documents (Theil-Sen, Mann-Kendall, Hamed-Rao, Benjamini-Hochberg, detrended and partial correlations, effective sample size, Fisher interval, standardised regression, temperature sensitivity) with the code that computes it |
| 8 | [Number map](08_number_map.md) | each number of the article and of the report, the file and key that hold it, and the step that writes it |
| 9 | [Figures and tables](09_figures_and_tables.md) | each figure and table of both documents, its file and its script |
| 10 | [The TROPOMI evaluation](10_tropomi.md) | the physics of the formaldehyde test, the screening, the four pre-specified tests, the diagnostics and the magnitude calculation worked out by hand |
| 11 | [Your own region](11_your_region.md) | the region wizard and how to adapt notebook 10 to another boundary |
| 12 | [Troubleshooting](12_troubleshooting.md) | messages you may see and what to do |
| 13 | [Provenance](13_provenance.md) | how outputs, code, data and results are identified: the `fixes.py` hash stamped in every output and checked by C0, the data checksums, the release assets, the comparison of `reproduce.py compare`, and the citation |
| | [References](references.md) | the literature cited in this tutorial |

## Reading paths

* **"I want to check a number in the article."** Chapter 8, then the file it names in `results/v10/`; chapter 5 if you
  want to recompute it.
* **"I want to reproduce everything on my laptop."** Chapters 3 and 5 (about an hour of reading, four minutes of
  computing).
* **"I want to understand what the model does."** Chapter 1, then `docs/ALGORITHM_SPEC.md` for the full transcription of
  the Fortran.
* **"I want to rerun the model from the raw data."** Chapters 3, 4 and 6.
* **"I want an inventory for my region."** Chapters 1 (section 1.8), 11 and 12.

## Conventions

* Paths are relative to the repository root. Commands are run from the root unless stated otherwise.
* "The run" is the production run of notebook 10 shipped in `data/run_v10/`.
* Numbers quoted in the tutorial come from the files of `results/v10/`; the chapter that quotes a number names its file.
  When a number is rounded here, the file holds the full value.
* BVOC means the 17 MEGAN classes other than soil NO and CO. "SOA precursors" means monoterpenes plus sesquiterpenes,
  as in both documents.
* Months are calendar months in UTC; "summer" means December to February.

## Glossary

| Term | Meaning |
|---|---|
| ER | emission activity: the non-dimensional factor computed by MEGVEA for each class, cell and time step |
| EF | emission factor (nmol m⁻² s⁻¹) under standard conditions, from the ecotype and growth-form tables |
| γ (gamma) | an activity factor: γTP (temperature and light, canopy-integrated), γage (leaf age), γSM (soil moisture), γCO2, γBD (bidirectional exchange of ethanol and acetaldehyde) |
| LAI | leaf area index (m² of leaf per m² of ground), from MODIS |
| LDF | light-dependent fraction of a class |
| PAR, PPFD | photosynthetically active radiation (W m⁻²) and photosynthetic photon flux density (µmol m⁻² s⁻¹) |
| growth form | tree, shrub, herb or crop; MapBiomas classes are mapped to these four |
| ecotype | MEGAN vegetation type that sets the species composition and therefore the emission factors |
| FIX-1 ... FIX-8 | the input and run procedures of `src/megan_sc/fixes.py` (chapter 1, section 1.9) |
| C0 ... C5 | the six automatic checks of notebook 11 |
| S1 ... S3 | the sensitivity tests of notebook 15: S1 time step, S2 LAI aggregation (planned, not run), S3 soil-moisture factor |
| ρ | canopy loss and production factor, 0.96 for isoprene (chapter 1, section 1.1) |
| ΔHCHO | TROPOMI formaldehyde column minus the Atlantic background |
| EECU | Earth Engine compute unit, the measure of Earth Engine batch compute |
