# MEGAN_SC tutorial: from the MEGAN equations to the biogenic inventory of any region

This tutorial is written for a reader who wants to understand, check, rerun or extend a biogenic-emissions inventory
made with this repository: the worked example of Santa Catarina (Brazil, 2003-2025) shipped with it, or a new
inventory for another region and period. It assumes some Python and some atmospheric chemistry, and explains the rest.
Each chapter can be read on its own; the reading paths below pick the chapters for a given purpose.

## Illustrated tutorial (PDF)

`TUTORIAL_MEGAN_SC_pt.pdf` (Portuguese) and `TUTORIAL_MEGAN_SC_en.pdf` (English), at the repository root, are the
illustrated tutorial for scientists and lay readers: MEGANv3.21 in plain words, the inputs, the repository, the Santa
Catarina example, Google Drive, Google Colab and Earth Engine, your own region, reading the results, limitations,
troubleshooting, glossary, references and a printable checklist. Their sources are in this folder: `tutorial.html`,
`tutorial_en.html`, `style.css` and `img/` (copies of figures of `results/v10/`). `python docs/tutorial/render.py`
renders both PDFs with Chromium through the Python package `playwright` (not part of the requirements files), and
`python docs/tutorial/check.py "$PWD/docs/tutorial/tutorial_en.html"` checks that every piece of text lies inside its
box.

## Chapters

| # | Chapter | What you will learn |
|---|---|---|
| 1 | [The science: MEGAN in one chapter](01_science.md) | what the model computes, every activity factor with its equation and constants, the canopy model, soil NO, speciation and units, the eight input and run procedures (FIX-1 to FIX-8) and the settings of the example run |
| 2 | [The data](02_data.md) | the inputs (ERA5-Land, MODIS, MapBiomas, MEGAN tables, IBGE boundaries), the shipped model output of the example, file formats, and how to open each file in Python or QGIS |
| 3 | [Setting up](03_setup.md) | the repository on a computer and in Google Drive, the Python environment, the unit tests, the release archives, and Google Earth Engine in Colab, on any computer and with a service account |
| 4 | [The pipeline](04_pipeline.md) | how the configuration, the notebooks and the scripts connect, from Earth Engine to the last figure, and what each one reads and writes |
| 5 | [Reproduce the Santa Catarina example](05_reproduce_the_example.md) | `reproduce.py` step by step, the expected output of each step, and how to read the comparison with the published files |
| 6 | [Run the model](06_run_the_model.md) | the run configuration, notebooks 10 (model), 11 (checks) and 25 (results) in Colab or Jupyter, and the other notebooks of the example (TROPOMI, sensitivity, LAI check, benchmark) |
| 7 | [The statistics](07_statistics.md) | every estimator used for the results (Theil-Sen, Mann-Kendall, Hamed-Rao, Benjamini-Hochberg, detrended and partial correlations, effective sample size, Fisher interval, standardised regression, temperature sensitivity) with the code that computes it |
| 8 | [Your own region and period](08_your_region.md) | the region wizard, the run of notebooks 10, 11 and 25 with its configuration, and the choices to review for another region |
| 9 | [Reading the results](09_reading_the_results.md) | units and conventions, annual totals, seasonal cycle, maps of intensity, drivers and correlations, trends, and what the results can and cannot support |
| 10 | [The TROPOMI evaluation](10_tropomi.md) | the physics of the formaldehyde test, the screening, the four pre-specified tests, the diagnostics and the magnitude calculation worked out by hand |
| 11 | [Troubleshooting](11_troubleshooting.md) | messages you may see and what to do |
| 12 | [Provenance](12_provenance.md) | how outputs, code, data and results are identified (the `fixes.py` hash, the data checksums, the release assets, the comparison of `reproduce.py compare`), what to publish with a run, and the citation |
| | [References](references.md) | the literature cited in this tutorial and in the notebooks |

## Reading paths

* **"I want an inventory for my region and period."** Chapters 3, 8 and 6, then 9.
* **"I want to understand what the model does."** Chapter 1, then `docs/ALGORITHM_SPEC.md` for the full transcription of
  the Fortran.
* **"I want to reproduce the Santa Catarina example on my computer."** Chapters 3 and 5 (about an hour of reading, a few
  minutes of computing).
* **"I want to rerun the example from the raw data."** Chapters 3, 4 and 6.
* **"I want to interpret or check a result."** Chapters 9 and 7, then the file named in the text.

## Conventions

* Paths are relative to the repository root. Commands are run from the root unless stated otherwise.
* "The example" is the Santa Catarina run of notebook 10 shipped in `data/run_v10/` (configuration
  `config/santa_catarina_2003_2025.yaml`), and `results/v10/` holds everything computed from it.
* Numbers quoted in the tutorial come from the files of `results/v10/`; the chapter that quotes a number names its file.
  When a number is rounded here, the file holds the full value.
* BVOC means the 17 MEGAN classes other than soil NO and CO. "SOA precursors" means monoterpenes plus sesquiterpenes.
* Months are calendar months in UTC; for Santa Catarina, "summer" means December to February.

## Glossary

| Term | Meaning |
|---|---|
| ER | emission activity: the non-dimensional factor computed by MEGVEA for each class, cell and time step |
| EF | emission factor (nmol m⁻² s⁻¹) under standard conditions, from the ecotype and growth-form tables |
| γ (gamma) | an activity factor: γTP (temperature and light, canopy-integrated), γage (leaf age), γSM (soil moisture), γCO2, γBD (bidirectional exchange of ethanol and acetaldehyde) |
| LAI | leaf area index (m² of leaf per m² of ground), from MODIS |
| LDF | light-dependent fraction of a class |
| PAR, PPFD | photosynthetically active radiation (W m⁻²) and photosynthetic photon flux density (µmol m⁻² s⁻¹) |
| growth form | tree, shrub, herb or crop; MapBiomas and MODIS classes are mapped to these four |
| ecotype | MEGAN vegetation type that sets the species composition and therefore the emission factors |
| run configuration | the YAML file of `config/` that defines a run: period, grid, region, inputs, output (chapter 6, section 6.2) |
| FIX-1 ... FIX-8 | the input and run procedures of `src/megan_sc/fixes.py` (chapter 1, section 1.9) |
| C0 ... C5 | the six automatic checks of notebook 11 |
| S1 ... S3 | the sensitivity tests of the example (notebook 15): S1 time step, S2 LAI aggregation (planned, not run), S3 soil-moisture factor |
| ρ | canopy loss and production factor, 0.96 for isoprene (chapter 1, section 1.1) |
| ΔHCHO | TROPOMI formaldehyde column minus the Atlantic background |
| EECU | Earth Engine compute unit, the measure of Earth Engine batch compute |
