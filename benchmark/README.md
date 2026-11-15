# Benchmark against the official MEGANv3.21 Fortran (`benchmark/`)

Scripts used to compare the Python engine with the executables built from `MEGANv3.21.tar.gz`. Notebook
`../notebooks/12_fortran_benchmark.ipynb` embeds them and runs the whole comparison in Colab; the files here are the
same code for running on a Linux machine.

| File | Purpose |
|---|---|
| `run_fortran_se_us.sh`, `run_fortran_se_us_soil_spec.sh` | build and run the Fortran on the test case shipped with the code (south-east US, 2 May 2018), with GAMSM/GAMBD off and on, and the soil-NO and speciation steps |
| `bench_se_us.py`, `bench_soil.py`, `bench_spec.py` | run the Python engine on the same inputs and compare stage by stage (DAYMET, MEGCAN, MEGVEA, end to end; soil NO; speciation and units) |
| `make_sc_inputs.py`, `run_fortran_sc.sh`, `bench_sc.py` | the same on the Santa Catarina grid for 15 January and 15 July 2019 |
| `ioapi.py` | minimal reader of the IOAPI/netCDF files written by the Fortran |
| `make_table.py` | prints the table of `../VALIDATION.md` from `results_colab/benchmark_summary.csv` |
| `results/` | CSV results of the local run |
| `results_colab/` | CSV results of the Colab run of notebook 12 (identical to `results/`) and `benchmark_summary.csv` |
| `colab_run/12_fortran_benchmark_ran.ipynb` | the executed notebook, with its printed output |

The shell scripts contain the paths of the machine where they ran; edit the variables at the top before reuse.
