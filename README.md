# HKGI Flagship Paper 2025

## Analysis

## Figures

### Steps to reproduce figures

1. Build an apptainer image containing Jupyter notebook with all R and python dependencies installed: `./figures/build/build.sh`
2. Run the created apptainer image. This will start a Jupyter notebook: `apptainer run r_hkgi_flagship_figures-v0.0.1.sif`
3. Open the http://localhost link that is printed to the terminal to connect to the Jupyter notebook.
4. Within Jupyter notebook, run the notebooks under `figures/src`
