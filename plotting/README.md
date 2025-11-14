# Plotting Utilities

Python scripts that convert the CSV/MAT outputs from the MATLAB pipeline into camera-ready figures. Each script assumes that `.../peb_results/<analysis>/<spec>/` contains the files produced by `peb/` (e.g., `connections.csv`, `class_connections.csv`, `driving_inputs_models.csv`, `M.mat`, and `LOO_*.mat`).

## Dependencies
Install via the top-level `requirements.txt`:

```
pip install -r requirements.txt
```

## Scripts
1. **`connections_plot.py`** – loads `connections.csv` (and `class_connections.csv` if present), draws directed graphs for intrinsic (A) and modulatory (B) connections above a posterior probability threshold, and writes PNGs to a `TNR/` subfolder.
2. **`models_plot.py`** – reads `driving_inputs_models.csv` (plus `driving_inputs_regions.csv` for completeness) and renders Bayesian model selection bar charts for the candidate driving-input configurations.
3. **`loo_violin_plot.py`** – ingests `LOO_*.mat` files alongside `M.mat`, recreates the STAT-T group splits, computes correlations, and outputs violin/strip plots summarising leave-one-out predictions.

## Typical Usage
Edit the `ROOT_RESULTS`/`BASE_DIR`, `ANALYSIS_SUBDIR`, `SPEC_NAME`, and file lists at the bottom of each script so they point at the folder created by the MATLAB runs, then execute:

```
python plotting/connections_plot.py
python plotting/models_plot.py
python plotting/loo_violin_plot.py
```

All figures are saved as 300 dpi PNGs in the same results directory (or its `TNR` child for network graphs).
