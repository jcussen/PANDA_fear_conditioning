# Group PEB & Cross-Validation

This folder hosts MATLAB code for the group-level Parametric Empirical Bayes (PEB) analyses, Bayesian model comparison, covariate modelling, and cross-validation built on top of the subject-level DCMs created in `dcm/`.

## Dependencies
- MATLAB with SPM25 on the path (`addpath(SPM_PATH); spm('Defaults','fMRI')`).
- Fitted group cell array `GCM_<version>_<task>_<spec>.mat` produced by `dcm/dcm_load_fit.m`.
- Writable results root (default `.../peb_results/<analysis>/<spec>/`).
- Optional: `stai_design.csv` (STAT-T latent class coding) for trait-anxiety analyses.

## Files & Recommended Order
1. **`PEB.m`** – loads a `GCM_*.mat`, computes variance explained, runs PEB on the A- and B-matrices, thresholds posterior probabilities, and saves `connections.csv` plus summary PNGs. Run this first after the GCM is ready.
2. **`driving_inputs.m`** – compares alternative C-matrix (driving input) configurations via Bayesian model comparison and outputs `driving_inputs_models.csv` + `driving_inputs_regions.csv`.
3. **`STAI.m`** – incorporates the STAT-T latent-class regressor (or any covariate in `design_csv`) into the PEB design, performs BMR/BMA focused on the B-matrix, and writes `class_connections.csv`. It also saves `PEB.mat`, `BMA_search.mat`, and `M.mat` (the design struct).
4. **`LOO.m`** – runs leave-one-out cross-validation for a single B-parameter specified by the `src`, `tgt`, and `condID` choices at the top of the script, using `M.mat` for the design matrix. Outputs paired `.png` and `.mat` files per connection.
5. **`stai_design.csv`** – subject IDs (EA###) with ±1 coding for the latent class regressor. Update or replace this file when using different covariates.

## Example Workflow
```matlab
spec = '.../peb_results/fear_task/GCM_base/GCM_base.mat';

% 1) Group summary
run('PEB.m');            % after editing the spec/rootDir parameters

% 2) Driving-input comparison
run('driving_inputs.m');

% 3) Trait-anxiety covariate analysis
run('STAI.m');

% 4) Leave-one-out evidence for a chosen B connection
run('LOO.m');
```

All scripts expect the `spec`, `rootDir`, and (where applicable) `design_csv` variables at the top of the file to point to valid locations before you execute them.
