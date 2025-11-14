# Effective Connectivity Reveals Hippocampal Disruptions in Adolescent Anxiety

## Abstract
Anxiety disorders peak during adolescence alongside rapid neural development, yet the neural mechanisms underlying adolescent threat learning remain unclear. Using fMRI and dynamic causal modeling, we mapped effective connectivity during threat conditioning and extinction in 87 adolescents (ages 11-17), focusing on the amygdala, hippocampus, and posterior ventromedial prefrontal cortex (p-vmPFC) circuit. The adolescent brain showed distinct circuit strategies for threat versus safety learning: threat learning (CS+) initially engaged only amygdala and p-vmPFC excitatory coupling without hippocampal involvement, while early safety learning (CS-) recruited hippocampal inhibition of threat-processing regions. By late conditioning, bidirectional inhibitory connectivity emerged between amygdala and p-vmPFC for both stimulus types, suggesting homeostatic regulation once associations consolidated. During extinction, the network transitioned from complex mutual amygdala-hippocampus inhibition during early CS+ trials to streamlined hippocampal-prefrontal pathways by late extinction, reflecting safety consolidation. Adolescents with higher trait anxiety demonstrated stronger hippocampal to p-vmPFC connectivity during early conditioning and altered hippocampal modulation during extinction, revealing temporally dynamic disruptions in directional signaling. This study provides the first effective connectivity map of adolescent threat learning and identifies neural mechanisms that may confer risk for anxiety disorders during this vulnerable developmental period.

## Project Overview
This repository packages the MATLAB + SPM25 code used to recreate the dynamic causal modeling (DCM) and Parametric Empirical Bayes (PEB) analyses for the paper above. The pipeline estimates subject-level GLMs for an fMRIPrep processed threat-conditioning task, extracts volumes of interest (VOIs), builds bilinear DCMs, aggregates them with PEB, and provides reproducible Python plotting scripts for the manuscript figures.

## Repository Layout
- `dcm/` – single-subject GLM specification, VOI extraction, DCM specification, and loading/fitting utilities written in MATLAB/SPM.
- `peb/` – MATLAB scripts for group-level inference (PEB/BMC/BMR), cross-validation, and trait-anxiety analyses, plus the accompanying design CSV.
- `plotting/` – lightweight Python utilities that transform CSV outputs into publication-ready network, model, and violin plots.

## Dependencies
- MATLAB R2022b+ with SPM25 (or SPM12 r7771+) on the MATLAB path.
- fMRIPrep outputs in MNI space for the fear task (`sub-EA###_task-fear#_space-MNI152...bold.nii`) plus 24-parameter motion-confound text files.
- Onset text files (`<condition>.txt`) that encode onset/duration pairs per task condition.
- Binary ROI masks (`hipp`, `postvmpfc`, `amygdala`, whole-brain mask) in MNI space.
- Python 3.10+ with the plotting stack listed in `requirements.txt` (install via `pip install -r requirements.txt`).

Set the `...` placeholders in each script to match your local directory layout before running anything.

## Analysis Workflow
1. **First-level GLMs:** Run `dcm/first_level_dcm.m` for every subject to build SPM designs, estimate GLMs, and create contrasts.
2. **VOI extraction:** Run `dcm/extract_vois.m` for each subject to derive VOIs from the contrasts + anatomical masks.
3. **DCM specification:** Use `dcm/dcm_spec.m` to define the bilinear three-node model per subject (one DCM per specification string).
4. **GCM assembly & fitting:** Execute `dcm/dcm_load_fit.m` once per analysis version/specification to load all subject DCMs into a cell array (`GCM_*.mat`) and fit them with subject-level PEB.
5. **Group PEB/BMC:** Switch to the `peb/` scripts:
   - `PEB.m` summarises intrinsic/modulatory connections and exports diagnostics.
   - `driving_inputs.m` compares alternative driving-input models (C-matrix) using BMC.
   - `STAI.m` fits group covariates (e.g., trait anxiety latent classes) and exports class-specific B-matrix effects.
   - `LOO.m` performs targeted leave-one-out model evidence tests on selected B parameters.
6. **Plotting:** After the MATLAB steps generate CSV/Mat files, use the Python scripts in `plotting/` (`connections_plot.py`, `models_plot.py`, `loo_violin_plot.py`) to recreate the manuscript figures. Each script expects the `peb_results/<analysis>/<spec>/` directory structure created by the MATLAB code.

## Data & Results
- Group-level results are written inside `.../peb_results/<analysis>/<spec>/`. You are expected to version-control derived CSVs/PNGs separately if needed.
- The `peb/stai_design.csv` file provides the STAT-T latent-class regressor used in the paper (EA### IDs, ±1 coding).

Questions or pull requests are welcome if you adopt or extend these scripts for related threat-learning datasets.
