# DCM Pipeline (Subject Level)

This folder contains all MATLAB/SPM25 scripts that prepare single-subject inputs for the effective connectivity analyses. Update the `...` path placeholders in each script before running.

## Dependencies
- MATLAB R2022b+ with SPM25 (or SPM12 r7771+) on the path.
- fMRIPrep outputs in MNI space (`sub-EA###_task-fear#_*_bold.nii` plus `motion24` files).
- Onset text files in `ONSETS_ROOT/fear#/` named exactly like the condition labels in `first_level_dcm.m`.
- Binary ROI masks (`hipp`, `postvmpfc`, `amygdala`, whole-brain mask) in MNI space.

## Script Inventory & Execution Order
1. **`first_level_dcm.m`** – runs a single-subject GLM for the fear task. It creates the subject output folder, specifies/estimates the model, and writes F- and T-contrasts needed later.
2. **`extract_vois.m`** – takes the GLM output for one subject, thresholds the task contrast, intersects it with anatomical masks, and saves VOI time series (`VOI_<roi>_<ID>.mat`).
3. **`dcm_spec.m`** – uses the subject’s VOIs and SPM.mat to create a three-node DCM (hippocampus, posterior vmPFC, amygdala). One `.mat` is produced per specification tag.
4. **`dcm_load_fit.m`** – loops over the subject list, loads all DCM files, fits them with empirical priors (`spm_dcm_peb_fit`), and saves a group cell array `GCM_<version>_<task>_<spec>.mat` for downstream PEB scripts.

## Typical Usage (Batch Script)
These functions are normally orchestrated from via batch bash scripting to run for multiple participants concurrently.
