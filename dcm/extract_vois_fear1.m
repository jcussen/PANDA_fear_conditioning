function extract_vois_fear1(subject_id, version)
% EXTRACT_VOIS_FEAR1
%
% Extracts subject-specific VOIs for a single contrast (e.g. task > baseline)
% for the fear task (task 1), using SPM25. VOIs are defined as the
% conjunction of:
%   (1) A thresholded SPM contrast map, and
%   (2) An anatomical ROI mask (e.g. hippocampus, vmPFC, amygdala).
%
% INPUTS
%   subject_id : numeric ID (e.g. 1, 2, ..., 999, 1000)
%   version    : string specifying analysis version (matches first-level code)
%
% This function operates on a single subject. In typical usage, it is called
% from a wrapper script that loops over all subjects.
%
% Dependencies:
%   - SPM25
%   - First-level SPM.mat and contrasts already estimated (see GLM script)
%   - Binary ROI masks in MNI space (NIfTI format)

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS
% -------------------------------------------------------------------------
% Root directory for first-level SPM output (must match GLM script)
OUTPUT_ROOT = '.../data/output';

% Directory containing anatomical ROI masks
MASKS_ROOT = '.../data/masks';

%% ------------------------------------------------------------------------
%  GENERAL SETUP
% -------------------------------------------------------------------------
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

% Contrast index in SPM.xCon to use for VOI definition
% (e.g. 2 = task-related contrast; must correspond to your first-level setup)
contrast_num = 2;

% Cluster-defining p-value threshold for the contrast image
p_threshold = 0.05;

% Task parameters (fear task 1)
task_num     = 1;
task_num_str = num2str(task_num);
fear_task    = ['fear_', task_num_str];
folder       = fear_task;

% ROI masks (binary images in MNI space)
masks_dir = MASKS_ROOT;
masks = struct( ...
    'hipp',       fullfile(masks_dir, 'hipp_bin_clean.nii'), ...
    'postvmpfc',  fullfile(masks_dir, 'postvmpfc_bin_clean.nii'), ...
    'amygdala',   fullfile(masks_dir, 'amygdala_bin_clean.nii'), ...
    );
roi_names = fieldnames(masks);

% Output directory for this analysis version
output_dir = fullfile(OUTPUT_ROOT, version, folder);

%% ------------------------------------------------------------------------
%  SUBJECT-SPECIFIC SETUP
% -------------------------------------------------------------------------
disp(['Processing subject ', num2str(subject_id), '...']);

% Subject ID string (EA1000 is a special case in this dataset)
if subject_id == 1000
    subject_id_str = '1000';
else
    subject_id_str = sprintf('%03d', subject_id);
end

% Directory containing SPM.mat and contrast images for this subject
spm_dir = fullfile(output_dir, subject_id_str);

clear matlabbatch

%% ------------------------------------------------------------------------
%  VOI DEFINITION PER ROI
% -------------------------------------------------------------------------
for r = 1:length(roi_names)  % Loop over anatomical ROIs

    % Link to subjects first-level model
    matlabbatch{r}.spm.util.voi.spmmat(1) = {fullfile(spm_dir, 'SPM.mat')};

    % Adjust data for F-contrast #1 (Effects of interest)
    % (0 = no adjustment, NaN = adjust for all effects)
    matlabbatch{r}.spm.util.voi.adjust    = 1;
    matlabbatch{r}.spm.util.voi.session   = 1;
    matlabbatch{r}.spm.util.voi.name      = [roi_names{r}, '_', subject_id_str];

    % ---------------------------------------------------------------------
    % First ROI component: thresholded SPM contrast map
    % ---------------------------------------------------------------------
    matlabbatch{r}.spm.util.voi.roi{1}.spm.spmmat    = {''};
    matlabbatch{r}.spm.util.voi.roi{1}.spm.contrast  = contrast_num;   % e.g. task contrast
    matlabbatch{r}.spm.util.voi.roi{1}.spm.conjunction = 1;
    matlabbatch{r}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
    matlabbatch{r}.spm.util.voi.roi{1}.spm.thresh    = p_threshold;    % e.g. p < .05
    matlabbatch{r}.spm.util.voi.roi{1}.spm.extent    = 0;              % no cluster extent threshold
    matlabbatch{r}.spm.util.voi.roi{1}.spm.mask      = struct('contrast', {}, 'thresh', {}, 'mtype', {});

    % ---------------------------------------------------------------------
    % Second ROI component: anatomical mask for the current ROI
    % ---------------------------------------------------------------------
    matlabbatch{r}.spm.util.voi.roi{2}.mask.image     = {fullfile(masks_dir, [roi_names{r}, '_bin_clean.nii'])};
    matlabbatch{r}.spm.util.voi.roi{2}.mask.threshold = 0.5;

    % Final VOI is the conjunction of (1) functional activation & (2) anatomical mask
    matlabbatch{r}.spm.util.voi.expression = 'i1 & i2';
end

% Run VOI extraction for all ROIs for this subject
spm_jobman('run', matlabbatch);

%% ------------------------------------------------------------------------
%  VOXEL COUNT SUMMARY (WRITTEN TO STDOUT)
% -------------------------------------------------------------------------
for r = 1:numel(roi_names)
    voi_nii = fullfile(spm_dir, ['VOI_', roi_names{r}, '_', subject_id_str, '_mask.nii']);
    if exist(voi_nii, 'file')
        nvox = nnz(spm_read_vols(spm_vol(voi_nii)));
        fprintf('Subject %s  ROI %-10s : %d voxels in final mask\n', ...
                subject_id_str, roi_names{r}, nvox);
    else
        fprintf(2, 'WARNING: VOI file %s not found – skipping voxel count\n', voi_nii);
    end
end

end