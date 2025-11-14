function dcm_load_fit(version, spec)
% DCM_LOAD_FIT
%
% Loads and fits subject-level DCMs for the fear task (task 1) using
% Parametric Empirical Bayes (PEB) in SPM25.
%
% For each subject in the hard-coded subject list below, this function:
%   1) Locates the subject-specific DCM file (DCM_*.mat)
%   2) Collects all DCMs into a GCM cell array
%   3) Runs spm_dcm_peb_fit to fit each DCM
%   4) Saves the fitted GCM to disk
%
% INPUTS
%   version : string specifying analysis version (matches first-level / DCM spec)
%   spec    : string tag identifying the DCM specification (used in filenames)
%
% OUTPUT
%   Saves a GCM_*.mat file in the group output directory.
%
% Dependencies:
%   - SPM25
%   - Subject-level DCM files already specified and saved (see dcm_spec_*)

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS
% -------------------------------------------------------------------------
% Root directory for first-level and DCM output (must match other scripts)
OUTPUT_ROOT = '.../data/output';

%% ------------------------------------------------------------------------
%  GENERAL SETUP
% -------------------------------------------------------------------------
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

% Task information
task_num     = 1;
task_num_str = num2str(task_num);
fear_task    = ['fear_', task_num_str];
folder       = fear_task;

% Group-level output directory for this analysis version
output_dir = fullfile(OUTPUT_ROOT, version, folder);

%% ------------------------------------------------------------------------
%  SUBJECT LIST
% -------------------------------------------------------------------------
% List of subject IDs included in the group DCM analysis.
% This is dataset-specific and should be adapted for other cohorts.
subject_ids = [ ...
     8,  9, 10, 16, 21, 26, 36, 37, 39, 44, 48, 49, 58, 92, 94, ...
   104,115,120,121,136,143,146,156,161,169,190,193,198,220,221, ...
   222,254,282,303,323,341,350,363,375,393,398,403,407,421,433, ...
   448,449,485,528,544,595,596,641,643,660,673,682,688,693,714, ...
   730,746,747,760,777,779,806,840,847,858,861,864,866,870,875, ...
   885,889,895,898,910,913,914,918,920,939,946,1000 ...
   ];
num_subjects = numel(subject_ids);

% Initialise container for file information returned by dir()
filelist = '';

%% ------------------------------------------------------------------------
%  LOCATE SUBJECT-LEVEL DCM FILES
% -------------------------------------------------------------------------
for i = 1:num_subjects
    subject_id = subject_ids(i);

    % Subject ID string (EA1000 is a special case in this dataset)
    if subject_id == 1000
        subject_id_str = '1000';
    else
        subject_id_str = sprintf('%03d', subject_id);
    end

    % Subject-specific directory containing DCM_*.mat
    sub_dir = fullfile(output_dir, subject_id_str);

    % Pattern for subject-specific DCM to load
    % (spm_dcm_load expects a list of full filenames)
    DCM_full    = fullfile(sub_dir, ['DCM_', spec, '_', subject_id_str, '.mat']);
    cur_filelist = dir(DCM_full);

    % Accumulate file information across subjects
    filelist = [filelist, cur_filelist]; %#ok<AGROW>
end

%% ------------------------------------------------------------------------
%  FIT DCMs USING PEB
% -------------------------------------------------------------------------
% Convert struct array from dir() into a cell array of full file paths
GCM = fullfile({filelist.folder}, {filelist.name})';

% Load all DCMs into a GCM cell array
GCM = spm_dcm_load(GCM);

% Fit each DCM with empirical priors (PEB)
GCM = spm_dcm_peb_fit(GCM);

% Save group GCM to disk (v7.3 for potentially large data)
save(fullfile(output_dir, ['GCM_', version, '_', folder, '_', spec, '.mat']), ...
     'GCM', '-v7.3');

end