function first_level_dcm(subject_id, version)
% FIRST_LEVEL_DCM
%
% Single-subject, single-session first-level GLM for the fear task
% (task 1) preprocessed with fMRIPrep and analysed in SPM25.
%
% This script was written for DCM preparation. It:
%   1) Specifies a first-level model
%   2) Estimates the GLM
%   3) Defines condition-wise and differential contrasts
%
% INPUTS
%   subject_id : numeric ID (e.g. 1, 2, ..., 999, 1000)
%   version    : string specifying analysis version (used in output path)
%
% Dependencies:
%   - SPM25 (https://www.fil.ion.ucl.ac.uk/spm/)
%   - fMRIPrep-preprocessed BOLD in MNI space (single 4D NIfTI per subject)
%   - Onset text files: one per condition with [onset duration] columns
%   - Motion confounds file: "motion24" 24-parameter regressors

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS
%  (Set these for your local environment before running.)
% -------------------------------------------------------------------------
% Root directory containing fMRIPrep outputs (BOLD files)
DATA_ROOT = '.../panda_fmri_prep';

% Root directory where first-level SPM output will be written
OUTPUT_ROOT = '.../data/output';

% Root directory for onset files
ONSETS_ROOT = '.../onsets';

% Whole-brain analysis mask (MNI space, 2 mm)
MASK_FILE = '.../data/masks/MNI152_T1_2mm_brain_mask_new_bin_clean.nii,1';

% Motion confounds file name pattern (24-parameter model)
MOTION_FILE_PATTERN = 'sub-EA%s_task-fear%s_motion24.txt';

%% ------------------------------------------------------------------------
%  GENERAL SETUP
% -------------------------------------------------------------------------
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

% Design parameters
task_num     = 1;                  % fear task number (fixed here to 1)
task_num_str = num2str(task_num);
TR           = 3;                  % Repetition time in seconds
hpf_cutoff   = 128;                % High-pass filter cut-off (seconds)

% Condition names (order defines regressor ordering in the design matrix)
conditions = { ...
    'Task_Question', ...
    'CS+_earlyconditioning', 'CS-_earlyconditioning', ...
    'CS+_earlyextinction',   'CS-_earlyextinction', ...
    'CS+_lateconditioning',  'CS-_lateconditioning', ...
    'CS+_lateextinction',    'CS-_lateextinction', ...
    'Fixation', 'Scream', ...
    };

% Subject ID as string (EA1000 is special case in this dataset)
if subject_id == 1000
    subject_id_str = '1000';
else
    subject_id_str = sprintf('%03d', subject_id);
end

% Input and output directories
sub_dir       = fullfile(DATA_ROOT, ['sub-EA' subject_id_str], 'func');
fear_task     = ['fear_', task_num_str];
folder        = fear_task;
output_dir    = fullfile(OUTPUT_ROOT, version, folder);
sub_output_dir = fullfile(output_dir, subject_id_str);

% Onset directory for the current task
onset_dir = fullfile(ONSETS_ROOT, ['fear' task_num_str '_onsets']);

%% ------------------------------------------------------------------------
% 1) CREATE SUBJECT-SPECIFIC OUTPUT DIRECTORY
% -------------------------------------------------------------------------
clear matlabbatch;
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {output_dir};
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name   = subject_id_str;

%% ------------------------------------------------------------------------
% 2) SPECIFY FIRST-LEVEL GLM
% -------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir              = {sub_output_dir};
matlabbatch{2}.spm.stats.fmri_spec.timing.units     = 'secs';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT        = TR;
matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t    = 16;
matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0   = 8;
matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];   % Canonical HRF only
matlabbatch{2}.spm.stats.fmri_spec.volt             = 1;
matlabbatch{2}.spm.stats.fmri_spec.global           = 'None';
matlabbatch{2}.spm.stats.fmri_spec.mthresh          = 0.01;
matlabbatch{2}.spm.stats.fmri_spec.mask             = {MASK_FILE};
matlabbatch{2}.spm.stats.fmri_spec.cvi              = 'FAST';

% fMRIPrep BOLD file pattern (single 4D NIfTI per subject, per task)
pattern = [ ...
    '^sub-EA' subject_id_str ...
    '_task-fear' task_num_str ...
    '_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold\.nii$'];

% Select 4D BOLD file (all volumes)
files = spm_select('ExtFPList', sub_dir, pattern, Inf);

if isempty(files)
    error('Subject %s – task %s missing, aborting single-session setup', ...
          subject_id_str, task_num_str);
end

% Single-session design (all volumes from the 4D file)
matlabbatch{2}.spm.stats.fmri_spec.sess(1).scans     = cellstr(files);
matlabbatch{2}.spm.stats.fmri_spec.sess(1).multi     = {''};

% -------------------------------------------------------------------------
% Motion confounds
% -------------------------------------------------------------------------
confounds_file = fullfile(sub_dir, ...
    sprintf(MOTION_FILE_PATTERN, subject_id_str, task_num_str));

% Read 24-parameter motion file and forward-fill any missing values
M       = readmatrix(confounds_file, 'Delimiter', '\t');
M_filled = fillmissing(M, 'next');
writematrix(M_filled, confounds_file, 'Delimiter', 'tab');

matlabbatch{2}.spm.stats.fmri_spec.sess(1).multi_reg = {confounds_file};
matlabbatch{2}.spm.stats.fmri_spec.sess(1).regress   = struct('name', {}, 'val', {});
matlabbatch{2}.spm.stats.fmri_spec.sess(1).hpf       = hpf_cutoff;

% -------------------------------------------------------------------------
% Condition onsets
% Each condition file: [onset duration] in seconds, one event per row
% -------------------------------------------------------------------------
for c = 1:numel(conditions)
    onsets = load(fullfile(onset_dir, [conditions{c}, '.txt']));

    matlabbatch{2}.spm.stats.fmri_spec.sess(1).cond(c).name     = conditions{c};
    matlabbatch{2}.spm.stats.fmri_spec.sess(1).cond(c).onset    = onsets(:, 1)';
    matlabbatch{2}.spm.stats.fmri_spec.sess(1).cond(c).duration = onsets(:, 2)';

    % No parametric or temporal modulation, and no orthogonalisation
    matlabbatch{2}.spm.stats.fmri_spec.sess(1).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{2}.spm.stats.fmri_spec.sess(1).cond(c).tmod = 0;
    matlabbatch{2}.spm.stats.fmri_spec.sess(1).cond(c).orth = 0;
end

%% ------------------------------------------------------------------------
% 3) ESTIMATE THE MODEL
% -------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat           = {fullfile(sub_output_dir, 'SPM.mat')};
matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;

%% ------------------------------------------------------------------------
% 4) DEFINE CONTRASTS
% -------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.spmmat = {fullfile(sub_output_dir, 'SPM.mat')};

% Indices of the nine task-related regressors used in DCM (ignore Fixation/Scream)
idx_dcm = [1 2 3 4 5 6 7 8 9];

% Effects of interest F-contrast (one regressor per condition of interest)
E = zeros(numel(idx_dcm), length(conditions));
for k = 1:numel(idx_dcm)
    E(k, idx_dcm(k)) = 1;
end

% F-contrast: all task-related conditions jointly
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name    = 'Effects of interest';
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = E;
matlabbatch{4}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

% -------------------------------------------------------------------------
% Simple T-contrasts: each condition of interest vs implicit baseline
% (Weights correspond to the first 9 task regressors; last 2 are Fixation/Scream.)
% -------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.consess{2}.tcon.name     = conditions{idx_dcm(1)};
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights  = [1 0 0 0 0 0 0 0 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{3}.tcon.name     = conditions{idx_dcm(2)};
matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights  = [0 1 0 0 0 0 0 0 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{4}.tcon.name     = conditions{idx_dcm(3)};
matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights  = [0 0 1 0 0 0 0 0 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{5}.tcon.name     = conditions{idx_dcm(4)};
matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights  = [0 0 0 1 0 0 0 0 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{5}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{6}.tcon.name     = conditions{idx_dcm(5)};
matlabbatch{4}.spm.stats.con.consess{6}.tcon.weights  = [0 0 0 0 1 0 0 0 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{6}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{7}.tcon.name     = conditions{idx_dcm(6)};
matlabbatch{4}.spm.stats.con.consess{7}.tcon.weights  = [0 0 0 0 0 1 0 0 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{7}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{8}.tcon.name     = conditions{idx_dcm(7)};
matlabbatch{4}.spm.stats.con.consess{8}.tcon.weights  = [0 0 0 0 0 0 1 0 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{8}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{9}.tcon.name     = conditions{idx_dcm(8)};
matlabbatch{4}.spm.stats.con.consess{9}.tcon.weights  = [0 0 0 0 0 0 0 1 0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{9}.tcon.sessrep  = 'none';

matlabbatch{4}.spm.stats.con.consess{10}.tcon.name    = conditions{idx_dcm(9)};
matlabbatch{4}.spm.stats.con.consess{10}.tcon.weights = [0 0 0 0 0 0 0 0 1 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

% -------------------------------------------------------------------------
% Differential T-contrasts (CS+ > CS-, CS- > CS+) within each phase
% -------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.consess{11}.tcon.name    = [conditions{idx_dcm(2)} '>' conditions{idx_dcm(3)}];
matlabbatch{4}.spm.stats.con.consess{11}.tcon.weights = [0  1 -1  0  0  0  0  0  0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{12}.tcon.name    = [conditions{idx_dcm(3)} '>' conditions{idx_dcm(2)}];
matlabbatch{4}.spm.stats.con.consess{12}.tcon.weights = [0 -1  1  0  0  0  0  0  0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{13}.tcon.name    = [conditions{idx_dcm(4)} '>' conditions{idx_dcm(5)}];
matlabbatch{4}.spm.stats.con.consess{13}.tcon.weights = [0  0  0  1 -1  0  0  0  0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{14}.tcon.name    = [conditions{idx_dcm(5)} '>' conditions{idx_dcm(4)}];
matlabbatch{4}.spm.stats.con.consess{14}.tcon.weights = [0  0  0 -1  1  0  0  0  0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{15}.tcon.name    = [conditions{idx_dcm(6)} '>' conditions{idx_dcm(7)}];
matlabbatch{4}.spm.stats.con.consess{15}.tcon.weights = [0  0  0  0  0  1 -1  0  0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{15}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{16}.tcon.name    = [conditions{idx_dcm(7)} '>' conditions{idx_dcm(6)}];
matlabbatch{4}.spm.stats.con.consess{16}.tcon.weights = [0  0  0  0  0 -1  1  0  0 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{16}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{17}.tcon.name    = [conditions{idx_dcm(8)} '>' conditions{idx_dcm(9)}];
matlabbatch{4}.spm.stats.con.consess{17}.tcon.weights = [0  0  0  0  0  0  0  1 -1 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{17}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{18}.tcon.name    = [conditions{idx_dcm(9)} '>' conditions{idx_dcm(8)}];
matlabbatch{4}.spm.stats.con.consess{18}.tcon.weights = [0  0  0  0  0  0  0 -1  1 zeros(1, 2)];
matlabbatch{4}.spm.stats.con.consess{18}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.delete = 0;

%% ------------------------------------------------------------------------
% 5) RUN PIPELINE
% -------------------------------------------------------------------------
spm_jobman('run', matlabbatch);

end