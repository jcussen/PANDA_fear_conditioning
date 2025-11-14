function dcm_spec_fear1_all(subject_id, version, spec)
% DCM_SPEC_FEAR1_ALL
%
% Specifies a bilinear task-driven DCM for the fear task (task 1) in SPM25.
% The model includes three regions (hippocampus, posterior vmPFC, amygdala),
% with fixed (A) and modulatory (B) connections and a single driving input (C).
%
% INPUTS
%   subject_id : numeric ID (e.g. 1, 2, ..., 999, 1000)
%   version    : string specifying analysis version (matches first-level code)
%   spec       : string tag identifying the DCM specification (used in filename)
%
% OUTPUT
%   Saves a DCM_*.mat file to the subjects first-level directory.
%
% Dependencies:
%   - SPM25
%   - First-level SPM.mat and VOI .mat files (see GLM and VOI extraction scripts)

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS
% -------------------------------------------------------------------------
% Root directory for first-level SPM output (must match GLM / VOI scripts)
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

% Condition inclusion vector:
% length = number of original conditions in the design (here 11),
% where 1 = include condition as input to DCM, 0 = ignore.
% (Here the first 9 task-related conditions are included; the last 2 are not.)
include = [1 1 1 1 1 1 1 1 1 0 0]';

% DCM analysis parameters
TR         = 3;        % Repetition time (s)
TE         = 0.04;     % Echo time (s)
nregions   = 3;        % Number of ROIs included in the model
nconditions = sum(include);    % Number of included experimental conditions

% Indexing for regions
hipp      = 1;
postvmpfc = 2;
amygdala  = 3;

% Subject ID string (EA1000 is a special case in this dataset)
if subject_id == 1000
    subject_id_str = '1000';
else
    subject_id_str = sprintf('%03d', subject_id);
end

% Subject-specific first-level directory (contains SPM.mat and VOIs)
output_dir = fullfile(OUTPUT_ROOT, version, folder);
spm_dir    = fullfile(output_dir, subject_id_str);

%% ------------------------------------------------------------------------
%  EFFECTIVE CONNECTIVITY MATRICES (A, B, C, D)
% -------------------------------------------------------------------------
% A-matrix: fixed (endogenous) connectivity
a = zeros(nregions, nregions);

% Limbic and cortical loops (bidirectional)
a(amygdala,  postvmpfc) = 1;
a(postvmpfc, amygdala)  = 1;

a(amygdala,  hipp)      = 1;
a(hipp,      amygdala)  = 1;

a(postvmpfc, hipp)      = 1;
a(hipp,      postvmpfc) = 1;

% C-matrix: driving inputs
% Here, the first included condition is treated as the main task regressor,
% driving the amygdala.
c = zeros(nregions, nconditions);
c(amygdala, 1) = 1;
% c(hipp,1)      = 1; % (alternative specification, not used)
% c(postvmpfc,1) = 1;

% D-matrix: nonlinear effects (not used in this model)
d = zeros(nregions, nregions, 0);

% B-matrix: modulatory effects (condition-specific changes in connectivity)
b = zeros(nregions, nregions, nconditions);

% For each condition from 2:nconditions, define symmetric modulatory effects
% between regions:
for cond = 2:nconditions
    % Amygdala ↔ Hippocampus
    b(hipp,     amygdala,  cond) = 1;
    b(amygdala, hipp,      cond) = 1;

    % Amygdala ↔ Posterior vmPFC
    b(amygdala,  postvmpfc, cond) = 1;
    b(postvmpfc, amygdala,  cond) = 1;

    % Hippocampus → Posterior vmPFC
    b(postvmpfc, hipp,      cond) = 1;

end

%% ------------------------------------------------------------------------
%  LOAD FIRST-LEVEL MODEL AND VOIs
% -------------------------------------------------------------------------
% Load first-level SPM structure
SPM = load(fullfile(spm_dir, 'SPM.mat'));
SPM = SPM.SPM;

% VOI .mat files (created by VOI extraction script)
% Each contains:
%   - xY: VOI structure (metadata)
%   - Y : extracted time series (single column)
voi_files = { ...
    fullfile(spm_dir, ['VOI_hipp_',      subject_id_str, '_1.mat']), ...
    fullfile(spm_dir, ['VOI_postvmpfc_', subject_id_str, '_1.mat']), ...
    fullfile(spm_dir, ['VOI_amygdala_',  subject_id_str, '_1.mat']) ...
    };

% Construct xY array for spm_dcm_specify
for r = 1:length(voi_files)
    XY      = load(voi_files{r});   % VOI.xY (struct), VOI.Y (time series)
    xY_temp = XY.xY;                % Copy metadata (name, Sess, etc.)
    xY_temp.xY = XY.Y;              % Store the time series in .xY field
    xY(r)      = xY_temp;           %#ok<AGROW>
end

%% ------------------------------------------------------------------------
%  SPECIFY AND SAVE DCM
% -------------------------------------------------------------------------
s             = struct();
s.name        = subject_id_str;
s.u           = include;                    % Included conditions (design inputs)
s.delays      = repmat(TR/2, 1, nregions); % Slice timing per region (seconds)
s.TE          = TE;
s.nonlinear   = false;   % Bilinear model
s.two_state   = false;   % One-state per region
s.stochastic  = false;   % Deterministic
s.centre      = true;    % Centre inputs
s.induced     = 0;
s.a           = a;
s.b           = b;
s.c           = c;
s.d           = d;

% Specify DCM using SPM
DCM = spm_dcm_specify(SPM, xY, s);

% Save DCM to subject directory
save(fullfile(spm_dir, ['DCM_', spec, '_', subject_id_str, '.mat']), 'DCM');

end