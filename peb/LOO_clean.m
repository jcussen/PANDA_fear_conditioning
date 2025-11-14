%% ========================================================================
%  LOO CROSS-VALIDATION FOR A SINGLE B-MATRIX CONNECTION (FEAR TASK)
%
%  This script performs leave-one-out (LOO) cross-validation for a single
%  connection in the DCM B-matrix, using SPM25’s spm_dcm_loo.
%
%  It:
%    1) Loads a group of subject-level DCMs (GCM_*.mat)
%    2) Loads the corresponding group-level design matrix (M.mat)
%    3) Specifies one B-matrix parameter (labelled 'B(tgt,src,cond)')
%    4) Runs LOO cross-validation for that parameter
%    5) Saves the LOO figure (PNG) and numerical output (MAT)
%
%  Intended as shareable analysis code to accompany a research paper.
% ========================================================================

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS / PARAMETERS
% -------------------------------------------------------------------------
SPM_PATH = '.../spm-main';

% Group-level DCM file (relative to rootDir), containing cell array GCM
spec    = '.../GCM.mat';

% Root directory where PEB / GCM results and M.mat live
rootDir = '.../peb_results';

% ROI names (must match ROI ordering in the DCMs)
roiName = {'hipp','postvmpfc','amygdala'};   % indices 1–3

% Indices of regions (for clarity when choosing src/tgt)
hipp      = 1;
postvmpfc = 2;
amygdala  = 3;

% Condition index to test (must correspond to B-matrix condition index)
condID   = 2;

% Condition names (for naming the output files only)
condName = { ...
    'Task_Question', ...
    'CS+_earlyconditioning', 'CS-_earlyconditioning', ...
    'CS+_earlyextinction',   'CS-_earlyextinction', ...
    'CS+_lateconditioning',  'CS-_lateconditioning', ...
    'CS+_lateextinction',    'CS-_lateextinction', ...
    'Fixation', 'Scream'};

% Connection to test (source → target in B-matrix)
src = hipp;
tgt = postvmpfc;

%% ------------------------------------------------------------------------
%  INITIALISE SPM
% -------------------------------------------------------------------------
addpath(SPM_PATH);
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
clearvars -except spm*  % keep SPM configuration; clear everything else

%% ------------------------------------------------------------------------
%  HOUSE-KEEPING / LOADING DATA
% -------------------------------------------------------------------------
outDir = fullfile(rootDir, erase(spec, '.mat'));
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created %s\n', outDir);
end

% Load group-level DCMs and design matrix
GCM = load(spec).GCM;                         % cell array of subject DCMs
M   = load(fullfile(outDir, 'M.mat')).M;     % design-matrix struct used for PEB

%% ------------------------------------------------------------------------
%  LEAVE-ONE-OUT CROSS-VALIDATION
% -------------------------------------------------------------------------
% Label of the B-parameter to test:
%   B(tgt, src, condID)   – i.e. source → target for a specific condition
label = sprintf('B(%d,%d,%d)', tgt, src, condID);

% Run LOO CV for this parameter
[qE, qC, Q] = spm_dcm_loo(GCM, M, {label});

% -------------------------------------------------------------------------
%  SAVE FIGURE AND NUMERICAL OUTPUT
% -------------------------------------------------------------------------
base = fullfile(outDir, sprintf('LOO_%s_to_%s_%s', ...
    roiName{src}, roiName{tgt}, condName{condID}));

% Save current figure (300 dpi PNG)
print(gcf, [base '.png'], '-dpng', '-r300');

% Save LOO results (posterior means, covariances, and predictions)
save([base '.mat'], 'qE', 'qC', 'Q');

fprintf('Saved %s.[png|mat]\n', base);