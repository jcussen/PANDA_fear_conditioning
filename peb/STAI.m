%% ========================================================================
%  PEB PIPELINE FOR CLASS EFFECTS ON B-MATRIX (FEAR TASK, JULY 2025)
%
%  This script:
%    1) Loads a group of subject-level DCMs (GCM_*.mat)
%    2) Builds a between-subjects design matrix including a class regressor
%    3) Runs PEB on the B-matrix only
%    4) Performs Bayesian Model Reduction (BMR) + Bayesian Model Averaging (BMA)
%    5) Exports review plots and a table of class-related connections
%
%  Requirements:
%    - SPM25 on the MATLAB path
%    - GCM .mat file containing a cell array GCM (one DCM per subject)
%    - CSV file with a 'class' (or equivalent) column as group covariate
%
%  Intended as shareable analysis code to accompany a research paper.
%  Adjust the paths in the USER-ADJUSTABLE section to match your setup.
% ========================================================================

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS / PARAMETERS
% -------------------------------------------------------------------------
SPM_PATH   = '.../spm-main';

% Group DCM file (relative to rootDir)
spec       = '.../GCM.mat';

% Root directory for PEB / BMA outputs
rootDir    = '.../peb_results';

% CSV file providing subject-level covariates (e.g. latent class assignment)
design_csv = '.../stai_design.csv';

% ROI labels (must correspond to the ROI order in your DCMs)
roi        = {'hipp','postvmpfc','amygdala'};  

% Posterior probability threshold for reporting class effects
thr        = 0.95;

%% ------------------------------------------------------------------------
%  SETUP
% -------------------------------------------------------------------------
addpath(SPM_PATH);

spm('Defaults','fMRI');
spm_jobman('initcfg');
clearvars -except spm*   % keep SPM configuration; clear other variables

% Output directory derived from the GCM filename
outDir = fullfile(rootDir, erase(spec, '.mat'));
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created %s\n', outDir);
end

% Load group of subject DCMs
GCM = load(spec).GCM;   % cell array of subject-level DCM structures

%% ------------------------------------------------------------------------
% 1. DESIGN MATRIX: INTERCEPT + CLASS REGRESSOR
% -------------------------------------------------------------------------
% The design matrix encodes group-level effects. Here:
%   - Column 1: intercept (group mean)
%   - Column 2: class regressor (e.g. latent class membership)
%
% Non-intercept columns are mean-centred.
T  = readtable(design_csv);

X  = [ones(height(T), 1), T{:, 2}];      % intercept + class column
X(:, 2:end) = X(:, 2:end) - mean(X(:, 2:end));  % mean-centre non-intercept

M  = struct('Q', 'all', ...
            'X', X, ...
            'Xnames', {{'mean', 'class'}});

%% ------------------------------------------------------------------------
% 2. PEB ESTIMATION (B-MATRIX ONLY) + BMR / BMA
% -------------------------------------------------------------------------
% Estimate PEB model over B-matrix parameters, then perform BMR/BMA.
[PEB, RCM] = spm_dcm_peb(GCM, M, {'B'});
BMA        = spm_dcm_peb_bmc(PEB);

% Save PEB, BMA, and the design matrix for reproducibility
save(fullfile(outDir, 'PEB.mat'), 'PEB', 'RCM');
save(fullfile(outDir, 'BMA_search.mat'), 'BMA');
save(fullfile(outDir, 'M.mat'), 'M');   % optional: design matrix

%% ------------------------------------------------------------------------
% 3. REVIEW & EXPORT PEB PLOTS
% -------------------------------------------------------------------------
% Launch PEB review interface for visual diagnostics, then export figures.
spm_dcm_peb_review(BMA, GCM);
export_open_figs(outDir, 'PEB_review_plot');

%% ------------------------------------------------------------------------
% 4. DECODE & SAVE CLASS-SPECIFIC B-MATRIX PARAMETERS
% -------------------------------------------------------------------------
% BMA.Pp is concatenated across regressors (mean, class, etc.).
% Here we:
%   - Assume two regressors (mean, class)
%   - Take only the second half of parameters (class-specific effects)
%   - Keep parameters with Pp ≥ thr
pp   = full(BMA.Pp);
half = numel(pp) / 2;

keep = find(pp >= thr & (1:numel(pp))' > half);  % Pp ≥ thr AND in 2nd half

out  = {};
for p = keep(:)'
    % Wrap index back to the corresponding mean-effect parameter name
    lbl  = BMA.Pnames{mod(p-1, half) + 1};      % 21→1, 22→2, …
    v    = sscanf(lbl, 'B(%d,%d,%d)');          % [tgt src cond]

    cond = sprintf('%d', v(3));                 % condition index
    src  = roi{v(2)};                           % source region
    tgt  = roi{v(1)};                           % target region

    Ep   = full(BMA.Ep(p));                     % posterior mean
    Cp   = full(BMA.Cp(p, p));                  % posterior variance

    fprintf('B%s: %s → %s   Ep=%.4f  Pp=%.3f  Cp=%.3f\n', ...
            cond, src, tgt, Ep, pp(p), Cp);

    out(end+1, :) = {'B', cond, src, tgt, Ep, pp(p), Cp}; %#ok<AGROW>
end

writetable(cell2table(out, ...
          'VariableNames', {'matrix', 'condition', 'source', 'target', 'Ep', 'Pp', 'Cp'}), ...
          fullfile(outDir, 'class_connections.csv'));

%% ========================================================================
%                            Helper function
%% ========================================================================
function export_open_figs(folder, basename)
% EXPORT_OPEN_FIGS Save every open figure as a 300-dpi PNG.
%
%   EXPORT_OPEN_FIGS(folder, basename) saves all currently open figures to
%   the specified folder, with filenames of the form:
%       basename1.png, basename2.png, ...
%
    figs = findall(groot, 'Type', 'figure');
    [~, idx] = sort([figs.Number]);
    figs = figs(idx);  % sort in numeric order: 1, 2, 3, ...

    for k = 1:numel(figs)
        print(figs(k), ...
              fullfile(folder, sprintf('%s%d.png', basename, k)), ...
              '-dpng', '-r300');   % 300 dpi PNG
    end
end