%% ========================================================================
%  PEB / BMC SUMMARY PIPELINE FOR FEAR TASK (TASK 1)
%
%  This script:
%    1) Loads a group of subject-level DCMs (GCM_*.mat)
%    2) Computes basic diagnostics (variance explained per subject)
%    3) Runs PEB for A- and B-matrices
%    4) Applies a posterior-probability threshold and writes a summary table
%    5) Saves a small set of summary figures
%
%  Assumptions:
%    - SPM25 is installed and accessible.
%    - The specified GCM_*.mat file contains a cell array GCM.
%    - ROIs are in the order: hippocampus, posterior vmPFC, amygdala.
%
%  This script is intended as shareable analysis code to accompany a
%  research paper. Edit the paths and file names in the block below to
%  match your local setup.
% ========================================================================

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS / PARAMETERS
% -------------------------------------------------------------------------
SPM_PATH = '.../spm-main';

rootDir  = '.../peb_results';
spec     = '.../GCM.mat';

% ROI labels and posterior probability threshold for reporting
roi      = {'hipp','postvmpfc','amygdala'};
thresh   = 0.10;                             % posterior-probability cut-off

%% ------------------------------------------------------------------------
%  SETUP
% -------------------------------------------------------------------------
addpath(SPM_PATH);

spm('Defaults','fMRI');
spm_jobman('initcfg');
clearvars -except spm*  % keep SPM config in workspace, clear everything else

% Output directory for this particular GCM
outDir = fullfile(rootDir, erase(spec, '.mat'));

if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created %s\n', outDir);
end

% Load the group DCM cell array GCM from file
GCM  = load(spec).GCM;                      % bring GCM into workspace
nROI = numel(roi); %#ok<NASGU>

%% ------------------------------------------------------------------------
% 1) BASIC DIAGNOSTIC: VARIANCE EXPLAINED PER SUBJECT
% -------------------------------------------------------------------------
% spm_dcm_fmri_check returns a struct per subject; we extract the
% 'diagnostics' field (first entry) if present, else NaN.
vals = cellfun(@(x) pickfield(x, 'diagnostics', @nan), ...
               spm_dcm_fmri_check(GCM)); %#ok<*NASGU>

writematrix(vals, fullfile(outDir, 'variance_explained.txt'));

%% ------------------------------------------------------------------------
% 2) PEB FOR A- AND B-MATRICES
% -------------------------------------------------------------------------
% We fit separate PEB models for the A- and B-matrices, then extract
% connections that exceed the specified posterior probability threshold.
results = {};   % will grow row-by-row

results = run_peb(GCM, 'A', roi, thresh, results);
results = run_peb(GCM, 'B', roi, thresh, results);

writetable(cell2table(results, ...
           'VariableNames', {'matrix','condition','source','target','Ep','Pp','Cp'}), ...
           fullfile(outDir, 'connections.csv'));
fprintf('Saved connection table.\n');

%% ------------------------------------------------------------------------
% 3) SAVE A SMALL SET OF SUMMARY FIGURES
% -------------------------------------------------------------------------
% Exports selected open figures (if present) at 300 dpi.
save_figs(outDir, {'posterior_probabilities','driving_input_BMC'});

% Tidy up current figure handle (optional)
close(gcf);

%% ========================================================================
%                       ----  HELPER FUNCTIONS  ----
%% ========================================================================

function val = pickfield(s, f, def)
% PICKFIELD Return s.(f)(1) if present/non-empty, otherwise def()
%
%   val = PICKFIELD(s, f, def) returns the first element of field f in
%   struct s if it exists and is non-empty; otherwise it evaluates and
%   returns def(), where def is a function handle.
    if isfield(s, f) && ~isempty(s.(f))
        val = s.(f)(1);
    else
        val = def();
    end
end

% -------------------------------------------------------------------------
function results = run_peb(GCM, mat, roi, thr, results)
% RUN_PEB Run PEB for a specified connectivity matrix (A or B)
%
%   results = RUN_PEB(GCM, mat, roi, thr, results) fits a PEB model for the
%   matrix indicated by mat ('A' or 'B'), performs Bayesian model
%   comparison, and appends connections with posterior probability >= thr
%   to the cell array results.

    % Simple group design: intercept-only model
    M.X      = ones(numel(GCM), 1);
    [PEB, ~] = spm_dcm_peb(GCM, M, {mat});
    [BMA, ~] = spm_dcm_peb_bmc(PEB);

    % Retain connections above the posterior probability threshold
    keep = find(full(BMA.Pp) >= thr);

    for p = keep(:)'
        [cond, src, tgt] = decode_label(BMA.Pnames{p}, mat, roi);
        Ep  = full(BMA.Ep(p));
        Pp  = full(BMA.Pp(p));
        Cp  = full(BMA.Cp(p, p));

        fprintf('%s%s: %s → %s   Ep=%.4f  Pp=%.3f  Cp=%.3f\n', ...
                mat, cond, src, tgt, Ep, Pp, Cp);

        results(end+1, :) = {mat, cond, src, tgt, Ep, Pp, Cp}; %#ok<AGROW>
    end
end

% -------------------------------------------------------------------------
function [cond, src, tgt] = decode_label(lbl, mat, roi)
% DECODE_LABEL Translate parameter names (e.g. 'A(2,1)' or 'B(3,2,1)')
%              into condition, source, and target labels.
%
%   [cond, src, tgt] = DECODE_LABEL(lbl, mat, roi) returns:
%       cond : condition index as a string ('' for A-matrix parameters)
%       src  : name of source ROI
%       tgt  : name of target ROI

    if mat == "A"
        v    = sscanf(lbl, 'A(%d,%d)');
        cond = '';                      %#ok<NASGU>
    else
        v    = sscanf(lbl, 'B(%d,%d,%d)');
        cond = sprintf('%d', v(3));
    end

    src = roi{v(2)};
    tgt = roi{v(1)};
end

% -------------------------------------------------------------------------
function save_figs(outDir, names)
% SAVE_FIGS Export selected open figures to PNG files
%
%   SAVE_FIGS(outDir, names) saves figures in numeric order, mapping
%   names{1} to the second figure, names{2} to the third figure, etc.
%   (The first figure is often an SPM window or an overview plot.)

    figs = findall(groot, 'Type', 'figure');
    [~, idx] = sort([figs.Number]);
    figs = figs(idx);

    for k = 1:numel(names)
        if k + 1 <= numel(figs)
            exportgraphics(figs(k+1), ...
                fullfile(outDir, [names{k} '.png']), ...
                'Resolution', 300, 'BackgroundColor', 'white');
        end
    end
end