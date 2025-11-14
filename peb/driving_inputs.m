%% DCM_PEB_DRIVING_INPUTS_FEAR1
%
% This script performs Bayesian model comparison (BMC) on alternative
% driving-input (C-matrix) configurations for a group DCM analysis of the
% fear task (task 1) using SPM25.
%
% Steps:
%   1) Load a group of subject-level DCMs (GCM_*.mat)
%   2) Fit a PEB model over subjects (effects on B and C matrices)
%   3) Define candidate models with different C-matrix structures
%   4) Perform BMC at the group level
%   5) Save model- and parameter-level posterior probabilities to CSV
%
% This script assumes that GCM_*.mat already exists (see dcm_load_fit_*).

%% ------------------------------------------------------------------------
%  USER-ADJUSTABLE PATHS
% -------------------------------------------------------------------------
% Path to SPM25
SPM_PATH = '.../spm-main';

% Root directory where PEB/DCM results are stored
rootDir  = '.../peb_results';

% Relative path to the GCM file (from rootDir)
spec = '.../GCM.mat';

%% ------------------------------------------------------------------------
%  GENERAL SETUP
% -------------------------------------------------------------------------
addpath(SPM_PATH);

% Initialise SPM
spm('Defaults', 'fMRI');
spm_jobman('initcfg');

% Output directory for this specific GCM
outDir = fullfile(rootDir, strrep(spec, '.mat', ''));

if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output folder: %s\n', outDir);
end

% Load group DCM cell array (GCM)
load(spec);   % expects variable GCM in the .mat file

%% ------------------------------------------------------------------------
%  PEB OVER B AND C MATRICES
% -------------------------------------------------------------------------
% Design matrix for PEB: intercept-only (group-average effects)
M_c      = [];
M_c.X    = ones(numel(GCM), 1);

% Run PEB on B and C parameters
[PEB_c, ~] = spm_dcm_peb(GCM, M_c, {'B', 'C'});

%% ------------------------------------------------------------------------
%  DEFINE CANDIDATE DRIVING-INPUT MODELS (C-MATRIX)
% -------------------------------------------------------------------------
% Region indices (must match the original DCMs)
hipp      = 1;
postvmpfc = 2;
amygdala  = 3;
striatum  = 4; %#ok<NASGU>  % defined for completeness, not used here

% Use the first subjects DCM as a template
DCM_full = GCM{1};

% If the DCM has already been estimated, remove previous priors
% so that changes to A/B/C are honoured
if isfield(DCM_full, 'M')
    DCM_full = rmfield(DCM_full, 'M');
end

% NOTE:
% Each model below modifies only the C-matrix (driving inputs),
% keeping A and B matrices identical across models.

%% Model 1: All regions driven (as in the original DCM_full)
% (DCM_full itself is used as the "all" model)

%% Model 2: Hippocampus only
DCM_model1      = DCM_full;
DCM_model1.c    = DCM_full.c;
DCM_model1.c(postvmpfc, 1) = 0;   % switch off postvmpfc in C
DCM_model1.c(amygdala,  1) = 0;   % switch off amygdala in C

%% Model 3: Posterior vmPFC only
DCM_model2      = DCM_full;
DCM_model2.c    = DCM_full.c;
DCM_model2.c(amygdala, 1)  = 0;   % switch off amygdala in C
DCM_model2.c(hipp,     1)  = 0;   % switch off hippocampus in C

%% Model 4: Amygdala only
DCM_model3      = DCM_full;
DCM_model3.c    = DCM_full.c;
DCM_model3.c(hipp,      1) = 0;   % switch off hippocampus in C
DCM_model3.c(postvmpfc, 1) = 0;   % switch off postvmpfc in C

%% Model 5: Hippocampus OFF (others as in DCM_full)
DCM_model4      = DCM_full;
DCM_model4.c    = DCM_full.c;
DCM_model4.c(hipp, 1)  = 0;

%% Model 6: Posterior vmPFC OFF
DCM_model5      = DCM_full;
DCM_model5.c    = DCM_full.c;
DCM_model5.c(postvmpfc, 1) = 0;

%% Model 7: Amygdala OFF
DCM_model6      = DCM_full;
DCM_model6.c    = DCM_full.c;
DCM_model6.c(amygdala, 1) = 0;

%% Collect candidate models into a GCM cell array
% One row per candidate model (not yet estimated).
GCM_c = { ...
    DCM_full, ...   % all regions driven
    DCM_model1, ... % hippocampus only
    DCM_model2, ... % postvmpfc only
    DCM_model3, ... % amygdala only
    DCM_model6, ... % hipp + postvmpfc (amygdala OFF)
    DCM_model5, ... % hipp + amygdala   (postvmpfc OFF)
    DCM_model4  ... % postvmpfc + amygdala (hipp OFF)
    };

%% ------------------------------------------------------------------------
%  BAYESIAN MODEL COMPARISON (BMC) AND OUTPUT
% -------------------------------------------------------------------------
% Run BMC on the PEB model using the candidate DCMs as templates
[BMA_c, ~] = spm_dcm_peb_bmc(PEB_c, GCM_c);   % BMA_c contains P and Pw

% Model labels and region names
models  = {'all', 'hipp', 'postvmpfc', 'amygdala', ...
           'hipp+postvmpfc', 'hipp+amygdala', 'postvmpfc+amygdala'};
regions = {'hipp', 'postvmpfc', 'amygdala'};

% 1. Model-level posterior probabilities
if isfield(BMA_c, 'P')
    Pm = BMA_c.P(:);      % column vector, one entry per candidate model
else
    error('BMA_c does not contain field ''P'' – cannot find model probabilities.');
end

% 2. Parameter-level posterior probabilities (driving input present?)
if isfield(BMA_c, 'Pw')
    Pp = BMA_c.Pw(:);     % one value per region (length 3)
else
    error('BMA_c does not contain field ''Pw'' – cannot find connection probabilities.');
end

% 3. Assemble & write tables
T1 = table(models(:), Pm, ...
           'VariableNames', {'model', 'model_probability'});

writetable(T1, fullfile(outDir, 'driving_inputs_models.csv'));
fprintf('Saved driving-input BMC model results to %s\n', ...
        fullfile(outDir, 'driving_inputs_models.csv'));

T2 = table(regions(:), Pp, ...
           'VariableNames', {'region', 'parameter_probability'});

writetable(T2, fullfile(outDir, 'driving_inputs_regions.csv'));
fprintf('Saved driving-input BMC parameter results to %s\n', ...
        fullfile(outDir, 'driving_inputs_regions.csv'));

