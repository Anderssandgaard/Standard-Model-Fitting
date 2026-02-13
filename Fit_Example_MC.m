%% Standard Model Fitting: Intra-Axonal Monte Carlo Simulation
% -------------------------------------------------------------------------
% Author:   Anders Dyhr Sandgaard
% Date:     Jan 2026
% Project:  Standard_Model_Fitting
% Purpose:  Performs SMPR and ISO fitting on dMRI signals using GPU 
%           acceleration and spherical harmonics.
% -------------------------------------------------------------------------

% 1. Cleanup and Persistent Reset
clear; clc; close all;
clear comps_SM      ; % Crucial to reset persistent variables in your function
clear comps_SMPR    ; % Crucial to reset persistent variables in your function
clear comps_SMPR_ISO; % Crucial to reset persistent variables in your function

%% 1. Path Configuration
rootPath = setpath;

fprintf('Project Root: %s\n', rootPath);

%% 2. Data & Experiment Definitions
NAMES = {'SHAM_ipsi_25', 'SHAM_ipsi_49', 'SHAM_contra_25', 'SHAM_contra_49', ...
         'TBI_24_ipsi', 'TBI_24_contra', 'TBI_2_contra', 'TBI_2_ipsi'};

NAMES_PLOT = strrep(NAMES, '_', '-'); 

%% 3. Main Processing Loop
FIT = struct;
numDatasets = length(NAMES);

for i_name = 1:numDatasets
    currentName = NAMES{i_name};
    fprintf('\n--> Processing [%d/%d]: %s\n', i_name, numDatasets, currentName);
    
    % Construct path to "Signals_for fitting" (sibling to "Run")
    dataFile = fullfile(rootPath, 'Signals_for fitting', currentName);
    
    if ~exist([dataFile, '.mat'], 'file')
        warning('Datafile not found: %s. Skipping.', dataFile);
        continue;
    end
    load(dataFile);
    
    % Initialize results storage
    pout_SMPR     = [];
    pout_SMPR_ISO = [];
    
    % Metadata Extraction
    B0_count    = size(MRI.params.B0, 2);
    Bdir_count  = size(MRI.params.B0_Dir, 2);
    valid_b     = find(MRI.dMRI.bvals < inf);
    
    % Parameter Vectorization (Preparing for GPU)
    bvals = reshape(permute(repmat(MRI.dMRI.bvals(valid_b)', [1 5 5]), [2 3 1]), 1, []);
    bvecs = reshape(permute(repmat(MRI.dMRI.bvecs(valid_b, :), [1 1 5 5]), [2 3 4 1]), 3, []);
    
    TEs  = repmat(MRI.params.TEs', [1, length(MRI.params.DTEs), length(valid_b)]);
    DTEs = repmat(MRI.params.DTEs, [length(MRI.params.TEs), 1, length(valid_b)]);
    TE_vec  = TEs(:)';
    DTE_vec = DTEs(:)';
    
    for i_b0 = 1:B0_count
        for i_dir = 1:Bdir_count
            
            % Signal Preparation & Normalization
            Sig_raw    = squeeze(Signal.S1.STEdDI(:, :, i_b0, i_dir, valid_b));
            N_norm     = MRI.params.volfrac_IA * MRI.params.N_particles;
            Sig_norm   = fun.complexcat(1, Sig_raw(:)) ./ N_norm;
            
            % Setup Fitting Information (GPU Optimized)
            lmax          = 12;
            leb           = getLebedevSphere(3074);
            B0_vec        = MRI.params.B0_Dir(:, i_dir)';
            
            infos         = struct;
            infos.gpu     = 1;
            infos.lmax    = gpuArray(single(lmax));
            infos.Ymat    = gpuArray(single(SH.sphericalHarmonicsMatrix([leb.x, leb.y, leb.z], lmax)));
            infos.nw      = gpuArray(single(leb.w));
            infos.TE      = gpuArray(single(TE_vec * 1e-4));
            infos.dTE     = gpuArray(single(DTE_vec / 100));
            infos.dTER    = gpuArray(single(sin(acos([leb.x, leb.y, leb.z] * B0_vec')).^2 * (DTE_vec / 100)));
            infos.b       = gpuArray(single(bvals));
            infos.bgn     = gpuArray(single(([leb.x, leb.y, leb.z] * bvecs).^2 .* bvals));
            infos.w       = gpuArray(single(ones(1, size(Sig_norm, 2))));
            infos.N       = 1;

            % --- Fit 1: SMPR ---
            pars_init_smpr = single([1; 0.5; 0; 0; 0; 0]);
            [res_p_smpr] = fit.Intra.comps_SMPR(gpuArray(single(Sig_norm)), ...
                                                pars_init_smpr, infos, fit.options('imax', 1e3));
            pout_SMPR(i_b0, i_dir, :) = gather(res_p_smpr);
            
            S_fit_smpr = fun.Intra.comps_SMPR(res_p_smpr, gpuArray(single(Sig_norm)), infos);
            FIT.(currentName).RMSE_Sig(i_b0, i_dir) = gather(sum((S_fit_smpr - gpuArray(single(Sig_norm))).^2, 1));
            
            % --- Fit 2: ISO ---
            pars_init_iso = single([1; 0.5; 0; 0; 0]);
            [res_p_iso] = fit.Intra.comps_SMPR_ISO(gpuArray(single(Sig_norm)), ...
                                                   pars_init_iso, infos, fit.options('imax', 1e3));
            pout_SMPR_ISO(i_b0, i_dir, :) = gather(res_p_iso);
            
            S_fit_iso = fun.Intra.comps_SMPR_ISO(res_p_iso, gpuArray(single(Sig_norm)), infos);
            FIT.(currentName).RMSE_Sig_ISO(i_b0, i_dir) = gather(sum((S_fit_iso - gpuArray(single(Sig_norm))).^2, 1));
            
            % --- Metrics and Storage ---
            FIT.(currentName).Params_fit(i_b0, i_dir, :)     = pout_SMPR(i_b0, i_dir, :);
            FIT.(currentName).Params_fit_ISO(i_b0, i_dir, :) = pout_SMPR_ISO(i_b0, i_dir, :); 
            
            % Theoretical estimates
            gamma = MRI.params.gamma;
            B0_val = MRI.params.B0(i_b0);
            chi = MRI.params.chi_Bulk;
            
            FIT.(currentName).Est_A_theory(i_b0, i_dir) = (0.5 * gamma * B0_val) * chi * 1000 / (2*pi);
            FIT.(currentName).Est_B_theory(i_b0, i_dir) = chi * (-1/3 * gamma * B0_val * 1000 / (2*pi));
          
            % Clear specific GPU variables to keep memory lean
            wait(gpuDevice);
        end
    end
end

% Final Save (Optional but recommended)
% save(fullfile(rootPath, 'Results_Fitting.mat'), 'FIT', '-v7.3');
fprintf('\nProcessing complete for all datasets.\n');