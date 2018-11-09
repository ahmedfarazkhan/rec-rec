% BuildAdjacency.m
%
% Build adjacency matrix A from density distribution data
%
% Methods for A
%   - | GENIE3 |
%   - GENIE3 x sign corr
%   - (GENIE3 x sign corr)  x Pvalue(<0.05)
%   - deconvoluted corr
%   - partial corr (Matlab partialcorr(x,y,z))
%   - Maximize similarity between species (e.g. Bayesian, GENIE methods)
%   - Kullback-Leibler/other distance measure
%
% Ahmed, October 2018
% 
% Human rec indexing:
% [AMPA NMDA kainate muscimol flum cgp pire oxo damp epib praz [rx/uk14] dpat keta
% sch]

% Same receptor order as human data, A matrix
rec_list = {'AMPA', 'MK80', 'KAIN', 'MUSC', 'FLUM', 'CGP5', 'PIRE', 'OXOT', 'DAMP', 'EPIB', 'PRAZ', 'UK14', 'KETA', 'DPAT', 'SCH2'};
reg_list = {'Au1_l', 'Au1_r', 'C_l', 'C_r', 'CM', 'CPu_l', 'CPu_r', 'Hip_l', 'Hip_r', 'M1_l', 'M1_r', 'RN', 'S1BF_l', 'S1BF_r', 'V1_l', 'V1_r', 'VPL_l', 'VPL_r', 'VPM_l', 'VPM_r'};

N_REGS = numel(reg_list);
N_RECS = numel(rec_list);

ctrl_dens = importdata('.\output\imputed_ctrl_densities.mat');
n_data = size(ctrl_dens);

% Combine data from all regions for more points
rec_data = reshape(ctrl_dens, n_data(1) * n_data(2), n_data(3));

%% Data Imputation
% % Load mice receptor densities
% % [mouse, region, receptor]
% ctrl_densities = importdata('.\data\Knockout mice_data/ctrl_densities.mat');
% ko_densities = importdata('.\data\Knockout mice_data/ko_densities.mat');
% 
% % Impute data using Trimmed Scores Regression (TSR)
% missing_data = cat(1, ctrl_densities, ko_densities);
% imputed_data = TSR(missing_data);
% 
% % Note: 70% explained variance because of low sample size 
% imputed_data = reshape(imputed_data, N_CTRL+N_KO, N_REGS, N_RECS);
% 
% % Split and reshape data - note Matlab reshapes using a different dimension
% % order compared to Python
% reshaped_data = permute(reshape(imputed_data, N_CTRL+N_KO, N_RECS, N_REGS), [1 3 2]);
% ctrl_dens = reshaped_data(1:N_CTRL,:,:);
% ko_dens = reshaped_data(N_CTRL+1:N_CTRL+N_KO,:,:);
% 
% save('.\output\imputed_ctrl_densities.mat','ctrl_dens');
% save('.\output\imputed_ko_densities.mat','ko_dens');

%% Human

disp('Building human adjacency matrices')

% Load human A matrix
human_A = importdata('.\data\Human\A_recp_vs_recp_matrix_human.mat');
A_h_unsign = human_A.A_human;
A_h_unsign = A_h_unsign + (eye(size(A_h_unsign)) - diag(sum(abs(A_h_unsign))));

A_h_sign = human_A.A_human_signs;
A_h_sign = A_h_sign + (eye(size(A_h_sign)) - diag(sum(abs(A_h_sign))));

% Only siginifcant p-vals
Pvalues = human_A.Pvalues;
A_h_pval = A_h_sign;
p_indices = find(Pvalues>=0.05);  
A_h_pval(p_indices) = 0;

%% 
% Load Julich human data to generate A matrices
human_data = xlsread('.\data\Human\receptor_data.xls');
data_size = size(human_data);
human_data = human_data(:, 2:data_size(2));

% Genie3
[A_h_corr, p] = corrcoef(human_data);
A_h_genie = genie3(human_data);
A_h_genie = A_h_genie + (eye(size(A_h_genie)) - diag(sum(abs(A_h_genie))));

A_h_genie_sign = A_h_genie .* sign(A_h_corr);
A_h_genie_sign = A_h_genie_sign + (eye(size(A_h_genie_sign)) - diag(sum(abs(A_h_genie_sign))));

% Only significant (p<0.05) correlations
A_h_genie_pval = A_h_genie_sign; 
p_indices = find(p>=0.05);  
A_h_genie_pval(p_indices) = 0; % Todo: add diagonal?

% Network deconvolution on A matrix
A_h_ND_corr = ND(A_h_corr);
A_h_ND_genie = ND(A_h_genie);
A_h_pc = partialcorr(human_data);

%% Mice 
disp('Generating mice adjacency matrices')

% GENIE3
[A_m_corr, p] = corrcoef(rec_data);
A_m_genie = genie3(rec_data);
A_m_genie = A_m_genie + (eye(size(A_m_genie)) - diag(sum(abs(A_m_genie))));

A_m_genie_sign = A_m_genie .* sign(A_m_corr);
A_m_genie_sign = A_m_genie_sign + (eye(size(A_m_genie_sign)) - diag(sum(abs(A_m_genie_sign))));

% Only significant (p<0.05) correlations
A_m_genie_pval = A_m_genie_sign; 
p_indices = find(p>=0.05);  
A_m_genie_pval(p_indices) = 0; % Todo: add diagonal?

% Network deconvolution on A matrix
A_m_ND_corr = ND(A_m_corr);
A_m_ND_genie = ND(A_m_genie);
%A_ND = A_ND + (eye(size(A_ND)) - diag(sum(abs(A_ND))));
A_m_pc = partialcorr(rec_data);


%%
disp('Saving A matrices')

Anames = {'A_{h,unsign}', 'A_{h,sign}', 'A_{h,pval}', 'A_{h,corr}', 'A_{h,genie,unsign}', ... 
    'A_{h,genie,sign}', 'A_{h,genie,p}', 'A_{h,ND,corr}', 'A_{h,ND,genie}', ...
    'A_{h,pc}', 'A_{m,corr}', 'A_{m,genie,unsign}', 'A_{m,genie,sign}', 'A_{m,genie,p}', ...
    'A_{m,ND,corr}', 'A_{m,ND,genie}', 'A_{m,pc}'};

As = cat(3, A_h_unsign, A_h_sign, A_h_pval, A_h_corr,  A_h_genie, A_h_genie_sign, A_h_genie_pval, A_h_ND_corr, A_h_ND_genie, A_h_pc, A_m_corr, A_m_genie, A_m_genie_sign, A_m_genie_pval, A_m_ND_corr, A_m_ND_genie, A_m_pc);

% 'A_genie', 'A_genie_sign', 'A_genie_pval', 'A_ND_corr', 'A_ND_genie', 'A_pc',
save('.\output\adjacency_matrices.mat', 'As', 'Anames');

