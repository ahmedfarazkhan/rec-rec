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
%
%
% Ahmed, October 2018
% 
% Human rec indexing:
% [AMPA NMDA kainate muscimol flum cgp pire oxo damp epib praz [rx/uk14] dpat keta
% sch]

% Same receptor order as human data, A matrix
rec_list = {'AMPA', 'MK80', 'KAIN', 'MUSC', 'FLUM', 'CGP5', 'PIRE', 'OXOT', 'DAMP', 'EPIB', 'PRAZ', 'UK14', 'KETA', 'DPAT', 'SCH2'};
reg_list = {'Au1_l', 'Au1_r', 'C_l', 'C_r', 'CM', 'CPu_l', 'CPu_r', 'Hip_l', 'Hip_r', 'M1_l', 'M1_r', 'RN', 'S1BF_l', 'S1BF_r', 'V1_l', 'V1_r', 'VPL_l', 'VPL_r', 'VPM_l', 'VPM_r'};

ctrl_dens = importdata('.\output\imputed_ctrl_densities.mat');
n_data = size(ctrl_dens);

% Combine data from all regions for more points
rec_data = reshape(ctrl_dens, n_data(1) * n_data(2), n_data(3));

%%
% GENIE3
[A_corr, p] = corrcoef(rec_data);
A_genie = genie3(rec_data);
A_genie = A_genie + (eye(size(A_genie)) - diag(sum(abs(A_genie))));

A_genie_sign = A_genie .* sign(A_corr);
A_genie_sign = A_genie_sign + (eye(size(A_genie_sign)) - diag(sum(abs(A_genie_sign))));

% Only significant (p<0.05) correlations
A_genie_pval = A_genie_sign; 
p_indices = find(p>=0.05);  
A_genie_pval(p_indices) = 0; % Todo: add diagonal?

% Network deconvolution on A matrix?
A_ND = ND(A_genie);
A_ND = A_ND + (eye(size(A_ND)) - diag(sum(abs(A_ND))));
A_pc = partialcorr(rec_data);

save('.\output\A_mice.mat', 'A_genie', 'A_genie_sign', 'A_genie_pval', 'A_ND', 'A_pc');