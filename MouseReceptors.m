% MouseReceptors.m
%
% Receptor-receptor interactions in mice 
% using human-derived receptor-receptor A matrix
%
% Level 1: 
%   - One region
%   - Receptor-receptor interactions
%
% Ahmed, October 2018
% 
% Human rec indexing:
% [AMPA NMDA kainate muscimol flum cgp pire oxo damp epib praz [rx/uk14] dpat keta
% sch]
% X_j^WT, X_j^T
% S(t) = X_j^T - X_j^WT(t_0)
% dS(t)/dt =  A^Human S + Bu
%
% Region i, receptor j MECS matrix
% To do: better imputation, A signed?, mouse A, close to singular gram

% Same receptor order as human data, A matrix
rec_list = {'AMPA', 'MK80', 'KAIN', 'MUSC', 'FLUM', 'CGP5', 'PIRE', 'OXOT', 'DAMP', 'EPIB', 'PRAZ', 'UK14', 'KETA', 'DPAT', 'SCH2'};
reg_list = {'Au1_l', 'Au1_r', 'C_l', 'C_r', 'CM', 'CPu_l', 'CPu_r', 'Hip_l', 'Hip_r', 'M1_l', 'M1_r', 'RN', 'S1BF_l', 'S1BF_r', 'V1_l', 'V1_r', 'VPL_l', 'VPL_r', 'VPM_l', 'VPM_r'};

% Load human A matrix
human_data = importdata('.\data\Human\A_recp_vs_recp_matrix_human.mat');
A_human = human_data.A_human;
A_human_signs = human_data.A_human_signs;
A_human_signs = A_human_signs + (eye(size(A_human_signs)) - diag(sum(abs(A_human_signs))));
Pvalues = human_data.Pvalues;

% Load mice receptor densities
% [mouse, region, receptor]
ctrl_densities = importdata('.\data\Knockout mice_data/ctrl_densities.mat');
ko_densities = importdata('.\data\Knockout mice_data/ko_densities.mat');
ko_dims = size(ko_densities);

% Impute data using Trimmed Scores Regression (TSR)
missing_data = cat(1, ctrl_densities, ko_densities);
imputed_data = TSR(missing_data);

MECS_matrix_signs = zeros(ko_dims(2), ko_dims(3));

 
% % Calculating average densities
% X_ctrl = squeeze(mean(ctrl_densities,1)); M = mean(squeeze(mean(ctrl_densities,1))); S = std(squeeze(mean(ctrl_densities,1)));
% X_ctrl = (X_ctrl - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);
% X_ko   = squeeze(mean(ko_densities,1));
% X_ko   = (X_ko - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);
% 
% % Iterate over knockout mouse brain regions
% for reg=1:ko_dims(2)
%     %     X_ctrl = squeeze(ctrl_densities(1,reg,:));
%     %     X_ko = squeeze(ko_densities(1,reg,:));
%     z_t0 = X_ctrl(reg,:)'-X_ctrl(reg,:)';
%     z_tf = X_ko(reg,:)'-X_ctrl(reg,:)';
%     [MECS,U_MECS,MECS_times,B] = TargetControl(A_human_signs,z_tf,z_t0,0.0,0.5);
%     MECS_matrix_signs(reg,:) = MECS;
%     disp(MECS)
% end
% 
% %disp(MECS_matrix_signs)
% save('MECS_matrix_signs.mat','MECS_matrix_signs');
