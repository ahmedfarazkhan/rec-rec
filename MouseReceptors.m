% MouseReceptors.m
%
% Receptor-receptor interactions in mice 
% using human-derived receptor-receptor A matrix
%
% Level 1: 
%   - One region
%   - Receptor-receptor interactions
%
% Level 2: Methods for A
%   - | genie |
%   - genie x sign corr
%   - (genie x sign corr)  x Pvalue(<0.05)
%   - deconvolutued corr
%   - partial corr (Matlab partialcorr(x,y,z))
% 
%
% Level 3:
%   - Control
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
% A: - species (mice), derivation methods
% Control: energy definition, time window

N_CTRL = 7;
N_KO = 6;
N_REGS = 20;
N_RECS = 15;

% Same receptor order as human data, A matrix
rec_list = {'AMPA', 'MK80', 'KAIN', 'MUSC', 'FLUM', 'CGP5', 'PIRE', 'OXOT', 'DAMP', 'EPIB', 'PRAZ', 'UK14', 'KETA', 'DPAT', 'SCH2'};
reg_list = {'Au1_l', 'Au1_r', 'C_l', 'C_r', 'CM', 'CPu_l', 'CPu_r', 'Hip_l', 'Hip_r', 'M1_l', 'M1_r', 'RN', 'S1BF_l', 'S1BF_r', 'V1_l', 'V1_r', 'VPL_l', 'VPL_r', 'VPM_l', 'VPM_r'};

% Load human A matrix
human_data = importdata('.\data\Human\A_recp_vs_recp_matrix_human.mat');
A_human = human_data.A_human;
A_human = A_human + (eye(size(A_human)) - diag(sum(abs(A_human))));

A_human_signs = human_data.A_human_signs;
A_human_signs = A_human_signs + (eye(size(A_human_signs)) - diag(sum(abs(A_human_signs))));

% Only siginifcant p-vals
Pvalues = human_data.Pvalues;
A_human_pval = A_human_signs;
p_indices = find(Pvalues>=0.05);  
A_human_pval(p_indices) = 0;

% Load mice receptor densities
% [mouse, region, receptor]
ctrl_densities = importdata('.\data\Knockout mice_data/ctrl_densities.mat');
ko_densities = importdata('.\data\Knockout mice_data/ko_densities.mat');

% Impute data using Trimmed Scores Regression (TSR)
missing_data = cat(1, ctrl_densities, ko_densities);
imputed_data = TSR(missing_data);

% Note: 70% explained variance because of low sample size 
imputed_data = reshape(imputed_data, N_CTRL+N_KO, N_REGS, N_RECS);

%% 

% Split and reshape data - note Matlab reshapes using a different dimension
% order compared to Python
reshaped_data = permute(reshape(imputed_data, N_CTRL+N_KO, N_RECS, N_REGS), [1 3 2]);
ctrl_dens = reshaped_data(1:N_CTRL,:,:);
ko_dens = reshaped_data(N_CTRL+1:N_CTRL+N_KO,:,:);

save('.\output\imputed_ctrl_densities.mat','ctrl_dens');
save('.\output\imputed_ko_densities.mat','ko_dens');

%%
% Calculate average densities from wild type (control) mice
X_ctrl = squeeze(mean(ctrl_dens,1));
M = mean(squeeze(mean(ctrl_dens,1)));
S = std(squeeze(mean(ctrl_dens,1)));
X_ctrl = (X_ctrl - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);

% Calculate average densities of knockout mice
X_ko   = squeeze(mean(ko_dens,1));
X_ko   = (X_ko - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);

%%
% MECS_matrix = zeros(ko_dims(2), ko_dims(3));
% 
% % Calculate average densities from wild type (control) mice
% X_ctrl = squeeze(mean(ctrl_dens,1)); 
% M = mean(squeeze(mean(ctrl_dens,1))); 
% S = std(squeeze(mean(ctrl_dens,1)));
% X_ctrl = (X_ctrl - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);
% 
% % Calculate average densities of knockout mice
% X_ko   = squeeze(mean(ko_dens,1));
% X_ko   = (X_ko - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);
% 
% % Iterate over knockout mouse brain regions
% for reg=1:ko_dims(2)
%     %     X_ctrl = squeeze(ctrl_dens(1,reg,:));
%     %     X_ko = squeeze(ko_dens(1,reg,:));
%     z_t0 = X_ctrl(reg,:)'-X_ctrl(reg,:)';
%     z_tf = X_ko(reg,:)'-X_ctrl(reg,:)';
%     [MECS,U_MECS,MECS_times,B] = TargetControl(A_human_signs,z_tf,z_t0,0.0,1.0);
%     MECS_matrix(reg,:) = MECS;
%     disp(MECS)
% end
% 
% save('.\output\MECS_matrix_signed.mat','MECS_matrix');


%%
% Using mice A matrices
load('.\output\A_mice.mat', 'A_genie', 'A_genie_sign', 'A_genie_pval', 'A_ND', 'A_pc');

Anames = {'A_{human}', 'A_{human,sign}', 'A_{human,pval}', 'A_{genie}', 'A_{genie,sign}', 'A_{genie,p}', 'A_{ND}', 'A_{pc}'};
As = cat(3, A_human, A_human_signs, A_human_pval, A_genie, A_genie_sign, A_genie_pval, A_ND, A_pc);

%%

% Initialize Minimum Control Energies matrix
ko_dims = size(ko_dens);
MECS_matrices = zeros(numel(Anames), ko_dims(2), ko_dims(3));

for atype=1:numel(Anames)
    
    A_curr = As(:,:,atype);
    %MECS_matrix = zeros(ko_dims(2), ko_dims(3));
    
    % Iterate over knockout mouse brain regions
    for reg=1:ko_dims(2)
        z_t0 = X_ctrl(reg,:)'-X_ctrl(reg,:)';
        z_tf = X_ko(reg,:)'-X_ctrl(reg,:)';
        [MECS,U_MECS,MECS_times,B] = TargetControl(A_curr, z_tf, z_t0, 0.0, 1.0);
        MECS_matrices(atype, reg, :) = MECS;
        %disp(MECS(atype, :, :))
        disp('A | Reg')
        disp([atype reg])
    end
    disp(atype)
%     imagesc(MECS_matrix)
%     title('Unsigned');
%     xlabel('Receptors');
%     ylabel('Regions');
%     set(gca, 'XTick', [1:1:15], 'XTickLabel', rec_list)
%     set(gca, 'YTick', [1:1:20], 'YTickLabel', reg_list)
    
end

save('.\output\MECS_matrices.mat','MECS_matrices', 'As', 'Anames', 'rec_list', 'reg_list');
    