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
%   - Bayesian methods
%
% Level 3:
%   - Control
%
%   - ND without diagonals and A matrix, human with PC/ND,  
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


% Same receptor order as human data, A matrix
rec_list = {'AMPA', 'MK80', 'KAIN', 'MUSC', 'FLUM', 'CGP5', 'PIRE', 'OXOT', 'DAMP', 'EPIB', 'PRAZ', 'UK14', 'KETA', 'DPAT', 'SCH2'};
reg_list = {'Au1_l', 'Au1_r', 'C_l', 'C_r', 'CM', 'CPu_l', 'CPu_r', 'Hip_l', 'Hip_r', 'M1_l', 'M1_r', 'RN', 'S1BF_l', 'S1BF_r', 'V1_l', 'V1_r', 'VPL_l', 'VPL_r', 'VPM_l', 'VPM_r'};
N_CTRL = 7;
N_KO = 6;
N_REGS = numel(reg_list);
N_RECS = numel(rec_list);

%% Load imputed data

load('.\output\imputed_ctrl_densities.mat','ctrl_dens');
load('.\output\imputed_ko_densities.mat','ko_dens');

%%
% Calculate average densities from wild type (control) mice
 X_ctrl = squeeze(mean(ctrl_dens,1));
% M = mean(squeeze(mean(ctrl_dens,1)));
% S = std(squeeze(mean(ctrl_dens,1)));
% X_ctrl = (X_ctrl - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);

% Calculate average densities of knockout mice
 X_ko   = squeeze(mean(ko_dens,1));
% X_ko   = (X_ko - repmat(M,[size(X_ctrl,1) 1]))./repmat(S,[size(X_ctrl,1) 1]);

%%
% Using mice A matrices
%load('.\output\adjacency_matrices.mat', 'As', 'Anames');

%Anames = cat(2, Anames, {'A_{h,unsign}', 'A_{h,sign}', 'A_{h,pval}'});
%As = cat(3, As, A_h_unsign, A_h_sign, A_h_pval);

%%

% Initialize Minimum Control Energies matrix
%load('.\output\A_{m,bn,ko}.mat','BootsAdjMatMice_ko')
load('.\output\A_{m,bn,ctrl}.mat','BootsAdjMatMice_ctrl')

Anames = {'A_{m,bn,ctrl}_nonnormed_MECS', 'A_{m,bn,ctrl_trans}_nonnormed_MECS'};

As = cat(3, BootsAdjMatMice_ctrl(1:N_RECS, 1:N_RECS), BootsAdjMatMice_ctrl(1:N_RECS, 1:N_RECS)');%BootsAdjMatMice_ko(1:N_RECS, 1:N_RECS));

%%

MECS_matrices = zeros(numel(Anames), N_REGS, N_RECS);
% Todo add atype
Us_MECS = zeros(N_REGS, N_RECS, 1001);
Bs = zeros(numel(Anames), N_REGS, N_RECS);

%%
for atype=1:1%numel(Anames)
    
    A_curr = As(:,:,atype);
    
    % For Bayesian network adjacency matrices
    % Ones along diagonal (proportion of structural connectivity
    % plus random noise for MECS
    if atype %>= 18 
        A_curr = A_curr + diag(diag(ones(N_RECS))) + (randn(size(A_curr)) * 0.001);
    end
    
    % Iterate over knockout mouse brain regions
    for reg=1:N_REGS
        % Reachability, not controllability 
        % Transgenic densities are the final state
        z_t0 = X_ctrl(reg,:)'-X_ctrl(reg,:)';
        z_tf = X_ko(reg,:)'-X_ctrl(reg,:)';
        [MECS,U_MECS,MECS_times,B] = TargetControl(A_curr, z_tf, z_t0, 0.0, 1.0);       
        Us_MECS(reg,:,:) = U_MECS;
        MECS_matrices(atype, reg, :) = MECS;
        Bs(atype, reg, :) = B;
        %disp(MECS(atype, :, :))
        disp(size(Us_MECS))
        disp('A | Reg')
        disp([atype reg])
    end
    disp(atype)
end


save('.\output\MECS_BN_U','MECS_matrices', 'Us_MECS', 'As', 'Anames', 'rec_list', 'reg_list');




%%
save('.\output\MECS_matrices','MECS_matrices', 'As', 'Anames', 'rec_list', 'reg_list');
    


%%
A = squeeze(As(:,:,1));

T_STEP = 1/1000;
final_Zs = zeros(N_REGS, N_RECS, N_RECS); %

% Numerical forward 
% dX/dt = Ax + Bu
for reg=1:N_REGS
    B = squeeze(Bs(1, reg, :));
    for rec=1:N_RECS  % Only target a single receptor at a time
        Z = zeros(N_RECS, 1); 
        
        for t=1:size(Us_MECS(3) - 1)
            dX = (A * Z) + (B * Us_MECS(reg, rec, t));
            Z = Z + T_STEP * dX;
        end
        
        final_Zs(reg, rec, :) = squeeze(Z);
    end
end

%% Compare predictions to expected final states

final_Zs(1,1,:)
zs_tf = zeros(N_REGS, N_RECS);

for reg=1:N_REGS
    zs_tf(reg, :) = squeeze(X_ko(reg,:)'-X_ctrl(reg,:)');
end

x_axis = repmat(reshape(zs_tf, [1 N_REGS*N_RECS]), [1 N_RECS] ) ;
y_axis = reshape(final_Zs, [1 N_REGS*N_RECS*N_RECS]);
scatter(x_axis, y_axis)


for reg=1:N_REGS
    disp(corrcoef( reshape(zs_tf, [1 N_REGS*N_RECS]), squeeze(final_Zs(:,:,reg)) ))
end


%%
save('.\output\control_fw','final_Xs', 'Us_MECS', 'As', 'Anames', 'rec_list', 'reg_list');
    