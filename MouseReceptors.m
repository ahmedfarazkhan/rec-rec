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

N_META_COLS = 3; % Metadata in Excel file (group, animal, region)
N_MICE_KO = 6;
N_MICE_CTRL = 7;
N_RECEPTORS = 18;
N_REGIONS = 20;

% Load human A matrix
human_data = importdata('.\data\Human\A_recp_vs_recp_matrix_human.mat');
A_human = human_data.A_human;
A_human_signs = human_data.A_human_signs;
Pvalues = human_data.Pvalues;

% Load mice receptor densities
% [mouse, region, receptor]
ctrl_densities = importdata('.\data\Knockout mice_data/ctrl_densities.mat');
ko_densities = importdata('.\data\Knockout mice_data/ko_densities.mat');

ko_dims = size(ko_densities);

MECS_matrix = zeros(ko_dims(2), ko_dims(3));

% Iterate over knockout mouse brain regions
for reg=1:ko_dims(2)
    X_ctrl = squeeze(ctrl_densities(1,reg,:));
    X_ko = squeeze(ko_densities(1,reg,:)); 
    [MECS,U_MECS,MECS_times,B] = TargetControl(A_human,X_ctrl,X_ko,0.0,1.0);
    MECS_matrix(reg,:) = MECS;
    disp(MECS)
end

disp(MECS_matrix)
save('MECS_matrix.mat','MECS_matrix');

%[MECS,U_MECS,MECS_times,B] = TargetControl(A_human,X_ctrl,X_ko,0.0,1.0);
