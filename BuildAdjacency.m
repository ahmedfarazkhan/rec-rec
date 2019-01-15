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
%   - NORMALIZE Mice/human - but 2D?
%   - Randomize data order? and check phenotype predictions
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
n_data_ctrl = size(ctrl_dens);
ko_dens = importdata('.\output\imputed_ko_densities.mat');
n_data_ko = size(ko_dens);

% Combine data from all regions for more points
ctrl_dens = reshape(ctrl_dens, [n_data_ctrl(1)*n_data_ctrl(2) n_data_ctrl(3)]);
ko_dens = reshape(ko_dens, [n_data_ko(1)*n_data_ko(2) n_data_ko(3)]);

% Load Julich human data to generate A matrices
human_data = xlsread('.\data\Human\receptor_data.xls');
data_size = size(human_data);
human_data = human_data(:, 2:data_size(2));

% Normalize control data
ctrl_mean = mean(ctrl_dens); 
ctrl_std = std(ctrl_dens);
ctrl_normed = (ctrl_dens - repmat(ctrl_mean,[size(ctrl_dens,1) 1]))./repmat(ctrl_std,[size(ctrl_dens,1) 1]);

human_mean = mean(human_data); 
human_std = std(human_data);
human_normed = (human_data - repmat(human_mean,[size(human_data,1) 1]))./repmat(human_std,[size(human_data,1) 1]);

% Normalize KO data using control normalization parameters
ko_normed = (ko_dens - repmat(ctrl_mean,[size(ko_dens,1) 1]))./repmat(ctrl_std,[size(ko_dens,1) 1]);

% Data not normalized
ctrl_normed = ctrl_dens;
ko_normed = ko_dens;
human_normed = human_data;

% %% Random permutations of data rows
%   randseed = 1;
%   
%   ko_perm = randperm(size(ko_normed, 1));
%   ko_normed = ko_normed(ko_perm, :);
%  
%   ctrl_perm = randperm(size(ctrl_normed, 1));
%   ctrl_normed = ctrl_normed(ctrl_perm, :);
%   
%   human_perm = randperm(size(human_normed, 1));
%   human_normed = human_normed(human_perm, :);

%% Human - existing A matrices

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

%% Human - generated

% Genie3
[A_h_corr, p] = corrcoef(human_normed);
A_h_genie = genie3(human_normed);
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
A_h_pc = partialcorr(human_normed);

%% Mice 
disp('Generating mice adjacency matrices')

% GENIE3
[A_m_corr, p] = corrcoef(ctrl_normed);
A_m_genie = genie3(ctrl_normed);
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
A_m_pc = partialcorr(ctrl_normed);


%% Bayesian networks - Mice
disp('Bayesian networks')
data = cat(1, ctrl_normed, ko_normed);
n_data_ctrl = size(ctrl_normed);
n_data_ko = size(ko_normed);

N_signs = numel(rec_list);
disc    = zeros(N_signs,1);
disc(N_signs + 1) = 1;
for i = 1:N_signs
    Signal_names{i,1} = char(rec_list(i)); %['receptor_' num2str(i)];
    Signal(:,i)       = data(:,i); %Z1(i,:);
end
Signal(:,N_signs+1) = cat(2, zeros(1, n_data_ctrl(1)), ones(1, n_data_ko(1))); %bp4_pff;
Signal_names{N_signs+1,1} = 'KO_or_WT';
pheno                 = 'KO_or_WT';
priorPrecision.nu     = 100;
priorPrecision.alpha  = 100;
priorPrecision.sigma2 = 1;
priorPrecision.maxParents = N_signs + 1;
BFTHRESH = 0;
nboots   = 1000;
verbose  = 1;
BootsAdjMatMice = BootstrapLearn(Signal, Signal_names, pheno, priorPrecision, nboots, 'ExhaustiveBN', verbose, BFTHRESH);
save('.\output\A_{m,bn}.mat','BootsAdjMatMice')

A_m_bn = BootsAdjMatMice(1:N_RECS, 1:N_RECS);
%A_m_bn = A_m_bn + (randn(size(A_m_bn)) * 0.001); 
% Frequency of appearance in BNs interpretation
A_m_bn = A_m_bn - diag(diag(A_m_bn)) + eye(size(A_m_bn));
%% Human-mice Bayesian network

data = cat(1, human_normed, ctrl_normed, ko_normed);
n_data_human = size(human_normed);
n_data_ctrl = size(ctrl_normed);
n_data_ko = size(ko_normed);

clear Signal
N_signs = numel(rec_list);
disc    = zeros(N_signs,1);
disc(N_signs + 1) = 1;
for i = 1:N_signs
    Signal_names{i,1} = char(rec_list(i)); %['receptor_' num2str(i)];
    Signal(:,i)       = data(:,i); %Z1(i,:);
end

Signal(:,N_signs+1) = cat(2, zeros(1, n_data_human(1) + n_data_ctrl(1)), ones(1, n_data_ko(1))); %bp4_pff;
Signal_names{N_signs+1,1} = 'KO_or_WT';
pheno                 = 'KO_or_WT';
priorPrecision.nu     = 100;
priorPrecision.alpha  = 100;
priorPrecision.sigma2 = 1;
priorPrecision.maxParents = N_signs + 1;
BFTHRESH = 0;
nboots   = 1000;
verbose  = 1;
BootsAdjMatHM = BootstrapLearn(Signal, Signal_names, pheno, priorPrecision, nboots, 'ExhaustiveBN', verbose, BFTHRESH);
save('.\output\A_{h+m,bn}.mat','BootsAdjMatHM')

A_hm_bn = BootsAdjMatHM(1:N_RECS, 1:N_RECS);
%A_hm_bn = A_hm_bn %+ (randn(size(A_hm_bn)) * 0.001); 
% Frequency of appearance in BNs interpretation
A_hm_bn = A_hm_bn - diag(diag(A_hm_bn)) + eye(size(A_hm_bn));



%%
disp('Saving A matrices')

Anames = {'A_{h,unsign}', 'A_{h,sign}', 'A_{h,pval}', 'A_{h,corr}', 'A_{h,genie,unsign}', ... 
    'A_{h,genie,sign}', 'A_{h,genie,p}', 'A_{h,ND,corr}', 'A_{h,ND,genie}', ...
    'A_{h,pc}', 'A_{m,corr}', 'A_{m,genie,unsign}', 'A_{m,genie,sign}', 'A_{m,genie,p}', ...
    'A_{m,ND,corr}', 'A_{m,ND,genie}', 'A_{m,pc}', 'A_{m,bn}', 'A_{h+m,bn}'};

As = cat(3, A_h_unsign, A_h_sign, A_h_pval, A_h_corr,  A_h_genie, A_h_genie_sign, A_h_genie_pval, A_h_ND_corr, A_h_ND_genie, A_h_pc, A_m_corr, A_m_genie, A_m_genie_sign, A_m_genie_pval, A_m_ND_corr, A_m_ND_genie, A_m_pc, A_m_bn, A_hm_bn);

% 'A_genie', 'A_genie_sign', 'A_genie_pval', 'A_ND_corr', 'A_ND_genie', 'A_pc',
%fname = sprintf('.\output\adjacency_matrices_randperm_seed=%d.mat', randseed);
save('.\output\adjacency_matrices', 'As', 'Anames');
    

%% Bayesian networks - KO or WT
disp('Bayesian networks - KO or WT')

data = ctrl_normed; %ko_normed
n_data = size(data);


N_signs = numel(rec_list);
disc    = zeros(N_signs,1);
disc(N_signs + 1) = 1;
for i = 1:N_signs
    Signal_names{i,1} = char(rec_list(i)); %['receptor_' num2str(i)];
    Signal(:,i)       = data(:,i); %Z1(i,:);
end

subjects = 1:(n_data(1) / N_REGS);
Signal(:,N_signs+1) = repelem(subjects, N_REGS);
Signal_names{N_signs+1,1} = 'KO_or_WT';
pheno                 = 'KO_or_WT';
priorPrecision.nu     = 100;
priorPrecision.alpha  = 100;
priorPrecision.sigma2 = 1;
priorPrecision.maxParents = N_signs + 1;
BFTHRESH = 0;
nboots   = 1000;
verbose  = 1;

analysis_title='PhenoPredict';
BootsAdjMatMice_clear = BootstrapLearn(Signal, Signal_names, pheno, priorPrecision, nboots, 'ExhaustiveBN', verbose, BFTHRESH);
save('.\output\A_{m,bn,ctrl}.mat','BootsAdjMatMice_ctrl')

% fprintf(1,'Learning Most Predictive Network Structure for %s\n', analysis_title);
% [MBNet, FullBN, outstats] = LearnStructure(Signal, Signal_names, pheno, priorPrecision, [analysis_title,'-net']);
% 
% for i=1:N_RECS+1
%     FullBN.nodes(i)
% end
