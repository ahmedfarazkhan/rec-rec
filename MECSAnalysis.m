% MECSAnalysis.m
%
% Compare MECS from different A matrices

load('.\output\MECS_matrices.mat','MECS_matrices', 'As', 'Anames', 'rec_list', 'reg_list');

N_REGS = 20;
N_RECS = 15;


%%
Rs = zeros(numel(Anames), numel(Anames));

for i=1:numel(Anames)
    for j=i:numel(Anames)
        v1 = As(:,:,i); v2 = As(:,:,j);
        ind = find(v1(:) & v2(:));
        temp = corrcoef(v1(ind),v2(ind));
        Rs(i,j) = temp(1,2); 
    end
end

Rs = Rs + Rs' - eye(size(Rs));

%%
imagesc(Rs)
title('A correlations');
xlabel('Method');
ylabel('Method');
set(gca, 'XTick', [1:1:numel(Anames)], 'XTickLabel', Anames)
set(gca, 'YTick', [1:1:numel(Anames)], 'YTickLabel', Anames)
colormap jet

%% Inter-species
human = As(:,:,10); mice = As(:,:,17); x = human(:); y = mice(:);
scatter(reshape(human, [1 15*15]), reshape(mice, [1 15*15]))
cftool

%% Intra-species
mice1 = As(:,:,18); mice2 = As(:,:,19); x1 = mice1(:); y1 = mice2(:);

scatter(reshape(mice1, [1 N_RECS*N_RECS]) , reshape(mice2, [1 N_RECS*N_RECS]))

cftool

%
%% MECS matrices

for atype=1:numel(Anames)
    figure
    
    MECS_curr = squeeze(MECS_matrices(atype, :, :));
    im = imagesc(MECS_curr);
    mecs_title = sprintf('MECS with %s', char(Anames(atype)));
    title(mecs_title);
    xlabel('Rec');
    ylabel('Region');
    set(gca, 'XTick', [1:1:15], 'XTickLabel', rec_list)
    set(gca, 'YTick', [1:1:20], 'YTickLabel', reg_list)
    colorbar
    f_name = sprintf('MECS_%s.png', char(Anames(atype)));
    %imwrite(MECS_curr, strcat('.\output\', f_name))
    saveas(im, strcat('.\output\', f_name))
end

%% A matrices

for atype=1:numel(Anames)
    figure
    A_curr = squeeze(As(:, :, atype));
    im = imagesc(A_curr);
    mecs_title = sprintf('Adjacency %s', char(Anames(atype)));
    title(mecs_title);
    xlabel('Rec');
    ylabel('Rec');
    set(gca, 'XTick', [1:1:15], 'XTickLabel', rec_list)
    set(gca, 'YTick', [1:1:15], 'YTickLabel', rec_list)
    colorbar
    f_name = sprintf('%s.png', char(Anames(atype)));
    saveas(im, strcat('.\output\', f_name))
end

%%


%% Pairwise differences in norms of A matrices

delta_A = zeros(numel(Anames), numel(Anames));

for i=1:numel(Anames)
    for j=i:numel(Anames)
        delta_A(i,j) = norm(As(:,:,i) - As(:,:,j));
    end
end

delta_A = delta_A + delta_A';

%% Receptor control energy rankings by A matrix and region
rec_ranks = zeros(numel(Anames), N_REGS, N_RECS, N_RECS);

for atype=1:numel(Anames)
    for reg=1:N_REGS
        %[X, rec_ranks(atype, reg, :)] =  sort( - MECS_matrices(atype, reg, :) ) ; % Descending order of E
    end
end

% 
%squeeze(rec_ranks(:, 1, :));

%% figure; stem(sum(mecs_pc))

%% Phenotype prediction

%  INPUT: 
%  BN : a BayesNet class object that represents the Bayes Net to predict
%       upon.  Must have these fields appropriate instantiated:
%    BN.pheno : string indicating which column to predict
%    BN.data : a matrix of data to predict
%    BN.cols : cell array of strings for column titles, or variable names.
%       must match BN.data and one element matches BN.pheno
%    BN.nodes : master list of NODES making up the bayes net, instantiated
%       by LearnParams()
%    BN.tree : ClusterSetTree instantiated by LearnParams()
%  VERBOSE : if true increases output. default = false;

A_m_bn = As(:,:,18);
%A_hm_bn = As(:,:,19); Predict human WT

%discrete = zeros(N_RECS); discrete(N_RECS) = 1;
%nodes = BNfromAdjMat(A_m_bn, discrete);

n_data_ctrl = size(ctrl_normed);
n_data_ko = size(ko_normed);
discvals = cat(2, zeros(1, n_data_ctrl(1)), ones(1, n_data_ko(1)));

data = cat(1, ctrl_normed, ko_normed);

disc    = zeros(N_RECS,1);
disc(N_RECS + 1) = 1;
for i = 1:N_RECS
    Signal_names{i,1} = char(rec_list(i)); 
    Signal(:,i)       = data(:,i); 
end
Signal(:,N_RECS+1) = cat(2, zeros(1, n_data_ctrl(1)), ones(1, n_data_ko(1))); %bp4_pff;
Signal_names{N_RECS+1,1} = 'KO_or_WT';
pheno                 = 'KO_or_WT';

priorPrecision.nu     = 100;
priorPrecision.alpha  = 100;
priorPrecision.sigma2 = 1;
priorPrecision.maxParents = N_RECS + 1;

analysis_title='PhenoPredict';

verbose = 1;

discrete = zeros(1, N_RECS+1); discrete(N_RECS+1) = 1;

%%


nodes_from_adj = BNfromAdjMat(A_m_bn, discrete, Signal_names);

[tree, nodes] = LearnParams(nodes_from_adj, '', Signal, Signal_names, priorPrecision);

%%
fprintf(1,'Learning Most Predictive Network Structure for %s\n', analysis_title);
[MBNet, FullBN, outstats] = LearnStructure(Signal, Signal_names, pheno, priorPrecision, [analysis_title,'-net']);
% fprintf(1,'Learning Network Parameters\n');
% 
% thresh = 0.5;
% A_thresh=zeros(N_RECS, N_RECS);
% low_ind = find(abs(A_m_bn)<thresh);
% high_ind = find(abs(A_m_bn)>=thresh);
% A_thresh(low_ind) = 0;
% A_thresh(high_ind) = 1;
% A_thresh = A_thresh - diag(diag(A_thresh));
% figure;
% imagesc(A_thresh);
% figure;
% imagesc(FullBN.adjmat)


nodes_from_FullBN_adj = BNfromAdjMat(FullBN.adjmat, discrete, Signal_names);

fprintf(1,'BN from Adj\n');
[tree, nodes] = LearnParams(nodes_from_FullBN_adj, '', Signal, Signal_names, priorPrecision);

FullBN.nodes = nodes;
FullBN.tree = tree;

fprintf(1,'Predicting on Training Data\n');
[acc, p, z] = PredictPheno(BN, verbose);

%[acc, p, z] = PredictPheno(tree, nodes, filename, pheno, data, colnames, verbose)
%[acc, p, z] = PredictPheno(BN, verbose);