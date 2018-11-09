% BNRecRec.m
%
% Adjacency matrix using Bayesian Networks (with bootstrapping)

rec_list = {'AMPA', 'MK80', 'KAIN', 'MUSC', 'FLUM', 'CGP5', 'PIRE', 'OXOT', 'DAMP', 'EPIB', 'PRAZ', 'UK14', 'KETA', 'DPAT', 'SCH2'};

mice_ctrl_dens = importdata('.\output\imputed_ctrl_densities.mat');
mice_ko_dens = importdata('.\output\imputed_ko_densities.mat');

n_data_ctrl = size(mice_ctrl_dens);
n_data_ko = size(mice_ko_dens);

mice_ctrl_dens = reshape(mice_ctrl_dens, [n_data_ctrl(1)*n_data_ctrl(2) n_data_ctrl(3)]);
mice_ko_dens = reshape(mice_ko_dens, [n_data_ko(1)*n_data_ko(2) n_data_ko(3)]);

data = cat(1, mice_ctrl_dens, mice_ko_dens);
n_data_ctrl = size(mice_ctrl_dens);
n_data_ko = size(mice_ko_dens);

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
nboots   = 10000;
verbose  = 1;
BootsAdjMat = BootstrapLearn(Signal, Signal_names, pheno, priorPrecision, nboots, 'ExhaustiveBN', verbose, BFTHRESH);

save('.\output\A_BN.mat','BootsAdjMat')

%%

figure
im = imagesc(BootsAdjMat);
title('Bayesian Network Adjacency');
xlabel('Rec');
ylabel('Rec');
set(gca, 'XTick', [1:1:15], 'XTickLabel', rec_list)
set(gca, 'YTick', [1:1:15], 'YTickLabel', rec_list)
colorbar
%f_name = sprintf('%s.png', char(Anames(atype)));
%saveas(im, strcat('.\output\', f_name))