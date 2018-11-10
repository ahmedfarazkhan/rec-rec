%% Data Imputation
%
% Ahmed

% Load mice receptor densities
% [mouse, region, receptor]
ctrl_densities = importdata('.\data\Knockout mice_data/ctrl_densities.mat');
ko_densities = importdata('.\data\Knockout mice_data/ko_densities.mat');

% Impute data using Trimmed Scores Regression (TSR)
missing_data = cat(1, ctrl_densities, ko_densities);
imputed_data = TSR(missing_data);

% Note: 70% explained variance because of low sample size 
imputed_data = reshape(imputed_data, N_CTRL+N_KO, N_REGS, N_RECS);

% Split and reshape data - note Matlab reshapes using a different dimension
% order compared to Python
reshaped_data = permute(reshape(imputed_data, N_CTRL+N_KO, N_RECS, N_REGS), [1 3 2]);
ctrl_dens = reshaped_data(1:N_CTRL,:,:);
ko_dens = reshaped_data(N_CTRL+1:N_CTRL+N_KO,:,:);

save('.\output\imputed_ctrl_densities.mat','ctrl_dens');
save('.\output\imputed_ko_densities.mat','ko_dens');