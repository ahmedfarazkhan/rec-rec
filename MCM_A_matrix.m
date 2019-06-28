% Dynamics matrix for receptor-factor multifactorial model
% Factor-factor interactions as a function of receptors
% Ahmed Faraz Khan - May 2019

cd('/export02/data/Work/rec-rec')
addpath('/export02/data/Work/conn18b')

addpath(genpath('/export02/data/Work/MATLAB/'))
% Conflicts with stats.lasso
%addpath(genpath('/export02/data/Work/SPAMtoolbox/lar.m'))

load('before_A.mat');

% Todo: regression for gender + educationthen subtract that
% y_c = y - b_control*x_control

DIAG_ID = 1;

%% Order healthy subjects by age 

% X0 is averaged and normalized data for each subject
X0_combined = zeros(N_facs * N_regs, N_subjects);
individual_diags = zeros(N_subjects,1);

% Subject and sample index of data
ind_S = [];
min_times = []; % Of all samples
min_ages = zeros(N_subjects, 1); % Only baseline for each subject
for subj=1:N_subjects
    subj_times = global_data(subj).times;
    for sample=1:nnz(sum(global_data(subj).S, 1))
        
        ind_S = [ind_S; subj sample]; 
        
        % Remove cognitive and 0-rows before computing column mean times
        min_times = [min_times; min(subj_times(find(subj_times(1:N_facs, sample)), sample))];
        
        if sample == 1
            min_ages(subj, 1) = min(subj_times(find(subj_times(1:N_facs, sample)), sample));
        end
    end
    X0_combined(:, subj) = global_data(subj).X0;
    individual_diags(subj, 1) = global_data(subj).diags;
end

[sorted_min_times, ind_sort] = sort(min_times);
sorted_ind_S = ind_S(ind_sort, :);

% prev = 0;
% for ind=1:size(ind_S,1)
%     next = ind_S(ind,1);
%     if next < prev
%         disp("[Warning] Non-consecutive sample index");
%     end
%     prev = next;
% end

%% Regression for receptor-factor interactions + factor diffusion
N_samples = numel(min_times);
S_combined = zeros(N_facs * N_regs, N_subjects);
subject_diags = zeros(N_subjects, 1);

for i=1:N_subjects
    subj = sorted_ind_S(i, 1);
    sample = sorted_ind_S(i, 2);
    temp = global_data(subj).S;
    S_combined(:, i) = temp(:, sample);
    subject_diags(i,1) = global_data(subj).diags;
end


ind_healthy = find(subject_diags == DIAG_ID);
healthy_S = S_combined(:, ind_healthy);
sorted_mean_times_healthy = sorted_min_times(ind_healthy);

ind_healthy_X0 = find(individual_diags == DIAG_ID);
healthy_X0 = X0_combined(:, ind_healthy_X0);
healthy_ages = min_ages(ind_healthy_X0);
[sorted_ages, sorted_healthy_ind_X0] = sort(healthy_ages);
sorted_healthy_X0 = healthy_X0(:, sorted_healthy_ind_X0);
N_healthy = numel(ind_healthy_X0);


%% Regression for receptor-factor interactions + factor diffusion

% Treat 0s as missing data
healthy_X_reg = reshape(sorted_healthy_X0, N_regs, N_facs, numel(ind_healthy_X0));
%healthy_X_reg(healthy_X_reg == 0) = NaN;

%dS_dt = zeros((N_healthy -1)* N_regs, N_facs);
dS_dt = zeros(N_healthy - 1, N_regs, N_facs);
for i=1:N_healthy - 1
    dt = sorted_ages(i+1) - sorted_ages(i);
    if dt~=0
        for j=1:N_regs
            %dS_dt(((i-1)*N_regs) + j, :) = squeeze(healthy_X_reg(j,:,i+1) - healthy_X_reg(j,:,i)) / dt;
            dS_dt(i, j, :) = squeeze(healthy_X_reg(j,:,i+1) - healthy_X_reg(j,:,i)) / dt;
        end
    end
end

%% Factor interaction matrices for each region
A_ffs_flat = zeros(N_regs, N_facs * N_facs);

for reg=1:N_regs
    dS_dt_reg = squeeze(dS_dt(:, reg, :));
    A_ff_reg = cov(dS_dt_reg);
    A_ffs_flat(reg, :) = reshape(A_ff_reg, N_facs * N_facs, 1);
end

%% How much does each receptor affect the interaction between factors?

% y = Xb 
% y = A_ff_flat (N_facs * N_facs)
% X = Z0 (N_recs * N_regs)

bs = zeros(N_facs * N_facs, N_recs_J + N_recs_5ht + 1);

bs_lasso = zeros(N_facs * N_facs, N_recs_J + N_recs_5ht + 1);
visual_bs = zeros(N_facs * N_facs, N_recs_J + N_recs_5ht);

for ff=1:N_facs*N_facs

    X = [ones(N_regs, 1) Z0'];
    y = A_ffs_flat(:, ff);
    
    bs(ff, :) = regress(y,X);
    b = lasso(X,y);
    bs_lasso(ff,:) = b(:, 98);
    
    
end

for ff=1:N_facs*N_facs
    visual_bs(ff,:) = abs(bs(ff,2:end)) ./ max(abs(bs(ff,2:end)));
    
end



imagesc(visual_bs);
xlabel("Receptor Beta");
ylabel("Factor-Factor Interaction");

%% Simple linear regression for each factor and modality

% dS/dt = [sum_(n) alpha_i^n->m S_i^n] + spreading + inputs 
% dS/dt = AS + Bu
% a_i^(n->m)  = f(GE_i, NT-R_i)
% = beta_0^n->m + sum_r beta_r^(n->m) NT-R_(i,r)
% spreading = D^m ( sum_j C_(j->i)S_j^m   - sum_j C_(i->j) S_i^m )
% = D^m *sum_j C_j,i (S_j^m - S_i^m)

% Simplify
% m = n, i=/=j [same factor, different regions]
% A^(m,n)_(i,j) = beta_0^(m,n) + sum_j^N_RE beta_j^(m,n) RE_(j,i) - sum_k^N_roi C_i->k S^m

% m=/=n, i=/=j [different factors, regions]
% A^(m,n)_(i,j) = 0

% m=/n, i=j [different factors, same region]
% A^(m,n)_(i,j) =  b0^(m,n) + sum_j^(N_RE) beta_j^(m,n) RE_(j,i)

% Homogeneous population regression across regions for same factor to get
% fn->m (RE)

% y = Xb 
% X = [1 Re C*S]

%S = healthy_X_reg;
% Receptor interactions, factor interactions, spreading, no inputs
N_params = N_recs_J + N_recs_5ht + (N_facs - 1) + (N_regs - 1);
bs = zeros(N_regs, N_facs, N_params);
bs_lasso = zeros(N_regs, N_facs, N_params);
bs_ridge = zeros(N_regs, N_facs, N_params);
bs_elastic = zeros(N_regs, N_facs, N_params);
ks = 1e-5;%5e-3; for ridge 
bs_lar = zeros(N_regs, N_facs, N_params);
bs_bhs = zeros(N_regs, N_facs, N_params);

 for reg=1:N_regs

     ind_other_regs = setdiff(1:N_regs, reg);
     conn_reg = squeeze(sc(reg,ind_other_regs));

    for fac=1:N_facs
        y = squeeze(dS_dt(fac, :));

        % Every other factor and receptor
        ind_other_facs = setdiff(1:N_facs, fac);
        conn_S = conn_reg .* squeeze(S(ind_other_regs, fac, :))';
        inter_fac = squeeze(S(reg, ind_other_facs, :))';

        X = [ones(N_healthy, 1) repmat(squeeze(Z0(:, reg)'), N_healthy, 1) inter_fac conn_S];

        % Normalize X first?

        b = regress(y, X);
        bs(reg, fac, :) = b(2:end);

        [b_lasso, fitinfo] = lasso(X, y);
        [~, ind_lambda] = min(fitinfo.MSE);
        bs_lasso(reg, fac, :) = b_lasso(2:end, ind_lambda);

        % No constants
        b_ridge = ridge(y, X(:, 2:end), ks);
        bs_ridge(reg, fac, 2:end) = b_ridge(2:end, ind_lambda);

        [b_elastic, fitinfo] = lasso(X, y, 'alpha', 0.1);% 'kfold', 10);
        [~, ind_lambda] = min(fitinfo.MSE);
        bs_elastic(reg, fac, :) = b_elastic(2:end, ind_lambda);

        % Do not store path of vars
        %[b_lar, fitinfo] = lar(X, y, -10, false);


        [beta, b0, s2, t2, l2] = bhs(X, y, 5, 10, 1);% nsamples, burnin, thin);
        bs_bhs(reg, fac, :) = beta(2:end, ind_lambda);
        %bs_bhs(reg, fac, :) = bhs(X, y, nsamples, burnin, thin);

    end
end

%%


%% Plot nonzero parameters

plot_nnz(bs, 'Linear');
plot_nnz(bs_lasso, 'Lasso');
plot_nnz(bs_ridge, "Ridge");
plot_nnz(bs_elastic, "ElasticNet");
plot_nnz(bs_lar, "LeastAngle");
plot_nnz(bs_bhs, "BayesianHorseshoe");

%% Check which parameters are nonzero

figure(1);
title('Regression Parameters')
x = 1:size(bs, 3);
nnz_regress = zeros(size(bs,3), 1);



%% Assemble into A matrix

A = bs;
w = rand(N_facs, 1);
N_steps = 3;
S0 = squeeze(S(:, :, 1)); % Sorted by age
S_final = MCM_Rec(sc, sc, sc, A, S0, Z0', w, N_steps);

%%% CHECK SELF TERM CORRECTLY IMPLEMENTED IN MCM_REC.M?



%% Control (reachability)

A_2D = reshape(A, N_regs * N_facs, N_regs * N_facs);
% Todo add atype
Us_MECS = zeros(N_regs, N_facs, 1001);
Bs = zeros(N_regs, N_facs);

% Observable state

% Iterate over brain regions + factors (flatten)
for reg=1:N_regs
    for fac=1:N_facs
        other_facs = squeeze(S(reg, setdiff(1:N_facs, fac)));
        other_regs = squeeze(S(setdiff(1:N_regs, reg), fac))';
        state = [squeeze(Z(reg,:)) other_facs other_regs];

        C = zeros(state);
        C(N_recs_J + N_recs_5ht + fac) = 1;

        z_t0 = S0;
        z_tf = S_final;
        [MECS,U_MECS,MECS_times,B] = TargetControl(A, z_tf, z_t0, 0.0, 1.0);
        Us_MECS(reg,:,:) = U_MECS;
        MECS_matrices(reg, :) = MECS;
        Bs(reg, :) = B;
        %disp(MECS(atype, :, :))
        disp(size(Us_MECS))
        sprintf('Region %d', reg);
    end
end
