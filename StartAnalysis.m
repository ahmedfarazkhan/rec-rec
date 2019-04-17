% rec-rec startup file
%
% Create JB hybrid atlas 
% Filled-in Brodmann areas indexed starting from last Julich data region
% for easier reverse mapping

% NeuroPM_JB Hybrid Atlas Regions
%   reg_centroids: [x,y,z] coordinates of region centroids
%   reg_indices_lr: indices of 144 left/right split regions in NeuroPM-JB
%   reg_list: names of regions (left/right) for Julich + Brodmann regions
%       with data

% Factors
%   [CSF, A-beta, f-conn, glucose metabolism, GM density, tau]
% Receptors
%   ["5-HT1A", "5-HT1B", "5-HT2A", "5-HTT", 
% AMPA NMDA kainate muscimol flum cgp pire oxo damp epib praz rx dpat keta sch]

cd('/export02/data/Work/rec-rec')
addpath('/export02/data/Work/MATLAB/spm12')

rec_list = ["AMPA", "NMDA", "kainate", "muscimol", "flum", "cgp", "pire", ...
    "oxo", "damp", "epib", "praz", "rx", "dpat", "keta", "sch", "5-HT1A", ...
    "5-HT1B", "5-HT2A", "5-HTT"];
rec_types = ["glut", "glut", "glut", "GABA", "GABA", "GABA", "ACh", "ACh", ...
    "ACh", "ACh", "nor", "??", "ser", "ser", "dopamine", "ser", "ser", "ser", "ser"];


% Serotonin map data (Med Uni Wien)
paths_serotonin_maps = ["data/SerotoninMaps/5-HT1A/rWAY_HC36_mean.hdr", ...
"data/SerotoninMaps/5-HT1B/rP943_HC22_mean.hdr", ...
"data/SerotoninMaps/5-HT2A/rALT_HC19_mean.hdr", ...
"data/SerotoninMaps/5-HTT/rDASB_HC30_mean.hdr"];

% Reference: DKT Atlas
V_DKT = niftiread("data/DKTAtlas/corrected_DKTatlas.nii");


%% Build matrix of regional receptor densities using DKT regions
% N_recs = size(paths_serotonin_maps, 2);
% N_regs = 83; % DKT atlas regions
% 
% % Regional densities
% M = zeros(N_recs, N_regs);
% 
% for rec=1:N_recs
%     V_rec = spm_vol(char(paths_serotonin_maps(rec)));
%     P_rec = spm_read_vols(V_rec);
%     
%     for reg=1:N_regs
%         ind = find(V_DKT == reg);
%         %[X,Y,Z] = ind2sub(size(V_DKT), ind);
%         
%         % Collect densities form all voxels in region
%         P_recregs = P_rec(ind);
%         
%         % Mean receptor density in region
%         M(rec, reg) = mean(P_recregs); 
%    
%     end   
%     
%     M(rec, :) = M(rec,:) / sum(M(rec,:));
% end
% 
% %M = normr(M);
% imagesc(Z0_5ht)
% yticks([0.5, 1.5, 2.5, 3.5]);
% yticklabels(rec_list(16:end));
% xlabel("Regions (DTK Atlas)");
% ylabel("Receptors");
% save('.\output\M_PET','M')
% 
% % Find top regions
% [~,p] = sort(Z0_5ht(4,:),'descend');


%% 3D image of hybrid Julich + Brodmann atlas

% Julich atlas registered to MNI space (with reference brain extracted without eye, neck)
% using parameters from white matter image then resliced to DKT dimensions
% (NN interpolation throughout)
V_Julich = niftiread('data/AtlasesMRI/rrJulich_Anatomy_v22c.nii');
%load('data/AtlasesMRI/Julich_Anatomy_v22c_MPM.mat');

% Brodmann atlas resliced to DKT dimensions
f_Brodmann = spm_vol(char("data/AtlasesMRI/rbrodmann.hdr"));
V_Brodmann = spm_read_vols(f_Brodmann);

% Area conversion
[~, ~, region_map] = xlsread('data/Human/receptor_data_area_atlas_correlation.xls');
jb_area_conversion = cell2mat(region_map(2:size(region_map, 1), 3:5)); % [Julich, Brodmann]

% Used for extracting Julich data 
jb_shorthand = rmmissing(string(region_map(2:size(region_map, 1),7)));
jb_indices = string(region_map(2:size(region_map, 1),5));
j_shorthand = string(region_map(2:size(region_map, 1),1));

%reg_list = cat(1, jb_shorthand, brod_missing);
reg_brodmann = unique(V_Brodmann);
brod_in_JB = rmmissing(string(region_map(2:size(region_map, 1), 4)));
brod_missing = []; % Brodmann regions without Julich receptor data
for reg=1:size(reg_brodmann)
    if ~ismember(string(reg), brod_in_JB)
        if reg > 0
            brod_missing = cat(1, brod_missing, strcat("B", string(reg)));
        end
    end
end

reg_list = cat(1, jb_shorthand, brod_missing);

%% Create map for missing regions from Julich ROI maps 

% Fill JB atlas
V_JB = zeros(size(V_DKT));

j_roi_dir = 'data/Julich_ROI/Julich_ROI_extracted/r';
% [JB index, filename]
roi_map = cat(2, region_map(2:end,5), region_map(2:end,6));

% Find entries for which a ROI map exists
ind_maps = find(~ismissing(string(roi_map(:, 2))));

for reg=1:numel(ind_maps)
    % Get list of ROI maps for region and populate V_JB
    maps = split(string(roi_map(ind_maps(reg), 2)), ', ');
    for m=1:numel(maps)
        if ~isempty(maps(m))
            path = strrep(strcat(j_roi_dir, maps(m)), ' ', '');
            V_map = niftiread(char(strrep(path, '.gz', '')));
            ind_vox = find(V_map==1);
            V_JB(ind_vox) = cell2mat(roi_map(ind_maps(reg)));
        end
    end  
end


%% Create hybrid Julich + Brodmann atlas

% For every area in Julich & Brodmann atlases, assign an index
% Todo: intersection case?

for x=1:size(V_DKT,1)
    for y=1:size(V_DKT,2)
        for z=1:size(V_DKT,3)
            % Check if Brodmann or Julich codes are in the conversion table
            julich_code = V_Julich(x,y,z);
            brodmann_code = V_Brodmann(x,y,z);
                       
            [r_j,~] = find(jb_area_conversion(:,1) == julich_code);
                        
            if ~isempty(r_j) && julich_code > 0
                % Julich area found
                V_JB(x,y,z) = jb_area_conversion(r_j,3);
            else            
                [r_b,~] = find(jb_area_conversion(:,2) == brodmann_code);
                
                if ~isempty(r_b)
                    % Brodmann area found
                    V_JB(x,y,z) = jb_area_conversion(r_b,3); 
                    
                elseif V_Brodmann(x,y,z) % && ~V_JB(x,y,z)
                    % Use Brodmann area offset by JB index
                    % Only do this if Julich ROI map not used
                    % max(...) not size(unique(...)) because of NaN
                    V_JB(x,y,z) = brodmann_code + max(jb_area_conversion(:,3));
                    
               end
            end   
            
        end
    end
end


V = f_Brodmann;
V.fname = ['data/AtlasesMRI/NeuroPM_JB.nii'];
spm_write_vol(V,V_JB)



%% Get region center coordinates and neuroreceptor values from Julich regions

% Julich receptor data
[~, ~, raw_rec_data] = xlsread('data/Human/receptor_data.xls');
rec_list = string(raw_rec_data(1, 2:size(raw_rec_data,2)));
reg_list_Julich = string(raw_rec_data(2:size(raw_rec_data,1),1));
data_reg_indices = unique(double(jb_indices));

% rec_data: [areas, receptors]
rec_data = cell2mat(raw_rec_data(2:end, 2:end));
j_reg_names = string( string( raw_rec_data(2:end, 1)) );

% Geometric centroids of regions in each hemisphere
reg_centroids = zeros(max(V_JB(:)) * 2, 3); 
reg_rec_densities = zeros(max(V_JB(:)), numel(rec_list)); % NT-r densities
reg_indices_lr = horzcat(1:max(V_JB(:)), 1:max(V_JB(:))); % prune missing, indices in JB atlas

% Fill Julich receptor densities
for i_tar=1:size(data_reg_indices, 1)
    % Find regions with JB index (not contains but exact string match)
    ind_jb = find(strcmp(jb_indices, string(data_reg_indices(i_tar))));
    
    % [Special case 1] Multiple Julich regions to one JB index
    mean_data = zeros(1, size(rec_data,2));
    for j=1:size(ind_jb, 1)
        % Name of a region associated with index i
        % Julich index, not J-JB conversion table index
        reg_j = find(strcmp(j_reg_names, j_shorthand(ind_jb(j))));
        
        % [Special case 2] Multiple differentiated regions, same data
        if numel(reg_j) > 1
           reg_j = reg_j(1); 
        end
        mean_data = mean_data + rec_data(reg_j, :);
    end
    mean_data = mean_data / size(ind_jb, 1);
    
    reg_rec_densities(i_tar, :) = mean_data;  
       
end

%% Split hemispheres

V_JB = niftiread('data/AtlasesMRI/NeuroPM_JB.nii');
V_hemi_mask = niftiread('data/AtlasesMRI/rrAnatMask.nii');

V_JB_r = zeros(size(V_DKT));
V_JB_l = zeros(size(V_DKT));
V_JB_lr = zeros(size(V_DKT));

% [0, right, left]
mask_vals = unique(V_hemi_mask);

ind_r = find(V_hemi_mask == mask_vals(2));
ind_l = find(V_hemi_mask == mask_vals(3));

V_JB_r(ind_r) = V_JB(ind_r);
V_JB_l(ind_l) = V_JB(ind_l);

V_JB_lr(ind_r) = V_JB(ind_r) + max(V_JB_l(:));
ind_0 = find(V_JB_lr==max(V_JB_l(:)));
V_JB_lr(ind_0) = 0;
V_JB_lr(ind_l) = V_JB(ind_l);

V = f_Brodmann;
V.fname = ['data/AtlasesMRI/NeuroPM_JB_r.nii'];
spm_write_vol(V,V_JB_r)

V = f_Brodmann;
V.fname = ['data/AtlasesMRI/NeuroPM_JB_l.nii'];
spm_write_vol(V,V_JB_l)

V = f_Brodmann;
V.fname = ['data/AtlasesMRI/NeuroPM_JB_lr.nii'];
spm_write_vol(V,V_JB_lr)

%% Fill remaining with NaN
ind = find(reg_rec_densities==0, 1, 'first');
reg_rec_densities_nan = reg_rec_densities;
reg_rec_densities_nan(ind:max(V_JB(:)), :) = NaN(max(V_JB(:)) - ind + 1, numel(rec_list));

reg_rec_densities = cat(1, reg_rec_densities, reg_rec_densities);
reg_rec_densities_nan = cat(1, reg_rec_densities_nan, reg_rec_densities_nan);

% Find region centroids
for i_tar=1:max(V_JB(:))
    ind_V = find(ismember(V_JB_l, i_tar));
    [X,Y,Z] = ind2sub(size(V_JB_l), ind_V);
    reg_centroids(i_tar, :) = mean([X,Y,Z]);   
    
    ind_V = find(ismember(V_JB_r, i_tar));
    [X,Y,Z] = ind2sub(size(V_JB_r), ind_V);
    reg_centroids(i_tar+max(V_JB(:)), :) = mean([X,Y,Z]);   
end

% First NaN removal - atlas regions
% Remove NaN from indices and data - repeated Brodmann-filled regions
ind_nan = find(isnan(reg_centroids));
[ind_del, ~] = ind2sub(size(reg_centroids), ind_nan);
ind_del = unique(ind_del);

for i_tar=fliplr(1:size(ind_del))
   reg_centroids(ind_del(i_tar), :) = [];
   reg_rec_densities(ind_del(i_tar), :) = [];
   reg_rec_densities_nan(ind_del(i_tar), :) = [];
   reg_indices_lr(ind_del(i_tar)) = [];
end

save("data/AtlasesMRI/JB_rec_reg.mat", 'reg_centroids', 'reg_rec_densities', 'reg_indices_lr');

%% Normalize values for interpolation?

% Second NaN removal -> interpolation vs. data regions
[reg_rec_data, ind_nan] = rmmissing([reg_centroids reg_rec_densities_nan]); 

ind_interp = find(ind_nan == 1);
ind_regs_data = find(ind_nan==0);

% Normalized
% [X_train, X_tr_mu, X_tr_sigma] = zscore(reg_rec_data(:, 1:3));
% [Z_train, Z_tr_mu, Z_tr_sigma] = zscore(reg_rec_data(:, 4:end));
% Z_test = (reg_rec_densities - Z_tr_mu) ./ Z_tr_sigma;
% X_test = (reg_centroids - X_tr_mu) ./ X_tr_sigma;

% Non-normalized
X_train = reg_rec_data(:, 1:3);
X_tr_mu = 0;
X_tr_sigma = 1;
Z_tr_mu = 0;
Z_tr_sigma = 1;
Z_train = reg_rec_data(:, 4:end);
Z_test = reg_rec_densities;
X_test = reg_centroids;

N_regs = size(Z_test, 1);
N_recs_J = size(Z_test, 2);
N_regs_data = size(Z_train, 1);

% Consider only regions with data
%Z0_LOO_all = zeros(N_regs, N_regs, N_recs_J);
error_LOO = zeros(4, N_regs_data, N_recs_J);
acc_LOO = zeros(4, N_recs_J);
acc_prediction = zeros(4, N_recs_J); 


%% Kernel regression

disp("1) Kernel Regression")

[~, ~, sigma] = zscore(X_test);

% Actual interpolation
Z0_Julich_regressed = KernelRegression(X_train, Z_train, X_test, sigma);

% Error calculation
Z0_LOO_datareg = zeros(N_regs_data, N_recs_J); % LOO pred for each data reg
Z0_Julich_datareg = zeros(N_regs_data, N_recs_J); % Fit at data regions

for i_tar=1:numel(ind_regs_data)
    reg = ind_regs_data(i_tar);
    X_train_LOO = X_train(setdiff(1:N_regs_data, i_tar), :); 
    Z_train_LOO = Z_train(setdiff(1:N_regs_data, i_tar), :); 
    
    Z0_LOO = KernelRegression(X_train_LOO, Z_train_LOO, X_train, sigma);
    

    Z0_LOO_datareg(i_tar, :) = Z0_LOO(i_tar, :);
    Z0_Julich_datareg(i_tar, :) = Z0_Julich_regressed(reg, :);
    error_LOO(1, i_tar, :) = Z0_Julich_regressed(reg, :) - Z0_LOO(i_tar);   
end

acc_LOO(1, :) = diag(corr(Z0_Julich_datareg, Z0_LOO_datareg) .^2);



%% Scattered interpolation
% Delauney triangulation followed by adjacent vertex interpolation

disp("2) Scattered interpolation")

Z0_Julich_scat = zeros(size(Z_test));

for rec=1:N_recs_J
    F_rec = scatteredInterpolant(X_train, Z_train(:,rec), 'natural');
    Z0_Julich_scat(:, rec) = F_rec(X_test);
end

% Error calculation
Z0_LOO_datareg = zeros(N_regs_data, N_recs_J); % LOO pred for each data reg
Z0_Julich_datareg = zeros(N_regs_data, N_recs_J); % Fit at data regions

for i_tar=1:numel(ind_regs_data)
    reg = ind_regs_data(i_tar);
    X_train_LOO = X_train(setdiff(1:N_regs_data, i_tar), :);
    Z_train_LOO = Z_train(setdiff(1:N_regs_data, i_tar), :);
    
    Z0_LOO_reg = zeros(1, N_recs_J);
    
    for rec=1:N_recs_J
        F_rec_LOO = scatteredInterpolant(X_train_LOO, Z_train_LOO(:,rec), 'natural');
        rec_interp = F_rec_LOO(X_test);
        
        Z0_LOO_reg(1, rec) = rec_interp(i_tar);
        error_LOO(2, i_tar, rec) = Z0_Julich_scat(reg, rec) - rec_interp(reg);  
    end
    
    Z0_LOO_datareg(i_tar, :) = Z0_LOO_reg;
    Z0_Julich_datareg(i_tar, :) = Z0_Julich_scat(reg, :);
     
end

acc_LOO(2, :) = diag(corr(Z0_Julich_datareg, Z0_LOO_datareg) .^2);


%% Interpolate NT-R based on SC

disp("3) Structural connectivity interpolation")

% Regional structural (anatomical) connectivity
hcp_conn = load('data/tractography/hcp_fib_conn.mat'); % Averaged
%ntu_conn = load('data/tractography/ntu_src_conn.mat'); % Individual
sc = zscore(hcp_conn.connectivity);

% Interpolation
% 0 entries on regions without data 
Z0_Julich_sc_interp = DistInterpolation(Z_test, sc, ind_interp, Z_tr_mu, Z_tr_sigma);

%Validation error
Z0_LOO_datareg = zeros(N_regs_data, N_recs_J); % LOO pred for each data reg
Z0_Julich_datareg = zeros(N_regs_data, N_recs_J); % Fit at data regions

for i_tar=1:numel(ind_regs_data)
    reg = ind_regs_data(i_tar);
    
    Z_train_LOO = Z_train;
    Z_train_LOO(i_tar, :) = 0;

    sc_LOO = sc(ind_regs_data, ind_regs_data);
    
    Z0_LOO = DistInterpolation(Z_train_LOO, sc_LOO, i_tar, Z_tr_mu, Z_tr_sigma);
    
    Z0_LOO_datareg(i_tar, :) = Z0_LOO(i_tar, :);
    Z0_Julich_datareg(i_tar, :) = Z0_Julich_sc_interp(reg, :);
     
    error_LOO(3, i_tar, :) = Z0_Julich_sc_interp(reg, :) - Z0_LOO(i_tar, :);
end

acc_LOO(3, :) = diag(corr(Z0_Julich_datareg, Z0_LOO_datareg) .^2);

%% Compare with Euclidean distance

disp("4) Euclidean distance interpolation")

% Don't consider regions without data
euclid_dists = zscore(pdist2(X_test, X_test));
euclid_dists = euclid_dists + diag(diag(1));
euclid_conn = 1 ./ euclid_dists;
euclid_conn = euclid_conn  - diag(diag(euclid_conn));

Z0_Julich_euclid_interp = DistInterpolation(Z_test, euclid_conn, ind_interp, Z_tr_mu, Z_tr_sigma);

% Validation error
Z0_LOO_datareg = zeros(N_regs_data, N_recs_J); % LOO pred for each data reg
Z0_Julich_datareg = zeros(N_regs_data, N_recs_J); % Fit at data regions

for i_tar=1:numel(ind_regs_data)
    reg = ind_regs_data(i_tar);
    
    Z_train_LOO = Z_train;
    Z_train_LOO(i_tar, :) = 0;
    
    conn_LOO = euclid_dists(ind_regs_data, ind_regs_data);
    conn_LOO(i_tar, :) = 0; % [source, target]
    
    Z0_LOO = DistInterpolation(Z_train_LOO, conn_LOO, i_tar, Z_tr_mu, Z_tr_sigma);

    Z0_LOO_datareg(i_tar, :) = Z0_LOO(i_tar, :);
    Z0_Julich_datareg(i_tar, :) = Z0_Julich_euclid_interp(reg, :);
     
    error_LOO(4, i_tar, :) = Z0_Julich_euclid_interp(reg, :) - Z0_LOO(i_tar, :);
end

acc_LOO(4, :) = diag(corr(Z0_Julich_datareg, Z0_LOO_datareg) .^2);


%% Distance-based interpolation
% 
% disp("3,4) Structural connectivity and euclidean distance interpolation")
% 
% % Regional structural (anatomical) connectivity
% hcp_conn = load('data/tractography/hcp_fib_conn.mat'); % Averaged
% %ntu_conn = load('data/tractography/ntu_src_conn.mat'); % Individual
% sc = zscore(hcp_conn.connectivity);
% 
% % Don't consider regions without data
% euclid_dists = pdist2(X_test, X_test);
% euclid_dists = euclid_dists + diag(diag(1));
% euclid_conn = 1 ./ euclid_dists;
% euclid_conn = euclid_conn  - diag(diag(euclid_conn));
% 
% methods = ["SC", "ED"];
% conn_mats = [sc, euclid_conn];
% 
% Z0_Julich_interpolated = zeros(2, size(Z_test));
% 
% for method=1:size(methods, 2)
%     conn_mat = squeeze(conn_mats(1,:,:));
%     Z0_Julich_interpolated(m,:,:) = DistInterpolation(Z_test, conn_mat, ind_interp, Z_tr_mu, Z_tr_sigma);
% 
%     % Validation error
%     Z0_LOO_datareg = zeros(N_regs_data, N_recs_J); % LOO pred for each data reg
%     Z0_Julich_datareg = zeros(N_regs_data, N_recs_J); % Fit at data regions
% 
%     for i=1:numel(ind_regs_data)
%         reg = ind_regs_data(i);
% 
%         Z_train_LOO = Z_train;
%         Z_train_LOO(i, :) = 0;
% 
%         conn_LOO = conn_mat(ind_regs_data, ind_regs_data);
%         conn_LOO(i, :) = 0; % [source, target]
% 
%         Z0_LOO = DistInterpolation(Z_train_LOO, conn_LOO, i, Z_tr_mu, Z_tr_sigma);
% 
%         Z0_LOO_datareg(i, :) = Z0_LOO(i, :);
%         Z0_Julich_datareg(i, :) = squeeze(Z0_Julich_interpolated(m, reg, :));
% 
%         error_LOO(2+m, i, :) = Z0_Julich_datareg(i, :) - Z0_LOO(i, :);
%     end
% 
%     acc_LOO(2+m, :) = diag(corr(Z0_Julich_datareg, Z0_LOO_datareg) .^2);
% 
% end

%% Visualize leave-one-out error

figure(1);
A = imagesc(squeeze(error_LOO(1, :, :)));
title('Gaussian Regression LOO Error');
xlabel('Receptors');
ylabel('Regions');
saveas(A,'Gaussian Reg LOO Error','png')

figure(2);
A = imagesc(squeeze(error_LOO(2, :, :)));
title('Scattered Interpolation LOO Error');
xlabel('Receptors');
ylabel('Regions');
saveas(A,'Scattered Interpolation LOO Error','png')

figure(3);
A = imagesc(squeeze(error_LOO(3, :, :)));
title('Structural Conn LOO Error');
xlabel('Receptors');
ylabel('Regions');
saveas(A,'SC Interpolation LOO Error','png')

figure(4);
A = imagesc(squeeze(error_LOO(4, :, :)));
title('Euclidean Distance LOO Error');
xlabel('Receptors');
ylabel('Regions');
saveas(A,'ED Interpolation LOO Error','png')


figure(10);
A = imagesc(acc_LOO);
title('LOO Accuracy');
xlabel('Receptors');
ylabel('Methods');

set(gca,'Ytick',1:4,'YTickLabel',{'GR', 'Scat.', 'SC', 'ED'})
saveas(A,'LOO Accuracy','png')

%% Interpolated data

figure(5);
A = imagesc(Z0_Julich_regressed(ind_regs_data,:));
title('Gaussian Regression');
xlabel('Receptors');
ylabel('Regions');
saveas(A,'Gaussian Reg','png')

figure(6);
A = imagesc(Z0_Julich_scat(ind_regs_data,:));
title('Interpolation - Scattered');
xlabel('Receptors');
ylabel('Regions');
saveas(A,'Scattered Interpolation','png')

figure(7);
A = imagesc(Z0_Julich_sc_interp(ind_regs_data,:));
title('Interpolation - SC');
xlabel('Receptors');
ylabel('Regions');

saveas(A,'SC Interpolation','png')

figure(8);
A = imagesc(Z0_Julich_euclid_interp(ind_regs_data,:));
title('Interpolation - ED');
xlabel('Receptors');
ylabel('Regions');
saveas(A,'ED Interpolation','png')

%% Extract serotonin data from maps 

N_recs_5ht = size(paths_serotonin_maps, 2);
N_regs = size(reg_centroids, 1);

% Regional densities
Z0_5ht = zeros(N_regs, N_recs_5ht);

for rec=1:N_recs_5ht
    V_rec = spm_vol(char(paths_serotonin_maps(rec)));
    P_rec = spm_read_vols(V_rec);
    
    for reg=1:N_regs
        % Convert to JB reg
        ind = find(V_JB_lr == reg_indices_lr(reg));
        %[X,Y,Z] = ind2sub(size(V_DKT), ind);
        
        % Collect densities form all voxels in region
        P_recregs = P_rec(ind);
        
        % Mean receptor density in region
        Z0_5ht(reg,rec) = mean(P_recregs); 
    end   
    
    Z0_5ht(:, rec) = Z0_5ht(:, rec) / sum(Z0_5ht(:, rec));
end
%% Combine multimodal data
% [CSF, A-beta, functional connectivity, glucose metabolism, GM density, tau]

load('data/global_data_ADNI_standarization_3_atlas_JulichB.mat');

N_subjects = size(global_data, 2);
N_regs = size(Z0_Julich_scat, 1);
N_facs = size(global_data(1).X0, 1) / N_regs;
X0 = zeros(N_subjects, N_facs, N_regs);
diagnoses = zeros(1, N_subjects);

for subj=1:N_subjects
   % Check correct dimensions order
   temp = reshape(global_data(subj).X0, N_facs, N_regs);
   % [N_subjects, N_regions, N_factors]
   X0(subj, :, :) = permute(temp, [1 3 2]);
   diagnoses(1, subj) = global_data(subj).diags;
end

%% Combine receptor types

Z0_t = cat(2, Z0_Julich_regressed, Z0_5ht);
Z0 = Z0_t';

%% Estimate interfactor interactions

% Interaction coefficients as correlations between factors for all subjects
% in each region with data
% [CSF, A-beta, fMRI, glucose metabolism, GM density, tau]

N_factors = size(X0, 2);
alphas_corr = zeros(N_regs_data, N_factors, N_factors);

ind_healthy = find(diagnoses == 1);


for i_reg=1:N_regs_data
    reg = ind_regs_data(i_reg);
    for i_tar=1:N_factors
        
        src = squeeze(X0(ind_healthy, i_tar, reg));
        
        % Symmetric matrix, only fill half
        for j=i_tar:N_factors
            tar = squeeze(X0(ind_healthy, j, reg));
            
	
            alphas_corr(i_reg, i_tar, j) = corr(src, tar);
        end
    end
    
    % Zero diagonals
    alphas_corr(i_reg, :, :) = squeeze(alphas_corr(i_reg, :, :)) + squeeze(alphas_corr(i_reg, :, :))'- diag(ones(N_factors, 1));
end


%% Alternative alpha - regression coefficients not corr

alphas_regress = zeros(N_regs_data, N_factors, N_factors);
for i_reg=1:N_regs_data
    reg = ind_regs_data(i_reg);
    
    for i_tar=1:N_factors
        
        tar = squeeze(X0(ind_healthy, i_tar, reg));
        ind_srcs = setdiff(1:N_factors, i_tar);
        srcs = squeeze(X0(ind_healthy, ind_srcs, reg));
        b = regress(tar, srcs);
        
        % Zero diagonals
        alphas_regress(i_reg, ind_srcs, i_tar) = b;
    end
    
end

%% Neurovascular coupling for each region

mean_cbf = squeeze(mean(X0(ind_healthy, 1, ind_regs_data), 1));
mean_fALFF = squeeze(mean(X0(ind_healthy, 3, ind_regs_data), 1));

nvc = mean_cbf ./ mean_fALFF;

%% Predictability of neurovascular coupling

% Todo: Just data 

% From NT-R
rs_nt_nvc = zeros(size(Z0, 1), 1);
for rec=1:size(Z0, 1)
    r = corrcoef(Z0(rec, ind_regs_data)', nvc);
    rs_nt_nvc(rec, 1) = r(2,1);
end

r2s = rs_nt_nvc .^ 2;
indmax = find(r2s == max(r2s));

% Multiple correlation
R_NT = corr(Z0');
R2_NT = rs_nt_nvc' * (R_NT \ rs_nt_nvc);

sprintf("Best predicting NT is %s with %f%%, linear combination predicts %f%%",...
    rec_list(indmax), r2s(indmax) * 100, R2_NT * 100)

% From CBF
r2 = corrcoef(mean_cbf, nvc);
sprintf("CBF vs. NVC R^2 = %f", r2(2,1) ^ 2)

% From glucose metabolism
mean_glucose = squeeze(mean(X0(ind_healthy, 4, ind_regs_data), 1));
r2 = corrcoef(mean_glucose, nvc);
sprintf("Glucose metabolism vs. NVC R^2 = %f", r2(2,1) ^ 2)

% fMRI
r2 = corrcoef(mean_fALFF, nvc) ^ 2;
sprintf("Functional activity vs. NVC R^2 = %f", r2(2,1) ^ 2)

% NT-NVC Regression
b = regress(nvc, Z0(:, ind_regs_data)');
nvc_pred = zeros(size(nvc));
for reg=1:N_regs_data
    
    nvc_pred(reg, 1) = Z0(:, ind_regs_data(reg))' * b;
end

figure(100);
title('Neurovascular Coupling')
x = 1:numel(nvc);
hold on
scatter(x,nvc)
scatter(x,nvc_pred)
xlabel("Region");
ylabel("NVC");
legend("Actual", "Predicted");
hold off


figure(101);
title('NVC Error')
x = 1:numel(nvc);
hold on
scatter(x,nvc-nvc_pred)
scatter(x, nvc_pred)
xlabel("Region");
ylabel("NVC");
legend("Error", "Predicted");
hold off


%% Receptor-factor interactions (not regional)

% Todo: non-linear model?
XZ0 = cat(1, squeeze(mean(X0,1)), Z0); 
XZ0_normed = normalize(XZ0, 1);

% Ignore X-X interactions from this matrix
N_factors_receptors = size(XZ0, 1);
alphas_XZ = zeros(N_factors_receptors, N_factors_receptors);

for i_reg=1:N_regs_data
    reg = ind_regs_data(i_reg);
    
    for i_tar=1:N_factors_receptors
        
        tar = squeeze(XZ0(ind_healthy, i_tar, reg));
        ind_srcs = setdiff(1:N_factors_receptors, i_tar);
        srcs = squeeze(XZ0(ind_healthy, ind_srcs, reg));
        b = regress(tar, srcs);
        
        % Zero diagonals
        alphas_XZ(i_reg, ind_srcs, i_tar) = b;
    end
    
end



%% Visualize interpolation

% 
% for rec=1:size(reg_rec_densities, 2)
%     V_rec_data = zeros(size(V_JB));
%     V_rec_interp = zeros(size(V_JB));
%     
%     
%     data = 255 * Y_train(:, rec) / max(Y_train(:, rec));
%     
%     interp = 255 * interp_rec(:, rec) / max(interp_rec(:, rec));
%     
%     for reg=1:size(reg_rec_densities, 1)
%         % Set region to interpolated or actual centroid value
%         
%         reg_JB = reg_indices_lr(reg);
%         ind = find(V_JB == reg_JB);
% 
%         % Same data for both hemispheres
%         if reg <= size(Y_train, 1) / 2
%             V_rec_data(ind) = data(reg);
%         end
%         
%         % Left hemisphere
%         if reg <= max(V_JB(:))
%             ind_split = ind(ismember(ind, ind_l));
%         % Right hemisphere
%         else
%             ind_split = ind(ismember(ind, ind_r));
%         end
%         
%         V_rec_interp(ind_split) = interp(reg);
%     end
%     
%     V = f_Brodmann;
%     V.fname = [sprintf('data/AtlasesMRI/[Data]Rec%d.nii', rec)];
%     spm_write_vol(V,V_rec_data)
%     
%     V = f_Brodmann;
%     V.fname = [sprintf('data/AtlasesMRI/[Interp]Rec%d.nii', rec)];
%     spm_write_vol(V,V_rec_interp)
% end

