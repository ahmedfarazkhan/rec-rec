% rec-rec startup file
%
% Create JB hybrid atlas 
% Filled-in Brodmann areas indexed starting from last Julich data region
% for easier reverse mapping

addpath('/export02/data/Work/MATLAB/spm12')

% Serotonin map data (Med Uni Wien)

% .hdr or .img?
paths_serotonin_maps = ["data/SerotoninMaps/5-HT1A/rWAY_HC36_mean.hdr", ...
"data/SerotoninMaps/5-HT1B/rP943_HC22_mean.hdr", ...
"data/SerotoninMaps/5-HT2A/rALT_HC19_mean.hdr", ...
"data/SerotoninMaps/5-HTT/rDASB_HC30_mean.hdr"];
names_serotonin_maps = ["5-HT1A", "5-HT1B", "5-HT2A", "5-HTT"];

% Reference: DKT Atlas
V_DKT = niftiread("data/DKTAtlas/corrected_DKTatlas.nii");


%% Build matrix of regional receptor densities using DKT regions
N_recs = size(names_serotonin_maps, 2);
N_regs = 83; % DKT atlas regions

% Regional densities
M = zeros(N_recs, N_regs);

for rec=1:N_recs
    V_rec = spm_vol(char(paths_serotonin_maps(rec)));
    P_rec = spm_read_vols(V_rec);
    
    for reg=1:N_regs
        ind = find(V_DKT == reg);
        %[X,Y,Z] = ind2sub(size(V_DKT), ind);
        
        % Collect densities form all voxels in region
        P_recregs = P_rec(ind);
        
        % Mean receptor density in region
        M(rec, reg) = mean(P_recregs); 
   
    end   
    
    M(rec, :) = M(rec,:) / sum(M(rec,:));
end

%M = normr(M);
imagesc(M)
yticks([0.5, 1.5, 2.5, 3.5]);
yticklabels(names_serotonin_maps);
xlabel("Regions (DTK Atlas)");
ylabel("Receptors");
save('.\output\M_PET','M')

% Find top regions
[~,p] = sort(M(4,:),'descend');


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
[~, ~, region_map] = xlsread('data/Human/receptor_data_area_atlas_correlation.xlsx');
jb_area_conversion = cell2mat(region_map(2:size(region_map, 1), 3:5)); % [Julich, Brodmann]

% Used for extracting Julich data 
jb_shorthand = string(region_map(2:size(region_map, 1),1));
jb_indices = string(region_map(2:size(region_map, 1),6));


%% Create map for missing regions from Julich ROI maps 

% Fill JB atlas
V_JB = zeros(size(V_DKT));

j_roi_dir = 'data/Julich_ROI/Julich_ROI_extracted/r';
% [JB index, filename]
roi_map = cat(2, region_map(2:end,5), region_map(2:end,7));

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


%% Create bybrid Julich + Brodmann atlas

%V_JB = zeros(size(V_DKT));
%V_B_full = zeros(size(V_DKT));

% For every area in Julich & Brodmann atlases, assign an index
% Todo: intersection case?
for x=1:size(V_DKT,1)
    for y=1:size(V_DKT,2)
        for z=1:size(V_DKT,3)
            % Check if Brodmann or Julich codes are in the conversion table
            julich_code = V_Julich(x,y,z);
            brodmann_code = V_Brodmann(x,y,z);
                       
            [r_j,~] = find(jb_area_conversion(:,1) == julich_code);
            
            if ~isempty(r_j)
                % Julich area found
                V_JB(x,y,z) = jb_area_conversion(r_j,3);
            else            
                [r_b,~] = find(jb_area_conversion(:,2) == brodmann_code);
                
                if ~isempty(r_b)
                    % Brodmann area found
                    V_JB(x,y,z) = jb_area_conversion(r_b,3); 
                    
                elseif V_Brodmann(x,y,z) & ~V_JB(x,y,z)
                    % Use Brodmann area offset by JB index
                    % Todo: fix gaps of Brodmann indices
                    % Only do this if Julich ROI map not used
                    V_JB(x,y,z) = brodmann_code + max(jb_area_conversion(:,3));
                end
            end   
            
        end
    end
end

V = f_Brodmann;
V.fname = ['data/AtlasesMRI/NeuroPM_JB.nii'];
spm_write_vol(V,V_JB)

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

%% Get region center coordinates and neuroreceptor values
% Repeat for each hemisphere
V_JB_split = [V_JB_l, V_JB_r];

% Julich receptor data
[~, ~, raw_rec_data] = xlsread('data/Human/receptor_data.xls');
rec_list = string(raw_rec_data(1, 2:size(raw_rec_data,2)));
combined_reg_list = string(raw_rec_data(2:size(raw_rec_data,1),1));
reg_list = cat(1, strcat(combined_reg_list, '_l'), strcat(combined_reg_list, '_r')); 

% rec_data: [areas, receptors]
rec_data = cell2mat(raw_rec_data(2:size(raw_rec_data,1), 2:size(raw_rec_data,2)));

% Geometric centroids of regions in each hemisphere
reg_centroids = zeros(max(V_JB(:)) * 2, 3); 
reg_rec_densities = zeros(max(V_JB(:)), numel(rec_list)); % NT-r densities
reg_indices_lr = horzcat(1:max(V_JB(:)), 1:max(V_JB(:))); % prune missing, indices in JB atlas

% Fill Julich receptor densities
for i=1:size(combined_reg_list,1)
    % Find region JB index (not contains but exact string match)
    % Index of conversion table
    ind_conv = find(strcmp(jb_shorthand, combined_reg_list(i)));
    
    % [Special case 1] Handle 1 region to multiple indices
    ind_reg_JB = double(unique(jb_indices(ind_conv)));
    
    % Region centroid
    % [Special case 2] Multiple regions to 1 index
    %ind_V = find(ismember(V_JB, ind_reg_JB));
    %[X,Y,Z] = ind2sub(size(V_JB), ind_V);
    %reg_centroids(ind_reg_JB, :) = mean([X,Y,Z]);
    reg_rec_densities(ind_reg_JB, :) = rec_data(i, :);  
end

%% Fill remaining with NaN
ind = find(reg_rec_densities==0, 1, 'first');
reg_rec_densities(ind:max(V_JB(:)), :) = NaN(max(V_JB(:)) - ind + 1, numel(rec_list));
reg_rec_densities = cat(1, reg_rec_densities, reg_rec_densities);

% Find region centroids
for i=1:max(V_JB(:))
    ind_V = find(ismember(V_JB_l, i));
    [X,Y,Z] = ind2sub(size(V_JB_l), ind_V);
    reg_centroids(i, :) = mean([X,Y,Z]);   
    
    ind_V = find(ismember(V_JB_r, i));
    [X,Y,Z] = ind2sub(size(V_JB_r), ind_V);
    reg_centroids(i+max(V_JB(:)), :) = mean([X,Y,Z]);   
end

% Remove NaN from indices and data - repeated regions
ind_nan = find(isnan(reg_centroids));
[ind_del, ~] = ind2sub(size(reg_centroids), ind_nan);
ind_del = unique(ind_del);

for i=fliplr(1:size(ind_del))
   reg_centroids(ind_del(i), :) = [];
   reg_rec_densities(ind_del(i), :) = [];
   reg_indices_lr(ind_del(i)) = [];
end

save("data/AtlasesMRI/JB_rec_reg.mat", 'reg_centroids', 'reg_rec_densities', 'reg_indices_lr');

%% Check region centroids
% 
% V_JB_centroids = V_JB;
% 
% for i=1:size(reg_centroids, 1)
%    ind_centroid = num2cell(round(reg_centroids(i,:)));
%    V_JB_centroids(ind_centroid{:}) = 0; 
% end
% 
% V = f_Brodmann;
% V.fname = ['data/AtlasesMRI/NeuroPM_region_centroids.nii'];
% spm_write_vol(V,V_JB_centroids)


%% Interpolate values

R = rmmissing([reg_centroids reg_rec_densities]); 
[Z_train, train_mu, train_sigma] = zscore(R(:, 1:3)) ; 
Y_train = R(:, 4:size(reg_rec_densities, 2) + 3); 

Z_pred = (reg_centroids - train_mu) ./ train_sigma;
Y_pred = zeros(size(reg_rec_densities));

%% Kernel regression
% Todo: kernel widths inversely proportionate to region size?
% Left/right
counts_lr = zeros(size(Y_train));
for reg=1:numel(counts_lr)
    counts_lr(reg) = numel(find(V_JB_lr==reg_indices_lr(reg)));
end
counts_lr = counts_lr / sum(counts_lr);

for rec=1:size(reg_rec_densities, 2)
%    [Md_rec, fit_info] = fitrkernel(Z_train, Y(:, rec),'KernelScale', 0.001);
%    Y_fit(:, rec) = predict(Md_rec, Z_test );
    Y_pred(:, rec) = KernelRegression(Z_train, Y_train(:, rec), Z_pred, 0.003);
end

% Compare known NT densities with regressed densities for the same regions
% Per hemisphere
N_reg_known = size(Z_train, 1);
N_reg_total = size(Z_pred, 1);
Y_compare = cat(1, Y_pred(1:N_reg_known/2, :), Y_pred((N_reg_total/2) + 1 : (N_reg_total/2) + (N_reg_known/2) , :));

errors = (Y_compare - Y_train) ./ Y_train;
error_sum = sum(errors, 2);
error_corr = corr(error_sum, counts_lr);

%% Scattered interpolation

% Create Delauney triangulation
%tri = delaunay(Z_train);
%trimesh(tri, Z_train(:,1), Z_train(:,2), Z_train(:,3))
%trisurf(tri, Z_train(:,1), Z_train(:,2), Z_train(:,3))

Y_interp = zeros(size(Y_pred));

for rec=1:size(Y_train,2)
    F_rec = scatteredInterpolant(Z_train, Y_train(:,rec), 'natural');
    Y_interp(:, rec) = F_rec(Z_pred);
end

%% Visualize interpolation


for rec=1:size(reg_rec_densities, 2)
    V_rec_data = zeros(size(V_JB));
    V_rec_interp = zeros(size(V_JB));
    
    
    data = 255 * Y_train(:, rec) / max(Y_train(:, rec));
    
    interp = 255 * Y_interp(:, rec) / max(Y_interp(:, rec));
    
    for reg=1:size(reg_rec_densities, 1)
        % Set region to interpolated or actual centroid value
        
        reg_JB = reg_indices_lr(reg);
        ind = find(V_JB == reg_JB);

        % Same data for both hemispheres
        if reg <= size(Y_train, 1) / 2
            V_rec_data(ind) = data(reg);
        end
        
        % Left hemisphere
        if reg <= max(V_JB(:))
            ind_split = ind(ismember(ind, ind_l));
        % Right hemisphere
        else
            ind_split = ind(ismember(ind, ind_r));
        end
        
        V_rec_interp(ind_split) = interp(reg);
    end
    
    V = f_Brodmann;
    V.fname = [sprintf('data/AtlasesMRI/[Data]Rec%d.nii', rec)];
    spm_write_vol(V,V_rec_data)
    
    V = f_Brodmann;
    V.fname = [sprintf('data/AtlasesMRI/[Interp]Rec%d.nii', rec)];
    spm_write_vol(V,V_rec_interp)
end

%% Combine multimodal data
