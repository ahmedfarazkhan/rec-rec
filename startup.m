% rec-rec startup file

%addpath('/home/bic/akhan/Documents/MATLAB/spm12')


% Serotonin map data (Med Uni Wien)

% .hdr or .img?
paths_serotonin_maps = ["data/SerotoninMaps/5-HT1A/rWAY_HC36_mean.hdr", ...
"data/SerotoninMaps/5-HT1B/rP943_HC22_mean.hdr", ...
"data/SerotoninMaps/5-HT2A/rALT_HC19_mean.hdr", ...
"data/SerotoninMaps/5-HTT/rDASB_HC30_mean.hdr"];
names_serotonin_maps = ["5-HT1A", "5-HT1B", "5-HT2A", "5-HTT"];

% Reference: DKT Atlas
V_DKT = niftiread("data/DKTAtlas/corrected_DKTatlas.nii");

% Build matrix of regional receptor densities
N_recs = size(names_serotonin_maps, 2);
N_regs = 83; % DKT atlas regions

% Regional densities
M = zeros(N_recs, N_regs);

for rec=1:N_recs
    V_rec = spm_vol(char(paths_serotonin_maps(rec)));
    P_rec = spm_read_vols(V_rec);
    
    for reg=1:N_regs
        ind = find(V_DKT == reg);
        [X,Y,Z] = ind2sub(size(V_DKT), ind);
        
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

% Julich atlas registered to MNI space then resliced to DKT dimensions
V_Julich = niftiread("data/AtlasesMRI/rmni_Julich_Anatomy_v22c.nii");

% Brodmann atlas resliced to DKT dimensions
f_Brodmann = spm_vol(char("data/AtlasesMRI/rbrodmann.hdr"));
V_Brodmann = spm_read_vols(f_Brodmann);

% Area conversion
julich_area_conversion = xlsread('.\data\Human\receptor_data_area_atlas_correlation.xlsx');

%% commit to git repo
% Hybrid atlas
JB = zeros(size(V_DKT));



