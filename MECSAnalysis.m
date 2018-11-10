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