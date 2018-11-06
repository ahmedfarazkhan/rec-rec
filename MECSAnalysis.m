% MECSAnalysis.m
%
% Compare MECS from different A matrices

load('.\output\MECS_matrices.mat','MECS_matrices', 'As', 'Anames', 'rec_list', 'reg_list');

Anames = {'A_{human}', 'A_{human,sign}', 'A_{human,pval}', 'A_{genie}', 'A_{genie,sign}', 'A_{genie,p}', 'A_{ND}', 'A_{pc}'};


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
    f_name = sprintf('MECS_%s.png', char(Anames(atype)));
    %imwrite(MECS_curr, strcat('.\output\', f_name))
    saveas(im, f_name)
end

%% A matrices

for atype=1:numel(Anames)
    figure
    A_curr = squeeze(As(:, :, atype));
    im = imagesc(As(atype));
    mecs_title = sprintf('MECS with %s', char(Anames(atype)));
    title(mecs_title);
    xlabel('Rec');
    ylabel('Rec');
    set(gca, 'XTick', [1:1:15], 'XTickLabel', rec_list)
    set(gca, 'YTick', [1:1:15], 'YTickLabel', rec_list)
    f_name = sprintf('%s.png', char(Anames(atype)));
    saveas(im, f_name)
end

