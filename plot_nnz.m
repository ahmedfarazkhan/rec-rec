function plot_nnz(params, fig_title)
%PLOT_NNZ Plot the number of nonzero parameters 
%   Inputs:
%       - params: 3D tensor
%       - title: 

figure(randi(100,1,1));
nnz_count = zeros(size(params,3), 1);

for p=1:size(params, 3)
    for fac=1:size(params(2))
        nnz_count(p) = nnz_count(p) + nnz(params(:, fac, p));
    end
end

scatter(1:size(params,3), nnz_count);
title(fig_title);

xlabel("Parameter");
ylabel("Nonzero Count");

filename = sprintf("[Regularizaton]%s.png", fig_title);
saveas(gcf, filename);
end

