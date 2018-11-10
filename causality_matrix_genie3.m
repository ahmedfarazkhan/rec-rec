function [A,Asigns,Pvalues_A,self_effects] = causality_matrix_genie3(data)
% in data, rows are observations, columns variables/nodes
% in A, Asigns & Pvalues, rows are target nodes, columns are direct
% modulators.
% in self_effects, each element i represents the fraction of the data in
% node i that couldn't be explained by the other nodes (we assume that this
% is representative of self effects).
%-----------------------------------------------------------------
% Yasser Iturria Medina, MNI, August. 2018.

% using Genie3:
[N_observ,N_nodes] = size(data);
for node_i = 1:N_nodes
    rec_names{node_i,1} = ['Node ' num2str(node_i)];
end
A      = genie3(data)'; % 
Asigns = A.*sign(corr(data)');
self_effects = 1 - sum(A,2);
% A = A + diag(self_effects); Asigns = Asigns + diag(self_effects);

% In the following graph, y-axis has target receptors, x-axis the direct modulators (excluding self-interactions for visual clarity)
figure; imagesc(A); colorbar; title('Local node-node causality (consistent across observations)')
xticklabel_rotate([1:length(rec_names)],45,rec_names,'interpreter','none'); colorbar; colormap Jet;
yticks([1:length(rec_names)]); yticklabels(rec_names);

% Creating Null distribution to test significance of causality
N_null = 5000;
A_perm = zeros([N_nodes N_nodes N_null]);
h = waitbar(0,'Performing permutations for significance ...');
for perm = 1:N_null
    waitbar(perm/N_null);
    rand_posc = round(1 + (N_observ - 1)*rand(N_observ,N_nodes));
    A_perm(:,:,perm) = genie3(data(rand_posc))';
end
close(h);

% Testing significance II
for i = 1:N_nodes
    for j = 1:N_nodes
        if i ~= j
            mean_null_ij = mean(A_perm(i,j,1:N_null));
            std_null_ij  = std(A_perm(i,j,1:N_null));
            Pvalues_A(i,j) = 1 - normcdf(A(i,j),mean_null_ij,std_null_ij);
        end
    end
end
alpha_forP = 0.05;
% In the following graph, y-axis has target nodes, x-axis the direct modulators (excluding self-interactions for visual clarity)
figure; imagesc((Pvalues_A < alpha_forP).*A); title('Significant node-node direct/causal links');
xticklabel_rotate([1:length(rec_names)],45,rec_names,'interpreter','none'); colorbar; colormap Jet;
yticks([1:length(rec_names)]); yticklabels(rec_names);

figure; stem(self_effects); title('Self effects')
xticklabel_rotate([1:length(rec_names)],45,rec_names,'interpreter','none'); colorbar; colormap Jet;