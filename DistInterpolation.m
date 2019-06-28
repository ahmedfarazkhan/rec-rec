function Y_interp = DistInterpolation(Y_data, conn, ind_targets, Y_mu, Y_sigma)
% DistInterpolation
% Interpolate based on distances between data points
% Inputs
%   - Y_data: for n_pts but only non ind_target indices filled
%   - conn: [n_pts, n_pts] for all data points
%   - ind_targets: indices of pts in dists matrix to interpolate

n_pts = size(conn, 2); 
n_dims = size(Y_data, 2);

Y_interp = zeros(n_pts, n_dims);
Y_data_norm = (Y_data .* Y_sigma) + Y_mu;

for i=1:n_pts
    % Data points, do not interpolate
    if ~ismember(i, ind_targets)
        Y_interp(i, :) = Y_data(i, :);
      
    else    
        % Weighted sum by distance (rows sum to one)
        % [n_dims] = [n_dims, n_pts] x [n_pts, 1]
        interp_val = Y_data_norm' * conn(:, i);
        
        % Connectivity sums to 1
        if sum(conn(:,i))
            interp_val = interp_val / sum(conn(:, 1));
            Y_interp(i, :) = ( interp_val - Y_mu) / (Y_sigma);  
            
        end
        
    end
end

