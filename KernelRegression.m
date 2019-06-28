function Y_pred = KernelRegression(X_train, Y_train, X_pred, sigma)
% Kernel regression (Gaussian kernel & Euclidean distance by default)
% Inputs:
%   X_train: [M, K] array of training features
%   Y_train: [M, D] array of training labels
%   X_pred: [N, K] array of test features
%   vars: [K]

    mu = zeros(1, size(X_train,2));
    Y_pred = zeros(size(X_pred,1), size(Y_train,2));
    
    for i=1:size(X_pred, 1)
        
        kernel_norm = 0.0;
        
        for j=1:size(X_train,1)
            
            %kernel_sample = mvnpdf(sqrt(sum((X_pred(i,:) - X_train(j,:)) .^2) ), mu, sigma);
            kernel_sample = mvnpdf(X_pred(i,:), X_train(j,:), sigma);            
            Y_pred(i,:) = Y_pred(i,:) + (kernel_sample * Y_train(j, :)); 
            kernel_norm = kernel_norm + kernel_sample;
                   
        end
        
        if kernel_norm
            Y_pred(i, :) = Y_pred(i,:) / kernel_norm;
        end
        
    end

end

