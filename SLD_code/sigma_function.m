

function [ sigma ] = sigma_function( t, X, Y, num_brow )

    % fsigma = (1s 1s);
    
    m = length(X);
    sigma = zeros(m,2,num_brow);
    epsilon = 0.05;
        
    sigma(:,1,1) = epsilon*ones(m,1);
    sigma(:,2,1) = zeros(m,1);
    
    sigma(:,1,2) = zeros(m,1);
    sigma(:,2,2) = epsilon*ones(m,1);

end

