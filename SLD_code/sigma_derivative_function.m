

function [ sigma_der_x, sigma_der_y ] = sigma_derivative_function( t, X, Y, num_brow )

    % fsigma_der = (0s 0s);
    
    m = length(X);
    sigma_der_x = zeros(m,2,num_brow);
    sigma_der_y = zeros(m,2,num_brow);
    epsilon = 0;
    
    sigma_der_x(:,1,1) = zeros(m,1);
    sigma_der_x(:,2,1) = epsilon*ones(m,1);
    sigma_der_x(:,1,2) = zeros(m,1);
    sigma_der_x(:,2,2) = zeros(m,1);
    
    sigma_der_y(:,1,1) = zeros(m,1);
    sigma_der_y(:,2,1) = zeros(m,1);
    sigma_der_y(:,1,2) = zeros(m,1);
    sigma_der_y(:,2,2) = zeros(m,1);
        
    %sigma_der(:,1,1) = zeros(m,1);
    %sigma_der(:,2,1) = zeros(m,1);
    
    %sigma_der(:,1,2) = zeros(m,1);
    %sigma_der(:,2,2) = zeros(m,1);

end

