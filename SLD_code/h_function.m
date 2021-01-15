

function [ XY ] = h_function( t, X, Y, h_temp, Delta_brow )

    % dX_t = X_0 + b(t,X,Y) dt + sigma(t,X,Y) dW_t
    m = length(X);
    num_brow = 2;
    suma_sigma_brow = zeros(m,2);
    suma_sigma_sigma_der = zeros(m,2);
        
    b = dgyre_function(t, X, Y);
    
    sigma = sigma_function(t, X, Y, num_brow);
      
    [sigma_der_x, sigma_der_y] = sigma_derivative_function(t, X, Y, num_brow);
    
    for k=1:num_brow
        suma_sigma_brow = suma_sigma_brow + sigma(:,:,k) * Delta_brow(1,k);
        suma_sigma_sigma_der(:,1) = suma_sigma_sigma_der(:,1) + sigma(:,1,k).*sigma_der_x(:,1,k) + sigma(:,2,k).*sigma_der_y(:,1,k);
        suma_sigma_sigma_der(:,2) = suma_sigma_sigma_der(:,2) + sigma(:,1,k).*sigma_der_x(:,2,k) + sigma(:,2,k).*sigma_der_y(:,2,k);
    end
    
    XY = ( b - 0.5*suma_sigma_sigma_der ) * h_temp + suma_sigma_brow;

end

