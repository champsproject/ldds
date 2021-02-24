

function [ B ] = euler_maruyama( t, X, Y, htemp, Delta_brow )

    pos_brow = round(2*t/htemp + 1);
    num_brow = 2;
    m = length(X);
    suma_sigma_brow = zeros(m,2);
    XY(:,1) = X;
    XY(:,2) = Y;
    
    b = saddle_function( t, X, Y );
    sigma = sigma_function( t, X, Y, num_brow );
    
    for k=1:num_brow        
        suma_sigma_brow = suma_sigma_brow + sigma(:,:,k) * Delta_brow(pos_brow,k);        
    end
    
    B = XY + b*htemp + suma_sigma_brow;  

end

