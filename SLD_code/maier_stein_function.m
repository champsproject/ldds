

function [ XY ] = maier_stein_function( t, X, Y )

    m = length(X);
    XY = zeros(m,2);
    beta = 10;
    
    XY(:,1) = X - X.^3 - beta*X.*Y.^2;
    XY(:,2) = -(1 + X.^2).*Y;

end

