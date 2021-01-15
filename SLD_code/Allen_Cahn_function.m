

function [ XY ] = Allen_Cahn_function( t, X, Y )

    m = length(X);
    XY = zeros(m,2);
    alpha = 1.12;
    
    XY(:,1) = 1/alpha*(X - Y - X.^3 + Y.^3) - X;
    XY(:,2) = 1/alpha*(Y - X - Y.^3 + X.^3) - Y;

end
