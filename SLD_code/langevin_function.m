

function [ XY ] = langevin_function( t, X, Y )

    m = length(X);
    XY = zeros(m,2);
    
    V0 = 0.1;
    a = 100;
    gamma = 0.001;
    V_prima_X = (V0/(a*sqrt(pi))) * exp(-(X.^2)/(2*a*a));

    XY(:,1) = Y;
    XY(:,2) = - V_prima_X - gamma*Y;


end
