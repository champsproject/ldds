

function [ XY ] = dgyre_function( t, X, Y )

    m = length(X);
    XY = zeros(m,2);

    A = 0.25;%10
    s = 1;%50
    mu = 0;%0.005;
    epsilon = 0.25;%0.01;
    omega = 2*pi;%1/pi;
    psi = 0;
    
    f = epsilon * sin(omega*t + psi) .* X.^2 + (1 - 2*epsilon*sin(omega*t + psi)) .* X;
    partial_f_x = 2*epsilon * sin(omega*t + psi) .* X + 1 - 2*epsilon*sin(omega*t + psi);
    
    XY(:,1) = -pi*A * sin(pi/s .* f) .* cos(pi/s .* Y) - mu*X;
    XY(:,2) = pi*A * cos(pi/s .* f) .* sin(pi/s .* Y) .* partial_f_x - mu*Y;


end

