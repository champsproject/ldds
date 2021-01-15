clear all;

% Setup Parameters

% Initial Time
to = 0;  
% Time Step
dt = 0.05;
% Integration Time for LDs
tau = 5;

xo_grid = 0;%-50;
xf_grid = 2;%50;
yo_grid = 0;%50;
yf_grid = 1;%0;
grid_points = 400   ;

potencia = 0.75;
num_exp = 30;
nfig = 1;

% Formamos la matriz de condiciones iniciales

[ X, Y, x, y ] = form_grid( xo_grid, xf_grid, yo_grid, yf_grid, grid_points );

MSp_value = 0;

for omega=1:num_exp

    [Lagrangiano, a, b] = SLD( X, Y, dt, to, tau, potencia, omega );
    MSp_value = MSp_value + Lagrangiano;
    
end

MSp_value = MSp_value/num_exp;


pintar( X, Y, MSp_value, x, y, a, b, tau, potencia, nfig );


