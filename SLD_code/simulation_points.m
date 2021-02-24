clear all;

% Simulation of experiment

num_simulations = 1000;
num_brow = 2;

% Spatial parameters

xo = 0;
yo = 0;
xy_forw(:,1) = xo*ones(num_simulations,1);
xy_forw(:,2) = yo*ones(num_simulations,1);
xy_back(:,1) = xo*ones(num_simulations,1);
xy_back(:,2) = yo*ones(num_simulations,1);
b_forw = zeros(num_simulations,2);
b_back = zeros(num_simulations,2);

%x_point = xo*ones(num_simulations);
%y_point = yo*ones(num_simulations);

% Temporal parameters

to = 0;
tf = 50;
htemp = 0.075;

% Axis parameters

xo_axis = -1.2;
xf_axis = 1.2;
yo_axis = -0.1;
yf_axis = 0.4;

Delta_brow_forw = zeros(1,num_brow);
Delta_brow_back = zeros(1,num_brow);

figure(1);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

for  t = to:htemp:tf
    
    for j = 1:num_simulations
        
        for i=1:num_brow
            Delta_brow_forw(1,i) = normrnd(0,sqrt(htemp));
            Delta_brow_back(1,i) = normrnd(0,sqrt(htemp));
        end
    
        b_forw(j,:) = rk4( t, xy_forw(j,1), xy_forw(j,2), htemp, Delta_brow_forw);        
        b_back(j,:) = rk4( -t, xy_back(j,1), xy_back(j,2), -htemp, Delta_brow_back);        
                        
    end
    
    xy_forw = b_forw;
    xy_back = b_back;
    
    % Plot
    
    plot( xy_forw(:,1), xy_forw(:,2),'m*','LineWidth',1);
    hold on;
    plot( xy_back(:,1), xy_back(:,2),'g*','LineWidth',1);    
    title( ['Magenta = forw iter; Green = backw iter, time = ', num2str(t) ] );
    xlabel('x','FontSize',16);
    ylabel('y','FontSize',16);
    axis([xo_axis, xf_axis, yo_axis, yf_axis]);
    pause(0.5);
    %shading flat;
    hold off;
    
end

