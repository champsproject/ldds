%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               TWO-DIMENSIONAL MAP DEFINITION                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ds = 'HenonMap';               % Name of Forward Iteration Function
ds_inv = [ds,'_inv'];    % Name of Backward Iteration Function

A = 0.298;
B = 1;

model_params = [A B];

%%%%%%%%%%%%%%%%%%%
%       DLD Method Setup Parameters     %
%%%%%%%%%%%%%%%%%%%
N = 12;  % Numer of Iterations of the Map (Forward and Backward)
flag_m = 1;
p_val = 1/2;

% Variable Iteration DLD Setup
flag_vt = 1;
bound_x1 = -5;
bound_x2 = 5;
bound_y1 = -5;
bound_y2 = 5;

DLD_params = [N flag_m p_val flag_vt bound_x1 bound_x2 bound_y1 bound_y2];

% Mesh to compute LDs
xi = -1.5;
xf = 1.5;
nx = 1200;
yi = -1.5;
yf = 1.5;
ny = 1200;

mesh_params = [xi xf nx yi yf ny];

% Computation of LDs
[ld_fw,ld_bw,~,~] = DLD_maps(ds,ds_inv,model_params,mesh_params,DLD_params);

% Draw Forward LDs
flag_type = 1;
flag_av = 0;
figure
draw_DLD_maps(mesh_params,DLD_params,flag_type,flag_av,ld_fw,ld_bw);

% Draw Backward LDs
flag_type = 2;
flag_av = 0;
figure
draw_DLD_maps(mesh_params,DLD_params,flag_type,flag_av,ld_fw,ld_bw);

% Draw LDs
flag_type = 3;
flag_av = 0;
figure
draw_DLD_maps(mesh_params,DLD_params,flag_type,flag_av,ld_fw,ld_bw);