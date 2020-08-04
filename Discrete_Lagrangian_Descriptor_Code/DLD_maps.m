function [ld_fw,ld_bw,N_it_fw,N_it_bw] = DLD_maps(ds,ds_inv,mod_params,mesh_params,DLD_params)

    xi = mesh_params(1);
    xf = mesh_params(2);
    nx = mesh_params(3);
    yi = mesh_params(4);
    yf = mesh_params(5);
    ny = mesh_params(6);

    [X,Y] = meshgrid(linspace(xi,xf,nx),linspace(yi,yf,ny));
   
    grid_size = nx * ny;
    mesh = [X(:) Y(:)];

    N = DLD_params(1);
    flag_m = DLD_params(2);
    pval = DLD_params(3);
    
    % Forward Iterations
    N_fw = zeros(grid_size,1);
    ld_fw = zeros(grid_size,1);
    
    p_fw = [mesh N_fw];           
    f = str2func(ds);
      
    for k = 1 : N
        p_fw_2 = f(p_fw,mod_params,DLD_params);
        if flag_m == 1
            if pval <= 1
                ld_fw = ld_fw + (abs(p_fw_2(:,1) - p_fw(:,1)).^pval + abs(p_fw_2(:,2) - p_fw(:,2)).^pval);
            end
        end
        p_fw = p_fw_2;
    end
    
    N_it_fw = p_fw(:,3);
    
    
    % Backward Iterations
    N_bw = zeros(grid_size,1);
    ld_bw = zeros(grid_size,1);
    
    f_inv = str2func(ds_inv);
    p_bw = [mesh N_bw];
    
    for k = 1 : N
        p_bw_2 = f_inv(p_bw,mod_params,DLD_params);
        
        if flag_m == 1
            if pval <= 1
                ld_bw = ld_bw + (abs(p_bw_2(:,1) - p_bw(:,1)).^pval + abs(p_bw_2(:,2) - p_bw(:,2)).^pval);
            end
        end
        p_bw = p_bw_2;
    end
    
    N_it_bw = p_bw(:,3);
    
end
