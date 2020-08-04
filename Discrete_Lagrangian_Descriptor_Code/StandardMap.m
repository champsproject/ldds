function p_new = StandardMap(p,model_params,DLD_params)

    p_new = p;
    
    x = p(:,1);
    y = p(:,2);
    N_it = p(:,3);

    K = model_params(1);
    
    flag_vt = DLD_params(4);
    b_x1 = DLD_params(5);
    b_x2 = DLD_params(6);
    b_y1 = DLD_params(7);
    b_y2 = DLD_params(8);
    
    if flag_vt
        % Escape Condition for Variable-Iterations
        idx = (x >= b_x1) & (x <= b_x2) & (y >= b_y1) & (y <= b_y2); 
    
        x_reg = x(idx);
        y_reg = y(idx);
        N_it_reg = N_it(idx);
    
        p_new(idx,1) = x_reg + y_reg - (K/(2*pi)) .* sin(2 .* pi .* x_reg);
        p_new(idx,2) = y_reg - (K/(2*pi)) .* sin(2 .* pi .* x_reg);
        p_new(idx,3) = N_it_reg + 1;
    else
        p_new(:,1) = x + y - (K/(2*pi)) .* sin(2 .* pi .* x);
        p_new(:,2) = y - (K/(2*pi)) .* sin(2.*pi .* x);
        p_new(:,3) = N_it + 1;
    end
    
end

