function p_new = HenonMap(p,model_params,DLD_params)
    p_new = p;
    
    x = p(:,1);
    y = p(:,2);
    N_it = p(:,3);
    
    A = model_params(1);
    B = model_params(2);
    
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
    
        p_new(idx,1) = A - x_reg.^2 + B .* y_reg;
        p_new(idx,2) = x_reg;
        p_new(idx,3) = N_it_reg + 1;
    else
        p_new(:,1) = A - x.^2 + B .* y;
        p_new(:,2) = x;
        p_new(:,3) = N_it + 1;
    end
end


