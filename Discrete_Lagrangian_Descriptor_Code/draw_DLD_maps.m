function draw_DLD_maps(mesh_params,DLD_params,flag_type,flag_av,ld_fw,ld_bw)
    
    xi = mesh_params(1);
    xf = mesh_params(2);
    nx = mesh_params(3);
    yi = mesh_params(4);
    yf = mesh_params(5);
    ny = mesh_params(6);
    
    xp = linspace(xi,xf,nx);
    yp = linspace(yi,yf,ny);
    
    N = DLD_params(1);
    flag_m = DLD_params(2);
    pval = DLD_params(3);
    
    if flag_m
        str_meth = ['p-norm (p = ',num2str(pval),') -- '];
    else
        str_meth = '';
    end
        
     if flag_type == 1
        LD = reshape(ld_fw,ny,nx);
        if flag_av
            LD = LD ./ N;
        end
        string_title = ['Forward LDs ',str_meth,'(','$N = $ ',num2str(N),')'];
    elseif flag_type == 2
        LD = reshape(ld_bw,ny,nx);
        if flag_av
            LD = LD ./ N;
        end
        string_title = ['Backward LDs ',str_meth,'(','$N = $ ',num2str(N),')'];
    else
        LD = reshape(ld_fw + ld_bw,ny,nx);
        if flag_av
            LD = LD ./ (2.*N);
        end
        string_title = ['LDs ',str_meth,'(','$N = $ ',num2str(N),')'];
     end
     
     LD = LD ./ max(max(LD));  % Scale LD output
     
     
    pcolor(xp,yp,LD);
    
    shading interp
    colormap bone
      
    set(gca,'FontSize',20)
    
    title(string_title,'FontSize',18,'Interpreter','latex');
   
    axis square
    
    xlabel('$x$','FontSize',20,'Interpreter','latex')
    ylabel('$y$','Interpreter','latex','FontSize',20,'Rotation',0);
    
    colorbar
    
    
end

