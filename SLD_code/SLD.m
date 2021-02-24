

function [ C, a, b ] = SLD( X, Y, htemp, to, tf, potencia, omega )

    a = 0;
    b = 0;
    
    [m,n] = size(X);    
    
    tamano_grid = m*n;
    
    A(:,1) = reshape(X, tamano_grid, 1);
    A(:,2) = reshape(Y, tamano_grid, 1);    

    A_forw = A; % Esta copia es para hacer RK4 hacia delante.
    A_back = A; % Esta copia es para hacer RK4 hacia atr'as.
    
    %B = zeros(tamano_grid,2);
    %C = zeros(n,n); % Matriz a pintar con los valores de la integral.
    
    Integral = zeros(tamano_grid,1);
    
    % Create brownian
    
    num_brow = 2;
    long_delta = round(4*abs(tf-to)/htemp + 1);
    Delta_brow_forw = zeros(long_delta,num_brow);
    Delta_brow_back = zeros(long_delta,num_brow);
    
    total_time = 0;
    
    for i=1:long_delta
        for j=1:num_brow
            Delta_brow_forw(i,j) = normrnd(0, sqrt(htemp/2));
            if i==1
                Delta_brow_back(i,j) = Delta_brow_forw(i,j);
            else
                Delta_brow_back(i,j) = normrnd(0, sqrt(htemp/2));
            end
        end
    end
    
    % Calculamos las trayectorias con RK4 y luego integramos hacia delante y hacia atr'as:
    
    for  t = to:htemp:tf
            
        tic        
        
        B_forw = rk4( t, A_forw(:,1), A_forw(:,2), htemp, Delta_brow_forw);
        B_back = rk4( -t, A_back(:,1), A_back(:,2), -htemp, Delta_brow_back);
        
        if potencia<=1
        
            Integral = Integral + (abs(B_forw(:,1) - A_forw(:,1))).^potencia + (abs(B_forw(:,2) - A_forw(:,2))).^potencia;% .* abs(htemp) .^ (1-potencia);
            Integral = Integral + (abs(B_back(:,1) - A_back(:,1))).^potencia + (abs(B_back(:,2) - A_back(:,2))).^potencia;% .* abs(htemp) .^ (1-potencia);
        
        else
            
            Integral = Integral + sqrt((B_forw(:,1) - A_forw(:,1)).^2 + (B_forw(:,2) - A_forw(:,2)).^2);
            Integral = Integral + sqrt((B_back(:,1) - A_back(:,1)).^2 + (B_back(:,2) - A_back(:,2)).^2);
        
        end
        
        pos_brow = round(2*t/htemp + 1);
        
        a = -exp(-t)*Delta_brow_forw(pos_brow,1) + a;
        b = -exp(-t-htemp)*Delta_brow_back(pos_brow,2) + b;
        
        time = toc;
        disp(['Omega ',num2str(omega),' time ',num2str(t),' finished in ',num2str(time),'seconds']);
        
        A_forw = B_forw;
        A_back = B_back;
            
        total_time = total_time + toc;
    
    end
        
    disp(['Elapsed time ',num2str(total_time),'seconds']);
    
    % Recolocamos la matriz integral como una matriz n*n para poder pintarla

    C = reshape(Integral,m,n);

end

