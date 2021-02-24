

% Perform the modulus of the gradient of a matrix

function [ edge_x, edge_y, gradient_mod ] = gradient_modulus( map, h_file, h_col )

    %[ydim, xdim] = size(map);

    [m,n] = size(map);
    
    edge_x = zeros(m,n);
    edge_y = zeros(m,n); 
    
    edge_x(2:m-1,2:n-1) = ( map(2:m-1,3:n) - map(2:m-1,1:n-2) ) / (2*h_col);
    edge_y(2:m-1,2:n-1) = ( map(3:m,2:n-1) - map(1:m-2,2:n-1) ) / (2*h_file);
    
    gradient_mod = sqrt(edge_x.*edge_x + edge_y.*edge_y);

end

