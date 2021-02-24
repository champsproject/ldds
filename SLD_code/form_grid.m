

function [ X, Y, x, y ] = form_grid( xo_grid, xf_grid, yo_grid, yf_grid, grid_points )

    xgrid = (xf_grid-xo_grid)/grid_points;
    ygrid = (yf_grid-yo_grid)/grid_points;

    x = xo_grid:xgrid:xf_grid;
    y = yo_grid:ygrid:yf_grid;
    
    [X, Y] = meshgrid(x,y);    

end

