import numpy as np
from operator import itemgetter
from pylds.base import generate_points, lagrangian_descriptor

def check_if_points_escape_box(u, box_boundaries):
    x, y = u.T
    # Escape condition
    box_x_min, box_x_max = box_boundaries[0]
    box_y_min, box_y_max = box_boundaries[1]
    u_indices = (x >= box_x_min) & (x <= box_x_max) & (y >= box_y_min) & (y <= box_y_max)
    return u_indices

def compute_lagrangian_descriptor(grid_parameters, discrete_map, N_iterations, p_value=0.5, box_boundaries=False):
    y0, mask = generate_points(grid_parameters)
    y0 = y0.reshape(-1,3)
    y0 = y0[:,:-1]
    
    if box_boundaries == False:
        box_boundaries = [(-np.infty, np.infty), (-np.infty, np.infty)]
    
    f = discrete_map

    LD_values = np.zeros(len(y0))
    for i in range(N_iterations):
        y = f(y0)
        y_inbox = check_if_points_escape_box(y, box_boundaries)

        y[y_inbox == True]  = f(y0[y_inbox])
        y[y_inbox == False] = y0[y_inbox == False]

        LD_values = LD_values + lagrangian_descriptor(y0, y-y0, p_value)        
        y0 = y

    N_points_slice_axes = list( map(itemgetter(-1), grid_parameters))
    LD = LD_values.reshape(*N_points_slice_axes) #reshape to 2-D array  

    if p_value<=1:
        return LD
    else:
        return LD**(1/p_value)
