import numpy as np
from operator import itemgetter
from pylds.base import generate_points, lagrangian_descriptor

def check_if_points_escape_box(u, box_boundaries):
    """
    Determine if points in 2D plane, u, have escaped box with user-defined dimensions.
    
    Parameters
    ----------
    u : array_like, shape(n, )
        points in plane to check if outside box boundaries
    
    box_boundaries : list of 2-tuples of floats
        box lower and upper limits along X and Y axes
        
    Returns
    -------
    u_indices : array_like, shape(n, )
        array of True/False bool values if points inside/outside the box
    """
    x, y = u.T
    # Escape condition
    box_x_min, box_x_max = box_boundaries[0]
    box_y_min, box_y_max = box_boundaries[1]
    u_indices = (x >= box_x_min) & (x <= box_x_max) & (y >= box_y_min) & (y <= box_y_max)
    return u_indices

def compute_lagrangian_descriptor(grid_parameters, discrete_map, N_iterations, p_value=0.5, box_boundaries=False, periodic_boundaries=False):
    """
    Returns the values of the LD function from trajectories from iterated initial conditions in plane by a map.
    
    Parameters
    ----------
    grid_parameters : list of 3-tuples of floats
        input parameters of limits and size of mesh per axis.
    
    discrete_map: function
        map of discrete 2D dynamical system.
        
    tau : float
        Upper limit of integration.
        
    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.
    
    box_boundaries : list of 2-tuples, optional
        Box boundaries for escape condition of variable time integration.
        Boundaries are infinite by default.
    
    Returns
    -------
    LD : ndarray, shape (Nx, Ny)
        Array of computed Lagrangian descriptor values for all initial conditions.
    """
    N_mesh_axes = len(grid_parameters)+1
    y0, mask = generate_points(grid_parameters)
    y0 = y0.reshape(-1,N_mesh_axes)
    y0 = y0[:,:-1] # exclude LD-axis
        
    f = discrete_map

    LD_values = np.zeros(len(y0))
    for i in range(N_iterations):
        y = f(y0)
        # Escape box condition
        if box_boundaries:
            y_inbox = check_if_points_escape_box(y, box_boundaries)
            y[y_inbox == False] = y0[y_inbox == False]
        
        # Periodic Boundary conditions
        dy = y-y0
        if periodic_boundaries:
            nint = lambda x: np.round(x).astype(int) #nearest integer
            L = np.asarray(periodic_boundaries)
            if not L==0:
                dy = dy - nint(dy/L) #minimum image criterion
        
        LD_values = LD_values + lagrangian_descriptor(y0, dy, p_value)        
        y0 = y

    N_points_slice_axes = [x[-1] for x in grid_parameters] #take number of points
    LD = LD_values.reshape(*N_points_slice_axes) #reshape to 2-D array  

    if p_value<=1:
        return LD
    else:
        return LD**(1/p_value)
