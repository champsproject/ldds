import numpy as np
from operator import itemgetter
from pylds.base import generate_points, lagrangian_descriptor

def check_if_points_escape_box(u, box_boundaries):
    """
    Determine if points, u, in 2D plane have escaped box with user-defined dimensions.
    
    Parameters
    ----------
    u : array_like, shape(n, 2)
        Points in plane to check if outside box boundaries
    
    box_boundaries : list of 2-tuples of floats
        Box lower and upper limits along X and Y axes
        
    Returns
    -------
    u_indices : array_like, shape(n, 2)
        array of True/False bool values if points inside/outside the box
    """
    x, y = u.T
    # Escape condition
    box_x_min, box_x_max = box_boundaries[0]
    box_y_min, box_y_max = box_boundaries[1]
    u_indices = (x >= box_x_min) & (x <= box_x_max) & (y >= box_y_min) & (y <= box_y_max)
    return u_indices

def pbc_correction_coords_single_axis(x, box_origin, box_length):
    """
    Correct single-axis coordinates by Periodic Boundary Conditions (PBCs).

    Parameters
    ----------
    x : array
        Points in axis to correct by PBCs
    
    box_origin : float or False(bool)
        If False, no PBC-correction applied
        
    box_length : float or False(bool)
        If False, no PBC-correction applied
        
    Returns
    -------
    x_pbc : array
        1D array of PBC-corrected coordinates
    """
    x_pbc = x
    x0 = box_origin
    L = box_length
    if not isinstance(x0, bool) or not isinstance(L, bool):
        #apply PBC correction
        x = x + L/2 - x0
        x = np.mod(x + 2*L, L) 
        x_pbc = x - L/2 + x0
        
    return x_pbc

def pbc_correction_coords(u, periodic_boundaries):
    """
    Correct coordinates in 2D plane by Periodic Boundary Conditions (PBCs).
    The periodic box can be a rectangle (torus) or an infinite band (cylinder) in the plane.
    
    Parameters
    ----------
    u : array_like, shape(n, 2)
        Points in 2D plane to correct by PBCs.
    
    periodic_boundaries : list of 2-tuples of floats or False (bools)
        The first tuple must be the box origin coordinates in the plane.
        The second tuple must be the lengths of the box edges.
        For a torus PBC, both tuples entries must be floats.
        For a cylinder PBC, tuple entries for the non-periodic axis must be False.
        If no PBCs, this must be a False (single bool).
        
    Returns
    -------
    u_pbc : array_like, shape(n, 2)
        array of PBC-corrected points for periodic box.
    """
    N_dims = len(periodic_boundaries)
    pbc = periodic_boundaries
    f = pbc_correction_coords_single_axis
    u_pbc = [f(u.T[i], *pbc[i]) for i in range(N_dims)]
    
    return np.column_stack(u_pbc)

def pbc_correction_distances_single_axis(dx, box_length):
    """
    Correct distances in a single axis by Periodic Boundary Conditions (PBCs).

    Parameters
    ----------
    dx : array_like, shape(n, )
        Distances in axis to correct by PBCs.
        
    box_length : float or bool
        If no PBCs, this must be False (bool).
        
    Returns
    -------
    dx_pbc : array_like, shape(n, )
        array of corrected distances for periodic box.
    """
    dx_pbc = dx
    L = box_length
    if not isinstance(L, bool):
        L = abs(L)
        nint = lambda x: np.round(x).astype(int) #nearest integer
        dx_pbc = dx - nint(dx/L) #minimum image criterion
    
    return dx_pbc
    
def pbc_correction_distances(du, periodic_boundaries):
    """
    Correct distances in 2D plane by Periodic Boundary Conditions (PBCs).
    The periodic box can be a rectangle (torus) or an infinite band (cylinder) in the plane.

    Parameters
    ----------
    du : array_like, shape(n, 2)
        Distances in plane to correct by PBCs.
    
    periodic_boundaries : list of 2-tuples of floats or False (bools)
        The first tuple must be the box origin coordinates in the plane.
        The second tuple must be the lengths of the box edges.
        For a torus PBC, both tuples entries must be floats.
        For a cylinder PBC, tuple entries for the non-periodic axis must be False.
        If no PBCs, this must be a False (single bool).
        
    Returns
    -------
    du_pbc : array_like, shape(n, 2)
        array of PB-corrected distances for periodic box.
    """
    N_dims = len(periodic_boundaries)
    pbc = periodic_boundaries
    f = pbc_correction_distances_single_axis
    du_pbc = [f(du.T[i], pbc[i][-1]) for i in range(N_dims)]
    
    return np.column_stack(du_pbc)

def compute_lagrangian_descriptor(grid_parameters, discrete_map, N_iterations, p_value=0.5, box_boundaries=False, periodic_boundaries=False):
    """
    Returns the values of the LD function from trajectories from iterated initial conditions in plane by a map.
    
    Parameters
    ----------
    grid_parameters : list of 3-tuples of floats
        Input parameters of limits and size of mesh per axis.
    
    discrete_map: function
        Map of discrete 2D dynamical system.
        
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
        
    perodic_boundaries: list of floats, optional
        Lenght values of periodic box axes (2D default).
        PBC are False by default.
    
    Returns
    -------
    LD : ndarray, shape (Nx, Ny)
        array of computed Lagrangian descriptor values for all initial conditions.
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
            dy = pbc_correction_distances(dy, periodic_boundaries)
            y0 = pbc_correction_coords(y0, periodic_boundaries)
            y  = pbc_correction_coords(y , periodic_boundaries)
                
        LD_values = LD_values + lagrangian_descriptor(y0, dy, p_value)
        y0 = y

    N_points_slice_axes = [x[-1] for x in grid_parameters] #take number of points
    LD = LD_values.reshape(*N_points_slice_axes) #reshape to 2-D array  

    if p_value<=1:
        return LD
    else:
        return LD**(1/p_value)
