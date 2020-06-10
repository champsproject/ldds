"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def generate_points(grid_parameters):
    """
    Returns a 1D array of all points from a on a uniform grid with dimensions and size defined by list of input parameters.
    An additional dimension initiallised with zeros is added for the calculation of Lagrangian Descriptors.
    
    Parameters
    ----------
    grid_parameters : list of 3-tuples of floats
        input parameters of limits and size of mesh per axis
        
    Returns
    -------
    mesh : 1d numpy array
        flattened array of initial conditions
    """
    x_min, x_max, Nx = grid_parameters[0]
    y_min, y_max, Ny = grid_parameters[1]
    points_x = np.linspace(x_min, x_max, Nx)
    points_y = np.linspace(y_min, y_max, Ny)    
    Y, X = np.meshgrid(points_y, points_x)  # Grid in phase-space
    mesh = np.transpose([X.flatten(), Y.flatten(), np.zeros(Nx*Ny)]) # 2D grid + a zero column for LDs
    return mesh.flatten()

def perturb_field(vector_field, perturbation):
    """
    Returns the vector field function with a linearly added pertubation
    Both input function should input (t, u), with t: float, and u: ndarray
    Also, the output of these funcs must be ndarrays of the same shape
    
    Parameters
    ----------
        vector_field: function
            unperturbed vector field
    
        perturbation: function
            forcing added to the vector field
    
    Returns
    -------
        perturbed function
    """
    return lambda t, u: vector_field(t, u) + perturbation(t, u)


def check_if_points_escape_box(u, box_boundaries):
    """
    Determine if points in phase-space u have scaped box with user-defined defined dimensions
    
    Parameters
    ----------
    u : array_like, shape(n, )
        points in phase-space to check if outside box boundaries
    
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

def lag_des(u, v, p_value = 0.5):
    """
    Vector field equation for Lagrangian descriptor.

    Parameters
    ----------
    v : ndarray, shape(n,2)
        Vector field at given point.
            
    p_value : float
        Exponent in Lagrangian descriptor definition.

    Returns
    -------
    LD : ndarray, shape(n,1)
        Vector field for Lagrangian descriptor dimension.
    """
    if p_value == 0:
        LD = np.abs(u[:,1]*v[:,0])
    elif p_value>0:
        LD = np.sum(np.abs(v)**p_value, axis=1)
    else:
        LD=np.zeros(len(u[:,0]))
    return LD

def vector_field_flat(t, points, vector_field, p_value, box_boundaries):
    """
    Returns vector field values for integration of flattened input array.
    
    Parameters
    ----------
    t : float
        time
    
    points : ndarray, shape(n,3)
    
    vector_field: function
        User defined vector field.
    
    p_value : float
        Exponent in Lagrangian descriptor definition.
    
    box_boundaries : list of 2-tuples, optional
        box boundaries for escape condition of variable time integration
        boundaries are infinite by default.
        
    Returns
    ------- 
    1d array
        y0 values for integrator 
    """
    u = points.reshape((-1,3))
    u = u[:,:-1] #remove LD values
    # Apply Escape condition
    u_inbox = check_if_points_escape_box(u, box_boundaries)
    # Define output vector field in combination with escape condition
    v = np.zeros(u.shape)
    v[u_inbox == True] = vector_field(t, u[u_inbox == True])
    # Calculate LD vector field
    LD_vec = np.zeros(len(u))
    LD_vec [u_inbox == True] = lag_des(u[u_inbox == True], v[u_inbox == True], p_value)
    # Add LD
    v_out=np.column_stack((v, LD_vec))
    return v_out.flatten()

def compute_lagrangian_descriptor(grid_parameters, vector_field, tau, p_value=0.5, box_boundaries = [(-np.infty, np.infty), (-np.infty, np.infty)]):
    """
    Returns the values of the LD function from integrated trajectories from initial conditions in phase-space.
    
    Parameters
    ----------
    grid_parameters : list of 3-tuples of floats
        input parameters of limits and size of mesh per axis
    
    vector_field: function
        vector field over phase-space
        
    tau : float
        Upper limit of integration.
        
    p_value : float, optional
        Exponent in Lagrangian descriptor definition. The default is 0.5.
    
    box_boundaries : list of 2-tuples, optional
        Box boundaries for escape condition of variable time integration.
        Boundaries are infinite by default.
    
    Returns
    -------
    LD : ndarray, shape (Nx, Ny)
        Array of computed Lagrangian descriptor values for all initial conditions.
    """
    
    f = lambda x, y: vector_field_flat(x, y, vector_field, p_value, box_boundaries)
    y0 = generate_points(grid_parameters)
    
    solution = solve_ivp(f, [0,tau], y0, t_eval=[tau], rtol=1.0e-4)
    
    n=len(grid_parameters)
    LD_values = solution.y[n::n+1] #values corresponding to LD
    
    LD=np.abs(LD_values).reshape((grid_parameters[0][2],grid_parameters[1][2])).T #reshape to 2D array
    if p_value<=1:
        return LD
    else:
        return LD**(1/p_value)

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
