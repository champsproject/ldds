"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp as integrator # RK45 (default)

def generate_points(GRID_PARAMETERS):
    """
    Returns a 1D array of all points from a meshgrid with dimensions and size defined by list of input parameters.
    
    Parameters
    ----------
    GRID_PARAMETERS : list of 3-tuples of floats
        input parameters of limits and size of mesh per axis
        
    Returns
    -------
    mesh : 1d numpy array
        points from meshgrid 
    """
    x_min, x_max, Nx = GRID_PARAMETERS[0]
    y_min, y_max, Ny = GRID_PARAMETERS[1]
    points_x = np.linspace(x_min, x_max, Nx)
    points_y = np.linspace(y_min, y_max, Ny)    
    Y, X = np.meshgrid(points_y, points_x)  # Grid in phase-space
    mesh = np.transpose([X.flatten(), Y.flatten()]) # 2D grid
    return mesh

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


def check_if_points_escape_box(u, box_BOUNDARIES):
    """
    Determine if points in phase-space u have scaped box with user-defined defined dimensions
    
    Parameters
    ----------
    u : array_like, shape(n, )
        points in phase-space to check if outside box boundaries
    
    box_BOUNDARIES : list of 2-tuples of floats
        box lower and upper limits along X and Y axes
        
    Returns
    -------
    u_indices : array_like, shape(n, )
        array of True/False bool values if points inside/outside the box
    """
    x, y = u.T
    # Escape condition
    box_x_min, box_x_max = box_BOUNDARIES[0]
    box_y_min, box_y_max = box_BOUNDARIES[1]
    u_indices = (x >= box_x_min) & (x <= box_x_max) & (y >= box_y_min) & (y <= box_y_max)
    
    return u_indices


def vector_field_flat(t, points, vector_field, box_BOUNDARIES):
    """
    Returns input vector field values (y0) for integrator as 1d array 
    obtained from applying vector_field to a set of initial conditions (points)
    
    Parameters
    ----------
    t : float
        time
    
    points : ndarray, shape(n,2)
    
    vector_field: function
    
    box_BOUNDARIES : list of 2-tuples, optional
        box boundaries for escape condition of variable time integration
        boundaries are infinite by default.
        
    Returns
    ------- 
    1d array
        y0 values for integrator 
    """
    u = points.reshape((-1,2))
    # Apply Escape condition
    u_inbox = check_if_points_escape_box(u, box_BOUNDARIES)
    # Define output vector field in combination with escape condition
    v = np.zeros(u.shape)
    v[u_inbox == True] = vector_field(t, u[u_inbox == True])
    
    return v.flatten()


def accumulate_lagrangian_descriptor(LD, points_initial, points_final, LD_PARAMETERS):
    """
    Returns the cumulative values of the LD function from trajectory segments in the evolution of a system.
    The trajectory segments run from points_initial to points_final.
    
    Parameters
    ----------
    LD : array_like, shape(n, )
        previous cummulative values of the LD function.
    
    points_initial : array_like, shape(n, )
        initial points in phase space in the evolution of trajectories.
        
    points_final : array_like, shape(n, )
        points ahead of initial points in the evolution of their dynamics.
        
    LD_PARAMETERS : list made of a 3-tuple and a float
        input parameters for LD computation
        3-tuple contains floats t_initial, t_final, dt (timestep)
        float is p-value of Lp-norm.

    Returns
    -------
    LD_accum : array_like, shape(n, )
        new cummulative values of the LD function from trajectory segment
    """
    t_initial, t_final, dt = LD_PARAMETERS[0]
    p_norm = LD_PARAMETERS[1]
    
    F = dt**(1 - p_norm) # Scaling factor from discrete integration
    points_difference = points_final - points_initial
    
    if p_norm <= 1: # p-norm LD
        LD_accum = LD + F*np.sum( np.abs(points_difference)**p_norm, axis=1)
    if p_norm == 2: # Arclength LD
        LD_accum = LD + np.linalg.norm(points_difference, axis=1)
    elif 1 < p_norm < 2: # Discrete arclength LD
        LD_accum = LD + np.sum( (points_difference**p_norm)**(1/p_norm), axis=1)
    
    return LD_accum


def extract_y_final(solution_object):
    """
    Returns the y-values from only at time `t + dt` from solution_object outputted by integrator.
    
    Parameters
    ----------
    solution_object: scipy.integrate.solve_ivp output object.
    
    Returns
    -------
    y_final: ndarray, shape (n, 2)
        extracted y-values
    """
    y_final = solution_object.y.T[-1]
    
    return y_final


def compute_lagrangian_descriptor(points_initial, vector_field, LD_PARAMETERS, box_BOUNDARIES = [(-np.infty, np.infty), (-np.infty, np.infty)]):
    """
    Returns the values of the LD function from integrated trajectories from initial conditions in phase-space.
    
    Parameters
    ----------
    points_initial : ndarray, shape(n, )
        initial conditions in phase-space.
    
    vector_field: function
        vector field over phase-space
        
    LD_PARAMETERS : list made of a 3-tuple and a float
        input parameters for LD computation
        3-tuple contains floats t_initial, t_final, dt (timestep)
        float is p-value of Lp-norm.
    
    box_BOUNDARIES : list of 2-tuples, optional
        box boundaries for escape condition of variable time integration
        boundaries are infinite by default.
    
    Returns
    -------
    LD : 1d array, shape (n, )
        array of computed LD values for all initial conditions.
    """
    f = vector_field_flat
    
    t_initial, t_final, dt = LD_PARAMETERS[0]
    p_norm = LD_PARAMETERS[1]
    
    if (t_final - t_initial)*dt < 0:
        dt = -dt
        
    shape = points_initial.shape
    LD = np.zeros(shape[0])  # Array to store forward LDs
    for t in np.arange(t_initial, t_final - dt, dt):
        ###################################
        # Inputs for integrator
        t_span = (t, t + dt)
        y0 = points_initial.flatten()
        ###################################
        # Outputs from integration
        solution_object = integrator(f, t_span, y0, args=(vector_field, box_BOUNDARIES))
        points_final = extract_y_final(solution_object)
        points_final = points_final.reshape(shape)
        ###################################
        # Compute LD
        LD = accumulate_lagrangian_descriptor(LD, points_initial, points_final, LD_PARAMETERS)
        ###################################
        # For next iteration
        points_initial = points_final
    
    return LD


__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido'
__status__ = 'Development'
