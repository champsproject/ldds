"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""

import sys
import numpy as np
import numpy.ma as ma
from operator import itemgetter
from scipy.integrate import solve_ivp
from functools import reduce

def energy_conservation_condition(phase_space_axes, H0, potential_energy, momentum_sign):
    N_dim = len(phase_space_axes)
    phase_space_positions = np.transpose(phase_space_axes[:int(N_dim/2)])
    phase_space_momenta = np.transpose(phase_space_axes[int(N_dim/2):])
    V = potential_energy(phase_space_positions)
    points_dof_remaining = momentum_sign * \
        np.sqrt(2*(H0 - V) - (phase_space_momenta**2).sum(axis=1))
    return points_dof_remaining

def generate_points(grid_parameters):
    """
    Returns a 1D array of all points from a on a uniform grid with dimensions and size defined by list of input parameters.
    An additional dimension initiallised with zeros is added for the calculation of Lagrangian Descriptors.

    Parameters
    ----------
    slice_parameters : list of 3-tuples of floats
        input parameters of limits and size of mesh per axis

    Returns
    -------
    mesh : 1d numpy array
        flattened array of initial conditions
    """
    if type(grid_parameters) == dict:
        # Unpack extra grid parameters
        slice_parameters = grid_parameters['slice_parameters']
        dof_slice = grid_parameters['dof_slice']
        dof_fixed, dof_fixed_values = grid_parameters['dof_fixed']
        momentum_sign = grid_parameters['momentum_sign']
        potential_energy = grid_parameters['potential_energy']
        H0 = grid_parameters['energy_level']

        N_dim = len(dof_slice)  # Total number of DoF
        # Check N DoF is even.
        if N_dim % 2 != 0:
            error_mssg = ("ERROR: Number of DoF not even. ",
                          "Check your extra grid parameters")
            print(error_mssg)
            sys.exit()
        # Determine number of DoF for Energy conservation.
        # There must be only one.
        dof_remaining = 1 - (np.array(dof_fixed) + np.array(dof_slice))
        if list(dof_remaining).count(1) > 1:
            error_mssg = ("ERROR: More than one remaing DoF. ",
                          "Cannot use Energy conservation to define high-dim grid.")
            print(error_mssg)
            sys.exit()

        # Axes of visualisation slice
        def pass_axis_parameters(parameters): return np.linspace(*parameters)
        points_slice_axes = list(map(pass_axis_parameters, slice_parameters))
        slice_mesh = np.meshgrid(*points_slice_axes)
        dof_slice_axes = [axis.flatten() for axis in slice_mesh]

        # Axes of fixed DoF
        dof_fixed_axes = dof_fixed_values

        N_points_slice = np.prod(list(map(itemgetter(-1), slice_parameters)))
        phase_space_axes = {}

        # Define and sort axes in phase space
        k = 0
        for i in range(N_dim):
            if dof_slice[i] == 1:
                phase_space_axes[i] = dof_slice_axes[k]
                k += 1
        k = 0
        for i in range(N_dim):
            if dof_fixed[i] == 1:
                phase_space_axes[i] = dof_fixed_axes[k] * \
                    np.zeros(N_points_slice)
                k += 1

        # Set axis to be determined by energy conservation
        idx_dof_H0 = list(set(range(N_dim))-set(phase_space_axes.keys()))[0]
        phase_space_axes[idx_dof_H0] = np.zeros(N_points_slice)

        # List of all phase space axes
        phase_space_axes = [phase_space_axes[i] for i in range(int(N_dim))]

        # Determine undefined axis via energy conservation
        phase_space_axes[idx_dof_H0] = energy_conservation_condition(
            phase_space_axes, H0, potential_energy, momentum_sign)
        
        mask = np.isnan(phase_space_axes[idx_dof_H0]) # Mask grid points
        phase_space_axes[idx_dof_H0] = np.nan_to_num(phase_space_axes[idx_dof_H0])
        
        # Return array of mesh points for integrator
        lagrangian_descriptor_axis = [np.zeros(N_points_slice)]
        mesh = np.transpose(phase_space_axes + lagrangian_descriptor_axis)

        return mesh.flatten(), mask
    
    else:
        if len(grid_parameters) > 2:
            error_mssg = ("ERROR: grid_parameters must be a list for 2D slices for 1DoF systems. ",
                          "Provide a dictionary for a higher-dimensional systems.")
        else:
            x_min, x_max, Nx = grid_parameters[0]
            y_min, y_max, Ny = grid_parameters[1]
            points_x = np.linspace(x_min, x_max, Nx)
            points_y = np.linspace(y_min, y_max, Ny)
            X, Y = np.meshgrid(points_x, points_y)  # Grid in phase-space
            # 2D grid + a zero column for LDs
            mesh = np.transpose([X.flatten(), Y.flatten(), np.zeros(Nx*Ny)])
            mask = False
            
            return mesh.flatten(), mask

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
    N_dim = u.shape[-1]
    points_positions = u.T[:int(N_dim/2)]
    
    if len(points_positions) == len(box_boundaries):
        check = lambda x, box_axis_limits: (box_axis_limits[0]<=x)&(x<=box_axis_limits[1])
        positions_within_box = [check(points_positions[i], box_boundaries[i]) for i in range(int(N_dim/2))]
        u_indices = reduce(lambda x, y: x&y, positions_within_box)
        return u_indices
    else:
        error_mssg = ("ERROR: Number of box axes and configuration space axes do not match. ",
                      "Check the defintion of your box boundaries.")
        sys.exit()

def lagrangian_descriptor(u, v, p_value = 0.5):
    """
    Vector field equation for Lagrangian descriptor.

    Parameters
    ----------
    v : ndarray, shape(n,2)
        Vector field at given point.
            
    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.

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
        LD = np.zeros(len(u[:,0]))
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
    
    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.
    
    box_boundaries : list of 2-tuples, optional
        box boundaries for escape condition of variable time integration
        boundaries are infinite by default.
        
    Returns
    ------- 
    1d array
        y0 values for integrator 
    """
    N_mesh_axes = 2*len(box_boundaries)+1
    u = points.reshape((-1,N_mesh_axes))
    u = u[:,:-1] #remove LD-values axis
    
    # Apply Escape condition
    u_inbox = check_if_points_escape_box(u, box_boundaries)
    
    # Define output vector field in combination with escape condition
    v = np.zeros(u.shape)
    v[u_inbox == True] = vector_field(t, u[u_inbox == True])
    
    # Calculate LD vector field
    LD_vec = np.zeros(len(u))
    LD_vec [u_inbox == True] = lagrangian_descriptor(u[u_inbox == True], v[u_inbox == True], p_value)
    
    # Add LD
    v_out=np.column_stack((v, LD_vec))
    return v_out.flatten()

def compute_lagrangian_descriptor(grid_parameters, vector_field, tau, p_value=0.5, box_boundaries=False):
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
    #get visualisation slice parameters and Number of DoF
    if type(grid_parameters) == dict:
        #n-DoF systems
        slice_parameters = grid_parameters['slice_parameters'] # 2n-D grid
        N_dim = len(grid_parameters['dof_slice'])
    else:
        #1-DoF systems
        slice_parameters = grid_parameters # 2-D grid
        N_dim = len(slice_parameters)
        
    #set boundaries for escape-box condition, if not defined
    if not box_boundaries:
        box_boundaries = int(N_dim/2)*[[-np.infty, np.infty]] #restricted to configuration space
    
    #solve initial value problem
    f = lambda t, y: vector_field_flat(t, y, vector_field, p_value, box_boundaries)
    y0, mask = generate_points(grid_parameters)
    
    solution = solve_ivp(f, [0,tau], y0, t_eval=[tau], rtol=1.0e-4)

    LD_values = solution.y[N_dim::N_dim+1] #values corresponding to LD
    
    N_points_slice_axes = list( map(itemgetter(-1), slice_parameters)) 
    LD = np.abs(LD_values).reshape(*N_points_slice_axes) #reshape to 2-D array    
    LD = ma.masked_array(LD, mask=mask)
    
    if p_value<=1:
        return LD
    else:
        return LD**(1/p_value)

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
