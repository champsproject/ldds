"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""

import sys
import numpy as np
import numpy.ma as ma
from scipy.integrate import solve_ivp
from functools import partial
from multiprocessing import Pool,cpu_count

def energy_conservation_condition(phase_space_axes, H0, potential, momentum_sign):
    N_dim = len(phase_space_axes)
    positions = np.transpose(phase_space_axes[:int(N_dim/2)])
    momenta = np.transpose(phase_space_axes[int(N_dim/2):])
    V = potential(positions)
    point_dim_remaining = momentum_sign * \
        np.sqrt(2*(H0 - V) - (momenta**2).sum(axis=1))
    return point_dim_remaining

def generate_grid(grid_parameters):
    """
    Returns a 1D array of all points from a on a uniform grid with dimensions and size defined by list of input parameters.
    An additional dimension initiallised with zeros is added for the calculation of Lagrangian Descriptors.
    NOTE: For n-DoF systems, currently energy conservation is only used to determine momenta dimensions.

    Parameters
    ----------
    grid_parameters : list (1-DoF systems) or dict (n-DoF systems)
        if 1-DoF, list should have two 3-tuples of floats
        entries are input parameters of limits and size of mesh per axis

        if n-DoF, dict should have the following keys
        * 'slice_parameters' : list, should have two 3-tuples of floats, for a 2D slice
        * 'dims_slice' : list of 0 and 1, ones indicate slice axes
        * 'dims_fixed' : list of 0 and 1, ones indicate fixed axes
        * 'momentum_sign' : int, -1 / 1, for negative/positive momentum for remaining axis
        * 'potential_energy' : func used by energy conservation condition to determine remaining momentum axis
        * 'energy_level' : float, energy value for energy conservation condition
    Returns
    -------
    mesh : 1d numpy array
        flattened array of initial conditions
    """
    if type(grid_parameters) == dict:
        # Unpack extra grid parameters
        slice_parameters = grid_parameters['slice_parameters']
        dims_slice = grid_parameters['dims_slice']
        dims_fixed, dims_fixed_values = grid_parameters['dims_fixed']
        momentum_sign = grid_parameters['momentum_sign']
        potential_energy = grid_parameters['potential_energy']
        H0 = grid_parameters['energy_level']

        N_dim = len(dims_slice)  # Phase space dimensions

        # Check N DoF is even.
        if N_dim % 2 != 0:
            error_mssg = ("Error: Number of phase space dimensions not even. ")
            print(error_mssg)
            sys.exit()

        # Determine number of dimensions for energy conservation.
        # There must be only one.
        dims_remaining = 1 - (np.array(dims_fixed) + np.array(dims_slice))
        if list(dims_remaining).count(1) > 1:
            error_mssg = ("Error: More than one remaing dimension. ",
                          "Cannot use energy conservation to define high dim grid.")
            print(error_mssg)
            sys.exit()

        # Axes of visualisation slice
        points_slice_axes = list(map(np.linspace, *slice_parameters.T))
        slice_mesh = np.meshgrid(*points_slice_axes)
        dims_slice_axes = [axis.flatten() for axis in slice_mesh]

        N_points_slice = np.prod(slice_parameters[:,-1])
        phase_space_axes = {}

        # Define and sort axes in phase space
        k = 0
        for i in range(N_dim):
            if dims_slice[i] == 1:
                phase_space_axes[i] = dims_slice_axes[k]
                k += 1
        k = 0
        for i in range(N_dim):
            if dims_fixed[i] == 1:
                phase_space_axes[i] = dims_fixed_values[k] * \
                    np.ones(N_points_slice)
                k += 1

        # Set axis to be determined by energy conservation
        idx_dims_H0 = list(set(range(N_dim))-set(phase_space_axes.keys()))[0]

        # Check if remaining dimension falls in configuration space
        if idx_dims_H0 < int(N_dim/2):
            error_mssg = ("Error: The remaining dimension fall in configuration space.",
                          "Cannot do fixed-momentum slices (yet).")
            print(error_mssg)
            sys.exit()

        phase_space_axes[idx_dims_H0] = np.zeros(N_points_slice)

        # List of all phase space axes
        phase_space_axes = [phase_space_axes[i] for i in range(N_dim)]

        # Determine undefined axis via energy conservation
        phase_space_axes[idx_dims_H0] = energy_conservation_condition(
            phase_space_axes, H0, potential_energy, momentum_sign)

        phase_space_axes[idx_dims_H0] = np.nan_to_num(phase_space_axes[idx_dims_H0])

        # Return array of mesh points for integrator
        lagrangian_descriptor_axis = [np.zeros(N_points_slice)]
        mesh = np.transpose(phase_space_axes + lagrangian_descriptor_axis)

        return mesh

    else:
        if len(grid_parameters) > 2:
            error_mssg = ("Error: Provide grid parameters in a dictionary for systems with 2 and more DoF.")
        else:
            x_min, x_max, Nx = grid_parameters[0]
            y_min, y_max, Ny = grid_parameters[1]
            points_x = np.linspace(x_min, x_max, Nx)
            points_y = np.linspace(y_min, y_max, Ny)
            X, Y = np.meshgrid(points_x, points_y)  # Grid in phase space
            # 2D grid + a zero column for LDs
            mesh = np.transpose([X.flatten(), Y.flatten(), np.zeros(Nx*Ny)])

            return mesh

def vector_field_ld(vector_field, p_value, t, u):
    """
    Combines vector_field with equation for Lagrangian descriptor.

    Parameters
    ----------
    t : float
        time

    u : ndarray, shape(n)
        point in phase space.

    vector_field: function
        User defined vector field.

    p_value : float, optional
        Exponent in Lagrangian descriptor definition.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.

    Returns
    -------
    dudt : ndarray, shape(n)
        vector field with LD
    """
    dudt = vector_field(t, u[:-1])
    if p_value == 0:
        n = 0.5*(len(u)-1)
        LD = np.abs(u[:n]*dudt[n:])
    elif 2 >= p_value > 0:
        LD = np.sum(np.abs(dudt)**p_value)
    return np.append(dudt, LD)

def compute_ld(y_in, vector_field, tau, p_value):
    """
    Returns the LD value for a single initial condition.

    Parameters
    ----------
    y_in : ndarray, shape(n)
        initial condition.

    vector_field: function(t,y)
        Vector field of the system.

    tau : float
        Integration time.

    p_value : float, optional
        Lagrangian descriptor parameter.

    Returns
    -------
    LD : float
        LD value.
    """
    # vec = partial(vector_field_ld,vector_field, p_value)

    def vec (t, y):
        return vector_field_ld(vector_field, p_value, t, y)

    #combine vector field with LDs
    solution = solve_ivp(vec, [0,tau], y_in, t_eval=[tau], rtol=1.0e-4)

    LD = np.abs(solution.y[-1]) #values corresponding to LD

    return LD

def compute_ld_box(y_in, vector_field, tau, p_value, box_boundaries):
    """
    Returns the LD value for a single initial condition.

    Parameters
    ----------
    y_in : ndarray, shape(n)
        initial condition.

    vector_field: function(t,y)
        Vector field of the system.

    tau : float
        Integration time.

    p_value : float, optional
        Lagrangian descriptor parameter.

    box_boundaries : list of 2-tuples, optional
        Box boundaries for escape condition of variable time integration.

    Returns
    -------
    LD : float
        LD value.
    """
    # vec = partial(vector_field_ld,vector_field, p_value)

    def vec (t, y):
        return vector_field_ld(vector_field, p_value, t, y)

    def event_function(t,y):
        positions = y[:len(y)/2]
        upper_lim = np.array(box_boundaries)[:,1]
        lower_lim = np.array(box_boundaries)[:,0]
        return np.all(upper_lim > positions) and np.all(positions > lower_lim)
    solution = solve_ivp(vec, [0,tau], y_in, t_eval=[tau], events=event_function, rtol=1.0e-4)


    LD = np.abs(solution.y[-1]) #values corresponding to LD

    return LD

def compute_lagrangian_descriptors(grid_parameters, vector_field, tau, p_value=0.5, box_boundaries=np.nan):
    """
    Returns Lagrangian descriptor values for system defined by vector_field on a uniform grid defined by _grid_parameters.

    Parameters
    ----------
    grid_parameters : list of 3-tuples of floats
        Input parameters for grid.

    vector_field: function(t,y)
        Vector field of the system.

    tau : float
        Integration time.

    p_value : float, optional
        Lagrangian descriptor parameter.
        0 is the acton-based LD,
        0 < p_value < 1 is the Lp quasinorm,
        1 <= p_value < 2 is the Lp norm LD,
        2 is the arclength LD.
        The default is 0.5.

    box_boundaries : list of 2-tuples, optional
        Box boundaries for escape condition of variable time integration.
        Boundaries are off by default.

    Returns
    -------
    LD : ndarray, shape (Nx, Ny)
        LD values.
    """

    if type(grid_parameters) == dict:
        #n-DoF systems
        slice_parameters = grid_parameters['slice_parameters'] # 2D grid
    else:
        #1-DoF systems
        slice_parameters = grid_parameters # 2D grid

    y_grid = generate_grid(grid_parameters)

    with Pool(cpu_count()) as p:
        if np.any(np.isnan(box_boundaries)):
            y_init = [[y, vector_field, tau, p_value] for y in y_grid]
            LD=p.starmap(compute_ld, y_init)
        else:
            y_init = [[y, vector_field, tau, p_value, box_boundaries] for y in y_grid]
            LD=p.starmap(compute_ld_box, y_init)

    LD_mat = np.abs(LD).reshape([s[-1] for s in slice_parameters])

    if p_value<=1:
        return LD_mat
    else:
        return LD_mat**(1/p_value)

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
