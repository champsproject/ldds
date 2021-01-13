import numpy as np

def StandardMap(u_initial, PARAMETERS=[0.3]):
    """
    2D Standard map for initial conditions in a unit square, centred at the origin (James Meiss).
    The map will return unwrapped trajectories for iterations of initial conditions, unlike when using PBCs.
    For computation of the Lagrangian Descriptor relative displacements are only needed. 
    
    Number of model parameters: 1 . PARAMETERS = [K]
    Functional form: 
    x_next  = x_initial + y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial
    y_next  = y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial))
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : array_like, shape(n,)
        initial points in unit square to determine their next iteration under the map
    
    PARAMETERS : list of floats
        map parameters
    
    Returns
    -------
    u_next : array_like, shape(n,)
        points u_next, in the 2D plane (not in the unit square necessarily).
    """
    x_initial, y_initial = u_initial.T
    # Map parameters
    K, = PARAMETERS
    
    # Map components 
    x_next = x_initial + y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial)
    y_next = y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial)
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])
    return u_next

def StandardMap_inverse(u_initial, PARAMETERS=[0.3]):
    """
    Inverse of 2D Standard map for initial conditions in a unit square, centred at the origin (James Meiss).
    The map will return unwrapped trajectories for iterations of initial conditions, unlike when using PBCs.
    For computation of the Lagrangian Descriptor relative displacements are only needed. 
    
    Number of model parameters: 1 . PARAMETERS = [K]
    Functional form: 
    x_next = x_initial - y_initial
    y_next = y_initial + (K/(2*np.pi))*np.sin(2*np.pi*(x_initial - y_initial))
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : array_like, shape(n,)
        initial points in unit square to determine their next iteration under the map
    
    PARAMETERS : list of floats
        map parameters
    
    Returns
    -------
    u_next : array_like, shape(n,)
        points u_next, in the 2D plane.
    """
    x_initial, y_initial = u_initial.T
    # Map parameters
    K, = PARAMETERS
    
    # Map components 
    x_next = x_initial - y_initial
    y_next = y_initial + (K/(2*np.pi))*np.sin(2*np.pi*(x_initial - y_initial))
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])

    return u_next

def HenonMap(u_initial, PARAMETERS=[0.298, 1]):
    """
    2D Henon map. 
    
    Number of model parameters: 2 . PARAMETERS = [a, b]
    Functional form: 
    x_next = a - x_initial^2 + b*y_initial
    y_next = x_initial
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : array_like, shape(n,)
        initial points in unit square to determine their next iteration under the map
    
    PARAMETERS : list of floats
        map parameters
    
    Returns
    -------
    u_next : array_like, shape(n,)
        points u_next, in the 2D plane.
    """
    x_initial, y_initial = u_initial.T
    # Map parameters
    a, b = PARAMETERS
    
    # Map components
    x_next = a - x_initial**2 + b*y_initial
    y_next = x_initial
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])
    
    return u_next

def HenonMap_inverse(u_initial, PARAMETERS=[0.298, 1]):
    """
    Inverse of 2D Henon map.
    
    Number of model parameters: 2 . PARAMETERS = [a, b]
    Functional form: 
    x_next = y_initial
    y_next = (x_initial - a + y_initial^2)/b
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : array_like, shape(n,)
        initial points in unit square to determine their next iteration under the map
    
    PARAMETERS : list of floats
        map parameters
    
    Returns
    -------
    u_next : array_like, shape(n,)
        points u_next, in the 2D plane (not in the unit square necessarily)
    """
    x_initial, y_initial = u_initial.T
    # Map parameters
    a, b = PARAMETERS
    
    # Map components
    x_next = y_initial
    y_next = (x_initial - a + y_initial**2)/b
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])
    
    return u_next
