
import numpy as np

def StandardMap(u_initial, PARAMETERS=[0.3]):
    """
    Returns 1D Hamilton-Centre vector field at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega]
    Functional form: v = (omega*y, - omega*x), with u = (x, y)
    
    Parameters
    ----------
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t
    
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    u_initial = np.mod(u_initial, 1)
    x_initial, y_initial = u_initial.T
    # Hamiltonian Model Parameter
    K, = PARAMETERS
    
    # Map components 
    x_next = x_initial + y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial)
    y_next = y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial)
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])
    return u_next

def StandardMap_inverse(u_initial, PARAMETERS=[0.3]):
    """
    Returns 1D Hamilton-Centre vector field at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega]
    Functional form: v = (omega*y, - omega*x), with u = (x, y)
    
    Parameters
    ----------
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t
    
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x_initial, y_initial = u_initial.T
    # Hamiltonian Model Parameter
    K, = PARAMETERS
    
    # Map components 
    x_next = x_initial - y_initial
    y_next = y_initial + (K/(2*np.pi))*np.sin(2*np.pi*(x_initial - y_initial))
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])

    return u_next

def HenonMap(u_initial, PARAMETERS=[0.298, 1]):
    """
    Returns 1D Hamilton-Centre vector field at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega]
    Functional form: v = (omega*y, - omega*x), with u = (x, y)
    
    Parameters
    ----------
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t
    
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x_initial, y_initial = u_initial.T
    # Hamiltonian Model Parameter
    a, b = PARAMETERS
    
    # Map components
    x_next = a - x_initial**2 + b*y_initial
    y_next = x_initial
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])
    
    return u_next

def HenonMap_inverse(u_initial, PARAMETERS=[0.298, 1]):
    """
    Returns 1D Hamilton-Centre vector field at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega]
    Functional form: v = (omega*y, - omega*x), with u = (x, y)
    
    Parameters
    ----------
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t
    
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x_initial, y_initial = u_initial.T
    # Hamiltonian Model Parameter
    a, b = PARAMETERS
    
    # Map components
    x_next = y_initial
    y_next = (x_initial - a + y_initial**2)/b
    
    # Map next iteration
    u_next = np.column_stack([ x_next, y_next])
    
    return u_next
