import numpy as np

def StandardMap(u_initial, PARAMETERS=[0.3]):
    """
    Chirikov standard map for initial conditions in a unit square, centred at the origin (as used in [1]_).
    The map is defined in a perioidic manner, but maps points outside the unit square (periodic domain).
    The Lagrangian descriptor calculations use a correction, so that points are mapped into the unit square. 
    
    Number of model parameters: 1 . PARAMETERS = [K]
    Functional form: 
    x_next  = x_initial + y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial
    y_next  = y_initial - (K/(2*np.pi))*np.sin(2*np.pi*x_initial))
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : ndarray, shape(n,),
        Points to be mapped.
    
    PARAMETERS : list of floats
        Map parameters.
    
    Returns
    -------
    u_next : ndarray, shape(n,)
        Array of forward images of u under Chirikov standard map (not in the unit square necessarily).
        
    References
    ----------
    .. [1] J.D. Meiss Visual explorations of dynamics: The standard map. Pramana - J Phys 70, 965–988 (2008). https://doi.org/10.1007/s12043-008-0103-3
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
    Inverse of Chirikov standard map for initial conditions in a unit square, centred at the origin (as used in [1]_).
    The map is defined in a perioidic manner, but maps points outside the unit square (periodic domain).
    The Lagrangian descriptor calculations use a correction, so that points are mapped into the unit square.
    
    Number of model parameters: 1 . PARAMETERS = [K]
    Functional form: 
    x_next = x_initial - y_initial
    y_next = y_initial + (K/(2*np.pi))*np.sin(2*np.pi*(x_initial - y_initial))
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : ndarray, shape(n,),
        Points to be mapped.
    
    PARAMETERS : list of floats
        Map parameters.
    
    Returns
    -------
    u_next : ndarray, shape(n,)
        Array of backward images of u under Chirikov standard map (not in the unit square necessarily).
        
    References
    ----------
    .. [1] J.D. Meiss Visual explorations of dynamics: The standard map. Pramana - J Phys 70, 965–988 (2008). https://doi.org/10.1007/s12043-008-0103-3
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
    Henon map. 
    
    Number of model parameters: 2 . PARAMETERS = [a, b]
    Functional form: 
    x_next = a - x_initial^2 + b*y_initial
    y_next = x_initial
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : ndarray, shape(n,),
        Points to be mapped.
    
    PARAMETERS : list of floats
        Map parameters.
    
    Returns
    -------
    u_next : ndarray, shape(n,)
        Array of forward images of u under Henon map.
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
    Inverse of Henon map.
    
    Number of model parameters: 2 . PARAMETERS = [a, b]
    Functional form: 
    x_next = y_initial
    y_next = (x_initial - a + y_initial^2)/b
    
    where u_initial = (x_initial, y_initial) 
    
    Parameters
    ----------
    u_initial : ndarray, shape(n,),
        Points to be mapped.
    
    PARAMETERS : list of floats
        Map parameters.
    
    Returns
    -------
    u_next : ndarray, shape(n,)
        Array of backward images of u under Henon map.
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
