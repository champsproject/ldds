import numpy as np

def HamCenter1D(t, u):
    """
    Returns 1D Hamilton-Centre vector field at time t, for an array of points in phase-space.
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase-space
        
    u : array_like, shape(n,)
        points in phase-space to determine vector field at common time t
        
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u in phase-space at time t
        functional form and parameters of the field is hard-coded.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    omega = 1
    v = np.array([ omega * y, - omega * x]).T
    return v

def HamSaddle1D(t, u):
    """
    Returns 1D Hamilton-Saddle vector field at time t, for an array of points in phase-space.
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase-space
        
    u : array_like, shape(n,)
        points in phase-space to determine vector field at common time t
        
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u in phase-space at time t
        functional form and parameters of the field is hard-coded.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    lamda = 1
    v = np.array([ lamda * y, lamda * x]).T
    return v

def HamDuffing1D(t, u):
    """
    Returns perturbed 1D Hamilton-Duffing vector field at time t, for an array of points in phase-space.
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase-space
        
    u : array_like, shape(n,)
        points in phase-space to determine vector field at common time t
        
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u in phase-space at time t
        functional form and parameters of the field is hard-coded.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    v = np.array([y, x - x**3]).T
    return v

def HamSN1D(t, u):
    """
    Returns 1D Hamilton-Saddle-Node vector field at time t, for an array of points in phase-space.
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase-space
        
    u : array_like, shape(n,)
        points in phase-space to determine vector field at common time t
        
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u in phase-space at time t
        functional form and parameters of the field is hard-coded.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    v = np.array([ y, -x -x**2]).T
    return v

def HamSN1D_inverted(t, u):
    """
    Returns 1D Inverted Hamilton-Duffing vector field at time t, for an array of points in phase-space.
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase-space
        
    u : array_like, shape(n,)
        points in phase-space to determine vector field at common time t
        
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u in phase-space at time t
        functional form and parameters of the field is hard-coded.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    # perturbation = forcing(t, u, flag_pert, perturbation_params)
    # v = np.array([y, - x + x**3 + perturbation]).T
    v = np.array([y, - x + x**3]).T
    return v

def forcing(t, u, perturbation_type = 1, perturbation_params = [0.15, 0.5]):
    x, y = u.T
    perturbation = np.zeros(u.shape)
    
    # Perturbation parameters
    A, freq = perturbation_params # Amplitude and Frequency
    
    if perturbation_type == 1:
        perturbation = perturbation + np.array([0, A * np.sin(freq*t)])
    elif perturbation_type == 2:
        perturbation = perturbation + np.array([0, A * np.sech(t) * np.sin(freq*t)])
    elif perturbation_type == None:
        perturbation = perturbation
    
    return perturbation

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido'
__status__ = 'Development'
