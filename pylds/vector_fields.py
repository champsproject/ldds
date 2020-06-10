"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""
import numpy as np

def HamCenter1D(t, u, PARAMETERS = [1]):
    """
    Returns 1D Hamilton-Centre vector field at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega]
    Functional form: v = (omega*y, - omega*x), with u = (x, y)
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase space
        
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t
    
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    omega, = PARAMETERS
    v = np.array([ omega * y, - omega * x]).T
    return v

def HamSaddle1D(t, u, PARAMETERS = [1]):
    """
    Returns 1D Hamilton-Saddle vector field at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [lamda]
    Functional form: v = (lamda*y, - lamda*x), with u = (x, y)
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase space.
        
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t.
    
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    lamda, = PARAMETERS
    v = np.array([ lamda * y, lamda * x]).T
    return v

def Duffing1D(t, u, PARAMETERS = [None]):
    """
    Returns 1D vector field of the Duffing oscillator, for an array of points in phase space.
    Number of model parameters: 0 . PARAMETERS = [None]
    Functional form: v = (y, x - x**3), with u = (x, y)
    
    Parameters
    ----------  
    t : float
        fixed time-point of vector field, for all points in phase space.
        
    u : array_like, shape(n,)
        Points in phase space.
        
    PARAMETERS : list of floats
        Vector field parameters.
    
    Returns
    -------
    v : array_like, shape(n,)
        Vector field corresponding to points u, in phase space at time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    v = np.array([y, x - x**3]).T
    return v

def HamSN1D(t, u, PARAMETERS = [None]):
    """
    Returns 1D Hamilton-Saddle-Node vector field at time t, for an array of points in phase space.
    Number of model parameters: 0 . PARAMETERS = [None]
    Functional form: v = (y, -x -x**2), with u = (x, y)
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase space.
        
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t.
        
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    v = np.array([ y, -x -x**2]).T
    return v

def HamSN1D_inverted(t, u, PARAMETERS = [None]):
    """
    Returns 1D Inverted Hamilton-Duffing vector field at time t, for an array of points in phase space.
    Number of model parameters: 0 . PARAMETERS = [None]
    Functional form: v = (y, - x + x**3), with u = (x, y)
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase space.
        
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t.
        
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    # perturbation = forcing(t, u, flag_pert, perturbation_params)
    # v = np.array([y, - x + x**3 + perturbation]).T
    v = np.array([y, - x + x**3]).T
    return v

def forcing(t, u, perturbation_params = [1, 0.15, 0.5]):
    """
    Returns 1D Inverted Hamilton-Duffing vector field at time t, for an array of points in phase space.
    Number of model parameters: 3. PARAMETERS = [perturbation_type, A, freq]
    Functional form: v = (, ), with u = (x, y)
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase space.
        
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t.
        
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x, y = u.T
    perturbation = np.zeros(u.shape)
    
    # Perturbation parameters
    perturbation_type, A, freq = perturbation_params # Amplitude and Frequency
    
    if perturbation_type == 1:
        perturbation = perturbation + np.array([0, A * np.sin(freq*t)])
    elif perturbation_type == 2:
        perturbation = perturbation + np.array([0, A * np.sech(t) * np.sin(freq*t)])
    elif perturbation_type == None:
        perturbation = perturbation
    
    return perturbation

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido'
__status__ = 'Development'
