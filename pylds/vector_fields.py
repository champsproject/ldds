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

def Duffing1D_inverted(t, u, PARAMETERS = [None]):
    """
    Returns 1D vector field of the inverted Duffing oscillator at time t, for an array of points in phase space.
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

def forcing(t, u, perturbation_params = [1, 0.15, 0.5]):
    """
    Returns vector field for a perturbation at time t, for an array of points in phase space.
    Number of model parameters: 3. perturbation_params = [perturbation_type, amplitude, frequency]
    Functional form: v = (, ), with u = (x, y)
    
    Parameters
    ----------
    t : float
        fixed time-point of vector field, for all points in phase space.
        
    u : array_like, shape(n,)
        points in phase space to determine vector field at time t.
        
    perturbation_params : list of floats, [perturbation_type, amplitude, frequency]
        vector field parameters
    
    Returns
    -------
    v : array_like, shape(n,)
        vector field corresponding to points u, in phase space at time t
    """
    x, y = u.T
    perturbation = np.zeros(u.shape)
    
    # Perturbation parameters
    perturbation_type, amplitude, freq = perturbation_params
    
    if perturbation_type == 1:
        perturbation = perturbation + np.array([0, amplitude * np.sin(freq*t)])
    elif perturbation_type == 2:
        perturbation = perturbation + np.array([0, amplitude * np.sech(t) * np.sin(freq*t)])
    
    return perturbation

def HenonHeiles_vector_field(t, u, PARAMETERS = None):
    """
    Returns 2D Henon-Heiles vector field at time t, for an array of points in phase space.
    Functional form: v = (p_x, p_y, -x - 2*x*y, -x**2 -y + y**2), with u = (x, y, p_x, p_y)
    
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
    N_dim = u.shape[-1]
    points_positions = u.T[:int(N_dim/2)]
    points_momenta = u.T[int(N_dim/2):]
    x, y = points_positions
    p_x, p_y = points_momenta 
    
    # Hamiltonian Model Parameter
    # None
    
    # Vector field defintion
    v_x   =  p_x
    v_y   =  p_y
    v_p_x = -x - 2*x*y
    v_p_y = -x**2 -y + y**2
    v = np.array([v_x, v_y, v_p_x, v_p_y]).T
    return v

def HenonHeiles_potential(positions, PARAMETERS = None):
    """
    Potential Energy Function (PEF) of 2DoF Henon-Heiles system.

    Parameters
    ----------
        
    positions : array_like, shape(n,)
        array of 2D points in configuration space.
        
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    V : array_like, shape(n,)
        potential energy corresponding to all positions
    
    """
    x, y = positions.T
    # Function parameters
    # None
    # Potential energy function
    V = (1/2)*(x**2 + y**2) + (y * x**2) - (1/3)*y**3
    return V

def NFSaddle_vector_field(t, u, PARAMETERS = None):
    """
    Returns 2D Index-1 Normal-Form Saddle vector field at time t, for an array of points in phase space.
    Functional form: v = (p_x, p_y, x, -y), with u = (x, y, p_x, p_y)
    
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
    N_dim = u.shape[-1]
    points_positions = u.T[:int(N_dim/2)]
    points_momenta = u.T[int(N_dim/2):]
    x, y = points_positions
    p_x, p_y = points_momenta 
    
    # Hamiltonian Model Parameter
    # None
    
    # Vector field defintion
    v_x   =  p_x
    v_y   =  p_y
    v_p_x = x
    v_p_y = -y
    v = np.array([v_x, v_y, v_p_x, v_p_y]).T
    return v

def NFSaddle_potential(positions, PARAMETERS = None):
    """
    Potential Energy Function (PEF) of 2DoF Index-1 Normal-Form Saddle system.

    Parameters
    ----------
        
    positions : array_like, shape(n,)
        array of 2D points in configuration space.
        
    PARAMETERS : list of floats
        vector field parameters
    
    Returns
    -------
    V : array_like, shape(n,)
        potential energy corresponding to all positions
    
    """
    x, y = positions.T
    # Function parameters
    # None
    
    # Potential energy function
    V = (1/2)*y**2 - (1/2)*x**2
    return V

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
