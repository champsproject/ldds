import numpy as np

def HamCenter1D(t, u, PARAMETERS = [1]):
    """
    Returns vector field for a 1DoF centre at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [omega]
    Functional form: v = (omega*y, - omega*x), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    omega, = PARAMETERS
    v = np.column_stack([ omega * y, - omega * x])
    return v

def Saddle1D(t, u, PARAMETERS = [1, 1]):
    """
    Returns 1D Hamilton-Saddle vector field at time t, for an array of points in phase space.
    Number of model parameters: 2 .PARAMETERS = [a1, a2]
    Functional form: v = (a1*x, -a2*y), with u = (x, y)

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
    a1, a2 = PARAMETERS
    v = np.column_stack([ a1*x, -a2*y])
    return v

def HamSaddle1D(t, u, PARAMETERS = [1]):
    """
    Returns vector field for a 1DoF saddle at time t, for an array of points in phase space.
    Number of model parameters: 1 . PARAMETERS = [lamda]
    Functional form: v = (lamda*y, - lamda*x), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    lamda, = PARAMETERS
    v = np.column_stack([ lamda * y, lamda * x])
    return v

def Duffing1D(t, u, PARAMETERS = [1, 1, -1]):
    """
    Returns vector field for the Duffing oscillator.
    Number of model parameters: 0 . PARAMETERS = [None]
    Functional form: v = (y, x - x**3), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t..
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    alpha, beta, gamma = PARAMETERS
    v = np.column_stack([ alpha*y, beta*x + gamma*x**3])
    return v

def Duffing1D_inverted(t, u, PARAMETERS = [None]):
    """
    Returns vector field for the inverted Duffing oscillator.
    Number of model parameters: 0 . PARAMETERS = [None]
    Functional form: v = (y, - x + x**3), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    # perturbation = forcing(t, u, flag_pert, perturbation_params)
    # v = np.array([y, - x + x**3 + perturbation]).T
    v = np.column_stack([y, - x + x**3])
    return v

def HamSN1D(t, u, PARAMETERS = [None]):
    """
    Returns vector field for the 1DoF saddle-node model.
    Number of model parameters: 0 . PARAMETERS = [None]
    Functional form: v = (y, -x -x**2), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    x, y = u.T
    # Hamiltonian Model Parameter
    v = np.column_stack([ y, -x -x**2])
    return v

def forcing(t, u, perturbation_params = [1, 0.15, 0.5]):
    """
    Returns vector field for a perturbation.
    Number of model parameters: 3. perturbation_params = [perturbation_type, amplitude, frequency]
    Functional form: v = (, ), with u = (x, y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    perturbation_params : list of floats, [perturbation_type, amplitude, frequency]
        Perturbation parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
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

def HenonHeiles_vector_field(t, u):
    """
    Returns Henon-Heiles vector field (2DoF).
    Functional form: v = (p_x, p_y, -x - 2*x*y, -x**2 -y + y**2), with u = (x, y, p_x, p_y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
    """
    points_positions = u.T[:2]
    points_momenta = u.T[2:4]
    x, y = points_positions
    p_x, p_y = points_momenta

    # Vector field defintion
    v_x   =  p_x
    v_y   =  p_y
    v_p_x = -x - 2*x*y
    v_p_y = -x**2 -y + y**2
    v = np.column_stack([v_x, v_y, v_p_x, v_p_y])
    return v

def NFSaddle_vector_field(t, u, PARAMETERS = None):
    """
    Returns vector field for a 2D index-1 saddle.
    Functional form: v = (p_x, p_y, x, -y), with u = (x, y, p_x, p_y)

    Parameters
    ----------
    t : float
        Time. (This vector field is independent of time.)

    u : ndarray, shape(n,)
        Points in phase space.

    PARAMETERS : list of floats
        Vector field parameters.

    Returns
    -------
    v : ndarray, shape(n,)
        Vector field at points u and time t.
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
    v = np.column_stack([v_x, v_y, v_p_x, v_p_y])
    return v


__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak, Shibabrat Naik'
__status__ = 'Development'
