"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""
import numpy as np

def HamCenter1(t, u, omega = 1):
    """
    Returns vector field of 1DoF Hamiltonian centre  at (t, u).
    Functional form: v = (omega*y, - omega*x), with u = (x, y).

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n)
        point in phase space.

    omega : float
        vector field parameter

    Returns
    -------
    dudt : ndarray, shape(n)
        vector field at (t, u)
    """
    dudt = np.zeros(len(u))

    dudt[0] = omega * u[1]
    dudt[1] = - omega * u[0]
    return dudt

def HamSaddle1(t, u, Lambda = 1):
    """
    Returns vector field of 1DoF Hamiltonian saddle at (t, u).
    Functional form: v = (Lambda*y, - Lambda*x), with u = (x, y).

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n)
        point in phase space.

    Lambda : float
        vector field parameter

    Returns
    -------
    dudt : ndarray, shape(n)
        vector field at (t, u)
    """
    dudt = np.zeros(len(u))

    dudt[0] = Lambda * u[1]
    dudt[1] = Lambda * u[0]
    return dudt
    # return Lambda * u[::-1]

def Duffing1(t, u):
    """
    Returns vector field of the 1DoF Duffing oscillator at (t,u).
    Functional form: v = (y, x - x**3), with u = (x, y).

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n)
        point in phase space.

    Returns
    -------
    dudt : ndarray, shape(n)
        vector field at (t, u).
    """
    dudt = np.zeros(len(u))

    dudt[0] = u[1]
    dudt[1] = u[0] - u[0]**3
    return dudt

def Duffing1_inverted(t, u):
    """
    Returns vector field of the inverted 1DoF Duffing oscillator at (t, u).
    Functional form: v = (y, - x + x**3), with u = (x, y).

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n)
        point in phase space.

    Returns
    -------
    dudt : ndarray, shape(n)
        vector field at (t, u)
    """
    dudt = np.zeros(len(u))

    dudt[0] = u[1]
    dudt[1] = - u[0] + u[0]**3
    return dudt

def HamSN1(t, u):
    """
    Returns vector field of the 1DoF Hamiltonian saddle-centre at (t, u).
    Functional form: v = (y, -x -x**2), with u = (x, y).

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n)
        point in phase space.

    Returns
    -------
    dudt : ndarray, shape(n)
        vector field at (t, u)
    """
    dudt = np.zeros(len(u))

    dudt[0] = u[1]
    dudt[1] = - u[0] - u[0]**2
    return dudt

def forcing1(vector_field, perturbation_params = [1, 0.15, 0.5], *args):
    """
    Returns perturbed 1DoF vector field at (t,u).
    Number of model parameters: 3. perturbation_params = [perturbation_type, amplitude, frequency].

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n)
        point in phase space.

    parameter : float
        vector field parameter

    vector_field : function
        vector field to be perturbed

    perturbation_params : list of floats, [perturbation_type, amplitude, frequency]

    Returns
    -------
    perturbation : ndarray, shape(n)
        perturbed vector field
    """

    # Perturbation parameters
    perturbation_type, amplitude, freq = perturbation_params

    if perturbation_type == 1:
        perturbation = lambda t: np.array([0, amplitude * np.sin(freq*t)])
    elif perturbation_type == 2:
        perturbation = lambda t: np.array([0, amplitude * np.sech(t) * np.sin(freq*t)])
    else:
        perturbation = lambda t: np.array([0,0])

    return lambda t, u :vector_field(t, u, *args) + perturbation(t)

def HenonHeiles_vector_field(t, u):
    """
    Returns vector field of the 2DoF Hénon-Heiles system at (t, u).
    Functional form: v = (p_x, p_y, -x - 2*x*y, -x**2 -y + y**2), with u = (x, y, p_x, p_y).

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n,)
        point in phase space.

    Returns
    -------
    v : ndarray, shape(n,)
        vector field at (t, u)
    """
    dudt = np.zeros(len(u))

    dudt[0] = u[2]
    dudt[1] = u[3]
    dudt[2] = - u[0] - 2*u[0]*u[1]
    dudt[3] = - u[0]**2 - u[1] + u[1]**2
    return dudt

def HenonHeiles_potential(x):
    """
    Potential Energy of the 2DoF Hénon-Heiles system.

    Parameters
    ----------

    x : ndarray, shape(n)
        array of positions only.

    Returns
    -------
    V : float
        potential energy at x

    """
    V = 0.5*(x[0]**2 + x[1]**2) + (x[1] * x[0]**2) - (1/3)*x[1]**3
    return V

def Saddle2_vector_field(t, u):
    """
    Returns vector field of the 2DoF index-1 saddle at (t, u).
    Functional form: v = (p_x, p_y, x, -y), with u = (x, y, p_x, p_y).

    Parameters
    ----------
    t : float
        time.

    u : ndarray, shape(n)
        point in phase space.

    Returns
    -------
    v : ndarray, shape(n)
        vector field at (t, u)
    """
    dudt = np.zeros(len(u))

    dudt[0] = u[2]
    dudt[1] = u[3]
    dudt[2] = u[0]
    dudt[3] = - u[1]
    return dudt

def Saddle2_potential(x):
    """
    Potential Energy of the 2DoF index-1 saddle.

    Parameters
    ----------

    x : ndarray, shape(n)
        array of positions only.

    Returns
    -------
    V : float
        potential energy at x

    """
    return 0.5*(x[1]**2 - x[0]**2)

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
