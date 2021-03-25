Examples
========

Maps
----

1. Standard map

The standard map (kicked rotator) is a two-dimensional map used in dynamical systems to study a periodically 
kicked pendulum. Its equations are given by the expressions:

.. math::
   \begin{align}
    x_{n+1} &= x_{n} + y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n))
    y_{n+1} &= y_{n} - \dfrac{K}{2\pi} \sin(2\pi x_{n))
   \end{align}
where :math:`K` is the parameter that controls the forcing strength of the perturbation.
   
The inverse map is described by:

.. math::
   \begin{align}
    x_{n} = x_{n+1} - y_{n+1}
    y_{n} = y_{n+1} + \dfrac{K}{2\pi} \sin(2\pi (x_{n+1} - y_{n+1}))
   \end{align}

   
2. Hénon map

The Hénon map was introduced by Michel Hénon as a simplified model of the Poincaré section 
of the Lorenz model. The map equations are as follows:

.. math::
   \begin{align}
    x_{n+1} = a - x_{n}^2 + b y_{n}
    y_{n+1} = x_{n}
   \end{align}

The inverse Hénon map is:

.. math::
   \begin{align}
    x_{n} = y_{n+1}
    y_{n} = \dfrac{x_{n+1} - a + y_{n+1}^2}{b}
   \end{align}


Two dimensional phase space
---------------------------

1. Hamiltonian center

2. `Hamiltonian saddle <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#one-degree-of-freedom-hyperbolic-equilibrium-point>`_

3. Duffing oscillator 

a. Time dependent Hamiltonian 

.. math::
   H(x,p_x,t) = \dfrac{1}{2}p_x^2 - \dfrac{\alpha}{2}x^2 + \dfrac{\beta}{4}x^4 - f(t) x \label{eqn:ham_duff_forced}


The corresponding vector field is

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = f_1(x,p_x) = p_x \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = f_2(x,p_x) = \alpha x - \beta x^3 + f(t)
   \end{align}

where :math:`f(t)` is the time dependent forcing function. The package offers two different options for the external forcing, :math:`f(t) = A\mathrm{sech}(t)\sin(\omega t)` and also :math:`f(t) = A\sin(\omega t)`. Other versions can be added manually by the user in the forcing function of the vector_fields.py file.

b. Unforced 

This model system corresponds to the undamped Duffing oscillator with model parameters in Hamiltonian~\eqref{eqn:eqn:ham_duff} chosen as :math:`f(t) = 0`.

c. Inverted 

This model system corresponds to the undamped unforced Duffing equation with Hamiltonian ~\eqref{eqn:eqn:ham_duff}, where the model parameters are chosen as :math:`\alpha = \beta = - 1`.

4. Saddle-node Hamiltonian 

This system is defined by the Hamiltonian:
.. math::
    H(x,p_x) = \dfrac{1}{2}p_x^2 + \dfrac{1}{2}x^2 + \dfrac{1}{3}x^3 \label{eqn:ham_saddnode}

and its associated vector field is:
.. math::
    \begin{align}
    \dot{x} = \dfrac{\partial H}{\partial p_x} = f_1(x,p_x) = p_x \\
    \dot{p}_x = -\dfrac{\partial H}{\partial x} = f_2(x,p_x) = -x - x^2
    \end{align} 

5. Non-autonomous double-gyre flow


Four dimensional phase space
----------------------------

1. `Saddle-center Hamiltonian <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#two-degrees-of-freedom-and-the-hyperbolic-periodic-orbit>`_ 

2. Hénon-Heiles Hamiltonian

Six dimensional phase space
---------------------------

1. `Saddle-center-center Hamiltonian <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#three-and-more-degrees-of-freedom-and-nhims>`_



.. automodule:: vector_fields
   :members:


