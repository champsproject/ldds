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
   
where :math:`a,b` are the model parameters.

The inverse Hénon map is:

.. math::
   \begin{align}
    x_{n} = y_{n+1}
    y_{n} = \dfrac{x_{n+1} - a + y_{n+1}^2}{b}
   \end{align}


Two dimensional phase space
---------------------------

1. Hamiltonian center

The Hamiltonian function:

.. math::
   H(x,p_x) = \dfrac{\omega}{2} \left( p_x^2 + x^2 \right) \label{eqn:ham_center1dof}

defines the normal form of a 1 DoF system with a center equilibrium at the origin. The associated dynamical system is:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = f_1(x,p_x) = \omega p_x \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = f_2(x,p_x) = -\omega x
   \end{align}
   

2. `Hamiltonian saddle <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#one-degree-of-freedom-hyperbolic-equilibrium-point>`_

The Hamiltonian function:

.. math::
   H(x,p_x) = \dfrac{\lambda}{2} \left( p_x^2 - x^2 \right) \label{eqn:ham_saddle1dof}

defines the normal form of a 1 DoF system with a saddle equilibrium at the origin. The associated dynamcical system is:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = f_1(x,p_x) = \lambda p_x \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = f_2(x,p_x) = \lambda x
   \end{align}

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

The Saddle-center Hamiltonian function is:

.. math::
   H(x,y,p_x,p_y) = \dfrac{1}{2} \left( p_x^2 + p_y^2 + y^2 - x^2) \right)  \label{eqn:ham_saddcen}

and Hamilton's equations of motion are in this case:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = f_1(x,y,p_x,p_y) = p_x \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y} = f_2(x,y,p_x,p_y) = p_y \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = f_3(x,y,p_x,p_y) = x \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y} = f_4(x,y,p_x,p_y) = - y
   \end{align}

2. Hénon-Heiles Hamiltonian

The Hamiltonian function for the Hénon-Heiles system is given:

.. math::
   H(x,y,p_x,p_y) = \dfrac{1}{2} \left( p_x^2 + p_y^2 \right) + \dfrac{1}{2} \left( x^2 + y^2 \right) + yx^2 - \dfrac{1}{3} y^3 \label{eqn:ham_henheil}

and Hamilton's equations of motion are in this case:

.. math::
   \begin{align}
   \dot{x} &= \dfrac{\partial H}{\partial p_x} = f_1(x,y,p_x,p_y) = p_x \\
   \dot{y} &= \dfrac{\partial H}{\partial p_y} = f_2(x,y,p_x,p_y) = p_y \\
   \dot{p}_x &= -\dfrac{\partial H}{\partial x} = f_3(x,y,p_x,p_y) = - x - 2xy \\
   \dot{p}_y &= -\dfrac{\partial H}{\partial y} = f_4(x,y,p_x,p_y) = - x^2 - y + y^2
   \end{align}

Six dimensional phase space
---------------------------

1. `Saddle-center-center Hamiltonian <https://champsproject.github.io/lagrangian_descriptors/content/chapter2_1.html#three-and-more-degrees-of-freedom-and-nhims>`_



.. automodule:: vector_fields
   :members:


