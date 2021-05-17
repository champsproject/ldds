========================
Glossary of Common Terms
========================

This glossary gives the mathematical definition of dynamical systems concepts used in the software for computing Lagrangian descriptors. More detailed description can be found in books listed in the references; Ref. [wiggins2003]_. 


.. glossary::

Vector field, map, phase space, trajectory
------------------------------------------

      The equation of the form

      .. math::
         \begin{equation}
         \dot{x} = f(x,t; \mu)
         \end{equation}

      with :math:`x \in U \subset \mathbb{R}^n, t \in \mathbb{R}^1`, and :math:`\mu \in V \subset \mathbb{R}^p` wheere :math:`U, V` are open sets in :math:`\mathbb{R}^n, \mathbb{R}^p`, respectively, as a **vector field or ordinary differential equation**.

      The equation of the form

      .. math::
         \begin{equation}
         x \mapsto g(x; \mu)
         \end{equation}

      as a **map or difference equation**. 

      We refer to the space of dependent variables as the **phase space**. Both of these form of equations are referred to as dynamical system. Solution of the ordinary differential equations, :math:`x(t,t0,x0)` will be referred to as the **trajectory or phase curve** through the point :math:`x_0` at :math:`t = t_0`.

      The goal of dynamical systems analysis is to understand the geometry of solution curves in phase space. 



Equilibrium points, linear stability
------------------------------------



Periodic orbits and invariant manifolds
---------------------------------------






==========
References
==========
   
.. [wiggins2003] Wiggins, S., (2003) Introduction to applied nonlinear dynamical systems and chaos (Vol. 2). Springer Science & Business Media.




