========================
Glossary of Common Terms
========================

This glossary gives the mathematical definition of dynamical systems concepts used in the software for computing Lagrangian descriptors. More detailed description can be found in books listed in the references; Ref. [wiggins2003]_. 


.. glossary::

Dynamical system
----------------

   A system that changes in time; usually described by differential equations (continuous time) or difference equations (sometimes called *maps*) (discrete time), or, possibly, some combination of the two.


Vector field, map, phase space, trajectory
------------------------------------------

   The equation of the form

   .. math::
      \begin{equation}
      \dot{x} = f(x,t; \mu)
      \end{equation}

   with :math:`x \in U \subset \mathbb{R}^n, t \in \mathbb{R}^1`, and :math:`\mu \in V \subset \mathbb{R}^p` wheere :math:`U, V` are open sets in :math:`\mathbb{R}^n, \mathbb{R}^p`, respectively, as a **vector field** or an **ordinary differential equation**. The above form is referred to as a non-autonomous or time-dependent ODE.

   The equation of the form

   .. math::
      \begin{equation}
      x \mapsto g(x; \mu)
      \end{equation}

   as a **map or difference equation**. 

   We refer to the space of dependent variables as the **phase space**. Both of these form of equations are referred to as dynamical system. Solution of the ordinary differential equations, :math:`x(t,t0,x0)` will be referred to as the **trajectory or phase curve** through the point :math:`x_0` at :math:`t = t_0`.

   The goal of dynamical systems analysis is to understand the geometry of solution curves in phase space. 



Equilibrium solutions, linear stability
---------------------------------------

   We give the definition of equilibrium solutions and their stability for an autonomous vector field of the form

   .. math::
      \begin{equation}
      \dot{x} = f(x; \mu), \qquad x \in \mathbb{R}^n
      \end{equation}

   The **equilibrium solution** is a point :math:`\bar{x} \in \mathbb{R}^n` such that 
   
   .. math::
      \begin{equation}
      f(\bar{x}) = 0.
      \end{equation}

   This constant solution of a vector field, that is a point in phase space where the vector field is zero, is also referred to as a "fixed point", "steady state", "stationary point", "rest point".

   Cautionary note: In case of non-autonomous vector fields, this definition of equilibrium solutions do not apply to the "frozen" (fixing :math:`t = \bar{t}`) vector field. Further explanation based on an example is given in the Chapter 1 of Ref. [wiggins2003]_ 

   Linear stability of an equilibrium solution (or point in phase space) of a vector field is given by the sign of the real part of the eigenvalues of the Jacobian of the vector field. More details on different type of stability is given in Ref. [wiggins2003]_ and [wiggins2017]_.




Periodic orbits and invariant manifolds
---------------------------------------






==========
References
==========
   
.. [wiggins2003] Wiggins, S., (2003) Introduction to applied nonlinear dynamical systems and chaos (Vol. 2). Springer Science & Business Media.

.. [wiggins2017] Wiggins, S. (2017) Ordinary Differential Equations (Ver. 2), https://figshare.com/articles/Ordinary_Differential_Equations/5311612,  Figshare.



