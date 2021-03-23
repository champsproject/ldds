Introduction
============

LDDS is a Python library for computing Lagrangian descriptors of dynamical systems.



Lagrangian descriptors
----------------------

The Lagrangian descriptor (LD) as presented in Ref.\cite{madrid2009} is the arc length of a trajectory calculated on a chosen initial time $t_0$ and measured for fixed forward and backward integration time, $\tau$. For continuous time dynamical systems, Ref.\cite{lopesino2017} gives an alternative definition of the LD which is useful for proving rigorous results and can be computed along with the trajectory. It provides a characterization of the notion of singular features of the LD that facilitates a proof for detecting invariant manifolds in certain model situations.  In addition, the ``additive nature'' of this new definition of LD provides 
an approach for assessing the influence of each degree-of-freedom separately on the Lagrangian descriptor.  This property was used in Ref.\cite{demian2017} which showed that a Lagrangian descriptor can be used to detect Lyapunov periodic orbits in the two degrees-of-freedom H{\'e}non-Heiles Hamiltonian system. We will describe this procedure for two and three degrees-of-freedom linear autonomous Hamiltonian systems. We begin by establishing notation in the general setting of a time-dependent vector field where 

.. math::
    \frac{d\mathbf{x}}{dt} = \mathbf{v}(\mathbf{x},t), \quad \mathbf{x} \in \mathbb{R}^n \;,\; t \in \mathbb{R}

where :math:`\mathbf{v}(\mathbf{x},t) \in C^r (r \geq 1)` in :math:`\mathbf{x}` and continuous in time. The definition of LDs depends on the initial condition :math:`\mathbf{x}_{0} = \mathbf{x}(t_0)`, on the initial time :math:`t_0` (trivial for autonomous systems) and the integration time :math:`\tau`, and the type of norm of the trajectory's components, and takes the form,


.. math::
    M_p(\mathbf{x}_{0},t_0,\tau) = \displaystyle{\int^{t_0+\tau}_{t_0-\tau} \sum_{i=1}^{n} |\dot{x}_{i}(t;\mathbf{x}_{0})|^p \; dt} \label{eqn:M_function}

where :math:`p \in (0,1]` and :math:`\tau \in \mathbb{R}^{+}` are freely chosen parameters,  and the overdot symbol represents the derivative with respect to time. It is to be noted here that there are three formulations of the function :math:`M_p` in the literature: the arc length of a trajectory in phase space~\cite{madrid2009}, the arc length of a trajectory projected on the configuration space~~\cite{junginger2016transition,junginger2016uncovering,junginger2017chemical,junginger2017variational}, and the sum of the :math:`p`-norm of the vector field components~\cite{lopesino_2015,lopesino2017}.
Although the latter formulation of the Lagrangian descriptor~\eqref{eqn:M_function} developed in Ref.~\cite{lopesino_2015,lopesino2017} does not resemble the arc length, the numerical results using either of these forms have been shown to be in agreement and promise of predictive capability in geophysical flows~\cite{amism11,mendoza2014,ggmwm15,ramos2018}. The formulation we adopt here is motivated by the fact that this allows for proving rigorous result, which we will discuss in the next section, connecting the singular features and minimum in the LD plots with NHIM and its stable and unstable manifolds. 
It follows from the result that 

.. math:: 
    \begin{align}
    \mathcal{W}^s(\mathbf{x}_0, t_0) & = \text{\rm argmin} \; \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) \\
    \mathcal{W}^u(\mathbf{x}_0, t_0) & = \text{\rm argmin} \; \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau)
    \end{align}

where the stable and unstable manifolds (:math:`\mathcal{W}^s(\mathbf{x}_0, t_0)` and :math:`\mathcal{W}^u(\mathbf{x}_0, t_0)`) denote the invariant manifolds at intial time :math:`t_0` and :math:`\text{\rm argmin} \; (\cdot)` denotes the argument that minimizes the function :math:`\mathcal{L}^{(\cdot)}(\mathbf{x}_0, t_0, \tau)` in forward and backward time, respectively. In addition, the coordinates of the NHIM at time :math:`t_0` is given by the intersection :math:`\mathcal{W}^s(\mathbf{x}_0, t_0)` and :math:`\mathcal{W}^u(\mathbf{x}_0, t_0)` of the stable and unstable manifolds, and thus given by

.. math::
    \begin{align}
    \mathcal{M}(\mathbf{x}_0, t_0) & = \text{\rm argmin} \; \left( \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) + \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau) \right) = \text{\rm argmin} \; \mathcal{L}(\mathbf{x}_0, t_0, \tau) \qquad \text{NHIM}
    \end{align}









