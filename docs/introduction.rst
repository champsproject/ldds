============
Introduction
============


The purpose of this guide is to give a brief overview of the theory, illustrate the computation and visualization of Lagrangian descriptors for dynamical systems, and automatically generated API documentation. We recommend that the user is familiar with nonlinear dynamical systems concepts mentioned in the :ref:`glossary`. Then, reading the brief introduction to the method of Lagrangian descriptor below and benchmark systems in :ref:`examples`.

Background and motivation
=========================

Understanding the geometry of the phase space of a dynamical system is a fundamental step in developing a complete picture of the dynamics. A trajectory diagnostic technique is a quick and practical approach to discovering the geometry of the structures that characterize phase space transport. This approach has become useful in analyzing a multitude of systems across geophysical fluid dynamics and chemical reactions. Lagrangian Descriptors for Dynamical Systems, or LDDS, is an open-source software written in Python represents a contribution in this direction using the method of Lagrangian descriptors. The basic idea behind this methodology is applicable to continuous or discrete dynamical system with or without periodic, aperiodic, stochastic forcing and with or without dissipation. The method encodes the geometry of phase space structures in the initial conditions on a two dimensional section by calculating a geometric property of the trajectories obtained from the initial conditions. 

Lagrangian descriptor is a scalar functional defined on a grid of initial conditions that assigns to every point a real number obtained by accumulating the values of a non-negative function along the trajectory starting from that initial condition in forward and backward time. There are many definitions of this non-negative function for computing Lagrangian descriptors (LDs) in the literature. Examples of such functions include the arclength of the trajectory, the p-norm of the components of the vector field defining the dynamical system under study, the Maupertuis classical action, etc. There is no definition of LDs that is better than the rest, and the choice depends on the dynamical characteristics of the system addressed. In this sense, there is always a trial and error stage where one would assess the different definitions, two-dimensional sections, integration time length for a given dynamical system. The LDDS package provides the user with an interface to perform such computations in a rapid prototyping manner. 

Lagrangian descriptors
======================

The Lagrangian descriptor (LD) as presented in Refs. [madrid2009]_, [mancho2013]_ is an arc-length of a trajectory calculated on a chosen initial time :math:`t_0` and measured for fixed forward and backward integration time, :math:`\tau`. For continuous time dynamical systems, Ref. [lopesino2017]_ gives an alternative definition of the LD which is useful for proving rigorous results and can be computed along with the trajectory. It provides a characterization of the notion of singular features of the LD that facilitates a proof for detecting invariant manifolds in certain model situations.  In addition, the "additive nature" of this new definition of LD provides 
an approach for assessing the influence of each degree-of-freedom separately on the Lagrangian descriptor.  This property was used in Ref. [demian2017]_ which showed that a Lagrangian descriptor can be used to detect Lyapunov periodic orbits in the two degrees-of-freedom Hénon-Heiles Hamiltonian system. We will describe this procedure for two and three degrees-of-freedom linear autonomous Hamiltonian systems. We begin by establishing notation in the general setting of a time-dependent vector field where 

.. math::
    \begin{equation}
    \frac{d\mathbf{x}}{dt} = \mathbf{v}(\mathbf{x},t), \quad \mathbf{x} \in \mathbb{R}^n \;,\; t \in \mathbb{R}
    \end{equation}

where :math:`\mathbf{v}(\mathbf{x},t) \in C^r (r \geq 1)` in :math:`\mathbf{x}` and continuous in time. The definition of LDs depends on the initial condition :math:`\mathbf{x}_{0} = \mathbf{x}(t_0)`, on the initial time :math:`t_0` (trivial for autonomous systems) and the integration time :math:`\tau`, and the type of norm of the trajectory's components, and takes the form,


.. math::
    \begin{equation}
    M_p(\mathbf{x}_{0},t_0,\tau) = \displaystyle{\int^{t_0+\tau}_{t_0-\tau} \sum_{i=1}^{n} |\dot{x}_{i}(t;\mathbf{x}_{0})|^p \; dt} \label{eqn:M_function}
    \end{equation}

where :math:`p \in (0,1]` and :math:`\tau \in \mathbb{R}^{+}` are freely chosen parameters,  and the overdot symbol represents the derivative with respect to time. It is to be noted here that there are three formulations of the function :math:`M_p` in the literature: the arc length of a trajectory in phase space [madrid2009]_, the arc length of a trajectory projected on the configuration space [junginger2016lagrangian]_, [junginger2016transition]_, [junginger2016uncovering]_, [junginger2017chemical]_ and the sum of the :math:`p`-norm of the vector field components [lopesino2015]_, [lopesino2017]_.
Although the latter formulation of the Lagrangian descriptor developed in Refs. [lopesino2015]_, [lopesino2017]_ does not resemble the arc length, the numerical results using either of these forms have been shown to be in agreement and promise of predictive capability in geophysical flows ([delacamara2012]_, [garciagarrido2015]_, [ramos2018]_, [mendoza2014lagrangian]_). The formulation we adopt here is motivated by the fact that this allows for proving rigorous result, which we will discuss in the next section, connecting the singular features and minimum in the LD plots with NHIM and its stable and unstable manifolds. 
It follows from the result that 

.. math:: 
    \begin{align}
    \mathcal{W}^s(\mathbf{x}_0, t_0) & = \text{argmin} \; \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) \\
    \mathcal{W}^u(\mathbf{x}_0, t_0) & = \text{argmin} \; \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau)
    \end{align}

where the stable and unstable manifolds (:math:`\mathcal{W}^s(\mathbf{x}_0, t_0)` and :math:`\mathcal{W}^u(\mathbf{x}_0, t_0)`) denote the invariant manifolds at intial time :math:`t_0` and :math:`\text{argmin} (\cdot)` denotes the argument that minimizes the function :math:`\mathcal{L}^{(\cdot)}(\mathbf{x}_0, t_0, \tau)` in forward and backward time, respectively. In addition, the coordinates of the NHIM at time :math:`t_0` is given by the intersection :math:`\mathcal{W}^s(\mathbf{x}_0, t_0)` and :math:`\mathcal{W}^u(\mathbf{x}_0, t_0)` of the stable and unstable manifolds, and thus given by

.. math::
    \begin{align}
    \mathcal{M}(\mathbf{x}_0, t_0) & = \text{argmin} \; \left( \mathcal{L}^{(f)}(\mathbf{x}_0, t_0, \tau) + \mathcal{L}^{(b)}(\mathbf{x}_0, t_0, \tau) \right) = \text{argmin} \; \mathcal{L}(\mathbf{x}_0, t_0, \tau) \qquad \text{NHIM}
    \end{align}


References
==========
   
   
.. [madrid2009] Madrid, J. A. J. and Mancho, A. M. (2009) Distinguished trajectories in time dependent vector fields. Chaos, 19, 013111.

.. [demian2017] Demian, A. S., and Wiggins, S. (2017). Detection of periodic orbits in Hamiltonian systems using Lagrangian descriptors. International Journal of Bifurcation and Chaos, 27(14), 1750225.

.. [lopesino2017] Lopesino, C., Balibrea-Iniesta, F., García-Garrido, V. J., Wiggins, S., and Mancho, A. M. (2017). A theoretical framework for Lagrangian descriptors. International Journal of Bifurcation and Chaos, 27(01), 1730001.

.. [lopesino2015] Lopesino, C., Balibrea, F., Wiggins, S., and Mancho, A. M. (2015). Lagrangian descriptors for two dimensional, area preserving, autonomous and nonautonomous maps. Communications in Nonlinear Science and Numerical Simulation, 27(1-3), 40–51.

.. [junginger2016lagrangian] Junginger, A. and Hernandez, R. (2016a). Lagrangian descriptors in dissipative systems. Physical Chemistry Chemical Physics, 18(44), 30282–30287. 

.. [junginger2016transition] Junginger, A., Craven, G. T., Bartsch, T., Revuelta, F., Borondo, F., Benito, R., and Hernandez, R. (2016). Transition state geometry of driven chemical reactions on time-dependent double- well potentials. Physical Chemistry Chemical Physics, 18(44), 30270–30281.

.. [junginger2016uncovering] Junginger, A. and Hernandez, R. (2016b). Uncovering the Geometry of Barrierless Reactions Using Lagrangian Descriptors. The Journal of Physical Chemistry B, 120(8), 1720–1725.


.. [junginger2017chemical] Junginger, A., Duvenbeck, L., Feldmaier, M., Main, J., Wunner, G., and Hernandez, R. (2017a). Chemical dynamics between wells across a time-dependent barrier: Self-similarity in the Lagrangian descriptor and reactive basins. The Journal of chemical physics, 147(6), 064101.


.. [mancho2013] Mancho, A. M., Wiggins, S., Curbelo, J., and Mendoza, C. (2013). Lagrangian Descriptors: A Method for Revealing Phase Space Structures of General Time Dependent Dynamical Systems. Communications in Nonlinear Science and Numerical, 18, 3530–3557.


.. [mendoza2014lagrangian] Mendoza, C., Mancho, A. M., and Wiggins, S. (2014). Lagrangian descriptors and the assessment of the predictive capacity of oceanic data sets. Nonlinear Processes in Geophysics, 21(3), 677–689.

.. [garciagarrido2015] García-Garrido, V. J., Mancho, A. M., and Wiggins, S. (2015). A dynamical systems approach to the surface search for debris associated with the disappearance of flight MH370. Nonlin. Proc. Geophys., 22, 701–712.

.. [ramos2018] Ramos, A. G., García-Garrido, V. J., Mancho, A. M., Wiggins, S., Coca, J., Glenn, S., Schofield, O., Kohut, J., Aragon, D., Kerfoot, J., Haskins, T., Miles, T., Haldeman, C., Strandskov, N., All- sup, B., Jones, C., and Shapiro., J. (2018). Lagrangian coherent structure assisted path planning for transoceanic autonomous underwater vehicle missions. Scientfic Reports, 4, 4575.

.. [delacamara2012] de la Cámara, A., Mancho, A. M., Ide, K., Serrano, E., and Mechoso, C. (2012). Routes of transport across the Antarctic polar vortex in the southern spring. J. Atmos. Sci., 69(2), 753–767.



