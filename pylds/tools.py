"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets

def ld_plot(LD, LD_gradient, grid_parameters, colormap, interactive, string_title, colormap_gradient):
    """
    Draws a Lagrangian descriptor contour plot and a contour plot showing the magnitude of its gradient field.

    Parameters
    ----------
    LD : ndarray, shape(n, )
        Array of Lagrangian Descriptor values.

    LD_gradient : ndarray, shape(n, )
        Array of Lagrangian Descriptor gradient values.

    grid_parameters : list of 3-tuples of floats
        Limits and size of mesh per axis.

    colormap : str, optional
        Name of matplotlib colormap for contouor plot.

    interactive : bool
        True allows interactively adjusting the gradient contouor plot minimum and maximum.

    string_title : string
        Plot title.

    colormap_gradient : string
        Name of matplotlib colormap for gradient contouor plot.

    Returns
    -------
        Nothing.
    """

    if type(grid_parameters) == dict:
        #n-DoF systems
        slice_parameters = grid_parameters['slice_parameters'] # 2n-D grid
        dims_slice = np.array(grid_parameters['dims_slice'])
        slice_axes_labels = np.array(['$x$','$y$','$p_x$','$p_y$'])
        slice_axes_labels = slice_axes_labels[dims_slice==1]
    else:
        #1-DoF systems
        slice_parameters = grid_parameters # 2-D grid
        slice_axes_labels = ['$x$', '$p_x$']

    ax1_min, ax1_max, N1 = slice_parameters[0]
    ax2_min, ax2_max, N2 = slice_parameters[1]

    points_ax1 = np.linspace(ax1_min, ax1_max, N1)
    points_ax2 = np.linspace(ax2_min, ax2_max, N2)

    LD = LD - np.nanmin(LD)
    LD = LD / np.nanmax(LD)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.5,3), dpi=130, sharex=True, sharey=True)
    fig.suptitle(string_title, fontsize=14, y=1.00)
    plt.subplots_adjust(bottom=0.13, top=0.85)

    # LD plot
    con0 = ax0.contourf(points_ax1, points_ax2, LD, cmap=colormap, vmin=0, vmax=1, levels=100)
    ax0.set_title('LD values')
    ax0.set_xlabel(slice_axes_labels[0])
    ax0.set_ylabel(slice_axes_labels[1])
    ticks_LD = np.linspace(0, 1, 11)
    fig.colorbar(con0, ax=ax0, ticks=ticks_LD, format='%.1f')

    #gradient plot
    vmin = np.nanmin(LD_gradient)
    vmax = np.nanmax(LD_gradient)
    con1 = ax1.contourf(points_ax1, points_ax2, LD_gradient, cmap=colormap_gradient, levels=100)
    if interactive:
        @widgets.interact(grad_min=(vmin, vmax-0.01, .01),grad_max=(vmin+0.01, vmax, .01))
        def update(grad_min=vmin,grad_max=vmax):
            grad_max = max(grad_min+0.01, grad_max)
            con1.set_clim(grad_min,grad_max)

    ticks_gradient = np.linspace(vmin,vmax, 11)
    fig.colorbar(con1, ax=ax1, ticks=ticks_gradient, format='%.1f')
    ax1.set_title('LD gradient magnitude')
    ax1.set_xlabel(slice_axes_labels[0])
    ax1.label_outer()

    plt.show()

def draw_all_lds(LD_forward, LD_backward, grid_parameters, tau, p_value, colormap='bone', interactive=False):
    """
    Draws the forward, backward and total Lagrangian descriptor contour plots and a contour plots showing the magnitude of its gradient field.

    Parameters
    ----------
    LD_forward : ndarray, shape(n, )
        Array of Lagrangian Descriptor values in forward time.

    LD_backward : ndarray, shape(n, )
        Array of Lagrangian Descriptor values in backward time.

    grid_parameters : list of 3-tuples of floats
        Limits and size of mesh per axis.

    tau : float
        Upper limit of integration.

    p_value : float
        Exponent in Lagrangian descriptor definition.

    norm : bool, optional
        True normalises LD values.

    colormap : str, optional
        Name of matplotlib colormap for plot.

    interactive : bool
        True allows interactively adjusting the gradient plot minimum and maximum.

    Returns
    -------
        Nothing.
    """

    # Prepare method name
    if p_value == 2:
        str_method = 'arclength LD'
    elif p_value >= 1:
        str_method = r'p-norm LD$(p={})$'.format(p_value)
    elif p_value == 0:
        str_method = 'action-based LD'
    elif p_value < 1:
        str_method = r'LD$_p$ $(p={})$'.format(p_value)
    t_final=abs(tau)

    # Prepare LD arrays
    def norm_and_grad(LD):
        gradient_x, gradient_y = np.gradient(LD)
        gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
        gradient_magnitude = gradient_magnitude - np.nanmin(gradient_magnitude)
        gradient_magnitude = gradient_magnitude / np.nanmax(gradient_magnitude)
        return gradient_magnitude

    # Plot LDs

    string_title = r'Forward {}, $\tau={}$'.format(str_method,t_final)
    LD_forward_gradient = norm_and_grad(LD_forward)
    ld_plot(LD_forward, LD_forward_gradient, grid_parameters, colormap, interactive, string_title, colormap_gradient='Reds')

    string_title = r'Backward {}, $\tau={}$'.format(str_method,t_final)
    LD_backward_gradient = -norm_and_grad(LD_backward)
    ld_plot(LD_backward, LD_backward_gradient, grid_parameters, colormap, interactive, string_title, colormap_gradient='Blues_r')

    string_title = r'Total {}, $\tau={}$'.format(str_method,t_final)
    ld_plot(LD_forward+LD_backward, LD_forward_gradient+LD_backward_gradient, grid_parameters, colormap, interactive, string_title, colormap_gradient='RdBu_r')

__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
