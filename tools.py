"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets

def draw_lagrangian_descriptor(LD, LD_type, grid_parameters, tau, p_value, colormap_name='bone', colormap_mode=1, interactive=False):
    """
    Draws a Lagrangian descriptor contour plot and a contour plot showing the magnitude of its gradient field.

    Parameters
    ----------
    LD : ndarray, shape(n, )
        Array of Lagrangian Descriptor values.
    
    LD_type : str
        Type of LD to plot. Options: 'forward', 'backward', 'total'.
    
    grid_parameters : list of 3-tuples of floats
        Limits and size of mesh per axis.
    
    tau : float
        Upper limit of integration.
        
    p_value : float
        Exponent in Lagrangian descriptor definition.
    
    norm : bool, optional
        True normalises LD values.
    
    colormap_name : str, optional
        Name of matplotlib colormap for plot.
        
    interactive : bool
        True allows interactively adjucting the gradient plot minimum.
    
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
        
    # normalise output
    LD = LD - np.nanmin(LD)
    LD = LD / np.nanmax(LD)

    # Plot LDs
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.5,3), dpi=130)
    
    points_ax1 = np.linspace(ax1_min, ax1_max, N1)
    points_ax2 = np.linspace(ax2_min, ax2_max, N2)
    
    if colormap_mode == 1:
        vmin, vmax = 0, 1
    elif colormap_mode == 2:
        vmin = np.nanmean(LD)-np.nanstd(LD)
        vmax = 1
    
    con0 = ax0.contourf(points_ax1, points_ax2, LD, cmap=colormap_name, vmin=vmin, vmax=vmax, levels=200)

    # Customise appearance
    if p_value == 2:
        str_method = 'arclength - '
    elif p_value >= 1:
        str_method = r'p-norm $(p={})$'.format(p_value)
    elif p_value == 0:        
        str_method = 'action-based'
    elif p_value < 1:
        str_method = r'LD$_p$ $(p={})$'.format(p_value)
    
    t_final=abs(tau)
    if LD_type == 'forward':
        string_title = r'Forward LD {}, $\tau={}$'.format(str_method,t_final)
    elif LD_type == 'backward':
        string_title = r'Backward LD {}, $\tau={}$'.format(str_method,t_final)
    elif LD_type == 'total':
        string_title = r'Total LD {}, $\tau={}$'.format(str_method,t_final)
    else: 
        string_title = ''
        print('Incorrect "LD_type". Valid options: forward, backward, total. Plot will appear without title')
    
    fig.suptitle(string_title, fontsize=14, y=1.00)
    plt.subplots_adjust(bottom=0.13, top=0.85)
    ax0.set_title('LD values')
    ax0.set_xlabel(slice_axes_labels[0])
    ax0.set_ylabel(slice_axes_labels[1])
    
    ticks_LD = np.linspace(0, 1, 11)
    fig.colorbar(con0, ax=ax0, ticks=ticks_LD, format='%.1f')
    
    gradient_x, gradient_y = np.gradient( LD, 0.05, 0.05)
    gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
    gradient_magnitude = gradient_magnitude - np.nanmin(gradient_magnitude)
    gradient_magnitude = gradient_magnitude / np.nanmax(gradient_magnitude)
    
    
    if interactive:
        @widgets.interact(grad_minimum=(0, 1, .1))
        def update(grad_minimum=0.0):
            ax1.clear()    
            con1 = ax1.contourf(points_ax1, points_ax2, gradient_magnitude, cmap='Reds', levels=np.linspace(grad_minimum,1,100))
            ax1.label_outer()
    con1 = ax1.contourf(points_ax1, points_ax2, gradient_magnitude, cmap='Reds', levels=np.linspace(0,1,100))
    ticks_gradient = np.linspace(0,1, 11)
    fig.colorbar(con1, ax=ax1, ticks=ticks_gradient, format='%.1f')
    ax1.set_title('LD gradient magnitude')
    ax1.set_xlabel(slice_axes_labels[0])
    ax1.label_outer()
    
    plt.show()
    
__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
