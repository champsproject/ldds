"""
Module description ...

Reference:
- Surename1, Forename1 Initials., Surename2, Forename2 Initials, YEAR. Publication/Book title
Publisher, Number(Volume No), pp.142-161.
"""
import numpy as np
import matplotlib.pyplot as plt

def draw_lagrangian_descriptor(LD, LD_type, grid_parameters, tau, p_value, norm = True, colormap_name='bone'):
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
    
    Returns
    -------
        Nothing.
    """
    x_min, x_max, Nx = grid_parameters[0]
    y_min, y_max, Ny = grid_parameters[1]
        
    if norm:
        LD = LD - np.nanmin(LD)  # Scale LD output
        LD = LD / np.nanmax(LD)  # Scale LD output
    
    # Plot LDs
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.5,3), dpi=130)
    
    points_x = np.linspace(x_min, x_max, Nx)
    points_y = np.linspace(y_min, y_max, Ny)    
    
    con0 = ax0.contourf(points_x, points_y, LD, cmap=colormap_name, levels=200)
    
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
    
    fig.suptitle(string_title, fontsize=14, y=1.04)
    ax0.set_title('LD values')
    ax0.set_xlabel('$x$')
    ax0.set_ylabel('$y$')
    
    ticks_LD = np.linspace(np.nanmin(LD), np.nanmax(LD), 11)
    fig.colorbar(con0, ax=ax0, ticks=ticks_LD, format='%.2f')
    
    gradient_x, gradient_y = np.gradient( LD, 0.05, 0.05)
    gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
    gradient_magnitude = gradient_magnitude/gradient_magnitude.max()
    
    con1 = ax1.contourf(points_x, points_y, gradient_magnitude, cmap='Reds', levels=200)
    ax1.set_title('LD gradient magnitude')
    ax1.set_xlabel('$x$')
    ax1.label_outer()
    
    ticks_gradient = np.linspace(np.nanmin(gradient_magnitude), np.nanmax(gradient_magnitude), 11)
    fig.colorbar(con1, ax=ax1, ticks=ticks_gradient, format='%.2f')
    
    plt.show()
    
__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
