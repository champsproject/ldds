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
    Returns ..

    Parameters
    ----------
    LD : ndarray, shape(n, )
        array of computed LD values for array of initial conditions
    
    LD_type : str
        type of LD to plot. Options: 'forward', 'backward', 'total'
    
    grid_parameters : list of 3-tuples of floats
        input parameters of limits and size of mesh per axis
    
    tau : float
        Upper limit of integration.
        
    p_value : float
        Exponent in Lagrangian descriptor definition.
    
    norm : bool, optional
        normalise LD values by maximum
    
    colormap_name : str, optional
        valid name of matplotlib color-map for plot
    
    Returns
    -------
        scatter plot of LD function in phase-space 
    """
    x_min, x_max, Nx = grid_parameters[0]
    y_min, y_max, Ny = grid_parameters[1]
        
    if norm:
        LD = LD / LD.max()  # Scale LD output
    
    # Plot LDs
    fig,ax = plt.subplots(1, 1, dpi = 100)
    
    points_x = np.linspace(x_min, x_max, Nx)
    points_y = np.linspace(y_min, y_max, Ny)    
    
    contour = plt.contourf(points_x,points_y,LD,cmap=colormap_name,levels=200)
    
    # Customise appearance
    if p_value == 2:
        str_method = 'arclength - '
    elif p_value >= 1:
        str_method = r'p-norm $(p={})$'.format(p_value)
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
    
    ax.set_title(string_title, fontsize=12)
    ax.set_xlabel('$x$', fontsize=18)
    ax.set_ylabel('$y$', fontsize=18)
    
    fig.colorbar(contour)

    plt.show()
    
__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
