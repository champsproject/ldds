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
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.5,3), dpi=130)
    
    points_x = np.linspace(x_min, x_max, Nx)
    points_y = np.linspace(y_min, y_max, Ny)    
    
    con0 = ax0.contourf(points_x,points_y,LD,cmap=colormap_name,levels=200)
    
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
    
    fig.colorbar(con0, ax=ax0, ticks=np.linspace(np.nanmin(LD),np.nanmax(LD),11),format='%.2f')
    
    gx,gy = np.gradient(LD, 0.05, 0.05)
    scalar = np.sqrt(gx**2 + gy**2)
    Z = scalar/scalar.max()
    
    con1 = ax1.contourf(points_x,points_y,Z,cmap='Reds',levels=200)
    ax1.set_title('LD gradient magnitude')
    ax1.set_xlabel('$x$')
    ax1.label_outer()
    fig.colorbar(con1, ax=ax1, ticks=np.linspace(np.nanmin(Z),np.nanmax(Z),11),format='%.2f')
    
    plt.show()
    
__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido, Vladimir Krajnak'
__status__ = 'Development'
