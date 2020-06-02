import numpy as np
import matplotlib.pyplot as plt

def draw_lagrangian_descriptor(LD, LD_type, GRID_PARAMETERS, LD_PARAMETERS, norm = True):
    """
    Returns ..

    Parameters
    ----------
    LD : ndarray, shape(n, )
        array of computed LD values for array of initial conditions
    
    LD_type : str
        Type of LD to plot. Options: 'forward', 'backward', 'total'
        
    p_norm = float
        p-value of Lp-norm used to calculated input LD
    
    Returns
    -------
        scatter plot of LD function in phase-space 
    """
    x_min, x_max, Nx = GRID_PARAMETERS[0]
    y_min, y_max, Ny = GRID_PARAMETERS[1]
    
    t_initial, t_final, dt = LD_PARAMETERS[0]
    p_norm = LD_PARAMETERS[1]
    
    tau = t_final - t_initial
    
    LD = LD.reshape(Nx, Ny).T # Reshape 1D array
    if norm:
        LD = LD / LD.max()  # Scale LD output
    ################################### 
    # Plot LDs
    resolution = 100 # in dpi
    fig,ax = plt.subplots(1, 1, dpi = resolution)
    
    points_x = np.linspace(x_min, x_max, Nx)
    points_y = np.linspace(y_min, y_max, Ny)    
    X, Y = np.meshgrid(points_x, points_y)  # Grid in phase-space
    
    colormap_name = 'bone'
    scatter = plt.scatter(X, Y, c = LD, cmap = colormap_name)
    ###################################
    # Customise appearance
    if p_norm != 2:
        str_meth = ' '.join(['p-norm (p=',str(p_norm),') - '])
    else:
        str_meth = 'arclength - '
    
    if LD_type == 'forward':
        string_title = ['Forward LD ', str_meth, '(','$\\tau=$ ',str(tau),' , ','$t_0=$',str(t_initial),')']
    elif LD_type == 'backward':
        string_title = ['Backward LD ', str_meth, '(','$\\tau=$ ',str(tau),' , ','$t_0=$',str(t_initial),')']
    elif LD_type == 'total':
        string_title = ['Total LD ', str_meth, '(','$\\tau=$ ',str(tau),' , ','$t_0=$',str(t_initial),')']
    else: 
        string_title = ['']
        print('Incorrect "LD_type". Valid options: forward, backward, total. Plot will appear without title')
    string_title = ' '.join(string_title)
    
    ax.set_title(string_title, fontsize=12)
    ax.set_xlabel('$x$', fontsize=18)
    ax.set_ylabel('$y$', fontsize=18)
    ax.set_aspect('auto')

    fig.colorbar(scatter) # Add color bar

    plt.show()
    
__author__ = 'Broncio Aguilar-Sanjuan, Victor-Jose Garcia-Garrido'
__status__ = 'Development'
