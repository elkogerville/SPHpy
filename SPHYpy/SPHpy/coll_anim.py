def SPH_coll_animation(position, density, N, M, start, end, **kwargs):
    """
    this function will create screenshots for the purpose of animating the SPH simulation. 
    PLEASE CREATE A FOLDER NAMED 'animation' IN THE SAME DIRECTORY AS THE SIMULATION CODE
    --------------------------------------------------------------------------------------
    position, density [numpy array]: N+Mx3 array of positions and densities
    N [integer]: number of particles in first star
    M [integer]: number of particles in second star
    start/end [integer]: start and end timesteps; must be within timestep range
    elev/azim [float]: sets camera angle
    width/height [float]: sets figure dimensions
    c1/c2 [matplotlib color]: sets color for first/second star
        ex: c1 = 'darkslateblue'
    s [float]: sets particle size
    lightmode [whatever u want]: if defined, sets background color to white; by default set to black
        ex: lightmode = 'yes please' or lightmode = True
    ------------------------------------------------------------------------------------
    OUTPUT: numbered screenshots corresponding to each timestep; use a video software to 
    string them into an animation
    """
    
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d

    # default fig and camera viewing angles
    width, height = 10, 10
    elev = 45 # z viewing angle (0 = edge on; 90 = Bird's eye view)
    azim = 90 # xy plane rotation angle
    c1 = 'mediumslateblue'
    c2 = 'darkorchid'
    s = 10
    
    if ('elev') in kwargs:
        elev = kwargs['elev']
    if ('azim') in kwargs:
        azim = kwargs['azim']
    if ('width') in kwargs:
        width = kwargs['width']
    if ('height') in kwargs:
        height = kwargs['height']
    if ('c1') in kwargs:
        c1 = kwargs['c1']
    if ('c2') in kwargs:
        c2 = kwargs['c2']
    if ('s') in kwargs:
        s = kwargs['s']
    
    # define figure
    fig = plt.figure(figsize=(width, height))
    ax = plt.axes(projection='3d')

    # formatting
    plt.rcParams['font.family'] = 'sans-serif' # set font
    # set background color to black
    ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
    ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
    ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
    
    if ('lightmode') in kwargs:
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    # set view angles
    ax.view_init(elev = elev, azim = azim)
    # set lables
    ax.set_xlabel('\u03A7')
    ax.set_ylabel('\u03A5')
    ax.set_zlabel('Z')
    
    Ntot = N + M
    
    for i in range(start, end):
        # clear previous plot
        ax.clear()
        # correct slicing index
        Slice = Ntot*i + Ntot
        SliceN = Slice - M

        # plot
        plotN = ax.scatter3D(position[Ntot*i:SliceN,0], position[Ntot*i:SliceN,1], position[Ntot*i:SliceN,2], s = s, c = c1) 
        plotM = ax.scatter3D(position[SliceN:Slice, 0], position[SliceN:Slice, 1], position[SliceN:Slice, 2], s = s, c = c2)
        
        plt.savefig('animation/animation{}.png'.format(i))
