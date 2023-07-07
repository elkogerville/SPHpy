def positions(N, mode = 'randn'):
    """
    this function creates initial position and velocity arrays
    by default uses random.randn but if mode = 'normal', will create positions from normal distribution
    centered around user input
    """
    import numpy as np
    if mode == 'randn':
        # Nx3 position and velocity arrays, starting at rest
        position = np.random.randn(N, 3)
        velocity = np.zeros((N, 3))
        return position, velocity, N
    if mode == 'normal':
        mean = input('mean: ')
        scale = input('scale: ')
        # Nx3 position and velocity arrays, starting at rest
        x = np.random.normal(loc = 0, scale = 5, size = N)
        y = np.random.normal(loc = 0, scale = 5, size = N)
        z = np.random.normal(loc = 0, scale = 5, size = N)
        # ensure correct shape
        position = np.array((x, y, z)).T
        velocity = np.zeros((N, 3))
        return position, velocity, N
    else:
        print('ERROR: please specify one of two modes: \n randn, normal')