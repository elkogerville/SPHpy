def main_SPH(position, velocity, N, M, R, k, n, nu, timesteps, dt):
    """
    smoothed particle hydrodynamical simulator for stellar structure based on Philip Mocz derivations
    source: https://pmocz.github.io/manuscripts/pmocz_sph.pdf
    --------------------------------------------------------------------------------------------------
    position [numpy array]: Nx3 array of x, y, z particle positions
    velocity [numpy array]: Nx3 array of x, y, z particle velocities
    N [integer]: number of particles in simulation
    M [float]: total mass of star
    R [float]: radius of star
    h [float]: smoothing length
    k [float]: state equation constant
    n [float]: polytropic index
    nu [float]: damping coefficient to achieve equilibrium state
    --------------------------------------------------------------
    OUTPUT [numpu array]: Mx3 position array and Mx3 density array
    """
    
    def smoothing_kernel(X, Y, Z, h):
        """
        evaluates the smoothing kernel function
        ---------------------------------------
        X, Y, Z [float]: x, y, z positions
        h [float]: smoothing length
        --------------------------------------------
        OUTPUT [float]: evaluated smoothing function 
        """
        import numpy as np
        r_sq = (X**2 + Y**2 + Z**2)
        r = np.sqrt(r_sq)
        pi_sq = np.sqrt(np.pi)

        return (1 / (h*pi_sq)**3) * np.exp(-r/h**2)

    def gradient_smooth(X, Y, Z, h):
        """
        evaluates the gradient of smoothing function
        --------------------------------------------
        X, Y, Z [float]: x, y, z positions
        h [float]: smoothing length
        ----------------------------------
        OUTPUT [float]: evaluated gradient
        """
        import numpy as np
        # calculate r and r^2 vector
        r_sq = (X**2 + Y**2 + Z**2)
        r = np.sqrt(r_sq)
        # calculate pi^(3/2)
        pi_sq3 = np.pi**(3/2)
        
        # evaluate gradient in 3 dimensions
        grad = -2*np.exp(-r/h**2) / (h**5*pi_sq3)
        delwx = grad * X
        delwy = grad * Y
        delwz = grad * Z
        
        # return x, y, z gradient
        return [delwx, delwy, delwz]
    
    def particle_separation(ri, rk):
        """
        evaluate particle-particle separations in 3 dimensions
        ------------------------------------------------------
        ri, rk [numpy array]: Nx3 and Nx3 arrays storing positions
        -------------------------------------------------------------
        OUTPUT [numpy array]: x, y, z particle - particle separations
        """
        # create Nx1 and Mx1 position arrays storing x, y, z positions
        rix, rkx = ri[:, 0:1], rk[:, 0:1]
        riy, rky = ri[:, 1:2], rk[:, 1:2]
        riz, rkz = ri[:, 2:3], rk[:, 2:3]
        
        # evaluate particle - particle separations with vectorized algorithm
        delx = rix - rkx.T
        dely = riy - rky.T
        delz = riz - rkz.T
        
        # return x, y, z particle separations
        return [delx, dely, delz]

    def density(r, pos, m, h):
        """
        evaluate mass density
        ---------------------
        r [numpy array]: Nx3 position array
        pos [numpy array]: Mx3 postion array
        m [float]: single particle mass
        h [float]: smoothing length
        ---------------------------------------
        OUTPUT [numpy array]: Nx1 density array
        """
        import numpy as np
        # evaluate particle separations
        dx, dy, dz = particle_separation(r, pos)
        # evaluate smoothing function using particle separations
        w = m * smoothing_kernel(dx, dy, dz, h)
        # calculate densities and store in Nx1 array
        rho = np.sum(w, 0).reshape((r.shape[0],1))
        
        return rho

    def pressure(rho, k, n):
        """
        evaluates the pressure from mass density
        ----------------------------------------
        rho [numpy array]: Nx1 pressure array
        k [float]: state equation constant
        n [float]: polytropic index
        ----------------------------------------
        OUTPUT [numpy array]: Nx1 pressure array
        """
        import numpy as np
        return k * rho*np.exp(1 + 1/n)
    
    def acceleration(position, velocity, mass, h, k, n, lamda, nu):
        """
        calculates acceleration onto particles from pressure gradient
        -------------------------------------------------------------
        position [numpy array]: Nx3 matrix of positions
        velocity [numpy array]: Nx3 matrix of velocities
        mass [float]: mass of single particle
        h [float]: smoothing length
        k [float]: state equation constant
        n [float]: polytropic index
        lamda [float]: gravitational term specified in Philip Mocz derivation
        nu [float]: damping coefficient to achieve equilibrium state
        --------------------------------------------
        OUTPUT [numpy array]: Nx3 acceleration array
        """
        # calculate density
        rho = density(position, position, mass, h)

        # calculate pressure
        P = pressure(rho, k, 1)

        # calculate pairwise separation and gradients
        dx, dy, dz = particle_separation(position, position)
        wx, wy, wz = gradient_smooth(dx, dy, dz, h)

        # compute acceleration components
        accelx = - np.sum(mass * ( (P/rho**2) + (P.T/rho.T**2) ) * wx, 1)
        accely = - np.sum(mass * ( (P/rho**2) + (P.T/rho.T**2) ) * wy, 1)
        accelz = - np.sum(mass * ( (P/rho**2) + (P.T/rho.T**2) ) * wz, 1)

        # create acceleration vector
        acceleration = np.vstack((accelx, accely, accelz)).T

        # add external force and damping term
        acceleration += - lamda*position - nu*velocity
        
        return acceleration

    import numpy as np
    # import gamma function
    from scipy.special import gamma
    
    ##########################
    # SIMULATION CONFIGURATION
    ##########################
    
    # define simulation constants; ensure quantities are integers as required
    N = int(N)
    timesteps = int(timesteps)
    # definition of h
    h = .04*(np.sqrt(N/1000))

    # definition of lambda from Philip Mocz, 2020
    lamda = 2*k*(1+n)*np.pi**(-3/(2*n)) * (M*gamma(5/2+n)/R**3/gamma(1+n))**(1/n) / R**2
    
    # mass of each particle (mass of star / number of particles)
    mass = M/N
    
    # initial acceleration
    accel = acceleration(position, velocity, mass, h, k, n, lamda, nu)
    # initial density
    rho = density(position, position, mass, h)

    # set up position and density arrays for storing
    pos_arr = position
    rho_arr = rho
    
    ######################
    # SIMULATION MAIN LOOP
    ######################
    
    print('simulation running....  /ᐠ –ꞈ –ᐟ\<[pls be patient]')
    
    for i in range(timesteps):
        # calculate timestep velocities
        velocity += accel * dt/2
        
        # calculate new position using timestep velocities
        position += velocity * dt
        
        # calculate timestep densities
        rho = density(position, position, mass, h)
        
        # save positions and densities
        pos_arr = np.append(pos_arr, position, 0)
        rho_arr = np.append(rho_arr, rho, 0)
        
        # update accelerations
        accel = acceleration(position, velocity, mass, h, k, n, lamda, nu)
        
        # update velocities
        velocity += accel * dt/2
        
        if i == int(timesteps/4):
            print('25% ...')
        if i == int(timesteps/2):
            print('50% ...')
        if i == int(timesteps * .75):
            print('75% ...')
    
    print('simulation complete [yay!!! (ﾐΦ ﻌ Φﾐ)✿ *ᵖᵘʳʳ*]')    
    
    # return positions and densities
    return pos_arr, rho_arr