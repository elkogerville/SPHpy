def collision(N, M, xpos, ypos, zpos, xvel, yvel, zvel):
    import numpy as np
    
    # initialize star at rest
    posN = np.random.randn(N, 3)
    velN = np.zeros((N, 3))
    
    # initialize merging star
    posM = np.random.randn(M, 3)
    posM[:, 0] = posM[:, 0] + xpos
    posM[:, 1] = posM[:, 1] + ypos
    posM[:, 2] = posM[:, 2] + zpos
    velM = np.ones((M, 3))
    velM[:, 0] = velM[:, 0] * xvel
    velM[:, 1] = velM[:, 1] * yvel
    velM[:, 2] = velM[:, 2] * zvel
    
    # append positions and velocities
    posN = np.append(posN, posM, 0)
    velN = np.append(velN, velM, 0)
    
    # total number of particles
    Ntot = N + M
    
    return posN, velN, Ntot