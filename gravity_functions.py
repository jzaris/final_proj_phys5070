import numpy as np
import matplotlib.pyplot as plt


def momentum_kick(pos, masses, n, N, dt, step, G = 6.67430*10**-11):
    """
    Calculated length 3 array containing momentum kick in (x,y,z) on particle n due to all other particles
    
    Arguments:
    ====
    pos: see 'evolve'
    masses: see 'evolve'
    n: index of mass whose momentum kick is calculated
    N: see 'evolve'
    dt: see 'evolve'
    step: see 'evolve'
    G: see 'evolve'
    
    Returns:
    =====
    Momentum kick on particle n in each direction (x,y,z)
    """
    nx = pos[step-1, 0, n]
    ny = pos[step-1, 1, n]
    nz = pos[step-1, 2, n]
    dpx=0
    dpy=0
    dpz=0
    for m in range(N):
        if m != n:
            mx = pos[step-1, 0, m]
            my = pos[step-1, 1, m]
            mz = pos[step-1, 2, m]
            rad2 = (mx-nx)**2+(my-ny)**2+(mz-nz)**2
            total_dp = dt*G*masses[m]*masses[n]/rad2
            dpx += total_dp * (mx-nx)/np.sqrt(rad2)
            dpy += total_dp * (my-ny)/np.sqrt(rad2)
            dpz += total_dp * (mz-nz)/np.sqrt(rad2)
    return dpx, dpy, dpz

def calc_energy(pos, vel, masses, N, step, G = 6.67430*10**-11):
    """
    Calculates the total energy of gravitational system at timestep given by 'step'
    
    Arguments:
    ====
    pos: see 'evolve'
    vel: see 'evolve'
    masses: see 'evolve'
    N: see 'evolve'
    step: timestep at which to calculate the energy of the system
    G: see 'evolve'
    
    Returns:
    =====
    Energy (kinetic plus potential) of gravitational system
    """
    PE=0
    KE=0
    for n in range(N):
        nx = pos[step, 0, n]
        ny = pos[step, 1, n]
        nz = pos[step, 2, n]
        KE += 0.5*masses[n]*(vel[0,n]**2+vel[1,n]**2+vel[2,n]**2)
        for m in range(N):
            if m != n:
                mx = pos[step, 0, m]
                my = pos[step, 1, m]
                mz = pos[step, 2, m]
                rad = np.sqrt((mx-nx)**2+(my-ny)**2+(mz-nz)**2)
                PE += -G*masses[m]*masses[n]/rad
    return KE+0.5*PE

def evolve(pos, vel, masses, N, steps, dt, G = 6.67430*10**-11):
    """
    Given initial parameters of gravitational system, returns array containing cooridnates of each mass at each timestep, as well as the energy at each timestep
    
    Arguments:
    ====
    pos: Array of positions of each mass - index 0: timestep, index 1: x, y, or z, index 2: mass index.  Passed array is empty except for the first timestep.
    vel: Array of velocities of each mass - ### index 0: Vx, Vy, or Vz, index 1: mass index 
    masses: Array of masses
    N: number of masses
    dt: timestep
    G: gravitational constant (set to one in certain examples)
    
    Returns:
    =====
    Filled position array, energy at each timestep
    """
    energy = np.empty(steps)
    energy[0] = calc_energy(pos, vel, masses, N, 0, G)
    for step in range(1, steps):
        dp = np.zeros([3, N])
        for n in range(N):
            dp[:,n] = momentum_kick(pos, masses, n, N, dt, step, G)
            vel[:,n] = vel[:,n] + dp[:,n]/masses[n]
            pos[step,:,n] = pos[step-1,:,n] + vel[:,n]*dt
        energy[step] = calc_energy(pos, vel, masses, N, step, G)
    return pos, energy
