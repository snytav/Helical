import numpy as np

def generate_initial_particle_distribution(N,rmax,zmax,dr,dth,dz):
    r = 0.5*np.ones(N) # rmax*np.random.random(N)
    theta = np.zeros(N) #2*np.pi*np.random.random(N)
    z = zmax * np.random.random(N)
    vr = np.zeros(N) # np.random.normal(0, dr, N)
    vz = np.zeros(N) #np.random.normal(0, dz, N)
    vth = 2*np.pi*np.random.normal(0, dth, N)
    vth = np.zeros(N) # np.multiply(r,vth)

    return [r,theta,z,vr,vth,vz]
