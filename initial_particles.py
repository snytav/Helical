import numpy as np

def generate_initial_particle_distribution(N,rmax,zmax,dr,dth,dz):
    r = rmax*np.random.random(N)
    theta = 2*np.pi*np.random.random(N)
    z = zmax * np.random.random(N)
    vr = np.random.normal(0, dr, N)
    vz = np.random.normal(0, dz, N)
    vth = 2*np.pi*np.random.normal(0, dth, N)
    vth = np.multiply(r,vth)

    return [r,theta,z,vr,vth,vz]
