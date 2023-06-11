import numpy as np
from monkey_boris import push_Boris
# from cyl_plot import  draw_particles_3D
import matplotlib.pyplot as plt

def get_theta(x,y):
    th = np.arctan2(y, x)
    if th < 0.0:
        th += 2 * np.pi
    if th > 2 * np.pi:
        th -= 2 * np.pi
    return th

def cart2pol(x, y):
    theta = get_theta(x,y)
    r     = np.hypot(x, y)
    return r,theta

def pol2cart( r,theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def vector_cart2pol(vx,vy,theta):
    vr  = np.cos(theta)*vx+ np.sin(theta)*vy
    vth = -np.sin(theta) * vx + np.cos(theta) * vy
    return vr,vth


def vector_pol2cart(vr, vth, theta):
    vx  = np.cos(theta) * vr - np.sin(theta) * vth
    vy =  np.sin(theta) * vr + np.cos(theta) * vth
    return vx,vy

def XtoL(x,x0,dh):
    lc = np.divide(x - x0,dh)

    return lc

def interpolate(data,i,j,k,di,dj,f1,f2):
    if j == data.shape[1]-1:
       j = 0

    val = (data[i][j][k] * (1 - di) * (1 - dj) * f2 +
           data[i + 1][j][k] * (di) * (1 - dj) * f2 +
           data[i + 1][j + 1][k] * (di) * (dj) * f1 +
           data[i][j + 1][k] * (1 - di) * (dj) * f1)
    return val


def get_polar_field_2D(lc,r0,dr,data):
    i = int(lc[0])
    di = lc[0] - i
    k = int(lc[2])
    dk = lc[2] - k

    j = int(lc[1])

    dj = lc[1] - j
    # compute correction factors
    rj = r0 + j * dr
    f1 = (rj + 0.5 * dj * dr) / (rj + 0.5 * dr)
    f2 = (rj + 0.5 * (dj + 1) * dr) / (rj + 0.5 * dr)

    #gather electric field onto particle position
    val = interpolate(data,i,j,k,di,dj,f1,f2)
    val1 = interpolate(data,i,j,k+1,di,dj,f1,f2)

    t = val*(1-dk)+ val1*dk
    return t


def reflect_particle(r,th,z,vr,vz,rmax,zmax):
    if th > 2*np.pi:
        th = th - 2*np.pi

    if th < 0:
        th = th + 2 * np.pi

    if r < 0:
       r = -r
       vr = -vr

    if r > rmax:
        r = 2*rmax - r
        vr = -vr

    if z < 0:
       z = -z
       vz = -vz

    if z > zmax:
        z = 2*zmax - z
        vz = -vz
    return [r,th,z,vr,vz]

def reflect(rr,theta,zz,vrr,vzz,rmax,zmax):

    for i in range(rr.shape[0]):
        r  = rr[i]
        th = theta[i]
        z  = zz[i]
        vr = vrr[i]
        vz = vzz[i]
        r,th,z,vr,vz = reflect_particle(r,th,z,vr,vz,rmax,zmax)
        rr[i]    = r
        theta[i] = th
        zz[i]    = z
        vrr[i]   = vr
        vzz[i]   = vz

    return [rr, theta, zz, vrr, vzz]

from polar_trace import draw_polar_trace


def cyl2cart_coordinates_fields(rr,theta,zz,vrr,vtheta,vzz,Er_spiral,Etheta_spiral,Ez_spiral,
             Br_spiral,Btheta_spiral,Bz_spiral,
             r_linspace,theta_linspace,z_linspace,dt,qm):

    dr     = r_linspace[1] - r_linspace[0]
    dtheta = theta_linspace[1] - theta_linspace[0]
    dz     = z_linspace[1] - z_linspace[0]
    dh = np.array([dr,dtheta,dz])
    x0 = np.zeros(3)

    rmax = np.max(r_linspace)
    zmax = np.max(z_linspace)

    xx = []
    vv = []
    EE = []
    BB = []
    for i in range(rr.shape[0]):
        r   = rr[i]
        th  = theta[i]
        z   = zz[i]
        vr  = vrr[i]
        vth = vtheta[i]
        vz  = vzz[i]
        xcyl = np.array([r,th,z])
        lc   = XtoL(xcyl,x0,dh)

        er = get_polar_field_2D(lc, r, dr, Er_spiral)
        et = get_polar_field_2D(lc, r, dr, Etheta_spiral)
        ez = get_polar_field_2D(lc, r, dr, Ez_spiral)
        br = get_polar_field_2D(lc, r, dr, Br_spiral)
        bt = get_polar_field_2D(lc, r, dr, Btheta_spiral)
        bz = get_polar_field_2D(lc, r, dr, Bz_spiral)

        x,y = pol2cart(r,th)
        vx,vy = vector_pol2cart(vr, vth,th)
        Ex,Ey = vector_pol2cart(er,et,th)
        Bx,By = vector_pol2cart(br,bt,th)

        x = np.array([x,y,z])
        v = np.array([vx, vy, vz])
        E = np.array([Ex, Ey, ez])
        B = np.array([Bx, By, bz])
        xx.append(x)
        vv.append(v)
        EE.append(E)
        BB.append(B)

    xx = np.array(xx)
    vv = np.array(vv)
    EE = np.array(EE)
    BB = np.array(BB)

    return xx,vv,EE,BB

def cart2cyl_coordinates_velocities(xx1,vv1,th):

    rr     = np.zeros(xx1.shape[0])
    theta  = np.zeros(xx1.shape[0])
    vrr    = np.zeros(xx1.shape[0])
    vtheta = np.zeros(xx1.shape[0])
    zz     = np.zeros(xx1.shape[0])
    vzz    = np.zeros(xx1.shape[0])

    for i in range(xx1.shape[0]):
        x1 = xx1[i]
        v1 = vv1[i]
        r1, th1 = cart2pol(x1[0], x1[1])
        vr1, vth1 = vector_cart2pol(v1[0], v1[1], th[i])
        rr[i] = r1
        theta[i] = th1
        zz[i] = x1[2]
        vrr[i] = vr1
        vtheta[i] = vth1
        vzz[i] = v1[2]

    return rr, theta, zz, vrr, vtheta, vzz


def push_cyl(rr,theta,zz,vrr,vtheta,vzz,Er_spiral,Etheta_spiral,Ez_spiral,
             Br_spiral,Btheta_spiral,Bz_spiral,
             r_linspace,theta_linspace,z_linspace,dt,qm):

    rmax = np.max(r_linspace)
    zmax = np.max(z_linspace)

    dr = r_linspace[1] - r_linspace[0]
    dtheta = theta_linspace[1] - theta_linspace[0]
    dz = z_linspace[1] - z_linspace[0]
    dh = np.array([dr, dtheta, dz])
    x0 = np.zeros(3)

    rmax = np.max(r_linspace)
    zmax = np.max(z_linspace)

    # for i in range(rr.shape[0]):
    #     r = rr[i]
    #     th = theta[i]
    #     z = zz[i]
    #     vr = vrr[i]
    #     vth = vtheta[i]
    #     vz = vzz[i]
    #     xcyl = np.array([r, th, z])
    #     lc = XtoL(xcyl, x0, dh)
    #
    #     er = get_polar_field_2D(lc, r, dr, Er_spiral)
    #     et = get_polar_field_2D(lc, r, dr, Etheta_spiral)
    #     ez = get_polar_field_2D(lc, r, dr, Ez_spiral)
    #     br = get_polar_field_2D(lc, r, dr, Er_spiral)
    #     bt = get_polar_field_2D(lc, r, dr, Etheta_spiral)
    #     bz = get_polar_field_2D(lc, r, dr, Ez_spiral)
    #
    #     x, y = pol2cart(r, th)
    #     vx, vy = vector_pol2cart(vr, vth, th)
    #     Ex, Ey = vector_pol2cart(er, et, th)
    #     Bx, By = vector_pol2cart(br, bt, th)
    #
    #     x = np.array([x, y, z])
    #     v = np.array([vx, vy, vz])
    #     E = np.array([Ex, Ey, ez])
    #     B = np.array([Bx, By, bz])

    x, v, E, B =     cyl2cart_coordinates_fields(rr, theta, zz,
                                 vrr, vtheta, vzz,
                                 Er_spiral, Etheta_spiral, Ez_spiral,
                                 Br_spiral, Btheta_spiral, Bz_spiral,
                                 r_linspace, theta_linspace, z_linspace, dt, qm)

    x1, v1 = push_Boris(x, v, qm, E, B, -dt*0.5)

    rr, th, zz, vrr, vtheta, vzz = cart2cyl_coordinates_velocities(x1, v1,theta)

    rr, theta, zz, vrr, vzz = reflect(rr, th, zz, vrr, vzz, rmax, zmax)

    return rr, theta, zz, vrr, vtheta, vzz,x1
