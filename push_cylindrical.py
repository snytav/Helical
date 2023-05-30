import numpy as np
from monkey_boris import push_Boris

def cart2pol(x, y):
    theta = np.arctan2(y, x)
    r     = np.hypot(x, y)
    return theta, r

def pol2cart(theta, r):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def vector_cart2pol(vx,vy,theta):
    vr  = np.cos(theta)*vx+ np.sin(theta)*vy
    vth = -np.sin(theta) * vx + np.cos(theta) * vy


def vector_pol2cart(vr, vth, theta):
    vx  = np.cos(theta) * vr - np.sin(theta) * vth
    vy =  np.sin(theta) * vr + np.cos(theta) * vth


def get_polar_field_2D(lc,r0,dr,data):
    i = int(lc[0])
    di = lc[0] - i

    j = int(lc[1])
    dj = lc[1] - j
    # compute correction factors
    rj = r0 + j * dr
    f1 = (rj + 0.5 * dj * dr) / (rj + 0.5 * dr)
    f2 = (rj + 0.5 * (dj + 1) * dr) / (rj + 0.5 * dr)

    #gather electric field onto particle position
    val = (data[i][j] * (1 - di) * (1 - dj) * f2 +
          data[i + 1][j] * (di) * (1 - dj) * f2 +
          data[i + 1][j + 1] * (di) * (dj) * f1 +
          data[i][j + 1] * (1 - di) * (dj) * f1)

    return val


def get_field_at_point(r,th,z,Er,Eth,Ez):


def push_cyl(rr,theta,zz,vrr,vtheta,vzz,Er_spiral,Etheta_spiral,Ez_spiral,
             Br_spiral,Btheta_spiral,Bz_spiral,
             r_linspace,theta_linspace,z_linspace):
    for r,th,z, vr,vth,vz in zip(rr,theta,zz,vrr,vtheta,vzz):
        x,y = pol2cart(th,r)
        Ex,Ey = vector_pol2cart(er,et)
        Bx,By = vector_pol2cart(br,bt)


