import numpy as np


def push_Boris(x,v,qm,E,B,dt):
    v_minus = v + qm * E * dt / 2
    t = qm * B * dt / 2
 #   t = t.reshape(1, v.shape[1])
   # t = np.repeat(t, v.shape[0], axis=0)
    v_prime = v_minus + np.cross(v_minus, t);
    cp = np.sum(np.multiply(t, t), axis=1)
    cp1 = np.add(np.ones_like(cp),cp)
    cp2 = cp1.reshape(1, cp1.shape[0])
    cp3 = np.repeat(cp2, t.shape[1], axis=0)
    cp3 = cp3.T
    s = 2.0 * np.divide(t,cp3)

    #s = np.repeat(s, v_prime.shape[0], axis=0)
    v_plus = v_minus + np.cross(v_prime, s)
    v = v_plus + qm * E * dt / 2

    x += - v * dt*2
    return x,v
