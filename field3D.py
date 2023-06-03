import numpy as np
from mayavi.mlab import *
import matplotlib.pyplot as plt

def volume_field_plot(s):
# u =    np.sin(np.pi*x) * np.cos(np.pi*z)
# v = -2*np.sin(np.pi*y) * np.cos(2*np.pi*z)
# w = np.cos(np.pi*x)*np.sin(np.pi*z) + np.cos(np.pi*y)*np.sin(2*np.pi*z)
    plt.figure()
    v = volume_slice(s, plane_orientation='x_axes', slice_index=10)
    outline()
    colorbar(v)
    axes()

    plt.show()
    qq = 0

plt.show()
qq = 0