import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings

REPORT_LOG = {
    'alphas': [],
    'norm_model': [],
    'norm_res': [],
}

def readraw(filename):
    f = open(settings.FILENAME, "r")
    return np.fromfile(f, dtype=np.float32)

def readmap(filename):
    a = readraw(filename)
    a = a.reshape(settings.NX, settings.NY)
    return a.T

if __name__ == "__main__":
    '''
    # plot real velocity
    raw_map = readmap(filename=settings.FILENAME)
    plt.imshow(
        raw_map, cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Real Velocity")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig("real_v.png")
    '''