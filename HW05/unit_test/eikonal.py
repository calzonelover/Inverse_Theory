import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import utility, settings


def main():
    v = utility.readraw(filename=settings.FILENAME)
    s = 1.0/v
    T = utility.get_travel_time(s, 6000.0, 2000.0, 2)
    plt.imshow(
        T.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
    )
    a = plt.colorbar()
    a.set_label('$T$')
    plt.title("Travel time")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.show()
    # plt.savefig("real_v.png")