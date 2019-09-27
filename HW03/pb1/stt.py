import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray

def readmap(filename):
    f = open(settings.FILENAME, "r")
    a = np.fromfile(f, dtype=np.float32)
    a = a.reshape(settings.NX, settings.NY)
    return a.T

if __name__ == "__main__":
    raw_map = readmap(filename=settings.FILENAME)
    for i_y in range(40):
        y = 40 + settings.DX * i_y
        one_ray_map = ray.ray_length(-10,210,1010,y+10).reshape(50, 50)
        plt.imshow(one_ray_map, extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
        plt.colorbar()
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.savefig("img/ray_{:02.0f}.png".format(i_y))
        plt.clf()
        # plt.show()

    # visualize
    plt.imshow(raw_map, cmap='jet', extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Real Velocity")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig("real_v.png")
    # plt.show()