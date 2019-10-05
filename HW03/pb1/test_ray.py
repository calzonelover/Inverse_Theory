import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray


if __name__ == "__main__":
    # Test ray tracing
    for i_y in range(settings.N_RECEIVER):
        y = settings.DX + settings.DX * i_y
        one_ray_map = ray.ray_length(-10,10,1010,y).reshape(50, 50)
        plt.imshow(one_ray_map, extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
        a = plt.colorbar()
        a.set_label('Ray length (m)')
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title("Ray Tracing (source: ({},{}), receiver({},{}))".format(-10, 510, 1010, y))
        plt.savefig("img/ray_{:02.0f}.png".format(i_y))
        plt.clf()