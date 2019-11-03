import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import utility, settings, ray

SOURCE_X = 1800.0
SOURCE_Y = 0.0

RECEIVER_X = 8000.0
RECEIVER_Y = 0.0

def main():
    v = utility.readraw(filename=settings.FILENAME)
    s = 1.0/v
    T = utility.get_travel_time(s, SOURCE_X, SOURCE_Y)
    L = ray.curved_ray(SOURCE_X, SOURCE_Y, RECEIVER_X, RECEIVER_Y, T)
    plt.imshow(
        L.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
    )
    a = plt.colorbar()
    a.set_label('$L$')
    plt.title("Ray [Source: ({}, {}), Receiver: ({}, {})]".format(SOURCE_X, SOURCE_Y, RECEIVER_X, RECEIVER_Y))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.show()
    # plt.savefig("unit_test/travel_time_x{}_y{}.png".format(int(SOURCE_X), int(SOURCE_Y)))