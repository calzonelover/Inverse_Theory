import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import utility, settings

SOURCE_X = 40.0
SOURCE_Y = 40.0

def main():
    v = utility.readraw(filename=settings.FILENAME)
    s = 1.0/v
    T = utility.get_travel_time(s, SOURCE_X, SOURCE_Y)
    plt.imshow(
        T.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
    )
    a = plt.colorbar()
    a.set_label('$T$')
    plt.title("Travel time [Source at (x,y)=({}, {})]".format(SOURCE_X, SOURCE_Y))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.show()
    # plt.savefig("unit_test/travel_time_x{}_y{}.png".format(int(SOURCE_X), int(SOURCE_Y)))