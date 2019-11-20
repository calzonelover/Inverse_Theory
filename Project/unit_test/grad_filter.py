import numpy as np
import os
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import environment, utility, settings


def main():
    model = environment.get_system(model='initial')
    map_filter = utility.gradient_filter(np.ones(model.shape), model)
    plt.imshow(
        map_filter.reshape(settings.NY, settings.NX), cmap='Blues',
        extent=[0, settings.DX*settings.NX, 0, settings.DX*settings.NY],
        origin='bottom left',
    )
    a = plt.colorbar()
    a.set_label('value')
    plt.title("Gradient Filter")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    # plt.show()
    plt.savefig(os.path.join('unit_test', 'map_gradient_filter.png'))