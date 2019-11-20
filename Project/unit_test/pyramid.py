import numpy as np
import os
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import environment, utility, settings

MODEL = 'initial' #  (initial, real)

def main():
    model = environment.get_system(model=MODEL)
    plt.imshow(
        model.reshape(settings.NY, settings.NX), cmap='summer',
        extent=[0, settings.DX*settings.NX, 0, settings.DX*settings.NY],
        origin='bottom left',
    )
    a = plt.colorbar()
    a.set_label('1 / Path length ($m^{-1}$)')
    name_title = "Real" if MODEL == 'real' else "Initial_Model"
    plt.title(name_title)
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    # plt.show()
    plt.savefig(os.path.join('unit_test', '{}.png'.format(MODEL)))