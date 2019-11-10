import os
import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, uniform_filter

import utility, settings

ALPHAS = np.logspace(-3.0, 3.0, num=20, base=10)

def main():
    v_real = utility.readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    L_real = utility.get_l(s_real, recalculate=True)
    t_obs = np.matmul(L_real, s_real)

    s_model = utility.smooth_map(s_real, mode='uniform')
    L = utility.get_l(s_model, recalculate=True)

    v_model = np.divide(1.0, s_model)
    plt.imshow(
        v_model.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Initial Model Velocity")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.show()
    plt.clf()

    for alpha in ALPHAS:
        print(alpha)
        r = utility.get_r(t_obs, s_model, L, magnitude=False)
        s_model = np.matmul(
            np.linalg.inv(
                np.add(
                    np.matmul(L.T, L),
                    np.multiply(alpha, np.identity(settings.NX*settings.NY, dtype=float))
                )
            ),
            np.matmul(L.T, r)
        )
        L = utility.get_l(s_model, recalculate=True)

        v_model = np.divide(1.0, s_model)
        plt.imshow(
            v_model.reshape(settings.NY, settings.NX), cmap='jet',
            extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
            vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
        )
        a = plt.colorbar()
        a.set_label('$v$')
        plt.title("Model Velocity (REG, alpha={})".format(alpha))
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.savefig(os.path.join('pb3'), 'model_v_alpha{0:.03d}.png'.format(alpha))
        plt.clf()