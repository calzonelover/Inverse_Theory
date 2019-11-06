import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, uniform_filter

import utility, settings

ALPHAS = np.logspace(-3, 3, num=50, base=10)

def main():
    v_real = utility.readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    L = utility.get_l(s_real, recalculate=False)
    t_obs = np.matmul(L, s_real)

    s_real_padded = np.pad(s_real.reshape(settings.NX, settings.NY).T, (100,100), 'edge')
    s_model = uniform_filter(s_real_padded, size=60)[100:-100, 100:-100].reshape(settings.NX*settings.NY)

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

    L = L.astype(np.float16)
    t_obs = t_obs.astype(np.float16)

    for alpha in ALPHAS:
        print(alpha)
        inn = np.linalg.inv(np.add(
            np.matmul(L.T, L),
            np.multiply(alpha, np.identity(s_model.shape, dtype=np.float16))
        ))
        s_model = np.matmul(
            inn,
            np.matmul(L.T, t_obs)
        )

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
        plt.show()
        plt.clf()