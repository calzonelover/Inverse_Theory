import numpy as np
import os
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import utility, settings

K_STOP = 10
LOG_RES = []

def main():
    v_real = utility.readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    L_real = utility.get_l(s_real, recalculate=True)
    t_obs = np.matmul(L_real, s_real)

    s_model = utility.smooth_map(s_real, mode='uniform')
    L = utility.get_l(s_model, recalculate=True)

    k = 0
    res = utility.get_r(t_obs, s_model, L)
    LOG_RES.append(res)

    v_model = np.divide(1.0, s_model)
    plt.imshow(
        v_model.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (SD, k={})".format(k))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig(os.path.join('pb1', 'model_v0.png'))
    plt.clf()
    
    try:
        while k < K_STOP and res > 1e-2:
            pk = np.negative(utility.grad(t_obs, s_model, L))
            pk = np.divide(pk, np.max(np.abs(pk)))
            pk = utility.smooth_map(pk, mode='uniform', kernel_size=(10, 40))
            # pk = utility.smooth_map(pk, mode='gaussian', kernel_size=(20, 20), pad_model='zero')

            plt.imshow(
                pk.reshape(settings.NY, settings.NX), cmap='jet',
                extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
            )
            a = plt.colorbar()
            a.set_label('$p_k$')
            plt.title("$p_k$ (SD, k={})".format(k))
            plt.xlabel("$x$")
            plt.ylabel("$y$")
            plt.savefig(os.path.join('pb1', 'pk_k%d.png'%k))
            plt.clf()

            alphak = utility.get_proper_alpha(t_obs, s_model, L, pk, method='backtrack')
            s_model = np.add(
                s_model,
                np.multiply(alphak, pk)
            )
            L = utility.get_l(s_model, recalculate=True)
            res = utility.get_r(t_obs, s_model, L)
            LOG_RES.append(res)
            print(k, res, alphak, np.min(s_model), np.max(s_model))
            k += 1
    except KeyboardInterrupt:
        pass

    plt.plot(LOG_RES)
    plt.xlabel('k')
    plt.ylabel('Residual')
    plt.title("Residual over iterations")
    plt.savefig(os.path.join('pb1', 'res.png'))
    plt.clf()

    v_model = np.divide(1.0, s_model)
    print('visual', k, res, alphak, np.min(v_model), np.max(s_model), np.max(v_model), np.min(s_model))
    plt.imshow(
        v_model.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (SD, k={})".format(k))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig(os.path.join('pb1', 'model_v.png'))
    print('The process is fully finish')