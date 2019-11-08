import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, uniform_filter

import utility, settings

K_STOP = 20
LOG_RES = []

def main():
    v_real = utility.readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    L = utility.get_l(s_real, recalculate=False)
    t_obs = np.matmul(L, s_real)

    s_model = utility.smooth_map(s_real)

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
    plt.title("Model Velocity (CG, k={})".format(k))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.show()
    plt.clf()
    
    try:
        pk = -1.0*utility.grad(t_obs, s_model, L)
        while k < K_STOP and res > 1e-2:
            alphak = utility.get_proper_alpha(t_obs, s_model, L, pk, method='cube_quad')
            s_model_new = np.add(
                s_model,
                np.multiply(alphak, pk)
            )
            s_model_new = utility.smooth_map(s_model_new, kernel_size=10)
            gradk = utility.grad(t_obs, s_model, L)
            gradk1 = utility.grad(t_obs, s_model_new, L)
            betak1 = np.matmul(gradk1.T, gradk1)/np.matmul(gradk.T, gradk)
            pk = np.add(np.negative(gradk1), np.multiply(betak1, pk))

            print('alphak:{}, min(pk):{}, max(pk):{}'.format(alphak, min(pk), max(pk)))
            
            res = utility.get_r(t_obs, s_model_new, L)
            LOG_RES.append(res)

            print('k:{}, res:{}, min_s:{}, max_s:{}'.format(k, res, min(s_model), max(s_model)))
            s_model = s_model_new
            k += 1
    except KeyboardInterrupt:
        pass
    

    plt.plot(LOG_RES)
    plt.xlabel('k')
    plt.ylabel('Residual')
    plt.title("Residual over iterations")
    plt.show()
    plt.clf()

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
    plt.show()