import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import utility, settings

K_STOP = 20
LOG_RES = []

def main():
    v_real = utility.readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    L_real = utility.get_l(s_real, recalculate=False)
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
        # vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (CG, k={})".format(k))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig(os.path.join('pb2', 'model_v0.png'))
    plt.clf()
    
    try:
        pk = np.negative(utility.grad(t_obs, s_model, L))
        while k < K_STOP and res > 1e-2:
            alphak = utility.get_proper_alpha(t_obs, s_model, L, pk, method='cube_quad')
            s_model_new = np.add(
                s_model,
                np.multiply(alphak, pk)
            )
            s_model_new = utility.smooth_map(s_model_new, kernel_size=50)
            L_new = utility.get_l(s_model_new, recalculate=True)
            gradk = utility.grad(t_obs, s_model, L)
            gradk1 = utility.grad(t_obs, s_model_new, L_new)
            betak1 = np.matmul(gradk1.T, gradk1)/np.matmul(gradk.T, gradk)

            pk = np.add(np.negative(gradk1), np.multiply(betak1, pk))
            pk = np.divide(pk, np.max(np.abs(pk)))
            pk = utility.smooth_map(pk, mode='uniform', kernel_size=(10, 10))

            plt.imshow(
                pk.reshape(settings.NY, settings.NX), cmap='jet',
                extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
            )
            a = plt.colorbar()
            a.set_label('$p_k$')
            plt.title("$p_k$ (SD, k={})".format(k))
            plt.xlabel("$x$")
            plt.ylabel("$y$")
            plt.savefig(os.path.join('pb2', 'pk_k%d.png'%k))
            plt.clf()
            
            res = utility.get_r(t_obs, s_model_new, L_new)
            LOG_RES.append(res)

            print(k, res, alphak, np.min(s_model), np.max(s_model))
            s_model = s_model_new
            L = L_new
            k += 1

    except KeyboardInterrupt:
        pass
    

    plt.plot(LOG_RES)
    plt.xlabel('k')
    plt.ylabel('Residual')
    plt.title("Residual over iterations")
    plt.savefig(os.path.join('pb2', 'res.png'))

    plt.clf()

    v_model = np.divide(1.0, s_model)
    plt.imshow(
        v_model.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        # vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (SD, k={})".format(k))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig(os.path.join('pb2', 'model_v.png'))