import os
import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, uniform_filter

import utility, settings

K_STOP = 10
ALPHAS = np.logspace(-3.0, 3.0, num=20, base=10)
REPORT_LOG = {
    'alphas': [],
    'norm_model': [],
    'norm_res': [],
}

def main():
    v_real = utility.readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    L_real = utility.get_l(s_real, recalculate=True)
    t_obs = np.matmul(L_real, s_real)

    v_real = np.divide(1.0, s_real)
    plt.imshow(
        v_real.reshape(settings.NX, settings.NY), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Resized Real Velocity")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig(os.path.join('pb3', 'model_v_resize_real.png'))
    plt.clf()

    s_model0 = utility.smooth_map(s_real, mode='uniform')
    L0 = utility.get_l(s_model0, recalculate=True)

    v_model0 = np.divide(1.0, s_model0)
    plt.imshow(
        v_model0.reshape(settings.NY, settings.NX), cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Initial Model Velocity")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig(os.path.join('pb3', 'model_v0.png'))
    plt.clf()

    for i_alpha, alpha in enumerate(ALPHAS):
        print("\n\n ALPHA: {}".format(alpha))
        s_model = s_model0
        L = L0
        LOG_RES = []
        dir_save = os.path.join('pb3', 'alpha%d'%i_alpha)
        if not os.path.exists(dir_save):
            os.makedirs(dir_save)
        for k in range(K_STOP):
            r = utility.get_r(t_obs, s_model, L, magnitude=False)
            pk = np.negative(np.matmul(
                np.linalg.inv(
                    np.add(
                        np.matmul(L.T, L),
                        np.multiply(alpha, np.identity(settings.NX*settings.NY, dtype=float))
                    )
                ),
                np.matmul(L.T, r)
            ))
            # pk = np.divide(pk, np.max(np.abs(pk)))
            pk = utility.smooth_map(pk, mode='uniform', kernel_size=(10, 40))

            s_model = np.add(
                s_model,
                pk
            )
            L = utility.get_l(s_model, recalculate=True)
            res = utility.get_r(t_obs, s_model, L)
            LOG_RES.append(res)

            plt.imshow(
                pk.reshape(settings.NY, settings.NX), cmap='jet',
                extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
            )
            a = plt.colorbar()
            a.set_label('$p_k$')
            plt.title("$p_k$ (LM, alpha={0:03f}, k={1})".format(alpha, k))
            plt.xlabel("$x$")
            plt.ylabel("$y$")
            plt.savefig(os.path.join(dir_save, 'pk_k%d.png'%k))
            plt.clf()

            print(alpha, k, res, np.min(s_model), np.max(s_model), np.min(pk), np.max(pk))
        ## end for k
        REPORT_LOG['alphas'].append(alpha)
        REPORT_LOG['norm_model'].append(np.linalg.norm(s_model))
        REPORT_LOG['norm_res'].append(res)
        ## visual res
        plt.plot(LOG_RES)
        plt.xlabel('k')
        plt.ylabel('Residual')
        plt.title("Residual over iterations")
        plt.savefig(os.path.join(dir_save, 'res.png'))
        plt.clf()
        ## visual model 
        v_model = np.divide(1.0, s_model)
        plt.imshow(
            v_model.reshape(settings.NY, settings.NX), cmap='jet',
            extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
            vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
        )
        a = plt.colorbar()
        a.set_label('$v$')
        plt.title("Model Velocity (LM, alpha={0:03f}, k={1})".format(alpha, K_STOP))
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.savefig(os.path.join(dir_save, 'model_v_alpha{}.png'.format(i_alpha)))
        plt.clf()
    # log l-curve
    f_curve = open(os.path.join("pb3","l_curve.log"), 'w')
    f_curve.write("alpha,norm_res,norm_model\n")
    for i in range(len(REPORT_LOG['norm_res'])):
        f_curve.write("{},{},{}\n".format(
            REPORT_LOG['alphas'],
            REPORT_LOG['norm_res'][i],
            REPORT_LOG['norm_model'][i])
        )
    f_curve.close()
    
    print('The process in LM is fully finish')