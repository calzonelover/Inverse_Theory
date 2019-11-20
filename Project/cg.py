import numpy as np
import os
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
else:
    import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import environment, utility, settings


def main():
    # settings 
    EPSILON = 1e-3
    MAX_ITER = 100
    ALPHA0 = 1.0
    ALPHA_DECAYRATE = 0.8

    real_lambda = environment.get_system(model='real')
    ray_paths = utility.get_ray_paths()
    zenith_angles = utility.get_source_receiver(is_separate=True)['zenith_angles']
    I_obs = utility.get_flux(zenith_angles, ray_paths, real_lambda)

    model0_lambda = environment.get_system(model='initial')
    model_lambda = environment.get_system(model='initial')
    err = utility.get_r(zenith_angles, ray_paths, model_lambda, I_obs)
    k = 0
    LOG_ERRS = []
    try:
        pk = utility.gradient(zenith_angles, ray_paths, model_lambda, I_obs)
        pk = np.negative(np.divide(pk, np.linalg.norm(pk)))
        if settings.SMOOTH_GRAD:
            pk = utility.smooth_map(pk)
        if settings.FILTER_GRAD:
            pk = utility.gradient_filter(pk, model0_lambda)
        while err > EPSILON and k < MAX_ITER:
            alphak = utility.get_proper_alpha(
                zenith_angles, ray_paths, model_lambda, I_obs, pk,
                ALPHA0=ALPHA0, ALPHA_DECAYRATE=ALPHA_DECAYRATE
            )
            model_lambda_new = utility.prevent_negative(np.add(model_lambda, np.multiply(alphak, pk)))
            gradk = utility.gradient(zenith_angles, ray_paths, model_lambda, I_obs)
            gradk1 = utility.gradient(zenith_angles, ray_paths, model_lambda_new, I_obs)
            betak1 = np.matmul(gradk1.T, gradk1)/np.matmul(gradk.T, gradk)

            pk = np.add(np.negative(gradk1), np.multiply(betak1, pk))
            pk = np.divide(pk, np.linalg.norm(pk))
            if settings.SMOOTH_GRAD:
                pk = utility.smooth_map(pk)
            if settings.FILTER_GRAD:
                pk = utility.gradient_filter(pk, model0_lambda)
            
            err = utility.get_r(zenith_angles, ray_paths, model_lambda, I_obs)
            LOG_ERRS.append(err)
            print(k, err)

            model_lambda = model_lambda_new
            k += 1
    except KeyboardInterrupt:
        pass

    # visualize
    file_name = 'k{}_filter{}_smooth{}'.format(
        k,
        1 if settings.FILTER_GRAD else 0,
        1 if settings.SMOOTH_GRAD else 0,
    )

    plt.plot(LOG_ERRS)
    plt.xlabel('k')
    plt.ylabel('Residual')
    plt.title("Residual over iterations (CG, k=%d, $\\alpha_0$=%.2f and $\\alpha_d$=%.2f)"%(k, ALPHA0, ALPHA_DECAYRATE))
    plt.savefig(os.path.join('cg', 'res_%s.png'%file_name))
    plt.clf()

    plt.imshow(
        model_lambda.reshape(settings.NY, settings.NX), cmap='summer',
        extent=[0, settings.DX*settings.NX, 0, settings.DX*settings.NY],
        origin='bottom left',
    )
    a = plt.colorbar()
    a.set_label('1 / Path length ($m^{-1}$)')
    plt.title("Result from CG (k=%d, $\\alpha_0$=%.2f and $\\alpha_d$=%.2f)"%(k, ALPHA0, ALPHA_DECAYRATE))
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    # plt.show()
    plt.savefig(os.path.join('cg', 'model_%s.png'%file_name))
    plt.clf()