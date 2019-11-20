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
    MAX_ITER = 400
    ALPHA0 = 0.8
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
        while err > EPSILON and k < MAX_ITER:
            pk = utility.gradient(zenith_angles, ray_paths, model_lambda, I_obs)
            pk = np.negative(np.divide(pk, np.linalg.norm(pk)))
            if settings.SMOOTH_GRAD:
                pk = utility.smooth_map(pk)
            if settings.FILTER_GRAD:
                pk = utility.gradient_filter(pk, model0_lambda)

            alpha = utility.get_proper_alpha(
                zenith_angles, ray_paths, model_lambda, I_obs, pk,
                ALPHA0=ALPHA0, ALPHA_DECAYRATE=ALPHA_DECAYRATE
            )
            model_lambda = utility.prevent_negative(np.add(model_lambda, np.multiply(alpha, pk)))

            err = utility.get_r(zenith_angles, ray_paths, model_lambda, I_obs)
            LOG_ERRS.append(err)
            print(k, err)
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
    plt.title("Residual over iterations (SD, k=%d, $\\alpha_0$=%.2f and $\\alpha_d$=%.2f)"%(k, ALPHA0, ALPHA_DECAYRATE))
    plt.savefig(os.path.join('sd', 'res_%s.png'%file_name))
    plt.clf()

    plt.imshow(
        model_lambda.reshape(settings.NY, settings.NX), cmap='summer',
        extent=[0, settings.DX*settings.NX, 0, settings.DX*settings.NY],
        origin='bottom left',
    )
    a = plt.colorbar()
    a.set_label('1 / Path length ($m^{-1}$)')
    plt.title("Result from SD (k=%d, $\\alpha_0$=%.2f and $\\alpha_d$=%.2f)"%(k, ALPHA0, ALPHA_DECAYRATE))
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    # plt.show()
    plt.savefig(os.path.join('sd', 'model_%s.png'%file_name))
    plt.clf()