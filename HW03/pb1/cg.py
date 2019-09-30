import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray

# settings 
ALPHA0 = 1e-2
EPSILON = 2e6

REPORT_LOG = {
    'alphas': [],
    'norm_model': [],
    'norm_res': [],
}

def readraw(filename):
    f = open(settings.FILENAME, "r")
    return np.fromfile(f, dtype=np.float32)

def readmap(filename):
    a = readraw(filename)
    a = a.reshape(settings.NX, settings.NY)
    return a.T

def get_l():
    dy_sr = (settings.NY * settings.DX)/settings.N_SOURCE
    l = []
    for s_i in range(settings.N_SOURCE):
        for r_j in range(settings.N_RECEIVER):
            y_s_i = settings.DX/2 + s_i * dy_sr
            y_r_i = settings.DX/2 + r_j * dy_sr
            l_i = ray.ray_length(-10,y_s_i,1010,y_r_i)
            l.append(l_i)
    return np.array(l)

def grad(l, s ,t):
    return 0.5*np.matmul(
        (
            np.subtract(
                np.matmul(s.T, l),
                t
            )
        ).T,
        l
    )

if __name__ == "__main__":
    s_real = readraw(filename=settings.FILENAME)
    l = get_l()
    t_obs = np.matmul(s_real.T, l)

    s_model = np.ones(shape=(settings.NX * settings.NY,), dtype=float)
    err = 0.5*np.linalg.norm(
        np.subtract(
            np.matmul(s_model.T, l),
            t_obs
        )
    )
    alpha = ALPHA0
    pk_old = -0.5*np.subtract(
            np.matmul(s_model.T, l),
            t_obs
    )
    while err > EPSILON:
        gradk = grad(l, s_model, t_obs)
        betak = np.divide(
            np.matmul(
                    np.matmul(s_model.T, pk_old).T,
                    np.matmul(s_model.T, gradk)
            ),
            np.sum(np.square(np.matmul(s_model.T, pk_old))),
        )
        pk = -gradk + np.multiply(betak, pk_old)
        pk_norm = pk/np.linalg.norm(pk)
        alpha = -np.divide(
            np.sum(np.matmul(gradk, pk)),
            np.sum(np.square(np.matmul(s_model.T, pk))),
        )

        s_model += np.multiply(alpha, pk_norm)
        err = 0.5*np.linalg.norm(
            np.subtract(
                np.matmul(s_model.T, l),
                t_obs
            )
        )
        pk_old = pk
        print(err)

    # visualize model
    map_model = s_model.reshape(settings.NX, settings.NY).T
    plt.imshow(map_model, cmap='jet', extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (SD)")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    # plt.savefig("img/model_v_alpha{}.png".format(alpha))
    plt.show()