import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray

# settings 
EPSILON = 8.5e-4

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
    return np.array(l).T

def grad(l, s ,t):
    return 0.5*np.matmul(
        np.subtract(
            np.matmul(l, s),
            t
        ).T,
        l
    )

if __name__ == "__main__":
    v_real = readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    l = get_l()
    t_obs = np.add(
        np.matmul(l, s_real),
        np.random.normal(loc=0.0, scale=1e-4, size=(settings.N_SOURCE*settings.N_RECEIVER))
    )

    s_model = 1e-3*np.ones(shape=s_real.shape)# np.divide(np.random.normal(size=s_real.shape), 1700)
    # s_model = np.matmul(
    #     np.linalg.inv(
    #         np.add(
    #             np.matmul(l.T, l),
    #             np.multiply(7e-2, np.identity(settings.NX*settings.NY, dtype=float))
    #         )
    #     ),
    #     np.matmul(l.T, t_obs)
    # )
    err = 0.5*np.linalg.norm(
        np.subtract(
            np.matmul(l, s_model),
            t_obs
        )
    )
    while err > EPSILON:
        pk = -grad(l, s_model, t_obs)
        alpha = np.divide(
            np.sum(np.square(pk)),
            np.sum(np.square(np.matmul(l, pk))),
        )

        s_model += np.multiply(alpha, pk)
        new_err = 0.5*np.linalg.norm(
            np.subtract(
                np.matmul(l, s_model),
                t_obs
            )
        )
        err = new_err
        print(err)

    # visualize model
    map_model = 1.0/s_model.reshape(settings.NX, settings.NY).T
    plt.imshow(map_model, cmap='jet', extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (SD)")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    # plt.savefig("img/model_v_alpha{}.png".format(alpha))
    plt.show()