import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray

# settings 
EPSILON = 8e-4

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

if __name__ == "__main__":
    v_real = readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    l = get_l()
    t_obs = np.add(
        np.matmul(l, s_real),
        np.random.normal(loc=0.0, scale=1e-4, size=(settings.N_SOURCE*settings.N_RECEIVER))
    )

    sk = 1e-3*np.ones(shape=s_real.shape)
    # sk = np.matmul(
    #     np.linalg.inv(
    #         np.add(
    #             np.matmul(l.T, l),
    #             np.multiply(7e-2, np.identity(settings.NX*settings.NY, dtype=float))
    #         )
    #     ),
    #     np.matmul(l.T, t_obs)
    # )
    rk = np.subtract(
        np.matmul(l, sk),
        t_obs,
    )
    pk = -rk
    err = 0.5*np.linalg.norm(rk)
    while err > EPSILON:
        alphak = np.divide(
            np.matmul(
                rk.T,
                rk
            ),
            np.matmul(
                pk.T,
                np.matmul(
                    l,
                    pk
                )
            )
        )
        sk1 = np.add(sk, np.multiply(alphak, pk))
        rk1 = np.add(rk, np.multiply(alphak,
            np.matmul(l, pk)
        ))
        betak1 = np.divide(
            np.matmul(
                rk1.T,
                rk1
            ),
            np.matmul(
                rk.T,
                rk
            ),
        )
        pk1 = np.add(-rk1, np.multiply(betak1, pk))

        sk = sk1
        rk = rk1
        pk = pk1
        err = np.linalg.norm(rk)
        print(err)
    s_model = sk
    # visualize model
    map_model = 1.0/s_model.reshape(settings.NX, settings.NY).T
    plt.imshow(map_model, cmap='jet', extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (CG)")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    # plt.savefig("img/model_v_alpha{}.png".format(alpha))
    plt.show()