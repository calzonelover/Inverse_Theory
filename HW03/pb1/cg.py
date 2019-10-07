import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray

# settings
EPSILON = 1e-3
LOG_ERRS = []

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

    mk = 1e-3*np.zeros(shape=s_real.shape)
    pk0 = np.zeros(shape=mk.shape)
    sk = np.subtract(
        t_obs,
        np.matmul(l, mk)
    )
    rk = np.matmul(l.T, sk)
    err = 0.5*np.linalg.norm(np.subtract(
        t_obs,
        np.matmul(l, mk)
    ))
    while err > EPSILON:
        try:
            betak = np.divide(
                np.matmul(rk.T, rk),
                np.matmul(rk0.T, rk0)
            )
        except NameError:
            betak = 0.0
        pk = np.add(
            rk,
            np.multiply(betak, pk0)
        )
        alphak = np.divide(
            np.linalg.norm(rk)**2,
            np.matmul(
                np.matmul(l, pk).T,
                np.matmul(l, pk)
            )
        )
        mk1 = np.add(
            mk,
            np.multiply(alphak, pk)
        )
        sk1 = np.subtract(
            sk,
            np.multiply(
                alphak,
                np.matmul(l, pk)
            )
        )
        rk1 = np.matmul(l.T, sk1)
        
        err = 0.5*np.linalg.norm(np.subtract(
            np.matmul(l, mk),
            t_obs,
        ))
        LOG_ERRS.append(err)
        print(err)
        pk0 = pk
        mk0 = mk
        mk = mk1
        sk0 = sk
        sk = sk1
        rk0 = rk
        rk = rk1

    s_model = mk
    # visualize model
    map_model = 1.0/s_model.reshape(settings.NX, settings.NY).T
    plt.imshow(
        map_model, cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (CGLS, $\epsilon$~{})".format(EPSILON))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig("v_cg.png")
    # plt.show()
    plt.clf()
    plt.plot(LOG_ERRS)
    plt.xlabel("k")
    plt.ylabel("$||t_{obs}-Ls_{model}||$")
    plt.yscale("log")
    plt.title("Decay of residue (CGLS)")
    plt.savefig("cg_r.png")