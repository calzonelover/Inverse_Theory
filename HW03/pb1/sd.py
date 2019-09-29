import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray

ALPHA = 0.1

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

if __name__ == "__main__":
    s_real = readraw(filename=settings.FILENAME)

    dy_sr = (settings.NY * settings.DX)/settings.N_SOURCE
    t = []
    l = []
    for s_i in range(settings.N_SOURCE):
        for r_j in range(settings.N_RECEIVER):
            y_s_i = settings.DX/2 + s_i * dy_sr
            y_r_i = settings.DX/2 + r_j * dy_sr

            l_i = ray.ray_length(-10,y_s_i,1010,y_r_i)
            t_i = np.add(
                np.matmul(s_real.T, l_i),
                np.random.normal(loc=0, scale=1e-4)
            )

            l.append(l_i)
            t.append(t_i)
    t = np.array(t)
    l = np.array(l)
    
    s_model = np.matmul(
        np.linalg.inv(
            np.add(
                np.matmul(l.T, l),
                np.multiply(alpha, np.identity(settings.NX*settings.NY, dtype=float))
            )
        ),
        np.matmul(l.T, t)
    )
    # visualize model
    map_model = s_model.reshape(settings.NX, settings.NY).T
    plt.imshow(map_model, cmap='jet', extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Model Velocity (alpha = {:6.6f})".format(alpha))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    # plt.savefig("img/model_v_alpha{}.png".format(alpha))
    plt.show()