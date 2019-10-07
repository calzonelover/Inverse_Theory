import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings, ray

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

if __name__ == "__main__":
    v_real = readraw(filename=settings.FILENAME)
    s_real = 1.0/v_real
    l = get_l()
    t_obs = np.add(
        np.matmul(l, s_real),
        np.random.normal(loc=0.0, scale=1e-4, size=(settings.N_SOURCE*settings.N_RECEIVER))
    )

    for i, alpha in enumerate(settings.ALPHAS):
        print(alpha)
        s_model = np.matmul(
            np.linalg.inv(
                np.add(
                    np.matmul(l.T, l),
                    np.multiply(alpha, np.identity(settings.NX*settings.NY, dtype=float))
                )
            ),
            np.matmul(l.T, t_obs)
        )
        REPORT_LOG['alphas'].append(alpha)
        REPORT_LOG['norm_model'].append(np.sqrt(np.sum(np.square(s_model))))
        REPORT_LOG['norm_res'].append(np.linalg.norm(
            np.subtract(
                t_obs,
                np.matmul(l, s_model)
        )))
        # visualize model
        map_model = 1.0/s_model.reshape(settings.NX, settings.NY).T
        plt.imshow(
            map_model, cmap='jet',
            extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
            vmin=settings.COLOR_VMIN, vmax=settings.COLOR_VMAX,    
        )
        a = plt.colorbar()
        a.set_label('$v$')
        plt.title("Model Velocity (alpha = {:010.05f})".format(alpha))
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.savefig("img/model_v_alpha{:03d}.png".format(i))
        plt.clf()
        # plt.show()

    # log l-curve
    f_curve = open("l_curve.log", 'w')
    f_curve.write("alpha,norm_res,norm_model\n")
    for i in range(len(REPORT_LOG['norm_res'])):
        f_curve.write("{},{},{}\n".format(
            settings.ALPHAS[i],
            REPORT_LOG['norm_res'][i],
            REPORT_LOG['norm_model'][i])
        )
    f_curve.close()

    # visualize raw
    raw_map = readmap(filename=settings.FILENAME)
    plt.imshow(
        raw_map, cmap='jet',
        extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
    )
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Real Velocity")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.savefig("real_v.png")
    # plt.show()

    '''
    # Test ray tracing
    for i_y in range(40):
        y = 10 + settings.DX * i_y
        one_ray_map = ray.ray_length(-10,510,1010,y).reshape(50, 50)
        plt.imshow(one_ray_map, extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0])
        a = plt.colorbar()
        a.set_label('Ray length (m)')
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title("Ray Tracing (source: ({},{}), receiver({},{}))".format(-10, 510, 1010, y))
        plt.savefig("img/ray_{:02.0f}.png".format(i_y))
        plt.clf()
    '''