# Reference
# https://www.doc.ic.ac.uk/~dfg/graphics/graphics2008/GraphicsLecture09.pdf
# web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/basic_algo.pdf

import math
import numpy as np
import settings

import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def ray_length(x1, y1, x2, y2, is_fast_tracing=True):
    r_s = np.array([x1, y1])
    r_r = np.array([x2, y2])
    D = np.subtract(r_r, r_s)
    length_D = np.linalg.norm(D)
    d = np.divide(D, length_D)
    d_angle = math.atan(d[1]/d[0])

    # reduce calculation time by scoping the possible square mesh
    if is_fast_tracing:
        i_x_min = math.floor(min(x1, x2)/settings.DX) - 1
        if i_x_min < 0: i_x_min = 0
        i_x_max = math.floor(max(x1, x2)/settings.DX) + 2
        if i_x_max > settings.NX-1: i_x_max = settings.NX
        i_y_min = math.floor(min(y1, y2)/settings.DX) - 1
        if i_y_min < 0: i_y_min = 0
        i_y_max = math.floor(max(y1, y2)/settings.DX) + 2
        if i_y_max > settings.NX-1: i_y_max = settings.NX
    else:
        i_x_min, i_x_max = 0, settings.NX
        i_y_min, i_y_max = 0, settings.NY

    s_map = np.zeros(shape=(settings.NX*settings.NY))
    for i_y in range(i_y_min ,i_y_max):
        for i_x in range(i_x_min ,i_x_max):
            g_i = i_x + settings.NX * i_y

            x_min = settings.DX * i_x
            x_max = settings.DX * (i_x + 1)
            y_min = settings.DX * i_y
            y_max = settings.DX * (i_y + 1)

            x_min_r = x_min - r_s[0]
            x_max_r = x_max - r_s[0]
            y_min_r = y_min - r_s[1]
            y_max_r = y_max - r_s[1]

            angles = np.array([
                math.atan(y_max_r/x_min_r), math.atan(y_max_r/x_max_r),
                math.atan(y_min_r/x_min_r), math.atan(y_min_r/x_max_r)
            ])
            if d_angle >= min(angles) and d_angle <= max(angles):
                t = np.array([
                    (x_min - r_s[0])/d[0], (x_max - r_s[0])/d[0],
                    (y_min - r_s[1])/d[1], (y_max - r_s[1])/d[1],
                ])
                t.sort()
                # same block
                if (r_s[0] > x_min and r_s[0] < x_max and r_s[1] > y_min and r_s[1] < y_max
                    and r_r[0] > x_min and r_r[0] < x_max and r_r[1] > y_min and r_r[1] < y_max ):
                # if t[0] < 0 and t[1] < 0 and t[2] > 0 and t[3] > 0:
                    s_map[g_i] = length_D
                # Not the same block
                elif t[1] > 0 and t[2] > 0 and t[3] > 0:
                    # far
                    if t[0] > 0 and t[1] < length_D and t[2] > length_D:
                        s_map[g_i] = t[2] - t[1]
                    # near
                    elif t[0] < 0 and t[2] < 2*settings.DX:
                        if length_D < t[2] and length_D - t[1] > 0:
                            s_map[g_i] = length_D - t[1]
                        else:
                            s_map[g_i] = t[2] - t[1]
    return s_map


def grad_mesh(t):
    gx = t[1,2] - t[1,0]
    gy = t[2,1] - t[0,1]
    norm = math.sqrt(gx*gx + gy*gy)
    grad = np.divide(np.array([gx, gy]), norm)
    return grad

def curved_ray(source_x, source_y, receiver_x, receiver_y, T):
    L = np.zeros(settings.NX*settings.NY)
    t = T.reshape(settings.NY, settings.NX)
    t = np.pad(t, (1,1), 'maximum')

    x_now, y_now = receiver_x, receiver_y
    distance_to_receiver = math.sqrt((source_x-x_now)*(source_x-x_now)+ (source_y-y_now)*(source_y-y_now))
    while distance_to_receiver > 1.5*settings.DX:
        i_receiver_x = math.floor(x_now/settings.DX) + 1
        i_receiver_y = math.floor(y_now/settings.DX) + 1
        grad_T = grad_mesh(t[i_receiver_y-1:i_receiver_y+2, i_receiver_x-1:i_receiver_x+2])

        # L = np.add(L, ray_length(2000, 2000, 6700, 1000))
        L = np.add(L, ray_length(x_now, y_now, x_now - settings.DX * grad_T[0], y_now - settings.DX * grad_T[1]))
        
        x_now = x_now - settings.DX * grad_T[0]
        y_now = y_now - settings.DX * grad_T[1]
        distance_to_receiver = math.sqrt((source_x-x_now)*(source_x-x_now)+ (source_y-y_now)*(source_y-y_now))
        # print(distance_to_receiver, x_now, y_now)

        # plt.imshow(
        #     L.reshape(settings.NY, settings.NX), cmap='jet',
        #     extent=[0, settings.DX*settings.NX, settings.DX*settings.NY, 0],
        # )
        # a = plt.colorbar()
        # a.set_label('$T$')
        # plt.title("Travel time")
        # plt.xlabel("$x$")
        # plt.ylabel("$y$")
        # plt.show()
        # plt.clf()
    return L