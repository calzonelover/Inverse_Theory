# Reference
# https://www.doc.ic.ac.uk/~dfg/graphics/graphics2008/GraphicsLecture09.pdf
# web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/basic_algo.pdf

import math
import numpy as np
import settings

def ray_length(x1, y1, x2, y2):
    r_s = np.array([x1, y1])
    r_r = np.array([x2, y2])
    D = np.subtract(r_r, r_s)
    d = np.divide(D, np.linalg.norm(D))
    d_angle = math.atan(d[1]/d[0])

    s_map = np.zeros(shape=(settings.NX*settings.NY))
    for i_y in range(settings.NY):
        for i_x in range(settings.NX):
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
                s_map[g_i] = t[2] - t[1]
    return s_map