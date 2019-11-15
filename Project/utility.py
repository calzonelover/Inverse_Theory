import multiprocessing as mp
import numpy as np
import math

import settings



def deg_to_rad(deg):
    return deg * math.pi/180.0
def rad_to_deg(rad):
    return rad * 180.0/math.pi

'''
    Ray
'''

def get_source_receiver(is_separate=False):
    source_positions = []
    r = settings.NY*settings.DX
    theta_min_rad = math.acos((settings.PYRAMID_H/2.0 + settings.SIDE_GAP)/r)
    dtheta_rad = (math.pi - 2.0*theta_min_rad)/(settings.N_SOURCE)
    for i_s in range(settings.N_SOURCE):
        theta_source_rad = theta_min_rad + i_s * dtheta_rad + dtheta_rad/2.0
        source_positions.append({
            'x': settings.G_W/2.0 + r*math.cos(theta_source_rad),
            'y': r*math.sin(theta_source_rad)
        })
    receiver_positions = []
    _x_offset_inside = settings.SIDE_GAP + settings.PYRAMID_W/2.0 - settings.CHAMBER_W/2.0
    _x_step_inside = settings.CHAMBER_W/(settings.N_INSIDE_PYRAMID + 1.0)
    for i_r in range(settings.N_INSIDE_PYRAMID):
        receiver_positions.append({
            'x': _x_offset_inside + i_r*_x_step_inside,
            'y': 0.0
        })
    _x_offset_right = settings.PYRAMID_W + settings.SIDE_GAP
    _x_step_outside = settings.SIDE_GAP/(settings.N_OUTSIDE_PYRAMID_PER_SIDE + 1.0)
    for i_r in range(settings.N_OUTSIDE_PYRAMID_PER_SIDE):
        receiver_positions.append({
            'x': i_r*_x_step_outside + _x_step_outside/2.0,
            'y': 0.0
        })
        receiver_positions.append({
            'x': _x_offset_right + i_r*_x_step_outside + _x_step_outside/2.0,
            'y': 0.0
        })
    if is_separate:
        return {
            's_positions': source_positions,
            'r_positions': receiver_positions
        }
    sr_pairs = []
    for r_pos in receiver_positions:
        for s_pos in source_positions:
            sr_pairs.append({
                's_pos': s_pos,
                'r_pos': r_pos
            })
    return sr_pairs
    

def ray_length(x1, y1, x2, y2, mode='circle', is_fast_tracing=True):
    s_map = np.zeros(shape=(settings.NX*settings.NY))
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
    if mode=='circle':
        m_sr = (y2 - y1)/(x2 - x1)
        c_sr = y1 - m_sr * x1
        y_sr = lambda x: x * m_sr + c_sr
        a = 1.0 + m_sr*m_sr
        r = settings.DX/2.0
        for i_y in range(i_y_min ,i_y_max):
            for i_x in range(i_x_min ,i_x_max):
                x_0 = settings.DX * i_x + r
                y_0 = settings.DX * i_y + r
                # print(x_0, y_0)
                b = 2.0*m_sr*(c_sr - y_0) - 2.0*x_0
                c = x_0*x_0 - r*r + (c_sr - y_0)*(c_sr - y_0)
                in_sqrt = b*b - 4.0*a*c
                if in_sqrt > 0:
                    x_c_2 = (-b + math.sqrt(in_sqrt)) / (2.0*a)
                    x_c_1 = (-b - math.sqrt(in_sqrt)) / (2.0*a)
                    y_c_2 = y_sr(x_c_2)
                    y_c_1 = y_sr(x_c_1)

                    g_i = i_x + settings.NX * i_y
                    _dx = x_c_2 - x_c_1
                    _dy = y_c_2 - y_c_1
                    s_map[g_i] = math.sqrt(_dx*_dx + _dy*_dy)
    elif mode=='cube':
        r_s = np.array([x1, y1])
        r_r = np.array([x2, y2])
        D = np.subtract(r_r, r_s)
        length_D = np.linalg.norm(D)
        d = np.divide(D, length_D)
        d_angle = math.atan(d[1]/d[0])

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
    else:
        raise Exception('Mode {} ray tracing that you request is not available yet'.format(mode))
    return s_map