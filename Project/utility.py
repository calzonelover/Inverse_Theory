from scipy.ndimage import gaussian_filter, uniform_filter
import multiprocessing as mp
import numpy as np
import math

import settings

'''
    Transformation
'''

def deg_to_rad(deg):
    return deg * math.pi/180.0
def rad_to_deg(rad):
    return rad * 180.0/math.pi

def smooth_map(map, kernel_size=(2, 2), mode='gaussian'):
    if mode == 'uniform':
        out = uniform_filter(map.reshape(settings.NY, settings.NX), size=kernel_size).reshape(settings.NX*settings.NY)
    elif mode == 'gaussian':
        out = gaussian_filter(map.reshape(settings.NY, settings.NX), sigma=kernel_size).reshape(settings.NX*settings.NY)
    return out

'''
    Governing equation
'''
def get_flux(zenith_angles, ray_paths, lamdba):
    return np.multiply(
        np.square(np.cos(zenith_angles)),
        np.exp(np.negative(
            np.matmul(ray_paths, lamdba)
        ))
    )

'''
    Optimization
'''

def prevent_negative(_model):
    _func = np.vectorize(lambda x: x if x > 0.0 else 0.0)
    return _func(_model)

def get_r(zenith_angles, ray_paths, lamdba, I_obs, magnitude=True):
    I_lambda = get_flux(zenith_angles, ray_paths, lamdba)
    r = np.subtract(I_obs, I_lambda)
    if magnitude:
        return np.linalg.norm(r)
    return r

def gradient(zenith_angles, ray_paths, lamdba, I_obs):
    I_lamdba = get_flux(zenith_angles, ray_paths, lamdba)
    return np.matmul(
        np.multiply(
            np.subtract(I_obs, I_lamdba),
            I_lamdba
        ).T,
        ray_paths
    )

def get_pk(zenith_angles, ray_paths, model_lambda, I_obs, model0_lambda):
    pk = gradient(zenith_angles, ray_paths, model_lambda, I_obs)
    pk = np.negative(np.divide(pk, np.linalg.norm(pk)))
    if settings.SMOOTH_GRAD:
        pk = smooth_map(pk)
    if settings.FILTER_GRAD:
        pk = gradient_filter(pk, model0_lambda)
    return pk

def gradient_filter(gradient, initial_lambda):
    _func = np.vectorize(lambda x: 1 if x == settings.LAMBDA_ROCK else 0)
    _kernel = _func(initial_lambda)
    return np.multiply(gradient, _kernel)

def get_proper_alpha(zenith_angles, ray_paths, lambda0, I_obs, pk, method='backtrack', ALPHA0=0.8, ALPHA_DECAYRATE=0.8):
    C = 0.1
    if method == 'backtrack':
        alphak = ALPHA0
        while True:
            # lambda1 = prevent_negative(np.add(lambda0, np.multiply(alphak, pk)))
            lambda1 = np.add(lambda0, np.multiply(alphak, pk))
            gradk = gradient(zenith_angles, ray_paths, lambda1, I_obs)
            if get_r(zenith_angles, ray_paths, lambda1, I_obs) <= get_r(zenith_angles, ray_paths, lambda0, I_obs) + np.multiply(C * alphak, np.matmul(gradk.T, pk)):
                break
            alphak *= ALPHA_DECAYRATE
            # print(alphak)
    else:
        print('Method {} to find the step length does not support yet'.format(method))
        exit()
    return alphak

'''
    Ray
'''

def get_source_receiver(is_separate=False):
    source_positions = []
    r = settings.NY*settings.DX
    theta_min_rad = 0.0 # math.acos((settings.PYRAMID_H/2.0 + settings.SIDE_GAP)/r)
    dtheta_rad = math.pi/settings.N_SOURCE# (math.pi - 2.0*theta_min_rad)/(settings.N_SOURCE)
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
    zenith_angles = [
        math.atan(abs((s_pos['x']-r_pos['x'])/(s_pos['y']-r_pos['y'])))
        for r_pos in receiver_positions
        for s_pos in source_positions
    ]
    if is_separate:
        return {
            's_positions': source_positions,
            'r_positions': receiver_positions,
            'zenith_angles': zenith_angles
        }
    sr_pairs = []
    for r_pos in receiver_positions:
        for s_pos in source_positions:
            sr_pairs.append({
                's_pos': s_pos,
                'r_pos': r_pos,
                'zenith_angle': math.atan(abs((s_pos['x']-r_pos['x'])/(s_pos['y']-r_pos['y'])))
            })
    return sr_pairs
    

def ray_length(x1, y1, x2, y2, mode='circle', is_fast_tracing=True):
    s_map = np.zeros(shape=(settings.NX*settings.NY))
    # reduce calculation time by scoping the possible square mesh
    if is_fast_tracing:
        i_x_min = math.floor(min(x1, x2)/settings.DX) - 2
        if i_x_min < 0: i_x_min = 0
        i_x_max = math.floor(max(x1, x2)/settings.DX) + 2
        if i_x_max >= settings.NX-1: i_x_max = settings.NX
        i_y_min = math.floor(min(y1, y2)/settings.DX) - 2
        if i_y_min < 0: i_y_min = 0
        i_y_max = math.floor(max(y1, y2)/settings.DX) + 2
        if i_y_max >= settings.NY-1: i_y_max = settings.NY
    else:
        i_x_min, i_x_max = 0, settings.NX
        i_y_min, i_y_max = 0, settings.NY
    if mode=='circle':
        # hor 
        if x1 == x2:
            i_x = math.floor(x1/settings.DX)
            for i_y in range(i_y_min ,i_y_max):
                y_0 = settings.DX * i_y + settings.DX/2.0
                if y_0 > min(y1, y2) and y_0 < max(y1, y2):
                    g_i = i_x + settings.NX * i_y
                    s_map[g_i] = settings.DX
        elif y1 == y2:
            i_y = math.floor(y1/settings.DX)
            for i_x in range(i_x_min ,i_x_max):
                x_0 = settings.DX * i_x + settings.DX/2.0
                if x_0 > min(x1, x2) and x_0 < max(x1, x2):
                    g_i = i_x + settings.NX * i_y
                    s_map[g_i] = settings.DX
        else:
            m_sr = (y2 - y1)/(x2 - x1)
            c_sr = y1 - m_sr * x1
            y_sr = lambda x: x * m_sr + c_sr
            a = 1.0 + m_sr*m_sr
            r = settings.DX/2.0
            for i_y in range(i_y_min ,i_y_max):
                for i_x in range(i_x_min ,i_x_max):
                    g_i = i_x + settings.NX * i_y
                    x_0 = settings.DX * i_x + r
                    y_0 = settings.DX * i_y + r
                    # if y_0 > min(y1, y2) and y_0 < max(y1, y2) and x_0 > min(x1, x2) and x_0 < max(x1, x2):
                    b = 2.0*m_sr*(c_sr - y_0) - 2.0*x_0
                    c = x_0*x_0 - r*r + (c_sr - y_0)*(c_sr - y_0)
                    in_sqrt = b*b - 4.0*a*c
                    if in_sqrt > 0:
                        x_c_2 = (-b + math.sqrt(in_sqrt)) / (2.0*a)
                        x_c_1 = (-b - math.sqrt(in_sqrt)) / (2.0*a)
                        y_c_2 = y_sr(x_c_2)
                        y_c_1 = y_sr(x_c_1)
                        if (x_c_2 >= min(x1,x2) and x_c_2 <= max(x1, x2) and x_c_1 >= min(x1,x2) and x_c_1 <= max(x1, x2)
                            and y_c_2 >= min(y1,y2) and y_c_2 <= max(y1, y2) and y_c_1 >= min(y1,y2) and y_c_1 <= max(y1, y2)):
                            _dx = x_c_2 - x_c_1
                            _dy = y_c_2 - y_c_1
                            s_map[g_i] = math.sqrt(_dx*_dx + _dy*_dy)
                    else:
                        x_min = settings.DX * i_x
                        x_max = settings.DX * (i_x + 1)
                        y_min = settings.DX * i_y
                        y_max = settings.DX * (i_y + 1)
                        if x1 > x_min and x1 < x_max and y1 > y_min and y1 < y_max:
                            r_s = np.array([x1, y1])
                            r_r = np.array([x2, y2])
                            D = np.subtract(r_r, r_s)
                            length_D = np.linalg.norm(D)
                            d = np.divide(D, length_D)
                            d_angle = math.atan(d[1]/d[0])
                            if x2 > x_min and x2 < x_max and y2 > y_min and y2 < y_max:
                                s_map[g_i] = length_D
                            else:
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
                                    dr_min = min([
                                        math.sqrt((x_c_1-x1)**2 + (y_c_1-y1)**2),
                                        math.sqrt((x_c_2-x1)**2 + (y_c_2-y1)**2)
                                    ])
                                    s_map[g_i] = math.sqrt(_dx*_dx + _dy*_dy) - dr_min
    elif mode=='cube':              
        ## old
        r_s = np.array([x1, y1])
        r_r = np.array([x2, y2])
        D = np.subtract(r_r, r_s)
        length_D = np.linalg.norm(D)
        d = np.divide(D, length_D)
        d_angle = math.atan2(d[1], d[0])

        for i_y in range(i_y_min ,i_y_max):
            for i_x in range(i_x_min ,i_x_max):
                g_i = i_x + settings.NX * i_y
                x_min = settings.DX * i_x
                x_max = settings.DX * (i_x + 1)
                y_min = settings.DX * i_y
                y_max = settings.DX * (i_y + 1)
                ## same block
                if (r_s[0] > x_min and r_s[0] < x_max and r_s[1] > y_min and r_s[1] < y_max
                    and r_r[0] > x_min and r_r[0] < x_max and r_r[1] > y_min and r_r[1] < y_max ):
                    s_map[g_i] = length_D
                ## not the same block
                else:
                    t_xmin = (x_min - r_s[0])/d[0]
                    t_xmax = (x_max - r_s[0])/d[0]
                    t_ymin = (y_min - r_s[1])/d[1]
                    t_ymax = (y_max - r_s[1])/d[1]
                    '''
                    # old
                    x_min_r = x_min - r_s[0]
                    x_max_r = x_max - r_s[0]
                    y_min_r = y_min - r_s[1]
                    y_max_r = y_max - r_s[1]                    
                    angles = np.array([
                        math.atan2(y_max_r, x_min_r), math.atan2(y_max_r, x_max_r),
                        math.atan2(y_min_r, x_min_r), math.atan2(y_min_r, x_max_r)
                    ])
                    t = [
                        t_xmin, t_xmax,
                        t_ymin, t_ymax,
                    ]
                    t.sort()
                    if d_angle <= max(angles) and d_angle >= min(angles):
                        ## Not the same block
                        # far
                        if t[0] > 0.0 and t[3] < length_D:
                            s_map[g_i] = t[2]  - t[1]
                        # stencil
                        elif t[0] < 0.0 and t[1] > 0.0 and t[1] < length_D and t[2] > length_D and length_D < 2.5*settings.DX:
                            s_map[g_i] = length_D  - t[1]
                    '''
                    # new
                    t_near, t_far = -1e10, 1e10
                    if t_ymin > t_ymax:
                        t_ymin, t_ymax = t_ymax, t_ymin
                    if t_xmin > t_xmax:
                        t_xmin, t_xmax = t_xmax, t_xmin
                    if max(t_xmin, t_ymin) > t_near:
                        t_near = max(t_xmin, t_ymin)
                    if min(t_xmax, t_ymax) < t_far:
                        t_far = min(t_xmax, t_ymax)
                    if t_near < t_far and t_far > 0:
                        s_map[g_i] = t_far - t_near
    else:
        raise Exception('Mode {} ray tracing that you request is not available yet'.format(mode))
    return s_map


'''
    instant ray
'''

def get_ray_paths():
    pair_srs = get_source_receiver()
    _l = []
    for pair_sr in pair_srs:
        s_pos = pair_sr['s_pos']
        r_pos = pair_sr['r_pos']
        _l.append(
            ray_length(
                s_pos['x'],s_pos['y'],r_pos['x'],r_pos['y'],
                mode=settings.TRACING_MODE
            )
        )
    return np.array(_l)