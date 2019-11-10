import math
import numpy as np
import platform
from scipy.ndimage import gaussian_filter, uniform_filter
import multiprocessing as mp

import settings, ray

def readraw(filename):
    f = open(settings.FILENAME, "r")
    return np.fromfile(f, dtype=np.float32)

def readmap(filename):
    a = readraw(filename)
    a = a.reshape(settings.NX, settings.NY)
    return a.T

'''
    Fast Sweeping Methods to find the 
'''
def get_mean_t(i_y, i_x, s, old_T):
    T_mean = None
    ### bc
    if i_x > 0 and i_x < settings.NX-1 and i_y > 0 and i_y < settings.NY-1:
        T_x_min = min(old_T[i_y*settings.NX+(i_x-1)], old_T[i_y*settings.NX+(i_x+1)])
        T_y_min = min(old_T[(i_y-1)*settings.NX+i_x], old_T[(i_y+1)*settings.NX+i_x])

    elif i_x == 0 and i_y > 0 and i_y < settings.NY-1: # left
        T_x_min = old_T[i_y*settings.NX + i_x+1]
        T_y_min = min(old_T[(i_y-1)*settings.NX+i_x], old_T[(i_y+1)*settings.NX+i_x])
    elif i_x == settings.NX-1 and i_y > 0 and i_y < settings.NY-1: # right
        T_x_min = old_T[i_y*settings.NX + i_x-1]
        T_y_min = min(old_T[(i_y-1)*settings.NX+i_x], old_T[(i_y+1)*settings.NX+i_x])

    elif i_y == 0 and i_x > 0 and i_x < settings.NX-1: # bottom
        T_x_min = min(old_T[i_y*settings.NX+(i_x-1)], old_T[i_y*settings.NX+(i_x+1)])
        T_y_min = old_T[(i_y+1)*settings.NX+i_x]
    elif i_y == settings.NY-1 and i_x > 0 and i_x < settings.NX-1: # top
        T_x_min = min(old_T[i_y*settings.NX+(i_x-1)], old_T[i_y*settings.NX+(i_x+1)])
        T_y_min = old_T[(i_y-1)*settings.NX+i_x]

    elif i_x == 0 and i_y == 0: # lower-left corner:
        T_x_min = old_T[i_y*settings.NX+i_x+1]
        T_y_min = old_T[(i_y+1)*settings.NX+i_x]
    elif i_x == settings.NX-1 and i_y == 0: # lower-right corner:
        T_x_min = old_T[i_y*settings.NX+i_x-1]
        T_y_min = old_T[(i_y+1)*settings.NX+i_x]
    elif i_x == 0 and i_y == settings.NY-1: # upper-left corner:
        T_x_min = old_T[i_y*settings.NX+i_x+1]
        T_y_min = old_T[(i_y-1)*settings.NX+i_x]
    elif i_x == settings.NX-1 and i_y == settings.NY-1: # upper-right corner:
        T_x_min = old_T[i_y*settings.NX+i_x-1]
        T_y_min = old_T[(i_y-1)*settings.NX+i_x]
    else:
        raise IndexError('Travel time indices i_y:{}, i_x:{} could not satisfy any boundary condition'.format(i_y, i_x))
        
    if abs(T_x_min - T_y_min) > s[i_y*settings.NX + i_x]*settings.DX:
        T_mean = min(T_x_min, T_y_min) + s[i_y*settings.NX + i_x]*settings.DX
    else:
        T_mean = 0.5*(
            T_x_min + T_y_min
            + math.sqrt(
                2.0*s[i_y*settings.NX + i_x]*s[i_y*settings.NX + i_x]*settings.DX*settings.DX
                - (T_x_min - T_y_min)*(T_x_min - T_y_min)
            )
        )
    return T_mean

def get_t_new(i_y, i_x, s, T):
    Tnew = get_mean_t(i_y, i_x, s, T)
    return min(T[i_y*settings.NX + i_x], Tnew)

def get_travel_time(s, source_x, source_y, n_sweep=settings.N_SWEEP):
    T = np.multiply(1e6, np.ones(settings.NY*settings.NX))
    i_source_x = math.floor(source_x/settings.DX)
    i_source_y = math.floor(source_y/settings.DX)
    T[i_source_y*settings.NX + i_source_x] = 0.0
    for i_sweep in range(n_sweep):
        old_T = T
        # quad 1
        # if not settings.N_PROCESSES:
        for i_y in range(settings.NY):
            for i_x in range(settings.NX):
                Tnew = get_mean_t(i_y, i_x, s, T)
                T[i_y*settings.NX + i_x] = min(T[i_y*settings.NX + i_x], Tnew)
        # else:
        #     pool = mp.Pool(processes=settings.N_PROCESSES)
        #     T = np.array([
        #             pool.apply(get_t_new, args=(i_y, i_x, s, T))
        #         for i_x in range(settings.NX) for i_y in range(settings.NY)
        #     ])
        # quad 2
        for i_y in range(settings.NY):
            for i_x in range(settings.NX-1, -1, -1):
                Tnew = get_mean_t(i_y, i_x, s, T)
                T[i_y*settings.NX + i_x] = min(T[i_y*settings.NX + i_x], Tnew)
        # quad 3
        for i_y in range(settings.NY-1, -1, -1):
            for i_x in range(settings.NX-1, -1, -1):
                Tnew = get_mean_t(i_y, i_x, s, T)
                T[i_y*settings.NX + i_x] = min(T[i_y*settings.NX + i_x], Tnew)
        # quad 4
        for i_y in range(settings.NY-1, -1, -1):
            for i_x in range(settings.NX):
                Tnew = get_mean_t(i_y, i_x, s, T)
                T[i_y*settings.NX + i_x] = min(T[i_y*settings.NX + i_x], Tnew)
        if abs(np.linalg.norm(old_T) - np.linalg.norm(T)) < settings.SWEEP_THRESHOLD:
            break
    return T

def get_quick_curved_ray(x_s, y_s, x_r, y_r, s):
    T = get_travel_time(s, x_s, y_s)
    return ray.curved_ray(x_s, y_s, x_r, y_r, T)

def get_l(s, recalculate=True, is_save=False):
    _x_min = 2.0*settings.DX
    _x_max = settings.NX * settings.DX - 2.0*settings.DX
    x_step = [ _x_min + i*settings.STEP_SR for i in range(int((_x_max-_x_min)/settings.STEP_SR))]

    if recalculate:
        # print("Calculating curve length")
        l = []
        if not settings.is_parallel:
            for x_s in x_step:
                T = get_travel_time(s, x_s, settings.FIX_S_Y)
                for x_r in x_step:
                    if x_s != x_r:
                        l.append(ray.curved_ray(x_s, settings.FIX_S_Y, x_r, settings.FIX_R_Y, T))
        else:
            pool = mp.Pool(mp.cpu_count())
            for x_s in x_step: 
                l.extend(
                    pool.starmap(
                        get_quick_curved_ray,
                        [(x_s, settings.FIX_S_Y, x_r, settings.FIX_R_Y, s) for x_r in x_step]
                    )
                )
            pool.close()
        _L = np.array(l)
        # print("Finished ray tracing")
        if is_save:
            np.savez('cache.npz', L=_L)
            print("save to cache file")
    else:
        _L = np.load('cache.npz')['L']
        print("Ray path (L) has been loaded from the cache file!")
    return _L

'''
Optimization
'''
def get_r(t_obs, s_model, L, magnitude=True):
    if magnitude:
        r = np.linalg.norm(np.subtract(
            np.matmul(L, s_model),
            t_obs
        ))
    else:
        r = np.subtract(
            np.matmul(L, s_model),
            t_obs
        )
    return r

def grad(t_obs, s_model, L, norm=True):
    _grad = np.matmul(
        np.subtract(np.matmul(L, s_model), t_obs).T,
        L,
    )
    return np.divide(_grad, np.linalg.norm(_grad)) if norm else _grad

def smooth_map(s_real, kernel_size=(100, 100), mode='gaussian', pad_model='edge'):
    if mode == 'uniform':
        s_model = uniform_filter(s_real.reshape(settings.NX, settings.NY).T, size=kernel_size).reshape(settings.NX*settings.NY)
    elif mode == 'gaussian':
        if pad_model == 'zero':
            s_real_padded = np.pad(s_real.reshape(settings.NX, settings.NY).T, kernel_size, 'constant', constant_values=(0.0, 0.0))
        else:
            s_real_padded = np.pad(s_real.reshape(settings.NX, settings.NY).T, kernel_size, pad_model)
        s_model = gaussian_filter(s_real_padded, sigma=kernel_size)[kernel_size[0]:-kernel_size[0], kernel_size[1]:-kernel_size[1]].reshape(settings.NX*settings.NY)
    return s_model

'''
line search
'''
def get_proper_alpha(t_obs, s0, L0, pk, method='backtrack'):
    ALPHA0 = 1e-4
    C = 0.1
    if method == 'backtrack':
        ALPHA_DECAYRATE = 0.5
        alphak = ALPHA0
        while True:
            s1 = np.add(s0, np.multiply(alphak, pk))
            L1 = get_l(s1, recalculate=True)
            gradk = grad(t_obs, s1, L1)
            if get_r(t_obs, s1, L1) <= get_r(t_obs, s0, L0) + np.multiply(C * alphak, np.matmul(gradk.T, pk)):
                break
            alphak *= ALPHA_DECAYRATE
    elif method == 'cube_quad':
        alpha0 = ALPHA0
        s_alpha0 = np.add(s0, np.multiply(alpha0, pk))
        L_alpha0 = get_l(s_alpha0, recalculate=True)
        phi_0 = get_r(t_obs, s0, L0)
        grad_0 = grad(t_obs, s0, L0)
        phi_d0 = np.matmul(grad_0.T, pk)
        phi_alpha0 = get_r(t_obs, s_alpha0, L_alpha0)
        # quad
        alpha1_numerator = -np.multiply(phi_d0, alpha0**2)
        alpha1_denominator = 2.0*(phi_alpha0 - phi_0 - alpha0*phi_d0)
        alpha1 = alpha1_numerator/alpha1_denominator
        s_alpha1 = s0 + np.multiply(alpha1, pk)
        L_alpha1 = get_l(s_alpha1, recalculate=True)
        phi_alpha1 = get_r(t_obs, s_alpha1, L_alpha1)
        if alpha1 > 1e-6 and phi_alpha1 <= phi_0 + C * alpha1 * phi_d0:
            alphak = alpha1
        else:
            while True:
                s_alpha0 = np.add(s0, np.multiply(alpha0, pk))
                L_alpha0 = get_l(s_alpha0, recalculate=True)
                s_alpha1 = np.add(s0, np.multiply(alpha1, pk))
                L_alpha1 = get_l(s_alpha1, recalculate=True)
                phi_alpha0 = get_r(t_obs, s_alpha0, L_alpha0)
                phi_alpha1 = get_r(t_obs, s_alpha1, L_alpha1)

                alpha_matrix = np.array([[alpha0**2, -alpha1**2], [-alpha0**3, alpha1**3]])
                phi_matrix = np.array([phi_alpha1 - phi_0 - alpha1*phi_d0, phi_alpha0 - phi_0 - alpha0*phi_d0])
                ab = (1.0/(alpha0**2*alpha1**2*(alpha1-alpha0)))*np.matmul(alpha_matrix, phi_matrix)
                a, b = ab[0], ab[1]

                alpha2 = (-b + math.sqrt(b**2-3*a*phi_d0))/(3*a)
                s_alpha2 = s0 + np.multiply(alpha2, pk)
                L_alpha2 = get_l(s_alpha2, recalculate=True)
                phi_alpha2 = get_r(t_obs, s_alpha2, L_alpha2)
                if phi_alpha2 <= phi_0 + C * alpha2 * phi_d0:
                    alphak = alpha2
                    break
                alpha0 = alpha1
                alpha1 = alpha2
                print('alpha_0:{}, alpha_1:{}'.format(alpha0, alpha1))
    else:
        print('Method {} to find the step length does not support yet'.format(method))
        exit()
    return alphak
'''
    LM algorithm
'''

def get_jacobian(t_obs, s_model, pk):
    _j = []
    for j in range(t_obs.shape):
        get_r(t_obs, s_model, L, magnitude=False)