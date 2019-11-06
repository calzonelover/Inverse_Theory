import math
import numpy as np
import platform
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

def get_travel_time(s, source_x, source_y, n_sweep=settings.N_SWEEP):
    T = np.multiply(1e6, np.ones(settings.NY*settings.NX))
    i_source_x = math.floor(source_x/settings.DX)
    i_source_y = math.floor(source_y/settings.DX)
    T[i_source_y*settings.NX + i_source_x] = 0.0
    for i_sweep in range(n_sweep):
        old_T = T
        # quad 1
        for i_y in range(settings.NY):
            for i_x in range(settings.NX):
                Tnew = get_mean_t(i_y, i_x, s, T)
                T[i_y*settings.NX + i_x] = min(T[i_y*settings.NX + i_x], Tnew)
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

def get_l(s, recalculate=False):
    _x_min = 2.0*settings.DX
    _x_max = settings.NX * settings.DX - 2.0*settings.DX
    x_step = [ _x_min + i*settings.STEP_SR for i in range(int((_x_max-_x_min)/settings.STEP_SR))]

    if recalculate:
        print("Calculating curve ray length")
        l = []
        for x_s in x_step:
            for x_r in x_step:
                if x_s != x_r:
                    print(x_s, x_r)
                    T = get_travel_time(s, x_s, 0.0)
                    l.append(ray.curved_ray(x_s, 0.0, x_r, 0.0, T))
        _L = np.array(l)
        print("Finished ray tracing")
        np.savez('cache.npz', L=_L)
        print("save to cache file")
    else:
        _L = np.load('cache.npz')['L']
        print("Ray path (L) has been loaded from the cache file!")
    return _L

'''
Optimization
'''
def get_r(t_obs, s_model, L):
    return np.square(np.linalg.norm(np.subtract(
        np.matmul(L, s_model),
        t_obs
    )))

def grad(t_obs, s_model, L, norm=True):
    _grad = np.matmul(
        np.subtract(np.matmul(L, s_model), t_obs).T,
        L,
    )
    return np.divide(_grad, np.linalg.norm(_grad)) if norm else _grad

'''
line search
'''
def get_proper_alpha(t_obs, s0, L, pk, method='backtrack'):
    if method == 'backtrack':
        ALPHA0 = 1.0
        C = 0.1
        ALPHA_DECAYRATE = 0.8

        alphak = ALPHA0
        while True:
            s1 = s0 + np.multiply(alphak, pk)
            gradk = grad(t_obs, s1, L)
            if get_r(t_obs, s1, L) <= get_r(t_obs, s0, L) + np.multiply(C * alphak, np.matmul(gradk.T, pk)):
                break
            alphak *= ALPHA_DECAYRATE
        # alphak = ALPHA0
        # while True:
        #     x1k1 = x1_0 + alphak * pk[0]
        #     x2k1 = x2_0 + alphak * pk[1]
        #     gradk = grad_func(x1_0, x2_0)
        #     if func(x1k1, x2k1) <= func(x1_0, x2_0) + C * alphak * np.matmul(gradk.T, pk):
        #         break
        #     alphak *= ALPHA_DECAYRATE
    # elif method == 'cube_quad':
    #     alpha0 = ALPHA0
    #     x1_alpha0 = x1_0 + alpha0 * pk[0]
    #     x2_alpha0 = x2_0 + alpha0 * pk[1]
    #     phi_0 = func(x1_0, x2_0)
    #     grad_0 = grad_func(x1_0, x2_0)
    #     phi_d0 = np.matmul(grad_0.T, pk)
    #     phi_alpha0 = func(x1_alpha0, x2_alpha0)
    #     # quad
    #     alpha1_numerator = -np.multiply(phi_d0, alpha0**2)
    #     alpha1_denominator = 2.0*(phi_alpha0 - phi_0 - alpha0*phi_d0)
    #     alpha1 = alpha1_numerator/alpha1_denominator
    #     x1_alpha1 = x1_0 + alpha1 * pk[0]
    #     x2_alpha1 = x2_0 + alpha1 * pk[1]
    #     phi_alpha1 = func(x1_alpha1, x2_alpha1)
    #     if alpha1 > 1e-5 and phi_alpha1 <= phi_0 + C * alpha1 * phi_d0:
    #         alphak = alpha1
    #     else:
    #         while True:
    #             x1_alpha0 = x1_0 + alpha0 * pk[0]
    #             x2_alpha0 = x2_0 + alpha0 * pk[1]
    #             x1_alpha1 = x1_0 + alpha1 * pk[0]
    #             x2_alpha1 = x2_0 + alpha1 * pk[1]
    #             phi_alpha0 = func(x1_alpha0, x2_alpha0)
    #             phi_alpha1 = func(x1_alpha1, x2_alpha1)

    #             alpha_matrix = np.array([[alpha0**2, -alpha1**2], [-alpha0**3, alpha1**3]])
    #             phi_matrix = np.array([phi_alpha1 - phi_0 - alpha1*phi_d0, phi_alpha0 - phi_0 - alpha0*phi_d0])
    #             ab = (1.0/(alpha0**2*alpha1**2*(alpha1-alpha0)))*np.matmul(alpha_matrix, phi_matrix)
    #             a, b = ab[0], ab[1]

    #             alpha2 = (-b + math.sqrt(b**2-3*a*phi_d0))/(3*a)
    #             x1_alpha2 = x1_0 + alpha2 * pk[0]
    #             x2_alpha2 = x2_0 + alpha2 * pk[1]
    #             phi_alpha2 = func(x1_alpha2, x2_alpha2)
    #             if phi_alpha2 <= phi_0 + C * alpha2 * phi_d0:
    #                 alphak = alpha2
    #                 break

    #             alpha0 = alpha1
    #             alpha1 = alpha2
    else:
        print('Please select the satisfiable method to find the step length')
        exit()
    return alphak