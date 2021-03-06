import math
import numpy as np
import platform
import settings

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
    try:
        T_x_min = min(old_T[i_y*settings.NX+(i_x-1)], old_T[i_y*settings.NX+(i_x+1)])
    except IndexError:
        try:
            T_x_min = old_T[i_y*settings.NX+(i_x-1)]
        except IndexError:
            T_x_min = old_T[i_y*settings.NX+(i_x+1)]
    try:
        T_y_min = min(old_T[(i_y-1)*settings.NX+i_x], old_T[(i_y+1)*settings.NX+i_x])
    except IndexError:
        try:
            T_y_min = old_T[(i_y-1)*settings.NX+i_x]
        except IndexError:
            T_y_min = old_T[(i_y+1)*settings.NX+i_x]
    if abs(T_x_min - T_y_min) >= s[i_y*settings.NX + i_x]*settings.DX:
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

def get_travel_time(s, source_x, source_y, n_sweep):
    T = np.multiply(1e6, np.ones(s.shape))
    i_source_x = math.floor(source_x/settings.DX) + 1
    i_source_y = math.floor(source_y/settings.DX) + 1
    T[i_source_y*settings.NX + i_source_x] = 0.0
    for i_sweep in range(n_sweep):
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
    return T