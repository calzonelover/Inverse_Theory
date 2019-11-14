import multiprocessing as mp
import numpy as np
import math

import settings

def is_inside_pyramid(x, y):
    slope = settings.PYRAMID_H/(settings.PYRAMID_W/2.0)
    y_left = lambda x: slope*x - slope*settings.SIDE_GAP
    y_right = lambda x: -slope*x + slope*(settings.PYRAMID_W + settings.SIDE_GAP)
    return True if  y < y_left(x) and y < y_right(x) else False

def is_main_chamber(x, y):
    x_min = settings.SIDE_GAP + settings.PYRAMID_W/2.0 - settings.CHAMBER_W/2.0
    x_max = settings.SIDE_GAP + settings.PYRAMID_W/2.0 + settings.CHAMBER_W/2.0
    y_min = 0.0
    y_max = settings.CHAMBER_H
    return True if (
        x > x_min and x < x_max and y > y_min and y < y_max
    ) else False

def is_hidden_void(x, y):
    x_min = settings.SIDE_GAP + settings.PYRAMID_W/2.0 - settings.VOID_W/2.0
    x_max = settings.SIDE_GAP + settings.PYRAMID_W/2.0 + settings.VOID_W/2.0
    y_min = settings.VOID_HEIGHT
    y_max = settings.VOID_HEIGHT + settings.VOID_H
    return True if (
        x > x_min and x < x_max and y > y_min and y < y_max
    ) else False

def get_element(x, y, model='initial'):
    if is_inside_pyramid(x, y):
        if is_main_chamber(x, y):
            _lambda = settings.LAMBDA_AIR
        elif model=='real' and is_hidden_void(x, y):
            _lambda = settings.LAMBDA_AIR
        else:
            _lambda = settings.LAMBDA_ROCK
    else:
        _lambda = settings.LAMBDA_AIR
    return _lambda

def get_system(model='initial'):
    _model = np.zeros(shape=(settings.NX*settings.NY))
    for i_y in range(settings.NY):
        for i_x in range(settings.NX):
            x = settings.DX * i_x + (settings.DX/2.0)
            y = settings.DX * i_y + (settings.DX/2.0)
            _model[i_x + settings.NX * i_y] = get_element(x, y, model)
    return _model