import numpy as np

# RAY TRACING
TRACING_MODE = 'cube' # circle, cube

# OPTIMIZATION
SMOOTH_GRAD = False
FILTER_GRAD = True

# GRID

NX = 220
NY = 160
DX = 2.0 # m

G_W = NX*DX
G_H = NY*DX

# DETECTOR
N_SOURCE = 90 # minumim 10
N_INSIDE_PYRAMID = 6 # min 3
N_OUTSIDE_PYRAMID_PER_SIDE = 3 # min 2
N_RECEIVER = N_INSIDE_PYRAMID + 2*N_OUTSIDE_PYRAMID_PER_SIDE

# PYRAMID
PYRAMID_W = 230
PYRAMID_H = 148
CHAMBER_W = 60
CHAMBER_H = 24

# void
VOID_W = 40
VOID_H = 40
VOID_HEIGHT = 60

# SIDE GAP
SIDE_GAP = (G_W - PYRAMID_W)/2.0

# MATERIALS
LAMBDA_ROCK = 1.0/12.18 # Fix E_muon = 2 GeV eq.4 [K.M. Tanaka 2016]
LAMBDA_AIR = 1.0/1e4 # https://arxiv.org/pdf/1208.1171.pdf