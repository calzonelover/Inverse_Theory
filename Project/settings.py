import numpy as np

# parallel
is_parallel = True

NX = 220
NY = 160
DX = 2.0 # m

G_W = NX*DX
G_H = NY*DX

# DETECTOR
N_SOURCE = 10
N_INSIDE_PYRAMID = 3
N_OUTSIDE_PYRAMID_PER_SIDE = 2

# PYRAMID
PYRAMID_W = 230
PYRAMID_H = 148
CHAMBER_W = 60
CHAMBER_H = 24

# void
VOID_W = 30
VOID_H = 30
VOID_HEIGHT = 60

# SIDE GAP
SIDE_GAP = (G_W - PYRAMID_W)/2.0

# MATERIALS
LAMBDA_ROCK = 10.0 # https://www.intechopen.com/books/cosmic-rays/cosmic-ray-muons-as-penetrating-probes-to-explore-the-world-around-us
LAMBDA_AIR = 1e-5