import numpy as np

# parallel
is_parallel = True

FILENAME = "vel_nx461_nz151_dx20.dat"
NX = 461
NY = 151
DX = 20
STEP_SR = 1600.0
FIX_S_Y = 40
FIX_R_Y = 60

is_reduce_size = True
FACTOR = 3
if is_reduce_size:
    NX = int(NX/FACTOR)
    NY = int(NY/FACTOR)
    DX = 20.0*FACTOR

# COLORBAR
COLOR_VMAX = 5800
COLOR_VMIN = 1500

# SWEEPING
N_SWEEP = 100
SWEEP_THRESHOLD = 1e-15