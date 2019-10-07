import numpy as np

FILENAME = "src/vel_nx50_nz50_dx20.dat"
NX = NY = 50
DX = 20

N_SOURCE = 50
N_RECEIVER = 50

# ALPHAS = [1e-3, ]
ALPHAS = np.logspace(-3, 3, num=50, base=10)

# COLORBAR
COLOR_VMAX = 2000
COLOR_VMIN = 1500