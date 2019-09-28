import numpy as np

FILENAME = "src/vel_nx50_nz50_dx20.dat"
NX = NY = 50
DX = 20

N_SOURCE = 50
N_RECEIVER = 50

# ALPHAS = [1e-3, ]
ALPHAS = np.logspace(-4, 4, num=20, base=10)