import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import utility, settings

def main():
    v = utility.readraw(filename=settings.FILENAME)
    s = 1.0/v
    L = utility.get_l(s, recalculate=True)
    print(L.shape)
    