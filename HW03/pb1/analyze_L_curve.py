import pandas as pd
import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

f = pd.read_csv("l_curve.log")

print(f['norm_res'], f['norm_model'])

plt.plot(f['norm_res'], f['norm_model'], 'o-')
plt.xscale('log')
plt.yscale('log')
plt.show()