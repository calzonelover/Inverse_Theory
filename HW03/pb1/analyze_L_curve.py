import pandas as pd
import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

f = pd.read_csv("l_curve.log")

row_corner = f[f.norm_res == f.norm_res.min()]

plt.plot(f['norm_res'], f['norm_model'], 'bo-', label="various alpha")
plt.plot(row_corner['norm_res'], row_corner['norm_model'], 'ro-', label="alpha = {:.4f}".format(float(row_corner['alpha'])))
plt.xlabel("$||t-ls||^2$")
plt.ylabel("$||m||^2$")
plt.xscale('log')
plt.yscale('log')
plt.title("L-curve for various alpha")
plt.legend()
plt.savefig("l_curve.png")
# plt.show()