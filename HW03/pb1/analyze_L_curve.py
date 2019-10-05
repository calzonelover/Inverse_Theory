import pandas as pd
import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

f = pd.read_csv("l_curve.log")

res_ref, model_ref = f.iloc[0]['norm_res'], f.iloc[-1]['norm_model']
f["distance"] = distance = (f['norm_res'] - res_ref)**2 + (f['norm_model'] - model_ref)**2 
print(f)
f_corner = f[f.distance == f.distance.min()]


plt.plot(f['norm_res'], f['norm_model'], 'bo-', label="various alpha")
plt.plot(f_corner['norm_res'], f_corner['norm_model'], 'ro-', label="alpha = {:.4f}".format(float(f_corner['alpha'])))
plt.xlabel("$||t-ls||^2$")
plt.ylabel("$||s||^2$")
plt.xscale('log')
plt.yscale('log')
plt.title("L-curve for various alpha")
plt.legend()
plt.tight_layout()
# plt.savefig("l_curve.png")
plt.show()