import numpy as np
import os
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import environment, utility, settings

TRACING_MODE = 'cube' # circle, cube

def main():
    pair_srs = utility.get_source_receiver()
    # Test ray tracing
    for k, sr_pos in enumerate(pair_srs):
        print(sr_pos)
        s_pos = sr_pos['s_pos']
        r_pos = sr_pos['r_pos']

        ray_map = utility.ray_length(s_pos['x'],s_pos['y'],r_pos['x'],r_pos['y'], mode=TRACING_MODE, is_fast_tracing=False)
        plt.imshow(
            ray_map.reshape(settings.NY, settings.NX), cmap='jet',
            extent=[0, settings.DX*settings.NX, 0, settings.DX*settings.NY],
            origin='bottom left',
        )
        a = plt.colorbar()
        a.set_label('Ray length (m)')
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title("{0} Ray Tracing (src: ({1:.1f},{2:.1f}), rec: ({3:.1f},{4:.1f}))".format(TRACING_MODE,s_pos['x'],s_pos['y'],r_pos['x'],r_pos['y']))
        plt.savefig(os.path.join('unit_test', 'log_muon_tracing', "{}_ray_{}.png".format(TRACING_MODE, k)))
        plt.clf()