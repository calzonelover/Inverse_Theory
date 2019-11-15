import numpy as np
import os
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import environment, utility, settings

def main():
    pair_srs = utility.get_source_receiver()
    # Test ray tracing
    for k, sr_pos in enumerate(pair_srs):
        print(sr_pos)
        s_pos = sr_pos['s_pos']
        r_pos = sr_pos['r_pos']

        ray_map = utility.ray_length(s_pos['x'],s_pos['y'],r_pos['x'],r_pos['y'], mode='circle', is_fast_tracing=False)
        plt.imshow(
            ray_map.reshape(settings.NY, settings.NX), cmap='jet',
            extent=[0, settings.DX*settings.NX, 0, settings.DX*settings.NY],
            origin='bottom left',
        )
        a = plt.colorbar()
        a.set_label('Ray length (m)')
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title("Ray Tracing (source: ({0:.2f},{1:.2f}), receiver({2:.2f},{3:.2f}))".format(s_pos['x'],s_pos['y'],r_pos['x'],r_pos['y']))
        plt.savefig(os.path.join('unit_test', 'log_muon_tracing', "ray_{}.png".format(k)))
        plt.clf()