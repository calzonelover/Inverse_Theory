import numpy as np
import platform
if platform.system() == "Darwin":
    import matplotlib
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# setting
FILENAME = "src/vel_nx50_nz50_dx20.dat"
NX = NY = 50
DX = 20

def readmap(filename):
    f = open(FILENAME, "r")
    a = np.fromfile(f, dtype=np.float32)
    a = a.reshape(NX, NY)
    return a.T

if __name__ == "__main__":
    raw_map = readmap(filename=FILENAME)

    # visualize
    plt.imshow(raw_map, cmap='jet', extent=[0, DX*NX, 0, DX*NY])
    a = plt.colorbar()
    a.set_label('$v$')
    plt.title("Real Velocity")
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")
    plt.savefig("real_v.png")
    # plt.show()