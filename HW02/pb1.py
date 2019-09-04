import numpy as np
import matplotlib.pyplot as plt
import scipy.io

BLUR_IMG_MAT_FILE = 'd.mat'
EPSILON = 1e-5

def load_blur_img(mat_file):
    return np.array(scipy.io.loadmat(mat_file)['d'], dtype=float).reshape(1,-1)

def vec_to_img(vec):
    return vec.reshape(100, 100).T

if __name__ == '__main__':
    blur_img_vec = load_blur_img(BLUR_IMG_MAT_FILE)
    blur_img = vec_to_img(blur_img_vec)

    # visual
    plt.imshow(blur_img, cmap=plt.cm.gray_r)
    plt.show()