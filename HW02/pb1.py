import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.io

WIDTH = 100
HEIGHT = 100
BLUR_IMG_MAT_FILE = 'd.mat'
EPSILON = 1e-3
D = 2


def load_blur_img(mat_file):
    return np.array(scipy.io.loadmat(mat_file)['d'], dtype=float).reshape(1,-1)

def vec_to_img(vec):
    return vec.reshape(WIDTH, HEIGHT)

def get_position(i):
    return 0.5 + float(i)

def get_identity_tensor(d, n=WIDTH):
    out = np.zeros( (n,) * d )
    out[ tuple([np.arange(n)] * d) ] = 1
    return out

if D == 2:
    # K as tensor rank 2 (Matrix)
    def kernel(i, ip):
        i, j = (i % WIDTH), math.floor(i/WIDTH)
        ip, jp = (ip % WIDTH), math.floor(ip/WIDTH)
        x, xp = get_position(i), get_position(j)
        y, yp = get_position(ip), get_position(jp)
        square_diff_x = (x-xp)**2 + (y-yp)**2
        return math.e**(-square_diff_x/2.0)/math.sqrt(math.pi)

    def get_k_tensor():
        init_k = np.zeros(shape=(WIDTH*HEIGHT, WIDTH*HEIGHT), dtype=float)
        for i in range(WIDTH*HEIGHT):
            for ip in range(WIDTH*HEIGHT):
                init_k[i, ip] = kernel(i, ip)
        return init_k
else:
    # K as tensor rank 4
    def kernel(i, j, ip, jp):
        x, xp = get_position(i), get_position(j)
        y, yp = get_position(ip), get_position(jp)
        square_diff_x = (x-xp)**2 + (y-yp)**2
        return math.e**(-square_diff_x/2.0)/math.sqrt(math.pi)

    def get_k_tensor():
        init_k = np.zeros(shape=(WIDTH, HEIGHT, WIDTH, HEIGHT), dtype=float)
        for i in range(WIDTH):
            for j in range(HEIGHT):
                for ip in range(WIDTH):
                    for jp in range(HEIGHT):
                        init_k[i, j, ip, jp] = kernel(i, j, ip, jp)
        return init_k

if __name__ == '__main__':
    blur_img_vec = load_blur_img(BLUR_IMG_MAT_FILE).reshape(-1, 1)
    blur_img = vec_to_img(blur_img_vec)

    A = get_k_tensor()
    print(A)

    if D == 2:
        inv_term = np.linalg.inv(np.add(
                        np.matmul(A.T, A),
                        np.multiply(np.identity(WIDTH*HEIGHT, dtype=float), EPSILON)
                    ))
        clear_img_vec = np.matmul(inv_term, np.matmul(A.T, blur_img_vec))
        clear_img = clear_img_vec.reshape(WIDTH, WIDTH)
    else:
        tensor_identity_ep = EPSILON*get_identity_tensor(4, n=WIDTH)
        inv_term = np.linalg.inv(np.add(
                        np.matmul(A.T, A),
                        tensor_identity_ep
                    ))
        clear_img = np.matmul(inv_term, np.matmul(A.T, blur_img))

    # visual
    plt.imshow(blur_img, cmap=plt.cm.gray_r)
    plt.show()

    plt.clf()

    plt.imshow(clear_img, cmap=plt.cm.gray_r)
    plt.show()