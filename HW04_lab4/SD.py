import math
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

A, B = 1, 100
C = 0.2
EPSILON = 1e-6
ALPHA0 = 1.0
ALPHA_DECAYRATE = 0.8

def func(x1, x2, a=A, b=B):
    return (a-x1)**2 + b*(x2-x1**2)**2
def grad_func(x1, x2, a=A, b=B):
    return np.array([
        2*(x1-a) + 4*b*x1*(x1**2-x2),
        2*b*(x2-x1**2)
    ])
def hess_func(x1, x2, a=A, b=B):
    return np.array([
        [2+4*b*(3*x1**2-x2), -4*b*x1],
        [-4*b*x1, 2*b],
    ])

def get_proper_alpha(x1k, x2k, pk):
    alphak = ALPHA0
    while True:
        x1k1 = x1k + alphak * pk[0]
        x2k1 = x2k + alphak * pk[1]
        gradk = grad_func(x1k, x2k)
        if func(x1k1, x2k1) <= func(x1k, x2k) + C * alphak * np.matmul(gradk.T, pk):
            break
        alphak *= ALPHA_DECAYRATE
    return alphak

if __name__ == '__main__':
    x1, x2 = 3.0, -2.0
    r = func(x1, x2)

    k = 0
    ps = []
    x1s = []
    x2s = []
    x1s.append(x1)
    x2s.append(x2)
    while r > EPSILON:
        pk = -1.0*grad_func(x1s[k],x2s[k])
        alphak = get_proper_alpha(x1s[k], x2s[k], pk)
        x1k1 = x1s[k] + alphak * pk[0]
        x2k1 = x2s[k] + alphak * pk[1]

        ps.append(pk)
        x1s.append(x1k1)
        x2s.append(x2k1)
        k +=1
        r = func(x1k1, x2k1)
    ## Visualize
    # contour
    x, y = np.meshgrid(np.linspace(-5, 5, 1000),
                   np.linspace(-5, 5, 1000))
    plt.contourf(x, y, func(x, y), locator=matplotlib.ticker.LogLocator(), cmap='gnuplot')
    plt.colorbar()
    # point
    plt.plot(x1s, x2s, 'g.-')
    plt.plot(x1s[0], x2s[0], 'bo',label='begin ({:.2f}, {:.2f})'.format(x1, x2))
    plt.plot(x1s[-1], x2s[-1], 'ro',label='end ({:.2f}, {:.2f})'.format(x1s[-1], x2s[-1]))
    plt.legend()
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title("Rosenbrock Potential")
    plt.show()