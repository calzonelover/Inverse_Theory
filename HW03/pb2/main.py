import math
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

A, B = 1, 100
EPSILON = 1e-2
ALPHA = 1e-5

def func(x1, x2, a=A, b=B):
    return (a-x1)**2 + b*(x2-x1**2)**2

if __name__ == '__main__':
    x, y = np.meshgrid(np.linspace(-5, 5, 1000),
                   np.linspace(-5, 5, 1000))
    plt.contourf(x, y, func(x, y), locator=matplotlib.ticker.LogLocator(), cmap='gnuplot')
    plt.colorbar()
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title("Rosenbrock Potential")
    plt.savefig("contour.png")