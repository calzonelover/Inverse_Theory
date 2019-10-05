import math
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

A, B = 1, 100
C = 0.1
EPSILON = 1e-2
ALPHA0 = 1.0

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

def get_proper_alpha(x1_0, x2_0, pk):
    alpha0 = ALPHA0
    x1_alpha0 = x1_0 + alpha0 * pk[0]
    x2_alpha0 = x2_0 + alpha0 * pk[1]
    phi_0 = func(x1_0, x2_0)
    grad_0 = grad_func(x1_0, x2_0)
    phi_d0 = np.matmul(grad_0.T, pk)
    phi_alpha0 = func(x1_alpha0, x2_alpha0)
    # quad
    alpha1_numerator = -np.multiply(phi_d0, alpha0**2)
    alpha1_denominator = 2.0*(phi_alpha0 - phi_0 - alpha0*phi_d0)
    alpha1 = alpha1_numerator/alpha1_denominator
    x1_alpha1 = x1_0 + alpha1 * pk[0]
    x2_alpha1 = x2_0 + alpha1 * pk[1]
    phi_alpha1 = func(x1_alpha1, x2_alpha1)
    if phi_alpha1 <= phi_0 + C * alpha1 * phi_d0:
        alphak = alpha1
    else:
        while True:
            x1_alpha0 = x1_0 + alpha0 * pk[0]
            x2_alpha0 = x2_0 + alpha0 * pk[1]
            x1_alpha1 = x1_0 + alpha1 * pk[0]
            x2_alpha1 = x2_0 + alpha1 * pk[1]
            phi_alpha0 = func(x1_alpha0, x2_alpha0)
            phi_alpha1 = func(x1_alpha1, x2_alpha1)

            alpha_matrix = np.array([[alpha0**2, -alpha1**2], [-alpha0**3, alpha1**3]])
            phi_matrix = np.array([phi_alpha1 - phi_0 - alpha1*phi_d0, phi_alpha0 - phi_0 - alpha0*phi_d0])
            ab = (1.0/(alpha0**2*alpha1**2*(alpha1-alpha0)))*np.matmul(alpha_matrix, phi_matrix)
            a, b = ab[0], ab[1]

            alpha2 = (-b + math.sqrt(b**2-3*a*phi_d0))/(3*a)
            x1_alpha2 = x1_0 + alpha2 * pk[0]
            x2_alpha2 = x2_0 + alpha2 * pk[1]
            phi_alpha2 = func(x1_alpha2, x2_alpha2)
            if phi_alpha2 <= phi_0 + C * alpha2 * phi_d0:
                alphak = alpha2
                break

            alpha0 = alpha1
            alpha1 = alpha2
    return alphak

if __name__ == '__main__':
    x1, x2 = 3.0, -2.0
    r = func(x1, x2)

    k = 0
    ps = []
    x1s = []
    x2s = []
    rs = []
    x1s.append(x1)
    x2s.append(x2)

    gradk = grad_func(x1,x2)
    pk = -1.0*gradk
    while np.linalg.norm(gradk) > EPSILON:
        alphak = get_proper_alpha(x1s[k], x2s[k], pk)
        x1k1 = x1s[k] + alphak * pk[0]
        x2k1 = x2s[k] + alphak * pk[1]
        
        gradk = grad_func(x1s[k],x2s[k])
        gradk1 = grad_func(x1k1,x2k1)
        betak1 = np.matmul(gradk1.T, gradk1)/np.matmul(gradk.T, gradk)
        pk = np.add(-1.0*gradk1, np.multiply(pk, betak1))

        ps.append(pk)
        x1s.append(x1k1)
        x2s.append(x2k1)
        k += 1
        r = func(x1k1, x2k1)
        rs.append(r)
        print(len(rs), r, np.linalg.norm(gradk))
    ## Visualize
    # contour
    x, y = np.meshgrid(
        np.linspace(-5, 5, 1000),
        np.linspace(-5, 5, 1000)
    )
    plt.contourf(x, y, func(x, y), locator=matplotlib.ticker.LogLocator(), cmap='gnuplot')
    plt.colorbar()
    # point
    plt.plot(x1s, x2s, 'g.-')
    plt.plot(x1s[0], x2s[0], 'bo',label='begin ({:.2f}, {:.2f})'.format(x1, x2))
    plt.plot(x1s[-1], x2s[-1], 'ro',label='end ({:.2f}, {:.2f})'.format(x1s[-1], x2s[-1]))
    plt.legend()
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title("Rosenbrock Potential (Fixed Alpha)")
    # plt.show()
    plt.savefig("cg.png")
    # decay
    plt.clf()
    plt.plot(rs)
    plt.yscale("log")
    plt.title("Decay of objective function")
    plt.savefig("cg_r.png")