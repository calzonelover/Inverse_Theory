import math
import numpy as np
A, B = 1, 100
EPSILON = 1e-7
ALPHA0 = 1.0
ALPHA_DECAYRATE = 0.7

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

if __name__ == '__main__':
    x1, x2 = 2.0, 2.0
    r = func(x1, x2)

    k = 0
    ps = alphas = x1s = x2s = []
    x1s.append(x1)
    x2s.append(x2)
    while r > EPSILON:
        ps[k] = grad_func(x1s[k],x2s[k])
        # find alpha
        alphak = ALPHA0
        while True:
            x1s[k+1] = x1s[k] + alphak * ps[k][0]
            x2s[k+1] = x2s[k] + alphak * ps[k][1]
            print(x1s[k+1], x2s[k+1])
            gradk = grad_func(x1s[k], x2s[k])
            if func(x1s[k+1], x2s[k+1]) < func(x1s[k], x2s[k]) + alphak * np.matmul(gradk.T, ps[k]):
                break
            alphak /= ALPHA_DECAYRATE
        # end find alpha
        k +=1
        r = func(x1, x2)
    print(x1s[-1], x2s[-1])