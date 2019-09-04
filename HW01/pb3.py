import numpy as np
import matplotlib.pyplot as plt
import random

N = 20
N_OUTLIER = 2
DISTANCE_BIAS = 3
EPSILON = 1e-6 # tolerance

a = 1.0
b = -0.5
SIGMA = 0.2

MAX_X = 5.0
MIN_X = -5.0

def linear_model(x):
	n = np.random.normal(loc=0, scale=SIGMA)
	return a*x + b + n

def get_sampling(n, n_outlier, distance_bias=DISTANCE_BIAS, max_x=MAX_X, min_x=MIN_X):
    step_x = (max_x - min_x) / n
    x_samples = np.array([i*step_x for i in range(n)])
    y_samples = np.array(list(map(lambda x: linear_model(x), x_samples)))

    index_outliers = np.random.randint(low=1, high=n, size=n_outlier)
    for index_outlier in index_outliers:
        y_samples[index_outlier] += random.choice([distance_bias, -distance_bias])
    return {'x': x_samples.reshape(-1,1), 'y': y_samples.reshape(-1,1)}

# SLQ
def get_model_by_slq(A, d):
	inverse_A = np.linalg.inv(np.matmul(A.T, A))
	return np.matmul(inverse_A.T, np.matmul(A.T, d))

# L1
def get_R(A, d, m_ki):
	r = np.subtract(d, np.matmul(A, m_ki))
	return np.linalg.inv(np.multiply(np.identity(N, dtype=float), r))

def get_m(A, d, R_ki):
	inv_term = np.linalg.inv(np.matmul(A.T, np.matmul(R_ki, A)))
	return np.matmul(inv_term, np.matmul(A.T, np.matmul(R_ki, d)))

def get_e(m_ki, m_kf):
	return np.linalg.norm(np.subtract(m_kf, m_ki))/(1+np.linalg.norm(m_kf))

def get_model_by_l1r(A, d, m_0, epsilon):
	m_ki = m_0
	while True:
		R_ki = get_R(A, d, m_ki)
		m_kf = get_m(A, d, R_ki)
		e_now = get_e(m_ki, m_kf)
		if e_now < epsilon:
			break
		m_ki = m_kf
		# print(e_now, m_kf)
	return m_kf

if __name__ == "__main__":
	samples = get_sampling(N, N_OUTLIER)
	array_ones = np.ones(N, dtype=float).reshape(-1,1)
	A = np.concatenate((samples['x'], array_ones), axis=1)
	d = samples['y']

	print("Standard least squares method")
	m_slq = get_model_by_slq(A, d)
	print("Parameter a = {}, b = {}".format(m_slq[0, 0], m_slq[1, 0]))

	print("L1 Regression method")
	m_l1r = get_model_by_l1r(A, d, m_slq, EPSILON)
	print("Parameter a = {}, b = {}".format(m_l1r[0, 0], m_l1r[1, 0]))

    # visual
	plt.plot(samples['x'], samples['y'], 'ro', label="y = ax + b + $\mathcal{N}$(0, 0.2)")
	plt.ylabel("y")
	plt.xlabel("x")
	plt.legend()
	plt.savefig("dist3.png")