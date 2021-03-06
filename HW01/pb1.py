import numpy as np
import matplotlib.pyplot as plt

N = 20
a = 1.0
b = -0.5
SIGMA = 0.2

MAX_X = 5.0
MIN_X = -5.0

def linear_model(x):
	n = np.random.normal(loc=0, scale=SIGMA)
	return a*x + b + n

def get_sampling(n, max_x=MAX_X, min_x=MIN_X):
	step_x = (max_x - min_x) / n
	x_samples = np.array([i*step_x for i in range(n)])
	y_samples = np.array(list(map(lambda x: linear_model(x), x_samples)))
	return {'x': x_samples.reshape(-1,1), 'y': y_samples.reshape(-1,1)}

def get_model_by_slq(A, d):
	inverse_A = np.linalg.inv(np.matmul(A.T, A))
	return np.matmul(inverse_A.T, np.matmul(A.T, d))

def get_model_by_wlq(A, d, sigma):
	W = np.divide(np.identity(N, dtype=float), sigma)
	A_W, d_W = np.matmul(W, A), np.matmul(W, d)
	inverse_A_W = np.linalg.inv(np.matmul(A_W.T, A_W))
	return np.matmul(inverse_A_W.T, np.matmul(A_W.T, d_W))

if __name__ == "__main__":
	samples = get_sampling(N)
	array_ones = np.ones(N, dtype=float).reshape(-1,1)
	A = np.concatenate((samples['x'], array_ones), axis=1)
	d = samples['y']

	print("Standard least squares method")
	m_slq = get_model_by_slq(A, d)
	print("Parameter a = {}, b = {}".format(m_slq[0, 0], m_slq[1, 0]))
	print("Weight least squares method")
	m_mlq = get_model_by_wlq(A, d, SIGMA)
	print("Parameter a = {}, b = {}".format(m_mlq[0, 0], m_mlq[1, 0]))

	# visual
	plt.plot(samples['x'], samples['y'], 'ro', label="y = ax + b + $\mathcal{N}$(0, 0.2)")
	plt.ylabel("y")
	plt.xlabel("x")
	plt.legend()
	plt.savefig("dist1.png")

