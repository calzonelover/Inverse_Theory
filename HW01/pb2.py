import numpy as np
import matplotlib.pyplot as plt

N = 20
N_R = 10
a = 1.0
b = -0.5
SIGMA = 0.2

MAX_X = 5.0
MIN_X = -5.0

def linear_model(x):
	n = np.random.normal(loc=0, scale=SIGMA)
	return a*x + b + n

def get_sampling(n, n_r, max_x=MAX_X, min_x=MIN_X):
	step_x = (max_x - min_x) / n
	x_samples = np.array([i*step_x for i in range(n)])
	x_sample_vecs = np.array(list(map(lambda x_sample: np.multiply(np.ones(n_r), x_sample), x_samples)))
	y_samples_vecs = np.array(list(map(
		lambda x_sample_vec:
			list(map(lambda x: linear_model(x), 
			x_sample_vec)),
		x_sample_vecs)))
	y_mean_samples = np.array(list(map(lambda x: np.mean(x), y_samples_vecs)))
	y_std_samples = np.array(list(map(lambda x: np.std(x), y_samples_vecs)))
	return {
		'x': x_samples.reshape(-1, 1),
		'x_vec': x_sample_vecs.reshape(-1,n_r),
		'x_all': x_sample_vecs.reshape(-1, 1),
		'y_vec': y_samples_vecs.reshape(-1,n_r),
		'y_all': y_samples_vecs.reshape(-1, 1),
		'y_mean': y_mean_samples.reshape(-1,1),
		'y_std': y_std_samples.reshape(-1,1),
	}

def get_model_by_slq(A, d):
	inverse_A = np.linalg.inv(np.matmul(A.T, A))
	return np.matmul(inverse_A.T, np.matmul(A.T, d))

def get_model_by_wlq(A, d, sigma):
	inv_sigma = np.divide(np.ones(sigma.shape), sigma)
	W = np.multiply(np.identity(N, dtype=float), inv_sigma)
	A_W, d_W = np.matmul(W, A), np.matmul(W, d)
	inverse_A_W = np.linalg.inv(np.matmul(A_W.T, A_W))
	return np.matmul(inverse_A_W.T, np.matmul(A_W.T, d_W))

if __name__ == "__main__":
	samples = get_sampling(N, N_R)

	print("Standard least squares method")
	array_ones = np.ones(N*N_R, dtype=float).reshape(-1,1)
	A = np.concatenate((samples['x_all'], array_ones), axis=1)
	d = samples['y_all']
	m_slq = get_model_by_slq(A, d)
	print("Parameter a = {}, b = {}".format(m_slq[0, 0], m_slq[1, 0]))

	print("Weight least squares method")
	array_ones = np.ones(N, dtype=float).reshape(-1,1)
	A = np.concatenate((samples['x'], array_ones), axis=1)
	d = samples['y_mean']
	sigmas = samples['y_std']
	m_mlq = get_model_by_wlq(A, d, sigmas)
	print("Parameter a = {}, b = {}".format(m_mlq[0, 0], m_mlq[1, 0]))

	# visual
	plt.plot(samples['x_all'], samples['y_all'], 'ro', label="y = ax + b + $\mathcal{N}$(0, 0.2); N_R = 5")
	plt.plot(samples['x'], samples['y_mean'], 'bo', label="y = <ax + b + $\mathcal{N}$(0, 0.2)>")
	plt.ylabel("y")
	plt.xlabel("x")
	plt.legend()
	plt.savefig("dist2.png")