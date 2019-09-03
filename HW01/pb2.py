import numpy as np

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
	for i in range(n_r):
		y_samples = np.array([list(map(lambda x: linear_model(x), x_samples)) for i in range(n_r)])
	# y_mean_samples = np.array([list(map(lambda x: linear_model(x), x_samples)) for i in range(N_R)])
	print(np.array(y_samples).reshape(-1, n_r).shape)
	print(np.array(y_samples).reshape(-1, n_r))
	return {'x': x_samples.reshape(-1,1), 'y': y_samples.reshape(-1,1)}

# def 

if __name__ == "__main__":
	samples = get_sampling(N, N_R)
	array_ones = np.ones(N, dtype=float).reshape(-1,1)
	A = np.concatenate((samples['x'], array_ones), axis=1)
	d = samples['y']
