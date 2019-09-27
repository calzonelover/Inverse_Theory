import math
import random
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

HOUSE = 20
T = 300
T_DECAY = 0.9
N_ACCEPTED_MAX = 10

KEEP_DISTANCE = []

REGION = np.array([
    [np.random.uniform(100,-100), np.random.uniform(100,-100)] for i in range(HOUSE)
])

def compute_path(sequence):
    distance = 0
    for i in range(len(sequence)-1):
        distance += np.sqrt(np.sum(np.square(np.subtract(
            REGION[sequence[i+1]], REGION[sequence[i]]
        ))))
    return distance

def swop_sequence(original_sequence, n_region):
    sequence = list(map(lambda x: x, original_sequence))
    i = np.random.randint(n_region)
    j = np.random.randint(n_region)
    sequence[i], sequence[j] = sequence[j], sequence[i]
    return sequence

if __name__ == "__main__":
    n_region = len(REGION)
    accepted_sequence = random.sample([i for i in range(n_region)], n_region)
    ##
    t = T
    n = 0
    distance_old = 100
    while t > 1e-3:
        while n < N_ACCEPTED_MAX:
            # x <- x + dx
            sequence = swop_sequence(accepted_sequence, n_region)
            # dE
            distance = compute_path(sequence)
            diff_distance = distance - compute_path(accepted_sequence)
            if diff_distance < 0:
                n += 1
                accepted_sequence = sequence
            else:
                prob = math.e**(-diff_distance/t)
                if np.random.uniform(0,1) < prob:
                    n += 1
                    accepted_sequence = sequence
        KEEP_DISTANCE.append(compute_path(accepted_sequence))
        n = 0
        t *= T_DECAY
    print("End t = {}".format(t))
    # visulize
    plt.plot(REGION[:, 0], REGION[:, 1], 'o')

    path_line = np.array([REGION[i] for i in accepted_sequence])
    plt.plot(path_line[:, 0], path_line[:, 1], '.-')

    plt.plot(REGION[accepted_sequence[0], 0], REGION[accepted_sequence[0], 1], 'ro', label="start")
    plt.plot(REGION[accepted_sequence[-1], 0], REGION[accepted_sequence[-1], 1], 'go', label="end")
    plt.legend()

    plt.title('Saleman walking path')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.savefig('walking_path.png')
    # plt.show()

    plt.clf()
    plt.plot(KEEP_DISTANCE)
    plt.title('Distance over iteration')
    plt.xlabel('k')
    plt.ylabel('Distance')
    plt.savefig('distance.png')
    # plt.show()