import math
import random
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

HOUSE = 6
T = 30000
T_DECAY = 0.9
N_ACCEPTED_MAX = 10

REGION = np.array([
    [np.random.uniform(100,-100), np.random.uniform(100,-100)] for i in range(HOUSE)
])

def compute_path(sequence):
    distance = 0
    for i in range(len(sequence)-1):
        distance += np.sqrt(np.sum(np.square(np.subtract(REGION[sequence[i+1]], REGION[sequence[i]]))))
    return distance

def swop_sequence(original_sequence, n_region):
    sequence = original_sequence
    i = np.random.randint(n_region)
    j = np.random.randint(n_region)
    sequence[i], sequence[j] = sequence[j], sequence[i]
    print(i, j)
    print(sequence, original_sequence)
    return sequence

if __name__ == "__main__":
    n_region = len(REGION)
    accepted_sequence = random.sample([i for i in range(n_region)], n_region)
    ##
    t = T
    n = 0
    distance_old = 100
    while t > 1.0:
        # print("! loop !")
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
                # print(diff_distance)
                if np.random.uniform(0,1) < prob:
                    n += 1
                    accepted_sequence = sequence
        # print(accepted_sequence)
        # print(compute_path(accepted_sequence))
        # print("! end loop !")
        n = 0
        t *= T_DECAY
        # print("!", distance)
    print("End t = {}".format(t))
    # visulize
    plt.plot(REGION[:, 0], REGION[:, 1], 'o')
    
    # print(accepted_sequence)
    # print(compute_path(accepted_sequence))

    path_line = np.array([REGION[i] for i in accepted_sequence])
    plt.plot(path_line[:, 0], path_line[:, 1], '.-')

    plt.plot(REGION[accepted_sequence[0], 0], REGION[accepted_sequence[0], 1], 'ro', label="start")
    plt.plot(REGION[accepted_sequence[-1], 0], REGION[accepted_sequence[-1], 1], 'go', label="end")
    plt.legend()

    plt.show()