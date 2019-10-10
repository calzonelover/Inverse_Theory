import math
import random
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

N_PARTICLE = 32
OMEGA = 0.7
C_b = 2.0
C_B = 2.0
STOP_VAL = 1e-4

XMIN, XMAX = -40, 40
YMIN, YMAX = -40, 40
A = 10
B = 0.2
C = math.pi/5

# log
MEAN_REWARDS = []
MAX_REWARDS = []
GLOBAL_BEST_REWARD = None
GLOBAL_BEST_POSITION = []

def ackley(x1, x2, a=A, b=B, c=C):
    first_term = np.multiply(
        np.negative(a),
        np.exp(
            np.multiply(
                np.negative(b)/math.sqrt(2),
                np.sqrt(np.add(
                    np.square(x1),
                    np.square(x2)
                ))
            )
        )
    )
    second_term = np.negative(np.exp(
        np.divide(
            np.add(
                np.cos(np.multiply(c, x1)),
                np.cos(np.multiply(c, x2))
            ),
            2
        )
    ))
    return np.add(
        np.add(
            first_term,
            second_term
        ), 
        a + np.exp(1)
    )

def get_fitness(position):
    global GLOBAL_BEST_REWARD, GLOBAL_BEST_POSITION
    fitness = -ackley(position[:, 0], position[:, 1])
    MAX_REWARDS.append(np.max(fitness))
    MEAN_REWARDS.append(np.mean(fitness))
    if not GLOBAL_BEST_REWARD or GLOBAL_BEST_REWARD < np.max(fitness):
        GLOBAL_BEST_REWARD = np.max(fitness)
        GLOBAL_BEST_POSITION = position[np.argmax(fitness)]
    return fitness

def check_bc(position):
    if position[0] < XMIN:
        position[0] = XMIN
    if position[0] > XMAX:
        position[0] = XMAX
    if position[1] < YMIN:
        position[1] = YMIN
    if position[1] > YMAX:
        position[1] = YMAX
    return position

if __name__ == "__main__":
    current_position = np.random.uniform(low=XMIN, high=XMAX, size=(N_PARTICLE, 2))
    current_velocity = np.zeros(shape=(N_PARTICLE, 2))
    fitness = get_fitness(current_position)
    while MEAN_REWARDS[-1] < -STOP_VAL:
        fitness = get_fitness(current_position)
        best_i_particle = current_position[np.argmax(fitness)]
        new_vecolity = np.add(
            np.multiply(OMEGA, current_velocity),
            np.add(
                np.multiply(
                    np.multiply(C_b, np.random.uniform()),
                    np.subtract(best_i_particle, current_position)
                ),
                np.multiply(
                    np.multiply(C_B, np.random.uniform()),
                    np.subtract(GLOBAL_BEST_POSITION, current_position)
                )
            )
        )
        current_position = np.add(current_position, new_vecolity)
        current_velocity = new_vecolity
        current_position = np.array(list(map(lambda x: check_bc(x), current_position)))
        print(MEAN_REWARDS[-1])
        # visulize
        x, y = np.meshgrid(np.linspace(XMIN, XMAX, 1000),
                    np.linspace(YMIN, YMAX, 1000))
        ackley_potential = ackley(x, y)
        plt.imshow(
            ackley_potential, cmap='viridis',
            extent=[XMIN, XMAX, YMIN, YMAX],
        )
        plt.colorbar()
        plt.xlabel('$x_1$')
        plt.ylabel('$x_2$')
        plt.title("Ackley Potential (k={:04d})".format(len(MAX_REWARDS)))
        plt.plot(
            current_position[:, 0], current_position[:, 1],
            'ro' ,label='Particle'
        )
        plt.legend()
        plt.savefig("log/particle_swarm{:04d}.png".format(len(MAX_REWARDS)))
        plt.clf()
    # decay
    plt.clf()
    plt.plot(np.negative(MAX_REWARDS), label="MAX_REWARDS")
    plt.plot(np.negative(MEAN_REWARDS), label="MEAN_REWARDS")
    plt.legend()
    plt.title('Reward Iteration')
    plt.xlabel('k')
    plt.ylabel('Ackley function value')
    plt.savefig('reward_iteration.png')