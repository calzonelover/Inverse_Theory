import math
import random
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

T = 300
T_DECAY = 0.95
T_STOP = 1e-2
N_ACCEPTED_MAX = 30

START_POSITION = [20, -20]
XMIN, XMAX = -40, 40
YMIN, YMAX = -40, 40
A = 10
B = 0.2
C = math.pi/5

# log
KEEP_DISTANCE = []
KEEP_POSITION_X1 = [START_POSITION[0], ]
KEEP_POSITION_X2 = [START_POSITION[1], ]

def ackley(x, a=A, b=B, c=C):
    x1 = x[0]
    x2 = x[1]
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
    accepted_position = START_POSITION
    old_potential = ackley(accepted_position)
    ##
    t = T
    n = 0
    while t > T_STOP:
        while n < N_ACCEPTED_MAX:
            # x <- x + dx
            new_position = accepted_position + np.random.normal(loc=0.0, scale=1, size=2)
            new_position = check_bc(new_position)
            # dE
            new_potential = ackley(new_position)
            diff_potential = new_potential - old_potential
            if diff_potential < 0:
                n += 1
                accepted_position = new_position
                old_potential = ackley(accepted_position)
                KEEP_DISTANCE.append(new_potential)
                KEEP_POSITION_X1.append(accepted_position[0])
                KEEP_POSITION_X2.append(accepted_position[1])
            else:
                prob = math.e**(-1.0*diff_potential/t)
                if np.random.uniform(0,1) < prob:
                    n += 1
                    accepted_position = new_position
                    old_potential = ackley(accepted_position)
                    KEEP_DISTANCE.append(new_potential)
                    KEEP_POSITION_X1.append(accepted_position[0])
                    KEEP_POSITION_X2.append(accepted_position[1])
        t *= T_DECAY
        n = 0
    print("End t = {}".format(t))
    # visulize
    x, y = np.meshgrid(np.linspace(XMIN, XMAX, 1000),
                   np.linspace(YMIN, YMAX, 1000))
    ackley_potential = ackley([x, y])
    plt.imshow(
        ackley_potential, cmap='viridis',
        extent=[XMIN, XMAX, YMIN, YMAX],
    )
    plt.colorbar()
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title("Ackley Potential")
    plt.savefig("contour.png")
    plt.plot(KEEP_POSITION_X1, KEEP_POSITION_X2, 'c.-')
    plt.plot(
        KEEP_POSITION_X1[0], KEEP_POSITION_X2[0],
        'ro' ,label='begin ({:.2f}, {:.2f})'.format(KEEP_POSITION_X1[0],
        KEEP_POSITION_X2[0])
    )
    plt.plot(
        KEEP_POSITION_X1[-1], KEEP_POSITION_X2[-1],
        'ko' ,label='end ({:.2f}, {:.2f})'.format(KEEP_POSITION_X1[-1],
        KEEP_POSITION_X2[-1])
    )
    plt.legend()
    plt.savefig("walking_contour.png")
    # decay
    plt.clf()
    plt.plot(KEEP_DISTANCE)
    plt.title('Potential over iteration')
    plt.xlabel('k')
    plt.ylabel('Ackley Potential')
    plt.savefig('potential_decay.png')