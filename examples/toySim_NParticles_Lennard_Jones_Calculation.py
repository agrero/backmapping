from sympy import *
# Import necessary libraries

init_printing()

from IPython.display import display, HTML

display(HTML("<style>div.output_scroll { height: 60em; }</style>"))

import sympy
from sympy import *
from sympy import lambdify

import numpy as np
import math as math
import matplotlib.pyplot as plt
import time

import pandas as pd


# record start time
start = time.time()

N = 6;
epsilon = 1;
sigma = 2;

# Define Variables

k_list = []
z_list = []
m_list = []

n = 0
while n <= N - 1:
    k_list.append(6)
    z_list.append(1)
    m_list.append(2)
    n += 1

k = np.array(k_list)
z = np.array(z_list)
m = np.array(m_list)

# k = np.array([6, 6, 14])
# z = np.array([1, 1, 2])
# m = np.array([1, 1, 1])

kbT = 1

T = 1
delta = 0.1


####################################################################
def simulateTrimer_gradient(k, z, m, kbT, T, delta):
    steps = int(T / delta)

    #     N = 3;

    zetas = z;
    masses = m;

    zeta3 = np.kron(np.diag(zetas[:]), np.eye(3))
    masses3 = np.kron(np.diag(masses[:]), np.eye(3))

    leq = 1.0;
    th0 = (math.pi) / 2;
    #     phi0 = (math.pi) / 2;
    phi0 = 0;

    a = np.zeros((6 * N, steps + 1))

    ainit = np.zeros(6 * N)

    # Define Initial Coordinates
    n = 0
    while n <= N - 1:
        ainit[3 * n:3 * (n + 1)] = np.array([0.25 * (-1) ** (n + 2), 1 * n + 0.1 * (-1) ** (n + 2), 0])
        n += 1

    #     ainit[0:3] = np.array([0, 0.8, 0])
    #     ainit[3:6] = np.array([0, 0, 0])
    #     ainit[6:9] = np.array([1.2, 0, 0])
    #     ainit[9:12] = np.array([1, 0, 0.8])

    tinit = 0.

    times = np.linspace(tinit, tinit + T, steps + 1)

    a[:, 0] = ainit[:]

    #     print("a={}".format(a))

    sd = np.kron(np.sqrt(2 * kbT * zetas[:] * delta), np.ones(N))

    f = np.zeros((3 * N, steps))

    for i in range(0, 3 * N):
        f[i, :] = np.random.normal(0, sd[i], (1, steps))

    # Perform simulation

    Minv = np.linalg.inv(masses3[:, :])
    # fricMat can no longer be diagonal matrix (row number restricted to 3)
    fricMat = zeta3[:, :] @ Minv

    #     fricMat = np.zeros((3, N))

    #     fricMat[0] = np.array([1, 0, 0, 0])
    #     fricMat[1] = np.array([0, 1, 0, 0])
    #     fricMat[2] = np.array([0, 0, 1, 0])

    #     print(fricMat.size)

    potForce = np.zeros((3 * N, steps + 1))
    LJForce = np.zeros((3 * N, steps + 1))

    # Set Up Gradient
    x = [sympy.symbols('x%d' % i) for i in range(3 * N)]

    #     x1,y1,z1,x2,y2,z2,x3,y3,z3 = symbols("x1,y1,z1,x2,y2,z2,x3,y3,z3")

    # Define Vectors to Atoms
    r_vector_list = []
    n = 0
    while n <= N - 1:
        r_vector = Matrix([[x[3 * n]], [x[3 * n + 1]], [x[3 * n + 2]]])
        r_vector_list.append(r_vector)
        n += 1

    #     print("r_vector_list={}".format(r_vector_list))

    #     r1 = Matrix([[x1],[y1],[z1]])
    #     r2 = Matrix([[x2],[y2],[z2]])
    #     r3 = Matrix([[x3],[y3],[z3]])

    # Define Vectors along Bonds
    distance_vector_list = []
    n = 0
    while n <= N - 2:
        distance_vector_x = -(r_vector_list[n][0] - r_vector_list[n + 1][0])
        distance_vector_y = -(r_vector_list[n][1] - r_vector_list[n + 1][1])
        distance_vector_z = -(r_vector_list[n][2] - r_vector_list[n + 1][2])
        distance_vector_list.append(Matrix([distance_vector_x, distance_vector_y, distance_vector_z]))
        n += 1

    #     print("distance_vector_list={}".format(distance_vector_list))

    #     r12 = r2-r1
    #     r23 = r3-r2

    # Define Magnitudes of Vectors along Bonds
    magnitude_distance_vector_list = []
    n = 0
    while n <= N - 2:
        magnitude_distance_vector = sqrt(distance_vector_list[n].dot(distance_vector_list[n]))
        magnitude_distance_vector_list.append(magnitude_distance_vector)
        n += 1

    #     print("magnitude_distance_vector_list={}".format(magnitude_distance_vector_list))

    #     rs12 = sqrt(r12.dot(r12))
    #     rs23 = sqrt(r23.dot(r23))

    # Define Angles
    theta_list = []
    n = 0
    while n <= N - 3:
        dot_r1_r2 = distance_vector_list[n].dot(distance_vector_list[n + 1])
        #         print("dot_r1_r2={}".format(dot_r1_r2))
        cos_theta = dot_r1_r2 / (magnitude_distance_vector_list[n] * magnitude_distance_vector_list[n + 1])
        theta = acos(cos_theta)
        theta_list.append(theta)
        n += 1

    #     print("theta_list={}".format(theta_list))

    #     cos_theta_123 = -r12.dot(r23)/(rs12*rs23)
    #     theta_123 = acos(cos_theta_123)

    # Define Potential Energy
    U = 0

    # Add Spring Potential Energy
    n = 0
    while n <= N - 2:
        U += (1 / 2) * k[0] * (magnitude_distance_vector_list[n] - leq) ** 2
        n += 1

    # Add Bend Angle Potential Energy
    n = 0
    while n <= N - 3:
        U += (1 / 2) * k[0] * (theta_list[n] - th0) ** 2
        n += 1

    U_LJ = 0
    # Add Lennard Jones Potential Energy
    n = 0
    while n <= N - 5:
        r_n_n4 = -(r_vector_list[n] - r_vector_list[n + 4])
        magnitude_r_n_n4 = sqrt(r_n_n4.dot(r_n_n4))
        U += 4 * epsilon * ((sigma / magnitude_r_n_n4) ** 12 - (sigma / magnitude_r_n_n4) ** 6)
        U_LJ += 4 * epsilon * ((sigma / magnitude_r_n_n4) ** 12 - (sigma / magnitude_r_n_n4) ** 6)
        n += 1

    # Write Parameters to Text File
    file1 = open("Parameters.txt", "w")
    file1.write("toySim_NParticle_Lennard_Jones_Calculation Parameters\n")
    file1.write("\n")
    file1.write("N={}\n".format(N))
    file1.write("epsilon={}\n".format(epsilon))
    file1.write("sigma={}\n".format(sigma))
    file1.write("k={}\n".format(k))
    file1.write("z={}\n".format(z))
    file1.write("m={}\n".format(m))
    file1.write("kbT={}\n".format(kbT))
    file1.write("T={}\n".format(T))
    file1.write("delta={}\n".format(delta))
    file1.write("leq={}\n".format(leq))
    file1.write("th0={}\n".format(th0))
    file1.write("U={}\n".format(U))
    #     file1.write("phi0={}\n".format(phi0))
    file1.close()  # to change file access modes

    # Define Potential Energy
    # No Dihedral Potential
    #     U = (1/2)*k12*(rs12 - l12)**2 + (1/2)*k23*(rs23 - l23)**2 + (1/2)*k34*(rs34 - l34)**2 + (1/2)*kt*(theta_123 - th0)**2 + (1/2)*kt*(theta_234 - th0)**2
    # No Bend and No Dihedral Potential
    #     U = (1/2)*k12*(rs12 - l12)**2 + (1/2)*k23*(rs23 - l23)**2 + (1/2)*k34*(rs34 - l34)**2
    # No Dihedral Potential (3 Particle)
    #     U = (1/2)*k12*(rs12 - l12)**2 + (1/2)*k23*(rs23 - l23)**2 + (1/2)*kt*(theta_123 - th0)**2
    #     #No Bend and No Dihedral Potential (3 Particle)
    #     U = (1/2)*k12*(rs12 - l12)**2 + (1/2)*k23*(rs23 - l23)**2

    #     print("U=")
    #     display(U)

    #     print("U_LJ=")
    #     display(U_LJ)

    coord_list = x
    print(coord_list)
    #     coord_list = [x1,y1,z1,x2,y2,z2,x3,y3,z3]
    force_list = []
    force_LJ_list = []

    # Take gradient of each internal coordinate
    n = 0
    while n <= len(coord_list) - 1:
        #         print(n/len(coord_list))
        F_coord = -U.diff(coord_list[n])
        F_LJ_coord = -U_LJ.diff(coord_list[n])
        #         print("F_{}=".format(coord_list[n]))
        #         display(-Derivative(U,coord_list[n]))
        #         print("F_{}=".format(coord_list[n]))
        #         display(F_coord)
        force_list.append(F_coord)
        force_LJ_list.append(F_LJ_coord)
        n += 1

    #     print("force_LJ_list[0]={}".format(force_LJ_list[0]))
    #     print("force_LJ_list[3]={}".format(force_LJ_list[3]))

    # Time Loop
    for i in range(1, steps + 1):
        print("i/steps={}".format(i / steps))

        #         problem_function = fricMat[:, :] @ a[3*N:6*N, i - 1]
        #         print(problem_function)

        # Get Force From Gradient
        coord_list_num = []
        m = 0
        while m <= len(coord_list) - 1:
            coord_list_num.append(a[m, i - 1])
            m += 1

        force_num_list = []
        force_LJ_num_list = []
        n = 0
        while n <= len(coord_list) - 1:
            F_coord_num = lambdify(coord_list, force_list[n], 'numpy')
            F_LJ_coord_num = lambdify(coord_list, force_LJ_list[n], 'numpy')
            F_coord_num_eval = F_coord_num(*coord_list_num)  # The * puts list in as coordinates in function
            F_LJ_coord_num_eval = F_LJ_coord_num(*coord_list_num)  # The * puts list in as coordinates in function
            force_num_list.append(F_coord_num_eval)
            force_LJ_num_list.append(F_LJ_coord_num_eval)
            n += 1

        #         print("fy1 = {}".format(force_num_list[1]))

        #         phi_1234_num = lambdify([x1,y1,z1,x2,y2,z2,x3,y3,z3, x4, y4, z4], phi_1234, 'numpy')
        #         phi_1234_num_eval = phi_1234_num(*coord_list_num)
        #         print("phi_1234_num_eval={}".format(phi_1234_num_eval))

        potForce[0:3 * N, i - 1] = force_num_list

        # Put LJ Force in Output
        LJForce[0:3 * N, i - 1] = force_LJ_num_list

        # Time Propagation Steps to Determine momentum and position of each particle
        a[3 * N:6 * N, i] = a[3 * N:6 * N, i - 1] + potForce[:, i - 1] * delta - fricMat[:, :] @ a[3 * N:6 * N,
                                                                                                 i - 1] * delta + f[:,
                                                                                                                  i - 1]
        a[0:3 * N, i] = a[0:3 * N, i - 1] + Minv @ a[3 * N:6 * N, i] * delta

    #         #Time Propagation Steps to Determine momentum and position of each particle (Without Random Force)
    #         a[3*N:6*N, i] = a[3*N:6*N, i - 1] + potForce[:, i - 1] * delta - fricMat[:, :] @ a[3*N:6*N, i - 1]
    #         a[0:3*N, i] = a[0:3*N, i - 1] + Minv @ a[3*N:6*N, i] * delta

    #         #Time Propagation Steps to Determine momentum and position of each particle (Without Friction, Without Random Force)
    #         a[3*N:6*N, i] = a[3*N:6*N, i - 1] + potForce[:, i - 1] * delta
    #         a[0:3*N, i] = a[0:3*N, i - 1] + Minv @ a[3*N:6*N, i] * delta

    #         Time Propagation Steps to Determine momentum and position of each particle (Without Friction, all masses = 1)
    #         a[3*N:6*N, i] = a[3*N:6*N, i - 1] + potForce[:, i - 1] * delta + f[:, i - 1]
    #         a[0:3*N, i] = a[0:3*N, i - 1] + a[3*N:6*N, i] * delta

    #                 #Time Propagation Steps to Determine momentum and position of each particle (Without Friction, all masses = 1, no random force)
    #         a[3*N:6*N, i] = a[3*N:6*N, i - 1] + potForce[:, i - 1] * delta
    #         a[0:3*N, i] = a[0:3*N, i - 1] + a[3*N:6*N, i] * delta

    # Store Data in Excel Sheet

    t_list = []
    time_end = T
    n = 0
    while n <= time_end:
        t_list.append(n)
        n += delta

    # create an Empty DataFrame object
    df = pd.DataFrame()

    df['t'] = t_list

    # append positions to dataframe
    n = 0
    while n <= len(coord_list) - 1:
        df['{}'.format(coord_list[n])] = a[n, 0:len(t_list)]
        n += 1

        # append momenta to dataframe
    n = 0
    while n <= len(coord_list) - 1:
        df['p{}'.format(coord_list[n])] = a[3 * N + n, 0:len(t_list)]
        n += 1

    # append LJ Forces to dataframe
    n = 0
    while n <= len(coord_list) - 1:
        df['LJ_Force_{}'.format(coord_list[n])] = LJForce[n, 0:len(t_list)]
        n += 1

        # append Potential Forces to dataframe
    n = 0
    while n <= len(coord_list) - 1:
        df['Pot_Force_{}'.format(coord_list[n])] = potForce[n, 0:len(t_list)]
        n += 1

    # saving the excel
    df.to_excel('Data.xlsx')
    print('DataFrame is written to Excel File successfully.')

    return a, LJForce


A = simulateTrimer_gradient(k, z, m, kbT, T, delta)
a = A[0]
# print(a)

# record end time
end = time.time()

# Append-adds at last
file1 = open("Parameters.txt", "a")  # append mode
file1.write("run_time = {} seconds \n".format(end - start))
file1.close()
