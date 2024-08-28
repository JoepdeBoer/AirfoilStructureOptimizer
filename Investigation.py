import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from structure import Airfoil
def objective(design):
    t = design[0]  # thickness
    x1 = design[1]  # chord position
    x2 = design[2]  # spar position
    airfoil = Airfoil("NACA012.txt", t, [x1, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.structural_area()

def lenghtofspar(design):
    t = design[0]  # thickness
    x1 = design[1]  # chord position
    x2 = design[2]  # spar position
    airfoil = Airfoil("NACA012.txt", t, [x1, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.spar1_lenght()

if __name__ == "__main__":
#     points = 50
#     thickness = np.linspace(0.0001, 0.01, points)
#     spar1 = np.linspace(0.05, 0.9, points)
#     spar2 = np.linspace(0.35, 0.9, points)
#     l1 = []
#     l2 = []
#     l3 = []
#     for t1 in thickness:
#         x11 = 0.3
#         x12 = 0.95
#         l1.append(objective([t1, x11, x12]))
#     for x21 in spar1:
#         x22 = 0.95
#         t2 = 0.005
#         l2.append(objective([t2, x21, x22]))
#     for x32 in spar2:
#         x31 = 0.3
#         t3 = 0.005
#         l3.append(objective([t3, x31, x32]))
#
# #subplots Created with Codeium AI Assistant. (Version 1.0) [Software]. (2021). Retrieved from https://www.codeium.com
#     fig, axs = plt.subplots(1, 3, figsize=(15, 5))
#     fig.text(0.02, 0.5, '$A$ $[m^2]$', ha='center', va='center', rotation='vertical', fontsize=12)
#
#     axs[0].set_xlabel(r'$t [m]$', fontsize=12)
#     axs[0].plot(thickness, l1)
#     axs[0].set_title(r'$x1 = 0.3$, $x2 = 0.95$')
#
#     axs[1].set_xlabel(r'$x1$ $[-]$', fontsize=12)
#     axs[1].plot(spar1, l2)
#     axs[1].set_title(r'$t = 0.005$, $x2 = 0.95$')
#
#     axs[2].set_xlabel(r'$x2$ $[-]$', fontsize=12)
#     axs[2].plot(spar2, l3)
#     axs[2].set_title(r'$t = 0.005$, $x1 = 0.3$')
#
#     plt.tight_layout()
#     plt.show()
    # Investigate stress constraint
    points = 30
    thickness = np.linspace(0.003, 0.005, points)
    spar1 = np.linspace(0.05, 0.7, points)
    spar2 = np.linspace(0.55, 0.95, points)
    bigM = np.zeros((points, points))
    l1 = []

    for i, t in enumerate(thickness):
        airfoil = Airfoil("NACA012.txt", t, [0.5, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
        shear_stress = airfoil.find_max_shear_flow() / t
        l1.append(shear_stress - 207 * 10 ** 6)

    # t = 0.005
    # for j, x1 in enumerate(spar1):
    #     for k, x2 in enumerate(spar2):
    #             if x2 > x1 + 0.01:
    #                 airfoil = Airfoil("NACA012.txt", t, [x1, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    #                 shear_stress = airfoil.find_max_shear_flow()/t
    #                 bigM[k, j] = shear_stress - 207 * 10 ** 6
    #             else:
    #                 bigM[k, j] = np.NaN




    # Create a mesh grid for spar1 and spar2
    Spar1, Spar2 = np.meshgrid(spar1, spar2)

    # Create a normal plot of thickness against l1
    plt.figure()
    plt.plot(thickness, l1)
    plt.xlabel('$t$ $[m]$')
    plt.ylabel('$g_{5}$ $[Pa]$')


    # # Create a 3D surface plot of spar1, spar2, and bigM
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(Spar1, Spar2, bigM)
    # ax.set_xlabel(r'$x_{1}$ $[-]$')
    # ax.set_ylabel('$x_{2}$ $[-]$')
    # ax.set_zlabel(r'$g_{5}$ $[Pa]$')


    plt.show()



