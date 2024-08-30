import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from structure import Airfoil

mpl.rcParams.update({'font.size': 18})
def objective(design):
    t = design[0]  # thickness
    x1 = design[1]  # chord position
    x2 = design[2]  # spar position
    airfoil = Airfoil("NACA012.txt", t, [x1, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.structural_area()


if __name__ == "__main__":
    step = 10**-15
    t0 = 0.0012
    x10 = 0.6
    x20 = 0.95
    thickness = np.arange(t0, t0 + 30 * step, step)
    spar1 = np.arange(x10, x10 + 30 * step, step)
    spar2 = np.arange(x20, x20 + 30 * step, step)
    l1 = []
    l2 = []
    l3 = []


    for i, t in enumerate(thickness):
        airfoil = Airfoil("NACA012.txt", t, [x10, x20], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
        shear_stress = airfoil.find_max_shear_flow() / t
        l1.append((abs(shear_stress) - 207 * 10 ** 6)/(207 * 10 ** 6))

    for i, x1 in enumerate(spar1):
        airfoil = Airfoil("NACA012.txt", t0, [x1, x20], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
        shear_stress = airfoil.find_max_shear_flow() / t0
        l2.append((abs(shear_stress) - 207 * 10 ** 6)/(207 * 10 ** 6))

    for i, x2 in enumerate(spar2):
        airfoil = Airfoil("NACA012.txt", t0, [x10, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
        shear_stress = airfoil.find_max_shear_flow() / t0
        l3.append((abs(shear_stress) - 207 * 10 ** 6)/(207 * 10 ** 6))



    # # Create a normal plot of thickness against l1
    # plt.figure()
    # plt.plot(thickness, l1)
    # plt.xlabel(r'$t$ $[m]$')
    # plt.ylabel(r'$g_{5}$ $[-]$')
    # plt.show()
    #
    # plt.figure()
    # plt.plot(spar1, l2)
    # plt.xlabel(r'$x_{1}$ $[m]$')
    # plt.ylabel(r'$g_{5}$ $[-]$')
    # plt.show()
    #
    # plt.figure()
    # plt.plot(spar1, l2)
    # plt.xlabel(r'$x_{1}$ $[m]$')
    # plt.ylabel(r'$g_{5}$ $[-]$')
    # plt.show()

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    fig.text(0.02, 0.5, r'$g_{5}$', ha='center', va='center', rotation='vertical')

    axs[0].set_xlabel(r'$t [mm]$')
    axs[0].scatter(thickness * 1000, l1)
    axs[0].set_title(f'$x1 = {round(x10,2)}$, $x2 = {round(x20,2)}$')

    axs[1].set_xlabel(r'$x1$ $[-]$')
    axs[1].scatter(spar1, l2)
    axs[1].set_title(f'$t = {round(t0*1000,2)} mm$, $x2 = {round(x20,2)}$')

    axs[2].set_xlabel(r'$x2$ $[-]$')
    axs[2].scatter(spar2, l3)
    axs[2].set_title(f'$t = {round(t0*1000,2)} mm$, $x1 = {round(x10,2)}$')

    plt.show()
