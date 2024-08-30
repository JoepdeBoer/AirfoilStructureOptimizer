import scipy as sp
import numpy as np
from structure import Airfoil

def objective(design):
    t = design[0]*0.01  # thickness
    x1 = design[1]  # chord position
    x2 = design[2]
    if x2 <= x1:
        return 100
    airfoil = Airfoil("NACA012.txt", t, [x1, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.structural_area()

def grad(design):
    return [0, 0, 0]


def g_5(design):
    t = design[0] * 0.01 # thickness
    x1 = design[1]  # chord position
    x2 = design[2]
    if x2 <= x1:
        return 100
    airfoil = Airfoil("NACA012.txt", t, [x1, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    shear_stress = airfoil.find_max_shear_flow() / t
    return (abs(shear_stress) - 207 * 10 ** 6)/(207 * 10 ** 6)

def g_3(design):
    return design[1] - design[2] + 0.01

x0 = [0.12674015, 0.53824447, 0.54884649]


grad1 = sp.optimize.check_grad(objective, grad, x0=np.array(x0))
grad2 = sp.optimize.check_grad(g_5, grad, x0=np.array(x0))
print(grad1)
print(grad2)