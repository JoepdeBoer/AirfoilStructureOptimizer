import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from structure import Airfoil


def objective(design):
    t = design[0]*0.01  # thickness
    x1 = design[1]  # chord position
    airfoil = Airfoil("NACA012.txt", t, [x1, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.structural_area()

# def constraint(design):
#     t = design[0] * 0.01  # thickness
#     x1 = design[1]  # chord position
#     airfoil = Airfoil("NACA012.txt", t, [x1, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
#     shear_stress = airfoil.find_max_shear_flow()/t
#     return shear_stress

def g_5(design):
    t = design[0] * 0.01 # thickness
    x1 = design[1]  # chord position
    airfoil = Airfoil("NACA012.txt", t, [x1, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    shear_stress = airfoil.find_max_shear_flow() / t
    return (abs(shear_stress) - 207 * 10 ** 6)/(207 * 10 ** 6)


# Bounds
lb = np.array([0.05, 0.1]) # lower bound 0.5 mm thickness and 10% chord
ub = np.array([1, 0.8]) # upper bound 1 cm thickness and 80% chord
bounds = sp.optimize.Bounds(lb, ub)

# nonlinear constraints
shearstrenght = 207 * 10 ** 6 # [N/m] Aluminum 6061-T6; 6061-T651 https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma6061t6
nl_lb = - np.inf
nl_ub = 0
cons = sp.optimize.NonlinearConstraint(g_5, nl_lb, nl_ub, keep_feasible=False)

res = sp.optimize.minimize(objective, [0.2, 0.3, 0.5], method='SLSQP', bounds=bounds, constraints=cons ,tol=1e-10, options={'disp': False, 'maxiter': 10000})
#res3 = sp.optimize.minimize(objective, [0.2, 0.594], method='SLSQP', bounds=bounds, constraints=cons, options={'disp': False, 'maxiter': 10000}, tol=1e-10)
# res4 = sp.optimize.minimize(objective, [0.4, 0.3], method='Nelder-Mead', bounds=bounds, constraints=cons, options={'disp': True}, tol=1e-7)
print(f'slsqp result: {res}\n')
print(f'g_5: {g_5(res.x)}')




