import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from structure import Airfoil


def objective(design):
    t = design[0]  # thickness
    x1 = design[1]  # chord position
    airfoil = Airfoil("NACA012.txt", t, [x1, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.structural_area()

def constraint(design):
    t = design[0]  # thickness
    x1 = design[1]  # chord position
    airfoil = Airfoil("NACA012.txt", t, [x1, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    shear_stress = airfoil.find_max_shear_flow()/t
    return shear_stress

# Bounds
lb = np.array([0.0005, 0.05]) # lower bound 1 mm thickness and 5% chord
ub = np.array([0.1, 0.8]) # upper bound  10 cm thickness and 80% chord
bounds = sp.optimize.Bounds(lb, ub, keep_feasible=np.array([False, False]))

# nonlinear constraints
shearstrenght = 207 * 10 ** 6 # [N/m] Aluminum 6061-T6; 6061-T651 https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma6061t6
nl_lb = - shearstrenght
nl_ub = shearstrenght
cons = sp.optimize.NonlinearConstraint(constraint, nl_lb, nl_ub, keep_feasible=False)

# res = sp.optimize.minimize(objective, [0.005, 0.5], method='COBYLA', bounds=bounds, constraints=cons, options={'disp': True}, )
# print(f'linear result: {res}')
#res2 = sp.optimize.minimize(objective, [0.005, 0.5], method='COBYQA', bounds=bounds, constraints=cons, options={'disp': True}, )
res3 = sp.optimize.minimize(objective, [0.008, 0.5], method='SLSQP', bounds=bounds, constraints=cons, options={'disp': True}, )
print(f'nonlinear result: {res3}')