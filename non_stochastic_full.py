from sko.PSO import PSO
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib as mpl
import scipy as sp
import numpy as np
from structure import Airfoil

mpl.rcParams.update({'font.size': 20})  # Set the font size to 14 for all elements


def objective(design):
    t = design[0]*0.01  # thickness
    x1 = design[1]  # chord position
    x2 = design[2]
    if x2 <= x1:
        return 100
    airfoil = Airfoil("NACA012.txt", t, [x1, x2], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.structural_area()

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

def grad(design):
    return np.array([0, 0, 0])


lb = [0.05, 0.1, 0.11]
ub = [1, 0.8, 0.95]


bounds = sp.optimize.Bounds(lb, ub)

# nonlinear constraints
shearstrenght = 207 * 10 ** 6 # [N/m] Aluminum 6061-T6; 6061-T651 https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma6061t6
nl_lb = - np.inf
nl_ub = 0
con1 = sp.optimize.NonlinearConstraint(g_5, nl_lb, nl_ub)
con2 = sp.optimize.NonlinearConstraint(g_3, nl_lb, nl_ub)
cons =(con1, con2)
res = sp.optimize.minimize(objective, [0.12674015, 0.53824447, 0.54884649], method='SLSQP', bounds=bounds, constraints=cons ,tol=1e-10, options={'disp': False, 'maxiter': 10000})
#res3 = sp.optimize.minimize(objective, [0.2, 0.594], method='SLSQP', bounds=bounds, constraints=cons, options={'disp': False, 'maxiter': 10000}, tol=1e-10)
# res4 = sp.optimize.minimize(objective, [0.4, 0.3], method='Nelder-Mead', bounds=bounds, constraints=cons, options={'disp': True}, tol=1e-7)
print(f'slsqp result: {res}\n')
print(f'g_5: {g_5(res.x)}')
grad1 = sp.optimize.check_grad(objective, grad, x0=np.array(res.x))
grad2 = sp.optimize.check_grad(g_5, grad, x0=np.array(res.x))
print(grad1)
print(grad2)