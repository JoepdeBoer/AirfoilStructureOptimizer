from sko.PSO import PSO
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib as mpl
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

lb = [0.05, 0.1, 0.11]
ub = [1, 0.8, 0.95]

cons = (g_5, g_3)

iter = 50
particals = 50
pso = PSO(func=objective, n_dim=3, pop=particals, max_iter=iter, lb=lb, ub=ub
          , constraint_ueq=cons)
pso.record_mode = True
pso.run()

print('best_x is ', pso.gbest_x, 'best_y is', pso.gbest_y)
print('g_5', g_5(pso.gbest_x))
print('g_3', g_3(pso.gbest_x))

plt.plot(pso.gbest_y_hist)
plt.xlabel('iteration')
plt.ylabel('objective [m^2]')
plt.show()

# code from example on scikit-opt https://github.com/guofei9987/scikit-opt/blob/master/examples/demo_pso_ani.py
record_value = pso.record_value
X_list, V_list = record_value['X'], record_value['V']

fig, ax = plt.subplots(1, 1)
ax.set_title('PSO', loc='center')

print('starting animation')
x1_grid, x2_grid = np.meshgrid(np.linspace(0.1, 0.8, 20), np.linspace(0.11, 0.95, 20))
Z_grid = np.zeros(x1_grid.shape)
for j, x1 in enumerate(x1_grid[:, 0]):
    for i, x2 in enumerate(x2_grid[:, 0]):
        Z_grid[i, j] = objective((0.13, x1, x2))
ax.contour(x1_grid, x2_grid, Z_grid, 20*20)
line = ax.plot([], [], 'r.')
ax.set_xlabel(r'$x_{1}$')
ax.set_ylabel(r'$x_{2}$')
# ax.set_xlim(0.009, 0.9)
# ax.set_ylim(0.009, 1)

plt.ion()
p = plt.show()


def update_scatter(frame):
    i, j = frame // 10, frame % 10
    ax.set_title('iter = ' + str(i))
    X_tmp = X_list[i] + V_list[i] * j / 10.0
    plt.setp(line, 'xdata', X_tmp[:, 1], 'ydata', X_tmp[:, 2])
    return line


ani = FuncAnimation(fig, update_scatter, blit=True, interval=25, frames=iter * 10)
plt.show()

ani.save('pos_full_report.gif', writer='pillow')