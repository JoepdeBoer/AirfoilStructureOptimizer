from sko.PSO import PSO
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from structure import Airfoil

def objective(design):
    t = design[0]*0.01  # thickness
    x1 = design[1]  # chord position
    airfoil = Airfoil("NACA012.txt", t, [x1, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    return airfoil.structural_area()


def g_5(design):
    t = design[0] * 0.01 # thickness
    x1 = design[1]  # chord position
    airfoil = Airfoil("NACA012.txt", t, [x1, 0.9], 2.5, [0, 75 * 10 ** 3, 300 * 10 ** 3, 1])
    shear_stress = airfoil.find_max_shear_flow() / t
    return (abs(shear_stress) - 207 * 10 ** 6)/(207 * 10 ** 6)

lb = [0.05, 0.1]
ub = [1, 0.8]

cons = (g_5,)

iter = 50
particals = 15
pso = PSO(func=objective, n_dim=2, pop=particals, max_iter=iter, lb=lb, ub=ub
          , constraint_ueq=cons)
pso.record_mode = True
pso.run()

print('best_x is ', pso.gbest_x, 'best_y is', pso.gbest_y)
print('g_5', g_5(pso.gbest_x))
print('rms V', np.sqrt(np.mean(pso.V**2)))
plt.plot(pso.gbest_y_hist)
plt.xlabel('iteration')
plt.ylabel('objective [m^2]')
plt.show()

# code from example on scikit-opt https://github.com/guofei9987/scikit-opt/blob/master/examples/demo_pso_ani.py
record_value = pso.record_value
X_list, V_list = record_value['X'], record_value['V']

fig, ax = plt.subplots(1, 1)
ax.set_title('PSO', loc='center')
# ax.set_xlabel('iterations')
# ax.set_ylabel('objective [m^2]')
line = ax.plot([], [], 'b.')
print('starting animation')
t_grid, x1_grid = np.meshgrid(np.linspace(0.05, 1, 50), np.linspace(0.1, 0.8, 50))
Z_grid = np.zeros(t_grid.shape)
for j, t in enumerate(t_grid[0, :]):
    for i, x1 in enumerate(x1_grid[0, :]):
        Z_grid[i, j] = objective((t, x1))
ax.contour(t_grid, x1_grid, Z_grid, 30)

ax.set_xlim(0.02, 1.01)
ax.set_ylim(0.05, 0.81)

plt.ion()
p = plt.show()


def update_scatter(frame):
    i, j = frame // 10, frame % 10
    ax.set_title('iter = ' + str(i))
    X_tmp = X_list[i] + V_list[i] * j / 10.0
    plt.setp(line, 'xdata', X_tmp[:, 0], 'ydata', X_tmp[:, 1])
    return line


ani = FuncAnimation(fig, update_scatter, blit=True, interval=25, frames=iter * 10)
plt.show()

ani.save('pso.gif', writer='pillow')
