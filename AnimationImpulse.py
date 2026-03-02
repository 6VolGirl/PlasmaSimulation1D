import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib

matplotlib.use("TkAgg")
data = np.loadtxt("C:\\Users\\6anna\\CLionProjects\\PlasmaSimulator1D\\cmake-build-debug\\ImpulsePlasma.cvs", delimiter=",", skiprows=1)
t = data[:, 0]
x = data[:, 1]
Ez = data[:, 2]

t_vals = np.unique(t)
x_vals = np.unique(x)

Nt = len(t_vals)
Nx = len(x_vals)

# собираем Ez в матрицу (Nt, Nx)
Ez_grid = Ez.reshape(Nt, Nx)

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=1.5)
ax.set_xlim(x_vals[0], x_vals[-1])
ax.set_ylim(Ez_grid.min()*1.1, Ez_grid.max()*1.1)
ax.grid(True)

def update(k):
    print("frame", k)
    line.set_data(x_vals, Ez_grid[k])
    return line,

ani = FuncAnimation(fig, update, frames=Nt, interval=20, blit=False)
plt.show()


# plt.figure()
# plt.plot(x_vals, Ez_grid[550])   # первый момент
# plt.figure()
# plt.plot(x_vals, Ez_grid[1000])  # последний момент
# plt.figure()
# plt.plot(x_vals, Ez_grid[1500])
# plt.figure()
# plt.plot(x_vals, Ez_grid[2000])
# plt.figure()
# plt.plot(x_vals, Ez_grid[2500])
# plt.figure()
# plt.plot(x_vals, Ez_grid[3000])
# plt.show()