import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib

matplotlib.use("TkAgg")


i0 = 100
L  = 100
source_pos = 50
resolution = 20
plasma_x0 = (i0 - source_pos) / resolution
plasma_x1 = (i0 + L - source_pos) / resolution

data = np.loadtxt(
    r"C:\Users\6anna\CLionProjects\PlasmaSimulator1D\cmake-build-debug\ImpulsePlasma.cvs",
    delimiter=",",
    skiprows=1
)

t = data[:, 0]
x = data[:, 1]
Ez = data[:, 2]

t_vals = np.unique(t)
x_vals = np.unique(x)

Nt = len(t_vals)
Nx = len(x_vals)

Ez_grid = Ez.reshape(Nt, Nx)

fig, ax = plt.subplots()

# подсветка плазмы
ax.axvspan(plasma_x0, plasma_x1, color="orange", alpha=0.25, zorder=0)  # xmin..xmax в данных [web:103]

line, = ax.plot([], [], lw=1.5, color="C0", zorder=2)

ax.set_xlim(x_vals[0], x_vals[-1])
ax.set_ylim(Ez_grid.min() * 1.1, Ez_grid.max() * 1.1)
ax.grid(True)

def update(k):
    line.set_data(x_vals, Ez_grid[k])
    ax.set_title(f"frame={k}, t={t_vals[k]:.6g}")
    return line,

ani = FuncAnimation(fig, update, frames=Nt, interval=20, blit=False)
plt.show()
