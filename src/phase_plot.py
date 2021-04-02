import numpy as np
import matplotlib.pyplot as plt

from consts import (b1, b2, c1, c2, t_init, t_fin, time_step, x0, y0, z0)
from .rk4 import rk4_vec

from src.utils.utils import (set_size, set_size_square_plot)

plt.style.use("science")
plt.rcParams.update({
  "font.size": 10,
  "xtick.labelsize": 8,
  "ytick.labelsize": 8
})

def hopf(t, v, gam):
  return np.array([
    b1*v[1] + b2*(gam - (v[0]**2 + v[1]**2))*v[0],
    c1*v[0] + c2*(gam - (v[0]**2 + v[1]**2))*v[1]
  ])

def solve(solver, edo, initial, dt, nt, gam):
  # vector mesh -- will receive the result
  v_mesh = np.ones((nt, 2)) #
  # set inital conditions
  v_mesh[0] = initial

  for t in range(0, nt - 1):
    v_mesh[t + 1] = v_mesh[t] + solver(edo, v_mesh[t], t * dt, dt, gam)

  return v_mesh

def phase_plot():
  dt = time_step
  nt = int((t_fin - t_init) / dt)
  gam = 1.

  initial = [y0, z0]

  fig, ax = plt.subplots(1, 1, figsize=(set_size_square_plot(width="column-size")))

  y_mesh = np.linspace(start=-3, stop=3, num=100)
  z_mesh = np.linspace(start=-3, stop=3, num=100)

  # meshgrid useful to evaluate function on a grid
  Y, Z = np.meshgrid(y_mesh, z_mesh)

  f = hopf(0, [Y, Z], gam)
  ax.streamplot(Y, Z, f[0,:], f[1,:], color="#BBBBBB")

  results = solve(rk4_vec, hopf, initial, dt, nt, gam)

  ax.plot(results[:,0], results[:,1], 'blue')
  ax.set_title("Portrait de phase --- $\gamma = 1$")
  ax.set_xlim(-3, 3)
  ax.set_ylim(-3, 3)
  ax.set_xlabel('x')
  ax.set_ylabel('y')

  #plt.savefig("article/figures/phase-plot.pdf", dpi=300)
  plt.show()
