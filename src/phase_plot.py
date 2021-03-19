import numpy as np
#import matplotlib.pyplot as plt

from consts import (b1, b2, c1, c2, t_init, t_fin, time_step, x0, y0, z0)
from .rk4 import rk4

def hopf(y, z, gam):
  return np.array([
    b1*z + b2*(gam - (y**2 + z**2))*y,
    c1*y + c2*(gam - (y**2 + z**2))*z
  ])

def solve(solver, edo, initial, dt, nt, gam):
  # vector mesh -- will receive the result
  v_mesh = np.ones((nt, 2)) #
  # set inital conditions
  v_mesh[0] = initial

  for t in range(0, nt - 1):
    v_mesh[t + 1] = v_mesh[t] + solver(edo, v_mesh[t], dt, gam)

  return v_mesh

def phase_plot():
  dt = time_step
  nt = int((t_fin - t_init) / dt)
  gam = 1.

  initial = [x0, y0, z0]

  y_mesh = np.linspace(start=-3, stop=3, num=nt)
  z_mesh = np.linspace(start=-3, stop=3, num=nt)

  # meshgrid useful to evaluate function on a grid
  Y, Z = np.meshgrid(y_mesh, z_mesh)
  print(Y, Z)

  f = hopf(Y, Z, gam)
  print(f)
  #plt.streamplot(Y, Z, f[0,:], f[1,:], color="#BBBBBB")

  #results = solve(rk4, hopf, initial, dt, nt, gam)
  #plt.plot(results[0,:], results[1,:], 'blue')
  #plt.title("Portrait de Phase - Hopf Bifurcation")
  #plt.xlabel('x')
  #plt.ylabel('y')
