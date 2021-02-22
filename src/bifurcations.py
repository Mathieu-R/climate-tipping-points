import matplotlib.pyplot as plt
import numpy as np

import matplotlib
colors = matplotlib.cm.get_cmap('plasma')

from tqdm import tqdm # progress bar
from scipy.optimize import fsolve # find roots of equation

from constantss import (x0, y0, z0, r0, t_init, t_fin, time_step)
from .edo import (fold, hopf)
from .rk4 import rk4

class bifurcations():
  def __init__(self):
    self.phi = 0
    self.nphi = 100
    self.phi_mesh = np.linspace(start=-2, stop=2, num=self.nphi)
    # number of time steps
    self.nt = int((t_fin - t_init) / time_step)
    self.time_mesh = np.linspace(start=t_init, stop=t_fin, num=self.nt)

  def fold(self):
    equilibria_mesh = np.zeros(self.nphi)

    # for each phi
    for phi_index in tqdm(range(0, self.nphi)):
      # # create mesh
      # x = np.zeros(self.nt)
      # # set initial conditions (for each phi value)
      # x[0] = x0

      # # simulate the system over time
      # for t in range(0, self.nt-1):
      #   x[t+1] = rk4(fold, time_step, t, x, self.phi_mesh[phi_index])

      # find the equilibria of our system
      equilibria = fsolve(func=fold, x0=np.array([0]), args=(self.phi_mesh[phi_index]))

      # add to the mesh
      equilibria_mesh[phi_index] = equilibria

    self.plot(equilibria_mesh)

  def hopf(self):
    equilibria_mesh = np.zeros(self.nphi)

    # for each phi
    for phi_index in tqdm(range(0, self.nphi)):
      # create mesh
      # x = np.zeros(self.nt)
      # # set initial conditions (for each phi value)
      # x[0] = x0

      # # simulate the system over time
      # for t in range(0, self.nt-1):
      #   x[t+1] = rk4(fold, time_step, t, x, self.phi_mesh[phi_index])

      # find the equilibria of our system
      equilibria = fsolve(func=hopf, x0=np.array([1]), args=(self.phi_mesh[phi_index]))

      # add to the mesh
      equilibria_mesh[phi_index] = equilibria

    self.plot(equilibria_mesh)

  def plot(self, dataset):
    plt.plot(self.phi_mesh, dataset)
    plt.xlabel("$\phi$")
    plt.ylabel("x")
    plt.xlim(-2,2)
    plt.ylim(-5,5)
    plt.show()
