import matplotlib.pyplot as plt
import numpy as np

import matplotlib
colors = matplotlib.cm.get_cmap('plasma')

from tqdm import tqdm # progress bar
from scipy.optimize import fsolve # find roots of equation

from consts import (x0, y0, z0, r0, t_init, t_fin, time_step)
from .edo import (fold, hopf)
from .rk4 import rk4
phi = 0
nphi = 100
nguesses = 3
phi_mesh = np.linspace(start=-2, stop=2, num=nphi)

# number of time steps
nt = int((t_fin - t_init) / time_step)
time_mesh = np.linspace(start=t_init, stop=t_fin, num=nt)

#def bifurcation(edo, phi_list, guess):
#    equilibria_mesh = np.zeros((nphi), nguesses)
#    try:



def fold_bifurcation(guesses):
    equilibria_mesh = np.zeros((nphi, nguesses))

    # for each phi
    for phi_index in range(0, nphi-1):

      # find the equilibria of our system

      equilibria = []
      new_guesses = []

      # look for some equilibria
      for guess in guesses:
          equilibrium = fsolve(func=fold, x0=[guess], args=(phi_mesh[phi_index]))
          equilibria.append(equilibrium[0])
          new_guesses.append(equilibrium)

      guesses = new_guesses
      np.array([equilibria])
      # add to the mesh
      #print(np.shape(equilibria_mesh[phi_index]), np.shape(equilibria))

      equilibria_mesh[phi_index] = np.array([equilibria])

    plot(dataset=equilibria_mesh.copy(), ylabel="x")

def fold_hopf_bifurcation():
    equilibria_mesh = np.zeros((nphi, nguesses))

    # for each phi
    for phi_index in range(0, nphi-1):
      # create mesh
      # x = np.zeros(self.nt)
      # # set initial conditions (for each phi value)
      # x[0] = x0

      # # simulate the system over time
      # for t in range(0, self.nt-1):
      #   x[t+1] = rk4(fold, time_step, t, x, self.phi_mesh[phi_index])
      # jacobian matrix df(p)
      #jac =

      guesses = np.linspace(start=-3, stop=3, num=nguesses)
      equilibria = []

      # look for some equilibria
      for guess in guesses:
          equilibrium_fold = fsolve(func=fold, x0=[guess], args=(phi_mesh[phi_index]))
          equilibrium_hopf_coupled = fsolve(func=hopf_coupled, x0=[equilibrium_fold, guess, guess], args=(gamma, phi_mesh[phi_index]))
          equilibria.append(equilibrium_fold[0])
          equilibria.append(equilibrium_hopf_coupled[0])

      # add to the mesh
      equilibria_mesh[phi_index] = equilibria

    plot(equilibria_mesh.copy())

def plot(dataset, ylabel):
    fig, ax = plt.subplots(figsize=(15, 7))
    ax.plot(phi_mesh, dataset)

    plt.xlabel("$\phi$")
    plt.ylabel(ylabel)
    plt.xlim(-2,2)
    plt.ylim(-5,5)
    plt.show()

guesses = np.linspace(start=-3, stop=3, num=nguesses)

