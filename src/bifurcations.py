import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
colors = matplotlib.cm.get_cmap('plasma')

from tqdm import tqdm # progress bar
from scipy.optimize import fsolve # find roots of equation

fig = plt.figure(figsize=(15,7))
mpl.rc('font', size = 14)

from consts import (x0, y0, z0, r0, a1, a2, b1, b2, c1, c2, t_init, t_fin, time_step, TRESHOLD)
from .edo import (fold, hopf)
from .rk4 import rk4

def fold_bifurcation(guesses):
  phi = 0
  nphi = 100
  nguesses = 3
  phi_mesh = np.linspace(start=-2, stop=2, num=nphi)

  # number of time steps
  nt = int((t_fin - t_init) / time_step)
  time_mesh = np.linspace(start=t_init, stop=t_fin, num=nt)

  unstable_equ = []
  stable_equ = []
  for phi in phi_mesh:
    fx_set = []
    for x in np.linspace(0,2,1000):
      fx = a1 * (x ** 3) + a2 * x + phi
      fx_set.append(fx)
      # confident interval
      if -TRESHOLD < fx < TRESHOLD:
        # compute derivative to check the nature of equilibrium
        dfdx = a1 * ((3 * x) ** 2) + a2
        if dfdx > 0:
          unstable_equ.append([x, phi])
        if dfdx < 0:
          stable_equ.append([x, phi])

    unstable_equ = list(zip(*sorted(unstable_equ, key=lambda x: x[0])))
    unstable_equ_phi = unstable_equ[1]
    unstable_equ_x = unstable_equ[0]

    stable_equ = list(zip(*sorted(stable_equ, key=lambda x: x[0])))
    stable_equ_phi = stable_equ[1]
    stable_equ_x = stable_equ[0]

    plt.subplot(2, 1, 1)
    plt.xlim(-2, 2)
    plt.ylim(-5, 5)
    plt.plot(unstable_equ_phi, unstable_equ_x, linestyle="dashdot", label="ligne instable")
    plt.plot(stable_equ_phi, stable_equ_x, label="ligne stable")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# def fold_hopf_bifurcation():
#     equilibria_mesh = np.zeros((nphi, nguesses))

#     # for each phi
#     for phi_index in range(0, nphi-1):
#       # create mesh
#       # x = np.zeros(self.nt)
#       # # set initial conditions (for each phi value)
#       # x[0] = x0

#       # # simulate the system over time
#       # for t in range(0, self.nt-1):
#       #   x[t+1] = rk4(fold, time_step, t, x, self.phi_mesh[phi_index])
#       # jacobian matrix df(p)
#       #jac =

#       guesses = np.linspace(start=-3, stop=3, num=nguesses)
#       equilibria = []

#       # look for some equilibria
#       for guess in guesses:
#           equilibrium_fold = fsolve(func=fold, x0=[guess], args=(phi_mesh[phi_index]))
#           equilibrium_hopf_coupled = fsolve(func=hopf_coupled, x0=[equilibrium_fold, guess, guess], args=(gamma, phi_mesh[phi_index]))
#           equilibria.append(equilibrium_fold[0])
#           equilibria.append(equilibrium_hopf_coupled[0])

#       # add to the mesh
#       equilibria_mesh[phi_index] = equilibria

#     plot(equilibria_mesh.copy())

# def plot(dataset, ylabel):
#     fig, ax = plt.subplots(figsize=(15, 7))
#     ax.plot(phi_mesh, dataset)

#     plt.xlabel("$\phi$")
#     plt.ylabel(ylabel)
#     plt.xlim(-2,2)
#     plt.ylim(-5,5)
#     plt.show()

guesses = np.linspace(start=-3, stop=3, num=nguesses)

