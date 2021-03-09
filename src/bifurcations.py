import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#from tqdm import tqdm # progress bar
#from scipy.optimize import fsolve # find roots of equation

from consts import (x0, y0, z0, r0, a1, a2, b1, b2, c1, c2, t_init, t_fin, time_step, TRESHOLD)
from .edo import (fold, fold_df, hopf_polar, hopf_polar_df, hopf_polar_coupled, hopf_polar_coupled_df)
from .rk4 import rk4

def bifurcation(fx, df):
  stable_equilibria = []
  unstable_equilibria = []

  nphi = 100
  # bifurcation from phi = -2 to phi = 2
  phi_mesh = np.linspace(start=-2, stop=2, num=nphi)
  # x from x = -2 to x = 2
  x_mesh = np.linspace(start=-2, stop=2, num=1000)

  for phi in phi_mesh:
    for x in x_mesh:
      if -TRESHOLD < fx([x], phi) < TRESHOLD:
        if df([x], phi) > 0:
          unstable_equilibria.append([x, phi])
        elif df([x], phi) < 0:
          stable_equilibria.append([x, phi])

  # sorting equilibria
  unstable_equilibria = list(zip(*sorted(unstable_equilibria, key=lambda x: x[0])))
  stable_equilibria = list(zip(*sorted(stable_equilibria, key=lambda x: x[0])))
  return stable_equilibria, unstable_equilibria

def coupled_bifurcation(fx, df, gx, dg):
  stable_equilibria_f = []
  unstable_equilibria_f = []
  stable_equilibria_g = []
  unstable_equilibria_g = []

  nphi = 100
  # bifurcation from phi = -2 to phi = 2
  phi_mesh = np.linspace(start=-2, stop=2, num=nphi)
  # x from x = -2 to x = 2
  x_mesh = np.linspace(start=-2, stop=2, num=1000)

  for phi in phi_mesh:
    for x in x_mesh:
      if -TRESHOLD < fx(x, phi) < TRESHOLD:
        if df(x, phi) > 0:
          unstable_equilibria_f.append([x, phi])
        elif df(x, phi) < 0:
          stable_equilibria_f.append([x, phi])

      if -TRESHOLD < gx([x, x], phi):
        if df([x, x], phi) > 0:
          unstable_equilibria_g.append([x, phi])
        elif df([x, x], phi) < 0:
          stable_equilibria_g.append([x, phi])

  # sorting equilibria
  unstable_equilibria_f = list(zip(*sorted(unstable_equilibria_f, key=lambda x: x[0])))
  stable_equilibria_f = list(zip(*sorted(stable_equilibria_f, key=lambda x: x[0])))

  unstable_equilibria_g = list(zip(*sorted(unstable_equilibria_g, key=lambda x: x[0])))
  stable_equilibria_g = list(zip(*sorted(stable_equilibria_g, key=lambda x: x[0])))

  return stable_equilibria_f, unstable_equilibria_f, stable_equilibria_g, unstable_equilibria_g

  def bifurcations():
    # fold
    fold_stable_equilibria, fold_unstable_equilibria = bifurcation(fold, fold_df)
    # hopf
    hopf_stable_equilibria,hopf_unstable_equilibria = bifurcation(hopf_polar, hopf_polar_df)

    # fold-hopf
    fold_coupled_stable_equilibria, fold_coupled_unstable_equilibria, hopf_coupled_stable_equilibria, hopf_coupled_unstable_equilibria = bifurcation(fold, fold_df, hopf_polar_coupled, hopf_polar_coupled_df)

    # plt.subplot(2, 1, 1)
    # plt.xlim(-2, 2)
    # plt.ylim(-5, 5)
    # plt.plot(unstable_equ_phi, unstable_equ_x, linestyle="dashdot", label="ligne instable")
    # plt.plot(stable_equ_phi, stable_equ_x, label="ligne stable")
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

