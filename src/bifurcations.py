import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-whitegrid")

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
    hopf_stable_equilibria, hopf_unstable_equilibria = bifurcation(hopf_polar, hopf_polar_df)

    # fold-hopf
    fold_coupled_stable_equilibria, fold_coupled_unstable_equilibria, hopf_coupled_stable_equilibria, hopf_coupled_unstable_equilibria = coupled_bifurcation(fold, fold_df, hopf_polar_coupled, hopf_polar_coupled_df)

    fig, (ax1, ax2, ax3) = plt.subplot(3, 1, figsize=(15, 7))
    fig.suptitle("Diagrammes de bifurcation")

    ax1.plot(fold_stable_equilibria[0], fold_stable_equilibria[1], color="black", label="ligne stable")
    ax1.plot(fold_unstable_equilibria[0], fold_unstable_equilibria[1], linestyle="dashdot", color="red", label="ligne instable")
    ax1.set_xlabel("$\phi$")
    ax1.set_ylabel("$x$")
    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-5, 5)
    ax1.set_title("Fold")
    ax1.set_legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax2.plot(hopf_stable_equilibria[0], hopf_stable_equilibria[1], color="black", label="ligne stable")
    ax2.plot(hopf_unstable_equilibria[0], hopf_unstable_equilibria[1], linestyle="dashdot", color="red", label="ligne instable")
    ax2.set_xlabel("$\phi$")
    ax2.set_ylabel("$r$")
    ax2.set_xlim(-2, 2)
    ax2.set_ylim(-5, 5)
    ax2.set_title("Hopf")
    ax2.set_legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax3.plot(fold_coupled_stable_equilibria[0], fold_coupled_stable_equilibria[1], color="black", label="ligne stable")
    ax3.plot(fold_coupled_unstable_equilibria[0], fold_coupled_unstable_equilibria[1], linestyle="dashdot", color="red", label="ligne instable")
    ax3.plot(hopf_coupled_stable_equilibria[0], hopf_coupled_stable_equilibria[1], color="black", label="ligne stable")
    ax3.plot(hopf_coupled_unstable_equilibria[0], hopf_coupled_unstable_equilibria[1], linestyle="dashdot", color="red", label="ligne instable")
    ax3.set_xlabel("$\phi$")
    ax3.set_ylabel("$x$, $r$")
    ax3.set_xlim(-2, 2)
    ax3.set_ylim(-5, 5)
    ax3.set_title("Fold-Hopf")
    ax3.set_legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.show()
