import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-whitegrid")

#from tqdm import tqdm # progress bar
#from scipy.optimize import fsolve # find roots of equation

from consts import (x0, y0, z0, r0, a1, a2, b1, b2, c1, c2, t_init, t_fin, time_step, TRESHOLD)
from .edo import (fold, fold_df, hopf_polar, hopf_polar_df, hopf_polar_coupled, hopf_polar_coupled_df)
from .rk4 import rk4

np.set_printoptions(precision=3, suppress=True)

"""Leading vs forcing \phi

"""
def fold_bifurcation(fx, df, ax):
  stable_equ_up = []
  stable_equ_down = []
  unstable_equ = []

  nphi = 1000
  nx = 10000

  phi_mesh = np.linspace(start=-1, stop=1, num=nphi)
  x_mesh = np.linspace(start=-1.5, stop=1.5, num=nx)

  # algorithm
  for phi in phi_mesh:
    for x in x_mesh:
      f_value = fx(x, phi)
      #print(f_value)
      # check if I am in trusting interval (i.e. near an equilibirum point such that \dot{x} = 0)
      if - TRESHOLD < f_value < TRESHOLD:
        # YEAH ! We hit an equilibrium point.
        # Check its stability (stable or unstable)
        df_value = df(x, phi)
        if df_value > 0:
          unstable_equ.append([phi, x])
        elif df_value < 0:
          if x > 0.:
            stable_equ_up.append([phi, x])
          elif x < 0.:
            stable_equ_down.append([phi, x])

  # [[x1, phi1],...,[xN, phiN]] => [[x1,...,xN], [phi1,...,phiN]]
  unstable_equ = list(zip(*unstable_equ))
  stable_equ_up = list(zip(*stable_equ_up))
  stable_equ_down = list(zip(*stable_equ_down))

  # plot
  ax.plot(stable_equ_up[0], stable_equ_up[1], color="black", label="ligne stable")
  ax.plot(stable_equ_down[0], stable_equ_down[1], color="black", label="ligne stable")
  ax.plot(unstable_equ[0], unstable_equ[1], linestyle="dashdot", color="red", label="ligne instable")
  ax.set_xlabel("$\phi$")
  ax.set_ylabel("$x$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("Fold")
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

"""Following vs coupling \gamma

"""
def hopf_bifurcation(fx, df, ax):
  stable_equ_middle = []
  stable_equ_up = []
  stable_equ_down = []
  unstable_equ = []

  ngamma = 1000
  nr = 10000

  gamma_mesh = np.linspace(start=-1, stop=1, num=ngamma)
  r_mesh = np.linspace(start=-1.5, stop=1.5, num=nr)

  # algorithm
  for gamma in gamma_mesh:
    for r in r_mesh:
      f_value = fx(r, gamma)
      #print(f_value)
      # check if I am in trusting interval (i.e. near an equilibirum point such that \dot{r} = 0)
      if - TRESHOLD < f_value < TRESHOLD:
        # YEAH ! We hit an equilibrium point.
        # Check its stability (stable or unstable)
        df_value = df(r, gamma)
        if df_value > 0:
          unstable_equ.append([gamma, r])
        elif df_value < 0:
          # decorralate points (otherwise matplotlib connect every points of each branch with a line)
          if r > 0.:
            stable_equ_up.append([gamma, r])
          elif r < 0.:
            stable_equ_down.append([gamma, r])
          else:
            stable_equ_middle.append([gamma, r])

  unstable_equ = list(zip(*unstable_equ))
  stable_equ_middle = list(zip(*stable_equ_middle))
  stable_equ_up = list(zip(*stable_equ_up))
  stable_equ_down = list(zip(*stable_equ_down))

  # plot
  #ax2.plot(stable_equ_middle[0], stable_equ_middle[0], color="black", label="ligne stable")
  ax.plot(stable_equ_up[0], stable_equ_up[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(stable_equ_down[0], stable_equ_down[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(unstable_equ[0], unstable_equ[1], linestyle="dashdot", color="red", label="ligne instable")
  ax.set_xlabel("$\gamma$")
  ax.set_ylabel("$r$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("Hopf")
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

"""Following vs forcing \phi
  # 1. loop each \phi
  # # 1.1. check roots of the leading system \dot{x} for that \phi => get x values
  # # 1.2. loop each r
  # # # 1.2.1 check roots of the "following system" (polar)
  # # # 1.2.1 that value of x will change value of \gamma in the "following system"

  # sorry for that ugly notation.
  :px: function -- primary system ODE
  :dp: function -- primary system ODE JAC
  :sr: function -- secondary system ODE
  :ds: function -- secondary system ODE JAC
  :ax: Axes -- axe to plot
"""
def fold_hopf_bifurcations(px, dp, sr, ds, ax):
  stable_equ_middle = []
  stable_equ_up = []
  stable_equ_down = []
  unstable_equ = []

  nphi = 1000
  nr = 10000
  nx = 10000

  phi_mesh = np.linspace(start=-1, stop=1, num=nphi)
  x_mesh = np.linspace(start=-1.5, stop=1.5, num=nx)

  r_mesh = np.linspace(start=-1.5, stop=1.5, num=nr)
  # loop each \phi
  for phi in phi_mesh:
    for x in x_mesh:
      p_value = px(x, phi)
      # looking for equilibria in the leading system \dot{x}
      # check if I am in trusting interval (i.e. near an equilibirum point such that \dot{x} = 0)
      if - TRESHOLD < p_value < TRESHOLD:
        # YEAH ! We hit an equilibrium point.
        # loop for each r
        for r in r_mesh:
          s_value = sr(r, x)
          # looking for equilibria in the secondary system \dot{r}
          if - TRESHOLD < s_value < TRESHOLD:
            # YEAH WE GOT IT.
            # check the stability
            ds_value = ds(r, x)
            if ds_value > 0:
              unstable_equ.append([phi, r])
            elif ds_value < 0:
              if r > 0:
                stable_equ_up.append([phi, r])
              elif r < 0:
                stable_equ_down.append([phi, r])
              else:
                stable_equ_middle.append([phi, r])

  unstable_equ = list(zip(*unstable_equ))
  stable_equ_middle = list(zip(*stable_equ_middle))
  stable_equ_up = list(zip(*stable_equ_up))
  stable_equ_down = list(zip(*stable_equ_down))

  # plot
  #ax3.plot(stable_equ_middle[0], stable_equ_middle[0], color="black", label="ligne stable")
  ax.plot(stable_equ_up[0], stable_equ_up[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(stable_equ_down[0], stable_equ_down[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(unstable_equ[0], unstable_equ[1], linestyle="dashdot", color="red", label="ligne instable")
  ax.set_xlabel("$\phi$")
  ax.set_ylabel("$r$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("Fold-Hopf")
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

def run_bifurcations():
  fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 10))
  fig.suptitle("Diagrammes de bifurcation")

  fold_bifurcation(fold, fold_df, ax1)
  hopf_bifurcation(hopf_polar, hopf_polar_df, ax2)
  fold_hopf_bifurcations(fold, fold_df, hopf_polar_coupled, hopf_polar_coupled_df, ax3)

  plt.savefig("../article/figures/bifurcations.pdf", dpi=300)

  plt.tight_layout()
  plt.show()

