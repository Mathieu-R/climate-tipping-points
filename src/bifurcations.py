import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as sciopt

from src.utils.utils import set_size

from consts import (x0, y0, z0, r0, a1, a2, b1, b2, c1, c2, t_init, t_fin, time_step, TRESHOLD)
from .edo import (fold, fold_df, hopf_polar, hopf_polar_df, hopf_polar_coupled, hopf_polar_coupled_df)
from .rk4 import rk4

plt.style.use("science")
#np.set_printoptions(precision=3, suppress=True)

"""Leading vs forcing \phi

"""
def fold_bifurcation(fx, df, ax):
  stable_equ_up = []
  stable_equ_down = []
  unstable_equ = []

  nphi = 1000
  nx = 10000

  x_min = -1.5
  x_max = 1.5

  phi_mesh = np.linspace(start=-1, stop=1, num=nphi)
  x_mesh = np.linspace(start=x_min, stop=x_max, num=nx)

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

  # annotations
  phi_down = stable_equ_up[0][0]
  phi_up = stable_equ_down[0][-1]

  ax.plot(phi_down, stable_equ_up[1][0], 'o', color="red", markersize=5)
  ax.plot([phi_down, phi_down], [x_min, stable_equ_up[1][0]], '--', color='Grey', linewidth=0.5, zorder=-1)

  ax.plot(phi_up, stable_equ_down[1][-1], 'o', color="red", markersize=5)
  ax.plot([phi_up, phi_up], [x_min, stable_equ_down[1][-1]], '--', color='Grey', linewidth=0.5, zorder=-1)

  ax.set_xlabel("$\phi$")
  ax.set_ylabel("$x$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("système primaire vs $\phi$")
  ax.legend(loc="best")

  return phi_down, phi_up

def fold_bifurcation_roots(fx, df, ax):
  stable_equ_up = []
  stable_equ_down = []
  unstable_equ = []

  nphi = 1000
  ntries = 5

  x_min = -1.5
  x_max = 1.5

  phi_mesh = np.linspace(start=-1, stop=1, num=nphi)
  x_mesh = np.linspace(start=x_min, stop=x_max, num=ntries)

  # algorithm
  for phi in phi_mesh:
    for x in x_mesh:
      sol = sciopt.root(fun=fx, jac=df, x0=x, method="hybr", tol=0.01, args=(phi))
      root = sol.x[0]
      # Check its stability (stable or unstable)
      df_value = df(root, phi)
      if df_value > 0:
        unstable_equ.append([phi, root])
      elif df_value < 0:
        if x > 0.:
          stable_equ_up.append([phi, root])
        elif x < 0.:
          stable_equ_down.append([phi, root])

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
          if r > 0. + 10**-2:
            stable_equ_up.append([gamma, r])
          elif r < 0. - 10**-2:
            stable_equ_down.append([gamma, r])
          else:
            #print([gamma, r])
            stable_equ_middle.append([gamma, r])

  unstable_equ = list(zip(*unstable_equ))
  stable_equ_middle = list(zip(*stable_equ_middle))
  stable_equ_up = list(zip(*stable_equ_up))
  stable_equ_down = list(zip(*stable_equ_down))

  # plot
  ax.plot(stable_equ_middle[0], stable_equ_middle[1], color="black", label="ligne stable")
  ax.plot(stable_equ_up[0], stable_equ_up[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(stable_equ_down[0], stable_equ_down[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(unstable_equ[0], unstable_equ[1], linestyle="dashdot", color="red", label="ligne instable")
  ax.set_xlabel("$\gamma$")
  ax.set_ylabel("$r$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("système secondaire vs $\gamma$")
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

"""Following vs coupling \gamma

"""
def hopf_bifurcation_roots(fx, df, ax):
  stable_equ_middle = []
  stable_equ_up = []
  stable_equ_down = []
  unstable_equ = []

  ngamma = 1000
  nr = 10

  r_min = -1.5
  r_max = 1.5

  ntries = 5

  gamma_mesh = np.linspace(start=-1, stop=1, num=ngamma)
  r_mesh = np.linspace(start=r_min, stop=r_max, num=ntries)

  # algorithm
  # start from r_min as guess
  guess = r_min
  for gamma in gamma_mesh:
    # look for some roots.
    for r in r_mesh:
      sol = sciopt.root(fun=fx, jac=df, x0=r, method="hybr", tol=0.01, args=(gamma))
      root = sol.x[0]
      #print(root)
      #guess = root
      #print(root)
      # Check its stability (stable or unstable)
      df_value = df(root, gamma)
      if df_value > 0:
        unstable_equ.append([gamma, root])
      elif df_value < 0:
        # decorralate points (otherwise matplotlib connect every points of each branch with a line)
        if root > 0.:
          stable_equ_up.append([gamma, root])
        elif root < 0.:
          stable_equ_down.append([gamma, root])
        else:
          stable_equ_middle.append([gamma, root])

  unstable_equ = list(zip(*unstable_equ))
  stable_equ_middle = list(zip(*stable_equ_middle))
  stable_equ_up = list(zip(*stable_equ_up))
  stable_equ_down = list(zip(*stable_equ_down))

  # plot
  ax.plot(stable_equ_middle[0], stable_equ_middle[1], color="black", label="ligne stable")
  ax.plot(stable_equ_up[0], stable_equ_up[1], linestyle="dashdot", linewidth=1, color="green", label="ligne stable - osc.")
  ax.plot(stable_equ_down[0], stable_equ_down[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(unstable_equ[0], unstable_equ[1], linestyle="dashed", color="red", label="ligne instable")

  # annotations
  ax.plot(unstable_equ[0][0], unstable_equ[1][0], 'o', color="orange", markersize=5)
  ax.plot([unstable_equ[0][0], unstable_equ[0][0]], [r_min, unstable_equ[1][0]], '--', color='Grey', linewidth=0.5, zorder=-1)

  ax.set_xlabel("$\gamma$")
  ax.set_ylabel("$r$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("système secondaire vs $\gamma$")
  ax.legend(loc="best")


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
              if r > 0. + 5*10**-2:
                stable_equ_up.append([phi, r])
              elif r < 0 - 5*10**-2:
                stable_equ_down.append([phi, r])
              else:
                stable_equ_middle.append([phi, r])

  unstable_equ = list(zip(*unstable_equ))
  stable_equ_middle = list(zip(*stable_equ_middle))
  stable_equ_up = list(zip(*stable_equ_up))
  stable_equ_down = list(zip(*stable_equ_down))

  # plot
  ax.plot(stable_equ_middle[0], stable_equ_middle[1], color="black", label="ligne stable")
  ax.plot(stable_equ_up[0], stable_equ_up[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(stable_equ_down[0], stable_equ_down[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(unstable_equ[0], unstable_equ[1], linestyle="dashdot", color="red", label="ligne instable")
  ax.set_xlabel("$\phi$")
  ax.set_ylabel("$r$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("Système secondaire vs $\phi$")
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

def fold_hopf_bifurcations_roots(px, dp, sr, ds, ax, phi_down, phi_up):
  stable_equ_middle = []
  stable_equ_up = []
  stable_equ_down = []
  unstable_equ = []

  nphi = 1000
  ntries = 5

  x_min = r_min = -1.5
  x_max = r_max = 1.5

  phi_mesh = np.linspace(start=-1, stop=1, num=nphi)

  x_mesh = np.linspace(start=x_min, stop=x_max, num=10000)
  r_mesh = np.linspace(start=r_min, stop=r_max, num=ntries)
  # loop each \phi
  for phi in phi_mesh:
    for x in x_mesh:
      p_value = px(x, phi)
      if - TRESHOLD < p_value < TRESHOLD:
        # loop for each r
        for r in r_mesh:
          sol_r = sciopt.root(fun=sr, jac=ds, x0=r, method="hybr", tol=0.01, args=(x))
          root_r = sol_r.x[0]
          #print(root_r)
          # check the stability
          #dp_value = dp(root_x, phi)
          ds_value = ds(root_r, x)
          if ds_value > 0:
            #print([phi, root_r])
            unstable_equ.append([phi, root_r])
          elif ds_value < 0:
            if root_r > 0.:
              stable_equ_up.append([phi, root_r])
            elif root_r < 0.:
              stable_equ_down.append([phi, root_r])
            else:
              stable_equ_middle.append([phi, root_r])

  unstable_equ = list(zip(*unstable_equ))
  stable_equ_middle = list(zip(*stable_equ_middle))
  stable_equ_up = list(zip(*stable_equ_up))
  stable_equ_down = list(zip(*stable_equ_down))

  # plot
  ax.plot(stable_equ_middle[0], stable_equ_middle[1], color="black", label="ligne stable")
  ax.plot(stable_equ_up[0], stable_equ_up[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(stable_equ_down[0], stable_equ_down[1], linestyle="dashdot", color="green", label="ligne stable - osc.")
  ax.plot(unstable_equ[0], unstable_equ[1], linestyle="dashdot", color="red", label="ligne instable", zorder=-1)

  # annotations
  ax.plot(unstable_equ[0][0], unstable_equ[1][0], 'o', color="orange", markersize=5)
  ax.plot([unstable_equ[0][0], unstable_equ[0][0]], [r_min, unstable_equ[1][0]], '--', color='Grey', linewidth=0.5, zorder=-1)

  ax.plot(phi_down, stable_equ_middle[1][0], 'o', color="red", markersize=5)
  ax.plot([phi_down, phi_down], [r_min, stable_equ_middle[1][0]], '--', color='Grey', linewidth=0.5, zorder=-1)

  ax.plot(phi_up, stable_equ_middle[1][0], 'o', color="red", markersize=5)
  ax.plot([phi_up, phi_up], [r_min, stable_equ_middle[1][-1]], '--', color='Grey', linewidth=0.5, zorder=-1)

  ax.set_xlabel("$\phi$")
  ax.set_ylabel("$r$")
  ax.set_xlim(-1, 1)
  ax.set_ylim(-1.5, 1.5)
  ax.set_title("système secondaire vs $\phi$")
  ax.legend(loc="best")

def run_bifurcations():
  fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(set_size(width="full-size", subplots=(1,3))), sharey=True)
  fig.suptitle("Diagrammes de bifurcation")

  phi_down, phi_up = fold_bifurcation(fold, fold_df, ax1)
  hopf_bifurcation_roots(hopf_polar, hopf_polar_df, ax2)
  fold_hopf_bifurcations_roots(fold, fold_df, hopf_polar_coupled, hopf_polar_coupled_df, ax3, phi_down, phi_up)

  ax1.text(0.5, 1.1, "(a)", ha="center", transform=ax1.transAxes, size=8)
  ax2.text(0.5, 1.1, "(b)", ha="center", transform=ax2.transAxes, size=8)
  ax3.text(0.5, 1.1, "(c)", ha="center", transform=ax3.transAxes, size=8)

  #plt.savefig("article/figures/bifurcations.pdf", dpi=300)

  # tight_layout() not recommended for figure that go in article.
  #plt.tight_layout()
  plt.show()

