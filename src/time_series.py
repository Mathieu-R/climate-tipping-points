import numpy as np
import matplotlib.pyplot as plt

# import pandas as pd
# import seaborn as sns

from tqdm import tqdm
#from scipy.stats import median_abs_deviation

from .edo import (fold_hopf, gamma, dx, dy, dz)
from .edo_stoch import (fold_hopf_stoch, dx_stoch, dy_stoch, dz_stoch)
from .rk4 import (rk4, rk4_vec)
from .euler import forward_euler_maruyama

from src.utils.utils import set_size

from consts import (x0, y0, z0, t_init, t_fin, time_step, time_step_stoch, phi_factor, phi_factor_stoch)

plt.style.use("science")
np.set_printoptions(precision=3, suppress=True)

class time_series():
  def __init__(self):
    # initial conditions
    self.initial_conditions = np.array([[x0, y0, z0]])
    # forcing parameter
    self.phi = -1
    # time
    self.t0 = t_init
    self.tN = t_fin
    # stochastic number of iterations
    self.niter = 100
    # legend
    self.legends = ["$x$ (leading)", "$y$ (following)", "$z$ (following)"]

  def phi_time(self, t):
    # stabilize phi parameter for the first 100 a.u so the bifurcation does not start directly.
    # t = [0, 100] : phi = - 0.7
    # t = [101, end] : phi(t) = 0.05 * t ;
    # we stop at phi = 0.4
    return min(-0.7 + 0.05 * max(t - 100, 0), 0.7)

  """ Plot a temporal serie

  :solver: function -- (edo, v, dt, *args) -> ndarray
  :edo: function -- (v, phi) -> ndarray
  :dt: float -- time step
  :nt: int -- number of time steps (N)

  :return: ndarray -- [[x0, y0, z0], ..., [x{N-1}, y{N-1}, z{N-1}]]
  """
  def solve(self, solver, edox, edoy, edoz, dt, nt, phi_factor):
    # [x0, y0, z0]
    #v = self.initial_conditions[0]

    # vector mesh -- will receive the result
    #v_mesh = np.ones((nt, 3)) #
    x_mesh = np.zeros((nt))
    y_mesh = np.zeros((nt))
    z_mesh = np.zeros((nt))
    # set inital conditions
    #v_mesh[0] = self.initial_conditions[0]
    x_mesh[0] = x0
    y_mesh[0] = y0
    z_mesh[0] = z0
    #print(edo, dt)

    # phi mesh -- stock the forcing parameter function of time for plotting
    phi_mesh = np.ones((nt, 1))

    for k in range(0, nt - 1):
      t = k * dt
      #v_mesh[k + 1] = v_mesh[k] + solver(edo, v_mesh[k], t, dt, self.phi)
      x_mesh[k+1] = x_mesh[k] + solver(edox, t, x_mesh[k], y_mesh[k], z_mesh[k], dt)
      y_mesh[k+1] = y_mesh[k] + solver(edoy, t, x_mesh[k], y_mesh[k], z_mesh[k], dt)
      z_mesh[k+1] = z_mesh[k] + solver(edoz, t, x_mesh[k], y_mesh[k], z_mesh[k], dt)
      #print(self.phi_time(t))
      phi_mesh[k + 1] = self.phi_time(t)
      self.phi = self.phi_time(t) #phi_factor * t

    # reset forcing parameter
    # self.phi = -1
    return x_mesh, y_mesh, z_mesh, phi_mesh

  def basic(self):
    dt = time_step
    nt = int((self.tN - self.t0) / dt)
    time_mesh_basic = np.linspace(start=self.t0, stop=self.tN, num=nt)

    # [ [x0, y0, z0], [x1, y1, z1], ..., [xN, yN, zN] ]
    x, y, z, phi_mesh = self.solve(rk4, dx, dy, dz, dt, nt, phi_factor)
    return time_mesh_basic, x, y, z, phi_mesh

  # def stochastic_seaborn(self):
  #   dt = time_step_stoch
  #   nt = int((self.tN - self.t0) / dt)
  #   time_mesh_stoch = np.linspace(start=self.t0, stop=self.tN, num=nt)
  #   reshaped_time_mesh_stoch = np.reshape(time_mesh_stoch, (1, -1))

  #   results = np.ones((3, nt))
  #   stochastic_results = np.ones((self.niter, nt, 4))
  #   # compute a lot of simulations
  #   for i in tqdm(range(0, self.niter)):
  #     results = self.solve(forward_euler_maruyama, fold_hopf_stoch, dt, nt, phi_factor_stoch)
  #     # concatenate (x, y, z) columns with (t) column
  #     results = np.concatenate((reshaped_time_mesh_stoch.T, results), axis=1)
  #     stochastic_results[i] = results

  #   # compute the mean
  #   # mean = np.mean(stochastic_results, axis=0)

  #   # data transformation for data visualisation
  #   # stack all the (t, x, y, z) column vectors of all the simulations
  #   stacked_stochastic_results = np.reshape(stochastic_results, (-1, 4))

  #   # create panda dataframe from it
  #   df_stochastic = pd.DataFrame(data=stacked_stochastic_results, columns=["t", "x", "y", "z"])

  #   # aggregate x, y, z
  #   df_stochastic = df_stochastic.reset_index()[["t", "x", "y", "z"]].melt(
  #     id_vars="t",
  #     value_vars=("x", "y", "z"),
  #     var_name="position",
  #     value_name="variables"
  #   )

  #   return time_mesh_stoch, df_stochastic

  def stochastic(self):
    dt = time_step_stoch
    nt = int((self.tN - self.t0) / dt)
    time_mesh_stoch = np.linspace(start=self.t0, stop=self.tN, num=nt)

    stochastic_results_x = np.ones((self.niter, nt))
    stochastic_results_y = np.ones((self.niter, nt))
    stochastic_results_z = np.ones((self.niter, nt))
    # compute a lot of simulations
    for i in tqdm(range(0, self.niter)):
      x, y, z, phi_mesh = self.solve(forward_euler_maruyama, dx_stoch, dy_stoch, dz_stoch, dt, nt, phi_factor_stoch)
      stochastic_results_x[i] = x
      stochastic_results_y[i] = y
      stochastic_results_z[i] = z

    # compute the mean
    mean_x = np.mean(stochastic_results_x, axis=0)
    mean_y = np.mean(stochastic_results_y, axis=0)
    mean_z = np.mean(stochastic_results_z, axis=0)

    # compute standard deviation
    deviation_x = np.std(stochastic_results_x, axis=0)
    deviation_y = np.std(stochastic_results_y, axis=0)
    deviation_z = np.std(stochastic_results_z, axis=0)

    return time_mesh_stoch, mean_x, mean_y, mean_z, deviation_x, deviation_y, deviation_z


  def plot(self):
    time_mesh_basic, x, y, z, phi_mesh = self.basic()
    time_mesh_stoch, mean_x, mean_y, mean_z, deviation_x, deviation_y, deviation_z = self.stochastic()

    # 2 lines, 1 column
    fig, ((ax1), (ax2)) = plt.subplots(2, 1, figsize=(483.69687 * 1 / 72.27, 2), sharex=True)
    #fig.suptitle("Série temporelle")

    ax1.plot(time_mesh_basic, x, color="black")
    ax1.plot(time_mesh_basic, y, color="red")
    ax1.plot(time_mesh_basic, z, color="gold")
    # phi parameter function of time
    ax1.plot(time_mesh_basic, phi_mesh, color="green", alpha=0.7)
    # gamma parameter function of time
    #print(gamma(basic_results[:,0]))
    ax1.plot(time_mesh_basic, gamma(x), color="blue", alpha=0.7, zorder=-1)

    #ax1.set_xlabel("$t$")
    ax1.set_ylabel("$x$, $y$, $z$, $\phi$, $\gamma$")
    ax1.set_xlim(0,500)
    ax1.set_ylim(-3,3)
    #ax1.legend(self.legends, loc="center left", bbox_to_anchor=(1,0.5))
    #ax1.set_title("Basique")

    ax2.plot(time_mesh_stoch, mean_x, color="black")
    ax2.plot(time_mesh_stoch, mean_y, color="red")
    ax2.plot(time_mesh_stoch, mean_z, color="gold")

    ax2.fill_between(time_mesh_stoch, mean_x - deviation_x, mean_x + deviation_x, alpha=0.2, color="black")
    ax2.fill_between(time_mesh_stoch, mean_y - deviation_y, mean_y + deviation_y, alpha=0.2, color="red")
    ax2.fill_between(time_mesh_stoch, mean_z - deviation_z, mean_z + deviation_z, alpha=0.2, color="gold")

    ax2.set_xlabel("temps [u.a]")
    ax2.set_ylabel("$x$, $y$, $z$")
    ax2.set_xlim(0,500)
    ax2.set_ylim(-3,3)
    #ax2.legend(self.legends, loc="center left", bbox_to_anchor=(1,0.5))
    #ax2.set_title("Stochastique")

    ax1.text(0.035, 0.85, "(a)", ha="center", transform=ax1.transAxes, size=8, fontweight="bold")
    ax2.text(0.035, 0.85, "(b)", ha="center", transform=ax2.transAxes, size=8, fontweight="bold")

    plt.savefig("article/figures/time-series.pdf", dpi=300)

    #plt.tight_layout()
    plt.show()


