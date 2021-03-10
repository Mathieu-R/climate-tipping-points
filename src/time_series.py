import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from .edo import (fold_hopf, gamma)
from consts import (x0, y0, z0, t_init, t_fin, time_step)

from .rk4 import rk4
from .euler import forward_euler_maruyama

class time_series():
  def __init__(self):
    # initial conditions
    self.initial_conditions = np.array([[x0, y0, z0]])
    # forcing parameter
    self.phi = -2
    # time
    self.t_init = t_init
    self.t_fin = t_fin
    self.dt = time_step
    # data
    self.dataset = self.initial_conditions
    self.time_range = np.arange(t_init, t_fin + time_step, time_step)
    self.legends = ["$x$ (leading)", "$y$ (following)", "$z$ (following)"]

  def solve(self, solver):
    # [x0, y0, z0]
    v = self.initial_conditions[0]

    dataset = self.initial_conditions

    for t in self.time_range:
      if t == self.t_fin: break # LOL

      # increase forcing parameter
      self.phi += 0.001

      v += solver(fold_hopf, v, self.dt, gamma, self.phi)
      dataset = np.append(self.dataset, [v.copy()], axis=0)

    return dataset

  def basic(self):
    results = self.solve(rk4)
    return results

  def stochastic(self):
    stochastic_results = []
    # compute a lot of simulations
    for i in range(0,100):
      results = self.solve(forward_euler_maruyama)
      stochastic_results.append(results)
    return stochastic_results

  def plot(self, dataset):
    basic_results = self.basic()
    stochastic_results = self.stochastic()

    # 3 lines, 1 column
    fig, ax = plt.subplots(3, 1, figsize=(15, 7))
    fig.suptitle("Série temporelle")

    ax[0].plot(self.time_range, basic_results[:,0])
    ax[0].plot(self.time_range, basic_results[:,1])
    ax[0].plot(self.time_range, basic_results[:,2])

    ax[0].set_xlabel("$t$")
    ax[0].set_ylabel("$x$, $y$, $z$")
    ax[0].set_xlim(0,500)
    ax[0].set_ylim(-10,10)
    ax[0].set_legend(self.legends, loc="upper right")
    ax[0].set_title("Basique")

    ax[1].plot(self.time_range, stochastic_results[:,0])
    ax[1].plot(self.time_range, stochastic_results[:,1])
    ax[1].plot(self.time_range, stochastic_results[:,2])

    ax[1].set_xlabel("$t$")
    ax[1].set_ylabel("$x$, $y$, $z$")
    ax[1].set_xlim(0,500)
    ax[1].set_ylim(-10,10)
    ax[1].set_legend(self.legends, loc="center left", bbox_to_anchor=(1,0.5))
    ax[1].set_title("Stochastique")

    plt.tight_layout()
    plt.show()


