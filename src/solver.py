import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from .edo import (fold_hopf)
from constantss import (x0, y0, z0, t_init, t_fin, time_step)

class solver():
  def __init__(self):
    # initial conditions
    self.initial_conditions = [x0, y0, z0]
    # forcing parameter
    self.phi = 0.38
    # time
    self.t_init = t_init
    self.t_fin = t_fin
    self.time_step = time_step
    # data
    self.dataset = []
    self.time_range = np.arange(t_init, t_fin + time_step, time_step)
    self.legends = ["x", "y", "z"]

  # linear coupling parameter
  # proposed by Dekker et al. article
  def gamma(self, x):
    return (-0.1 + 0.12*x)

  def solve(self):
    results = solve_ivp(
      fun=fold_hopf,
      t_span=(self.t_init, self.t_fin),
      t_eval=self.time_range,
      y0=self.initial_conditions,
      first_step=self.time_step,
      method="RK45",
      args=(self.phi, self.gamma,)
    )

    self.dataset = results.y
    print(len(self.time_range))
    print(len(self.dataset[0]))
    self.plot()

  def plot(self):
    plt.plot(self.time_range, self.dataset[0,:], self.dataset[1,:], self.dataset[2,:])
    plt.xlabel("t")
    plt.ylabel("variables")
    plt.legend(self.legends, loc="upper right")
    plt.title("Fold-Hopf Bifurcation - Time serie")
    plt.show()
