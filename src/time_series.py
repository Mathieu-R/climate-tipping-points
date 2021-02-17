import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from .edo import fold_hopf
from constantss import (x0, y0, z0, t_init, t_fin, time_step)

class time_series():
  def __init__(self):
    # initial conditions
    self.initial_conditions = np.array([[x0, y0, z0]])
    # forcing parameter
    self.phi = -2
    # time
    self.t_init = t_init
    self.t_fin = t_fin
    self.time_step = time_step
    # data
    self.dataset = self.initial_conditions
    self.time_range = np.arange(t_init, t_fin + time_step, time_step)
    self.legends = ["$x$ (leading)", "$y$ (following)", "$z$ (following)"]

  # linear coupling parameter
  # proposed by Dekker et al. article
  def gamma(self, x):
    return (-0.1 + 0.12*x)

  def derivatives(self, tn, v):
    k1 = fold_hopf(self, tn, v)
    k2 = fold_hopf(self, tn + (self.time_step / 2), v + ((self.time_step / 2) * k1))
    k3 = fold_hopf(self, tn + (self.time_step / 2), v + ((self.time_step / 2) * k2))
    k4 = fold_hopf(self, tn + self.time_step, v + (self.time_step * k3))
    return (k1 + 2*k2 + 2*k3 + k4)

  def rk4(self):
    # [x0, y0, z0]
    v = self.initial_conditions[0]

    for t in self.time_range:
      if t == self.t_fin: break # LOL

      # increase forcing parameter
      self.phi += 0.01

      v += (self.time_step / 6) * self.derivatives(t, v)
      self.dataset = np.append(self.dataset, [v.copy()], axis=0)

    # results = solve_ivp(
    #   fun=fold_hopf,
    #   t_span=(self.t_init, self.t_fin),
    #   t_eval=self.time_range,
    #   y0=self.initial_conditions,
    #   first_step=self.time_step,
    #   method="RK45",
    #   args=(self.phi, self.gamma,)
    # )

    self.plot()

  def plot(self):
    fig, ax = plt.subplots()

    l1, = ax.plot(self.time_range, self.dataset[:,0])
    l2, = ax.plot(self.time_range, self.dataset[:,1])
    l3, = ax.plot(self.time_range, self.dataset[:,2])

    plt.xlabel("t")
    plt.ylabel("variables")
    plt.legend(self.legends, loc="upper right")
    plt.title("Fold-Hopf Bifurcation - Time serie")
    plt.show()
