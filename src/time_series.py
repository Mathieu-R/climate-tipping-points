import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from .edo import fold_hopf
from consts import (x0, y0, z0, t_init, t_fin, time_step)

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


  def rk4(self):
    # [x0, y0, z0]
    v = self.initial_conditions[0]

    for t in self.time_range:
      if t == self.t_fin: break # LOL

      # increase forcing parameter
      self.phi += 0.001

      v += rk4(fold_hopf, self.time_step, v, gamma, self.phi)
      self.dataset = np.append(self.dataset, [v.copy()], axis=0)

    self.plot()

  def plot(self):
    fig, ax = plt.subplots(figsize=(15, 7))

    l1, = ax.plot(self.time_range, self.dataset[:,0])
    l2, = ax.plot(self.time_range, self.dataset[:,1])
    l3, = ax.plot(self.time_range, self.dataset[:,2])

    plt.xlabel("t")
    plt.ylabel("variables")
    plt.xlim(0,500)
    plt.ylim(-10,10)
    plt.legend(self.legends, loc="upper right")
    plt.title("Fold-Hopf Bifurcation - Time serie")
    plt.show()

time_serie = time_series()
time_serie.rk4()
