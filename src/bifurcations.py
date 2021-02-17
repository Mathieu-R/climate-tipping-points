import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib
colors = matplotlib.cm.get_cmap('plasma')

from tqdm import tqdm # progress bar
from constantss import (x0, y0, z0, r0, t_init, t_fin, time_step)
from .edo import (fold, hopf)
from .rk4 import rk4

class bifurcations():
  def __init__(self):
    self.phi = 0
    self.phi_range = np.linspace(-2, 2, 100)
    self.time_range = np.arange(t_init, t_fin, time_step)

  def fold(self):
    initial_conditions = [x0]
    dataset = pd.DataFrame()

    for phi in tqdm(self.phi_range):
      time_serie = []
      v = initial_conditions
      for t in self.time_range:
        time_serie.append(v)
        v += rk4(fold, time_step, t, v, phi)

      dataset[phi] = pd.Series(time_serie)

    self.plot(dataset)

  def hopf(self):
    initial_conditions = [r0]
    dataset = pd.DataFrame()

    for phi in tqdm(self.phi_range):
      time_serie = []
      v = initial_conditions
      for t in self.time_range:
        time_serie.append(v)
        v += rk4(fold, time_step, t, v, phi)

      dataset[phi] = pd.Series(time_serie)

    self.plot(dataset)

  def plot(self, dataset):
    print(dataset)
    dataset.T.plot(
      figsize = (16,6),
      ylim = (-5,5),
      legend = False,
      colormap=colors ,
      alpha = 0.3,
      title = "Fold Bifurcation Diagram"
    )
    plt.xlabel("$\phi$")
    plt.ylabel("x")
