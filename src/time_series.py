import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

from .edo import (fold_hopf)
from .edo_stoch import (fold_hopf_stoch)
from .rk4 import rk4
from .euler import forward_euler_maruyama

from consts import (x0, y0, z0, t_init, t_fin, time_step, time_step_stoch)

plt.style.use("seaborn-whitegrid")
np.set_printoptions(precision=3, suppress=True)

class time_series():
  def __init__(self):
    # initial conditions
    self.initial_conditions = np.array([[x0, y0, z0]])
    # forcing parameter
    self.phi = -2
    # time
    self.t0 = t_init
    self.tN = t_fin
    # stochastic number of iterations
    self.niter = 10
    # legend
    self.legends = ["$x$ (leading)", "$y$ (following)", "$z$ (following)"]

  def solve(self, solver, edo, dt, nt):
    # [x0, y0, z0]
    #v = self.initial_conditions[0]

    # vector mesh -- will receive the result
    v_mesh = np.ones((nt, 3)) #
    # set inital conditions
    v_mesh[0] = self.initial_conditions[0]
    print(edo, dt)

    for t in tqdm(range(0, nt - 1)):
      v_mesh[t + 1] = v_mesh[t] + solver(edo, v_mesh[t], dt, self.phi)

      # increase forcing parameter
      self.phi += 0.001

    return v_mesh

  def basic(self):
    dt = time_step
    nt = int((self.tN - self.t0) / dt)
    time_mesh_basic = np.linspace(start=self.t0, stop=self.tN, num=nt)

    # [ [x0, y0, z0], [x1, y1, z1], ..., [xN, yN, zN] ]
    results = self.solve(rk4, fold_hopf, dt, nt)
    return time_mesh_basic, results

  def stochastic(self):
    dt = time_step_stoch
    nt = int((self.tN - self.t0) / dt)
    time_mesh_stoch = np.linspace(start=self.t0, stop=self.tN, num=nt)


    stochastic_results = np.ones((self.niter, nt, 3))
    # compute a lot of simulations
    for i in tqdm(range(0, self.niter)):
      results = self.solve(forward_euler_maruyama, fold_hopf_stoch, dt, nt)
      stochastic_results[i] = results
      return
    return time_mesh_stoch, stochastic_results

  def plot(self):
    #time_mesh_basic, basic_results = self.basic()
    time_mesh_stoch, stochastic_results = self.stochastic()
    return

    # 3 lines, 1 column
    fig, ((ax1), (ax2)) = plt.subplots(2, 1, figsize=(15, 7))
    fig.suptitle("Série temporelle")

    ax1.plot(time_mesh_basic, basic_results[:,0])
    ax1.plot(time_mesh_basic, basic_results[:,1])
    ax1.plot(time_mesh_basic, basic_results[:,2])

    ax1.set_xlabel("$t$")
    ax1.set_ylabel("$x$, $y$, $z$")
    ax1.set_xlim(0,500)
    ax1.set_ylim(-1.5,1.5)
    ax1.legend(self.legends, loc="upper right")
    ax1.set_title("Basique")

    print(stochastic_results)

    ax2.stackplot(time_mesh_stoch, stochastic_results[:,:,0])
    ax2.stackplot(time_mesh_stoch, stochastic_results[:,:,1])
    ax2.stackplot(time_mesh_stoch, stochastic_results[:,:,2])

    ax2.set_xlabel("$t$")
    ax2.set_ylabel("$x$, $y$, $z$")
    ax2.set_xlim(0,500)
    ax2.set_ylim(-1.5,1.5)
    ax2.legend(self.legends, loc="center left", bbox_to_anchor=(1,0.5))
    ax2.set_title("Stochastique")

    plt.tight_layout()
    plt.show()


