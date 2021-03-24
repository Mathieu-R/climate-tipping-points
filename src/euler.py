import numpy as np
from consts import (mean, variance)

"""
forward euler method with a stochastical term
x_{i+1} = x_i * dt + zeta + sqrt(dt)
(forward-euler + gaussian noise * sqrt(dt))
"""
def forward_euler_maruyama(edo, t, x, y, z, dt, *args):
  zeta = np.random.normal(loc=mean, scale=np.sqrt(variance))
  return edo(t, x, y, z, *args) * dt + zeta * np.sqrt(dt)

