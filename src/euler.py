import numpy as np
from consts import (mean, variance)

"""
forward euler method with a stochastical term
x_{i+1} = x_i * dt + zeta + sqrt(dt)
(forward-euler + gaussian noise * sqrt(dt))
"""
def forward_euler_maruyama(edo, v, dt, *args):
  zeta = np.random.normal(loc=mean, scale=np.sqrt(variance))
  print(zeta, edo(v, *args) * dt)
  return edo(v, *args) * dt + zeta * np.sqrt(dt)

