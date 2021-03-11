import numpy as np

"""
forward euler method with a stochastical term
"""
def forward_euler_maruyama(edo, v, dt, *args):
  eta = np.random.normal(loc=0.0, scale=np.sqrt(dt))
  return edo(v, *args) * dt + eta * np.sqrt(dt)

