import numpy as np

"""
forward euler method with a stochastical term
"""
def forward_euler_maruyama(edo, v, dt):
  eta = np.random.normal(loc=0.0, scale=np.sqrt(dt))
  return edo(v) * dt + eta * np.sqrt(dt)

