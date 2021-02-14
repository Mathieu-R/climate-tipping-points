import numpy as np

from constants import (a1, a2, b1, b2, c1, c2)

def fold(x, phi, gamma):
  return np.array([
    a1*(x**3) + a2*x + phi
  ])
  
def hopf(y, z, gamma):
  return np.array([
    b1*z + b2*(gamma(x) - (y**2 + z**2))*y,
    c1*y + c2*(gamma(x) - (y**2 + z**2))*z
  ])