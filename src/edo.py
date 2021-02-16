import numpy as np

from constantss import (a1, a2, b1, b2, c1, c2)

# x is a vector \vec{x}: [x, y, z]
def fold_hopf(t, x, phi, gamma):
  return np.array([
    a1 * (x[0] ** 3) + a2 * x[0] + phi,
    b1*x[2] + b2*(gamma(x[0]) - (x[1]**2 + x[2]**2))*x[1],
    c1*x[1] + c2*(gamma(x[0]) - (x[1]**2 + x[2]**2))*x[2]
  ])
