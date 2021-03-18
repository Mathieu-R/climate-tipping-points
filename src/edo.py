import numpy as np
from consts import (a1, a2, b1, b2, c1, c2, gamma1, gamma2)

# linear coupling parameter
# proposed by Dekker et al. article
def gamma(x):
  return gamma1 + (gamma2 * x)

"""Leading vs forcing \phi
  :x: float
"""
def fold(x, phi):
  return (a1 * (x ** 3)) + (a2 * x) + phi

def fold_df(x, phi):
  return (3 * a1 * (x ** 2)) + a2

"""Following vs coupling \gamma
    \gamma act like a forcing parameter
  :r: float
"""
def hopf_polar(r, gam):
  return (gam * r) - (r ** 3)

def hopf_polar_df(r, gam):
  return gam - (3 * (r ** 2))

"""Following vs forcing \phi
    \gamma act like a coupled parameter through x and so \phi
  :r: float
"""
def hopf_polar_coupled(r, x):
  return (gamma(x) * r) - (r ** 3)

def hopf_polar_coupled_df(r, x):
  return gamma(x) - (3 * (r ** 2))

# v is a vector \vec{v}: [x, y, z]
def fold_hopf(v, phi):
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi,
    b1*v[2] + b2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[1],
    c1*v[1] + c2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[2]
  ])
