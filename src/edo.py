import numpy as np
from consts import (a1, a2, b1, b2, c1, c2, gamma1, gamma2)

# phi time dependent for time serie
def phi_time(t):
  # stabilize phi parameter for the first 100 a.u so the bifurcation does not start directly.
  # t = [0, 100] : phi = - 0.7
  # t = [101, end] : phi(t) = 0.05 * t ;
  # we stop at phi = 0.4
  return min(-0.7 + 0.05 * max(t - 100, 0), 0.7)

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

def dx(t, x, y, z):
  return a1 * (x ** 3) + a2 * x + phi_time(t)

def dy(t, x, y, z):
  return b1*z + b2*(gamma(x) - (y**2 + z**2))*y

def dz(t, x, y, z):
  return c1*y + c2*(gamma(x) - (y**2 + z**2))*z

# v is a vector \vec{v}: [x, y, z]
def fold_hopf(t, v, phi):
  print(a1, a2, b1, b2, c1, c2)
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi_time(t),
    b1*v[2] + b2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[1],
    c1*v[1] + c2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[2]
  ])
