import numpy as np
from consts import (a1, a2, b1, b2, c1, c2)

# note: v is a vector. example: v = [x, y, z]

# linear coupling parameter
# proposed by Dekker et al. article
def gamma(x):
    return (-0.1 + 0.12*x)

"""
@param v: [x]
"""
def fold(v, phi):
  return a1 * (v[0] ** 3) + a2 * v[0] + phi

def fold_df(v, phi):
  return 3 * a1 * (v[0] ** 2) + a2

"""
@param v: [r]
"""
def hopf_polar(v, phi):
  return phi * v[0] - (v[0] ** 3),

def hopf_polar_df(v, phi):
  return phi - 3 * (v[0] ** 2)

"""
@param v: [x, r]
"""
def hopf_polar_coupled(v, phi):
  return gamma(v[0]) * v[1] - (v[1] ** 3)

def hopf_polar_coupled_df(v, phi):
  return gamma(v[0]) - 3 * (v[1] ** 2)

# def hopf_coupled(v, gamma, phi):
#   rsquared = v[1] ** 2 + v[2] ** 2
#   return np.array([
#     gamma(v[0]) * rsquared - (rsquared ** 3)
#   ])

# v is a vector \vec{v}: [x, y, z]
# def fold_hopf_polar(v, gamma, phi):
#   rsquared = v[1] ** 2 + v[2] ** 2
#   return np.array([
#     a1 * (v[0] ** 3) + a2 * v[0] + phi,
#     gamma(v[0]) * rsquared - (rsquared ** 3)
#   ])

# v is a vector \vec{v}: [x, y, z]
def fold_hopf(v, gamma, phi):
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi,
    b1*v[2] + b2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[1],
    c1*v[1] + c2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[2]
  ])
