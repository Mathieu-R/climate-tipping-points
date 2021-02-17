import numpy as np

from constantss import (a1, a2, b1, b2, c1, c2)

def fold(t, v, phi):
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi
  ])

def hopf(t, v, phi):
  return np.array([
    phi * v[1] - (v[1] ** 3),
  ])

# v is a vector \vec{v}: [x, y, z]
def fold_hopf(self, t, v, gamma, phi):
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi,
    b1*v[2] + b2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[1],
    c1*v[1] + c2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[2]
  ])

# # v is a vector \vec{v}: [x, r, \theta]
# def fold_hopf_polar(t, v, phi, gamma):
#   return np.array([
#     a1 * (v[0] ** 3) + a2 * v[0] + phi,
#     phi * v[1] - (v[1] ** 3),
#     -1
#   ])
