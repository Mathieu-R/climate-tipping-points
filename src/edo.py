import numpy as np

# linear coupling parameter
# proposed by Dekker et al. article
def gamma(x):
    return (-0.1 + 0.12*x)

def fold(v, phi):
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi
  ])

def hopf(v, phi):
  return np.array([
    phi * v[0] - (v[0] ** 3),
  ])

def hopf_coupled(v, gamma, phi):
  rsquared = v[1] ** 2 + v[2] ** 2
  return np.array([
    gamma(v[0]) * rsquared - (rsquared ** 3)
  ])

# v is a vector \vec{v}: [x, y, z]
def fold_hopf_polar(v, gamma, phi):
  rsquared = v[1] ** 2 + v[2] ** 2
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi,
    gamma(v[0]) * rsquared - (rsquared ** 3)
  ])

# v is a vector \vec{v}: [x, y, z]
def fold_hopf(v, gamma, phi):
  return np.array([
    a1 * (v[0] ** 3) + a2 * v[0] + phi,
    b1*v[2] + b2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[1],
    c1*v[1] + c2*(gamma(v[0]) - (v[1]**2 + v[2]**2))*v[2]
  ])
