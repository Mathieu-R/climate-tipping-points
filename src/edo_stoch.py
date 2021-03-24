import numpy as np
from consts import (a1_stoch, a2_stoch, b1_stoch, b2_stoch, c1_stoch, c2_stoch, gamma1_stoch, gamma2_stoch)

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
  return gamma1_stoch + (gamma2_stoch * x)

def dx_stoch(t, x, y, z):
  return a1_stoch * (x ** 3) + a2_stoch * x + phi_time(t)

def dy_stoch(t, x, y, z):
  return b1_stoch*z + b2_stoch*(gamma(x) - (y**2 + z**2))*y

def dz_stoch(t, x, y, z):
  return c1_stoch*y + c2_stoch*(gamma(x) - (y**2 + z**2))*z

def fold_hopf_stoch(t, v, phi):
  #print(v, phi)
  return np.array([
    a1_stoch * (v[0] ** 3) + a2_stoch * v[0] + phi_time(t),
    b1_stoch * v[2] + b2_stoch * (gamma(v[0]) - (v[1]**2 + v[2]**2))*v[1],
    c1_stoch * v[1] + c2_stoch * (gamma(v[0]) - (v[1]**2 + v[2]**2))*v[2]
  ])
