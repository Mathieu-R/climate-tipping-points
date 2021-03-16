import numpy as np
from consts import (a1_stoch, a2_stoch, b1_stoch, b2_stoch, c1_stoch, c2_stoch, gamma1_stoch, gamma2_stoch)


# linear coupling parameter
# proposed by Dekker et al. article
def gamma(x):
  return gamma1_stoch + (gamma2_stoch * x)

def fold_hopf_stoch(v, phi):
  #print(v, phi)
  return np.array([
    a1_stoch * (v[0] ** 3) + a2_stoch * v[0] + phi,
    b1_stoch * v[2] + b2_stoch * (gamma(v[0]) - (v[1]**2 + v[2]**2))*v[1],
    c1_stoch * v[1] + c2_stoch * (gamma(v[0]) - (v[1]**2 + v[2]**2))*v[2]
  ])
