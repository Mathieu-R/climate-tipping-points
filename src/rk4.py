def rk4_derivatives_vec(edo, v, t, dt, *args):
  k1 = dt * edo(t, v, *args)
  k2 = dt * edo(t + (dt / 2), v + ((dt / 2) * k1), *args)
  k3 = dt * edo(t + (dt / 2), v + ((dt / 2) * k2), *args)
  k4 = dt * edo(t + dt, v + (dt * k3), *args)
  return (k1 + 2*k2 + 2*k3 + k4)

def rk4_vec(edo, v, t, dt, *args):
  return rk4_derivatives_vec(edo, v, t, dt, *args) / 6

def rk4(edo, t, x, y, z, dt, *args):
  k1 = dt * edo(t, x, y, z, *args)
  k2 = dt * edo(t + (dt / 2), x + ((dt / 2) * k1), y + ((dt / 2) * k1), z + ((dt / 2) * k1), *args)
  k3 = dt * edo(t + (dt / 2), x + ((dt / 2) * k2), y + ((dt / 2) * k2), z + ((dt / 2) * k2), *args)
  k4 = dt * edo(t + dt, x + (dt * k3), y + (dt * k3), z + (dt * k3), *args)
  return (k1 + 2*k2 + 2*k3 + k4) / 6
