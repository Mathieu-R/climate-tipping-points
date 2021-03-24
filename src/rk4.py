def rk4_derivatives(edo, v, t, dt, *args):
  k1 = edo(t, v, *args)
  k2 = edo(t + (dt / 2), v + ((dt / 2) * k1), *args)
  k3 = edo(t + (dt / 2), v + ((dt / 2) * k2), *args)
  k4 = edo(t + dt, v + (dt * k3), *args)
  return (k1 + 2*k2 + 2*k3 + k4)

def rk4(edo, v, t, dt, *args):
  return (dt / 6) * rk4_derivatives(edo, v, t, dt, *args)
