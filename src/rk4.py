class rk4():
  def derivatives(self, edo, time_step, tn, v, *args):
    k1 = edo(tn, v, args)
    k2 = edo(tn + (time_step / 2), v + ((time_step / 2) * k1), args)
    k3 = edo(tn + (time_step / 2), v + ((time_step / 2) * k2), args)
    k4 = edo(tn + time_step, v + (time_step * k3), args)
    return (k1 + 2*k2 + 2*k3 + k4)

  @staticmethod
  def solve(self, edo, time_step, t, v, *args):
    return (time_step / 6) * self.derivatives(edo, time_step, t, v, args)
