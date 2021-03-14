# 1D EDO for Fold bifurcation paramaters
a1 = -1
a2 = 1

# 2D EDO for Hopf bifurcation parameters
b1 = b2 = 1
c1 = -1
c2 = 1

# coupling parameters
gamma1 = -0.1
gamma2 = 0.12

# initial conditions
x0 = -0.5
y0 = 1.0
z0 = -1.0

r0 = 5.0

# time range
t_init = 0
t_fin = 500
time_step = 0.01

# gaussian noise parameters
mean = 0.0
variance = 0.1

# parameters for stochastic plot
a1_stoch = -1
a2_stoch = 1
b1_stoch = 0.1
b2_stoch = 1
c1_stoch = -0.5
c2_stoch = 1
gamma1_stoch = -0.2
gamma2_stoch = 0.3
time_step_stoch = 0.5

TRESHOLD = 10**-4
