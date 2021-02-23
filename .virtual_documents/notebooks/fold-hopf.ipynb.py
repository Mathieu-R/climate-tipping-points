import numpy as np
from sympy import Matrix, symbols, init_printing, solve, lambdify
import matplotlib.pyplot as plt

init_printing()


a1, a2, b1, b2, c1, c2, x, y, z, gamma1, gamma2, phi = symbols("a1 a2 b1 b2 c1 c2 x y z phi")


def gamma(x):
    return gamma1 * x + gamma2

def f1(x):
    return a1 * (x ** 3) + a2 * x + phi

def f2(x, y, z):
    return b1 * z + b2 * (gamma(x) - (y**2 + z**2)) * y

def f3(x, y, z):
    return c1 * y + c2 * (gamma(x) - (y**2 + z**2)) * z


M = Matrix([f1(x), f2(x, y, z), f3(x, y, z)])


fixed_points = solve(M, [x, y, z])[2]
fixed_points


M.jacobian([x, y, z])


Jac = M.jacobian([x,y,z]).subs({x: fixed_points[0], y: fixed_points[1], z: fixed_points[2]})
Jac


Jac.charpoly()


Jac.eigenvals()


result = lambdify([a1, a2, b1, b2, c1, c2], J)



