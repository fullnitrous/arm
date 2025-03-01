from sympy import *

v2, a2, j_max, x0, x = symbols("v2 a2 j x0 x")

exp = v2 + a2*(x - x0) + 0.5*j_max*(x - x0)**2

print(collect(expand(exp), x0))
