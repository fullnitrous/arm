from sympy import *

x0, v0, a0, x1, v1, a1 = symbols("x0 v0 a0 x1 v1 a1")
a, b, c, d = symbols("a b c d")

eq0 = Eq(a*x0**3 + b*x0**2 + c*x0 + d - v0, 0)
eq1 = Eq(a*x1**3 + b*x1**2 + c*x1 + d - v1, 0)
eq2 = Eq(3*a*x0**2 + 2*b*x0 + c - a0, 0)
eq3 = Eq(3*a*x1**2 + 2*b*x1 + c - a1, 0)

sols = solve((eq0, eq1, eq2, eq3), (a, b, c, d))

print(sols)
