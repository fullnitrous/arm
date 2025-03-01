from sympy import *

x0, x1, x2 = symbols("x0 x1 x2")
v0, v1, v2 = symbols("v0 v1 v2")
a,  b,  c  = symbols("a b c")

eq0 = Eq((a*x0**2 + b*x0 + c) - v0, 0)
eq1 = Eq((a*x1**2 + b*x1 + c) - v1, 0)
eq2 = Eq((a*x2**2 + b*x2 + c) - v2, 0)

sols = solve((eq0, eq1, eq2), (a, b, c))

print(simplify(sols[a]))
print(simplify(sols[b]))
print(simplify(sols[c]))
