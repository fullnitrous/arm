from sympy import *

x0, x1, x2, x3, y0, y1, y2, y3 = symbols("x0 x1 x2 x3 y0 y1 y2 y3")
a, b, c, d = symbols("a b c d")

eq0 = Eq(a*x0**3 + b*x0**2 + c*x0 + d, y0)
eq1 = Eq(a*x1**3 + b*x1**2 + c*x1 + d, y1)
eq2 = Eq(a*x2**3 + b*x2**2 + c*x2 + d, y2)
eq3 = Eq(a*x3**3 + b*x3**2 + c*x3 + d, y3)

sols = solve((eq0, eq1, eq2, eq3), (a, b, c, d))

print(sols)
