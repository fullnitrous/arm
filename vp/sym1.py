from sympy import *

a0, v0, x0, a1, v1, x1, a2, v2, j, d = symbols("a0 v0 x0 a1 v1 x1 a2 v2 j d")

# - -> +
a1 = a0 + (-j)*x0
v1 = v0 + a0*x0 + 0.5*(-j)*x0**2

eq0 = Eq(d - (x0 + x1), 0)
eq1 = Eq(a2 - (a0 + (-j)*x0 + (j)*x1), 0)
eq2 = Eq(v2 - (v1 + a1*x1 + 0.5*j*x1**2), 0)

sols = solve((eq0, eq1, eq2), (x0, x1, j))

for sol in sols:
	print("---")
	print(sol)
	print("---")
