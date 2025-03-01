from sympy import *

x1, v1, a1 = symbols("x1 v1 a1")
a, b, c = symbols("a b c")
j_max = symbols("j_max")

eq0 = Eq(a*x1**2+b*x1+c-v1, 0)
eq1 = Eq(2*a*x1+b-a1, 0)
eq2 = Eq(2*a-j_max, 0)

sols = solve((eq0, eq1, eq2), (a, b, c))

print(sols)

#expa = (a0*x0 - a0*x1 + a1*x0 - a1*x1 - 2*v0 + 2*v1)/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)
#expb = (-a0*x0**2 - a0*x0*x1 + 2*a0*x1**2 - 2*a1*x0**2 + a1*x0*x1 + a1*x1**2 + 3*v0*x0 + 3*v0*x1 - 3*v1*x0 - 3*v1*x1)/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)
#expc = (2*a0*x0**2*x1 - a0*x0*x1**2 - a0*x1**3 + a1*x0**3 + a1*x0**2*x1 - 2*a1*x0*x1**2 - 6*v0*x0*x1 + 6*v1*x0*x1)/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)
#expd = (-a0*x0**2*x1**2 + a0*x0*x1**3 - a1*x0**3*x1 + a1*x0**2*x1**2 + 3*v0*x0*x1**2 - v0*x1**3 + v1*x0**3 - 3*v1*x0**2*x1)/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)

#print(simplify(expa))
#print(simplify(expb))
#print(simplify(expc))
#print(simplify(expd))
