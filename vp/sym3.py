from sympy import *

d, a0, a2, v0, v2 = symbols("d a0 a2 v0 v2")

exp = 0.25*d**2*(a0**2+a2**2)+0.5*(d*((v0-v2)*(a0+a2))+v0**2+v2**2)-v0*v2

print(expand(exp))
