from sympy import *

r, x, k = symbols('r x k')
a, b, c = symbols('a b c')
delta, zeta, mu = symbols('delta zeta mu')
b               = -(k**2+x**2+r**2)
a               = k**2
c               = x**2
delta           = b**2 - 4 * a*c
zeta           = sqrt( (- b + sqrt(delta)) / (2*a) )
mu = x/(k*zeta)

dmudx   = diff(mu,x)
dzetadx = diff(zeta,x)
dmudr   = diff(mu,r)
dzetadr = diff(zeta,r)


print('delta')
print(delta)
print('Mu')
print(mu)
print('zeta')
print(zeta)

print('dmu_dx=')
print(dmudx)
print('dzeta_dx=')
print(dzetadx)

print('dmu_dr=')
print(dmudr)
print('dzeta_dr=')
print(dzetadr)

print('test')
print(diff(sqrt(x),x))
