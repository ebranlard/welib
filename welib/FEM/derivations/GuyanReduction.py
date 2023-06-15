""" 
Show that Guyan-Reduction of a single element result in rigid body modes

"""

import numpy as np
import sympy
from sympy import Symbol
from sympy import Matrix
from sympy.abc import *
from fem.frame3d import frame3d_KeMe

display=lambda x: sympy.pprint(x, use_unicode=False,wrap_line=False)

Kv = Symbol('Kv')
Ix = Symbol('Ix')
Iy = Symbol('Iy')
Iz = Symbol('Iz')
Mass = Symbol('M')

Ke,Me = frame3d_KeMe(E,G,Kv,E*A,E*Ix,E*Iy,E*Iz,L,A,Mass)
Ke = Matrix(Ke)
Me = Matrix(Me)
Kmm = Ke[:6,:6]
Ksm = Ke[6:,:6]
Kss = Ke[6:,6:]
Kssm1 = Kss.inv()

print('-----------------------------------------------------------------------------------------')
print('--- Ke ')
print('-----------------------------------------------------------------------------------------')
display(Ke) 
#     display(Me) 
print('-----------------------------------------------------------------------------------------')
print('--- Kss ')
print('-----------------------------------------------------------------------------------------')
display(Kss)
print('-----------------------------------------------------------------------------------------')
print('--- Ksm ')
print('-----------------------------------------------------------------------------------------')
display(Ksm)
print('-----------------------------------------------------------------------------------------')
print('--- Kss^-1 ')
print('-----------------------------------------------------------------------------------------')
display(Kssm1)
print('-----------------------------------------------------------------------------------------')
print('--- Phi1 = -Kss^-1 Ksm ')
print('-----------------------------------------------------------------------------------------')
Phi1 = -Kssm1 * Ksm
display(Phi1)
