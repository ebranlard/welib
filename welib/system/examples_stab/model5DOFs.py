import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio



def systemMatrices(M, mb, l, kb, kx, ky, Omega, psi1, psi2, psi3, symb=False):
    MM = np.zeros((5,5))
    DD = np.zeros((5,5))
    KK = np.zeros((5,5))
    if symb:
        from sympy import Matrix, sin, cos
        MM = Matrix(MM)
        DD = Matrix(DD)
        KK = Matrix(KK)
    else:
        from np import sin, cos
    MM[0,0] = mb*l**2
    MM[1,1] = mb*l**2
    MM[2,2] = mb*l**2
    MM[3,3] = M + 3*mb
    MM[4,4] = M + 3*mb
    MM[0,3] = -mb*l*sin(psi1)
    MM[1,3] = -mb*l*sin(psi2)
    MM[2,3] = -mb*l*sin(psi3)
    MM[0,4] =  mb*l*cos(psi1)
    MM[1,4] =  mb*l*cos(psi2)
    MM[2,4] =  mb*l*cos(psi3)
    MM[3,0] =  MM[0,3]
    MM[3,1] =  MM[1,3]
    MM[3,2] =  MM[2,3]
    MM[4,0] =  MM[0,4]
    MM[4,1] =  MM[1,4]
    MM[4,2] =  MM[2,4]

    DD[3,0] = -2*mb *l *Omega * cos(psi1)
    DD[3,1] = -2*mb *l *Omega * cos(psi2)
    DD[3,2] = -2*mb *l *Omega * cos(psi3)
    DD[4,0] = -2*mb *l *Omega * sin(psi1)
    DD[4,1] = -2*mb *l *Omega * sin(psi2)
    DD[4,2] = -2*mb *l *Omega * sin(psi3)

    KK[0,0] =  l**2 * kb
    KK[1,1] =  l**2 * kb
    KK[2,2] =  l**2 * kb
    KK[3,3] =  kx
    KK[4,4] =  ky
    KK[3,0] =  mb *l *Omega**2 * sin(psi1)
    KK[3,1] =  mb *l *Omega**2 * sin(psi2)
    KK[3,2] =  mb *l *Omega**2 * sin(psi3)
    KK[4,0] = -mb *l *Omega**2 * cos(psi1)
    KK[4,1] = -mb *l *Omega**2 * cos(psi2)
    KK[4,2] = -mb *l *Omega**2 * cos(psi3)

    return MM, DD, KK


if __name__ == '__main__':
    pass
