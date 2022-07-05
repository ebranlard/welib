from sympy import Matrix, Symbol
import numpy as np


def stiffness6DOF(DOFs, frame, label='K', bDOFs=None, IKeep=None):
    """
    Return stiffness loads (force and moment)
     - DOFS: list of degrees of freedom, e.g. DOFs=[x, y, z, phi_x, phi_y, phi_z]
     - frame: frame to use for force and moment components
               e=[frame.x, frame.y, frame.z]
     - label: base for the symbols of the stiffnes matrix
     - bDOFs: which DOFs are active (true/false array of size 6)
     - IKeep: list of tuple pairs indicating the (i,j) index of the stiffness matrix to keep
              None implies all indices
    """
    # Index to keep in stiffness matrix
    if IKeep is None:
        IKeep = [(i,j) for i in range(6) for j in range(6)]

    # Creating symbols for stiffness coefficients
    KM=Matrix(np.zeros((6,6)))
    for i in range(6):
        for j in range(6):
            if j>=i:
                if (i,j) in IKeep:
                    s= Symbol('{}_{}{}'.format(label,i,j))
                    KM[i,j] = s
                    KM[j,i] = s
    # 
    fr=0; Mr=0
    e   =[frame.x, frame.y, frame.z]
    for i in range(3):
        for j in range(6):
            if bDOFs[j]:
                fr+= - KM[i,j  ] * DOFs[j] * e[i] 
                Mr+= - KM[i+3,j] * DOFs[j] * e[i]
    #print('IKeep:',IKeep)
    #print('bDOFs:',bDOFs)
    #print('KM   :\n',KM)
    #print('fr   :',fr)
    #print('Mr   :',Mr)
    #K_Mx, K_My, K_Mz          = symbols('K_x_M, K_y_M, K_z_M') # Mooring restoring
    #K_Mphix, K_Mphiy, K_Mphiz = symbols('K_phi_x_M, K_phi_y_M, K_phi_z_M') # Mooring restoring
    #C_Mx, C_My, C_Mz          = symbols('C_x_M, C_y_M, C_z_M') # Mooring restoring
    #C_Mphix, C_Mphiy, C_Mphiz = symbols('C_phi_x_M, C_phi_y_M, C_phi_z_M') # Mooring restoring
    #fr += -K_Mx * x *ref.frame.x #if bDOFs[0] else 0
    #fr += -K_My * y *ref.frame.y #if bDOFs[1] else 0
    #fr += -K_Mz * z *ref.frame.z #if bDOFs[2] else 0
    #Mr += -K_Mphix * phi_x *ref.frame.x #if bDOFs[3] else 0
    #Mr += -K_Mphiy * phi_y *ref.frame.y #if bDOFs[4] else 0
    #Mr += -K_Mphiz * phi_z *ref.frame.z #if bDOFs[5] else 0
    #fr += -C_Mx * x.diff(dynamicsymbols._t) *ref.frame.x #if bDOFs[0] else 0
    #fr += -C_My * y.diff(dynamicsymbols._t) *ref.frame.y #if bDOFs[1] else 0
    #fr += -C_Mz * z.diff(dynamicsymbols._t) *ref.frame.z #if bDOFs[2] else 0
    #Mr += -C_Mphix * phi_x.diff(dynamicsymbols._t) *ref.frame.x #if bDOFs[3] else 0
    #Mr += -C_Mphiy * phi_y.diff(dynamicsymbols._t) *ref.frame.y #if bDOFs[4] else 0
    #Mr += -C_Mphiz * phi_z.diff(dynamicsymbols._t) *ref.frame.z #if bDOFs[5] else 0
    return fr, Mr, KM
