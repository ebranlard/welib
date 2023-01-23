import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.system.wtmodels.model5CS import *
from welib.FEM.utils import deleteRowsCols
from welib.tools.eva import eigMCK, eigMK
from welib.tools.strings import printMat, printVec
from welib.system.mbc import Coleman2Comp


kt  = 11369784 # [Nm/rad]
f0 = 0.8      # [Hz] Blade edgewise natural frequency
om0 = 0.8*2*np.pi # [rad/s] Blade edgewise natural frequency
kx  = 200000   # [N/m] Lateral spring constant for M
ky  = 20000    # Vertical spring constant for M
l   = 30       # [m] Distance between m and M
m   = 500      # [kg] Mass of each rotating mass
M   = 50000    # [kg] Nacelle/generalized tower mass

kb = m*l**2*om0**2
plane='XYneg'; ordering='decreasing'



Omega = 0
# Omega = 0.8 
psi1  = 0
MM,CC,KK = systemMatrices(M, m, l, kb, kx, ky, Omega, psi1, plane=plane, ordering=ordering)
# MNR, CNR, KNR = systemMatricesNR(M, m, l, kb, kx, ky, Omega, plane=plane, ordering=ordering, method='numerical')
MNR, CNR, KNR = systemMatricesNR(M, m, l, kb, kx, ky, Omega, plane=plane, ordering=ordering, method='analytical')
MNR2 = deleteRowsCols(MNR, [4])
CNR2 = deleteRowsCols(CNR, [4])
KNR2 = deleteRowsCols(KNR, [4])
# printMat(MNR2, 'MNR2')
# printMat(CNR2, 'CNR2')
# printMat(KNR2, 'KNR2')

# Should return the same
# fd, zeta, Q, f0 = eigMCK(MNR2, CNR2, KNR2, method='full_matrix') #, normQ='byMax')
fd, zeta, Q, f0 = eigMCK(MNR2, CNR2, KNR2, method='state_space') #, normQ='byMax')
# fd, zeta, Q, f0 = eigMCK(MNR2, CNR2, KNR2, method='state_space_gen') #, normQ='byMax')
# fd, zeta, Q, f0 = eigMK(MNR2, KNR2)
printVec(f0,'f0', digits=3)
printVec(zeta,'zeta')
printMat(Q,'Q', digits=3)






print('------------- OMEGA =0.8rad/')
Omega = 0.8 
MNR, CNR, KNR = systemMatricesNR(M, m, l, kb, kx, ky, Omega, plane=plane, ordering=ordering, method='analytical')
MNR2 = deleteRowsCols(MNR, [4])
CNR2 = deleteRowsCols(CNR, [4])
KNR2 = deleteRowsCols(KNR, [4])
fd, zeta, Q, f0 = eigMCK(MNR2, CNR2, KNR2, method='state_space')
printVec(f0 ,'f0   ', digits=3)
printVec(zeta,'zeta', digits=3)
# print(Q[:,0])
print(Q[:,1])
Qprint = np.zeros(Q.shape)
for iMode in range(Q.shape[1]):
    a0,a1,b1 = Q[:3,iMode]
    x        = Q[3,iMode]
    A0,ABW,AFW,phi0,phiBW,phiFW = Coleman2Comp(a0, a1, b1)
    Qprint[0,iMode] = abs(x)
    Qprint[1:4,iMode] = [A0,ABW,AFW]

printMat(Qprint,'QSBW', digits=3)




# ---- Coleman
vOmega = np.linspace(0,np.pi,21)

f_NR = np.empty( (len(vOmega),  MNR2.shape[0]) )

for i, Omega in enumerate(vOmega):
    MNR, CNR, KNR = systemMatricesNR(M, m, l, kb, kx, ky, Omega, plane=plane, ordering=ordering, method='analytical')
    MNR2 = deleteRowsCols(MNR, [4])
    CNR2 = deleteRowsCols(CNR, [4])
    KNR2 = deleteRowsCols(KNR, [4])
    fd, zeta, Q, f0 = eigMCK(MNR2, CNR2, KNR2, method='state_space')
#     printVec(f0 ,'f0   ', digits=3)
#     printVec(zeta,'zeta', digits=3)
    # print(Q[:,0])
    f_NR[i,:] = f0

    Qprint = np.zeros(Q.shape)
    for iMode in range(Q.shape[1]):
        a0,a1,b1 = Q[:3,iMode]
        x        = Q[3,iMode]
        A0,ABW,AFW,phi0,phiBW,phiFW = Coleman2Comp(a0, a1, b1)
        Qprint[0,iMode] = abs(x)
        Qprint[1:4,iMode] = [A0,ABW,AFW]

# make the campbell diagram
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(7,3.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(vOmega*30/np.pi, f_NR[:,0], 'o', label='Tower side-side') 
ax.plot(vOmega*30/np.pi, f_NR[:,1], 'o', label='Backward whirling') 
ax.plot(vOmega*30/np.pi, f_NR[:,2], 'o', label='Symmetric') 
ax.plot(vOmega*30/np.pi, f_NR[:,3], 'o', label='Forward whirling') 
ax.plot(vOmega*30/np.pi,   vOmega/2/np.pi, ':', c='0.2', label='1P')
ax.plot(vOmega*30/np.pi, 3*vOmega/2/np.pi, '--', c='0.2', label='3P')
# ax.plot(vOmega*30/np.pi, om0/2/np.pi*np.ones_like(vOmega), '-.', c='0.5', label='Blade $f_n$')
ax.legend( bbox_to_anchor=(0., 1.04, 1., .102), loc='lower left',
           ncol=3, mode='expand', borderaxespad=0., fontsize=10)
ax.grid('on')
ax.set_xlabel(f'Rotational frequency [RPM]')
ax.set_ylabel('Frequencies [Hz]')
fig.tight_layout()



if __name__ == '__main__':

    from scipy.linalg import eig
    from welib.system.wtmodels.examples._ronnie_utils import get_nr_matrices, print_modal_components, animate_mode_shape

    # keyword arguments
    kwargs = dict(dofs=[1, 1, 1, 1, 0],  # flags to enable degrees of freedom [Bld1, Bld2, Bld3, x, y]
                  M=50_000,  # equivalent rotor/nacelle/tower mass [N/m]
                  kx=200_000,  # equivalent tower stiffness [N/m]
                  m=500,  # equivalent blade mass [kg]
                  l=30,  # equivalent blade length [m]
                  omega0=2*np.pi*0.8)  # blade natural frequency [rad/s]

    # get the system matrices in the non-rotating frame
    M_nr, C_nr, K_nr = get_nr_matrices(Omega, **kwargs)

    # assemble the state-space matrices
    n = M_nr.shape[0]  # number of dofs
    I, O = np.eye(n), np.zeros((n, n))
    A_nr = np.block([[O, I],
                     [-K_nr, -C_nr]])
    B_nr = np.block([[I, O],
                     [O, M_nr]])

    # solve the generalized eigenvalue problem
    [L, P] = eig(A_nr, B_nr)
    # print('P\n',P)

    # isolate the states that correspond to displacements
    L = L[::2]
    P = P[:n, ::2]
    # print('L\n',L)
    # print('P\n',P)

    # sort in increasing order
    argsort = np.argsort(np.abs(L))
    L = L[argsort]
    Q = P[:, argsort]
    # print('L\n',L)

    # calculate natural frequencies and damping
    omegas = np.abs(L)  # rad/s
    f0 = omegas / 2 / np.pi  # Hz
    zeta = -np.cos(np.angle(L)) * 100  # % critical
    printVec(f0, 'f0', digits=3)
    printVec(zeta, 'zeta', digits=3)
    printMat(Q, 'Q', digits=3)

    # extract the modal components in the non-rotating frame
    a0, a1, b1 = P[:3, :]  # extract the blade components from eigenvectors
    x = P[3, :]  # center mass displacement
    # print('a0',a0)

    plt.show()
