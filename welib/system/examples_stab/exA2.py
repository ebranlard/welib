import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .model15DOFs import *

from welib.system.statespace import StateMatrix
from welib.system.eva import eig, eigA, eigMCK
from pyFAST.linearization.campbell import campbellModeStyles

#       % 1  , 2   , 3   , 4   , 5   , 6   , 7 , 8 , 9 , 10 , 11
LIM_TOP=[0.6 , 0.4 , 2.4 , 0.3 , 0.3 , 0.2 , 3 , 3 , 3,  3,  3]
LIM_BOT=[0.6 , 0.4 , 2.4 , 0.3 , 0.3 , 0.2 , 3 , 3 , 3,  3,  3]
LIM_TOP=[0.6 , 0.4 , 2.4 , 0.3 , 0.3 , 0.2 , .3 , .3 , .3,  .3,  .3]
LIM_BOT=[0.6 , 0.4 , 2.4 , 0.3 , 0.3 , 0.2 , .3 , .3 , .3,  .3,  .3]


np.set_printoptions(linewidth=300, precision =5)

# vOmega=np.arange(0,1.55,0.05) #omegas up to rate omega Rotor Speed
vOmega=np.arange(0,1.55,0.15) #omegas up to rate omega Rotor Speed


p=parameters() # Turbine parameters
lambda1 = 1.8751 # NOTE: constant for analytical bending mode 1
c1      = 0.7341 # NOTE: constant for analytical bending mode 1, see bendingMode!



psi1 = 0
nDOF = 15         

Freq  = np.zeros((len(vOmega), nDOF))
Modes = np.zeros((len(vOmega), nDOF, nDOF), dtype=np.complex128)


for i,Omega in enumerate(vOmega):
    M,C,K        = MCKmat( p, Omega, psi1 )
    B,mu,Bdot,R  = Bmat( psi1, Omega )
    MB,CB,KB     = MbCbKbmat( B,M,C,K,mu )
    MBT, CBT,KBT = MCKTildemat( B,MB,CB,KB,mu,R )

    A=StateMatrix(M=MBT, K=KBT, C=CBT)
    # Q,freq =eig(K=KBT, M=MBT, freq_out=True, sort=True)
    freq_d1, zeta1, Q1, freq0_1 = eigMCK(MBT, CBT, KBT, method='full_matrix')
    freq_d2, zeta2, Q2, freq0_2 = eigA(A, nq=MBT.shape[0])
#     print('{:5.2f} {:5.2f}   {:5.4f} {:5.4f}'.format(freq0_2[4],  freq0_1[4], zeta2[4], zeta1[4]))
#     print('{:5.2f} {:5.2f}   {:5.4f} {:5.4f}'.format(freq0_2[-1], freq0_1[-1], zeta2[-1], zeta1[-1]))
    #print(np.around(freq0_1-freq0_2,4))

    Freq[i,:]    = freq0_2
    Modes[i,:,:] = Q2

print('nDOFs',A.shape, MBT.shape)

cols=['Omega_[rad/s]']+['f{:d}'.format(i+1) for i in np.arange(11)]
print(cols)
M = np.column_stack((vOmega, Freq[:,:11]))
print(M.shape)
df=pd.DataFrame(columns=cols, data=M)
df.to_csv('_outputs/A2_freq.csv', index=False, sep='\t')

# --- Plot Campbell
ModeID= {
        0:'1st tower ss?',
        1:'1st tower fa?',
        2:'Drivetrain Torsion?',
        3:'1st blade Flap BW',
        4:'1st blade Flap coll',
        5:'1st blade Flap FW',
        6:'1st blade Edge BW',
        7:'1st blade Edge FW',
        8:'2nd blade Flap BW',
        9:'2nd blade Flap coll',
       10:'2nd blade Flap FW',
}

fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.75, bottom=0.11, hspace=0.20, wspace=0.20)
for j in np.arange(11):
    try:
        lbl=ModeID[j]
    except:
        lbl='Mode {}'.format(j+1)

    c, ls, ms, mk = campbellModeStyles(j, lbl)


    ax.plot(vOmega, Freq[:, j], ls, marker=mk, color=c, label=lbl)

# ax.legend()
ax.legend(bbox_to_anchor=(0., 1.02, 1.00, .802), loc='lower left', ncol=3, mode="expand", borderaxespad=0.)




# ---

# for iMode in np.arange(0,11):
for iMode in np.arange(0):
    Mode= Modes[:,:,iMode]
    Amplitude=abs(Mode)


    a0f=Mode[:,0]
    a0e=Mode[:,1]
    a0t=Mode[:,2]
    a1f=Mode[:,3]
    a1e=Mode[:,4]
    a1t=Mode[:,5]
    b1f=Mode[:,6]
    b1e=Mode[:,7]
    b1t=Mode[:,8]
#     Result(iOm, iAcc+  (2:7)  )= [A0 , phi0, ABW, phiBW, AFW, phiFW ]; 
#     # idem edgewise component 
#     [A0,ABW,AFW,phi0,phiBW,phiFW]=fAss(a0e,a1e,b1e);
#     Result(iOm, iAcc+  (8:13)  )= [A0 , phi0, ABW, phiBW, AFW, phiFW ];
#     % idem torsional component
#     [A0,ABW,AFW,phi0,phiBW,phiFW]=fAss(a0t,a1t,b1t);




    fig,axes = plt.subplots(2, 3, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)


    # Amplitude and phase of the flapwise rotor component
    A0,ABW,AFW,phi0,phiBW,phiFW = Coleman2Comp(a0f, a1f, b1f)
    ax=axes[0,0]
    ax.plot(vOmega, A0 *2, 'r')
    ax.plot(vOmega, ABW*2, 'g')
    ax.plot(vOmega, AFW*2, 'b')
    ax.set_xlim([0,1.6])
    ax.set_ylim([0,LIM_TOP[iMode]])

    # Amplitude and phase of the edgewise rotor component
    A0,ABW,AFW,phi0,phiBW,phiFW = Coleman2Comp(a0e, a1e, b1e)
    ax=axes[0,1]
    ax.plot(vOmega, A0 *2, 'r')
    ax.plot(vOmega, ABW*2, 'g')
    ax.plot(vOmega, AFW*2, 'b')
    ax.set_xlim([0,1.6])
    ax.set_ylim([0,LIM_TOP[iMode]])

    # Amplitude and phase of the torsional rotor component
    A0,ABW,AFW,phi0,phiBW,phiFW = Coleman2Comp(a0t, a1t, b1t)
    ax=axes[0,2]
    ax.plot(vOmega, A0, 'r')
    ax.plot(vOmega, ABW, 'g')
    ax.plot(vOmega, AFW, 'b')
    ax.set_xlim([0,1.6])
    ax.set_ylim([0,LIM_TOP[iMode]])


    ax=axes[1,0]
    ax.plot(vOmega, Amplitude[:,9], 'r')
    ax.plot(vOmega, Amplitude[:,10], 'g')
    ax.set_xlim([0,1.6])
    ax.set_ylim([0,LIM_BOT[iMode]])

    ax=axes[1,1]
    ax.plot(vOmega, Amplitude[:,11]*180/np.pi, 'r')
    ax.plot(vOmega, Amplitude[:,12]*180/np.pi, 'g')
    ax.set_xlim([0,1.6])
    ax.set_ylim([0,LIM_BOT[iMode]])

    ax=axes[1,2]
    ax.plot(vOmega, Amplitude[:,13]*180/np.pi, 'r')
    ax.plot(vOmega, Amplitude[:,14]*180/np.pi, 'g')
    ax.set_xlim([0,1.6])
    ax.set_ylim([0,LIM_BOT[iMode]])


#     ax.set_xlabel('')
#     ax.set_ylabel('')



    fig.suptitle('Mode {}'.format(iMode+1))

# ax.legend()
plt.show()
