import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .model15DOFs import *

from welib.system.statespace import StateMatrix
from welib.system.eva import eig, eigA
from pyFAST.linearization.campbell import campbellModeStyles



np.set_printoptions(linewidth=300, precision =5)

vOmega=np.arange(0,1.55,0.05) #omegas up to rate omega Rotor Speed
# vOmega=np.arange(0,10.55,0.05) #omegas up to rate omega Rotor Speed


p=parameters() # Turbine parameters
lambda1 = 1.8751 # NOTE: constant for analytical bending mode 1
c1      = 0.7341 # NOTE: constant for analytical bending mode 1, see bendingMode!

df=pd.DataFrame()
df['Omega_[rad/s]']  = vOmega
df['f_flap_[Hz]']    = vOmega*0
df['f_edge_[Hz]']     = vOmega*0
df['f_tors_[Hz]'] = vOmega*0

for i,Omega in enumerate(vOmega):
    Mb,Kb,Gb,Mt,Gt,Kt = diagonalStructuralMatrices(p,Omega,lambda1,c1)
    A=StateMatrix(M=Mb,K=Kb)
    Q,freq =eig(K=Kb,M=Mb, freq_out=True, sort=True)
    freq2, zeta, Q2, freq_0  = eigA(A, nq=Mb.shape[0])
    df.loc[i,'f_flap_[Hz]'] = freq[0]
    df.loc[i,'f_edge_[Hz]'] = freq[1]
    df.loc[i,'f_tors_[Hz]'] = freq[2]

print(df)

df.to_csv('_outputs/A_freq.csv', sep='\t', index=False)


# --- Plot 
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.75, bottom=0.11, hspace=0.20, wspace=0.20)
lbl='1st blade flap' # TODO Mode 1
c, ls, ms, mk = campbellModeStyles(0, lbl)
ax.plot(vOmega, df['f_flap_[Hz]'], ls, marker=mk, color=c, label=lbl)

lbl='1st blade edge' # TODO Mode 2
c, ls, ms, mk = campbellModeStyles(1, lbl)
ax.plot(vOmega, df['f_edge_[Hz]'], ls, marker=mk, color=c, label=lbl)

lbl='1st blade torsion' # TODO TODO
c, ls, ms, mk = campbellModeStyles(2, lbl)
ax.plot(vOmega, df['f_tors_[Hz]'], ls, marker=mk, color=c, label=lbl)

ax.set_xlabel(r'$\Omega$ [rad/s]')
ax.set_ylabel(r'Frequency [Hz]')
ax.legend()


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
x=np.linspace(0,50,51)
nDOF= 3
phi_e = bendingMode(x,1)  
phi_t = torsionalMode(x,1)
phi_f = bendingMode(x,1)  
Phi = np.column_stack((phi_f, phi_e, phi_t))
# for j in np.arange(nDOF):
#     Phi[:,j] /= Phi[-1,j]

# print(Phi)

# --- Eigen vector 
def getShape(Omega):
    Mb,Kb,Gb,Mt,Gt,Kt = diagonalStructuralMatrices(p,Omega,lambda1,c1)
    A = StateMatrix(M=Mb,K=Kb)
    freq2, zeta, Q, freq_0  = eigA(A, nq=Mb.shape[0], normQ='byMax')

    Q0,freq0 =eig(K=Kb,M=Mb, freq_out=True, sort=True, normQ='byMax')
    phase = np.angle(Q)
    ampl  = np.abs(Q)
    #print(phase)
    #print(ampl)
    Modes = np.zeros((Q.shape[1], Phi.shape[0], Phi.shape[1]))
    for j in np.arange(Q.shape[1]):
        q_j = Q[:,j]
        #a_j = np.abs(q_j) # TODO figure out how to incorpoorate sign with phase
        phi_j = np.angle(q_j) 
        a_j = Q0[:,j] 
        Mode = np.multiply(Phi, a_j)
        Mode_end=Mode[-1,:]
        scale = np.max(Mode_end)
        Modes[j][:,:] = Mode/scale
    return Modes, phase

Modes0,phase0 = getShape(Omega=0)   
ModesR,phaseR = getShape(Omega=1.5)

for j in np.arange(nDOF):
    Mode0=Modes0[j]
    ModeR=ModesR[j]
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(x, Mode0[:,0], label='flap')
    ax.plot(x, Mode0[:,1], label='edge')
    ax.plot(x, Mode0[:,2], label='torsion')
    ax.plot(x, ModeR[:,0], 'kx')
    ax.plot(x, ModeR[:,1], 'kx')
    ax.plot(x, ModeR[:,2], 'kx')
    ax.set_xlabel('Blade span [m]')
    ax.set_ylabel('Mode content')
    ax.set_title('Mode {}'.format(j+1))
    ax.legend()
plt.show()




# plt.plot(df)
# plt.show()



# print(Omega)
# print('Mb\n',Mb)
# print('Kb\n',Kb)
# print('Mt\n',Mt)
# print('Gt\n',Gt)
# print('Kt\n',Kt)
# 
# print('Q\n',Q)
# print('freq\n',freq)
# 
# print('Q\n',np.imag(Q2))
# print('freq_d\n',freq_d)



# 
# fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
# fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
# ax.plot(x, phi_e    , label='')
# ax.set_xlabel('')
# ax.set_ylabel('')
# ax.legend()
# plt.show()
