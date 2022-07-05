import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .model15DOFs import *

from welib.system.statespace import StateMatrix
from welib.system.eva import eig, eigA, eigMCK
from pyFAST.linearization.campbell import campbellModeStyles



np.set_printoptions(linewidth=300, precision =5)



# --- Polar
alpha = np.linspace(-np.pi,np.pi,100)
Cl,Cd,Clp,Cdp = liftDrag(alpha)
alpha=alpha*180/np.pi
M=np.column_stack((alpha,Cl,Cd,Clp,Cdp))
df = pd.DataFrame(data=M,columns=['alpha_[deg]','Cl_[-]','Cd_[-]','Clp_[-]','Cdp_[-]'])
df.to_csv('_outputs/B1_polars.csv',index=False,sep=',')



aero=True

p=parameters() # Turbine parameters
lambda1 = 1.8751 # NOTE: constant for analytical bending mode 1
c1      = 0.7341 # NOTE: constant for analytical bending mode 1, see bendingMode!


Omega=1.5
vWS=np.arange(5,26,1)
nDOF=3


if aero:
    flag=''
else:
    flag='_NoAero'



x=np.linspace(0,50,51)
phi_e = bendingMode(x,1)  
phi_t = torsionalMode(x,1)
phi_f = bendingMode(x,1)  
Phi = np.column_stack((phi_f, phi_e, phi_t))


FreqDamp=np.zeros((len(vWS),nDOF*2));
ModesAmplitudes = np.zeros((len(vWS),nDOF, nDOF));
ModesPhases     = np.zeros((len(vWS),nDOF, nDOF));
ue = np.zeros((len(x), nDOF))
uf = np.zeros((len(x), nDOF))
ut = np.zeros((len(x), nDOF))
phie = np.zeros((len(x), nDOF))
phif = np.zeros((len(x), nDOF))
phit = np.zeros((len(x), nDOF))

for iWS,WS in enumerate(vWS):
    print('--------------------------',WS)
    Mb,Kb,Gb,Mt,Gt,Kt = diagonalStructuralMatrices(p,Omega,lambda1,c1)

    aeroMat = aerodynMatrix(p,Omega,WS)
    Cab=aeroMat['Cab'];
    Kab=aeroMat['Kab'];
    if aero:
        Kb+= Kab
        Cb = Cab
    else:
        Cb  = Kb*0
    A=StateMatrix(M=Mb, K=Kb, C=Cb)
    #Q,freq =eig(K=Kb, M=Mb, C=Cb, freq_out=True, sort=True)
    freq1, zeta1, Q1, freq_01  = eigMCK(Mb, Cb, Kb)
    freq2, zeta2, Q , freq_02  = eigA(A, nq=Mb.shape[0])
    np.testing.assert_almost_equal(freq1,freq2)
    np.testing.assert_almost_equal(zeta1,zeta2)

    FreqDamp[iWS,0:nDOF*2:2]=freq2 
    FreqDamp[iWS,1:nDOF*2:2]=zeta2

    # For each mode, extract "content"

    M=x

    for iMod in np.arange(nDOF):
        ModeVector = Q[:,iMod]
        imax       = np.argmax(np.abs(ModeVector))
        zm         = ModeVector[imax]
        ModeVectorScaled=ModeVector*np.conj(zm)/abs(zm)**2 # complex, but normalized to unity for max amplitude 
        ModesPhases[iWS,iMod,:]=np.angle(ModeVectorScaled)/(2*np.pi)

        PhaseMod = np.angle(ModeVectorScaled);
        AmpMod   = np.abs(ModeVectorScaled)  ;
#        
        if (imax==1 or imax ==2):
            Amp=AmpMod/2  # because shape functions are not normalized...
        else:
            Amp=AmpMod
        uf[:,iMod]=phi_f*Amp[0]
        ue[:,iMod]=phi_e*Amp[1]
        ut[:,iMod]=phi_t*Amp[2]
        phif[:,iMod]=[PhaseMod[0]]*len(x)
        phie[:,iMod]=[PhaseMod[1]]*len(x)
        phit[:,iMod]=[PhaseMod[2]]*len(x)

        M=np.column_stack((M,uf[:,iMod], phif[:,iMod], ue[:,iMod], phie[:,iMod], ut[:,iMod], phit[:,iMod]))
    if  WS==8 or WS == 16:
        filename='_outputs/B1_modes_{:}ms{}.csv'.format(WS,flag)
        print(filename)
        headers=['r_[m]','m1_a_flap_[-]','m1_phi_flap_[rad]','m1_a_edge_[-]','m1_phi_edge_[rad]','m1_a_tors_[-]','m1_phi_tors_[rad]','m2_a_flap_[-]','m2_phi_flap_[rad]','m2_a_edge_[-]','m2_phi_edge_[rad]','m2_a_tors_[-]','m2_phi_tors_[rad]','m3_a_flap_[-]','m3_phi_flap_[rad]','m3_a_edge_[-]','m3_phi_edge_[rad]','m3_a_tors_[-]','m3_phi_tors_[rad]'];
        df      = pd.DataFrame(data=M, columns = headers)
        df.to_csv(filename, index=False, sep=',')
        if WS==8:
            Res1=M
        else:
            Res2=M

# --- Export Cambell Data
M       = np.column_stack((vWS,FreqDamp))
headers = ['WS_[m/s]','F1_[Hz]','D1_[-]','F2_[Hz]','D2_[-]','F3_[Hz]','D3_[-]'];
df      = pd.DataFrame(data=M, columns = headers)
df.to_csv('_outputs/B1_freq{}.csv'.format(flag),index=False,sep=',')

# --- Plot Blade Campbell
fig,axes = plt.subplots(1, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.75, bottom=0.11, hspace=0.20, wspace=0.20)
lbl='Mode 1 (1st blade flap)' # TODO Mode 1
c, ls, ms, mk = campbellModeStyles(0, lbl)
axes[0].plot(vWS, df['F1_[Hz]'], ls, marker=mk, color=c, label=lbl)
axes[1].plot(vWS, df['D1_[-]'] , ls, marker=mk, color=c, label=lbl)

lbl='Mode 2 (1st blade edge)' # TODO Mode 2
c, ls, ms, mk = campbellModeStyles(1, lbl)
axes[0].plot(vWS, df['F2_[Hz]'], ls, marker=mk, color=c, label=lbl)
axes[1].plot(vWS, df['D2_[-]'],  ls, marker=mk, color=c, label=lbl)

lbl='Mode 3 (1st blade torsion)' # TODO TODO
c, ls, ms, mk = campbellModeStyles(2, lbl)
axes[0].plot(vWS, df['F3_[Hz]'], ls, marker=mk, color=c, label=lbl)
axes[1].plot(vWS, df['F3_[Hz]'], ls, marker=mk, color=c, label=lbl)

axes[0].set_xlabel(r'WS [m/s]')
axes[0].set_ylabel(r'Frequency [Hz]')
axes[1].set_ylabel(r'Damping ratio [-]')
axes[0].legend()



# --- Plot mode shapes at 8 and 16m/s
for iMod in np.arange(3):
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(x,Res1[:,iMod*2*nDOF+1:(iMod*2*nDOF+nDOF*2):2]  )
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('Mode {} - WS 8'.format(iMod)) 
    #ax.legend('Flap','Edge','Torsion',0)
# for iMod in np.arange(3):
#     figure(iMod+100)
#     plot(x,Res2(:,1+(iMod-1)*2*ndof+(1:2:(ndof*2)) ))
#     title(sprintf('Mode %d - WS 16',iMod)); 
#     legend('Flap','Edge','Torsion',0)
# end



plt.show()

