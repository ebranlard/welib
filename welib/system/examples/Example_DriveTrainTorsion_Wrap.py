"""
Integrate drive train torsion equation
Use as much wrapping functions as possible from welib/mech_system
  TODO: use yams simulator as well
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.system.mech_system import MechSystem
from welib.system.eva import eig, eigA
import welib.weio as weio
from welib.yams.windturbine import FASTWindTurbine

def main():
    # --- Parameters
    InFile = '../../../data/NREL5MW/Main_Onshore_DriveTrainTorsion.fst'

    # --- Turbine parameters
    OutFile = InFile.replace('.fst','.outb')
    WT = FASTWindTurbine(InFile)
    nGear  = WT.ED['GBRatio']
    K_DT   = WT.ED['DTTorSpr']
    D_DT   = WT.ED['DTTorDmp']
    Jr_LSS = WT.rot.inertia[0,0]*1.000    # bld + hub, LSS
    Jg_LSS = WT.gen.inertia[0,0]          # gen, LSS
    print('Jr',Jr_LSS)
    print('N ',nGear)


    # --- Read output from openfast
    df=weio.read(OutFile).toDataFrame()
    time  = df['Time_[s]'].values
    Q_r   = df['RtAeroMxh_[N-m]'].values
    Q_g   = df['GenTq_[kN-m]'].values*1000
    sDOF = ['Q_GeAz_[rad]', 'Q_DrTr_[rad]']
    sVel = ['QD_GeAz_[rad/s]', 'QD_DrTr_[rad/s]' ]
    sAcc = ['QD2_GeAz_[rad/s^2]','QD2_DrTr_[rad/s^2]']
    sF   = ['F_GeAz_[Nm]','F_DrTr_[Nm]']
    df[sF[0]] = Q_r - Q_g*nGear
    df[sF[1]] = Q_g*nGear

    # --- Setup Model
    # Everything on LSS: theta, and nu = theta_r-theta_g_LSS
    # [Jr + Jg N^2, -Jg N^2 ][theta_r]_dd + [0, 0][theta_r] + [0, 0][theta_r]_d =  Q_r - Q_g N
    # [  -Jg N^2  ,   Jg N^2][ nu    ]_dd + [0, K][ nu    ] + [0, B][ nu    ]_d =  Q_g N
    F = np.column_stack((Q_r-Q_g*nGear, Q_g*nGear)).T
    M = np.array([[Jr_LSS+Jg_LSS,-Jg_LSS],[-Jg_LSS,Jg_LSS]])
    K = K_DT*np.array([[0, 0],[0, 1]])
    D = D_DT*np.array([[0, 0],[0, 1]])

    # --- Define a system and perform time integration
    q0     = df.loc[0,sDOF].values
    qd0    = df.loc[0,sVel].values

    sys=MechSystem(M, D, K, F=(time,F), x0=q0     , xdot0=qd0 )
    print(sys)

    #res=sys.integrate(time, method='LSODA') # **options):
    res=sys.integrate(time, method='RK45') # **options):

    dfNL=sys.toDataFrame(sDOF+sVel, acc=True, sAcc=sAcc, forcing=True, sForcing=sF)
    print(dfNL.columns)
    dfNL['Q_GeAz_[rad]'] = dfNL['Q_GeAz_[rad]'].mod(2*np.pi)

    # sys.plot()

    # --- Plot states
    fig,axes = plt.subplots( int(sys.nStates/2),4, sharex=True, sharey=False, figsize=(12.8,8.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08, hspace=0.30, wspace=0.30)

    for i,chan in enumerate(sDOF):
        axes[i,0].plot(dfNL['Time_[s]'], dfNL[chan], label='Python')
        axes[i,0].plot(df  ['Time_[s]'], df  [chan], 'k:', label='OpenFAST')
        axes[i,0].set_xlabel('Time [s]')
        axes[i,0].set_ylabel(chan)
    for i,chan in enumerate(sVel):
        axes[i,1].plot(dfNL['Time_[s]'], dfNL[chan], label='Python')
        axes[i,1].plot(df  ['Time_[s]'], df  [chan], 'k:', label='OpenFAST')
        axes[i,1].set_xlabel('Time [s]')
        axes[i,1].set_ylabel(chan)
    for i,chan in enumerate(sAcc):
        axes[i,2].plot(dfNL['Time_[s]'], dfNL[chan], label='Python')
        axes[i,2].plot(df  ['Time_[s]'], df  [chan], 'k:', label='OpenFAST')
        axes[i,2].set_xlabel('Time [s]')
        axes[i,2].set_ylabel(chan)
    for i,chan in enumerate(sF):
        axes[i,3].plot(dfNL['Time_[s]'], dfNL[chan], label='Python')
        axes[i,3].plot(df  ['Time_[s]'], df  [chan], 'k:', label='OpenFAST')
        axes[i,3].set_xlabel('Time [s]')
        axes[i,3].set_ylabel(chan)

#     axes[0].set_ylabel(r'$\theta_r$ [rad]')
    axes[0,0].legend()

    # --- Frequencies
    #Q,Lambda = eig(K, M)# , freq_out=False, sort=True, normQ=None, discardIm=False):
    #freq_02 = np.sqrt(np.diag(Lambda))/(2*np.pi) # frequencies [Hz]
    #freq_d, zeta, Q, freq_0 = eigA(sys.A) #, nq=None, nq1=None, fullEV=False, normQ=None)
    #print('EVA  (zeta)    :',zeta   )
    #print('f_EVA          :',freq_d , '(damped)' )
    #print('f_EVA          :',freq_02, '(natural)' )
    #print('f_EVA          :',freq_0 , '(natural)' )
    #print('f_free-free    :',np.sqrt(K_DT/Jr_LSS + K_DT/Jg_LSS)/(2*np.pi))
    #print('f_fixed-free   :',np.sqrt(K_DT/(Jr_LSS))/(2*np.pi))



if __name__=="__main__":
    main()
    plt.show()
if __name__=="__test__":
    pass
