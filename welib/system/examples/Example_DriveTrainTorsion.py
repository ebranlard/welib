"""
Integrate drive train torsion equation
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
    model=2 # 1=1DOF (theta_r), 2= 2DOF (theta_r, nu), 22: (theta_r, theta_g_LSS), 3: (theta_r, theta_g)

    # --- Turbine parameters
    OutFile = InFile.replace('.fst','.outb')
    WT = FASTWindTurbine(InFile)
    nGear  = WT.ED['GBRatio']
    K_DT   = WT.ED['DTTorSpr']
    D_DT   = WT.ED['DTTorDmp']
    Jr_LSS = WT.rot.inertia[0,0]*1.000    # bld + hub, LSS
    Jg_LSS = WT.gen.inertia[0,0]          # gen, LSS
    Jg_HSS = WT.gen.inertia[0,0]/nGear**2
    print('Jr',Jr_LSS)
    print('Jg',Jg_HSS)
    print('N ',nGear)


    # --- Read output from openfast
    df=weio.read(OutFile).toDataFrame()
    time  = df['Time_[s]'].values
    Q_r   = df['RtAeroMxh_[N-m]'].values
    Q_g   = df['GenTq_[kN-m]'].values*1000
    Q_LSS = df['RotTorq_[kN-m]'].values*1000
    sDOF=['Q_GeAz_[rad]', 'Q_DrTr_[rad]']
    sVel=['QD_GeAz_[rad/s]', 'QD_DrTr_[rad/s]' ]
    sAcc=['QD2_GeAz_[rad/s^2]','QD2_DrTr_[rad/s^2]']
    q_r    = df[sDOF[0]] # NOTE: weird initial condition
    q_DT   = df[sDOF[1]]
    qd_r   = df[sVel[0]]
    qd_DT  = df[sVel[1]]
    qdd_r  = df[sAcc[0]]
    qdd_DT = df[sAcc[1]]
    GenSpeed     = df['GenSpeed_[rpm]']*2*np.pi/60
    GenSpeed_LSS = GenSpeed/nGear
    RotSpeed     = df['RotSpeed_[rpm]']*2*np.pi/60
    NuSpeed      = qd_DT
    # Generator from rotor + drivetrain torsion (note: jumps due to modulo)
    q_g      = (q_r  + q_DT)*nGear
    qd_g     = (qd_r + qd_DT)*nGear
    q_g_LSS  = (q_r  + q_DT)
    qd_g_LSS = (qd_r  + qd_DT)
    # Drive train torsion torque
    Q_DT = K_DT * q_DT + D_DT * qd_DT

    # --- Drivetrain speed
    # theta_gen = theta_rot - Nu
    # fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    # fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # ax.plot(GenSpeed_LSS    , label='Gen')
    # ax.plot(RotSpeed        , label='Rot')
    # # ax.plot(NuSpeed              , label='Nu')
    # ax.plot(RotSpeed-NuSpeed ,'--'    , label='Gen=Nu+Rot')
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    # ax.legend()


    # --- Setup Model
    if model==1:
        # One degree of freedom on LSS (torsion omitted): theta
        F = Q_r-Q_g*nGear
        F= F.reshape((1,-1))
        M = np.array([[Jr_LSS+Jg_LSS]])
        K = K_DT*np.array([[0]])
        D = D_DT*np.array([[0]])
    elif model==2:
        # Everything on LSS: theta, and nu = theta_r-theta_g_LSS
        # [Jr + Jg N^2, -Jg N^2 ][theta_r]_dd + [0, 0][theta_r] + [0, 0][theta_r]_d =  Q_r - Q_g N
        # [  -Jg N^2  ,   Jg N^2][ nu    ]_dd + [0, K][ nu    ] + [0, B][ nu    ]_d =  Q_g N
        F = np.column_stack((Q_r-Q_g*nGear, Q_g*nGear)).T
        M = np.array([[Jr_LSS+Jg_LSS,-Jg_LSS],[-Jg_LSS,Jg_LSS]])
        K = K_DT*np.array([[0, 0],[0, 1]])
        D = D_DT*np.array([[0, 0],[0, 1]])
    elif model==22:
        # Everything on LSS: theta_r, and theta_g_LSS
        # [Jr , 0    ][theta_r    ]_dd + [ K,-K][theta_r    ] + [ B,-B][theta_r    ]_d =  Q_r 
        # [0  ,Jg N^2][theta_g_LSS]_dd + [-K, K][theta_g_LSS] + [-B, B][theta_g_LSS]_d =  Q_g N
        F = np.column_stack((Q_r, -Q_g*nGear)).T
        M = np.array([[Jr_LSS, 0],[0,Jg_LSS]])
        K = K_DT*np.array([[1,-1],[-1, 1]])
        D = D_DT*np.array([[1,-1],[-1, 1]])

    elif model==3:
        # theta_r on LSS, theta_g on HSS
        # [Jr , 0 ][theta_r]_dd + [ K  ,-K/N  ][theta_r] + [ B,  -B/N  ][theta_r]_d =  Q_r 
        # [0,  Jg ][theta_g]_dd + [-K/N, K/N^2][theta_g] + [-B/N, B/N^2][theta_g]_d =  Q_g 
        F = np.column_stack((Q_r, -Q_g)).T
        M = np.array([[Jr_LSS, 0],[0,Jg_HSS]])
        K = K_DT*np.array([[1,-1/nGear],[-1/nGear, 1/nGear**2]])
        D = D_DT*np.array([[1,-1/nGear],[-1/nGear, 1/nGear**2]])

    # --- Define a system and perform time integration
    q0     = np.zeros(2)
    qd0    = np.zeros(2)
    q0[0]  = q_r[0]
    qd0[0] = qd_r[0]
    if model==22:
        q0[1]  = q_g_LSS[0]
        qd0[1] = qd_g_LSS[0]
    if model==3:
        q0[1]  = q_g[0]
        qd0[1] = qd_g[0]

    if model==1:
        sys=MechSystem(M, D, K, F=(time,F), x0=q0[0:1], xdot0=qd0[0:1] )
    else:
        sys=MechSystem(M, D, K, F=(time,F), x0=q0     , xdot0=qd0 )
    print(sys)

    #res=sys.integrate(time, method='LSODA') # **options):
    res=sys.integrate(time, method='RK45') # **options):

    #dfLI = sysLI.toDataFrame(self.channels, self.FASTDOFScales, x0=qop, xd0=qdop)
    dfNL=sys.toDataFrame(sDOF+sVel, acc=True, sAcc=sAcc)
    print(dfNL.columns)

    # sys.plot()

    # --- Plot states
    fig,axes = plt.subplots( sys.nStates,1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    axes[0].plot(sys.res.t, np.mod(sys.res.y[0,:],2*np.pi), label='Python')
    axes[0].plot(sys.res.t, q_r, label='OpenFAST')
    axes[0].set_ylabel(r'$\theta_r$ [rad]')
    axes[0].legend()

    if model==1:
        axes[1].plot(sys.res.t, sys.res.y[1,:])
        axes[1].plot(sys.res.t, qd_r)
        axes[1].set_ylabel(r'$\dot{\theta}_r$ [rad/s]')

    else:
        axes[2].plot(sys.res.t, sys.res.y[2,:])
        axes[2].plot(sys.res.t, RotSpeed,'--') # qd_r is not good here...
        axes[2].set_ylabel(r'$\dot{\theta}_r$ [rad/s]')


    if model==2:

        axes[1].plot(sys.res.t, sys.res.y[1,:])
        axes[1].plot(sys.res.t, q_DT)
        axes[1].set_ylabel(r'$\nu$ [rad]')

        axes[3].plot(sys.res.t, sys.res.y[3,:])
        axes[3].plot(sys.res.t, qd_DT)
        axes[3].set_ylabel(r'$\dot{\nu}$ [rad/s]')

    elif model==22:

        axes[1].plot(sys.res.t, np.mod(sys.res.y[1,:], 2*np.pi))
        axes[1].plot(sys.res.t, q_g_LSS)

        axes[3].plot(sys.res.t, sys.res.y[3,:])
        axes[3].plot(sys.res.t, GenSpeed_LSS)

    elif model==3:

        axes[1].plot(sys.res.t, np.mod(sys.res.y[1,:], 2*np.pi))
        axes[1].plot(sys.res.t, q_g)

        axes[3].plot(sys.res.t, sys.res.y[3,:])
        axes[3].plot(sys.res.t, GenSpeed)

    axes[-1].set_xlabel('Time [s]')



    # --- Frequencies
    Q,Lambda = eig(K, M)# , freq_out=False, sort=True, normQ=None, discardIm=False):
    freq_02 = np.sqrt(np.diag(Lambda))/(2*np.pi) # frequencies [Hz]
    freq_d, zeta, Q, freq_0 = eigA(sys.A) #, nq=None, nq1=None, fullEV=False, normQ=None)
    print('EVA  (zeta)    :',zeta   )
    print('f_EVA          :',freq_d , '(damped)' )
    print('f_EVA          :',freq_02, '(natural)' )
    print('f_EVA          :',freq_0 , '(natural)' )
    print('f_free-free    :',np.sqrt(K_DT/Jr_LSS + K_DT/Jg_LSS)/(2*np.pi))
    print('f_fixed-free   :',np.sqrt(K_DT/(Jr_LSS))/(2*np.pi))



if __name__=="__main__":
    main()
    plt.show()
if __name__=="__test__":
    pass
