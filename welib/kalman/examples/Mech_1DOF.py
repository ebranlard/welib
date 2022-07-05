""" 
This scripts uses a Kalamn filter to estimate the position/speed and wave loading on a monopile based on the measurement of the acceleration.


"""
import numpy as np
import matplotlib.pyplot as plt
import os
from  welib.kalman.kalman import BuildSystem_Linear_MechOnly 
from  welib.kalman.kalmanfilter import KalmanFilter

MyDir=os.path.dirname(__file__)

def main():
    # --- Main parameters
    tRange=[0,100]  # Time range for simulation [s]
    nUnderSamp=1    # 1: use same time steps as measurements, >1: undersample
    NoiseRFactor=0  # Add noise to measurements. 0: no noise
    algo     = '3'  # Define algorithm kind for splitting qh into states and inputs

    # --- Parameters for the simulation that we will use as "measurements"
    simFile=os.path.join(MyDir,'../../../data/Monopile/MT100_LoadsMotions.csv')
    # Mappign to rename columns present in simulation file, to be uses as measurements and clean states (for comparison)
    colMap={
            'q'    : 'GX_[m]' ,
            'qdot' : 'GV_[m/s]' ,
            'TTacc': 'GA_[m/s^2]',
            'GF'   : 'GF_[N]', 
            'eta'  : '{Eta_[m]}-50', 
            'M_sb' : 'M_sb_[Nm]', 
            'qh'   : '{Eta_[m]}-50', 
            }

    # --- Parameters that are a function of the structure and ocean conditions
    # TODO put me in a function, based on D,CD,CM,m,EI,h etc
    # Hydro inputs
    k_h    = -231924.27 # Hydrodynamic factor GF_hydro= k_h   * qh
    #  Mechanical matrices
    MM = np.array([[4.936E+05]]) # Struct+AddedMass
    KK = np.array([[2.395E+07]])
    CC = KK*0.000
    k_h_mm1 = k_h/MM[0]
    #print('k_h/m',k_h_mm1)


    # --- Define names of physical states, augmented states, measurements, and inputs
    sStates     = np.array(['q','qdot'] )
    sAug        = np.array(['qh'] )  # Hydrodynamic 
    sMeas       = np.array(['TTacc'])
    sInp        = np.array(['qh_mean'])
    sStore      = np.array(['eta','GF','M_sb','qh_avg','qh'])

    # --- Setup state matrices, problem specific!
    # Generalized force from augmented states (qddot to p)
    GF_hydro_z = np.array([[k_h]]) # GF_h_z = \int phi(z) p_h(z) dz , should be (nDOF x nP)
    # State matrix A, and empty inputs/outputs B,C,D
    A,B,C,D = BuildSystem_Linear_MechOnly(MM, CC, KK, nP=len(sAug), nU=len(sInp), nY=len(sMeas), Fp=GF_hydro_z)
    # Setting B,C,D
    B[1,0]=k_h_mm1 
    C[:,:]=A[1,:] # Acceleration only
    D[0,0]=k_h_mm1

    # --- Initialize an empty Kalman Filter and set the state matrices
    KF = KalmanFilter(sX0=sStates, sXa=sAug, sU=sInp, sY=sMeas, sS=sStore)
    KF.setMat(A,B,C,D)

    # --- Loading "Measurements"
    # - Reference file is opened
    # - Measurements are extracted from it
    # - Other signals are extracted from the file, for comparison with estimates. These are referred as "clean" values
    # - Estimate sigmas from measurements (overriden in next section)
    KF.initFromSimulation(simFile, nUnderSamp=nUnderSamp, tRange=tRange, colMap=colMap, timeCol='Time_[s]')

    # --- Process and measurement uncertainties (standard deviation sigma)
    # Important parameters defnining uncertainties on the signals
    KF.sigX['qh']    = 40.2  # 2.1
    KF.sigX['q']     = 0.021*2  # 0.021
    KF.sigX['qdot']  = 0.064*2 # 0.064
    KF.sigY['TTacc'] = 0.440*2  # 0.440
    # 
    # --- Storage for plot, convert sigmas to covariance matrices (KF.R and KF.Q)
    KF.prepareTimeStepping()
    # --- Prepare measurements - Create noisy measurements
    KF.setYFromClean(R=KF.R, NoiseRFactor=NoiseRFactor)

    # --- Initial conditions
    x = KF.initFromClean()

    # --- Time loop
    qh_DC=0 # Part of qh that is put in the input (u) instead of the augmented states
    T_avg=5 # Averaging every n seconds for Option 1a
    n_avg=int(T_avg/KF.dt)
    for it in range(0,KF.nt-1):    
        # --- "Measurements"
        y  = KF.Y.iloc[it,:]
        # --- Inputs
        u = KF.U_clean.iloc[it,:]
        if it>1:
            if algo=='1a':
                iStart = max(it-n_avg,0) # Option 1a, average over n previous seconds
            else: # 1b
                iStart = 0               # Option 1b, average from t=0 to current time
            iEnd   = max(it-1,0)
            if iEnd-iStart>1:
                if algo=='1a' or algo=='1b':
                    qh_DC = np.mean(KF.X_hat['qh'][iStart:iEnd]) # Option 1: take mean from previous time
                elif algo=='2':
                    qh_DC = x[KF.iX['qh']] # Option 2: use previous value
                else:
                    qh_DC = 0  # Option 3: do not use an input
                u[0] =  qh_DC

        # --- Predictions of next time step based on current time step
        x,KF.P,_ = KF.estimateTimeStep(u,y,x,KF.P,KF.Q,KF.R)

        # --- Estimate derived states
        qh = x[KF.iX['qh']]
        eta  = qh+qh_DC         # NOTE: not really eta, a phase shift is missing
        GF   = (qh+qh_DC) * k_h
        M_sb = 0                # TODO compute M_sb based on qddot and estimate of Fhydro 

        # --- Store states and measurements
        KF.Y_hat.iloc[it+1,:]   = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)
        KF.X_hat.iloc[it+1,:]   = x
        # --- Store extra info
        KF.S_hat['eta'   ][it+1] = eta
        KF.S_hat['M_sb'  ][it+1] = M_sb
        KF.S_hat['GF'    ][it+1] = GF
        KF.S_hat['qh_avg'][it+1] = qh_DC
        KF.S_hat['qh'    ][it+1] = qh
        # --- Propagation to next time step
        if np.mod(it,500) == 0:
            print('Time step %8.0f t=%10.3f ' % (it,KF.time[it]))

    # --- 
    KF.plot_X()
    # KF.plot_Y()
    KF.plot_S()

    return KF

if __name__ == '__main__':
    main()
#     KF.print_sigmas()
#     print(KF)
    plt.show()

if __name__ == '__test__':
    main()
