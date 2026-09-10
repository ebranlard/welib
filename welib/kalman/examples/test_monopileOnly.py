"""
Monopile with no top mass excited by waves.


Std(P), error covariance::
    - std(P structure) = 10^-1: low, filter is confident in structural states
    - std(P wave)      = 10:   high compared to structure, forces
     
Gains (K)

Reduce damping zeta:
    will allow for more resonance, larger oscillation of qh for the same amount of push from the nmeasurements

"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
# Welib
from welib.essentials import *
from welib.kalman.KF_M import KalmanFilterMonopile


import pytest

scriptDir = os.path.dirname(__file__)

def main(fstFile=None, hydroShape=None, Tp=None, tRange=None, tRangeStats=None):

    # --- Main parameters
    export = False
    nUnderSamp   = 1                          # 1: use same time steps as measurements, >1: undersample
    NoiseRFactor = 0                          # Add noise to measurements. 0: no noise

    # --- Monopile Jonswap Hs=8.1 Tp=12.7
    if fstFile is None:
        tRange      = [0,12] # Time range for simulation [s]
        tRangeStats = [0,12] # Time range for stats [s]
        Tp = 12.7
        fstFile    = os.path.join(scriptDir, '../../../data/Monopile/Main_MT100_JONSWAP_UserDef.fst')
        hydroShape = os.path.join(scriptDir, '../../../data/Monopile/MT100_HydroShapeFunction_Hs=8.1_Tp=12.7_h=50.csv')

    # --- Script derived parameters
    simFile = fstFile.replace('.fst','.outb')              # Measurements

    KF = KalmanFilterMonopile()
    KF.setup_matrices(fstFilename=fstFile, hydroShape=hydroShape, Tp=Tp)

    # --- Loading "Measurements"
    # - Reference file is opened
    # - Measurements are extracted from it
    # - Other signals are extracted from the file, for comparison with estimates. These are referred as "clean" values
    # - Estimate sigmas from measurements (overriden in next section)
    KF.loadMeasurements(simFile, nUnderSamp=nUnderSamp, tRange=tRange, colMap=KF.colMap, timeCol='Time_[s]', raiseIfAbsent=True)
    KF.X_clean['qd_h'] = np.gradient(KF.X_clean['q_h'], KF.dt)

    # --- Process and measurement uncertainties (standard deviation sigma)
    # Important parameters defnining uncertainties on the signals
    # --- Storage for plot, convert sigmas to covariance matrices (KF.R and KF.Q)
    KF.prepareTimeStepping() 

    # KF.Xx.iloc[0,0] = -0.001 # Centering force
    # Based on influence in $C$.
    # Since $q_p$ is so sensitive, it needs a smaller $Q$ than $q_s$ to prevent it from dominating the filter's attention.
    # --- Process and measurement covariances
    dt_ref = 0.01 # NOTE: Q change with dt
    sigs = {'x':{}, 'y':{}, 'Q':{}}
    sigs['y']['TTacc'] = np.sqrt(1e-3)
    sigs['Q']['q_s']   = np.sqrt(KF.dt/dt_ref * 1e-6)
    sigs['Q']['q_p']   = np.sqrt(KF.dt/dt_ref * 1e-6)
    sigs['Q']['qd_s']  = np.sqrt(KF.dt/dt_ref * 1e-3)
    sigs['Q']['qd_p']  = np.sqrt(KF.dt/dt_ref * 1e-6)
    sigs['Q']['q_h']   = np.sqrt(KF.dt/dt_ref * 1e-6)
    sigs['Q']['qd_h']  = np.sqrt(KF.dt/dt_ref * KF.Sw)
    KF.setupCovariances(
            sigs=sigs,
            useDt=False, Pidentity=True, verbose=True)

    # --- Prepare measurements - Create noisy measurements
    KF.setYFromClean(R=KF.R, NoiseRFactor=NoiseRFactor)

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    # KF.print_sigmas()
    print(KF)
    print('>>> dt', KF.dt)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

    KF.timeLoop()

    # --- 
    print('-------------------------------------------------------')
    print('fstFile: ', fstFile)
    print('Qdiag  : ', np.diag(KF.Q))
    print('Rdiag  : ', np.diag(KF.R))
    print('Cmat   : ', KF.C.values)
    print(f'Tuning: zeta={KF.zeta}, qdhScale={KF.qdhScale}')

    statsDict = {}
    fig = KF.plot_X( printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
    fig = KF.plot_Y()
    # KF.plot_U()
    fig = KF.plot_S(printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
    # KF.plot_P()
    # KF.plot_K()
    # KF.plot_innovation()
    print('-------------------------------------------------------')

    KFpkl = fstFile.replace('.fst', '_KF302.pkl')
    if export:
        KF.save(KFpkl)
    # 
    return statsDict

def test_monopile():
    # --- Monopile Jonswap Hs=8.1 Tp=12.7, Default
    stats = main()
    np.testing.assert_array_less(stats['q_s']['eps'],  1.9)
    np.testing.assert_array_less(stats['M_sb']['eps'], 5.9)
    np.testing.assert_array_less(stats['F_sb']['eps'], 5.9)
    np.testing.assert_array_less(stats['eta']['eps'],  2.4)
    np.testing.assert_array_less(stats['Fhx']['eps'],  5.5)


if __name__ == '__main__':

    # --- Monopile Jonswap Hs=8.1 Tp=12.7, Default
    test_monopile()

    # --- Monopile Jonswap Hs=2.5 Tp=10
#     tRange=[0,600]; 
#     tRangeStats=[35,600]; 
#     Tp = 10.0
#     #tRange      = [0,12]  # Time range for simulation [s]
#     #tRangeStats = [0,12] # Time range for stats [s]
#     fstFile    = f'simulations/MT100/06_Jonswap/OF_Long_Hs=2.5_Tp=10_h=50.fst'; 
#     hydroShape = f'simulations/MT100/06_Jonswap/MT100_HydroShapeFunction_Hs=2.5_Tp=10_h=50.csv'
#     stats = main(fstFile =fstFile, hydroShape=hydroShape, Tp=Tp, tRange=tRange, tRangeStats=tRangeStats)

#     # --- Monopile Jonswap Hs=8.1 Tp=12.7, Long
#     tRange=[0,100]; 
#     tRangeStats=[35,600]; 
#     Tp = 12.7
#     #tRange      = [0,12]  # Time range for simulation [s]
#     #tRangeStats = [0,12] # Time range for stats [s]
#     fstFile    = os.path.join(scriptDir, '../../../data/Monopile/Main_MT100_JONSWAP_UserDef_Long.fst')
#     hydroShape = os.path.join(scriptDir, '../../../data/Monopile/MT100_HydroShapeFunction_Hs=8.1_Tp=12.7_h=50.csv')
#     stats = main(fstFile =fstFile, hydroShape=hydroShape, Tp=Tp, tRange=tRange, tRangeStats=tRangeStats)


    # --- Monopile Reg Wave
    # fstFile = '../simulations/MT100/06_Jonswap/OF.fst'
    # SSCase='_Hs=8.1_Tp=12.7_h=50'; fstFile = f'../simulations/MT100/06_Jonswap/OF_Long{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600]; 




    plt.show()
