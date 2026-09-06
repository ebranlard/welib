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
from welib.weio.fast_linearization_file import FASTLinearizationFile
from welib.kalman.kalman import BuildSystem_Linear_MechOnly 
from welib.kalman.kalmanfilter import KalmanFilter
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions
from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.yams.windturbine import monopileSetupFromOpenFAST

from welib.kalman.kalman_model import AugmentedLinModel

from numpy import trapezoid

import pytest

scriptDir = os.path.dirname(__file__)


# --------------------------------------------------------------------------------}
# -- Augmented Linear Model 
# --------------------------------------------------------------------------------{
class KalmanModelMonopile(AugmentedLinModel):

    def __init__(KM, *args, **kwargs):
        """ 
        See AugmentedLinModel for main data (e.g. sX, A)
        Additional data specific to this model:
          - WT
          - zDepth
          - pHD
          - zeta, qdhScale
        """
        AugmentedLinModel.__init__(KM)

        # --- Define names of physical states, augmented states, measurements, and inputs
        KM.sQ  = ['q_s','q_p','qd_s', 'qd_p']
        KM.sQa = ['q_h', 'qd_h']
        KM.sY  = ['TTacc']
        KM.sY += ['q_p']
        KM.sU  = ['w'] # White noise
        KM.sS  = ['M_sb','F_sb', 'eta', 'Fhx']

        KM.setup(*args, **kwargs)

    def setup(KM, fstFilename=None, shapes_sub=None, hydroShape=None, Tp=None,
              zeta =0.12, qdhScale=1,#Tuning
              ):
        # --- Default arguments
        if shapes_sub is None:
            shapes_sub =[0,4]
        
        # --- Method 1: using MTNSB
        WT = FASTmodel2MTNSB(fstFilename, shapes_sub=shapes_sub, shapes_twr=[], shapes_bld=[], bStiffening=True, main_axis='z', fixedShaft=True, algo='OpenFAST').WT
        if WT.pSS is not None:
            #NOTE('Setting Components', compFile)
            #WT.SS_setComponents(compFile)
            #NOTE('Setting Compute Eta')
            #WT.SS_computeEta(dfRef['Time_[s]'])
            if hydroShape is not None:
                WT.HD_setShapeFunction(hydroShape)
        GM_hydro = WT.pHD['GM_hydro']
        print('GM_hydro:\n', GM_hydro)
        for i, phi in enumerate(WT.fnd.PhiU):
            # NOTE: diagonal only?
            WT.MM[i,i]+=GM_hydro[i,i]

        pHD = WT.pHD
        # --- Method 2: using Legacy TNSB
        # pST, pSS, pHD, Sys, WT, ref = monopileSetupFromOpenFAST(fstFilename, shapes_sub=shapes_sub, TMIN=tRange[0], TMAX=tRange[1], nSubSample=nUnderSamp, hydroShape=hydroShape, reconHydro=True)
        # WT.fnd = WT.twr

        KM.WT = WT
        KM.pHD = pHD

        # --- Parameters that are a function of the structure and ocean conditions
        KM.zDepth = WT.fnd.s_span - WT.WtrDpth

        # --- Setup state matrices, problem specific!
        # State matrix A from MCK
        # Empty inputs/outputs B,C,D
        A,B,C,D = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK, nP=len(KM.sQa), nU=len(KM.sU), nY=len(KM.sY), Fp=None)

        # print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> USING A LIN')
        # linFile = '../simulations/MT100/00_EVA/OF_NoHydro.2.lin'
        # linFileH = '../simulations/MT100/00_EVA/OF.3.lin'
        # SX = ['PtfmSurge_[m]', 'PtfmPitch_[rad]', 'd_PtfmSurge_[m/s]', 'd_PtfmPitch_[rad/s]']
        # dfA = FASTLinearizationFile(linFile).toDataFrame()['A']
        # Alin = dfA.loc[SX,SX]
        # dfAH= FASTLinearizationFile(linFileH).toDataFrame()['A']
        # AlinH = dfAH.loc[SX,SX]
        # A[:4,:4] = Alin.values[:4,:4]
        # A[:4,:4] = AlinH.values[:4,:4]
        KM.qdhScale = qdhScale

        Minv = np.linalg.inv(WT.MM)
        IQD   =[KM.iX['qd_s'], KM.iX['qd_p']]
        A[IQD, KM.iX['qd_h']]  = Minv @ (pHD['k_h'][0], pHD['k_h'][1])*KM.qdhScale  # qd_h influence in mech DOF

        B[KM.iX['qd_h'], KM.iU['w']] = 1 # White noise

        C[KM.iY['TTacc'], :] = A[KM.iX['qd_s'],:] # TTacc is assumed to be qdd_s
        D[KM.iY['TTacc'], :] = B[KM.iX['qd_s'],:] # TTacc is assumed to be qdd_s
        C[KM.iY['q_p'], KM.iX['q_p']] = 1


        # --- HYDRO STATE EQ - See Script 180
        if Tp==12.7:
            KM.Sw= 2.3835e-01
        elif Tp==10.0:
            KM.Sw= 2.3835e-01/2
        else:
            raise NotImplementedError()
        KM.omega_p = 2*np.pi/Tp
        KM.zeta = zeta
        print('omega_p^2', KM.omega_p**2, '2 zeta omega_p', 2*KM.zeta*KM.omega_p)
        A[KM.iX['q_h'], KM.iX['qd_h']]  = 1
        A[KM.iX['qd_h'], KM.iX['q_h']]  = -KM.omega_p**2
        A[KM.iX['qd_h'], KM.iX['qd_h']]  = -2*KM.zeta*KM.omega_p

        KM.A, KM.B, KM.C, KM.D = A, B, C, D

        KM.colMap={
                'q_s'    : 'Q_Sg_[m]' ,
                'qd_s'   : 'QD_Sg_[m/s]' ,
                'TTacc ' : 'NcIMUTAxs_[m/s^2]' ,
                'eta'    : 'Wave1Elev_[m]', 
                'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                'Fhx'    : 'HydroFxi_[N]',
                'F_sb'   : '-ReactFXss_[N]',
                'M_sb'   : '-ReactMYss_[N*m]',
                'q_p'    : 'Q_P_[rad]' ,
                'qd_p'   : 'QD_P_[rad/s]'
            }




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

    KM = KalmanModelMonopile(fstFilename = fstFile, hydroShape=hydroShape, Tp=Tp)

    # --- Initialize an empty Kalman Filter 
    KF = KalmanFilter(KM=KM)

    # --- Loading "Measurements"
    # - Reference file is opened
    # - Measurements are extracted from it
    # - Other signals are extracted from the file, for comparison with estimates. These are referred as "clean" values
    # - Estimate sigmas from measurements (overriden in next section)

    KF.initFromSimulation(simFile, nUnderSamp=nUnderSamp, tRange=tRange, colMap=KM.colMap, timeCol='Time_[s]', raiseIfAbsent=True)
    KF.X_clean['qd_h'] = np.gradient(KF.X_clean['q_h'], KF.dt)



    # --- Process and measurement uncertainties (standard deviation sigma)
    # Important parameters defnining uncertainties on the signals
    # --- Storage for plot, convert sigmas to covariance matrices (KF.R and KF.Q)
    KF.prepareTimeStepping() 

    # KF.Xx.iloc[0,0] = -0.001 # Centering force
    # Based on influence in $C$.
    # Since $q_p$ is so sensitive, it needs a smaller $Q$ than $q_s$ to prevent it from dominating the filter's attention.
    # --- Process and measurement covariances
    KF.sigmasFromClean(factor=1, dt=None)

    KF.setupCovariances(useDt=False, Pidentity=True)
    dt_ref = 0.01 # NOTE: Q change with dt
    KF.R[0,0] =                1e-3   # TTacc
    KF.Q[0,0] = KF.dt/dt_ref * 1e-6   # q_p
    KF.Q[1,1] = KF.dt/dt_ref * 1e-6   # q_s
    KF.Q[2,2] = KF.dt/dt_ref * 1e-3   # qd_p
    KF.Q[3,3] = KF.dt/dt_ref * 1e-6   # qd_s
    KF.Q[4,4] = KF.dt/dt_ref * 1e-6   # eta
    KF.Q[5,5] = KF.dt/dt_ref * KM.Sw  # eta_dot

    KF.print_sigmas()

    # -------------------------------------------------------
    # fstFile:  ../simulations/MT100/06_Jonswap/OF.fst
    # Qdiag  :  [1.0e-06 1.0e-06 1.0e-03 1.0e-05 1.0e-03 2.3e-01]
    # Rdiag  :  [1.000e-03 2.704e-07]
    # Cmat   :  [[ 1.0346e+01 -1.7264e+03  6.2078e-02 -1.0358e+01  0.0000e+00  4.2735e-01]
    #  [ 0.0000e+00  1.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00]]
    # Tuning: zeta=0.12, qdhScale=1.0
    # q_s        σ_est/σ_ref = 0.999 - ε=0.2% - R²=1.000
    # q_p        σ_est/σ_ref = 0.997 - ε=0.2% - R²=1.000
    # qd_s       σ_est/σ_ref = 1.027 - ε=0.5% - R²=0.997
    # qd_p       σ_est/σ_ref = 1.475 - ε=2.8% - R²=0.810
    # q_h        σ_est/σ_ref = 0.933 - ε=1.8% - R²=0.974
    # qd_h       σ_est/σ_ref = 0.909 - ε=1.7% - R²=0.975
    # M_sb       σ_est/σ_ref = 1.033 - ε=0.9% - R²=0.993
    # F_sb       σ_est/σ_ref = 1.021 - ε=2.4% - R²=0.952
    # eta        σ_est/σ_ref = 0.933 - ε=1.8% - R²=0.974
    # Fhx        σ_est/σ_ref = 1.015 - ε=2.4% - R²=0.951
    # -------------------------------------------------------
    # -------------------------------------------------------
    # fstFile:  ../simulations/MT100/06_Jonswap/OF.fst
    # Qdiag  :  [1.0e-06 1.0e-06 1.0e-06 1.0e-06 1.0e-06 2.3e-01]
    # Rdiag  :  [1.000e-03 2.704e-07]
    # Cmat   :  [[ 1.0346e+01 -1.7264e+03  6.2078e-02 -1.0358e+01  0.0000e+00  4.2735e-01]
    #  [ 0.0000e+00  1.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00]]
    # Tuning: zeta=0.12, qdhScale=1.0, inclin=True
    # q_s        σ_est/σ_ref = 0.996 - ε=0.2% - R²=1.000
    # q_p        σ_est/σ_ref = 0.997 - ε=0.2% - R²=1.000
    # qd_s       σ_est/σ_ref = 1.000 - ε=0.2% - R²=0.999
    # qd_p       σ_est/σ_ref = 1.422 - ε=3.6% - R²=0.357
    # q_h        σ_est/σ_ref = 0.940 - ε=2.1% - R²=0.967
    # qd_h       σ_est/σ_ref = 0.911 - ε=1.7% - R²=0.975
    # M_sb       σ_est/σ_ref = 1.027 - ε=0.8% - R²=0.995
    # F_sb       σ_est/σ_ref = 1.019 - ε=2.3% - R²=0.953
    # eta        σ_est/σ_ref = 0.940 - ε=2.1% - R²=0.967
    # Fhx        σ_est/σ_ref = 1.016 - ε=2.4% - R²=0.951
    # -------------------------------------------------------


    # --- Prepare measurements - Create noisy measurements
    KF.setYFromClean(R=KF.R, NoiseRFactor=NoiseRFactor)

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    # KF.print_sigmas()
    print(KF)
    print('>>> dt', KF.dt)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')


    # --- Initial conditions
    x = KF.initFromClean()
    KF.X_hat.iloc[0,:] = x

    WT = KM.WT

    # --- Time loop
    for it in range(0,KF.nt-1):    
        # --- "Measurements"
        y  = KF.Y.iloc[it,:].values
        # --- Inputs
        u = KF.U_clean.iloc[it,:].values

        # --- Predictions of next time step based on current time step
        x, KF.P, _ = KF.estimateTimeStep(u, y, x, KF.P)

        # --- Estimate Generalized hydro force and bending moment (calc output)
        q_h     = x[KF.iX['q_h']]  # eta
        qd_h    = x[KF.iX['qd_h']] # eta_dot
        eta     = q_h
        eta_dot = qd_h
        p_hydro = KM.pHD['phi'] * eta_dot # p_h = k_h(z) q_h(t)
        p_hydro[KM.zDepth>0] = 0 # safety, shoudn't be necessray

        x_dot = np.dot(KF.A, x) + np.dot(KF.B, u)
        p_ext      = np.zeros((3,len(KM.zDepth)))
        p_ext[0,:] = p_hydro

        x_q   = np.array([x[0],x[1]])
        xd_q  = np.array([x_dot[0], x_dot[1]])
        xdd_q = np.array([x_dot[2], x_dot[3]])

        ## Top loads
        F_top = np.array((0.,0.,0.))
        M_top = np.array((0.,0.,0.))
        a_ext = np.array((0.,0.,-WT.gravity)) # external acceleration (gravity/earthquake)
        F_sec, M_sec, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext=a_ext, corrections=0, PhiK=WT.fnd.PhiK)

        xdd_q *=0
        F_sec_h, M_sec_h, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext = a_ext, PhiK=WT.fnd.PhiK)

        
        # --- Store extra info
        KF.S_hat.at[it+1, 'M_sb']   = M_sec[1,0]
        KF.S_hat.at[it+1, 'F_sb']   = F_sec[0,0]
        KF.S_hat.at[it+1, 'eta' ]   = eta
        KF.S_hat.at[it+1, 'Fhx']    = F_sec_h[0,0]
        
        # --- Propagation to next time step
        if np.mod(it,500) == 0:
            print('Time step %8.0f t=%10.3f ' % (it,KF.time[it]))

    # --- 
    print('-------------------------------------------------------')
    print('fstFile: ', fstFile)
    print('Qdiag  : ', np.diag(KF.Q))
    print('Rdiag  : ', np.diag(KF.R))
    print('Cmat   : ', KF.C.values)
    print(f'Tuning: zeta={KM.zeta}, qdhScale={KM.qdhScale}')

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
