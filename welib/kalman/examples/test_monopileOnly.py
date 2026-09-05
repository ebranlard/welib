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

from numpy import trapezoid

import pytest

scriptDir = os.path.dirname(__file__)


# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def monopileSetupFromOpenFAST(fstFile, shapes_sub=[0,4], TMIN=0.0, TMAX=None, nSubSample=1, compFile=None, tuneM=False, hydroShape=None, reconHydro=False):
    import welib.weio as weio
    #from welib.yams.models.FNSB_FAST import FASTmodel2FNSB
    from welib.yams.models.MNSB_FAST import FASTmodel2MNSB
    from welib.yams.windturbine import FASTWindTurbine
    from welib.system.mech_system import MechSystem
    from welib.hydro.wavekin import wavenumber, elevation2d

    # --- Loading a "YAMS" model
    WT = FASTmodel2MNSB(fstFile, shapes_sub=shapes_sub, shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z', assembly='manual', fixedShaft=True).WT

    print(WT)
    zBeam = np.array(sorted(WT.twr.SD.pointsMN['z'].unique()))
    zDepth = WT.twr.s_span - WT.WtrDpth
    zDepth = zBeam
    if not np.array_equal(zBeam, zDepth):
        raise Exception('Problem in model, both z are different')

    # --- Structure
    pST         = dict()
    pST['PhiU']   = WT.twr.PhiU
    pST['PhiV']   = WT.twr.PhiV
    pST['PhiK']   = WT.twr.PhiK
    pST['m']      = WT.twr.m
    pST['s_span'] = WT.twr.s_span
    pST['gravity'] = WT.gravity

    pST['z']    = WT.twr.s_span - WT.WtrDpth
    # --- Sea state and hydrodynamics
    pSS = None
    pHD = None
    if WT.Hydro:
        # --- Sea state, wave kinematics
        pSS = dict()
        pSS['rho']        = WT.WtrDens      # [kg/m3]
        pSS['WaterDepth'] = WT.WtrDpth      # [m]
        if compFile is not None:
            dfComp = weio.read(compFile).toDataFrame()
            pSS['ap']   = dfComp['Amplitude_[m]']
            pSS['fp']   = dfComp['Frequency_[Hz]']
            pSS['epsp'] = dfComp['Phase_[rad]']
        else:
            print('>>> Using Wave of amplitude 3 and period 12 for now')
            pSS['ap']   = np.array([3])     # Amplitudes
            pSS['fp']   = np.array([1/12])  # frequencies [Hz]
            pSS['epsp'] = np.array([np.pi]) # Deterministic phas

        # --- Wave Kinematics
        pSS['kp'] = wavenumber(pSS['fp'], pSS['WaterDepth'], pST['gravity']) # Wave numbers
        pSS['compFile'] = compFile

        # ---  Hydrodynamics 
        HD = WT.HD
        try:
            cprop = HD.getTab('SectionPropCyl')
        except:
            cprop = HD.getTab('SectionProp')
        D_ = cprop['PropD'].values[0]
        sprop = HD.getTab('SmplPropCyl')
        Cp_ = sprop['SimplCp'].values[0]
        Cd_ = sprop['SimplCd'].values[0]
        Ca_ = sprop['SimplCa'].values[0]
        print(f'HD props: Cp={Cp_} Cd={Cd_} Ca={Ca_} D={D_} ')
        PlaceHolderOnes = np.ones([len(zDepth),1])
        pHD = dict()
        pHD['zDepth'] = zDepth
        pHD['D']      = D_  * PlaceHolderOnes
        pHD['Cd']     = Cd_ * PlaceHolderOnes
        pHD['CM']     = (Ca_+Cp_) *PlaceHolderOnes
        pHD['Ca']     = Ca_ * PlaceHolderOnes
        pHD['Cp']     = Cp_ * PlaceHolderOnes
        pHD['m_hydro']= pSS['rho'] * np.pi / 4 * pHD['D'].ravel()**2 * pHD['Ca'].ravel()
            
        # Generalized hydro mass matrix 
        GM_hydro = np.zeros((len(WT.twr.PhiU),len(WT.twr.PhiU)))
        bWet = zDepth<=0
        for i, phi in enumerate(WT.twr.PhiU):
            phi_x = phi[0,:]
            GM_hydro[i,i] = np.trapezoid(pHD['m_hydro'][bWet] * phi_x[bWet]**2, zDepth[bWet])

        pHD['GM_hydro']= GM_hydro
            
        print('GM_hydro:\n', GM_hydro)
        for i, phi in enumerate(WT.twr.PhiU):
            # NOTE: diagonal only?
            WT.MM[i,i]+=GM_hydro[i,i]

    if tuneM:
#         raise Exception('Removed ?')
        # TODO remove this, it's application specific
        # Tuning
        factM1=1.009
        factM2=1.003
        WT.MM[0,0]*=factM1
        WT.MM[1,1]*=factM2
        WT.KK[0,0]*=factM1
        WT.KK[1,1]*=factM2


    # --- Load hydroShape
    if hydroShape is not None:
        dfH = weio.read(hydroShape).toDataFrame()
        if not np.array_equal(dfH['z_[m]'],pHD['zDepth']):
            raise Exception('z depth different with shape function and hydrodyn')
        pHD['phi']  = dfH['phi_[Ns/m^2]'].values # Used to be called k_h_z
        pHD['phit'] = dfH['phit_[-]'].values
        bWet = zDepth<=0
        k_h = np.zeros(len(WT.twr.PhiU))
        for i, phi in enumerate(WT.twr.PhiU):
            phi_x = phi[0,:]
            k_h[i] = np.trapezoid(pHD['phi'][bWet] * phi_x[bWet], zDepth[bWet])
        print('k_h     : ', k_h)
        pHD['k_h'] = k_h


    # --- OPENFAST as a reference simulation
    outFile = fstFile.replace('.fst','.outb')
    if not os.path.exists(outFile):
        outFile = fstFile.replace('.fst','.out')
    if not os.path.exists(outFile):
        WARN('Cannot setup reference, out file does not exist: '+outFile)
        ref = None
    else:
        ref = dict()
        # --- Reading input loads
        df = weio.read(outFile).toDataFrame()
        df = df.iloc[::nSubSample] # SubSampling for shorter comp time
        df = df[df['Time_[s]']<TMAX]   # Limiting time
        df = df[df['Time_[s]']>=TMIN]   # Limiting time
        vTime = df['Time_[s]'].values
        print('nSteps  : ',len(vTime))
        # --- State
        if shapes_sub==[0,4]:
            ref['state']   = df[ ['Q_Sg_[m]', 'Q_P_[rad]', 'QD_Sg_[m/s]', 'QD_P_[rad/s]']].values.T
            ref['state_d'] = df[ ['QD_Sg_[m/s]', 'QD_P_[rad/s]', 'QD2_Sg_[m/s^2]', 'QD2_P_[rad/s^2]']].values.T
        elif shapes_sub==[0]:
            ref['state']   = df[ ['Q_Sg_[m]',    'QD_Sg_[m/s]']].values.T
            ref['state_d'] = df[ ['QD_Sg_[m/s]', 'QD2_Sg_[m/s^2]']].values.T
        else:
            raise NotImplementedError()
        ## --- Nicknames
        if 'Wave1Elev_[m]' in df.columns:
            ref['eta']         = df['Wave1Elev_[m]']
        else:
            ref['eta']         = df['Time_[s]'] * 0.0
        if 'HydroFxi_[N]' in df.columns:
            ref['HydroFx']     = df['HydroFxi_[N]'].values
        else:
            ref['HydroFx']     = df['Time_[s]'].values * 0.0
        # --- Section loads Ref
        # SubDyn secton outputs
        zBeam, F_sec, r_sec =  WT.twr.SD.beamSecOutputs(df, verbose=False)
        ref['df']    = df
        ref['z']     = zBeam
        ref['F_sec'] = F_sec
        ref['r_sec'] = r_sec

    if ref is not None and WT.Hydro:
        #print('Wave freq {}  amplitude {} '.format(np.max(pSS['fp']),np.max(pSS['ap'])))
        eta_sim = elevation2d(pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], vTime, x=0)
        dt = vTime[1]-vTime[0]
        eta_dot = np.gradient(eta_sim, dt) # Velocity state
        pSS['eta_time'] = vTime
        pSS['eta']      = eta_sim
        pSS['eta_dot']  = eta_dot

    # --- Setup Sys
    Sys = MechSystem(WT.MM, WT.DD, WT.KK)
    Sys.setStateInitialConditions(WT.z0.values)
    # --- Set external loads time series
    ForceFunction = lambda t,q,qd : monopileGF(t, q, qd, pST=pST, pSS=pSS, pHD=pHD, reconHydro=reconHydro)[0]
    Sys.setForceFunction( ForceFunction )

    return pST, pSS, pHD, Sys, WT, ref



def main(fstFile=None, hydroShape=None, Tp=None, tRange=None, tRangeStats=None):

    # --- Main parameters
    export = False
    nUnderSamp   = 1                          # 1: use same time steps as measurements, >1: undersample
    NoiseRFactor = 0                          # Add noise to measurements. 0: no noise
    shapes_sub =[0,4]

    # --- Monopile Jonswap Hs=8.1 Tp=12.7
    if fstFile is None:
        tRange      = [0,12] # Time range for simulation [s]
        tRangeStats = [0,12] # Time range for stats [s]
        Tp = 12.7
        fstFile    = os.path.join(scriptDir, '../../../data/Monopile/Main_MT100_JONSWAP_UserDef.fst')
        hydroShape = os.path.join(scriptDir, '../../../data/Monopile/MT100_HydroShapeFunction_Hs=8.1_Tp=12.7_h=50.csv')

    # --- Script derived parameters
    simFile = fstFile.replace('.fst','.outb')              # Measurements

    # --- Loading Lin model
    pST, pSS, pHD, Sys, WT, ref = monopileSetupFromOpenFAST(fstFile, shapes_sub=shapes_sub, TMIN=tRange[0], TMAX=tRange[1], nSubSample=nUnderSamp, hydroShape=hydroShape, reconHydro=True)
    # To avoid any kind of cheating
    del pSS['kp']
    del pSS['fp']
    del pSS['ap']

    # --- Parameters that are a function of the structure and ocean conditions
    zDepth =pST['z']

    # --- Define names of physical states, augmented states, measurements, and inputs
    sStates = ['q_s','q_p','qd_s', 'qd_p']
    sAug    = ['q_h', 'qd_h']
    sMeas   = ['TTacc']
    sMeas  += ['q_p']
    sInp    = ['w'] # White noise
    sStore  = ['M_sb','F_sb', 'eta', 'Fhx']

    # --- Initialize an empty Kalman Filter 
    KF = KalmanFilter(sX0=sStates, sXa=sAug, sU=sInp, sY=sMeas, sS=sStore)

    # --- Setup state matrices, problem specific!
    # State matrix A from MCK
    # Empty inputs/outputs B,C,D
    A,B,C,D = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK, nP=len(sAug), nU=len(sInp), nY=len(sMeas), Fp=None)

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
    qdhScale =1.0

    Minv = np.linalg.inv(WT.MM)
    IQD   =[KF.iX['qd_s'], KF.iX['qd_p']]
    A[IQD, KF.iX['qd_h']]  = Minv @ (pHD['k_h'][0], pHD['k_h'][1])*qdhScale  # qd_h influence in mech DOF

    B[KF.iX['qd_h'], KF.iU['w']] = 1 # White noise

    C[KF.iY['TTacc'], :] = A[KF.iX['qd_s'],:] # TTacc is assumed to be qdd_s
    D[KF.iY['TTacc'], :] = B[KF.iX['qd_s'],:] # TTacc is assumed to be qdd_s
    C[KF.iY['q_p'], KF.iX['q_p']] = 1


    # --- HYDRO STATE EQ - See Script 180
    if Tp==12.7:
        Sw= 2.3835e-01
    elif Tp==10.0:
        Sw= 2.3835e-01/2
    else:
        raise NotImplementedError()
    omega_p = 2*np.pi/Tp
    zeta = 0.12 
    print('omega_p^2', omega_p**2, '2 zeta omega_p', 2*zeta*omega_p)
    A[KF.iX['q_h'], KF.iX['qd_h']]  = 1
    A[KF.iX['qd_h'], KF.iX['q_h']]  = -omega_p**2
    A[KF.iX['qd_h'], KF.iX['qd_h']]  = -2*zeta*omega_p


    KF.setMat(A,B,C,D)
    # import pdb; pdb.set_trace()

    # --- Loading "Measurements"
    # - Reference file is opened
    # - Measurements are extracted from it
    # - Other signals are extracted from the file, for comparison with estimates. These are referred as "clean" values
    # - Estimate sigmas from measurements (overriden in next section)
    colMap={
            'q_s'  : 'Q_Sg_[m]' ,
            'qd_s' : 'QD_Sg_[m/s]' ,
            'TTacc ' : 'NcIMUTAxs_[m/s^2]' ,
            'eta'    : 'Wave1Elev_[m]', 
            'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
            'Fhx'    : 'HydroFxi_[N]',
            'F_sb'   : '-ReactFXss_[N]',
            'M_sb'   : '-ReactMYss_[N*m]',
            'q_p'    : 'Q_P_[rad]' ,
            'qd_p'   : 'QD_P_[rad/s]'
        }

    KF.initFromSimulation(simFile, nUnderSamp=nUnderSamp, tRange=tRange, colMap=colMap, timeCol='Time_[s]', raiseIfAbsent=True)
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
    KF.Q[5,5] = KF.dt/dt_ref * Sw     # eta_dot

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
        p_hydro = pHD['phi'] * eta_dot # p_h = k_h(z) q_h(t)
        p_hydro[zDepth>0] = 0 # safety, shoudn't be necessray

        x_dot = np.dot(KF.A, x) + np.dot(KF.B, u)
        p_ext      = np.zeros((3,len(zDepth)))
        p_ext[0,:] = p_hydro

        x_q   = np.array([x[0],x[1]])
        xd_q  = np.array([x_dot[0], x_dot[1]])
        xdd_q = np.array([x_dot[2], x_dot[3]])

        ## Top loads
        #GF, outH = monopileGF(t, q=x_q, qd=xd_q, pST, pSS, pHD, qdd=xdd_q, reconHydro=True, nOut=[0], eta=None, eta_dot=eta_dot)
        F_top = np.array((0.,0.,0.))
        M_top = np.array((0.,0.,0.))
        a_ext = np.array((0.,0.,-WT.gravity)) # external acceleration (gravity/earthquake)
        F_sec, M_sec, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.twr.s_span, WT.twr.PhiU, WT.twr.PhiV, WT.twr.m, a_ext=a_ext, corrections=0, PhiK=WT.twr.PhiK)

        xdd_q *=0
        F_sec_h, M_sec_h, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.twr.s_span, WT.twr.PhiU, WT.twr.PhiV, WT.twr.m, a_ext = a_ext, PhiK=WT.twr.PhiK)

        
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
    print(f'Tuning: zeta={zeta}, qdhScale={qdhScale}')
    # KF.print_sigmas()
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

    # --- Monopile Jonswap Hs=8.1 Tp=12.7, Long
    tRange=[0,100]; 
    tRangeStats=[35,600]; 
    Tp = 12.7
    #tRange      = [0,12]  # Time range for simulation [s]
    #tRangeStats = [0,12] # Time range for stats [s]
    fstFile    = os.path.join(scriptDir, '../../../data/Monopile/Main_MT100_JONSWAP_UserDef_Long.fst')
    hydroShape = os.path.join(scriptDir, '../../../data/Monopile/MT100_HydroShapeFunction_Hs=8.1_Tp=12.7_h=50.csv')
    stats = main(fstFile =fstFile, hydroShape=hydroShape, Tp=Tp, tRange=tRange, tRangeStats=tRangeStats)


    # --- Monopile Reg Wave
    # fstFile = '../simulations/MT100/06_Jonswap/OF.fst'
    # SSCase='_Hs=8.1_Tp=12.7_h=50'; fstFile = f'../simulations/MT100/06_Jonswap/OF_Long{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600]; 




    plt.show()
