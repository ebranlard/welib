""" Documentation
OpenFAST Lin uses:
  - 2 mechanical DOF:    'u' (tower bending)   'psi' (shaft rotation)
  - 4 measurements:               'TT acc',  'omega_rotor' ,    'Mgen' , 'Pitch' 
  - 5 states:             'u',  'azimuth', 'udot', 'omega_rotor', 'Qaero'
  - 3 inputs:             'T',  'Qgen', 'pitch'
 The estimated states are compared to the simulation at the end
 
YAMS uses:
  - 2 mechanical DOF:    'u' (tower bending)   'psi' (shaft rotation)
  - 3 measurements:      'TT acc',  'omega_rotor' ,    'Mgen'
  - 7 states:            'u',  'azimuth', 'udot', 'omega_rotor', 'T','Qaero' 'Qgen'
 Used to be called onshore_YAMS/300_Kalman_2DOF_7States
 Used to be called onshore_FASTlin/301_Kalman_2DOF_5States
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from welib.essentials import *
from welib.kalman.KF_TN    import KalmanFilterTNSim 
from welib.kalman.KF_TNLin import KalmanFilterTNLinSim
from welib.fast.FASTLin import FASTLin

from welib.tools.fatigue import eq_load

import pytest

scriptDir = os.path.dirname(__file__)

def main(bYAMS=True, StateModel='nt1_nx5', test=False):
    # Options for 7 states
    if test:
        tRange=[200,300]
    else:
        tRange=[0,700]
    nUnderSamp=5
    bExport=False
    # bExport=True
    # bNoise=True
    bNoise=False
    bMoreNoise=False
    bFilterAcc=False  # FILTER ACC IMPROVES SPECTRA BUT INCREASE REL ERR OF My
    # bFilterAcc=True
    nFilt=15
    Qgen_LSS = True

    # --- Parameters

    sPref=''
    NoiseRFactor=0
    if bNoise:
        sPref+='_Noise'
        NoiseRFactor=1/10
    if bMoreNoise:
        sPref+='More'
        NoiseRFactor=1/5
    if bFilterAcc:
        sPref+='_FilterAcc'+str(nFilt)

    OutDir       = os.path.join(scriptDir, './../../data/NREL5MW/onshore/_kalman/')
    FstFile      = os.path.join(scriptDir, '../../../data/NREL5MW/onshore/Hat.fst')
    linFile      = os.path.join(scriptDir, '../../../data/NREL5MW/onshore/ws_5_lin.1.lin')
    linStateFile = os.path.join(scriptDir, '../../../data/NREL5MW/onshore/ws_5_lin_FASTLin_2DOF.pkl') # Will be generated from linFile if not existing
    aeroMapFile  = os.path.join(scriptDir, '../../../data/NREL5MW/NREL5MW_CPCTCQ.txt')


    MeasFile   = FstFile.replace('.fst','.outb')
    Case       = os.path.basename(FstFile.replace('.fst',''))
    if bYAMS:
        OutputFile = OutDir+Case+'_NREL5MW_{:s}'.format('YAMS')+sPref+'.csv'
    else:
        OutputFile = OutDir+Case+'_NREL5MW_{:s}'.format(StateModel)+sPref+'.csv'
    if not test:
        INFO('OutputFile: ' + OutputFile)

    # --- Sigmas
    useStdFromMeas = False  # if True, the sigma are estimated based on the std from meas
    if not useStdFromMeas:
        # States
        sigs = {'x':{}, 'y':{}, 'Q':{}}
        sigs['x']['ut1']    = 1.0
        sigs['x']['psi']    = 0.1
        sigs['x']['ut1dot'] = 0.1
        sigs['x']['omega']  = 0.1
        sigs['x']['Thrust'] = 1000000
        sigs['x']['Qaero']  = 8*10**6*1.0
        sigs['x']['Qgen']   = 1.0*10**6
        sigs['x']['WS']     = 1.0
        sigs['Q'] = sigs['x'].copy()
        # Measurements - more or less half the std
        sigs['y']['TTacc'] = 0.08  # m/s^2
        sigs['y']['omega'] = 0.05 # rad/s
        sigs['y']['Qgen']  = 1*10**6
        sigs['y']['pitch'] = 2.00

    # --------------------------------------------------------------------------------}
    # --- Kalman filter estimation 
    # --------------------------------------------------------------------------------{
    with Timer('Simulation Loop'):
        if bYAMS:
            bThrustInStates=True
            KF= KalmanFilterTNSim(FstFile, MeasFile, OutputFile, aeroMapFile, bThrustInStates, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigs=sigs, bExport=bExport)
        else:

            if not os.path.exists(linStateFile):
                # We create the state file
                # sX_sel  = ['qt1FA_[m]', 'psi_rot_[rad]', 'd_qt1FA_[m/s]', 'd_psi_rot_[rad/s]']
                sED_DOF = ['7_TwFADOF1' ,'13_GeAz']
                sU_DOF  = ['HubFxN1_[N]','Qgen_[Nm]','PitchColl_[rad]']
                sY_DOF  = ['NcIMUTAxs_[m/s^2]', 'RotSpeed_[rpm]','SvDGenTq_[kNm]', 'BPitch1_[deg]']
                FL = FASTLin(linfiles=[linFile])# sX, sU, sY, sED = FL.xdescr, FL.udescr, FL.ydescr, FL.EDdescr
                Ar, Br, Cr, Dr, Mr = FL.average_subset(sU_sel=sU_DOF, sY_sel=sY_DOF, sE_sel=sED_DOF, exportFile=linStateFile, baseDict={'model':'TNSB'})
                print('A:\n',Ar)
                print('D:\n',Dr)
                print('M:\n',Mr)


            KF= KalmanFilterTNLinSim(FstFile, MeasFile, OutputFile, aeroMapFile, linStateFile, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigs=sigs, bExport=bExport, StateModel=StateModel, Qgen_LSS=Qgen_LSS, ThrustHack=True)
    # --------------------------------------------------------------------------------}
    # --- PostPro  
    # --------------------------------------------------------------------------------{
    def Leq(t,y,m=5):
        return eq_load(y, m=m, neq=t[-1]-t[0])[0][0]
    def SNR(y):
        return np.mean(y**2)/np.std(y)**2
    # --- Leq
    for j in [2,5]:
        print('Leq  ref: {:.2f} - est: {:.2f}'.format(Leq(KF.time,KF.M_ref[j]),Leq(KF.time,KF.M_sim[j])))
    # --- Signal-to-noise ratio
#     for j,s in enumerate(KF.sY):
#         print('SNR {:s} - clean: {:.2f} - meas {:.2f} - est {:.2f}'.format(s, SNR(KF.Y_clean[j,:]), SNR(KF.Y[j,:]), SNR(KF.Y_hat[j,:])))
    if not test:
        # --- Compare Measurements
        KF.plot_Y()
        # --- Compare Intermediate values
        KF.plot_S()
        # --- Compare States
        KF.plot_X()
        #
        KF.plot_moments()
    fig, stats = KF.plot_summary()
    #if test:
    #    plt.close(fig)
    #                                                  
    return stats

def test_onshore_TNS_YAMS(test=True):
    stats = main(bYAMS=False, StateModel='nt1_nx5', test=test)
    np.testing.assert_array_less(stats['Qaero']['eps'] , 3.85)
    np.testing.assert_array_less(stats['WS']['eps']    , 3.65)
    np.testing.assert_array_less(stats['Thrust']['eps'], 3.15)
    np.testing.assert_array_less(stats['M1']['eps']    , 5.05)
    np.testing.assert_array_less(stats['M6']['eps']    , 2.75)
    np.testing.assert_array_less(stats['M9']['eps']    ,21.55)

def test_onshore_TNS_OFLin(test=True):
    stats = main(bYAMS=True, StateModel='nt1_nx5' , test=test)
    np.testing.assert_array_less(stats['Qaero']['eps'] , 3.85)
    np.testing.assert_array_less(stats['WS']['eps']    , 3.65)
    np.testing.assert_array_less(stats['Thrust']['eps'], 3.15)
    np.testing.assert_array_less(stats['M1']['eps']    , 4.85)
    np.testing.assert_array_less(stats['M6']['eps']    , 2.65)
    np.testing.assert_array_less(stats['M9']['eps']    ,21.55)

if __name__ == '__main__':
    test_onshore_TNS_YAMS (test=True)
    test_onshore_TNS_OFLin(test=True)

    plt.show()
