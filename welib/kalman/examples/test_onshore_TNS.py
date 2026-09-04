""" Documentation
 This scripts uses:
  - 2 mechanical DOF:    'u' (tower bending)   'psi' (shaft rotation)
  - 4 measurements:               'TT acc',  'omega_rotor' ,    'Mgen' , 'Pitch' 
  - 5 states:             'u',  'azimuth', 'udot', 'omega_rotor', 'Qaero'
  - 3 inputs:             'T',  'Qgen', 'pitch'
 The estimated states are compared to the simulation at the end


 Used to be called 301_Kalman_2DOF_5States
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from welib.essentials import *
from welib.kalman.TN    import KalmanFilterTNSim 
from welib.kalman.TNLin import KalmanFilterTNLinSim

import pytest

scriptDir = os.path.dirname(__file__)


class KalmanModel():
    def __init__(self,StateModel='nt1_nx7',Qgen_LSS=True):
        self.StateModel=StateModel
        self.Qgen_LSS=Qgen_LSS
        self.ThrustHack=False

        if Qgen_LSS:
            self.ColMap={
              ' ut1    ' : ' TTDspFA_[m]                   ' ,
              ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
              ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
              ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
              ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
              ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
              ' Qgen   ' : ' 97*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
              ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
              ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
              ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
            }
        else:
            self.ColMap={
              ' ut1    ' : ' TTDspFA_[m]                   ' ,
              ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
              ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
              ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
              ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
              ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
              ' Qgen   ' : ' {GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
              ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
              ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
              ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
            }

        if self.StateModel=='nt1_nx8':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Thrust','Qaero','Qgen','WS'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['pitch'])
            self.sStor       = np.array(['WS'])
            self.bWSInStates     = True
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx7':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Thrust','Qaero','Qgen'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['pitch'])
            self.sStor       = np.array(['WS'])
            self.bWSInStates     = False
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx6':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Thrust','Qaero'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['Qgen','pitch'])
            self.sStor       = np.array(['WS'])
            self.bWSInStates     = False
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx5':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Qaero'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['Thrust','Qgen','pitch'])
            self.sStor       = np.array(['Thrust','WS'])
            self.bWSInStates     = False
            self.bThrustInStates = False


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

    sPref='_Base'
    NoiseRFactor=0
    if bNoise:
        sPref+='_Noise'
        NoiseRFactor=1/10
    if bMoreNoise:
        sPref+='More'
        NoiseRFactor=1/5
    if bFilterAcc:
        sPref+='_FilterAcc'+str(nFilt)

    OutDir    = os.path.join(scriptDir, './_data_onshore_TNS/')
    FstFile   = os.path.join(scriptDir, './_data_onshore_TNS/Hat.fst')
    StateFile = os.path.join(scriptDir, './_data_onshore_TNS/NREL5MW_2DOF_ABCD_mean_NEW.dat')
    base      = os.path.join(scriptDir, './_data_onshore_TNS/NREL5MW')
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
        sigX=dict()
        sigX['ut1']    = 1.0
        sigX['psi']    = 0.1
        sigX['ut1dot'] = 0.1
        sigX['omega']  = 0.1
        sigX['Thrust'] = 1000000
        sigX['Qaero']  = 8*10**6*1.0
        sigX['Qgen']   = 1.0*10**6
        sigX['WS']     = 1.0
        # Measurements - more or less half the std
        sigY=dict()
        sigY['TTacc'] = 0.08  # m/s^2
        sigY['omega'] = 0.05 # rad/s
        sigY['Qgen']  = 1*10**6
        sigY['pitch'] = 2.00

    # --------------------------------------------------------------------------------}
    # --- Kalman filter estimation 
    # --------------------------------------------------------------------------------{
    with Timer('Simulation Loop'):
        if bYAMS:
            bThrustInStates=True
            KF= KalmanFilterTNSim(FstFile, MeasFile, OutputFile, base, bThrustInStates, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX, sigY, bExport)
        else:
            KM = KalmanModel(StateModel=StateModel, Qgen_LSS=Qgen_LSS)
            KM.ThrustHack=True
            KF= KalmanFilterTNLinSim(KM, FstFile, MeasFile, OutputFile, base, StateFile, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX, sigY, bExport)
    # --------------------------------------------------------------------------------}
    # --- PostPro  
    # --------------------------------------------------------------------------------{
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
    return
    stats = main(bYAMS=False, StateModel='nt1_nx5', test=test)
    np.testing.assert_array_less(stats['Qaero']['eps'] , 3.61)
    np.testing.assert_array_less(stats['WS']['eps']    , 3.61)
    np.testing.assert_array_less(stats['Thrust']['eps'], 3.01)
    np.testing.assert_array_less(stats['M2']['eps']    , 3.01)
    np.testing.assert_array_less(stats['M7']['eps']    , 8.21)

def test_onshore_TNS_OFLin(test=True):
    return
    stats = main(bYAMS=True, StateModel='nt1_nx5' , test=test)
    np.testing.assert_array_less(stats['Qaero']['eps'] , 3.61)
    np.testing.assert_array_less(stats['WS']['eps']    , 3.61)
    np.testing.assert_array_less(stats['Thrust']['eps'], 3.01)
    np.testing.assert_array_less(stats['M2']['eps']    , 3.01)
    np.testing.assert_array_less(stats['M7']['eps']    , 8.21)

if __name__ == '__main__':
    test_onshore_TNS_YAMS (test=True)
    test_onshore_TNS_OFLin(test=True)

    plt.show()
