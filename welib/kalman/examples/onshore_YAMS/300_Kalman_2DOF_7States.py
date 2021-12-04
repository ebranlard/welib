""" Documentation
 This scripts uses:
  - 2 mechanical DOF:    'u' (tower bending)   'psi' (shaft rotation)
  - 3 measurements:               'TT acc',  'omega_rotor' ,    'Mgen'
  - 7 states:             'u',  'azimuth', 'udot', 'omega_rotor', 'T','Qaero' 'Qgen'
 The estimated states are compared to the simulation at the end
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import welib.weio # https://github.com/ebranlard/weio

from wtDigiTwin.kalman.TN import KalmanFilterTN, KalmanFilterTNSim
from wtDigiTwin.tools.colors import *
from wtDigiTwin.tools.figure import *
from wtDigiTwin.tools.tictoc import Timer
from wtDigiTwin.tools.clean_exceptions import *
from wtDigiTwin.tools.fatigue import eq_load

# --- General Parameters and State specific parameters
bThrustInStates=True
OutDir = '../../_data/kalman/'
tRange=[0,700]
nUnderSamp=5
bExport=False
bExport=True
# bNoise=True
bNoise=False
bMoreNoise=False
bFilterAcc=False
# bFilterAcc=True
nFilt=15

# --- Parameters
base    = '../../_data/NREL5MW'
FstFile = '../../_data/NREL5MW_SimpleSimulations/TurbWindStep_AllDOF.fst'
sPref='_Base'
NoiseRFactor=0
if bNoise:
    sPref+='_Noise'
    NoiseRFactor=1/10
if bMoreNoise:
    sPref+='More'
    NoiseRFactor=1/5
if bFilterAcc:
    sPref+='_FilterAcc'

MeasFile   = FstFile.replace('.fst','.outb')
Case       = os.path.basename(FstFile.replace('.fst',''))
OutputFile = OutDir+'300_7states_'+Case+sPref+'.csv'
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
    sigX['Qaero']  = 8*10**6
    sigX['Qgen']   = 1.5*10**6
    sigX['WS']     = 1.0
    # Measurements - more or less half the std
    sigY=dict()
    sigY['TTacc'] = 0.08  # m/s^2
    sigY['omega'] = 0.05 # rad/s
    sigY['Qgen']  = 1*10**6
    sigY['pitch'] = 2.00

try:
    os.mkdir(OutDir)
except:
    pass

# --------------------------------------------------------------------------------}
# --- Kalman filter estimation 
# --------------------------------------------------------------------------------{
with Timer('Simulation Loop'):
    KF= KalmanFilterTNSim(FstFile, MeasFile, OutputFile, base, bThrustInStates, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX, sigY, bExport)

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
for j,s in enumerate(KF.sY):
    print('SNR {:s} - clean: {:.2f} - meas {:.2f} - est {:.2f}'.format(s, SNR(KF.Y_clean[j,:]), SNR(KF.Y[j,:]), SNR(KF.Y_hat[j,:])))

# --- Compare Measurements
KF.plot_Y()
# --- Compare States
KF.plot_X()
#
KF.plot_summary()
KF.plot_moments()
                                                 
plt.show()
