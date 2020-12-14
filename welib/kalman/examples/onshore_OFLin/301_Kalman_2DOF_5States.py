""" Documentation
 This scripts uses:
  - 2 mechanical DOF:    'u' (tower bending)   'psi' (shaft rotation)
  - 4 measurements:               'TT acc',  'omega_rotor' ,    'Mgen' , 'Pitch' 
  - 5 states:             'u',  'azimuth', 'udot', 'omega_rotor', 'Qaero'
  - 3 inputs:             'T',  'Qgen', 'pitch'
 The estimated states are compared to the simulation at the end
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from kalman_models import KalmanModel
from wtDigiTwin.kalman.TN    import KalmanFilterTN, KalmanFilterTNSim
from wtDigiTwin.kalman.TNLin import KalmanFilterTNLin, KalmanFilterTNLinSim
from wtDigiTwin.tools.fatigue import eq_load
from wtDigiTwin.tools.tictoc import Timer


# --- General Parameters and State specific parameters
sWT='NREL5MW'
StateModel='nt1_nx8' # aug = Ta, Qa, Qg, WS
StateModel='nt1_nx5' # aug = Qa
# bYAMS    = True # nt1_nx8 with YAMS matrices
bYAMS    = False
Qgen_LSS = True


# Options for 7 states
OutDir = '../../_data/kalman/'
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

# --- Parameters
FstFile = '../../_data/{:s}_SimpleSimulations/TurbWindStep_AllDOF.fst'.format(sWT)

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

MeasFile   = FstFile.replace('.fst','.outb')
Case       = os.path.basename(FstFile.replace('.fst',''))
if bYAMS:
    OutputFile = OutDir+Case+'_{:s}_{:s}'.format(sWT, 'YAMS')+sPref+'.csv'
else:
    OutputFile = OutDir+Case+'_{:s}_{:s}'.format(sWT, StateModel)+sPref+'.csv'
StateFile = '../../_data/{}_2DOF_ABCD_mean.dat'.format(sWT)
# StateFile = '../data/turbines/{}_2DOF_ABCD_ws12.dat'.format(sWT)
# StateFile = '../data/turbines/{}_2DOF_ABCD_ws5.dat'.format(sWT)
base      = '../../_data/{}'.format(sWT)
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

try:
    os.mkdir(OutDir)
except:
    pass

# --------------------------------------------------------------------------------}
# --- Kalman filter estimation 
# --------------------------------------------------------------------------------{
with Timer('Simulation Loop'):
    if bYAMS:
        bThrustInStates=True
        KF= KalmanFilterTNSim(FstFile, MeasFile, OutputFile, base, bThrustInStates, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX, sigY, bExport)
    else:
        KM = KalmanModel(StateModel, Qgen_LSS)
        KM.ThrustHack=True
        KF= KalmanFilterTNLinSim(KM, FstFile, MeasFile, OutputFile, base, StateFile, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX, sigY, bExport)
# --------------------------------------------------------------------------------}
# --- PostPro  
# --------------------------------------------------------------------------------{
def Leq(t,y,m=5):
    return eq_load(y, m=m, neq=t[-1]-t[0])[0][0]
def SNR(y):
    return np.mean(y**2)/np.std(y)**2
# --- Leq
try:
    for j in [2,5]:
        print('Leq  ref: {:.2f} - est: {:.2f}'.format(Leq(KF.time,KF.M_ref[j]),Leq(KF.time,KF.M_sim[j])))
    # --- Signal-to-noise ratio
    for j,s in enumerate(KF.sY):
        print('SNR {:s} - clean: {:.2f} - meas {:.2f} - est {:.2f}'.format(s, SNR(KF.Y_clean[j,:]), SNR(KF.Y[j,:]), SNR(KF.Y_hat[j,:])))
except:
    pass
# --- Compare Measurements
KF.plot_Y()
# --- Compare Intermediate values
KF.plot_S()
# --- Compare States
KF.plot_X()
#
KF.plot_moments()
KF.plot_summary()
#                                                  
print(OutputFile)
plt.show()
