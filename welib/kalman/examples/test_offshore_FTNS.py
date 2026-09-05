""" 
Perform a Digital twin simulation of the NREL 5-MW Spar
The Wind Speed Estimator is not used in this script.

"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from welib.weio.fast_output_file import FASTOutputFile
from welib.tools.stats import *
from welib.essentials import *
from welib.tools.colors import MathematicaBlue, ManuDarkOrange
COLRS =[MathematicaBlue, ManuDarkOrange]

from welib.kalman.FTNS_DigitalTwin import DigitalTwin

import pytest

scriptDir=os.path.dirname(__file__)

np.set_printoptions(linewidth=300, precision=2)
pd.set_option('display.max_rows', 50, 'display.max_columns', 50,'display.width', 400, 'display.precision',2)

def main(sWT='FTNS',test=False):
    # --------------------------------------------------------------------------------
    # --- Parameters for digital twin simulation
    # --------------------------------------------------------------------------------
    tRange = None
    # ---- Script parameters
    if sWT=='FTNS':
        nGear = 1  # gear ratio
        fstLin          = os.path.join(scriptDir, '_data/FTNS_lin_NoAero_OP2/Main.fst'); labelLin='NewOpNoAero'
        fstSim          = os.path.join(scriptDir, '_data/FTNS_IrregWave_Turb/Main.fst'); 
        WSE_pklFilename = os.path.join(scriptDir, '_data/FTNS_CPCT.pkl')
        tRange = [300,400]
        virtualSensing = True

    elif sWT=='Spar':
        nGear = 97  # gear ratio
        fstLin          = os.path.join(scriptDir,'_data/Spar_lin/Main.fst');
        fstSim          = os.path.join(scriptDir,'_data/Spar_sim/Main.fst');
        WSE_pklFilename = None                                            # TODO provide a pickle file to be able to estimate thrust and wind speed
        virtualSensing = False



    nUnderSamp = 0     # Undersampling of measurement inputs
    NoiseRFactor=0.   # Signal to noise ratio to add to measurements
    bFilterAcc=False  # Filter the accelration measurements FILTER ACC IMPROVES SPECTRA BUT INCREASE REL ERR OF My
    nFilt=15
    Isec = [0,4,7] # Section indices

    # --- Derived Params
    sLabel=''
    if NoiseRFactor>0:
        sLabel+='_Noise{}'.format(NoiseRFactor)
    if bFilterAcc>0:
        sLabel+='_FilterAcc'

    # --- DEFAULT
    labelLin=''
    sFramework='OpenFAST'
    modelName=None
#     modelName = 'F6T1N0S1_fnd_moorO_hydroO';


    # --- Measurement 
    MeasFile = fstSim.replace('.fst','.outb') # Here we use an OpenFAST simulation as "measurements"

    # --- State estimator settings
    # State-space model (Q: states, Qa: augmented states, Y: outputs, U: inputs, S: storage)
    sQ  = 'x,y,z,phi_x,phi_y,phi_z,q_FA1,psi,dx,dy,dz,dphi_x,dphi_y,dphi_z,dq_FA1,dpsi'
    sQa = 'Qaero'
    sY  = 'x,y,phi_x,phi_y,dpsi,NcIMUAx,NcIMUAy,NcIMUAz,Qgen'
    sU  = 'Qgen,pitch,Thrust'
    if sWT=='FTNS':
        sS  = 'WS,'
        sS += ','.join(['TwHt{}MLyt_[kN-m]'.format(i+1) for i in Isec])
    else:
        sS  = '' # e.g. WS

    # --- Tuning option for state estimator
    # DEFAULT
    tuning={}
    tuning['zero_threshold'] = 1e-9          # Watch out for 1/J ~ e-8
    tuning['sThrust']        = 'NacFxN1_[N]'
    tuning['sThrust']        = 'Thrust'
    tuning['kIMUz_z']        = 1             # Tuning factor z into IMUz
    tuning['fullColumns']    = False         # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< IMPORTANT NOTE: breaks psi (debug)
    # Case specific
    tuning['kThrust']     = 0.3 # Tuning factor for thrust measured acceleration feedback
    tuning['kThrustA']    = 1.  # Tuning factor for thrust state accelerations
    tuning['kIMUz_z']     = 3   # Tuning factor z into IMUz
    tuning['fullColumns'] = True # Use full columns of lin matrices or key values of it <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< IMPORTANT
    useDtForCov=True
    Pidentity=False
    Pidentity=True
    useDtForCov=False
#     tuning['kSigQaero']   = 1.0  # Tuning of covariance for aero torque
#     tuning['kSigPsi']     = 1  # Tuning of covariance for psi
#     tuning['kSigQaero']   = 10  # Tuning of covariance for aero torque

    sigXDict={}
    sigQDict={}
#     sigXDict['Qaero'] = 750000.000
    dt_for_sigQ = 0.05
    sigQDict['Qaero'] = 700000.000
    sigQDict['x']        =       0.270*10  #       0.270
    sigQDict['y']       =        0.020*1   #       0.020
    sigQDict['z']       =        0.065*1   #       0.065
    sigQDict['phi_x']   =        0.000     #       0.000
    sigQDict['phi_y']   =        0.002     #       0.002
    sigQDict['phi_z']   =        0.002     #       0.002
    sigQDict['q_FA1']   =        0.007*1   #       0.007
    sigQDict['psi']     =        0.085     #       0.085
    sigQDict['dx']      =        0.007     #       0.007
    sigQDict['dy']      =        0.000     #       0.000
    sigQDict['dz']      =        0.006     #       0.006
    sigQDict['dphi_x']  =        0.000     #       0.000
    sigQDict['dphi_y']  =        0.000     #       0.000
    sigQDict['dphi_z']  =        0.000     #       0.000
    sigQDict['dq_FA1']  =        0.019     #       0.019
    sigQDict['dpsi']    =        0.001     #       0.001
    sigQDict['Qaero']   =   700000.000     #  700000.000





#     sigQDict['Qaero'] = 600000.000
#     14000000.000

    # Combining state Estimator Options into one dictionary
    SEOpts={'modelName':modelName, 'fstLin':fstLin, 'sQ':sQ, 'sU':sU, 'sQa':sQa, 'sS':sS, 'sY':sY, 
            'sFramework':sFramework, 'nGear':nGear, 'tuning':tuning, 'fstFilename':fstSim}


    # --------------------------------------------------------------------------------
    # --- Digital twin simulation
    # --------------------------------------------------------------------------------
    # --- Initialization of the digital twin with aero estimator, state estimator and virtual sensing
    DT = DigitalTwin()
    DT.setupAeroEstimator(fstSim, WSE_pklFilename) #, df=MeasFile)
    DT.setupStateEstimator(**SEOpts)
    DT.setupVirtualSensing(vsType='SL_YAMS', fstFile=fstSim)
    # --- Load measurement data (given time series)
    DT.setupMeasurementData(MeasFile, tRange=tRange, nUnderSamp=nUnderSamp, bFilterAcc=bFilterAcc, nFilt=nFilt, NoiseRFactor=NoiseRFactor)

    DT.setupCovariances(tuning=tuning, useDt=useDtForCov, Pidentity=Pidentity, sigXDict=sigXDict, sigQDict=sigQDict, dt_for_sigQ=dt_for_sigQ, verbose=True)


    # --- Perform digital twin simulation for that measurement timeseries
    KF = DT.timeLoop(virtualSensing=True)

    # --- Export results to file
    dfAll=KF.toDataFrame()
    resOut = fstSim.replace('.fst','_KF_{}_tmax{}.outb'.format(labelLin+sLabel, str(tRange[1])))
    dfAll = KF.saveOutputs(resOut)

    print('-------------------------------------------------------')
    print('fstFile: ', fstSim)
    print('Qdiag  : ', np.diag(KF.Q))
    print('Rdiag  : ', np.diag(KF.R))
    print('Cmat   : ', KF.C.values)

    # --- Simple plots and stats
    statsDict = {}
    figNames=[]
    tRangeStats=None
    labelSim = os.path.dirname(fstSim)
    figX = KF.plot_X(nPlotCols=2, figSize=(12.8,8.2), 
                     printStats=True, tRangeStats=tRangeStats, statsDict=statsDict,
                     title='States - LinFile:{} Sim:{}{}'.format(labelLin,labelSim,sLabel.replace('_',' ')), COLRS=COLRS)
    # fig.subplots_adjust(left=0.12, right=0.98, top=0.955, bottom=0.12, hspace=0.20, wspace=0.20)
    figName = fstSim.replace('.fst','_KF_{}.png'.format(labelLin+sLabel))
    figX.savefig(figName)
    figNames.append(figName)


    #figY = KF.plot_Y()
    #figU = KF.plot_U()
    figS = KF.plot_S(printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
    # KF.plot_P()
    # KF.plot_K()
    # KF.plot_innovation()


    if figS is not None:
        wse=KF.wse
        if wse is not None:
            Ylim1 = [0,20] # WS
            dfY = KF.Y_clean
            dfU = KF.U_clean
            dfX = KF.X_clean

            # Where the data is invalid
            bInv = np.logical_or.reduce((dfY['dpsi'] <   wse.omega[0],   dfY['dpsi']  >wse.omega[-1]))
            bInv = np.logical_or.reduce((dfU['pitch']   <wse.pitch[0],   dfU['pitch'] >wse.pitch[-1], bInv))
            bInv = np.logical_or.reduce((dfX['phi_y']*180/np.pi<wse.phiy[0],   dfX['phi_y']*180/np.pi>wse.phiy[-1] , bInv))
        #     bInv = np.logical_or.reduce((df['WS_ref_[m/s]']   <wse.WS[0],     df['WS_ref_[m/s]']   >wse.WS[-1]   , bInv))
            bInv = np.logical_or.reduce((dfU['Qgen']<100, bInv))
    #         bInv =dfU['Qgen']<100
            print('>> nInvalid:',sum(bInv))
        #     t = df['Time_[s]'].values 
        # 
        #     #with Timer('ValidValues'):
        #     #    bInv2 =  ~wse.validValues(df['WS_ref_[m/s]'], df['Omega_[rad/s]'], df['Pitch_[deg]']  , df['PtfmPitch_[deg]'])
        #     bInv3 = np.logical_or(bInv, bInv2)
        # 
        #     bInv = bInv3
        #     b    = ~bInv3
        # 
        #     # WS plot TODO find which one
            axes=figS.axes
            ax=axes[0]
            ax.fill_between(KF.time, Ylim1[0], Ylim1[1], where=bInv, alpha=0.1, color=(0.5,0.5,0.5))
            ax.set_ylim(Ylim1)


        figName = fstSim.replace('.fst','_KF_{}_WS.png'.format(labelLin+sLabel))
        figS.savefig(figName)
        figNames.append(figName)


    if KF.dfExtra is not None:
        from welib.tools.strings import latexStrip
        df = KF.df 
        dfSL = KF.dfExtra
        # --------------------------------------------------------------------------------}
        # --- PLOT MOMENTS 
        # --------------------------------------------------------------------------------{
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(12.4,12.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        for iiED,iED in enumerate(Isec):
            sT = 'TwHt{}MLyt_[kN-m]'.format(iED+1)
            t1 = df['Time_[s]'].values
            y1 = df[sT].values
            t2 = dfSL['Time_[s]'].values
            y2 = dfSL[sT].values
            stats, sStats =  comparison_stats(t1, y1, t2, y2, stats='sigRatio,eps,R2', method='1-2', latex=True)
            sT = sT.split('_')[0]
            statsDict[sT] = stats
            print(f"{sT:10s} "+latexStrip(sStats))
            ax.plot(t1, y1, 'k-'                        , label='OpenFAST' if iiED==0 else None)
            ax.plot(t2, y2,  '--' , color=fColrs(iiED)  , label='Sim Ht{}'.format(iED+1))
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.legend()

        figName = fstSim.replace('.fst','_KF_{}_Moments.png'.format(labelLin+sLabel))
        fig.savefig(figName)
        figNames.append(figName)
    print('-------------------------------------------------------')
    for figName in figNames:
        print('>>> FIG', figName);

    return statsDict

def test_offshore_FTNS():
    if os.getenv('GITHUB_ACTIONS') == 'true':
        NOTE('test Offshore FTNS not ready yet')
        pytest.skip("Skipping local-only test on GitHub Actions")

    stats = main(sWT='FTNS', test=False)
    np.testing.assert_array_less(stats['x']['eps'],     2.2)
    np.testing.assert_array_less(stats['y']['eps'],     4.1)
    np.testing.assert_array_less(stats['phi_y']['eps'], 6.1)
    np.testing.assert_array_less(stats['Qaero']['eps'], 4.5)
    np.testing.assert_array_less(stats['WS']['eps'],    5.0)
    np.testing.assert_array_less(stats['TwHt1MLyt_[kN-m]']['eps'],    5.2)
    np.testing.assert_array_less(stats['TwHt5MLyt_[kN-m]']['eps'],    5.0)
    np.testing.assert_array_less(stats['TwHt8MLyt_[kN-m]']['eps'],    6.1)

if __name__ == '__main__':
    test_offshore_FTNS()
    #statsMoments = main(sWT='FTNS', test=False)
    #statsMoments = main(sWT='Spar', test=False)

    plt.show()
