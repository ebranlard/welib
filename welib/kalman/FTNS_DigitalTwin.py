
import pandas as pd
import os
import numpy as np    
import matplotlib.pyplot as plt
import importlib

import welib.weio as weio
# WELIB
from welib.essentials import *
from welib.ws_estimator.tabulated_floating import TabulatedWSEstimatorFloating
from welib.kalman.KF_FTNS import KalmanFilterFTNSLin
from welib.kalman.FTNS_SectionLoadsCalc import YAMSSectionLoadCalculatorOptimized

# --- Kalman filter model
class DigitalTwin():
    """ 
    Main Data:
     - SE: state estimator (KalmanFilter)
     - AE: aerodynamic estimator 
     - VS: Virtual sensing algorithm
     - MD: Measurement data

    """
    def __init__(self):
        # --- Main Data
        self.SE = None # State estimator
        self.AE = None # Aerodynamic estimator
        self.VS = None # Virtual sensing algorithm
        self.MD = None # Measurement data

    def setupAeroEstimator(self, fstFile, pickleFile, df=None):
        """ NOTE: fst file is mostly for rotor diameter and airdensity"""
        if pickleFile is None:
            self.AE = None
            print('[WARN] DigiTwin: No aerodynamic estimator')
        else:
            self.AE = TabulatedWSEstimatorFloating(fstFile=fstFile, pickleFile=pickleFile)
            if df is not None:
                if isinstance(df, str):
                    df = weio.read(df).toDataFrame()
                # Instead of using the pklFile data for P,T, we'll use the dataframe
                NOTE('DigitTwin: Aero estimator using time series!')
                self.AE.setFromTimeSeries(df)

    def setupStateEstimator(self, **opts):
        """ Prepare a Kalman filter based on the model"""
        print('---------------------------DIGITAL TWIN SETUP STATE EST ----------------------------')
        self.SE = KalmanFilterFTNSLin(AE=self.AE, **opts)

    def setupVirtualSensing(self, vsType='SL_YAMS', **opts):
        if vsType=='SL_YAMS':
            self.VS = YAMSSectionLoadCalculatorOptimized(fstFile=opts['fstFile'])
        else:
            raise NotImplementedError()

    def setupMeasurementData(self, MeasFile, tRange=None, nUnderSamp=1, bFilterPhi=False ,bFilterAcc=False, bFilterOm=False, nFilt=15, NoiseRFactor=0, colMap=None):
        # TODO TODO We need to rething the KF and split the 
        #  - loadMeasurements
        #  - init time stepping
        #  - computation of sigma
        #  - introducing of noise
        #  differently. We should be able to add noise to the measurements in this setup step
        print('------------------- DIGITAL TWIN SETUP MEASUREMENT TIME SERIES ---------------------')
        INFO('FTNS_DigitalTwin: setupMeasurementData')

        KF = self.SE
        # --- Loading "Measurements" (Defining "clean" values, estimate sigmas from measurements)
        if colMap is None:
            colMap = self.SE.colMap
        KF.loadMeasurements(MeasFile, nUnderSamp=nUnderSamp, tRange=tRange, colMap=colMap)

        # --- Storage for plot
        KF.prepareTimeStepping()  

        # --- Set Initial conditions
        x = KF.initFromClean()

        # --- Creating noise measuremnts
        KF.prepareMeasurements(NoiseRFactor=NoiseRFactor, bFilterAcc=bFilterAcc, nFilt=nFilt, bFilterPhi=bFilterPhi, bFilterOm=bFilterOm)




    def timeLoop(self, virtualSensing=True):
        print('----------------------------- DIGITAL TWIN TIME LOOP -------------------------------')
        KF = self.SE
        # --- Time integrations
        with Timer('Kalman filter time loop'):
            KF.timeLoop()
        # --- virtualSensing 
        if virtualSensing:
            with Timer('Extrapolation'):
                KF.dfExtra = self.virtualSensing()
        else:
            KF.dfExtra = None

        return KF


    def virtualSensing(self):
        """ """
        print('--------------------------- DIGITAL TWIN VIRTUAL SENSING ---------------------------')
        #DOFNames = []
        KF = self.SE
        X_clean = KF.X_clean
        S_clean = KF.S_clean
        X       = KF.X_hat
        XD      = KF.XD_hat

        dfIn = self.VS.emptyInputDF(nt=len(self.SE.time), inputFrame='R_xs')
        dfIn['Time'] = KF.time

        MAPQ = {'Sg':'x', 'Sw':'y', 'Hv':'z' ,'R':'phi_x', 'P':'phi_y', 'Y':'phi_z', 'TFA1':'q_FA1', 'TSS1':'q_SS1', 'Yaw':'q_yaw'}

        for sDOF,sShort in MAPQ.items():
            sq   = sShort
            if sq in X.keys():
                if sq in ['phi_z']:
                    #print('[ OK ] Using zero for', sq)
                    dfIn['Q_'+sDOF]   = 0
                    dfIn['QD_'+sDOF]  = 0
                    dfIn['QD2_'+sDOF] = 0
#                 elif sq in ['y','phi_x']:
#                     print('[ OK ] Using vel/acc from signal for', sq)
#                     x_smooth = moving_average(df[sq].values, n=15)
#                     vel = ddt(x_smooth, df['Time'].values)
#                     acc = ddt(vel, df['Time'].values)
#                     dfIn['Q_'+sDOF]   = df[sq]
#                     dfIn['QD_'+sDOF]  = vel
#                     dfIn['QD2_'+sDOF] = 0
                else:
                    #print('[ OK ] Using Hat for', sq)
                    dfIn['Q_'+sDOF]   = X[sq]
                    dfIn['QD_'+sDOF]  = X['d'+sq]
                    dfIn['QD2_'+sDOF] = XD['dd'+sq]  - np.mean(XD['dd'+sq])
# #                 dfIn['Q_'+sDOF]   = X[sq]
# #                 dfIn['QD_'+sDOF]  = X['d'+sq]
#                 dfIn['Q_'+sDOF]   = X_clean[sq]
#                 dfIn['QD_'+sDOF]  = X_clean['d'+sq]
#                 dfIn['QD2_'+sDOF] = XD['dd'+sq]
            else:
                print('[INFO] DigiTwin: Virtual sensing: state no present  ', sShort)
        if 'Qaero' in X.keys():
            dfIn['Madd_R_xs'] = X['Qaero']
        if 'Thrust' in X.keys():
            print('[INFO] DigiTwin: Virtual Sensing: using Thrust from X')
            dfIn['Fadd_R_xs'] = X['Thrust']
        elif 'Thrust' in KF.U_hat.keys():
            print('[INFO] DigiTwin: Virtual Sensing: using Thrust from U_hat')
            dfIn['Fadd_R_xs'] = KF.U_hat['Thrust']
        elif 'Thrust' in KF.S_hat.keys():
            print('[INFO] DigiTwin: Virtual Sensing: using Thrust from S_hat')
            dfIn['Fadd_R_xs'] = KF.S_hat['Thrust']

        dfSL = self.VS.fromDF(dfIn, useTopLoadsFromDF=False)

        # Store in KF
        for s in KF.sS:
            if s not in dfSL.columns:
                print('[INFO] DigiTwin: Virtual Sensing: storage column not computed: ',s)
            else:
                KF.S_hat[s] = dfSL[s]

        return dfSL
#         WT=KF.WT2
#         z_test = fastlib.ED_TwrGag(WT.ED) - WT.ED['TowerBsHt']
#         EI     = np.interp(z_test, WT.Twr.s_span, WT.Twr.EI[0,:])
#         kappa  = np.interp(z_test, WT.Twr.s_span, WT.Twr.PhiK[0][0,:])
#         qx    = KF.X_hat[KF.iX['ut1']]
#         KF.M_sim = [qx*EI[i]*kappa[i]/1000 for i in range(len(z_test))]                 # in [kNm]
#         KF.M_ref=[]
#         for i in range(len(z_test)):
#             try:
#                 val=KF.df['TwHt{:d}MLyt_[kN-m]'.format(i+1)].values
#             except:
#                 try:
#                     val=KF.df['TwHt{:d}MLyt'.format(i+1)].values
#                 except:
#                    val=KF.time*0
#             KF.M_ref.append(val)


