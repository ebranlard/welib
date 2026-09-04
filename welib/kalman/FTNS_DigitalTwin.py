
import pandas as pd
import os
import numpy as np    
import matplotlib.pyplot as plt
import importlib
# import dill as pickle
import dill as pickle

import welib.weio as weio
# WELIB
from welib.kalman.kalman import EmptyStateMat, EmptyStateDF
from welib.ws_estimator.tabulated_floating import TabulatedWSEstimatorFloating
from welib.tools.clean_exceptions import *
from welib.tools.tictoc import Timer
from welib.yams.windturbine import FASTWindTurbine
from welib.yams.models.simulator import SimulatorFromOF , hydroMatToSysMat
from welib.yams.models.generator_oneRigidBody import generateOneRigidBodyModel
from welib.yams.models.generator import generateModel
from welib.fast.hydrodyn import HydroDyn
from welib.FEM.utils import rigidTransformationTwoPoints, rigidTransformationTwoPoints_Loads

# For YAMS
from welib.yams.models.packman import IMUjacobian
from welib.fast.extract import extractIMUAccFromLinFile, mainLinInputs

# Open OpenFAST lin
from welib.fast.linmodel import DEFAULT_COL_MAP_LIN
from welib.fast.linmodel import FASTLinModelFTNSB

# For both
from welib.fast.tools.lin import subMat, matSimpleStateLabels, matToSIunits, renameList



# Local
from welib.kalman.FTNS_SectionLoadsCalc import YAMSSectionLoadCalculatorOptimized
from welib.kalman.FTNS_KalmanFilter import KalmanFilterFTNSLin


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
                # Instead of using the pklFile data for P,T, we'll use the dataframe
                print('[WARN] DigitTwin: Aero estimator using time series!')
                self.AE.setFromTimeSeries(df)

    def setupStateEstimator(self, **opts):
        """ Prepare a Kalman filter based on the model"""
        print('---------------------------DIGITAL TWIN SETUP STATE EST ----------------------------')
        self.KM = KalmanModel(**opts)
        self.SE = KalmanFilterFTNSLin(self.KM, AE=self.AE)

    def setupVirtualSensing(self, vsType='SL_YAMS', **opts):
        if vsType=='SL_YAMS':
            self.VS = YAMSSectionLoadCalculatorOptimized(fstFile=opts['fstFile'])
        else:
            raise NotImplementedError()

    def setupMeasurementData(self, MeasFile, tRange=None, nUnderSamp=1, bFilterPhi=False ,bFilterAcc=False, bFilterOm=False, nFilt=15, NoiseRFactor=0, tuning=None, colMap=None, sigXDict=None):
        # TODO TODO We need to rething the KF and split the 
        #  - loadMeasurements
        #  - init time stepping
        #  - computation of sigma
        #  - introducing of noise
        #  differently. We should be able to add noise to the measurements in this setup step
        print('------------------- DIGITAL TWIN SETUP MEASUREMENT TIME SERIES ---------------------')

        KF = self.SE
        # --- Loading "Measurements" (Defining "clean" values, estimate sigmas from measurements)
        if colMap is None:
            colMap = self.KM.ColMap
        KF.loadMeasurements(MeasFile, nUnderSamp=nUnderSamp, tRange=tRange, ColMap=colMap)
        # --- Default arguments
        if tuning is None:
            tuning['kSigQaero'] = 1
        # --- Tuning of Sigmas
        # KF.sigY['z']/=10000
        # KF.sigY['NcIMUAz']*=100
        # KF.sigX['z']=1
        # KF.sigY=sigY
        # for k,v in KF.sigX.items():
        #     KF.sigX[k]=v/10
        # for k,v in KF.sigY.items():
        #     KF.sigY[k]=v/10
        if sigXDict is not None:
            for k,v in sigXDict.items():
                if k in KF.sigX:
                    print('[INFO] DigiTwin: Overiding Sigma x',k,v)
                    KF.sigX[k] = v

        if 'Qaero' in KF.sigX:
            KF.sigX['Qaero']*=tuning['kSigQaero']
        # if 'x' in KF.sigX:
        #     KF.sigX['x']*=0.0000001
        #KF.print_sigmas()

        # --- Storage for plot and setting up covariances from sigmas
        KF.prepareTimeStepping()  

        # --- Creating noise measuremnts
        KF.prepareMeasurements(NoiseRFactor=NoiseRFactor, bFilterAcc=bFilterAcc, nFilt=nFilt, bFilterPhi=bFilterPhi, bFilterOm=bFilterOm)

        # --- Set Initial conditions
        x = KF.initFromClean()

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


# --- Kalman filter model
class KalmanModel():
    """
    Prepare the state matrices to be given to a Kalman Filter algorithm

    Open up a linear physical model, convert it to an augmented system

    Main Data:
        KM.sQ : list of states
        KM.sQa: list of augmented states
        KM.sY : list of measurements
        KM.sU : list of inputs
        KM.sS : list of additional storage
        KM.Xx = Xx
        KM.Xu = Xu
        KM.Yx = Yx
        KM.Yu = Yu
    """
    def __init__(KM, modelName=None, fstLin=None, usePickle=True, sQ='', sY='', sU='', sQa='', sS='', qop=None, qdop=None, 
            sFramework='OpenFAST',
            tuning=None,
            nGear=1, # TODO get this from WT
            ):
        KM.StateModel=''

        # --- Main Data
        KM.Xx =  None
        KM.Xu =  None
        KM.Yx =  None
        KM.Yu =  None
        KM.sQ  = None
        KM.sQa = None
        KM.sY  = None
        KM.sU  = None
        KM.sS  = None

        # --- Default arguments
        if tuning is None:
            tuning={}
            tuning['zero_threshold'] = 1e-16         # Watch out for 1/J ~ e-8
            tuning['ksThrust']        = 'Thrust'
            tuning['kThrust']         = 10            # Tuning factor for thrust measured acceleration feedback
            tuning['kThrustA']        = 1             # Tuning factor for thrust state accelerations
            tuning['kIMUz_z']         = 1             # Tuning factor z into IMUz
            tuning['kSigQaero']       = 1             # Tuning of covariance for aero torque
            tuning['fullColumns']    = True          # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< IMPORTANT
        fullColumns=tuning['fullColumns']



        # Col MAP for OpenFAST OutFile "Measurements" used for "clean" values
        sIMU=['NcIMUAx','NcIMUAy','NcIMUAz']
        sIMU2=['NcIMUAx','NcIMUAy','NcIMUAz','NcIMUVx','NcIMUVy','NcIMUVz']
        KM.ColMap={
          ' x      ' : ' PtfmSurge_[m]                   '              ,
          ' y      ' : ' PtfmSway_[m]                   '               ,
          ' z      ' : ' {PtfmHeave_[m]}                '              ,
          ' phi_x  ' : ' {PtfmRoll_[deg]}   * np.pi/180               ' , # SI [deg] -> [rad]
          ' phi_y  ' : ' {PtfmPitch_[deg]}  * np.pi/180                ', # SI [deg] -> [rad]
          ' phi_z  ' : ' {PtfmYaw_[deg]}    * np.pi/180              '  , # SI [deg] -> [rad]
          #' q_FA1  ' : ' TTDspFA_[m]                   '                ,
          ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   '                , # SI [deg] -> [rad]
          ' q_FA1  ' : ' Q_TFA1_[m]                   '                ,
          ' q_SS1  ' : ' Q_TSS1_[m]                   '                ,
          ' dpsi  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 '                , # SI [rpm] -> [rad/s]
          ' dq_FA1 ' : ' QD_TFA1_[m/s]               '                ,
          ' dq_SS1 ' : ' QD_TSS1_[m/s]               '                ,
          ' dx     ' : ' QD_Sg_[m/s]               '              ,
          ' dy     ' : ' QD_Sw_[m/s]                '               ,
          ' dz     ' : ' QD_Hv_[m/s]               '              ,
          ' dphi_x ' : ' QD_R_[rad/s]                             ' ,
          ' dphi_y ' : ' QD_P_[rad/s]                             ',
          ' dphi_z ' : ' QD_Y_[rad/s]                             '  ,
          ' ddpsi  ' : 'QD2_GeAz_[rad/s^2]'  ,
          ' ddq_FA1' : 'QD2_TFA1_[m/s^2]               '                ,
          ' ddq_SS1' : 'QD2_TSS1_[m/s^2]               '                ,
          ' ddx    ' : 'QD2_Sg_[m/s^2]               '              ,
          ' ddy    ' : 'QD2_Sw_[m/s^2]                '               ,
          ' ddz    ' : 'QD2_Hv_[m/s^2]               '              ,
          ' ddphi_x' : 'QD2_R_[rad/s^2]                             ' ,
          ' ddphi_y' : 'QD2_P_[rad/s^2]                             ',
          ' ddphi_z' : 'QD2_Y_[rad/s^2]                             '  ,
          ' Thrust ' : ' RtFldFxh_[N]                 '                ,
          ' Qaero  ' : ' RtFldMxh_[N-m]               '                ,
#           ' Qgen   ' : ' {GenTq_[kN-m]}  *1000         '             , # [kNm] -> [Nm]
          ' Qgen   ' : ' {GenTq_[kN-m]}  *1000 '+'*{}'.format(nGear)             , # [kNm] -> [Nm]
          ' WS     ' : ' RtVAvgxh_[m/s]                '                ,
          ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 '                , # SI [deg]->[rad]
          ' NcIMUAx ' : ' NcIMUTAxs_[m/s^2]             ',
          ' NcIMUAy ' : ' NcIMUTAys_[m/s^2]             ',
          ' NcIMUAz ' : ' NcIMUTAzs_[m/s^2]             ',
          ' NcIMUVx ' : ' NcIMUTVxs_[m/s]             ',
          ' NcIMUVy ' : ' NcIMUTVys_[m/s]             ',
          ' NcIMUVz ' : ' NcIMUTVzs_[m/s]             ',
          # Extrapolations
          }


        ColMapLinFile = DEFAULT_COL_MAP_LIN


        # --------------------------------------------------------------------------------}
        # ---  Linear Physical model
        # --------------------------------------------------------------------------------{

        frameworks = [s.strip() for s in sFramework.split(',')]
        if 'YAMS' in frameworks:
            # --- YAMS
            WT = FASTWindTurbine(fstFilename, twrShapes=[0], algo='OpenFAST')
            sysLI, sim = get_physical_model(WT, modelName, fstLin, qop=qop, qdop=qdop, usePickle=usePickle, noBlin=True)
            sX0= list(sim.WT.DOFname)
            sX = sX0 + ['d'+sx for sx in sX0]
            sXd =['d'+sx for sx in sX]
            A_YAMS = pd.DataFrame(data=sysLI.A, index=sXd, columns=sX)
            #dq  = ((np.max(dfFS[sq]) -np.min(dfFS[sq]))/100).values
            #dqd = ((np.max(dfFS[sqd])-np.min(dfFS[sqd]))/100).values
            #Kacc_fd, Cacc_fd = IMUjacobian(pkg, q0, qd0, p, 'finiteDifferences', dq, dqd)
            uop = sim.uop
            u=dict()
            for key in sim.pkg.info()['su']:
                u[key]= lambda t,q=None,qd=None: 0
            if qop is None:
                qop = sim.qop
            if qdop is None:
                qdop = sim.qdop
            Kacc, Cacc, acc0 = IMUjacobian(sim.pkg, q0=qop, qd0=qdop, p=sim.p, u=u, uop=uop, method='packageJacobians', sDOFs=sX0)
            Kacc2, Cacc2, acc02 = IMUjacobian(sim.pkg, q0=qop, qd0=qdop, p=sim.p, u=u, uop=uop, method='finiteDifferences', sDOFs=sX0, dq=qop*0+0.01, dqd=qdop*0+0.01)
            CIMU_YAMS = pd.concat((Kacc,Cacc),axis=1)
            CIMU_YAMS.index=sIMU

        if 'OpenFAST' in frameworks:
            # --- OpenFAST
            # --- FASTLinModelFTNSB is an instance of FASTLinModel, instance of LinearStateSpace
            # Used to handle one or several lin files
            linmodel = FASTLinModelFTNSB(fstFilename=fstLin, usePickle=usePickle)
            linmodel.rename(verbose=False)
            #print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> lin model')
            #print(linmodel)
            A_OF, B, C, D = linmodel.toDataFrames()
            #linmodel.extract(sX=sX, sU=sU, sY=sY, check=False)

        # Chose between OF or YAMS
        if frameworks[0]=='OpenFAST':
            A=A_OF
        else:
            A=A_YAMS

        # --- 
        A[abs(A)<tuning['zero_threshold']]=0
        B[abs(B)<tuning['zero_threshold']]=0
        C[abs(C)<tuning['zero_threshold']]=0
        D[abs(D)<tuning['zero_threshold']]=0


        # --------------------------------------------------------------------------------}
        # ---  Kalman model (Augmented/modified physical model)
        # --------------------------------------------------------------------------------{
        KM.sQ   = [s.strip() for s in sQ.split(',') if len(s)>0] # Assumed to include derivatives
        KM.sU   = [s.strip() for s in sU.split(',') if len(s)>0]   if len(sU)>0 else  []
        KM.sY   = [s.strip() for s in sY.split(',') if len(s)>0]   if len(sY)>0 else  []
        KM.sS   = [s.strip() for s in sS.split(',') if len(s)>0]   if len(sS)>0 else  []
        KM.sQa  = [s.strip() for s in sQa.split(',') if len(s)>0]  if len(sQa)>0 else []
        KM.sQD  = ['d'+s for s in KM.sQ] # All states derivatives

        # --- Build linear system
        nX = len(KM.sQ)+len(KM.sQa)
        nU = len(KM.sU   )
        nY = len(KM.sY  )
        nq = len(KM.sQ)
        sX0= list(KM.sQ)
        sX = list(KM.sQ)+list(KM.sQa)
        sQd = KM.sQD
        sU = KM.sU
        sY = KM.sY

        Xx, Xu, Yx, Yu = EmptyStateDF(nX,nU,nY,sX,sU,sY)

        # --------------------------------------------------------------------------------}
        # --- Filling state matrix Xx
        # --------------------------------------------------------------------------------{

        # basic A matrix
        sXd = ['d'+s for s in sX]
        for sqd in sQd:
            if sqd not in A.index:
                raise Exception('{} not present in Xx ({})'.format(sqd, A.index))
            for sq in KM.sQ:
                Xx.loc[sqd,sq] = A.loc[sqd,sq]

        # --- Hard Coding
        def setter(MM, sM, srow, scol, value, verbose=True):
            if srow not in MM.index:
                print('[WARN] KalmanModel: Matrix {}: Row {} is not present in matrix {}'.format(sM, srow))
                return
            if scol not in MM.columns:
                print('[WARN] KalmanModel: Matrix {}: Col {} is not present in matrix {}'.format(sM, srow))
                return
            if verbose:
                print('[INFO] KalmanModel: Matrix {}: Replacing [{:10s} x  {:10s}] from {:} to {} '.format(sM, srow, scol, MM.loc[srow, scol], value))
            MM.loc[srow,scol]  = value

        # --- Who influences omega
        if not fullColumns:
            for s in KM.sQ:
                if 'dpsi' in Xx.columns:
                    setter(Xx, 'Xx', 'dpsi' , s, 0) 
                if 'ddpsi' in Xx.columns:
                    setter(Xx, 'Xx', 'ddpsi', s, 0) 
        else:
            if 'psi' in Xx.columns:
                setter(Xx, 'Xx', 'ddpsi', 'psi', 0) 

        
        # --- Useful channels from lin file
        colAugForce = mainLinInputs(hub=2, nac=1, ptfm=2, gen=1, pitch=1)
        colAugForce2 = renameList(colAugForce, ColMapLinFile)
        colAugForce3 = [c for c in colAugForce2 if c in KM.sQa or c in KM.sU or c in KM.sQ]

        # --- Main B matrix
        # Thrust fay faz Qaero may maz
        B = subMat(B, rows=None, cols=colAugForce2, check=True)

        # --- Rotor Inertia
        if 'ddpsi' in A.index:
            J_LSS_YAMS = linmodel.WT.rot.inertia[0,0] 
            J_LSS_OF_Qgen   = -1/B.loc['ddpsi','Qgen']
            J_LSS_OF_Qaero  =  1/B.loc['ddpsi','Qaero']
            J_LSS = J_LSS_OF_Qgen # Selection
            print('[INFO] KalmanModel: Rotor Inertia seleted: {}'.format(J_LSS))

            if 'Qaero' in KM.sQa:
                if fullColumns:
                    Xx.loc[sQd, 'Qaero'] = B.loc[sQd, 'Qaero'] # <<<<<
                setter(Xx, 'Xx', 'ddpsi', 'Qaero', 1/J_LSS)

            if 'Qgen' in KM.sQa:
                if fullColumns:
                    Xx.loc[sQd, 'Qgen'] = B.loc[sQd, 'Qgen'] # <<<<<
                setter(Xx, 'Xx', 'ddpsi', 'Qgen' , -1/J_LSS)
            if 'Qgen' in KM.sU:
                if fullColumns:
                    Xu.loc[sQd, 'Qgen'] = B.loc[sQd, 'Qgen'] # <<<<<
                setter(Xu, 'Xu', 'ddpsi', 'Qgen' ,-1/J_LSS)

        # --- Thrust
        sThrust = tuning['sThrust']
        if 'Thrust' in KM.sQa or 'Thrust' in KM.sU:
            BFHx = B.loc[sQd, 'Thrust'] # Hub x force
            BFNx = B.loc[sQd, 'NacFxN1_[N]'] # Nacelle x force
            BFx_selected = B.loc[sQd, sThrust]*tuning['kThrustA']
            print('[INFO] KalmanModel: Thrust ddq relation: {}'.format(BFx_selected.loc['ddq_FA1']))
            if 'Thrust' in KM.sQa:
                if fullColumns:
                    Xx.loc[sQd, 'Thrust'] = BFx_selected.loc[sQd]
                Xx.loc['ddq_FA1','Thrust'] = BFx_selected.loc['ddq_FA1']
            if 'Thrust' in KM.sU:
                if fullColumns:
                    Xu.loc[sQd, 'Thrust'] = BFx_selected.loc[sQd]
                Xu.loc['ddq_FA1','Thrust'] = BFx_selected.loc['ddq_FA1']

        # --------------------------------------------------------------------------------}
        # --- Filling output matrix Yx
        # --------------------------------------------------------------------------------{
        # States directly measured (States that are in Y directly)
        sYX = [sx for sx in sX if sx in sY] # States that are in Y
        for sxy in sYX:
            sy=sxy
            Yx.loc[sy,sxy] = 1
        # IMU
        _,_,CIMU, DIMU = linmodel.extract(sX=KM.sQ, sU=colAugForce3, sY=sIMU2, verbose=False, check=False, inPlace=False)
        sYIMU = [sy for sy in sIMU2 if sy in sY]
        for sy in sYIMU:
            for sx in KM.sQ:
                Yx.loc[sy,sx] = CIMU.loc[sy,sx]
            for sx in KM.sQa:
                if sx in Yx.columns: # Augmented states
                    Yx.loc[sy,sx] = DIMU.loc[sy,sx]
        if 'Thrust' in Yx.columns:
            if fullColumns:
                Yx.loc[sYIMU,'Thrust'] = DIMU.loc[sYIMU,sThrust]*tuning['kThrust']
            else:
                Yx.loc[sYIMU,'Thrust'] = 0
                Yx.loc['NcIMUAx','Thrust'] = DIMU.loc['NcIMUAx',sThrust]*tuning['kThrust']


        # --------------------------------------------------------------------------------}
        # ---  Yu matrix
        # --------------------------------------------------------------------------------{
        # --- Inputs directly measured (Inputs that are in Y directly)
        sUY = [su for su in sU if su in sY] # Inputs that are in Y
        for su in sUY:
            Yu.loc[su,su] = 1
        # --- IMU
        for sy in sYIMU:
            for sx in DIMU.columns:
                if sx in Yu.columns:
                    Yu.loc[sy,sx] = DIMU.loc[sy,sx]

        if 'Thrust' in Yu.columns:
            if fullColumns:
                Yu.loc[sYIMU,'Thrust'] = DIMU.loc[sYIMU,sThrust]*tuning['kThrust']
            else:
                Yu.loc[sYIMU,'Thrust'] = 0
                Yu.loc['NcIMUAx','Thrust'] = DIMU.loc['NcIMUAx',sThrust]*tuning['kThrust']

        # --- HACK Attempts
        # Heave is a bit too crazy
        if 'z' in Yx.columns and 'NcIMUAz' in Yx.index:
            Yx.loc['NcIMUAz','z'] *= tuning['kIMUz_z']

        #printMat('Xx',Xx, xmin=1e-8)
        #printMat('Yx',Yx, xmin=1e-8)
        #print('A --------------------------------------------------------\n',Xx)
        #print('C --------------------------------------------------------\n',Yx)
        #print('B --------------------------------------------------------\n',Xu)
        #print('D --------------------------------------------------------\n',Yu)
        #print('   --------------------------------------------------------\n')

        KM.Xx = Xx
        KM.Xu = Xu
        KM.Yx = Yx
        KM.Yu = Yu

        hasControl=True
        try:
            import control
        except:
            hasControl=False
        if hasControl:
            O = control.obsv(Xx,Yx)
            try:
                sys = control.StateSpace(Xx,Xu,Yx,Yu)
            except:
                print('[FAIL] State space')
                pass
            try:
                Wc = control.gram(sys, 'c')
            except:
                print('[FAIL] gramian Controlability')
                pass
            try:
                print('[FAIL] gramian Observability')
                Wo  = control.gram(sys, 'o')
            except:
                pass





def get_physical_model(WT, modelName, fstFilename, qop=None, qdop=None, usePickle=True, qopFst=False, noBlin=True, MCKh=None):
    """ 
    Return YAMS physical model.
    Less and less use
    """
    #
    pickleFilename = os.path.splitext(fstFilename)[0]+'_linModelYAMS.pkl'

    if usePickle:
        # If a pickle exist, we load it, and then return
        if os.path.exists(pickleFilename):
            sysLI, sim = pickle.load(open(pickleFilename,'rb'))
            sim.reloadPackage()
            return sysLI, sim
        else:
            print('[FAIL] Pickle file not found:',pickleFilename)
    tMax=0
    # --- Setup Sim
    print('----------------------- SETUP SIMULATION -----------------------------------------')
    sim = SimulatorFromOF(WT, modelName=modelName, packageDir='py')
    if modelName[0]=='B':
        time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
        zRef = -sim.p['z_B0']
    else:
        time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)
        zRef =  sim.p['z_OT']
    su = sim.pkg.info()['su']
    sq = sim.WT.DOFname
    sqd = sim.WT.dDOFname

    # --- uop
    print('----------------------- OPERATING POINT ------------------------------------------')
    # --- Q0
    qop_ = pd.Series(data=np.zeros(len(sq)), index=sq)
    if qop is not None:
        for i,s in enumerate(sq): 
            if s in qop_.index:
                qop_.loc[s] =qop[i]
            else:
                print('[WARN] {} not found in qop'.format(s))
        print('[INFO] Setting qop to:', dict(qop_))
    else:
        if qopFst:
            q0=WT.q0
            for s in sq:
                if s not in q0:
                    raise Exception('DOF {} is not found in fst simulation (available:{})'.format(s,dict(q0)))
                else:
                    qop_[s] = q0[s]
            # sanity check
            for s in q0.keys():
                if s not in qop_.index:
                    print('[WARN] DOF {} present in fst simulation but not used: '.format(s))
            print('[INFO] Using q0 from FST:', dict(qop_))
        else:
            print('[INFO] Using q0 is zero:', dict(qop_))
    # --- QD0
    qdop_ = pd.Series(data=np.zeros(len(sqd)), index=sqd)
    if qdop is not None:
        for i,s in enumerate(sqd): 
            if s in qdop_.index:
                qdop_.loc[s] =qdop[i]
            else:
                print('[WARN] {} not found in qdop'.format(s))
        print('[INFO] Setting qdop to:', dict(qdop_))
    else:
        if qopFst:
            qd0=WT.qd0
            for s in sqd:
                if s not in qd0:
                    raise Exception('DOF {} is not found in fst simulation (available:{})'.format(s,dict(qd0)))
                else:
                    qdop_[s] = qd0[s]
            # sanity check
            for s in qd0.keys():
                if s not in qdop_.index:
                    print('[WARN] DOF {} present in fst simulation but not used: '.format(s))
            print('[INFO] Using qd0 from FST:', dict(qdop_))
        else:
            print('[INFO] Using qd0 is zero:', dict(qdop_))

    uop = sim.uop
    sim.qop  = qop_.values.flatten()
    sim.qdop = qdop_.values.flatten()


    # --- Linear Hydro
    print('----------------------- LINEAR HYDRO  --------------------------------------------')
    q0h_ = pd.Series(data=np.zeros(6), index=['x','y','z','phi_x','phi_y','phi_z'])
    for s in enumerate(q0h_.index): 
        if s in qop_.index:
            q0h_.loc[s] = qop_[s]
        #
    q0h = q0h_.values.flatten()
    hd = HydroDyn(fstFilename)
    if MCKh == 0:
        Mh=np.zeros((6,6))
        Ch=np.zeros((6,6))
        Kh=np.zeros((6,6))
    if MCKh is None:
        if 'hydroO' in modelName:
            MCKFh = hd.linearize_RigidMotion2Loads(q0h, RefPointMotion=(0,0,zRef), RefPointMapping=(0,0,zRef) ) # <<< Good if hydroO model
        else:
            MCKFh = hd.linearize_RigidMotion2Loads(q0h, RefPointMotion=(0,0,zRef), RefPointMapping=(0,0,0) ) # <<< Good if hydro0 model
    #       MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,0), RefPointMapping=(0,0,0) ) # OLD and BAD
        Mh,Ch,Kh,Fh0=MCKFh
    #print('Ch\n',Ch)

    Mh_=hydroMatToSysMat(Mh, su, sq)
    Ch_=hydroMatToSysMat(Ch, su, sq)
    Kh_=hydroMatToSysMat(Kh, su, sq)
    Fh_=hydroMatToSysMat(Fh0, su)
    #print('>>> Ch_\n',Ch_)
    #print('>>> Mh_\n',Mh_)
    #print('>>> Kh_\n',Kh_)
    #print('>>> Fh_\n',Fh_)
    MCKu = Mh_, Ch_, Kh_

    if WT.MAP is not None:
        print('----------------------- LINEAR MOOR ----------------------------------------------')
        print("Mooring stiffness matrix (0,0,zRef={})".format(sim.p['z_OT']))
        print(WT.MAP._K_lin) # TODO might depend on qop

    # --- Simulation
    sysLI = sim.linmodel(MCKextra=None, MCKu=MCKu, noBlin=noBlin)
    print(sysLI)

    #dfNL = sysNL.toDataFrame(self.channels, self.FASTDOFScales, acc=acc, forcing=forcing, sAcc=self.acc_channels)

    if usePickle:
        WT.MAP=None # Can't output C
        sim.unloadPackage() # Can't store imported module
        pickle.dump((sysLI, sim), open(pickleFilename,'wb'))
        print('>>> Pickle file written:',pickleFilename)
        sysLI, sim = pickle.load(open(pickleFilename,'rb'))
        sim.reloadPackage()

    return sysLI, sim




