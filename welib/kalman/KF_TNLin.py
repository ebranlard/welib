""" 
Kalman filter model for "Tower Nacelle Shaft" (based on yams TNSB)"

Uses Lin file from OpenFAST

"""

import numpy as np
from .kalman import *
from .kalmanfilter import KalmanFilter
from .kalman_model import AugmentedLinModel
from .KF_TN import KalmanFilterTN
from .filters import moving_average
from welib.ws_estimator.tabulated import TabulatedWSEstimator
from welib.fast.linmodel import FASTLinModel, FASTLinModelTNSB
from welib.yams.models.TNSB_FAST import FASTmodel2TNSB
from welib.tools.stats import comparison_stats
import welib.fast.fastlib as fastlib
import welib.weio as weio

# --------------------------------------------------------------------------------}
# -- Augmented Linear Model 
# --------------------------------------------------------------------------------{
# This is the complicated step, setting up the state matrices based on various inputs 
class KalmanModelTNLin(AugmentedLinModel):
    def __init__(self, FstFile, StateFile, StateModel='nt1_nx7', Qgen_LSS=True, ThrustHack=False):
        AugmentedLinModel.__init__(self)

        self.StateModel=StateModel
        self.Qgen_LSS=Qgen_LSS
        self.ThrustHack=ThrustHack

        WT2= FASTmodel2TNSB(FstFile , shapes_twr=[0],shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z').WT
        #WT2.DD      = WT2.DD*3.5 # increased damping to account for aero damping
        self.WT2 = WT2
        
        nGear = WT2.ED['GBRatio']

        if Qgen_LSS:
            self.ColMap={
              ' ut1    ' : ' TTDspFA_[m]                   ' ,
              ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
              ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
              ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
              ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
              ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
              ' Qgen   ' : f'{nGear}'+'*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
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
            self.sQ  = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sQa = np.array(['Thrust','Qaero','Qgen','WS'])
            self.sY  = np.array(['TTacc','omega','Qgen','pitch'])
            self.sU  = np.array(['pitch'])
            self.sS  = np.array(['WS'])
            self.bWSInStates     = True
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx7':
            self.sQ  = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sQa = np.array(['Thrust','Qaero','Qgen'])
            self.sY  = np.array(['TTacc','omega','Qgen','pitch'])
            self.sU  = np.array(['pitch'])
            self.sS  = np.array(['WS'])
            self.bWSInStates     = False
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx6':
            self.sQ  = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sQa = np.array(['Thrust','Qaero'])
            self.sY  = np.array(['TTacc','omega','Qgen','pitch'])
            self.sU  = np.array(['Qgen','pitch'])
            self.sS  = np.array(['WS'])
            self.bWSInStates     = False
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx5':
            self.sQ  = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sQa = np.array(['Qaero'])
            self.sY  = np.array(['TTacc','omega','Qgen','pitch'])
            self.sU  = np.array(['Thrust','Qgen','pitch'])
            self.sS  = np.array(['Thrust','WS'])
            self.bWSInStates     = False
            self.bThrustInStates = False


        
        iX = self.iX
        iY = self.iY
        iU = self.iU

        # --- Mechanical system and turbine data
  

        WT  = FASTLinModelTNSB(FstFile, StateFile=StateFile, DEBUG=False)
        self.WT=WT
        A,B,C,D,M = WT.A, WT.B, WT.C, WT.D, WT.M # To Shorten notations

        # --- Build linear system
        nX = len(self.sQ)+len(self.sQa)
        nU = len(self.sU   )
        nY = len(self.sY  )
        nq = len(self.sQ)
        #
        nGear = WT.nGear
        Mqt      =  1/B.iloc[2,0]
        J_LSSnG  = -1/B.iloc[3,1]  # This is JLSS * nGear (scaled by nGea since it takes Torque at HSS and return influences at LSS). It is not J_HSS !!

        Mqt_ED   = M.iloc[0,0]
        J_LSS_ED = M.iloc[1,1]
        
        Xx, Xu, Yx, Yu = EmptyStateMat(nX, nU, nY)
        # --- Filling extended state matrices
        if self.StateModel=='nt1_nx8' or self.StateModel=='nt1_nx7': # sQa =  ['Thrust','Qaero','Qgen','WS']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values[:,2:]
            Yu[:  ,:   ] = D.values[:,2:]
            #----
            Xx[iX['omega'],iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            Xx[:nq,iX['Thrust']]  = B.values[:,0]
            Yx[:,  iX['Thrust']]  = D.values[:,0]
            Yx[iY['Qgen'],iX['Qgen']] = 1
            # --- Value Hack
            if self.ThrustHack:
                Xx[iX['ut1dot'], iX['Thrust']] =  2.285e-06  # Thrust
#             Xx[iX['omega'],  iX['Qaero']]  =  2.345e-08  # Qa
#             Xx[2,0:4] =[ -6.132e+00,      0,   -5.730e-02,      0]
#             Xx[3,4  ] =  0
#             Xu[3,0  ] =  0
            # --- Consistency
            if self.Qgen_LSS:
                Xx[iX['omega'],  iX['Qgen']]   =-Xx[iX['omega'],iX['Qaero']]
            else:
                Xx[iX['omega'],  iX['Qgen']]   =-Xx[iX['omega'],iX['Qaero']]*nGear
            Yx[0,0:] =Xx[2,0:]  # <<<< Important


        elif self.StateModel=='nt1_nx6': # sQa = ['Qaero','Thrust']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values[:,1:]
            Yu[:  ,:   ] = D.values[:,1:]
            #----
            Xx[iX['omega'],  iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            Xx[:nq,iX['Thrust']] = B.values[:,0]
            Yx[:,  iX['Thrust']] = D.values[:,0]
            #  Value Hack
#             Xx[2,0:4] =[ -6.132e+00,      0,   -5.730e-02,      0]
            if self.ThrustHack:
                Xx[2,iX['Thrust']] =  2.285e-06  # Thrust
#             Xx[3,iX['Qaero' ]] =  2.345e-08  # Torque
#             Xx[3,4  ] =  0
#             Xu[2,0  ] =  0
#             Yu[0,0  ] =  0
            # Consistency
            if self.Qgen_LSS:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]
            else:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]*nGear
            Yx[0,0:6] =Xx[2,0:6]  # <<<< Important

        elif self.StateModel=='nt1_nx5':  # sQa = ['Qaero']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values
            Yu[:  ,:   ] = D.values
            #----
            Xx[iX['omega'],  iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            #  Value Hack
#             Xx[2,0:4] =[ -6.132e+00,      0,   -5.730e-02,      0]
            if self.ThrustHack:
                Xu[2,0  ] =  2.285e-06  # Thrust
#             Xx[3,4]   =  2.345e-08  # Torque
            # Consistency
            if self.Qgen_LSS:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]
            else:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]*nGear
            Yx[0,0:4] = Xx[2,0:4]
            Yu[0,0]   = Xu[2,0]

        self.A = Xx
        self.B = Xu
        self.C = Yx
        self.D = Yu
        
# --------------------------------------------------------------------------------}
# -- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterTNLin(KalmanFilterTN):
    def __init__(KF, KM, WSE=None, debug=False):
        """

        """
        KalmanFilterTN.__init__(KF, KM, WSE=WSE)
        KF.WT2 = KM.WT2
        KF.WT  = KM.WT2




    # --- Methods Common between TN and TNLin
    # loadMeasurements, prepareMeasurements


    def timeLoop(KF):
        # --- Initial conditions
        x = KF.initFromClean()
        P = KF.P        

        # --- WSE
        WS_last  = KF.S_clean.loc[0,'WS']
        KF.S_hat.loc[0,'WS']= WS_last

        if 'Thrust' not in KF.sX:
            Thrust_last = KF.S_clean.loc[0,'Thrust']
            KF.S_hat.loc[0,'Thrust']= Thrust_last
        
        KF.X_hat.iloc[0,:]   = x


        WSavg      = np.zeros((50,1))
        WSavg[:]=WS_last

        for it in range(0,KF.nt-1):    
            t = it*KF.dt
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values

            # --- KF predictions
            u=KF.U_clean.iloc[it,:].values.copy()
            if 'Thrust' not in KF.sX:
                u[0] = Thrust_last # (we don't know the thrust)
            x,P,_ = KF.estimateTimeStep(u,y,x,P,KF.Q,KF.R)

            # --- Estimate thrust and WS - Non generic code
            if 'WS' in KF.sX:
                WS_last=x[KF.iX['WS']]
            pitch     = y[KF.iY['pitch']]*180/np.pi # deg
            Qaero_hat = x[KF.iX['Qaero']]
            omega     = x[KF.iX['omega']]
            WS_hat, _ = KF.wse.estimate(Qaero_hat, pitch, omega, WS_last, relaxation = 0, WSavg=np.mean(WSavg))
            Qaero_hat = np.max(Qaero_hat,0)
            Thrust = KF.wse.Thrust(WS_hat, pitch, omega)

            GF = Thrust
            GF = KF.WT2.GF_lin(Thrust,x,bFull=True)

            # --- Store
            if 'Thrust' in KF.sX:
                x[KF.iX['Thrust']] = GF
            else:
                KF.S_hat.loc[it+1, 'Thrust']= GF
            if 'WS' in KF.sX:
                x[KF.iX['WS']] = WS_hat
            KF.S_hat.loc[it+1, 'WS'    ]= WS_hat
            x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)
            KF.X_hat.iloc[it+1,:]   = x
            KF.Y_hat.iloc[it+1,:]   = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)
            # --- Propagation to next time step
            Thrust_last = GF
            WS_last     = WS_hat
            WSavg[1:] = WSavg[0:-1]
            WSavg[0]  = WS_hat

            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it],WS_hat,Thrust))
        KF.P = P

    # --- Methods Common between TN and TNLin
    # moments, export, plot_summary, plot_moments
    # 


# --------------------------------------------------------------------------------}
# --- Wrapper For Simulation 
# --------------------------------------------------------------------------------{
def KalmanFilterTNLinSim(KM, FstFile, MeasFile, OutputFile, aeroMapFile, StateFile, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX=None, sigY=None, sigQ=None, bExport=False, ColMap=None, debug=False):

    # ---
    # --- Creating a wind speed estimator (reads tabulated aerodynamic data)    
    wse = TabulatedWSEstimator(fstFile=FstFile, aeroMapFile=aeroMapFile)
    # ---
    KF = KalmanFilterTNLin(KM, WSE=wse)
    if debug:
        print(KF.wse)
        print(KF.WT)
        print(KF)
    # --- Loading "Measurements"
    # Defining "clean" values 
    # Estimate sigmas from measurements
    KF.loadMeasurements(MeasFile, nUnderSamp=nUnderSamp, tRange=tRange, ColMap=ColMap)
    KF.sigX=sigX
    KF.sigY=sigY
    KF.sigQ=sigQ
    if debug:
        KF.print_sigmas()

    # --- Process and measurement covariances
    KF.setupCovariances(useDt=False, Pidentity=True)
    # --- Storage for plot
    KF.prepareTimeStepping()
    # --- Creating noise measuremnts
    KF.prepareMeasurements(NoiseRFactor=NoiseRFactor, bFilterAcc=bFilterAcc, nFilt=nFilt)

    # --- Time loop
    if debug:
        print(OutputFile)
    KF.timeLoop()
    KF.moments()

    if bExport:
        KF.export(OutputFile)
    return KF


