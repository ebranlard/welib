""" 
Kalman filter model for "Tower Nacelle Shaft" (based on yams TNSB)"

Uses Lin file from OpenFAST

"""

import numpy as np
from .kalman import *
from .kalmanfilter import KalmanFilter
from .KF_TN import KalmanFilterTN
from .filters import moving_average
from welib.ws_estimator.tabulated import TabulatedWSEstimator
from welib.fast.linmodel import FASTLinModel, FASTLinModelTNSB
from welib.yams.models.TNSB_FAST import FASTmodel2TNSB
from welib.tools.stats import comparison_stats
import welib.fast.fastlib as fastlib
import welib.weio as weio

# --------------------------------------------------------------------------------}
# -- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that change from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterTNLin(KalmanFilterTN):
    def __init__(KF, StateModel='nt1_nx7', WSE=None, debug=False):
        KF.StateModel = StateModel
        if StateModel=='nt1_nx8':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Thrust','Qaero','Qgen','WS'])
            sY  = np.array(['TTacc','dpsi','Qgen','pitch'])
            sU  = np.array(['pitch'])
            sS  = np.array(['WS'])
            KF.bWSInStates     = True
            KF.bThrustInStates = True
        elif StateModel=='nt1_nx7':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Thrust','Qaero','Qgen'])
            sY  = np.array(['TTacc','dpsi','Qgen','pitch'])
            sU  = np.array(['pitch'])
            sS  = np.array(['WS'])
            KF.bWSInStates     = False
            KF.bThrustInStates = True
        elif StateModel=='nt1_nx6':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Thrust','Qaero'])
            sY  = np.array(['TTacc','dpsi','Qgen','pitch'])
            sU  = np.array(['Qgen','pitch'])
            sS  = np.array(['WS'])
            KF.bWSInStates     = False
            KF.bThrustInStates = True
        elif StateModel=='nt1_nx5':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Qaero'])
            sY  = np.array(['TTacc','dpsi','Qgen','pitch'])
            sU  = np.array(['Thrust','Qgen','pitch'])
            sS  = np.array(['Thrust','WS'])
            KF.bWSInStates     = False
            KF.bThrustInStates = False
        else:
            raise ValueError('Unknown StateModel: {}'.format(StateModel))

        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)
        KF.wse = WSE
        KF.debug = debug

    def setup_matrices(KF, FstFile, StateFile, Qgen_LSS=True, ThrustHack=False):
        KF.Qgen_LSS = Qgen_LSS
        KF.ThrustHack = ThrustHack

        WT2= FASTmodel2TNSB(FstFile , shapes_twr=[0],shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z').WT
        #WT2.DD      = WT2.DD*3.5 # increased damping to account for aero damping
        KF.WT2 = WT2
        KF.WT  = WT2
        
        nGear = WT2.ED['GBRatio']

        if Qgen_LSS:
            KF.colMap={
              ' q_FA1    ' : ' TTDspFA_[m]                   ' ,
              ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
              ' dq_FA1 ' : ' NcIMUTVxs_[m/s]               ' ,
              ' dpsi   ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
              ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
              ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
              ' Qgen   ' : f'{nGear}'+'*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
              ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
              ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
              ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
            }
        else:
            KF.colMap={
              ' q_FA1    ' : ' TTDspFA_[m]                   ' ,
              ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
              ' dq_FA1 ' : ' NcIMUTVxs_[m/s]               ' ,
              ' dpsi   ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
              ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
              ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
              ' Qgen   ' : ' {GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
              ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
              ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
              ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
            }

        iX = KF.iX
        iY = KF.iY
        iU = KF.iU

        # --- Mechanical system and turbine data
        WT  = FASTLinModelTNSB(FstFile, StateFile=StateFile, DEBUG=False)
        KF.linWT = WT
        A,B,C,D,M = WT.A, WT.B, WT.C, WT.D, WT.M # To Shorten notations

        # --- Build linear system
        nX = KF.nX
        nU = KF.nU
        nY = KF.nY
        nq = KF.nX0
        #
        nGear = WT.nGear
        Mqt      =  1/B.iloc[2,0]
        J_LSSnG  = -1/B.iloc[3,1]  # This is JLSS * nGear (scaled by nGea since it takes Torque at HSS and return influences at LSS). It is not J_HSS !!

        Mqt_ED   = M.iloc[0,0]
        J_LSS_ED = M.iloc[1,1]
        
        Xx, Xu, Yx, Yu = EmptyStateMat(nX, nU, nY)
        # --- Filling extended state matrices
        if KF.StateModel=='nt1_nx8' or KF.StateModel=='nt1_nx7': # sQa =  ['Thrust','Qaero','Qgen','WS']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values[:,2:]
            Yu[:  ,:   ] = D.values[:,2:]
            #----
            Xx[iX['dpsi'],iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            Xx[:nq,iX['Thrust']]  = B.values[:,0]
            Yx[:,  iX['Thrust']]  = D.values[:,0]
            Yx[iY['Qgen'],iX['Qgen']] = 1
            # --- Value Hack
            if KF.ThrustHack:
                Xx[iX['dq_FA1'], iX['Thrust']] =  2.285e-06  # Thrust
            # --- Consistency
            if KF.Qgen_LSS:
                Xx[iX['dpsi'],  iX['Qgen']]   =-Xx[iX['dpsi'],iX['Qaero']]
            else:
                Xx[iX['dpsi'],  iX['Qgen']]   =-Xx[iX['dpsi'],iX['Qaero']]*nGear
            Yx[0,0:] =Xx[2,0:]  # <<<< Important


        elif KF.StateModel=='nt1_nx6': # sQa = ['Qaero','Thrust']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values[:,1:]
            Yu[:  ,:   ] = D.values[:,1:]
            #----
            Xx[iX['dpsi'],  iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            Xx[:nq,iX['Thrust']] = B.values[:,0]
            Yx[:,  iX['Thrust']] = D.values[:,0]
            if KF.ThrustHack:
                Xx[2,iX['Thrust']] =  2.285e-06  # Thrust
            # Consistency
            if KF.Qgen_LSS:
                Xu[iX['dpsi'],iU['Qgen']]  =-Xx[iX['dpsi'],iX['Qaero']]
            else:
                Xu[iX['dpsi'],iU['Qgen']]  =-Xx[iX['dpsi'],iX['Qaero']]*nGear
            Yx[0,0:6] =Xx[2,0:6]  # <<<< Important

        elif KF.StateModel=='nt1_nx5':  # sQa = ['Qaero']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values
            Yu[:  ,:   ] = D.values
            #----
            Xx[iX['dpsi'],  iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            if KF.ThrustHack:
                Xu[2,0  ] =  2.285e-06  # Thrust
            # Consistency
            if KF.Qgen_LSS:
                Xu[iX['dpsi'],iU['Qgen']]  =-Xx[iX['dpsi'],iX['Qaero']]
            else:
                Xu[iX['dpsi'],iU['Qgen']]  =-Xx[iX['dpsi'],iX['Qaero']]*nGear
            Yx[0,0:4] = Xx[2,0:4]
            Yu[0,0]   = Xu[2,0]

        KF.setMat(Xx, Xu, Yx, Yu)




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
            omega     = x[KF.iX['dpsi']]
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
def KalmanFilterTNLinSim(FstFile, MeasFile, OutputFile, aeroMapFile, StateFile, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigs=None, bExport=False, colMap=None, debug=False, StateModel='nt1_nx7', Qgen_LSS=True, ThrustHack=False):

    # --- Creating a wind speed estimator (reads tabulated aerodynamic data)    
    wse = TabulatedWSEstimator(fstFile=FstFile, aeroMapFile=aeroMapFile)
    KF = KalmanFilterTNLin(StateModel=StateModel, WSE=wse, debug=debug)
    KF.setup_matrices(FstFile, StateFile, Qgen_LSS=Qgen_LSS, ThrustHack=ThrustHack)
    if debug:
        print(KF.wse)
        print(KF.WT)
        print(KF)
    # --- Loading "Measurements"
    # Defining "clean" values 
    # Estimate sigmas from measurements
    if colMap is None:
        colMap = KF.colMap
    KF.loadMeasurements(MeasFile, nUnderSamp=nUnderSamp, tRange=tRange, colMap=colMap)
    # --- Process and measurement covariances
    KF.setupCovariances(useDt=False, Pidentity=True, sigs=sigs, verbose=debug)
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


