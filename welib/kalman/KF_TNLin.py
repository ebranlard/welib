""" 
Kalman filter model for "Tower Nacelle Shaft" (based on yams TNSB)"

Uses Lin file from OpenFAST

"""
import os
import numpy as np
# Welib
from welib.essentials import *
from welib.tools.stats import comparison_stats
from welib.ws_estimator.tabulated import TabulatedWSEstimator
import welib.fast.fastlib as fastlib
import welib.weio as weio
from welib.fast.linmodel import FASTLinModel, FASTLinModelTNSB
# Kalman
from welib.kalman.kalman import *
from welib.kalman.kalmanfilter import KalmanFilter
from .KF_TN import KalmanFilterTN
from .filters import moving_average






# YAMS
from welib.yams.models.TNSB_FAST import FASTmodel2TNSB

# --------------------------------------------------------------------------------}
# --- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterTNLin(KalmanFilterTN):
    def __init__(KF, StateModel='nt1_nx7', WSE=None, debug=False, hacks=None):
        KF.StateModel = StateModel
        if StateModel=='nt1_nx8':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Thrust','Qaero','Qgen','WS'])
            sY  = np.array(['NcIMUAx','dpsi','Qgen','pitch'])
            sU  = np.array(['pitch'])
            sS  = np.array(['WS'])
            KF.bWSInStates     = True
            KF.bThrustInStates = True
        elif StateModel=='nt1_nx7':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Thrust','Qaero','Qgen'])
            sY  = np.array(['NcIMUAx','dpsi','Qgen','pitch'])
            sU  = np.array(['pitch'])
            sS  = np.array(['WS'])
            KF.bWSInStates     = False
            KF.bThrustInStates = True
        elif StateModel=='nt1_nx6':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Thrust','Qaero'])
            sY  = np.array(['NcIMUAx','dpsi','Qgen','pitch'])
            sU  = np.array(['Qgen','pitch'])
            sS  = np.array(['WS'])
            KF.bWSInStates     = False
            KF.bThrustInStates = True
        elif StateModel=='nt1_nx5':
            sQ  = np.array(['q_FA1'  ,'psi'  ,'dq_FA1','dpsi'] )
            sQa = np.array(['Qaero'])
            sY  = np.array(['NcIMUAx','dpsi','Qgen','pitch'])
            sU  = np.array(['Thrust','Qgen','pitch'])
            sS  = np.array(['Thrust','WS'])
            KF.bWSInStates     = False
            KF.bThrustInStates = False
        else:
            raise ValueError('Unknown StateModel: {}'.format(StateModel))

        # --- Parent init
        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)
        # --- Storing wind speed estimator (based on tabulated aerodynamic data)
        KF.wse = WSE
        KF.debug = debug
        # Hacks
        hacks_def = {'thrust':None, 'WSE':None, 'SL_cleanQ':False, 'SL_cleanFtop':False}
        if hacks is None:
            KF.hacks = hacks_def
        else:
            hacks_def.update(hacks)
            KF.hacks = hacks_def
        if KF.hacks['thrust']=='clean':
            WARN('HACKING, using thrust from measurements for DEBUG ONLY!')
        if KF.hacks['WSE']=='clean_inputs':
            WARN('HACKING, using WSE inputs from measurements for DEBUG ONLY!')
        if KF.hacks['SL_cleanQ']:
            WARN('HACKING, using clean Q for section loads.')
        if KF.hacks['SL_cleanFtop']:
            WARN('HACKING, using clean F for section loads.')

    def __repr__(self):
        s = KalmanFilter.__repr__(self)
        s+=' - hacks  : {} \n'.format(self.hacks)
        return s

    def setup_matrices(KF, 
                       fstFile,
                       StateFile, Qgen_LSS=True, ThrustHack=False
                       ):
        KF.Qgen_LSS = Qgen_LSS
        KF.ThrustHack = ThrustHack
        shapes_twr =[0]   # TODO detemine this based on sQ
        # --- Windturbine model
        WT = FASTmodel2TNSB(fstFile, shapes_twr=shapes_twr, shapes_bld=[], 
                            DEBUG=False, bStiffening=True, main_axis='z'
                            ).WT
        KF.WT2 = WT
        KF.WT  = WT
        
        nGear = WT.ED['GBRatio']
        
        # --- ColMap
        if Qgen_LSS:
            KF.colMap={
                'q_FA1'  : 'Q_TFA1_[m]',
                'psi'    : '{Azimuth_[deg]} * np.pi/180',   # [deg] -> [rad]
                'dq_FA1' : ' NcIMUTVxs_[m/s]               ' ,
                'dpsi'   : '{RotSpeed_[rpm]} * 2*np.pi/60', # [rpm] -> [rad/s]
                'NcIMUAx': 'NcIMUTAxs_[m/s^2]',
                'Qgen'   : f'{nGear}'+'*{GenTq_[kN-m]}*1000', # [kNm] -> [Nm]  # NOTE: nGear
                'pitch'  : '{BldPitch1_[deg]} * np.pi/180',   # [deg] -> [rad]
                'Thrust' : 'RtAeroFxh_[N]',
                'Qaero'  : 'RtAeroMxh_[N-m]',
                'WS'     : 'RtVAvgxh_[m/s]',


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
              ' NcIMUAx' : ' NcIMUTAxs_[m/s^2]             ' 
            }

        iX = KF.iX
        iY = KF.iY
        iU = KF.iU

        # --- Mechanical system and turbine data
        WT  = FASTLinModelTNSB(fstFile, StateFile=StateFile, DEBUG=False)
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
        # --- Aliases to shorten notations
        # --- Initial conditions
        x = KF.initFromClean(var='x,y,u')
  

        # --- WSE
        WS_last = KF.S_clean['WS'].iloc[0]
        KF.S_hat.loc[0,'WS']= WS_last

        if 'Thrust' not in KF.sX:
            Thrust_last = KF.S_clean.loc[0,'Thrust']
            KF.S_hat.loc[0,'Thrust']= Thrust_last
        
        KF.X_hat.iloc[0,:]   = x


        WSavg      = np.zeros((50,1))
        WSavg[:]=WS_last

        # --- Time loop
        for it in range(0, KF.nt-1):
            t = KF.time[it]
            # --- "Measurements"
            y = KF.Y.iloc[it,:].values
            # --- Inputs
            u = KF.U_clean.iloc[it,:].values.copy()
            if 'Thrust' not in KF.sX:
                u[0] = Thrust_last # (we don't know the thrust)

            # --- Predictions of next time step based on current time step
            t = KF.time[it+1]
            x, KF.P, _ = KF.estimateTimeStep(u, y, x, KF.P, KF.Q, KF.R)

            # --- Estimate Wind Speed
            if KF.hacks['WSE'] == 'clean_inputs':
                Qaero_hat = KF.X_clean['Qaero'].iloc[it]
                omega = KF.X_clean['dpsi'].iloc[it]
            else:
                if 'WS' in KF.sX:
                    WS_last=x[KF.iX['WS']]

                Qaero_hat = x[KF.iX['Qaero']]
                omega     = x[KF.iX['dpsi']]
            pitch = y[KF.iY['pitch']] * 180 / np.pi # deg
            WS_hat, _ = KF.wse.estimate(Qaero_hat, pitch=pitch, omega=omega, WS0=WS_last, relaxation = 0, WSavg=np.mean(WSavg))
            Qaero_hat = np.max(Qaero_hat,0)
            
            # --- Estimate Thrust
            if KF.hacks['thrust'] == 'clean':
                Thrust = KF.U_clean['Thrust'].iloc[it]
            else:
                Thrust = KF.wse.Thrust(WS_hat, pitch=pitch, omega=omega)

            GF = Thrust
            GF = KF.WT2.GF_lin(Thrust,x,bFull=True)

            # --- Store state (partially done by KF.estimateTimeStep)
            if 'Thrust' in KF.sX:
                x[KF.iX['Thrust']] = GF
            if 'WS' in KF.sX:
                x[KF.iX['WS']] = WS_hat
            x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)
            KF.X_hat.iloc[it+1,:]   = x
            KF.Y_hat.iloc[it+1,:]   = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)
            # --- Store extra info
            # Environment
            KF.S_hat.at[it+1, 'WS']     = WS_hat
            # Loads
            if 'Thrust' not in KF.sX:
                KF.S_hat.loc[it+1, 'Thrust']= GF
                
            # --- Propagation to next time step
            Thrust_last = GF
            WS_last   = WS_hat
            WSavg[1:] = WSavg[0:-1]
            WSavg[0]  = WS_hat

            # --- Print status to screen
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it, KF.time[it], WS_hat, Thrust))

    # --- Methods Common between TN and TNLin
    # moments, export, plot_summary, plot_moments
    # 


# --------------------------------------------------------------------------------}
# --- Wrapper For Simulation 
# --------------------------------------------------------------------------------{
def KalmanFilterTNLinSim(fstFile, MeasFile, OutputFile, aeroMapFile, StateFile, 
                         nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigs=None, bExport=False, colMap=None, debug=False, 
                         StateModel='nt1_nx7', Qgen_LSS=True, ThrustHack=False,
                         operFile=None
                         ):

    # --- Creating a wind speed estimator (reads tabulated aerodynamic data)    
    wse = TabulatedWSEstimator(fstFile=fstFile, aeroMapFile=aeroMapFile, operFile=operFile)
    KF = KalmanFilterTNLin(StateModel=StateModel, WSE=wse, debug=debug)
    KF.setup_matrices(fstFile, StateFile, Qgen_LSS=Qgen_LSS, ThrustHack=ThrustHack)
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


