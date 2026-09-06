import os
import numpy as np

from welib.kalman.kalman import *
from .kalmanfilter import KalmanFilter
from .kalman_model import AugmentedLinModel

from welib.essentials import *
from welib.kalman.filters import moving_average

# --- External dependencies!
import welib.fast.fastlib as fastlib
import welib.weio as weio

from welib.yams.windturbine import FASTWindTurbine
from welib.yams.models.simulator import SimulatorFromOF , hydroMatToSysMat
from welib.fast.hydrodyn import HydroDyn

# For YAMS
from welib.yams.models.packman import IMUjacobian
from welib.fast.extract import mainLinInputs

# Open OpenFAST lin
from welib.fast.linmodel import DEFAULT_COL_MAP_LIN
from welib.fast.linmodel import FASTLinModelFTNSB

# For both
from welib.fast.tools.lin import subMat, matSimpleStateLabels, matToSIunits, renameList

# Local
from welib.kalman.kalman import EmptyStateMat, EmptyStateDF

# 
# #          'WS':'Wind1VelX', 'pitch':'BldPitch1','TTacc':'NcIMUTAxs'}
# #          'Thrust':'RotThrust','Qaero':'RtAeroMxh','Qgen':'GenTq',
# # NOTE: RotThrust contain gravity and inertia
# DEFAULT_COL_MAP={
# #   ' ut1    ' : ' TTDspFA_[m]                   ' ,
# #   ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
# #   ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
# #   ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
# #   ' thrust ' : ' rtaerofxh_[n]                 ' ,
# #   ' qaero  ' : ' rtaeromxh_[n-m]               ' ,
# #   ' thrust ' : ' rtfldfxh_[n]                 ' ,
# #   ' qaero  ' : ' rtfldmxh_[n-m]               ' ,
# # #   ' qgen   ' : ' {gentq_[kn-m]}  *1000         ' , # [knm] -> [nm]
# #   ' qgen   ' : ' 97*{gentq_[kn-m]}  *1000         ' , # [knm] -> [nm]
# #   ' ws     ' : ' rtvavgxh_[m/s]                ' ,
# #   ' pitch  ' : ' {bldpitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
# #   ' ttacc  ' : ' ncimutaxs_[m/s^2]             ' 
#   ' x      ' : ' ptfmsurge_[m]                   '              ,
#   ' y      ' : ' ptfmsway_[m]                   '               ,
#   ' z      ' : ' {ptfmheave_[m]}                '              ,
#   ' phi_x  ' : ' {ptfmroll_[deg]}   * np.pi/180               ' , # si [deg] -> [rad]
#   ' phi_y  ' : ' {ptfmpitch_[deg]}  * np.pi/180                ', # si [deg] -> [rad]
#   ' phi_z  ' : ' {ptfmyaw_[deg]}    * np.pi/180              '  , # si [deg] -> [rad]
#   #' q_fa1  ' : ' ttdspfa_[m]                   '                ,
#   ' psi    ' : ' {azimuth_[deg]} * np.pi/180   '                , # si [deg] -> [rad]
#   ' q_fa1  ' : ' q_tfa1_[m]                   '                ,
#   ' q_ss1  ' : ' q_tss1_[m]                   '                ,
#   ' dpsi  ' : ' {rotspeed_[rpm]} * 2*np.pi/60 '                , # si [rpm] -> [rad/s]
#   ' dq_fa1 ' : ' qd_tfa1_[m/s]               '                ,
#   ' dq_ss1 ' : ' qd_tss1_[m/s]               '                ,
#   ' dx     ' : ' qd_sg_[m/s]               '              ,
#   ' dy     ' : ' qd_sw_[m/s]                '               ,
#   ' dz     ' : ' qd_hv_[m/s]               '              ,
#   ' dphi_x ' : ' qd_r_[rad/s]                             ' ,
#   ' dphi_y ' : ' qd_p_[rad/s]                             ',
#   ' dphi_z ' : ' qd_y_[rad/s]                             '  ,
#   ' ddpsi  ' : 'qd2_geaz_[rad/s^2]'  ,
#   ' ddq_fa1' : 'qd2_tfa1_[m/s^2]               '                ,
#   ' ddq_ss1' : 'qd2_tss1_[m/s^2]               '                ,
#   ' ddx    ' : 'qd2_sg_[m/s^2]               '              ,
#   ' ddy    ' : 'qd2_sw_[m/s^2]                '               ,
#   ' ddz    ' : 'qd2_hv_[m/s^2]               '              ,
#   ' ddphi_x' : 'qd2_r_[rad/s^2]                             ' ,
#   ' ddphi_y' : 'qd2_p_[rad/s^2]                             ',
#   ' ddphi_z' : 'qd2_y_[rad/s^2]                             '  ,
#   ' thrust ' : ' rtfldfxh_[n]                 '                ,
#   ' qaero  ' : ' rtfldmxh_[n-m]               '                ,
# #           ' qgen   ' : ' {gentq_[kn-m]}  *1000         '             , # [knm] -> [nm]
#   ' qgen   ' : ' {gentq_[kn-m]}  *1000 '+'*{}'.format(1)             , # [knm] -> [nm] # <<<<<<< todo todo todo ngear
#   ' power  ' : ' {genpwr_[kw]}  *1000 '            , # [knm] -> [nm] # <<<<<<< todo todo todo ngear
#   ' ws     ' : ' rtvavgxh_[m/s]                '                ,
#   ' pitch  ' : ' {bldpitch1_[deg]} * np.pi/180 '                , # si [deg]->[rad]
#   ' ncimuax ' : ' ncimutaxs_[m/s^2]             ',
#   ' ncimuay ' : ' ncimutays_[m/s^2]             ',
#   ' ncimuaz ' : ' ncimutazs_[m/s^2]             ',
#   ' ncimuvx ' : ' ncimutvxs_[m/s]             ',
#   ' ncimuvy ' : ' ncimutvys_[m/s]             ',
#   ' ncimuvz ' : ' ncimutvzs_[m/s]             ',
# }

# --------------------------------------------------------------------------------}
# -- Augmented Linear Model 
# --------------------------------------------------------------------------------{
class KalmanModelFTNS(AugmentedLinModel):

    def __init__(KM, modelName=None, fstLin=None, usePickle=True, fstFilename=None,
                 sQ='', sY='', sU='', sQa='', sS='', qop=None, qdop=None, 
            sFramework='OpenFAST',
            tuning=None,
            nGear=1, # TODO get this from WT
            ):
        AugmentedLinModel.__init__(KM)


        KM.StateModel=''

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
        KM.colMap={
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


        colMapLinFile = DEFAULT_COL_MAP_LIN


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
        colAugForce2 = renameList(colAugForce, colMapLinFile)
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

        KM.A = Xx
        KM.B = Xu
        KM.C = Yx
        KM.D = Yu

        hasControl=True
        try:
            import control
        except:
            hasControl=False
        if hasControl:
            O = control.obsv(A,C)
            try:
                sys = control.StateSpace(A, B, C, D)
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
    import dill as pickle
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








# --------------------------------------------------------------------------------}
# -- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing
class KalmanFilterFTNSLin(KalmanFilter):
    def __init__(KF, KM=None, AE=None, debug=False):
        """

        """
        # --- Initialize Kalman Filter, variables names (e.g. sX) and matrices (Xx=A)
        KalmanFilter.__init__(KF, KM=KM)

        # --- Creating a wind speed estimator (reads tabulated aerodynamic data)
        if AE:
            # --- Wind Speed estimator
            if 'phi_y' not in KF.iX:
                print('[WARN] WSE: phi_y not in X, will assume phi_y=0')
            if 'pitch' not in KF.iU:
                print('[WARN] WSE: pitch not in U, will assume pitch=0')
            if not ('Qaero' in KF.iX or 'Qaero' in KF.iU):
                print('[WARN] WSE: not running WSE becasue Qaero is not in X or U')
                pickleFile=None
                KF.wse=None
            else:
                KF.wse=AE
        else:
            KF.wse = None

    def prepareMeasurements(KF, NoiseRFactor=0, bFilterAcc=False, bFilterOm=False, nFilt=15, bFilterPhi=False):
        # --- Creating noise measuremnts
        KF.setYFromClean(R=KF.R_c, NoiseRFactor=NoiseRFactor)
        if bFilterAcc:
            if 'NcIMUAx' in KF.Y:
                KF.Y['NcIMUAx'] = moving_average(KF.Y['NcIMUAx'],n=nFilt) 
            if 'NcIMUAy' in KF.Y:
                KF.Y['NcIMUAy'] = moving_average(KF.Y['NcIMUAy'],n=nFilt) 
            if 'NcIMUAz' in KF.Y:
                KF.Y['NcIMUAz'] = moving_average(KF.Y['NcIMUAz'],n=nFilt) 
        if bFilterPhi:
            if 'phi_y' in KF.Y:
                print('>>>>> FILTERING PHI_y')
                KF.Y['phi_y'] = moving_average(KF.Y['phi_y'],n=nFilt) 
        if bFilterOm:
            if 'dpsi' in KF.Y:
                print('>>>>> FILTERING OMEGA')
                KF.Y['dpsi'] = moving_average(KF.Y['dpsi'],n=nFilt) 

    def timeLoop(KF):
        # --- Initial conditions
        x = KF.initFromClean()
        P = KF.P
        KF.U_hat.iloc[0,:] = KF.U_clean.iloc[0,:]

        # --- WSE
        if KF.wse:
            WS_last     = KF.S_clean['WS'].values[0].copy()
            KF.S_hat.loc[0, 'WS'] = WS_last
            WSavg      = np.zeros((50,1))
            WSavg[:]=WS_last


        if 'Thrust' in KF.sU:
            Thrust_last = KF.U_clean['Thrust'][0]
            iThrust = list(KF.sU).index('Thrust')

        for it in range(0,KF.nt-1):    
            t = it*KF.dt
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values

            # --- KF predictions
            u=KF.U_clean.iloc[it,:].values.copy()
            if 'Thrust' in KF.sU:
                # We use previous estimated thrust as input.
                u[iThrust] = Thrust_last  
            x,P,_ = KF.estimateTimeStep(u,y,x,P,KF.Q,KF.R)

            # --- Estimate thrust and WS - Non generic code
            if KF.wse:
                if 'WS' in KF.iX:
                    WS_last=x[KF.iX['WS']]
                if 'dpsi' in KF.iX:
                    omega     = x[KF.iX['dpsi']] # in rad/s for WSE
                if 'phi_y' in KF.iX:
                    phiy     = x[KF.iX['phi_y']] * 180/np.pi # in deg for WSE
                elif 'phi_y' in KF.iU:
                    phiy     = u[KF.iU['phi_y']] * 180/np.pi # in deg for WSE
                else:
                    raise Exception('Cannot run WSE if phi_y not in X') # Relax this later
                if 'pitch' in KF.iU:
                    pitch     = u[KF.iU['pitch']]*180/np.pi # in deg for WSE
                if 'Qaero' in KF.iX:
                    Qaero_hat = x[KF.iX['Qaero']]
                elif 'Qaero' in KF.iU:
                    Qaero_hat = u[KF.iU['Qaero']]
                else:
                    raise Exception('Cannot run WSE if Qaero not in X or U')

                #def estimate(self, Qa, omega, pitch,  phiy , WS0, relaxation=0, method='crossing', deltaWSMax=1, verbose=False, debug=False, t=0, WSref=np.nan): 
                WS_hat, _ = KF.wse.estimate(Qaero_hat, omega=omega, pitch=pitch, phiy=phiy, WS0=WS_last, relaxation=0.5, method='oper-crossing', t=t)
                Qaero_hat = np.max(Qaero_hat,0)
                Thrust = KF.wse.Thrust(WS_hat, omega=omega, pitch=pitch, phiy=phiy)
                GF = Thrust
                # GF = KF.WT2.GF_lin(Thrust,x,bFull=True) # TODO TODO
            else:
                WS_hat=0
                GF     = 0
                omega  = 0
                phiy   = 0
                omega  = 0
                pitch  = 0
                Thrust = 0

            # --- Store
            # TODO TODO WHY U IS NOT STORED?
            if 'Thrust' in KF.iX:
                x[KF.iX['Thrust']] = GF
            elif 'Thrust' in KF.iU:
                pass
                #u[KF.iU['Thrust']] = GF
            elif 'Thrust' in KF.iS:
                KF.S_hat.loc[it+1, 'Thrust']= GF

            if 'WS' in KF.iX:
                x[KF.iX['WS']] = WS_hat
            elif 'WS' in KF.iS:
                KF.S_hat.loc[it+1, 'WS'    ]= WS_hat
            if 'phi_y' in KF.iS:
                KF.S_hat.loc[it+1, 'phi_y'    ]= phiy*np.pi/180
            if 'Qaero' in KF.iS:
                KF.S_hat.loc[it+1, 'Qaero'    ]= Qaero_hat
            if 'dpsi' in KF.iS:
                KF.S_hat.loc[it+1, 'dpsi'    ]= omega
            if 'pitch' in KF.iS:
                KF.S_hat.loc[it+1, 'pitch'    ]= pitch*np.pi/180
            if 'WS0' in KF.iS:
                KF.S_hat.loc[it+1, 'WS0'    ]= WS_last

            if 'psi' in KF.iX:
                x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)

            KF.U_hat.iloc[it+1,:]   = u
            KF.X_hat.iloc[it+1,:]   = x
            KF.Y_hat.iloc[it+1,:]   = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)
            KF.XD_hat.iloc[it+1,:]  = np.dot(KF.Xx,x) + np.dot(KF.Xu,u) # Accelerations
            # --- Propagation to next time step
            if not np.isnan(GF):
                Thrust_last = GF
            WS_last     = WS_hat
            #WSavg[1:] = WSavg[0:-1]
            #WSavg[0]  = WS_hat
# 
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it],WS_hat,Thrust))
        KF.P = P

    # --------------------------------------------------------------------------------}
    # --- Extrapolation (calculation of moments  
    # --------------------------------------------------------------------------------{
    def virtualSensing(KF):
        """ """
        pass


    def moments(KF):
        WT=KF.WT2
        z_test = fastlib.ED_TwrGag(WT.ED) - WT.ED['TowerBsHt']
        EI     = np.interp(z_test, WT.Twr.s_span, WT.Twr.EI[0,:])
        kappa  = np.interp(z_test, WT.Twr.s_span, WT.Twr.PhiK[0][0,:])
        qx    = KF.X_hat[KF.iX['ut1']]
        KF.M_sim = [qx*EI[i]*kappa[i]/1000 for i in range(len(z_test))]                 # in [kNm]
        KF.M_ref=[]
        for i in range(len(z_test)):
            try:
                val=KF.df['TwHt{:d}MLyt_[kN-m]'.format(i+1)].values
            except:
                try:
                    val=KF.df['TwHt{:d}MLyt'.format(i+1)].values
                except:
                   val=KF.time*0
            KF.M_ref.append(val)
        return KF.M_sim, KF.M_ref

    def export(KF,OutputFile):
        M=np.column_stack([KF.time]+[KF.X_clean[j,:] for j,_ in enumerate(KF.sX)])
        M=np.column_stack([M]+[KF.X_hat  [j,:] for j,_ in enumerate(KF.sX)])
        M=np.column_stack([M]+[KF.Y      [j,:] for j,_ in enumerate(KF.sY)])
        M=np.column_stack([M]+[KF.Y_hat  [j,:] for j,_ in enumerate(KF.sY)])
        if len(KF.sS)>0:
           M=np.column_stack([M]+[KF.S_clean[j,:] for j,_ in enumerate(KF.sS)])
           M=np.column_stack([M]+[KF.S_hat  [j,:] for j,_ in enumerate(KF.sS)])
        M=np.column_stack([M]+KF.M_ref)
        M=np.column_stack([M]+KF.M_sim)
        header='time'+','
        header+=','.join(['x_'+s+'_ref' for s in KF.sX])+','
        header+=','.join(['x_'+s+'_est' for s in KF.sX])+','
        header+=','.join(['y_'+s+'_ref' for s in KF.sY])+','
        header+=','.join(['y_'+s+'_est' for s in KF.sY])+','
        if len(KF.sS)>0:
            header+=','.join([s+'_ref' for s in KF.sS])+','
            header+=','.join([s+'_est' for s in KF.sS])+','
        header+=','.join(['My_ref{:d}'.format(j) for j,_ in enumerate(KF.M_ref)])+','
        header+=','.join(['My_est{:d}'.format(j) for j,_ in enumerate(KF.M_sim)])
        np.savetxt(OutputFile,M,delimiter=',',header=header)

    def plot_summary(KF):
        import matplotlib
        import matplotlib.pyplot as plt
        cmap = matplotlib.cm.get_cmap('viridis')
        COLRS = [(cmap(v)[0],cmap(v)[1],cmap(v)[2]) for v in np.linspace(0,1,3+1)]

        def spec_plot(ax,t,ref,sim):
            try:
                from pybra.spectral import fft_wrap
            except:
                return
            f1,S1,Info = fft_wrap(t,ref,output_type = 'PSD',averaging = 'Welch', nExp=10, detrend=True)
            f2,S2,Info = fft_wrap(t,sim,output_type = 'PSD',averaging = 'Welch', nExp=10, detrend=True)
            ax.plot(f1,S1,'-' , color=COLRS[0],label='Reference')
            ax.plot(f2,S2,'--', color=COLRS[1],label='simulation')
            ax.set_xlim([0,4])
            ax.set_xlabel('Frequency [Hz]')
            ax.set_yscale('log')
            
        def mean_rel_err(t1,y1,t2,y2):
            if len(y1)!=len(y2):
                y2=np.interp(t1,t2,y2)
            # Method 1 relative to mean
            ref_val = np.mean(y1)
            meanrelerr0=np.mean(np.abs(y1-y2)/ref_val)*100 
            print('Mean rel error {:7.2f} %'.format( meanrelerr0))
            # Method 2 scaling signals
            Min=min(np.min(y1), np.min(y2))
            Max=max(np.max(y1), np.max(y2))
            y1=(y1-Min)/(Max-Min)+0.001
            y2=(y2-Min)/(Max-Min)+0.001
            meanrelerr=np.mean(np.abs(y1-y2)/np.abs(y1))*100 
            print('Mean rel error {:7.2f} %'.format( meanrelerr))
            return meanrelerr,meanrelerr0

        def time_plot(ax,t,ref,sim):
            t=t[1:]
            ref=ref[0:-1]
            sim=sim[1:]

            eps=mean_rel_err(t,ref,t,sim)[1]
            sig_ref=np.std(ref)
            sig_sim=np.std(sim)
            ax.plot(t,ref,'-' , color=COLRS[0])
            ax.plot(t,sim,'--', color=COLRS[1])
            Ylim=ax.get_ylim()
            Xlim=ax.get_xlim()
            ax.text(Xlim[0],Ylim[0]+(Ylim[1]-Ylim[0])*0.8,r'$\epsilon=$'+r'{:.1f}%'.format(eps)+r' - $\sigma_\mathrm{est}/\sigma_\mathrm{ref} = $'+r'{:.3f}'.format(sig_sim/sig_ref), fontsize=11 )

        # Aliases to shorten notations
        iX, iY, iS = KF.iX, KF.iY, KF.iS
        X_clean, X_hat = KF.X_clean, KF.X_hat
        S_clean, S_hat = KF.S_clean, KF.S_hat
        time = KF.time

        ##
        fig=plt.figure()
        # fig.set_size_inches(13.8,4.8,forward=True) # default is (6.4,4.8)
        fig.set_size_inches(13.8,8.8,forward=True) # default is (6.4,4.8)
        ax=fig.add_subplot(6,2,1)
        time_plot(ax,time,X_clean[iX['Qaero'],:]/ 1000, X_hat[iX['Qaero'],:]/ 1000)
        ax.set_ylabel('Aerodynamic Torque [kNm]')

        ax=fig.add_subplot(6,2,2)
        spec_plot(ax,time,X_clean[iX['Qaero'],:]/ 1000, X_hat[iX['Qaero'],:]/ 1000)
        # ax.set_ylabel('Power Spectral Density (Welch Avg.)') 


        ax=fig.add_subplot(6,2,3)
        try:
            time_plot(ax,time,X_clean[iX['WS'],:], X_hat[iX['WS'],:])
        except:
            time_plot(ax,time,S_clean[iS['WS'],:], S_hat[iS['WS'],:])
        ax.set_ylabel('WS [m/s]')

        ax=fig.add_subplot(6,2,4)
        try:
            spec_plot(ax,time,X_clean[iX['WS'],:], X_hat[iX['WS'],:])
        except:
            spec_plot(ax,time,S_clean[iS['WS'],:], S_hat[iS['WS'],:])

        ax=fig.add_subplot(6,2,5)
        time_plot(ax,time,X_clean[iX['omega'],:], X_hat[iX['omega'],:])
        ax.set_ylabel('Omega [RPM]')

        ax=fig.add_subplot(6,2,6)
        spec_plot(ax,time,X_clean[iX['omega'],:], X_hat[iX['omega'],:])

        ax=fig.add_subplot(6,2,7)
        try:
            time_plot(ax,time,X_clean[iX['Thrust'],:]/1000, X_hat[iX['Thrust'],:]/1000)
        except:
            time_plot(ax,time,S_clean[iS['Thrust'],:]/1000, S_hat[iS['Thrust'],:]/1000)
        ax.set_ylabel('Thrust [kN]')

        ax=fig.add_subplot(6,2,8)
        try:
            spec_plot(ax,time,X_clean[iX['Thrust'],:]/1000, X_hat[iX['Thrust'],:]/1000)
        except:
            spec_plot(ax,time,S_clean[iS['Thrust'],:]/1000, S_hat[iS['Thrust'],:]/1000)

        ax=fig.add_subplot(6,2,9)
        time_plot(ax,time,X_clean[iX['ut1'],:], X_hat[iX['ut1'],:])
        ax.set_ylabel('TT position [m]')
        ax=fig.add_subplot(6,2,10)
        spec_plot(ax,time,X_clean[iX['ut1'],:], X_hat[iX['ut1'],:])

        #                
#         z_test = list(fastlib.ED_TwrGag(KF.WT.ED) - KF.WT.ED['TowerBsHt'])
#         try:
#             for i,z in enumerate(z_test):
#                 if np.mean(np.abs(KF.M_ref[i] ))>1:
#                     ax=fig.add_subplot(6,2,11)
#                     time_plot(ax,time,KF.M_ref[i], KF.M_sim[i])
#                     ax.set_ylabel('My [kNm] - z={:.1f}'.format(z))
#                     ax=fig.add_subplot(6,2,12)
#                     spec_plot(ax,time,KF.M_ref[i], KF.M_sim[i])
#                     break
#         except:
#             pass
        try:
            ax=fig.add_subplot(6,2,11)
            time_plot(ax,time,KF.M_ref[2], KF.M_sim[2])
            ax.set_ylabel('My [kNm]')
            ax=fig.add_subplot(6,2,12)
            spec_plot(ax,time,KF.M_ref[2], KF.M_sim[2])
        except:
            pass
#
        #                                         
    def plot_moments(KF,fig=None,scaleByMean=False):
        import matplotlib
        import matplotlib.pyplot as plt

        z_test = list(fastlib.ED_TwrGag(KF.WT.ED) - KF.WT.ED['TowerBsHt'])
        print('z test:',z_test)
        n=len(z_test)
#         z_test.reverse()
        # --- Compare measurements
        cmap = matplotlib.cm.get_cmap('viridis')
        COLRS = [(cmap(v)[0],cmap(v)[1],cmap(v)[2]) for v in np.linspace(0,1,n+1)]
        if fig is None:
            fig=plt.figure()
        fig.set_size_inches(6.4,15.0,forward=True) # default is (6.4,4.8)
        for i,z in enumerate(z_test):
            ax = fig.add_subplot(n,1,i+1)
            M_sim =KF.M_sim[i]
            if scaleByMean:
                M_sim+=-np.mean(KF.M_sim[i])+np.mean(KF.M_ref[i])
            
            ax.plot (KF.time, KF.M_ref[i], 'k-', color='k',       label='Reference' , lw=1)
            ax.plot (KF.time,    M_sim   , '--', color=COLRS[i],label='Estimation', lw=0.8)
            ax.set_ylabel('My z={:.1f}'.format(z))
            ax.tick_params(direction='in')
#             if ii<2:
            if i<n-1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('Time [s]')
                ax.legend()
#             # plt.ylim(0.05*10**8,0.8*10**8)
        ax.set_title('KalmanLoads')

