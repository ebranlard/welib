import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
# Welib
from welib.essentials import *
from welib.kalman.kalman import BuildSystem_Linear_MechOnly 
from welib.kalman.kalmanfilter import KalmanFilter
from welib.kalman.kalman_model import AugmentedLinModel

from welib.weio.fast_linearization_file import FASTLinearizationFile
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions
from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.yams.windturbine import monopileSetupFromOpenFAST

# --------------------------------------------------------------------------------}
# -- Augmented Linear Model 
# --------------------------------------------------------------------------------{
class KalmanModelMonopile(AugmentedLinModel):

    def __init__(KM, *args, **kwargs):
        """ 
        See AugmentedLinModel for main data (e.g. sX, A)
        Additional data specific to this model:
          - WT
          - zDepth
          - pHD
          - zeta, qdhScale
        """
        AugmentedLinModel.__init__(KM)

        # --- Define names of physical states, augmented states, measurements, and inputs
        KM.sQ  = ['q_s','q_p','qd_s', 'qd_p']
        KM.sQa = ['q_h', 'qd_h']
        KM.sY  = ['TTacc']
        KM.sY += ['q_p']
        KM.sU  = ['w'] # White noise
        KM.sS  = ['M_sb','F_sb', 'eta', 'Fhx']

        KM.setup(*args, **kwargs)

    def setup(KM, fstFilename=None, shapes_sub=None, hydroShape=None, Tp=None,
              zeta =0.12, qdhScale=1,#Tuning
              ):
        # --- Default arguments
        if shapes_sub is None:
            shapes_sub =[0,4]
        
        # --- Method 1: using MTNSB
        WT = FASTmodel2MTNSB(fstFilename, shapes_sub=shapes_sub, shapes_twr=[], shapes_bld=[], bStiffening=True, main_axis='z', fixedShaft=True, algo='OpenFAST').WT
        if WT.pSS is not None:
            #NOTE('Setting Components', compFile)
            #WT.SS_setComponents(compFile)
            #NOTE('Setting Compute Eta')
            #WT.SS_computeEta(dfRef['Time_[s]'])
            if hydroShape is not None:
                WT.HD_setShapeFunction(hydroShape)
        GM_hydro = WT.pHD['GM_hydro']
        print('GM_hydro:\n', GM_hydro)
        for i, phi in enumerate(WT.fnd.PhiU):
            # NOTE: diagonal only?
            WT.MM[i,i]+=GM_hydro[i,i]

        pHD = WT.pHD
        # --- Method 2: using Legacy TNSB
        # pST, pSS, pHD, Sys, WT, ref = monopileSetupFromOpenFAST(fstFilename, shapes_sub=shapes_sub, TMIN=tRange[0], TMAX=tRange[1], nSubSample=nUnderSamp, hydroShape=hydroShape, reconHydro=True)
        # WT.fnd = WT.twr

        KM.WT = WT
        KM.pHD = pHD

        # --- Parameters that are a function of the structure and ocean conditions
        KM.zDepth = WT.fnd.s_span - WT.WtrDpth

        # --- Setup state matrices, problem specific!
        # State matrix A from MCK
        # Empty inputs/outputs B,C,D
        A,B,C,D = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK, nP=len(KM.sQa), nU=len(KM.sU), nY=len(KM.sY), Fp=None)

        # print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> USING A LIN')
        # linFile = '../simulations/MT100/00_EVA/OF_NoHydro.2.lin'
        # linFileH = '../simulations/MT100/00_EVA/OF.3.lin'
        # SX = ['PtfmSurge_[m]', 'PtfmPitch_[rad]', 'd_PtfmSurge_[m/s]', 'd_PtfmPitch_[rad/s]']
        # dfA = FASTLinearizationFile(linFile).toDataFrame()['A']
        # Alin = dfA.loc[SX,SX]
        # dfAH= FASTLinearizationFile(linFileH).toDataFrame()['A']
        # AlinH = dfAH.loc[SX,SX]
        # A[:4,:4] = Alin.values[:4,:4]
        # A[:4,:4] = AlinH.values[:4,:4]
        KM.qdhScale = qdhScale

        Minv = np.linalg.inv(WT.MM)
        IQD   =[KM.iX['qd_s'], KM.iX['qd_p']]
        A[IQD, KM.iX['qd_h']]  = Minv @ (pHD['k_h'][0], pHD['k_h'][1])*KM.qdhScale  # qd_h influence in mech DOF

        B[KM.iX['qd_h'], KM.iU['w']] = 1 # White noise

        C[KM.iY['TTacc'], :] = A[KM.iX['qd_s'],:] # TTacc is assumed to be qdd_s
        D[KM.iY['TTacc'], :] = B[KM.iX['qd_s'],:] # TTacc is assumed to be qdd_s
        C[KM.iY['q_p'], KM.iX['q_p']] = 1


        # --- HYDRO STATE EQ - See Script 180
        if Tp==12.7:
            KM.Sw= 2.3835e-01
        elif Tp==10.0:
            KM.Sw= 2.3835e-01/2
        else:
            raise NotImplementedError()
        KM.omega_p = 2*np.pi/Tp
        KM.zeta = zeta
        print('omega_p^2', KM.omega_p**2, '2 zeta omega_p', 2*KM.zeta*KM.omega_p)
        A[KM.iX['q_h'], KM.iX['qd_h']]  = 1
        A[KM.iX['qd_h'], KM.iX['q_h']]  = -KM.omega_p**2
        A[KM.iX['qd_h'], KM.iX['qd_h']]  = -2*KM.zeta*KM.omega_p

        KM.A, KM.B, KM.C, KM.D = A, B, C, D

        KM.colMap={
                'q_s'    : 'Q_Sg_[m]' ,
                'qd_s'   : 'QD_Sg_[m/s]' ,
                'TTacc ' : 'NcIMUTAxs_[m/s^2]' ,
                'eta'    : 'Wave1Elev_[m]', 
                'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                'Fhx'    : 'HydroFxi_[N]',
                'F_sb'   : '-ReactFXss_[N]',
                'M_sb'   : '-ReactMYss_[N*m]',
                'q_p'    : 'Q_P_[rad]' ,
                'qd_p'   : 'QD_P_[rad/s]'
            }


# --------------------------------------------------------------------------------}
# -- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterMonopile(KalmanFilter):

    def __init__(KF, KM=None, debug=False):
        """

        """
        # --- Initialize Kalman Filter, variables names (e.g. sX) and matrices (Xx=A)
        KalmanFilter.__init__(KF, KM=KM)
        KF.WT = KM.WT
        KF.ColMap = KM.colMap


    def timeLoop(KF):
        # --- Aliases to shorten notations
        KM = KF.KM
        WT = KM.WT

        # --- Initial conditions
        x = KF.initFromClean()
        KF.X_hat.iloc[0,:] = x


        # --- Time loop
        for it in range(0,KF.nt-1):    
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values
            # --- Inputs
            u = KF.U_clean.iloc[it,:].values

            # --- Predictions of next time step based on current time step
            x, KF.P, _ = KF.estimateTimeStep(u, y, x, KF.P)

            # --- Estimate Generalized hydro force and bending moment (calc output)
            q_h     = x[KF.iX['q_h']]  # eta
            qd_h    = x[KF.iX['qd_h']] # eta_dot
            eta     = q_h
            eta_dot = qd_h
            p_hydro = KM.pHD['phi'] * eta_dot # p_h = k_h(z) q_h(t)
            p_hydro[KM.zDepth>0] = 0 # safety, shoudn't be necessray

            x_dot = np.dot(KF.A, x) + np.dot(KF.B, u)
            p_ext      = np.zeros((3,len(KM.zDepth)))
            p_ext[0,:] = p_hydro

            x_q   = np.array([x[0],x[1]])
            xd_q  = np.array([x_dot[0], x_dot[1]])
            xdd_q = np.array([x_dot[2], x_dot[3]])

            ## Top loads
            F_top = np.array((0.,0.,0.))
            M_top = np.array((0.,0.,0.))
            a_ext = np.array((0.,0.,-WT.gravity)) # external acceleration (gravity/earthquake)
            F_sec, M_sec, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext=a_ext, corrections=0, PhiK=WT.fnd.PhiK)

            xdd_q *=0
            F_sec_h, M_sec_h, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext = a_ext, PhiK=WT.fnd.PhiK)

            
            # --- Store extra info
            KF.S_hat.at[it+1, 'M_sb']   = M_sec[1,0]
            KF.S_hat.at[it+1, 'F_sb']   = F_sec[0,0]
            KF.S_hat.at[it+1, 'eta' ]   = eta
            KF.S_hat.at[it+1, 'Fhx']    = F_sec_h[0,0]
            
            # --- Propagation to next time step
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f ' % (it,KF.time[it]))

