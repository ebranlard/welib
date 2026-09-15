"""
Kalman filter model for "Monopile Tower Nacelle Shaft" (based on yams MTNSB)

"""
import os
import numpy as np
import pandas as pd
import argparse
import sys
import matplotlib.pyplot as plt

# Welib
from welib.essentials import *
from welib.weio.fast_linearization_file import FASTLinearizationFile
from welib.fast.FASTLin import FASTLin
from welib.ws_estimator.tabulated import TabulatedWSEstimator
import welib.fast.fastlib as fastlib
import welib.weio as weio
# Kalman
from welib.kalman.kalman import *
from welib.kalman.kalmanfilter import KalmanFilter

# YAMS
from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions


# --------------------------------------------------------------------------------}
# --- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterMTNS(KalmanFilter):

    def __init__(KF, WSE=None, debug=False, hacks=None):
        """

        """    
        # --- Initialize Kalman Filter, variables names (e.g. sX) and matrices (Xx=A)
        sQ  = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi'] # Mechanical states
        sQa = ['q_h', 'dq_h', 'Qaero']                                         # Augmented states
        sU  = ['Qgen', 'pitch', 'Thrust', 'Fx_i', 'My_i', 'w']
        sY  = ['PtfmIMUAx', 'PtfmIncly', 'dpsi',  'NcIMUAx', 'Qgen']
        sS  = ['My_sb', 'Fx_sb', 'eta', 'Fx_h', 'WS', 'Thrust']
        # --- Parent init
        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)
        # --- Storing wind speed estimator (based on tabulated aerodynamic data)
        KF.wse = WSE # wind speed estimator
        KF.debug = debug
        # Hacks
        hacks_def = {'thrust':None, 'WSE':None, 'SL_cleanQ':False, 'SL_cleanFtop':False, 'SL_cleanEtaDot':False, 'SL_cleanP':False}
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
        if KF.hacks['SL_cleanEtaDot']:
            WARN('HACKING, using clean eta dot for section loads.')
        if KF.hacks['SL_cleanP']:
            WARN('HACKING, using clean p hydro dot for section loads.')



    def __repr__(self):
        s = KalmanFilter.__repr__(self)
        s+=' - hacks  : {} \n'.format(self.hacks)
        return s

    def setup_matrices(KF, 
                       fstFile, dfTime=None,
                       hydroShapeFile=None, compFile=None, Tp=None, zeta =0.12, qdhScale=1, # Hydro params
                       method='OpenFAST',
                       linFile=None, # For method=='OpenFAST
                       fullColumns=False,
                       ):
        """ Build WT model (sea state, hydro) and state matrices A,B,C,D """
        # --- Default arguments
        shapes_sub =[0,4] # TODO detemine this based on sQ
        shapes_twr =[0]   # TODO detemine this based on sQ
        # --- Windturbine model
        WT = FASTmodel2MTNSB(fstFile, shapes_sub=shapes_sub, shapes_twr=shapes_twr, shapes_bld=[],
                             DEBUG=False, bStiffening=True, main_axis='z', fixedShaft=False,
                             algo='OpenFAST').WT
        KF.WT = WT

        nGear = WT.ED['GBRatio']
        
        # --- ColMap
        # Col MAP for OpenFAST OutFile "Measurements" used for "clean" values
        KF.colMap={
                'x'      : 'Q_Sg_[m]' ,
                'phi_y'  : 'Q_P_[rad]' ,
                'q_FA1'  : 'Q_TFA1_[m]',
                'psi'    : '{Azimuth_[deg]} * np.pi/180', # SI [deg] -> [rad]
                'dx'     : 'QD_Sg_[m/s]',
                'dphi_y' : 'QD_P_[rad/s]',
                'dq_FA1' : 'QD_TFA1_[m/s]',
                'dpsi'   : '{RotSpeed_[rpm]} * 2*np.pi/60', # SI [rpm] -> [rad/s]
                'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                'ddx'    : 'QD2_Sg_[m/s^2]',
                'ddphi_y': 'QD2_P_[rad/s^2]',
                'ddq_FA1': 'QD2_TFA1_[m/s^2]',
                'ddpsi'  : 'QD2_GeAz_[rad/s^2]',
                'PtfmIMUAx': '{QD2_Sg_[m/s^2]}',
                'PtfmIncly': '{Q_P_[rad]}',
                'NcIMUAx': 'NcIMUTAxs_[m/s^2]',
                #'NcIMUAy': 'NcIMUTAys_[m/s^2]',
                #'NcIMUAz': 'NcIMUTAzs_[m/s^2]',
                'Qgen'   : f'{nGear}'+'*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]  # NOTE: nGear
                'pitch'  : '{BldPitch1_[deg]} * np.pi/180', # SI [deg]->[rad]
                'Thrust' : 'RtAeroFxh_[N]',
                'Qaero'  : 'RtAeroMxh_[N-m]',
                'WS'     : 'RtVAvgxh_[m/s]',
                'Fx_h'   : '{HydroFxi_[N]}',   # We use a trick  
                'Fx_sb'  : '-ReactFXss_[N]',
                'My_sb'  : '-ReactMYss_[N*m]',

        }
        if WT.fnd is not None:
            KF.colMap.update({
                'Fx_i'   : 'IntfFXss_[N]',   
                'My_i'   : 'IntfMYss_[N*m]',
                }
            )
        else:
            KF.colMap.update({
                'Fx_i'   : '{TwrBsFxt_[kN]}*1000',   
                'My_i'   : '{TwrBsMyt_[kN-m]}*1000',
                }
            )



        # --- Configure Sea State components and wave elevation
        if WT.pSS is not None:
            print('[INFO] Setting Components', compFile)
            WT.SS_setComponents(compFile)
            print('[INFO] Setting Compute Eta')
            WT.SS_computeEta(dfTime)
            if hydroShapeFile is not None:
                WT.HD_setShapeFunction(hydroShapeFile)
        # Ensure WT.MM contains the hydrodynamic mass if not already added
        pHD = WT.pHD
        if not getattr(WT, '_GM_hydro_added', False) and 'GM_hydro' in pHD:
            for i in range(min(pHD['GM_hydro'].shape[0], WT.MM.shape[0])):
                WT.MM[i, i] += pHD['GM_hydro'][i, i]
            WT._GM_hydro_added = True

        # --- Store turbine and hydro data in Kalman filter object
        KF.WT  = WT
        KF.pHD = pHD
        KF.zDepth = WT.fnd.s_span - WT.WtrDpth

        # --- Setup state matrices, problem specific!
        nX, nU, nY = len(KF.sX), len(KF.sU), len(KF.sY)
        # Empty inputs/outputs B,C,D
        MM_sub = WT.MM[:2,:2] # We only keep surge and pitch
        M_inv     = np.linalg.inv(WT.MM)
        M_inv_sub = np.linalg.inv(MM_sub)


        # --- Matrices from OpenFAST lin file
        if linFile is None or not os.path.exists(linFile):
            raise FileNotFoundError('An OpenFAST .lin file is required for now')
        A_OF = np.zeros((nX, nX))
        B_OF = np.zeros((nX, nU))
        C_OF = np.zeros((nY, nX))
        D_OF = np.zeros((nY, nU))
        OF_lin = KF._set_openfast_submatrix(A_OF, B_OF, C_OF, linFile)

        # --- Matrices from YAMS (partial
        A_YS = np.zeros((nX, nX))
        B_YS = np.zeros((nX, nU))
        C_YS = np.zeros((nY, nX))
        D_YS = np.zeros((nY, nU))
        As, Bs, Cs, Ds = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK)
        # --- Structural DOFs
        names = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi']
        for row, row_name in enumerate(names):
            for col, col_name in enumerate(names):
                A_YS[KF.iX[row_name], KF.iX[col_name]] = As[row, col]

        if method=='YAMS':
            A, B, C, D = A_YS, B_YS, C_YS, D_YS
        elif method=='OpenFAST':
            A, B, C, D = A_OF, B_OF, C_OF, D_OF

        # --------------------------------------------------------------------------------}
        # ---  Code common to OpenFAST and YAMS
        # --------------------------------------------------------------------------------{
        # After applying the "generic" A, we introduce the augmented states equations
        # and some additional tweaks

        # --- Rotor Inertia / Shaft equation
        A[KF.iX['psi'], KF.iX['dpsi']] = 1
        J_LSS_YAMS      = WT.rot.inertia[0,0]
        J_LSS_OF_Qgen   = -1/OF_lin['B'].loc['d_psi_rot_[rad/s]','Qgen_[Nm]']
        J_LSS = J_LSS_YAMS
        print('[INFO] KalmanModel: Rotor Inertia seleted: {:.1f} (YAMS: {:.1f} OF: {:.1f}'.format(J_LSS, J_LSS_YAMS, J_LSS_OF_Qgen))
        if 'Qaero' in KF.sXa:
            A[KF.iX['dpsi'], KF.iX['Qaero']] = 1 / J_LSS
        if 'Qgen' in KF.sU:
            B[KF.iX['dpsi'], KF.iU['Qgen']] = -1 / J_LSS


        # --- Thrust influence (physical derivation)
        try:
            r_TN = WT.r_TN_inT # if hasattr(WT, 'r_TN_inT') else [0,0,198.386]
            r_NS = WT.r_NS_inN # if hasattr(WT, 'r_NS_inN') else [0,0,4.143]
        except:
            raise NotImplementedError()
        h_hub = r_TN[2] + r_NS[2]

        sQd_OF = ['d_PtfmSurge_[m/s]', 'd_PtfmHeave_[m/s]', 'd_PtfmPitch_[rad/s]', 'd_qt1FA_[m/s]', 'd_psi_rot_[rad/s]']
        BFHx = OF_lin['B'].loc[sQd_OF, 'HubFxN1_[N]'] # Hub x force
        BFNx = OF_lin['B'].loc[sQd_OF, 'NacFxN1_[N]'] # Nacelle x force
        #BFx_selected = OF_lin['B'].loc[sQd, sThrust]*tuning['kThrustA']
        #print('[INFO] KalmanModel: Thrust ddq relation: {}'.format(BFx_selected.loc['ddq_FA1']))
        if 'Thrust' in KF.sU:
#             if fullColumns:
#                 B.loc[sQd, 'Thrust'] = BFx_selected.loc[sQd]

            # For q_FA1
            m_twr = WT.twr.mass
            GM_twr = KF.WT.twr.MM[6,6] # Generalized mass of tower
            m_rna = WT.RNA.mass
            M_modal = GM_twr + m_rna
            #M_modal_OF = OF_lin['B'].loc['d_qt1FA_[m/s]', 'NacFxN1_[N]'] # Nacelle x force
            B[KF.iX['dq_FA1'], KF.iU['Thrust']] = 1 / M_modal
            #B[KF.iX['dq_FA1'], KF.iU['Thrust']] = BFx_selected.loc['ddq_FA1']

            # B matrix columns for Thrust (applicable to both YAMS and OpenFAST)
            B[KF.iX['dx']    , KF.iU['Thrust']] = M_inv[0, 0] * 1.0 + M_inv[0, 1] * h_hub
            B[KF.iX['dphi_y'], KF.iU['Thrust']] = M_inv[1, 0] * 1.0 + M_inv[1, 1] * h_hub

        # --- Generalized hydro force 
        IQD   =[KF.iX['dx'], KF.iX['dphi_y']]
        if 'k_h' in pHD:
            A[KF.iX['dx'],     KF.iX['dq_h']] = M_inv_sub[0, :] @ pHD['k_h']
            A[KF.iX['dphi_y'], KF.iX['dq_h']] = M_inv_sub[1, :] @ pHD['k_h']
        else:
            FAIL('k_h not present')

        # --- Outputs 
        # Monopile top acceleration
        C[KF.iY['PtfmIMUAx'], :] = A[KF.iX['dx'],:] # PtfmIMUAx is assumed to be qdd_s
        D[KF.iY['PtfmIMUAx'], :] = B[KF.iX['dx'],:] # PtfmIMUAx is assumed to be qdd_s
        # Inclination
        C[KF.iY['PtfmIncly'], KF.iX['phi_y']] = 1   # We measure inclination        
        # Rotational speed
        C[KF.iY['dpsi'], KF.iX['dpsi']] = 1  # We measure rotational speed

        # Nacelle acceleration in x direction (including pitch coupling)
        h_nac = r_TN[2]
        if method == 'YAMS':
            C[KF.iY['NcIMUAx'], :] = A[KF.iX['dx'], :] + h_nac * A[KF.iX['dphi_y'], :] + A[KF.iX['dq_FA1'], :]
            D[KF.iY['NcIMUAx'], :] = B[KF.iX['dx'], :] + h_nac * B[KF.iX['dphi_y'], :] + B[KF.iX['dq_FA1'], :]

        #C[KF.iY['NcIMUAz'], KF.iX['phi_y']] = -9.81
        D[KF.iY['Qgen'], KF.iU['Qgen']] = 1

        # --- Shaping filter, Hydro state equation
        if Tp==12.7:
            KF.Sw= 2.3835e-01
        elif Tp==10.0:
            KF.Sw= 2.3835e-01/2
        else:
            raise NotImplementedError(f'Tp={Tp}')
        KF.omega_p = 2*np.pi/Tp
        KF.zeta = zeta
        print('omega_p^2', KF.omega_p**2, '2 zeta omega_p', 2*KF.zeta*KF.omega_p)
        A[KF.iX['q_h'], KF.iX['dq_h']]  = 1
        A[KF.iX['dq_h'], KF.iX['q_h']]  = -KF.omega_p**2
        A[KF.iX['dq_h'], KF.iX['dq_h']] = -2 * KF.zeta * KF.omega_p
        B[KF.iX['dq_h'], KF.iU['w']] = 1 # White noise

        # --- Finally, we set the matrices
        KF.setMat(A, B, C, D)

    def _set_openfast_submatrix(KF, A, B, C, linFile):
        """Insert measured OpenFAST state and IMU couplings when available."""
        from welib.fast.FASTLin import FASTLin
        sX_sel  =['PtfmSurge_[m]', 'PtfmHeave_[m]', 'PtfmPitch_[rad]', 'qt1FA_[m]', 'psi_rot_[rad]']
        sX_sel +=['d_PtfmSurge_[m/s]', 'd_PtfmHeave_[m/s]', 'd_PtfmPitch_[rad/s]', 'd_qt1FA_[m/s]', 'd_psi_rot_[rad/s]']

        sU_sel=[]
        sU_sel+=['WS_[m/s]', 'alpha_[-]', 'WD_[rad]', 'SEAWaveElevRefPoint_[m]']
        sU_sel+=['PtfmFxN1_[N]', 'PtfmFyN1_[N]', 'PtfmFzN1_[N]', 'PtfmMxN1_[Nm]', 'PtfmMyN1_[Nm]', 'PtfmMzN1_[Nm]']
        sU_sel+=['TwrFxN1_[N]']
        sU_sel+=['TwrMyN1_[N]']
        sU_sel+=['TwrFxN20_[N]']
        sU_sel+=['TwrMyN20_[N]']
        sU_sel+=['HubFxN1_[N]', 'HubFyN1_[N]', 'HubFzN1_[N]', 'HubMxN1_[Nm]', 'HubMyN1_[Nm]', 'HubMzN1_[Nm]'] 
        sU_sel+=['NacFxN1_[N]', 'NacFyN1_[N]', 'NacFzN1_[N]', 'NacMxN1_[Nm]', 'NacMyN1_[Nm]', 'NacMzN1_[Nm]'] 
        sU_sel+=['B1pitch_[rad]', 'B2pitch_[rad]', 'B3pitch_[rad]']
        sU_sel+=['Qgen_[Nm]'] 
        sU_sel+=['PitchColl_[rad]']
        sU_sel+=['ADNacTxN1_[m]', 'ADNacTyN1_[m]', 'ADNacTzN1_[m]', 'ADNacRxN1_[rad]', 'ADNacRyN1_[rad]', 'ADNacRzN1_[rad]']
        sU_sel+=['ADHubTxN1_[m]', 'ADHubTyN1_[m]', 'ADHubTzN1_[m]', 'ADHubRxN1_[rad]', 'ADHubRyN1_[rad]', 'ADHubRzN1_[rad]']
        sU_sel+=['ADHubRVxN1_[rad/s]', 'ADHubRVyN1_[rad/s]', 'ADHubRVzN1_[rad/s]', 'ADTwrTxN1_[m]']
        sU_sel+=['ADWS_[m/s]', 'ADalpha_[-]', 'ADWD_[rad]']
#         sU_sel+=['HDPtfm-RefPtTxN1_[m]'      , 'HDPtfm-RefPtTyN1_[m]'      , 'HDPtfm-RefPtTzN1_[m]'      , 'HDPtfm-RefPtRxN1_[rad]'      , 'HDPtfm-RefPtRyN1_[rad]'      , 'HDPtfm-RefPtRzN1_[rad]']
#         sU_sel+=['HDPtfm-RefPtTVxN1_[m/s]'   , 'HDPtfm-RefPtTVyN1_[m/s]'   , 'HDPtfm-RefPtTVzN1_[m/s]'   , 'HDPtfm-RefPtRVxN1_[rad/s]'   , 'HDPtfm-RefPtRVyN1_[rad/s]'   , 'HDPtfm-RefPtRVzN1_[rad/s]']
#         sU_sel+=['HDPtfm-RefPtTAxN1_[m/s^2]' , 'HDPtfm-RefPtTAyN1_[m/s^2]' , 'HDPtfm-RefPtTAzN1_[m/s^2]' , 'HDPtfm-RefPtRAxN1_[rad/s^2]' , 'HDPtfm-RefPtRAyN1_[rad/s^2]' , 'HDPtfm-RefPtRAzN1_[rad/s^2]']
        sU_sel+=['HDWaveElevRefPoint_[m]', 'HDhorizontalcurrentspeed_[m/s]', 'HDalpha_[-]', 'HDWD_[rad]']
        sU_sel+=['SDTPMeshTxN1_[m]'      , 'SDTPMeshTyN1_[m]'      , 'SDTPMeshTzN1_[m]'      , 'SDTPMeshRxN1_[rad]'      , 'SDTPMeshRyN1_[rad]'      , 'SDTPMeshRzN1_[rad]']
        sU_sel+=['SDTPMeshTVxN1_[m/s]'   , 'SDTPMeshTVyN1_[m/s]'   , 'SDTPMeshTVzN1_[m/s]'   , 'SDTPMeshRVxN1_[rad/s]'   , 'SDTPMeshRVyN1_[rad/s]'   , 'SDTPMeshRVzN1_[rad/s]']
        sU_sel+=['SDTPMeshTAxN1_[m/s^2]' , 'SDTPMeshTAyN1_[m/s^2]' , 'SDTPMeshTAzN1_[m/s^2]' , 'SDTPMeshRAxN1_[rad/s^2]' , 'SDTPMeshRAyN1_[rad/s^2]' , 'SDTPMeshRAzN1_[rad/s^2]']
        sU_sel+=['SDLMeshFxN1_[N]'  , 'SDLMeshFyN1_[N]'  , 'SDLMeshFzN1_[N]'  , 'SDLMeshMxN1_[Nm]'  , 'SDLMeshMyN1_[Nm]'  , 'SDLMeshMzN1_[Nm]']
        sU_sel+=['SDLMeshFxN50_[N]' , 'SDLMeshFyN50_[N]' , 'SDLMeshFzN50_[N]' , 'SDLMeshMxN50_[Nm]' , 'SDLMeshMyN50_[Nm]' , 'SDLMeshMzN50_[Nm]']
        
        sY_sel=[]
        sY_sel+=['Extendedoutput:WS_[m/s]', 'Extendedoutput:alpha_[-]', 'Extendedoutput:WD_[rad]']
        sY_sel+=['Wind1VelX_[m/s]', 'Wind1VelY_[m/s]', 'Wind1VelZ_[m/s]']
        sY_sel+=['SEAExtendedoutput:WaveElevRefPoint_[m]', 'SEAWave1Elev_[m]']
        sY_sel+=['PtfmTxN1_[m]'      , 'PtfmTyN1_[m]'      , 'PtfmTzN1_[m]'      , 'PtfmRxN1_[rad]'      , 'PtfmRyN1_[rad]'      , 'PtfmRzN1_[rad]']
        sY_sel+=['PtfmTVxN1_[m/s]'   , 'PtfmTVyN1_[m/s]'   , 'PtfmTVzN1_[m/s]'   , 'PtfmRVxN1_[rad/s]'   , 'PtfmRVyN1_[rad/s]'   , 'PtfmRVzN1_[rad/s]']
        sY_sel+=['PtfmTAxN1_[m/s^2]' , 'PtfmTAyN1_[m/s^2]' , 'PtfmTAzN1_[m/s^2]' , 'PtfmRAxN1_[rad/s^2]' , 'PtfmRAyN1_[rad/s^2]' , 'PtfmRAzN1_[rad/s^2]']
        sY_sel+=['TwrTAxN1_[m/s^2]'  , 'TwrTAyN1_[m/s^2]'  , 'TwrTAzN1_[m/s^2]'  , 'TwrRAyN1_[rad/s^2]'  , 'TwrRAzN1_[rad/s^2]']
        sY_sel+=['TwrTAxN22_[m/s^2]' , 'TwrTAyN22_[m/s^2]' , 'TwrTAzN22_[m/s^2]' , 'TwrRAyN22_[rad/s^2]' , 'TwrRAzN22_[rad/s^2]']
        sY_sel+=['HubTxN1_[m]'      , 'HubTyN1_[m]'      , 'HubTzN1_[m]'       , 'HubRxN1_[rad]' , 'HubRyN1_[rad]' , 'HubRzN1_[rad]']
        sY_sel+=['HubRVxN1_[rad/s]' , 'HubRVyN1_[rad/s]' , 'HubRVzN1_[rad/s]']
        sY_sel+=['NacTxN1_[m]'      , 'NacTyN1_[m]'      , 'NacTzN1_[m]'      , 'NacRxN1_[rad]'      , 'NacRyN1_[rad]'      , 'NacRzN1_[rad]']
        sY_sel+=['NacTVxN1_[m/s]'   , 'NacTVyN1_[m/s]'   , 'NacTVzN1_[m/s]'   , 'NacRVxN1_[rad/s]'   , 'NacRVyN1_[rad/s]'   , 'NacRVzN1_[rad/s]']
        sY_sel+=['NacTAxN1_[m/s^2]' , 'NacTAyN1_[m/s^2]' , 'NacTAzN1_[m/s^2]' , 'NacRAxN1_[rad/s^2]' , 'NacRAyN1_[rad/s^2]' , 'NacRAzN1_[rad/s^2]']

        sY_sel +=['HSSSpd_[rad/s]', 'Azimuth_[deg]', 'BPitch1_[deg]', 'GenSpeed_[rpm]', 'RotSpeed_[rpm]']
        sY_sel +=['RotThrust_[kN]', 'RotTorq_[kNm]', 'RotPwr_[kW]']
        sY_sel +=['LSSTipMya_[kNm]', 'LSSTipMza_[kNm]']
        sY_sel +=['LSSGagMya_[kNm]', 'LSSGagMza_[kNm]']
        sY_sel +=['NcIMUTAxs_[m/s^2]', 'NcIMUTAys_[m/s^2]', 'NcIMUTAzs_[m/s^2]', 'NcIMUTVxs_[m/s]', 'NcIMUTVys_[m/s]', 'NcIMUTVzs_[m/s]']
        sY_sel +=['TTDspFA_[m]', 'TTDspSS_[m]'] 
        sY_sel +=['NacYaw_[deg]']
        sY_sel +=['YawBrFxp_[kN]', 'YawBrFyp_[kN]', 'YawBrFzp_[kN]', 'YawBrMxp_[kNm]', 'YawBrMyp_[kNm]', 'YawBrMzp_[kNm]']
        sY_sel +=['YawBrTDxt_[m]', 'YawBrTDyt_[m]'] 
        sY_sel +=['TwrBsFxt_[kN]', 'TwrBsFyt_[kN]', 'TwrBsFzt_[kN]', 'TwrBsMxt_[kNm]', 'TwrBsMyt_[kNm]', 'TwrBsMzt_[kNm]']
        sY_sel +=['PtfmSurge_[m]', 'PtfmSway_[m]', 'PtfmHeave_[m]', 'PtfmRoll_[deg]', 'PtfmPitch_[deg]', 'PtfmYaw_[deg]']
        sY_sel +=['Q_GeAz_[rad]'    , 'Q_TFA1_[m]'    , 'Q_TSS1_[m]'    , 'Q_TFA2_[m]'    , 'Q_TSS2_[m]'    , 'Q_Sg_[m]'    , 'Q_Hv_[m]'    , 'Q_P_[rad]']
        sY_sel +=['QD_GeAz_[rad/s]' , 'QD_TFA1_[m/s]' , 'QD_TSS1_[m/s]' , 'QD_TFA2_[m/s]' , 'QD_TSS2_[m/s]' , 'QD_Sg_[m/s]' , 'QD_Hv_[m/s]' , 'QD_P_[rad/s]']
        sY_sel +=['QD2_GeAz_[rad/s^2]', 'QD2_TFA1_[m/s^2]', 'QD2_Sg_[m/s^2]', 'QD2_Hv_[m/s^2]', 'QD2_P_[rad/s^2]']
        sY_sel +=['TwHt1ALxt_[m/s^2]' , 'TwHt1ALyt_[m/s^2]' , 'TwHt1ALzt_[m/s^2]']
        sY_sel +=['TwHt1MLxt_[kNm]'   , 'TwHt1MLyt_[kNm]'   , 'TwHt1MLzt_[kNm]']
        sY_sel +=['TwHt1FLxt_[kN]'    , 'TwHt1FLyt_[kN]'    , 'TwHt1FLzt_[kN]']
        sY_sel +=['ADNacFxN1_[N]', 'ADNacFyN1_[N]', 'ADNacFzN1_[N]', 'ADNacMxN1_[Nm]', 'ADNacMyN1_[Nm]', 'ADNacMzN1_[Nm]']
        sY_sel +=['ADHubFxN1_[N]', 'ADHubFyN1_[N]', 'ADHubFzN1_[N]', 'ADHubMxN1_[Nm]', 'ADHubMyN1_[Nm]', 'ADHubMzN1_[Nm]']
        sY_sel +=['ADRtAeroCp', 'ADRtAeroCq', 'ADRtAeroCt', 'ADRtAeroPwr', 'ADRtArea']
        sY_sel +=['ADRtSkew_[deg]', 'ADRtSpeed_[rpm]', 'ADRtTSR'] 
        sY_sel +=['ADRtAeroFxh_[N]', 'ADRtAeroMxh_[Nm]']
        sY_sel +=['ADRtVAvgxh_[m/s]', 'ADRtVAvgyh_[m/s]', 'ADRtVAvgzh_[m/s]']
        sY_sel +=['HDMorisonLoadsFxN1_[N]', 'HDMorisonLoadsFyN1_[N]', 'HDMorisonLoadsFzN1_[N]', 'HDMorisonLoadsMxN1_[Nm]', 'HDMorisonLoadsMyN1_[Nm]', 'HDMorisonLoadsMzN1_[Nm]']
        sY_sel +=['HDMorisonLoadsFxN99_[N]', 'HDMorisonLoadsFyN99_[N]', 'HDMorisonLoadsFzN99_[N]','HDMorisonLoadsMyN99_[Nm]', 'HDMorisonLoadsMzN99_[Nm]']
        sY_sel +=['HDHydroFxi_[N]', 'HDHydroFyi_[N]', 'HDHydroFzi_[N]', 'HDHydroMxi_[Nm]', 'HDHydroMyi_[Nm]', 'HDHydroMzi_[Nm]']
        sY_sel +=['SDInterfacedisplacementFxN1_[N]', 'SDInterfacedisplacementFyN1_[N]', 'SDInterfacedisplacementFzN1_[N]', 'SDInterfacedisplacementMxN1_[Nm]', 'SDInterfacedisplacementMyN1_[Nm]', 'SDInterfacedisplacementMzN1_[Nm]'] 
        #sY_sel +=# 'SDSSQM01', 'SDSSQMD01', 'SDSSQMDD01', 'SDSSQM02', 'SDSSQMD02', 'SDSSQMDD02',
        sY_sel +=['SDIntfFxss_[N]'   , 'SDIntfFzss_[N]'   , 'SDIntfMyss'       , 'SDIntfTDXss_[m]' ]
        sY_sel +=['SDIntfTDZss_[m]' , 'SDIntfRDYss_[rad]']
        sY_sel +=['SDM1N1FKxe_[N]'   , 'SDM1N1FKze_[N]'   , 'SDM1N1MKye'       , 'SDM2N1FKxe_[N]'  , 'SDM2N1FKze_[N]']
        sY_sel +=['SDM49N1FKxe_[N]'  , 'SDM49N1FKze_[N]'  , 'SDM49N1MKye'      , 'SDM49N2FKxe_[N]' , 'SDM49N2FKze_[N]' , 'SDM49N2MKye']
        sY_sel +=['SD-ReactFxss_[N]' , 'SD-ReactFyss_[N]' , 'SD-ReactFzss_[N]' , 'SD-ReactMxss'    , 'SD-ReactMyss'    , 'SD-ReactMzss']

        pklFile = linFile.replace('.lin', '.pkl')
        if os.path.exists(pklFile):
            INFO('Loading lin file pickle: ', pklFile)
            FL = FASTLin.from_pickle(pklFile)
        else:
            FL = FASTLin(linfiles=[linFile], prefix='', verbose=False, sX_sel=sX_sel, sU_sel=sU_sel, sY_sel=sY_sel)
            FL.save(pklFile)
        data = FL.OP_Data[0].Data[0] # Instance of FASTLinearizationFile with keys A, B, C, D, x, u, y, x_info
        sX = [str(label) for label in FL.xdescr] # State names
        sY = [str(label) for label in FL.ydescr] # Output names
        sU = [str(label) for label in FL.udescr] # Input names
        state_map = {
            'x'      : 'PtfmSurge_[m]',
            'phi_y'  : 'PtfmPitch_[rad]',
            'q_FA1'  : 'qt1FA_[m]',
            'psi'    : 'psi_rot_[rad]',
            'dx'     : 'd_PtfmSurge_[m/s]',
            'dphi_y' : 'd_PtfmPitch_[rad/s]',
            'dq_FA1' : 'd_qt1FA_[m/s]',
            'dpsi'   : 'd_psi_rot_[rad/s]'}
        indices = {name: sX.index(label) for name, label in state_map.items() if label in sX}
        for row_name, row in indices.items():
            for col_name, col in indices.items():
                A[KF.iX[row_name], KF.iX[col_name]] = data.A.values[row, col]
        output_map = {
                'PtfmIMUAx': 'QD2_Sg_[m/s^2]',
                'PtfmIncly': 'Q_P_[rad]',
                'NcIMUAx': 'NcIMUTAxs_[m/s^2]',
                      }
        for output_name, label in output_map.items():
            if label in sY:
                row = sY.index(label)
                for state_name, col in indices.items():
                    C[KF.iY[output_name], KF.iX[state_name]] = data.C.values[row, col]
            else:
                WARN('Label missing from sY in lin model: ', label)

        input_map = {
                'Qgen': 'Qgen_[Nm]',
                # 'eta':  'HDWaveElevRefPoint_[m]'
        }
        for input_name, label in input_map.items():
            if label in sU:
                col = sU.index(label)
                for state_name, row in indices.items():
                    B[KF.iX[state_name], KF.iU[input_name]] = data.B.values[row, col]
            else:
                WARN('Label missing from sU in lin model: ', label)

        OF_lin = data.toDataFrame() 
        return OF_lin

        # --- OUTPUTS that we could use
#               1    0.00000000000E+00                                                 F               0         SEA Wave1Elev, (m)
#               7    3.57159790039E+02                                                 F               0         ED RotThrust, (kN)
#               8   -1.33902830157E-07                                                 F               0         ED RotTorq, (kN-m)
#               9   -1.44830606878E-02                                                 F               0         ED NcIMUTAxs, (m/s^2)
#              23    1.16558103561E+01                                                 F               0         ED TwrBsFxt, (kN)
#              25   -2.74127949219E+04                                                 F               0         ED TwrBsFzt, (kN)
#              27   -1.10646890625E+05                                                 F               0         ED TwrBsMyt, (kN-m)
#              35    4.71238899231E+00                                                 F               0         ED Q_GeAz, (rad)
#              36   -4.67881500721E-01                                                 F               0         ED Q_TFA1, (m)
#              40   -2.05030087382E-02                                                 F               0         ED Q_Sg, (m)
#              42   -8.86160356458E-04                                                 F               0         ED Q_P, (rad)
#              51   -5.10188657790E-03                                                 F               0         ED QD2_Sg, (m/s^2)
#              54    5.57291064453E+03                                                 F               0         HD HydroFxi, (N)
#              58   -5.41483320312E+04                                                 F               0         HD HydroMyi, (N-m)
#              66    3.59479296875E+04                                                 F               0         SD IntfFXss, (N)
#              68   -1.10646888000E+08                                                 F               0         SD IntfMYss, (N*m)
#              69   -2.05030087382E-02                                                 F               0         SD IntfTDXss, (m)
#              71   -8.86160589289E-04                                                 F               0         SD IntfRDYss, (rad)
#              72    4.39176328125E+04                                                 F               0         SD M1N1FKXe, (N)
#              74   -1.09302344000E+08                                                 F               0         SD M1N1MKYe, (N*m)
#             219    3.72363359375E+04                                                 F               0         SD M49N1FKXe, (N)
#             220   -2.84774700000E+07                                                 F               0         SD M49N1FKZe, (N)
#             221   -1.11201280000E+08                                                 F               0         SD M49N1MKYe, (N*m)
#             222    3.72363359375E+04                                                 F               0         SD M49N2FKXe, (N)
#             223   -2.84774700000E+07                                                 F               0         SD M49N2FKZe, (N)
#             224   -1.11238520000E+08                                                 F               0         SD M49N2MKYe, (N*m)
#             225    4.39180976562E+04                                                 F               0         SD -ReactFXss, (N)
#             227   -3.67519840000E+07                                                 F               0         SD -ReactFZss, (N)
#             229   -1.09302344000E+08                                                 F               0         SD -ReactMYss, (N*m)

    # --- Methods From Parent Class
    # loadMeasurements 
    # prepareTimeStepping 
    def prepareTimeStepping(KF, *args, **kwargs):
        KalmanFilter.prepareTimeStepping(KF, *args, **kwargs)
        KF.dInfo = KF.WT.calcOutputs_init(time=KF.time)

    # setupCovariances
        
    # --- Methods Common between TN and TNLin
    # prepareMeasurements


    def get_OF_DOFs(KF, x, x_dot=None):
        q   = KF.dInfo['q_default'].copy()
        qd  = KF.dInfo['q_default'].copy()
        qdd = KF.dInfo['q_default'].copy()
        # Model specific
        #        0      1         2      3      4      5       6          7
        #sQ  = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi'] # Mechanical states
        fnd_x_q   = np.array([x[0],x[1]])
        fnd_xd_q  = np.array([x_dot[0], x_dot[1]])
        fnd_xdd_q = np.array([x_dot[4], x_dot[5]])

        q  ['Sg']   = x[0]
        q  ['P']    = x[1]
        q  ['TFA1'] = x[2]
        q  ['Psi']  = x[3]

        qd ['Sg']   = x[4]
        qd ['P']    = x[5]
        qd ['TFA1'] = x[6]
        qd ['Psi']  = x[7]
        if x_dot is not None:
            qdd['Sg']   = x_dot[4]
            qdd['P']    = x_dot[5]
            qdd['TFA1'] = x_dot[6]
            qdd['Psi']  = x_dot[7]
        return q, qd, qdd, fnd_x_q, fnd_xd_q, fnd_xdd_q

    def get_OF_DOFs_clean(KF, it):
        #q   = dInfo['Q'].iloc[it,:].copy() 
        #qd  = dInfo['QD'].iloc[it,:].copy()
        #qdd = dInfo['QDD'].iloc[it,:].copy()
        q   = KF.dInfo['q_default'].copy()
        qd  = KF.dInfo['q_default'].copy()
        qdd = KF.dInfo['q_default'].copy()
        q['Sg']     = KF.df['x'].iloc[it]
        q['P']      = KF.df['phi_y'].iloc[it]
        q['TFA1']   = KF.df['q_FA1'].iloc[it]
        q['Psi']    = KF.df['psi'].iloc[it]
        qd['Sg']    = KF.df['dx'].iloc[it]
        qd['P']     = KF.df['dphi_y'].iloc[it]
        qd['TFA1']  = KF.df['dq_FA1'].iloc[it]
        qd['Psi']   = KF.df['dpsi'].iloc[it]
        qdd['Sg']   = KF.df['ddx'].iloc[it]
        qdd['P']    = KF.df['ddphi_y'].iloc[it]
        qdd['TFA1'] = KF.df['ddq_FA1'].iloc[it]
        qdd['Psi']  = KF.df['ddpsi'].iloc[it]
        return q, qd, qdd


    def computeSectionLoads(KF, t, x, x_dot, p_hydro, Thrust, Qaero, it=None):
        WT = KF.WT
        dInfo = KF.dInfo
        # ---
        p_ext      = np.zeros((3,len(KF.zDepth)))
        p_ext[0,:] = p_hydro

        # --- DOFs in the way expected by calcOutputs_step
        q, qd, qdd, x_q, xd_q, xdd_q = KF.get_OF_DOFs(x, x_dot=x_dot)

        # ---  Section Loads
        F_top = np.array((0.,0.,0.))
        M_top = np.array((0.,0.,0.))
        a_ext = np.array((0.,0.,-WT.gravity)) # external acceleration (gravity/earthquake)

        Qaero = 0
        ser_Loads = pd.Series({'Fadd_R_xs':Thrust, 'Fadd_R_ys':0, 'Fadd_R_zs':0, 'Madd_R_xs':Qaero, 'Madd_R_ys':0, 'Madd_R_zs':0})

        if KF.hacks['SL_cleanQ']:
            q, qd, qdd = KF.get_OF_DOFs_clean(it+1)

        if KF.hacks['SL_cleanFtop']:
            ser_Loads_Add =ser_Loads
            ser_Loads = KF.df.loc[it+1] # will pick up the YawBr loads
            dInfo['useTopLoadsFromDF'] = True
            # We add the "Fadd" just so that YawBr looks better, even though it's overwritten
            for key, val in ser_Loads_Add.items():
                ser_Loads[key] = val

        if KF.hacks['SL_cleanEtaDot']:
            eta_dot_true = WT.pSS['eta_dot'][it+1]
            p_hydro = KF.pHD['phi'] * eta_dot_true # p_h = k_h(z) q_h(t)
            p_ext[0,:] = p_hydro
            #p_ext = None

        p_ext_for_h = p_ext
        if KF.hacks['SL_cleanP']:
            p_ext = None  # if p_ext is None, WT will compute the p_ext based on the sea state, it's cheating

        rowOut, twr_sec, mnp_sec = WT.calcOutputs_step(q, qd, qdd, dInfo, t=t, ser_Loads=ser_Loads, mnp_p_ext=p_ext)

        # No acceleration # TODO get it from calcOutputs 
        xdd_q *=0
        F_sec_h, M_sec_h, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext_for_h, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext = a_ext, PhiK=WT.fnd.PhiK)
        Fx_h_est = F_sec_h[0,0]
        #wet_nodes = KF.zDepth <= 0
        #Fx_h_est = np.trapezoid(p_hydro[wet_nodes], KF.zDepth[wet_nodes])

        return rowOut, twr_sec, mnp_sec, Fx_h_est

    def timeLoop(KF):
        # --- Aliases to shorten notations
        WT = KF.WT
        dInfo = KF.dInfo

        # Prepare section output calculation
        KF.dfOut = dInfo['dfOut']
        
        # --- Initial conditions
        x = KF.initFromClean(var='x,y,u')
        
        Thrust = KF.U_clean['Thrust'].iloc[0]
        Fx_i   = KF.U_clean['Fx_i'].iloc[0]
        My_i   = KF.U_clean['My_i'].iloc[0]
        # --- WSE
        WS_last = KF.S_clean['WS'].iloc[0]

        # --- Section loads at t=0
        #x= x.values
        #x_dot = x*0
        #p_hydro = KF.pHD['phi'] * 0
        #print(x)
        #print(x_dot)
        #rowOut, twr_sec, mnp_sec, Fx_h_est = KF.computeSectionLoads(KF.time[0], x, x_dot, p_hydro, Thrust, Qaero=0, it=0)

        # --- Time loop
        for it in range(0, KF.nt-1):
            t = KF.time[it]
            # --- "Measurements"
            y = KF.Y.iloc[it,:].values
            # --- Inputs
            u = KF.U_clean.iloc[it,:].values.copy()
            u[KF.iU['Thrust']] = Thrust # We use previous estimated thrust as input.
            u[KF.iU['Fx_i']] = Fx_i # We use previous estimated interface loads
            u[KF.iU['My_i']] = My_i # We use previous estimated interface loads
            
            # --- Predictions of next time step based on current time step
            t = KF.time[it+1]
            x, KF.P, _ = KF.estimateTimeStep(u, y, x, KF.P)

            # --- Estimate Wind Speed
            if KF.hacks['WSE'] == 'clean_inputs':
                Qaero_hat = KF.X_clean['Qaero'].iloc[it]
                omega = KF.X_clean['dpsi'].iloc[it]
            else:
                Qaero_hat = x[KF.iX['Qaero']]
                omega     = x[KF.iX['dpsi']]
            pitch = u[KF.iU['pitch']] * 180 / np.pi # deg
            WS_hat, _ = KF.wse.estimate(Qaero_hat, pitch=pitch, omega=omega, WS0=WS_last, relaxation=0, method='oper-crossing', t=KF.time[it])
            Qaero_hat = np.max(Qaero_hat,0)
            WS_last = float(WS_hat)
            
            # --- Estimate Thrust
            if KF.hacks['thrust'] == 'clean':
                Thrust = KF.U_clean['Thrust'].iloc[it]
            else:
                Thrust = KF.wse.Thrust(WS_hat, pitch=pitch, omega=omega)

            # --- Estimate Generalized hydro force and bending moment (calc output)
            q_h     = x[KF.iX['q_h']]  # eta
            dq_h    = x[KF.iX['dq_h']] # eta_dot
            eta     = q_h
            eta_dot = dq_h
            p_hydro = KF.pHD['phi'] * eta_dot # p_h = k_h(z) q_h(t)
            p_hydro[KF.zDepth>0] = 0 # safety, shoudn't be necessary

            # --- Acceleration
            x_dot = np.dot(KF.A, x) + np.dot(KF.B, u)
            
            # --- Section Loads
            rowOut, twr_sec, mnp_sec, Fx_h_est = KF.computeSectionLoads(KF.time[it+1], x, x_dot, p_hydro, Thrust, Qaero_hat, it=it)
            KF.dfOut.loc[it+1] = rowOut
            # Section loads at interface for next time step
            Fx_i = twr_sec[0, 0]  # rowOut['TwrBsFxt_[kN]']
            My_i = twr_sec[4, 0]  # rowOut['TwrBsMyt_[kN-m]']
            
            
            # --- Sotre "updated"/"hacked" states and inputs
            if 'psi' in KF.iX:
                x[KF.iX['psi']] = np.mod(x[KF.iX['psi']], 2*np.pi)
            KF.U_hat.iloc[it+1,:]   = u
            KF.X_hat.iloc[it+1,:]   = x
            KF.XD_hat.iloc[it+1,:]  = x_dot

            # --- Store extra info
            # Environment
            KF.S_hat.at[it+1, 'WS']     = WS_hat
            KF.S_hat.at[it+1, 'eta']    = q_h
            # Loads
            KF.S_hat.at[it+1, 'Fx_sb']  = mnp_sec[0,0]
            KF.S_hat.at[it+1, 'My_sb']  = mnp_sec[4,0]
            KF.S_hat.at[it+1, 'Fx_h']   = Fx_h_est
            KF.S_hat.at[it+1, 'Thrust'] = Thrust

            # --- Propagation to next time step
            # --- Print status to screen
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it, KF.time[it], WS_hat, Thrust))

        # TODO evaluate at t=0, for now we just replicate the value
        index = KF.dfOut.index
        cols = KF.dfOut.columns.difference(['Time_[s]'])
        KF.dfOut.loc[0, cols] = KF.dfOut.loc[1, cols]
        KF.dfOut.loc[index[-1], cols] = KF.dfOut.loc[index[-2], cols]



    def calc_sectionLoads(KF, clean=True):
        # --- Aliases to shorten notations
        WT = KF.WT
        dInfo = KF.dInfo

        # Prepare section output calculation
        KF.dfOut = dInfo['dfOut']
        
        # --- Time loop
        for it in range(0, KF.nt):    
            # --- DOFs in the way expected by calcOutputs_step
            if not clean:
                raise Exception()
                #q, qd, qdd, x_q, xd_q, xdd_q = KF.get_OF_DOFs(x, x_dot=x_dot)
            else:
                Thrust = KF.U_clean['Thrust'].iloc[it]
                Qaero = 0
                ser_Loads_Add = pd.Series({'Fadd_R_xs':Thrust, 'Fadd_R_ys':0, 'Fadd_R_zs':0, 'Madd_R_xs':Qaero, 'Madd_R_ys':0, 'Madd_R_zs':0})

                q, qd, qdd = KF.get_OF_DOFs_clean(it)
                ser_Loads = KF.df.loc[it] # will pick up the YawBr loads

                # We add the "Fadd" just so that YawBr looks better, even though it's overwritten
                for key, val in ser_Loads_Add.items():
                    ser_Loads[key] = val
                dInfo['useTopLoadsFromDF'] = True

            p_ext = None

            rowOut, twr_it, mnp_it = WT.calcOutputs_step(q, qd, qdd, dInfo, t=KF.time[it], ser_Loads=ser_Loads, mnp_p_ext=p_ext)
            KF.dfOut.loc[it] = rowOut
        return KF.dfOut



"""Monopile/turbine digital twin using augmented Kalman estimation."""



scriptDir = os.path.dirname(__file__)

def main(fstFile, tmin=0, tmax=20, show=False,
         compFile=None, hydroShapeFile=None, 
         aeroMapFile=None, operFile=None,
         linFile=None, 
         method='YAMS', 
         hacks=None,
         Tp=None,
         nUnderSamp=10,
         tRangeStats=None):
    if tRangeStats is None:
         tRangeStats = [tmin, tmax]

    base = os.path.splitext(fstFile)[0] + '_DigitalTwin'

    # --- Read reference output DataFrame 
    outFile = fstFile.replace('.fst', '.outb')
    if not os.path.exists(outFile):
        raise FileNotFoundError(outFile)
    df_ref = weio.read(outFile).toDataFrame()
    df_ref = df_ref[(df_ref['Time_[s]'] >= tmin) & (df_ref['Time_[s]'] <= tmax)]

    df_ref=df_ref.iloc[::nUnderSamp,:] # reducing sampling
    df_ref.reset_index(inplace=True)

    # Read output file and ensure all mapped columns exist (e.g. for onshore case)
    missing_cols = ['Wave1Elev_[m]', 'HydroFxi_[N]', 'HydroMyi_[N-m]', '-ReactMYss_[N*m]', '-ReactFXss_[N]']
    for col in missing_cols:
        if col not in df_ref.columns:
            WARN('Missing column '+col)
            df_ref[col] = 0.0

    # --------------------------------------------------------------------------------}
    # --- Kalman filter estimation 
    # --------------------------------------------------------------------------------{
    # --- Wind speed estimator (reads tabulated aerodynamic data)
    wse = TabulatedWSEstimator(fstFile=fstFile, operFile=operFile, aeroMapFile=aeroMapFile)
    KF = KalmanFilterMTNS(WSE=wse, hacks=hacks)
    KF.setup_matrices(fstFile, 
                      compFile=compFile, hydroShapeFile=hydroShapeFile, Tp=Tp,
                      dfTime=df_ref['Time_[s]'].values, method=method, linFile=linFile,
                      )

    # --- Loading "Measurements"
    # - Reference file is opened
    # - Measurements are extracted from it
    # - Other signals are extracted from the file, for comparison with estimates. These are referred as "clean" values
    # - Estimate sigmas from measurements (overriden in next section)
    KF.loadMeasurements(measFile=df_ref, tRange=[tmin,tmax], colMap=KF.colMap, timeCol='Time_[s]', raiseIfAbsent=True)

    KF.X_clean['dq_h'] = np.gradient(KF.X_clean['q_h'], KF.dt)
    # --- Storage for plot
    KF.prepareTimeStepping()
    # --- Process and measurement covariances
    dt_ref = 0.02 # NOTE: Q change with dt
    # NOTE: sigs will be squared for P, Q, R
    sigs = {'x':{}, 'y':{}, 'Q':{}}
    sigs['y']['PtfmIMUAx'] = np.sqrt(1e-3)    
    sigs['y']['PtfmIncly'] = np.sqrt(2.7e-7)
    sigs['y']['dpsi']      = np.sqrt(1e-5)
    sigs['y']['NcIMUAx']  = np.sqrt(1e-2)    
    sigs['y']['Qgen']      = np.sqrt(1e-5)
#     sQ  = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi'] # Mechanical states
#     sQa = ['q_h', 'dq_h', 'Qaero']                                         # Augmented states
    sigs['Q']['x']      = np.sqrt(KF.dt/dt_ref * 2e-6)
    sigs['Q']['phi_y']  = np.sqrt(KF.dt/dt_ref * 2e-6)
    sigs['Q']['q_FA1']  = np.sqrt(KF.dt/dt_ref * 2e-5)
    sigs['Q']['psi']    = np.sqrt(KF.dt/dt_ref * 2e-6)

    sigs['Q']['dx']     = np.sqrt(KF.dt/dt_ref * 2e-3)
    sigs['Q']['dphi_y'] = np.sqrt(KF.dt/dt_ref * 2e-6)
    sigs['Q']['dq_FA1'] = np.sqrt(KF.dt/dt_ref * 2e-5)
    sigs['Q']['dpsi']   = np.sqrt(KF.dt/dt_ref * 2e-5)

    sigs['Q']['q_h']    = np.sqrt(KF.dt/dt_ref * 1e-6)
    sigs['Q']['dq_h']   = np.sqrt(KF.dt/dt_ref * 1e-4) #KF.Sw
    sigs['Q']['Qaero']  = np.sqrt(KF.dt/dt_ref * 1e10)
    sigs['x'] = sigs['Q'].copy()

    KF.setupCovariances(
            sigs=sigs,
            useDt=False, Pidentity=True, verbose=False)
    # TODO use sigs above instead
#     KF.R[:] = np.diag([1e-3, 2.7e-7, 1e-5, 1e-2, 1e4])
#     KF.Q[:] = np.diag([1e-6, 1e-6, 1e-5, 1e-6, 1e-5, 1e-6, 1e-5, 1e-5,
#                        1e-6, 0.0001, 1e10])

    # --- Prepare measurements - Create noisy measurements
    KF.setYFromClean(R=KF.R, NoiseRFactor=0)

    print('>>>>>>>>>>>>>>>> KALMAN FILTER >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    print(KF)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    # --------------------------------------------------------------------------------}
    # --- Time Loop 
    # --------------------------------------------------------------------------------{
    with Timer('KF time loop'):
        KF.timeLoop()
    df_sl = KF.dfOut
    file_sl = base + '_SectionLoads_KF_timeloop.outb'
    df_sl.to_outb(file_sl)
    print('Export:', file_sl)

    statsDict = {}    
    try:
        fig = KF.plot_X( printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        plt.savefig(base + '_KF_X.png')
        KF.plot_Y(printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        plt.savefig(base + '_KF_Y.png')
        fig = KF.plot_S(printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        plt.savefig(base + '_KF_S.png')
#         KF.plot_U()
    except Exception as e:
        FAIL('Plotting using KF plot functions failed:'+str(e))
    # KF.plot_P()
    # KF.plot_K()
    # KF.plot_innovation()
    
    if show:
        plt.show()
    return KF, df_ref, df_sl


if __name__ == '__main__':
    fstFile        = os.path.join(scriptDir, 'examples/_simulations/06_Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst')
#     linFile        = os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_OnlyWriteOutputs.1.lin')
    linFile        = os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1.1.lin')
    compFile       = os.path.join(scriptDir, 'examples/_simulations/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv')
    aeroMapFile    = os.path.join(scriptDir, 'examples/_simulations/IEA-22-280-RWT/IEA-22-280-RWT_Cp_Ct_Cq.rpf')
    operFile       = os.path.join(scriptDir, 'examples/_simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv')
    hydroShapeFile = os.path.join(scriptDir, 'examples/_data/IEAMonoPile_HydroShapeFunction_Hs=8.1_Tp=12.7.csv')
    hacks = {}

    # Super hack
#     hacks = {'thrust':'clean', 'WSE':'clean_inputs', 'SL_cleanQ':True, 
#              'SL_cleanFtop':True, 'SL_cleanEtaDot':True, 'SL_cleanP':True}
# 
#     # Intermediate hack: States are exact - Hydro loads are exact
#     hacks = {'SL_cleanQ':True, 'SL_cleanEtaDot':True, 'SL_cleanP':True} 

    # Intermediate hack
    #hacks = {'SL_cleanQ':True, 'SL_cleanEtaDot':True} # <<<< EXAMPLE

    method='YAMS'
    method='OpenFAST'
    show=True
    tRange = [190, 200]
    nUnderSamp=10
    #tRange = [150, 270]

    main(fstFile=fstFile, linFile=linFile, 
         compFile=compFile,  hydroShapeFile=hydroShapeFile, Tp=12.7,
         aeroMapFile=aeroMapFile, operFile=operFile,
         hacks=hacks, show=show,
         nUnderSamp=nUnderSamp,
         tmin=tRange[0], tmax=tRange[1],
         method=method
         )
