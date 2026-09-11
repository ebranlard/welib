"""
Kalman filter model for "Monopile Tower Nacelle Shaft" (based on yams MTNSB)

"""
import os
import numpy as np
from welib.kalman.kalman import *
from welib.kalman.kalmanfilter import KalmanFilter
from welib.ws_estimator.tabulated import TabulatedWSEstimator
from welib.yams.models.MTNSB import FASTmodel2MTNSB
import welib.fast.fastlib as fastlib
import welib.weio as weio

# --------------------------------------------------------------------------------}
# --- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterMTNS(KalmanFilter):

    def __init__(KF, WSE=None, debug=False):
        # 11 States
        sQ  = ['q_s', 'q_p', 'q_FA1', 'psi', 'qd_s', 'qd_p', 'qd_FA1', 'dpsi'] # Mechanical states
        sQa = ['q_h', 'qd_h', 'Qaero']                                         # Augmented states
        sU  = ['Qgen', 'pitch', 'Thrust', 'Fx_i', 'My_i', 'w']
        sY  = ['ddq_s', 'q_p', 'dpsi', 'NcIMUAx', 'NcIMUAy', 'NcIMUAz', 'Qgen']
        sS  = ['M_sb', 'F_sb', 'eta', 'Fx_h', 'WS', 'Thrust']
        # --- Parent init
        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)
        # --- Storing wind speed estimator (based on tabulated aerodynamic data)
        KF.wse = WSE # wind speed estimator
        KF.debug = debug
        # Hacks
        KF.hacks={'thrust':None, 'WSE':None}


    def setup_matrices(KF, 
          fstFilename, comp_file=None, hydro_shape_file=None, dfTime=None,
                       Tp=12.7, method='YAMS', lin_file=None):
        """ Build WT model (sea state, hydro) and state matrices A,B,C,D """
        from welib.kalman.kalman import BuildSystem_Linear_MechOnly
        from welib.fast.FASTLin import FASTLin
        from welib.tools.strings import FAIL
        # --- Default arguments

        # --- Windturbine model
        WT = FASTmodel2MTNSB(fstFilename, shapes_sub=[0, 4], shapes_twr=[], shapes_bld=[],
                             DEBUG=False, bStiffening=True, main_axis='z', fixedShaft=True,
                             algo='OpenFAST').WT

		# --- ColMap
        KF.colMap={
                'q_s'    : 'Q_Sg_[m]' ,
                'q_p'    : 'Q_P_[rad]' ,
                'q_FA1'  : 'Q_TFA1_[m]',
                'psi'    : '{Azimuth_[deg]} * np.pi/180',
                'qd_s'   : 'QD_Sg_[m/s]' ,
                'qd_p'   : 'QD_P_[rad/s]',
                'qd_FA1' : 'QD_TFA1_[m/s]',
                'dpsi'   : '{RotSpeed_[rpm]} * 2*np.pi/60',
                'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                'ddq_s'  : 'QD2_Sg_[m/s^2]',
                'NcIMUAx': 'NcIMUTAxs_[m/s^2]',
                'NcIMUAy': 'NcIMUTAys_[m/s^2]',
                'NcIMUAz': 'NcIMUTAzs_[m/s^2]',
                'Qgen'   : '{GenTq_[kN-m]} * 1000',
                'pitch'  : '{BldPitch1_[deg]} * np.pi/180',
                'Thrust' : 'RtAeroFxh_[N]',
                'Fx_i'   : 'HydroFxi_[N]',     # TODO TODO TODO THIS IS WRONG
                'My_i'   : 'HydroMyi_[N-m]',   # TODO TODO TODO THIS IS WRONG
                'WS'     : 'RtVAvgxh_[m/s]',
                'Fx_h'   : '{HydroFxi_[N]}',   # We use a trick  
                'F_sb'   : '-ReactFXss_[N]',
                'M_sb'   : '-ReactMYss_[N*m]',

            }



        # --- Configure Sea State components and wave elevation
        if WT.pSS is not None:
            print('[INFO] Setting Components', comp_file)
            WT.SS_setComponents(comp_file)
            print('[INFO] Setting Compute Eta')
            WT.SS_computeEta(dfTime)
            if hydro_shape_file is not None:
                WT.HD_setShapeFunction(hydro_shape_file)

        # --- Structural data used for postprocessing (section loads)
        pST = {}
        pST['PhiU'] = WT.fnd.PhiU
        pST['PhiV'] = WT.fnd.PhiV
        pST['PhiK'] = WT.fnd.PhiK
        pST['m'] = WT.fnd.m
        pST['s_span'] = WT.fnd.s_span
        pST['gravity'] = WT.gravity
        if hasattr(WT, 'WtrDpth'):
            pST['z'] = WT.fnd.s_span - WT.WtrDpth
        else:
            pST['z'] = WT.fnd.s_span

        # --- Hydrodynamic parameters (mock values for onshore case)
        if WT.pSS is None:
            zDepth = WT.fnd.s_span if hasattr(WT, 'fnd') else WT.twr.s_span
            pHD_mock = {
                'zDepth': zDepth,
                'D': np.zeros_like(zDepth), 'Cd': np.zeros_like(zDepth), 'CM': np.zeros_like(zDepth),
                'Ca': np.zeros_like(zDepth), 'Cp': np.zeros_like(zDepth), 'm_hydro': np.zeros_like(zDepth),
                'GM_hydro': np.zeros((2, 2)), 'phi': np.zeros_like(zDepth), 'phit': np.zeros_like(zDepth)
            }
            WT.pHD = pHD_mock
        pHD = WT.pHD

        # Ensure WT.MM contains the hydrodynamic mass if not already added
        if not getattr(WT, '_GM_hydro_added', False) and 'GM_hydro' in pHD:
            for i in range(min(pHD['GM_hydro'].shape[0], WT.MM.shape[0])):
                WT.MM[i, i] += pHD['GM_hydro'][i, i]
            WT._GM_hydro_added = True

        # --- Parameters that are a function of the structure and ocean conditions
        KF.WT = WT
        KF.pHD = pHD
        KF.pST = pST
        KF.zDepth = pST['z']

        # --- Setup state matrices, problem specific!
        nX, nU, nY = len(KF.sX), len(KF.sU), len(KF.sY)
        A = np.zeros((nX, nX))
        B = np.zeros((nX, nU))
        C = np.zeros((nY, nX))
        D = np.zeros((nY, nU))
        M_inv = np.linalg.inv(WT.MM)

        if method=='YAMS':
            As, _, _, _ = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK)

            # --- Structural DOFs
            names = ['q_s', 'q_p', 'qd_s', 'qd_p']
            for row, row_name in enumerate(names):
                for col, col_name in enumerate(names):
                    A[KF.iX[row_name], KF.iX[col_name]] = As[row, col]

            # Tower fore-aft first mode (1DOF oscillator)
            omega_t = 2 * np.pi * 0.32
            A[KF.iX['q_FA1'], KF.iX['qd_FA1']] = 1
            A[KF.iX['qd_FA1'], KF.iX['q_FA1']] = -omega_t**2
            A[KF.iX['qd_FA1'], KF.iX['qd_FA1']] = -2 * 0.03 * omega_t

        elif method=='OpenFAST':
            if lin_file is None or not os.path.exists(lin_file):
                raise FileNotFoundError('An OpenFAST .lin file is required for method=OpenFAST')
            KF._set_openfast_submatrix(A, B, C, lin_file)

        # --------------------------------------------------------------------------------}
        # ---  Code common to OpenFAST and YAMS
        # --------------------------------------------------------------------------------{
        # After applying the "generic" A, we introduce the augmented states equations
        # and some additional tweaks

        # --- Thrust influence (physical derivation)
        r_TN = WT.r_TN_inT if hasattr(WT, 'r_TN_inT') else [0,0,198.386]
        r_NS = WT.r_NS_inN if hasattr(WT, 'r_NS_inN') else [0,0,4.143]
        h_hub = r_TN[2] + r_NS[2]

        # B matrix columns for Thrust (applicable to both YAMS and OpenFAST)
        B[KF.iX['qd_s'], KF.iU['Thrust']] = M_inv[0, 0] * 1.0 + M_inv[0, 1] * h_hub
        B[KF.iX['qd_p'], KF.iU['Thrust']] = M_inv[1, 0] * 1.0 + M_inv[1, 1] * h_hub

        # For q_FA1, estimate modal mass: 0.25 * m_twr + m_rna
        m_twr = np.trapezoid(WT.twr.m, WT.twr.s_span) if (hasattr(WT, 'twr') and hasattr(WT.twr, 'm')) else 9.6e5
        m_rna = WT.M_RNA if hasattr(WT, 'M_RNA') else 1.189e6
        M_modal = 0.25 * m_twr + m_rna
        B[KF.iX['qd_FA1'], KF.iU['Thrust']] = 1 / M_modal

        # --- Shaft equation
        A[KF.iX['psi'], KF.iX['dpsi']] = 1
        inertia = float(np.asarray(WT.rot.inertia).ravel()[0])
        A[KF.iX['dpsi'], KF.iX['Qaero']] = 1 / inertia
        B[KF.iX['dpsi'], KF.iU['Qgen']] = -1 / inertia

        # --- Generalized hydro force 
        if 'k_h' in pHD:
            A[KF.iX['qd_s'], KF.iX['qd_h']] = M_inv[0, :] @ pHD['k_h']
            A[KF.iX['qd_p'], KF.iX['qd_h']] = M_inv[1, :] @ pHD['k_h']
        else:
            FAIL('k_h not present')

        # --- Output equation for monopile top acceleration
        C[KF.iY['ddq_s'], :] = A[KF.iX['qd_s'], :]
        D[KF.iY['ddq_s'], :] = B[KF.iX['qd_s'], :]
        C[KF.iY['q_p'],  KF.iX['q_p']] = 1   # We measure inclination
        C[KF.iY['dpsi'], KF.iX['dpsi']] = 1  # We measure rotational speed 

        # Nacelle acceleration in x direction (including pitch coupling)
        h_nac = r_TN[2]
        if method == 'YAMS':
            C[KF.iY['NcIMUAx'], :] = A[KF.iX['qd_s'], :] + h_nac * A[KF.iX['qd_p'], :] + A[KF.iX['qd_FA1'], :]
            D[KF.iY['NcIMUAx'], :] = B[KF.iX['qd_s'], :] + h_nac * B[KF.iX['qd_p'], :] + B[KF.iX['qd_FA1'], :]

        C[KF.iY['NcIMUAz'], KF.iX['q_p']] = -9.81
        D[KF.iY['Qgen'], KF.iU['Qgen']] = 1

        # --- Shaping filter, Hydro state equation
        KF.Sw      = 2.3835e-01
        KF.omega_p = 2 * np.pi / Tp
        KF.zeta    = 0.12
        A[KF.iX['q_h'], KF.iX['qd_h']]  = 1
        A[KF.iX['qd_h'], KF.iX['q_h']]  = -KF.omega_p**2
        A[KF.iX['qd_h'], KF.iX['qd_h']] = -2 * KF.zeta * KF.omega_p
        B[KF.iX['qd_h'], KF.iU['w']] = 1 # White noise

        KF.setMat(A, B, C, D)




    def _set_openfast_submatrix(KF, A, B, C, lin_file):
        """Insert measured OpenFAST state and IMU couplings when available."""
        from welib.fast.FASTLin import FASTLin
        FL = FASTLin(linfiles=[lin_file], prefix='', verbose=False)
        data = FL.OP_Data[0].Data[0]
        state_labels = [str(label) for label in FL.xdescr]
        state_map = {
            'q_s': 'PtfmSurge_[m]', 'q_p': 'PtfmPitch_[rad]',
            'q_FA1': 'qt1FA_[m]', 'psi': 'psi_rot_[rad]',
            'qd_s': 'd_PtfmSurge_[m/s]', 'qd_p': 'd_PtfmPitch_[rad/s]',
            'qd_FA1': 'd_qt1FA_[m/s]', 'dpsi': 'd_psi_rot_[rad/s]'}
        indices = {name: state_labels.index(label) for name, label in state_map.items()
                   if label in state_labels}
        for row_name, row in indices.items():
            for col_name, col in indices.items():
                A[KF.iX[row_name], KF.iX[col_name]] = data.A.values[row, col]
        output_map = {'NcIMUAx': 'NcIMUTAxs_[m/s^2]',
                      'NcIMUAy': 'NcIMUTAys_[m/s^2]',
                      'NcIMUAz': 'NcIMUTAzs_[m/s^2]'}
        for output_name, label in output_map.items():
            if label in FL.ydescr:
                row = list(FL.ydescr).index(label)
                for state_name, col in indices.items():
                    C[KF.iY[output_name], KF.iX[state_name]] = data.C.values[row, col]
        if 'Qgen_[Nm]' in FL.udescr:
            col = list(FL.udescr).index('Qgen_[Nm]')
            for state_name, row in indices.items():
                B[KF.iX[state_name], KF.iU['Qgen']] = data.B.values[row, col]
        return FL

    # --- Methods From Parent Class
    # loadMeasurements 
    # prepareTimeStepping 
    # setupCovariances
        
    # --- Methods Common between TN and TNLin
    # prepareMeasurements   

    def timeLoop(KF):
        # --- Aliases to shorten notations
        WT = KF.WT

        # Prepare section output calculation
        dInfo = WT.calcOutputs_init(time=KF.time)
        
        # --- Initial conditions
        x = KF.initFromClean(var='x,y,u')
        
        thrust = KF.U_clean['Thrust'].iloc[0]
        # --- WSE
        ws_last = KF.S_clean['WS'].iloc[0]

        # --- Time loop
        for it in range(0,KF.nt-1):    
            t = it*KF.dt
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values
            # --- Inputs
            u = KF.U_clean.iloc[it,:].values.copy()
            u[KF.iU['Thrust']] = thrust # We use previous estimated thrust as input.
            
            # --- KF predictions
            x, KF.P, _ = KF.estimateTimeStep(u, y, x, KF.P)

            # --- WSE hack
            if KF.hacks['WSE'] == 'clean_inputs':
                qaero_val = KF.X_clean['Qaero'].iloc[it]
                dpsi_val = KF.X_clean['dpsi'].iloc[it]
            else:
                qaero_val = x[KF.iX['Qaero']]
                dpsi_val = x[KF.iX['dpsi']]

            ws, _ = KF.wse.estimate(qaero_val, pitch=u[KF.iU['pitch']] * 180 / np.pi, omega=dpsi_val, WS0=ws_last, relaxation=0, method='oper-crossing', t=KF.time[it])
            ws_last = float(ws)

            # --- Thrust hack
            if KF.hacks['thrust'] == 'clean':
                thrust = KF.U_clean['Thrust'].iloc[it]
            else:
                thrust = float(np.asarray(KF.wse.Thrust(ws_last, pitch=u[KF.iU['pitch']] * 180 / np.pi, omega=dpsi_val)))

            # --- Estimate integrated hydro force (Fx_h)
            q_h  = x[KF.iX['q_h']]
            qd_h = x[KF.iX['qd_h']]
            p_hydro = KF.pHD['phi'] * qd_h
            p_hydro[KF.pST['z'] > 0] = 0.0
            wet_nodes = KF.pST['z'] <= 0
            Fx_h_est = np.trapezoid(p_hydro[wet_nodes], KF.pST['z'][wet_nodes])

            # --- Store extra info
            KF.S_hat.at[it + 1, 'WS']     = ws_last
            KF.S_hat.at[it + 1, 'eta']    = q_h
            KF.S_hat.at[it + 1, 'Fx_h']   = Fx_h_est
            KF.S_hat.at[it + 1, 'Thrust'] = thrust

            # --- Propagation to next time step
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it], ws, thrust))


