"""
Kalman filter model for "Monopile Tower Nacelle Shaft" (based on yams MTNSB)

"""
import os
import numpy as np
from welib.essentials import *
from welib.kalman.kalman import *
from welib.kalman.kalmanfilter import KalmanFilter
from welib.ws_estimator.tabulated import TabulatedWSEstimator
from welib.yams.models.MTNSB import FASTmodel2MTNSB

from welib.weio.fast_linearization_file import FASTLinearizationFile
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions

import welib.fast.fastlib as fastlib
import welib.weio as weio
from welib.kalman.kalman import BuildSystem_Linear_MechOnly
from welib.fast.FASTLin import FASTLin
from welib.tools.strings import FAIL

# --------------------------------------------------------------------------------}
# --- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterMTNS(KalmanFilter):

    def __init__(KF, WSE=None, debug=False, hacks=None):
        # 11 States
        sQ  = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi'] # Mechanical states
        sQa = ['q_h', 'dq_h', 'Qaero']                                         # Augmented states
        sU  = ['Qgen', 'pitch', 'Thrust', 'Fx_i', 'My_i', 'w']
        sY  = ['ddx', 'phi_y', 'dpsi', 'NcIMUAx', 'NcIMUAy', 'NcIMUAz', 'Qgen']
        sS  = ['My_sb', 'Fx_sb', 'eta', 'Fx_h', 'WS', 'Thrust']
        # --- Parent init
        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)
        # --- Storing wind speed estimator (based on tabulated aerodynamic data)
        KF.wse = WSE # wind speed estimator
        KF.debug = debug
        # Hacks
        hacks_def = {'thrust':None, 'WSE':None, 'useTopLoadsFromDF':False}
        if hacks is None:
            KF.hacks = hacks_def
        else:
            hacks_def.update(hacks)
            KF.hacks = hacks_def
        if KF.hacks['thrust']=='clean':
            WARN('HACKING, using thrust from measurements for DEBUG ONLY!')
        if KF.hacks['WSE']=='clean_inputs':
            WARN('HACKING, using WSE inputs from measurements for DEBUG ONLY!')


    def __repr__(self):
        s = KalmanFilter.__repr__(self)
        s+=' - hacks  : {} \n'.format(self.hacks)
        return s

    def setup_matrices(KF, 
                       fstFilename, dfTime=None,
                       hydro_shape_file=None, comp_file=None, Tp=None, zeta =0.12, qdhScale=1, # Hydro params
                       method='YAMS', 
                       lin_file=None, # For method=='OpenFAST
                       ):
        """ Build WT model (sea state, hydro) and state matrices A,B,C,D """
        # --- Default arguments
        shapes_sub =[0,4] # TODO detemine this based on sQ
        shapes_twr =[0]   # TODO detemine this based on sQ
        # --- Windturbine model
        WT = FASTmodel2MTNSB(fstFilename, shapes_sub=shapes_sub, shapes_twr=shapes_twr, shapes_bld=[],
                             DEBUG=False, bStiffening=True, main_axis='z', fixedShaft=False,
                             algo='OpenFAST').WT

		# --- ColMap
        KF.colMap={
                'x'      : 'Q_Sg_[m]' ,
                'phi_y'  : 'Q_P_[rad]' ,
                'q_FA1'  : 'Q_TFA1_[m]',
                'psi'    : '{Azimuth_[deg]} * np.pi/180',
                'dx'     : 'QD_Sg_[m/s]',
                'dphi_y' : 'QD_P_[rad/s]',
                'dq_FA1' : 'QD_TFA1_[m/s]',
                'dpsi'   : '{RotSpeed_[rpm]} * 2*np.pi/60',
                'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                'ddx'    : 'QD2_Sg_[m/s^2]',
                'NcIMUAx': 'NcIMUTAxs_[m/s^2]',
                'NcIMUAy': 'NcIMUTAys_[m/s^2]',
                'NcIMUAz': 'NcIMUTAzs_[m/s^2]',
                'Qgen'   : '{GenTq_[kN-m]} * 1000',
                'pitch'  : '{BldPitch1_[deg]} * np.pi/180',
                'Thrust' : 'RtAeroFxh_[N]',
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
            print('[INFO] Setting Components', comp_file)
            WT.SS_setComponents(comp_file)
            print('[INFO] Setting Compute Eta')
            WT.SS_computeEta(dfTime)
            if hydro_shape_file is not None:
                WT.HD_setShapeFunction(hydro_shape_file)
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
        A = np.zeros((nX, nX))
        B = np.zeros((nX, nU))
        C = np.zeros((nY, nX))
        D = np.zeros((nY, nU))
        # Empty inputs/outputs B,C,D
        MM_sub = WT.MM[:2,:2] # We only keep surge and pitch
        M_inv     = np.linalg.inv(WT.MM)
        M_inv_sub = np.linalg.inv(MM_sub)

        if method=='YAMS':
            As, _, _, _ = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK)

            # --- Structural DOFs
            names = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi']
            for row, row_name in enumerate(names):
                for col, col_name in enumerate(names):
                    A[KF.iX[row_name], KF.iX[col_name]] = As[row, col]

            # Tower fore-aft first mode (1DOF oscillator)
            #omega_t = 2 * np.pi * 0.32
            #A[KF.iX['q_FA1'],  KF.iX['dq_FA1']] = 1
            #A[KF.iX['dq_FA1'], KF.iX['q_FA1']]  = -omega_t**2
            #A[KF.iX['dq_FA1'], KF.iX['dq_FA1']] = -2 * 0.03 * omega_t

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
        try:
            r_TN = WT.r_TN_inT # if hasattr(WT, 'r_TN_inT') else [0,0,198.386]
            r_NS = WT.r_NS_inN # if hasattr(WT, 'r_NS_inN') else [0,0,4.143]
        except:
            raise NotImplementedError()
        h_hub = r_TN[2] + r_NS[2]

        # B matrix columns for Thrust (applicable to both YAMS and OpenFAST)
        B[KF.iX['dx']    , KF.iU['Thrust']] = M_inv[0, 0] * 1.0 + M_inv[0, 1] * h_hub
        B[KF.iX['dphi_y'], KF.iU['Thrust']] = M_inv[1, 0] * 1.0 + M_inv[1, 1] * h_hub

        # For q_FA1, estimate modal mass: 0.25 * m_twr + m_rna
        m_twr = np.trapezoid(WT.twr.m, WT.twr.s_span) if (hasattr(WT, 'twr') and hasattr(WT.twr, 'm')) else 9.6e5
        m_rna = WT.M_RNA if hasattr(WT, 'M_RNA') else 1.189e6
        M_modal = 0.25 * m_twr + m_rna
        B[KF.iX['dq_FA1'], KF.iU['Thrust']] = 1 / M_modal

        # --- Shaft equation
        A[KF.iX['psi'], KF.iX['dpsi']] = 1
        inertia = float(np.asarray(WT.rot.inertia).ravel()[0])
        A[KF.iX['dpsi'], KF.iX['Qaero']] = 1 / inertia
        B[KF.iX['dpsi'], KF.iU['Qgen']] = -1 / inertia

        # --- Generalized hydro force 
        if 'k_h' in pHD:
            A[KF.iX['dx'],     KF.iX['dq_h']] = M_inv_sub[0, :] @ pHD['k_h']
            A[KF.iX['dphi_y'], KF.iX['dq_h']] = M_inv_sub[1, :] @ pHD['k_h']
        else:
            FAIL('k_h not present')

        # --- Output equation for monopile top acceleration
        C[KF.iY['ddx'], :] = A[KF.iX['dx'], :]
        D[KF.iY['ddx'], :] = B[KF.iX['dx'], :]
        C[KF.iY['phi_y'],  KF.iX['phi_y']] = 1   # We measure inclination
        C[KF.iY['dpsi'], KF.iX['dpsi']] = 1  # We measure rotational speed 

        # Nacelle acceleration in x direction (including pitch coupling)
        h_nac = r_TN[2]
        if method == 'YAMS':
            C[KF.iY['NcIMUAx'], :] = A[KF.iX['dx'], :] + h_nac * A[KF.iX['dphi_y'], :] + A[KF.iX['dq_FA1'], :]
            D[KF.iY['NcIMUAx'], :] = B[KF.iX['dx'], :] + h_nac * B[KF.iX['dphi_y'], :] + B[KF.iX['dq_FA1'], :]

        C[KF.iY['NcIMUAz'], KF.iX['phi_y']] = -9.81
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




    def _set_openfast_submatrix(KF, A, B, C, lin_file):
        """Insert measured OpenFAST state and IMU couplings when available."""
        from welib.fast.FASTLin import FASTLin
        FL = FASTLin(linfiles=[lin_file], prefix='', verbose=False)
        data = FL.OP_Data[0].Data[0]
        state_labels = [str(label) for label in FL.xdescr]
        state_map = {
            'x': 'PtfmSurge_[m]', 'phi_y': 'PtfmPitch_[rad]',
            'q_FA1': 'qt1FA_[m]', 'psi': 'psi_rot_[rad]',
            'dx': 'd_PtfmSurge_[m/s]', 'dphi_y': 'd_PtfmPitch_[rad/s]',
            'dq_FA1': 'd_qt1FA_[m/s]', 'dpsi': 'd_psi_rot_[rad/s]'}
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
        for it in range(0, KF.nt-1):    
            t = it * KF.dt
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values
            # --- Inputs
            u = KF.U_clean.iloc[it,:].values.copy()
            u[KF.iU['Thrust']] = thrust # We use previous estimated thrust as input.
            
            # --- Predictions of next time step based on current time step
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

            # --- Estimate Generalized hydro force and bending moment (calc output)
            q_h     = x[KF.iX['q_h']]  # eta
            dq_h    = x[KF.iX['dq_h']] # eta_dot
            eta     = q_h
            eta_dot = dq_h
            p_hydro = KF.pHD['phi'] * eta_dot # p_h = k_h(z) q_h(t)
            p_hydro[KF.zDepth>0] = 0 # safety, shoudn't be necessary

            # --- Acceleration
            x_dot = np.dot(KF.A, x) + np.dot(KF.B, u)
            
            # ---
            p_ext      = np.zeros((3,len(KF.zDepth)))
            p_ext[0,:] = p_hydro

            x_q   = np.array([x[0],x[1]])
            xd_q  = np.array([x_dot[0], x_dot[1]])
            xdd_q = np.array([x_dot[2], x_dot[3]])

            # --- DOFs in the way expected by calcOutputs_step
            q   = dInfo['q_default'].copy()
            qd  = dInfo['q_default'].copy()
            qdd = dInfo['q_default'].copy()
            q  ['Sg']  = x_q[0]
            q  ['P']   = x_q[1]
            qd ['Sg']  = xd_q[0]
            qd ['P']   = xd_q[1]
            qdd['Sg']  = xdd_q[0]
            qdd['P']   = xdd_q[1]

            # ---  Top loads
            F_top = np.array((0.,0.,0.))
            M_top = np.array((0.,0.,0.))
            a_ext = np.array((0.,0.,-WT.gravity)) # external acceleration (gravity/earthquake)

            Thrust = thrust
            ser_Loads = pd.Series({'Fadd_R_xs':Thrust, 'Fadd_R_ys':0, 'Fadd_R_zs':0, 'Madd_R_xs':0, 'Madd_R_ys':0, 'Madd_R_zs':0})
            # TEMPORARY for backward compatibility
#             ser_Loads['TwrBsFxt_[kN]'] = 0
#             ser_Loads['TwrBsFyt_[kN]'] = 0
#             ser_Loads['TwrBsFzt_[kN]'] = 0
#             ser_Loads['TwrBsMxt_[kN-m]'] = 0
#             ser_Loads['TwrBsMyt_[kN-m]'] = 0
#             ser_Loads['TwrBsMzt_[kN-m]'] = 0
#             dInfo['useInterfaceLoadsFromDF'] = True

            rowOut, twr_it, mnp_it = WT.calcOutputs_step(q, qd, qdd, dInfo, t=KF.time[it], ser_Loads=ser_Loads, mnp_p_ext=p_ext)
            F_sec = mnp_it[0:3,:]
            M_sec = mnp_it[3:6,:]

            #F_sec, M_sec, outD = beamSectionLoadsFromShapeFunctions    (x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext=a_ext, corrections=0, PhiK=WT.fnd.PhiK)

            # No acceleration # TODO get it from calcOutputs 
            xdd_q *=0
            F_sec_h, M_sec_h, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext = a_ext, PhiK=WT.fnd.PhiK)

            wet_nodes = KF.zDepth <= 0
            Fx_h_est = np.trapezoid(p_hydro[wet_nodes], KF.zDepth[wet_nodes])
            
            
            # --- Store extra info
            # Environment            
            KF.S_hat.at[it+1, 'WS']     = ws_last
            KF.S_hat.at[it+1, 'eta']    = q_h
            # Loads            
            KF.S_hat.at[it+1, 'My_sb']  = M_sec[1,0]
            KF.S_hat.at[it+1, 'Fx_sb']  = F_sec[0,0]
            KF.S_hat.at[it+1, 'Fx_h']   = Fx_h_est
            KF.S_hat.at[it+1, 'Thrust'] = thrust

            # --- Propagation to next time step
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it], ws, thrust))


