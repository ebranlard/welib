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
from welib.yams.models.TNSB_FAST import FASTmodel2TNSB
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions


# --------------------------------------------------------------------------------}
# --- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterMTNS(KalmanFilter):

    def __init__(KF, WSE=None, debug=False, hacks=None, setup_opts=None):
        """

        """    
        # --- Default arguments
        setup_def = {'hydro_states':True, 'monopileDOFs':True}
        hacks_def = {'thrust':None, 'WSE':None, 'SL_cleanQ':False, 'SL_cleanFtop':False, 'SL_cleanEtaDot':False, 'SL_cleanP':False}
        if setup_opts is None:
            KF.setup_opts = setup_def
        else:
            setup_def.update(setup_opts)
            KF.setup_opts = setup_def
        if hacks is None:
            KF.hacks = hacks_def
        else:
            hacks_def.update(hacks)
            KF.hacks = hacks_def


        # --- Initialize Kalman Filter, variables names (e.g. sX) and matrices (Xx=A)
        sQ  = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi'] # Mechanical states
        sQa = ['q_h', 'dq_h', 'Qaero']              # Augmented states
        sU  = ['w', 'Qgen', 'pitch', 'GenThrust'] #, 'Fx_i', 'My_i']
        sY  = ['PtfmIMUAx', 'PtfmIncly',  'NcIMUAx', 'dpsi',  'Qgen']
        sS  = ['Fx_sb', 'My_sb', 'eta', 'Fx_h', 'Fx_i', 'My_i', 'WS', 'GenThrust', 'Thrust']
        if not KF.setup_opts['hydro_states']:
            sQa.remove('q_h'); sQa.remove('dq_h'); 
            sU.remove('w'); 
            sS.remove('eta'); sS.remove('Fx_h'); 
        if not KF.setup_opts['monopileDOFs']:
            sQ.remove('x'); sQ.remove('phi_y'); sQ.remove('dx'); sQ.remove('dphi_y'); 
            sY.remove('PtfmIMUAx'); sY.remove('PtfmIncly'); 
        if not KF.setup_opts['aero_est']:
            sQ.remove( 'psi'); sQ.remove( 'dpsi'); 
            sU.remove( 'Qgen'); sU.remove( 'pitch'); sU.remove('GenThrust'); 
            sY.remove( 'dpsi'); sY.remove( 'Qgen');
            sS.remove( 'WS'); sS.remove( 'GenThrust'); sS.remove('Thrust')
            sQa.remove('Qaero');
        if not KF.setup_opts['q_FA1']:
            sQ.remove( 'q_FA1'); sQ.remove( 'dq_FA1'); 
            sY.remove('NcIMUAx')
        # Constant input (value 1) carrying the static generalized forces (gravity overhang) of the full structure
        if KF.setup_opts['aero_est'] and KF.setup_opts['q_FA1']:
            sU.append('Static')
        if KF.setup_opts['monopileDOFs'] and KF.setup_opts['aero_est'] and KF.setup_opts['q_FA1']:
            if KF.setup_opts.get('bias_states', True):
                sQa += ['Bx', 'Bphi'] # Slowly varying generalized force biases (model static errors) on x and phi_y
        elif not KF.setup_opts['monopileDOFs'] and KF.setup_opts['q_FA1'] and KF.setup_opts['aero_est']:
            # Tower-only case (Case 2): bias on q_FA1 (slowly-varying model error on tower mode)
            if KF.setup_opts.get('bias_states', True):
                sQa += ['Bq']  # Tower mode bias

        # --- Parent init
        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)
        # --- Storing wind speed estimator (based on tabulated aerodynamic data)
        KF.wse = WSE # wind speed estimator
        KF.debug = debug
        # Hacks
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

    # --------------------------------------------------------------------------------}
    # --- Setup of the models
    # --------------------------------------------------------------------------------{
    # Two models are always built, and stored in KF.mats:
    #   - 'YS': YAMS model, physics from the structural matrices (M,D,K) and the geometry of the turbine
    #   - 'OF': OpenFAST model, built from the (averaged) linearized OpenFAST model
    # They share the same states/inputs/outputs and the same code for what is not method specific (_set_common_submatrix)
    # `method` only selects which one is used by the filter. See compareMatrices to compare them element by element.
    def setup_matrices(KF, fstFile, time_ref=None, hydroShapeFile=None, compFile=None, eta_ref=None, Tp=None, zeta=0.12,
                       method='YAMS', linFiles=None):
        """ Build WT model (sea state, hydro) and state matrices A,B,C,D """
        if linFiles is None:
            raise FileNotFoundError('An OpenFAST .lin file is required for now')
        elif not isinstance(linFiles, list):
            linFiles = [linFiles]
        # --- Models, maps to measurements, hydro
        KF._setup_WT(fstFile)
        KF._setup_colMap()
        KF._setup_hydro(time_ref, eta_ref, compFile, hydroShapeFile, Tp, zeta)
        KF.M_inv = np.linalg.inv(KF.WT.MM)
        OF_lin   = KF._load_openfast_lin(linFiles)
        # --- Matrices: initialize zero matrices for YAMS (YS) and OpenFAST (OF) models
        # Both models share the same states/inputs/outputs; only the physics differs
        nX, nU, nY = KF.nX, KF.nU, KF.nY
        YS = [np.zeros((nX, nX)), np.zeros((nX, nU)), np.zeros((nY, nX)), np.zeros((nY, nU))]
        OF = [np.zeros((nX, nX)), np.zeros((nX, nU)), np.zeros((nY, nX)), np.zeros((nY, nU))]
        # Rotor inertia: YAMS from the rotor model, OpenFAST from the Qgen column of the lin. model
        J_LSS_YS = KF.WT.rot.inertia[0,0] if 'psi' in KF.iX else None
        J_LSS_OF = -1/OF_lin['B'].loc['d_psi_rot_[rad/s]','Qgen_[Nm]'] if 'psi' in KF.iX else None
        if 'psi' in KF.iX:
            print('[INFO] KalmanModel: Rotor Inertia: YAMS: {:.1f} OF: {:.1f}'.format(J_LSS_YS, J_LSS_OF))

        # --- Fill matrix blocks: common terms (identical for both models)
        # Sets: rotor equation (psi,dpsi), wave shaping filter (q_h,dq_h), simple outputs (PtfmIncly, dpsi, Qgen)
        KF._set_common_submatrix(*YS, J_LSS_YS)
        KF._set_common_submatrix(*OF, J_LSS_OF)

        # --- Fill matrix blocks: model-specific terms
        # YAMS uses mechanical matrices (M,D,K) from YAMS; OpenFAST uses lin. model from OpenFAST
        KF._set_yams_submatrix(*YS, method)
        KF._set_openfast_submatrix(*OF, OF_lin)

        # --- Hybrid model: YAMS with selected blocks replaced by OpenFAST (if user wants it)
        HY = KF._set_hybrid_submatrix(YS, OF, OF_lin)

        # --- Selection of the model used by the filter: YAMS, OpenFAST, or Hybrid
        KF.method = method
        KF.mats   = dict(YS=tuple(YS), OF=tuple(OF), HY=tuple(HY), OF_lin=OF_lin, M_inv=KF.M_inv)
        KF.setMat(*dict(YAMS=YS, OpenFAST=OF, Hybrid=HY)[method])

    def _setup_WT(KF, fstFile):
        """ Windturbine models: full structure (WT, monopile and tower) and tower only (WTTN)"""
        # TODO shapes should be determined based on sQ
        KF.WT = FASTmodel2MTNSB(fstFile, shapes_sub=[0,4], shapes_twr=[0], shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z', fixedShaft=False, algo='OpenFAST').WT
        KF.WTTN = FASTmodel2TNSB(fstFile, shapes_twr=[0], shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z').WT
        KF.zDepth = KF.WT.fnd.s_span - KF.WT.WtrDpth

    def _setup_colMap(KF):
        """ Col MAP for OpenFAST OutFile "Measurements" used for "clean" values """
        WT = KF.WT
        nGear = WT.ED['GBRatio']
        KF.colMap={
                'x'      : 'Q_Sg_[m]' ,
                'phi_y'  : 'Q_P_[rad]' ,
                'psi'    : '{Azimuth_[deg]} * np.pi/180', # SI [deg] -> [rad]
                'dx'     : 'QD_Sg_[m/s]',
                'dphi_y' : 'QD_P_[rad/s]',
                'dpsi'   : '{RotSpeed_[rpm]} * 2*np.pi/60', # SI [rpm] -> [rad/s]
                'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                'ddx'    : 'QD2_Sg_[m/s^2]',
                'ddphi_y': 'QD2_P_[rad/s^2]',
                'PtfmIMUAx': '{QD2_Sg_[m/s^2]}',
                'PtfmIncly': '{Q_P_[rad]}',
                'NcIMUAx': 'NcIMUTAxs_[m/s^2]',
                'eta'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                #'NcIMUAy': 'NcIMUTAys_[m/s^2]',
                #'NcIMUAz': 'NcIMUTAzs_[m/s^2]',
                'pitch'  : '{BldPitch1_[deg]} * np.pi/180', # SI [deg]->[rad]
                'Fx_h'   : '{HydroFxi_[N]}',   # We use a trick  
                'Fx_sb'  : '-ReactFXss_[N]',
                'My_sb'  : '-ReactMYss_[N*m]',

        }
        if WT.fnd is not None:
            KF.colMap.update({
                'Fx_i'   : 'IntfFXss_[N]',   
                'My_i'   : 'IntfMYss_[N*m]',
                })
        else:
            KF.colMap.update({
                'Fx_i'   : '{TwrBsFxt_[kN]}*1000',   
                'My_i'   : '{TwrBsMyt_[kN-m]}*1000',
                })
        if 'psi' in KF.sX:
            KF.colMap.update({
                'Qgen'   : f'{nGear}'+'*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]  # NOTE: nGear
                'ddpsi'  : 'QD2_GeAz_[rad/s^2]',
                })
        if 'GenThrust' in KF.sX or 'GenThrust' in KF.sU:
            KF.colMap.update({
                'Thrust' : '{RtAeroFxh_[N]}',
                'GenThrust':'{RtAeroFxh_[N]}', # TODO need a proper map
                'Qaero'  : 'RtAeroMxh_[N-m]',
                'WS'     : 'RtVAvgxh_[m/s]',
                })
        if 'Static' in KF.sU:
            KF.colMap['Static'] = '{Time_[s]}*0 + 1'
        for nm in ('Bx', 'Bphi'):
            if nm in KF.sX:
                KF.colMap[nm] = '{Time_[s]}*0'
        if 'q_FA1' in KF.sX:
            KF.colMap.update({
                'q_FA1'  : 'Q_TFA1_[m]',
                'dq_FA1' : 'QD_TFA1_[m/s]',
                'ddq_FA1': 'QD2_TFA1_[m/s^2]',
                })

    def _setup_hydro(KF, time_ref, eta_ref, compFile, hydroShapeFile, Tp, zeta):
        """ Sea state components, wave elevation, hydro shape function, and shaping filter of the waves"""
        WT = KF.WT
        if 'q_h' not in KF.sX:
            return
        # --- Configure Sea State components and wave elevation
        if 'q_h' in KF.sX:
            if WT.pSS is not None:
                if compFile is not None:
                    print('[INFO] Setting Components from file:', compFile)
                    WT.SS_setComponents(compFile)
                elif eta_ref is not None:
                    NOTE('Setting Components from eta ref')
                    WT.SS_setComponents(time=time_ref, eta=eta_ref)
                print('[INFO] Compute Eta')
                WT.SS_computeEta(time_ref)
                if hydroShapeFile is not None:
                    WT.HD_setShapeFunction(hydroShapeFile)
            # Ensure WT.MM contains the hydrodynamic mass if not already added
            pHD = WT.pHD
            KF.pHD = pHD
            if not getattr(WT, '_GM_hydro_added', False) and 'GM_hydro' in pHD:
                for i in range(min(pHD['GM_hydro'].shape[0], WT.MM.shape[0])):
                    WT.MM[i, i] += pHD['GM_hydro'][i, i]
                WT._GM_hydro_added = True

        # Shaping filter parameters, state equation of the waves is in _set_common_submatrix
        if Tp==12.7:
            KF.Sw= 2.3835e-01
        elif Tp==10.0:
            KF.Sw= 2.3835e-01/2
        else:
            raise NotImplementedError(f'Tp={Tp}')
        KF.omega_p = 2*np.pi/Tp
        KF.zeta    = zeta
        print('omega_p^2', KF.omega_p**2, '2 zeta omega_p', 2*KF.zeta*KF.omega_p)

    def _k_h(KF):
        """ Generalized hydro force: k_h * dq_h, for x and phi_y (including user scaling)"""
        if 'k_h' not in KF.pHD:
            FAIL('k_h not present')
        return np.asarray(KF.pHD['k_h']) * np.asarray(KF.setup_opts.get('kh_scale', [1.0, 1.0]))

    def _set_nacelle_IMU(KF, A, B, C, D, w):
        """ Nacelle IMU acceleration as a combination of the generalized accelerations (rows of A and B), weights w for (x, phi_y, q_FA1) """
        if 'NcIMUAx' not in KF.iY: return
        iy = KF.iY['NcIMUAx']; C[iy,:] = 0; D[iy,:] = 0
        for wk, nm in zip(w, ['dx', 'dphi_y', 'dq_FA1']):
            if nm in KF.iX:
                C[iy,:] += wk * A[KF.iX[nm],:]
                D[iy,:] += wk * B[KF.iX[nm],:]

    # --------------------------------------------------------------------------------}
    # --- Common to YAMS and OpenFAST
    # --------------------------------------------------------------------------------{
    def _set_common_submatrix(KF, A, B, C, D, J_LSS):
        """ Augmented states equations and simple outputs, identical for both models. J_LSS: rotor inertia"""
        iX, iU, iY = KF.iX, KF.iU, KF.iY
        # --- Rotor Inertia / Shaft equation
        if 'psi' in iX:
            A[iX['psi'], iX['dpsi']] = 1
            if 'Qaero' in KF.sXa: A[iX['dpsi'], iX['Qaero']] =  1 / J_LSS
            if 'Qgen'  in iU    : B[iX['dpsi'], iU['Qgen']]  = -1 / J_LSS
            C[iY['dpsi'], iX['dpsi']] = 1  # We measure rotational speed
        # --- Outputs
        if 'x' in iX:
            C[iY['PtfmIncly'], iX['phi_y']] = 1   # We measure inclination
        if 'Qgen' in iY:
            D[iY['Qgen'], iU['Qgen']] = 1
        # --- Shaping filter, Hydro state equation
        if 'q_h' in iX:
            A[iX['q_h'],  iX['dq_h']]  = 1
            A[iX['dq_h'], iX['q_h']]   = -KF.omega_p**2
            A[iX['dq_h'], iX['dq_h']]  = -2 * KF.zeta * KF.omega_p
            B[iX['dq_h'], iU['w']]     = 1 # White noise

    # --------------------------------------------------------------------------------}
    # --- YAMS
    # --------------------------------------------------------------------------------{
    def _set_yams_submatrix(KF, A, B, C, D, method='YAMS'):
        """ Rows of the generalized accelerations from YAMS: mechanical matrices, thrust, static force, biases, hydro, and IMU outputs"""
        iX, iU, iY = KF.iX, KF.iU, KF.iY
        WT, M_inv = KF.WT, KF.M_inv
        dofs = [('dx',0), ('dphi_y',1), ('dq_FA1',2)] # Generalized accelerations and their index in M
        # --- Mechanical system
        As, _, _, _ = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK)
        names = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi']
        for row, row_name in enumerate(names):
            for col, col_name in enumerate(names):
                if row_name in iX and col_name in iX:
                    A[iX[row_name], iX[col_name]] = As[row, col]
        # --- Thrust and static force
        if 'GenThrust' in iU and 'x' in iX:
            # Input is the physical rotor thrust T along the shaft, applied at the rotor center R.
            # Generalized forces: Q = J*T + b, with J the Jacobian of the point R w.r.t. [x, phi_y, q_FA1]
            W = KF.WTTN
            r_TN = WT.r_TN_inT
            cs, sn = np.cos(W.shaft_tilt), np.sin(W.shaft_tilt)
            rNR  = np.asarray(W.r_NR_inN).flatten()
            rNG  = np.asarray(W.r_NGrna_inN).flatten()
            vy1c = W.twr.Bhat_t_bc[1,0]
            J = np.zeros(M_inv.shape[0]); bs = np.zeros(M_inv.shape[0])
            J[0] = cs
            J[1] = (r_TN[2] + rNR[2])*cs + rNR[0]*sn
            J[2] = cs + vy1c*(rNR[0]*sn + rNR[2]*cs)
            bs[1] = rNG[0]*W.M_RNA*W.gravity           # overhang moment of the RNA weight
            bs[2] = vy1c*rNG[0]*W.M_RNA*W.gravity
            sb = np.asarray(KF.setup_opts.get('static_bias', [0,0,0])); bs[:len(sb)] += sb
            KF.J_thrust, KF.b_static = J, bs
            for nm, i in dofs:
                if nm in iX:
                    B[iX[nm], iU['GenThrust']] = M_inv[i,:] @ J
                    if 'Static' in iU: B[iX[nm], iU['Static']] = M_inv[i,:] @ bs
        elif 'GenThrust' in iU:
            # Tower only: thrust at the top of the tower, generalized mass of the tower+RNA
            B[iX['dq_FA1'], iU['GenThrust']] = 1 / (WT.twr.MM[6,6] + WT.RNA.mass)
        # --- Generalized force biases
        for nmB, j in (('Bx',0), ('Bphi',1)):
            if nmB in iX:
                for nm, i in dofs:
                    if nm in iX: A[iX[nm], iX[nmB]] = M_inv[i, j]
        if 'Bq' in iX:  # Tower-only bias (Case 2)
            A[iX['dq_FA1'], iX['Bq']] = 1.0 / (WT.twr.MM[6,6] + WT.RNA.mass) if 'GenThrust' in iU else 1.0
        # --- Generalized hydro force
        if 'q_h' in iX:
            k_h = KF._k_h()
            for nm, i in dofs:
                if nm in iX: A[iX[nm], iX['dq_h']] = M_inv[i, :2] @ k_h
        # --- HACK for tower only: use the tower-only model (K/M) instead of the coupled model of the full structure
        # TODO TODO TODO the full structure model is not valid when monopile DOFs are removed (condensation needed)
        if 'q_FA1' in iX and 'x' not in iX:
            if method=='YAMS': FAIL('Using a super hack: YAMS tower-only -K/M instead of coupled model')
            A[iX['dq_FA1'], iX['q_FA1']] = -KF.WTTN.KK[0,0]/KF.WTTN.MM[0,0] # -1.71 instead of -55.16
        # --- Outputs
        if 'x' in iX:
            C[iY['PtfmIMUAx'], :] = A[iX['dx'],:] # PtfmIMUAx is assumed to be qdd_s
            D[iY['PtfmIMUAx'], :] = B[iX['dx'],:]
        KF._set_nacelle_IMU(A, B, C, D, (1, WT.r_TN_inT[2], 1))

    # --------------------------------------------------------------------------------}
    # --- OpenFAST
    # --------------------------------------------------------------------------------{
    def _of_inverse_mass(KF, Bm, dofs):
        """ Inverse mass matrix of the DOFs `dofs` (x, phi_y, q_FA1), deduced from the input columns of the lin. model.
        Generalized forces (x: PtfmFx, phi_y: PtfmMy) give the x and phi_y columns. The q_FA1 column comes from a nacelle force NacFx,
        applied at height h with tower-top mode shape phi: Minv @ [1, h, phi]; (h, phi) are solved from the x and phi_y rows."""
        if dofs != ['x', 'phi_y', 'q_FA1']: raise NotImplementedError('Inverse mass matrix for DOFs: '+str(dofs))
        rows = ['d_PtfmSurge_[m/s]', 'd_PtfmPitch_[rad/s]', 'd_qt1FA_[m/s]']
        Mi = np.zeros((3,3))
        Mi[:,0] = Bm.loc[rows, 'PtfmFxN1_[N]'].values;  Mi[:,1] = Bm.loc[rows, 'PtfmMyN1_[Nm]'].values
        Mi[:,2] = Mi[2,:]; Mi[2,2] = 0 # symmetry
        b = Bm.loc[rows, 'NacFxN1_[N]'].values
        # b[:2] = Mi[:2,0] + h Mi[:2,1] + phi Mi[:2,2]
        h, phi = np.linalg.solve(np.column_stack([Mi[:2,1], Mi[:2,2]]), b[:2] - Mi[:2,0])
        Mi[2,2] = (b[2] - Mi[2,0] - h*Mi[2,1]) / phi
        print('[INFO] OpenFAST lin. inverse mass: nacelle height {:.2f}, tower-top mode shape {:.3f}, M_qq {:.4e}'.format(h, phi, np.linalg.inv(Mi)[2,2]))
        return Mi

    def _set_openfast_submatrix(KF, A, B, C, D, OF_lin):
        """ Rows of the generalized accelerations from the linearized OpenFAST model (Am, Bm, Cm, Dm) for the DOFs present in the filter.
        Structure (A), thrust (HubFx/HubFz columns), generalized force biases and hydro force (PtfmFx/PtfmMy columns),
        IMU outputs (Cm/Dm rows), and static offset (operating point).  Heave and rotor-azimuth couplings are dropped.
        When the filter has fewer DOFs than the lin. model (e.g. tower only), the model is condensed: Mss^-1 (M a)_s """
        iX, iU, iY = KF.iX, KF.iU, KF.iY
        Am, Bm, Cm, Dm, FL = OF_lin['A'], OF_lin['B'], OF_lin['C'], OF_lin['D'], OF_lin['FL']
        opts = KF.setup_opts.get('of_opts', {})
        cs, sn = np.cos(KF.WTTN.shaft_tilt), np.sin(KF.WTTN.shaft_tilt)
        lq = {'x':'PtfmSurge_[m]', 'phi_y':'PtfmPitch_[rad]', 'q_FA1':'qt1FA_[m]'}
        lv = {'x':'d_PtfmSurge_[m/s]', 'phi_y':'d_PtfmPitch_[rad/s]', 'q_FA1':'d_qt1FA_[m/s]'}
        full = [k for k in lq if lq[k] in Am.columns]  # DOFs of the lin. model
        sub  = [k for k in full if k in iX]            # DOFs of the filter
        if len(sub)==0: return
        idx  = [full.index(k) for k in sub]
        sStates = [lq[k] for k in sub] + [lv[k] for k in sub]    # states of the filter (lin. labels)
        iStates = [iX[k] for k in sub] + [iX['d'+k] for k in sub]
        rowsF = [lv[k] for k in full]  # acceleration rows of the lin. model
        # --- Condensation, project(col): column of accelerations of the lin. model -> accelerations of the filter DOFs
        if sub == full:
            project = lambda col: np.asarray(col)
        else:
            M = np.linalg.inv(KF._of_inverse_mass(Bm, full))
            Mss_inv = np.linalg.inv(M[np.ix_(idx, idx)])
            project = lambda col: Mss_inv @ (M @ np.asarray(col))[idx]
        def colB(label):   # Input column of the lin. model, projected on the filter DOFs
            return project(Bm.loc[rowsF, label].values)
        def colT(M, rows): # thrust along shaft: [cos, 0, -sin] in global frame
            return cs * M.loc[rows, 'HubFxN1_[N]'].values - sn * M.loc[rows, 'HubFzN1_[N]'].values
        W = KF.WTTN; rNG = np.asarray(W.r_NGrna_inN).flatten(); Mg = W.M_RNA * W.gravity  # RNA weight, for static_mode 'gravity'
        # Operating point (mean over the OPs)
        xdescr = [str(l) for l in FL.xdescr]; xop = np.asarray(FL.xop_mean).ravel()
        xop_s = np.array([xop[xdescr.index(l)] for l in sStates])
        ydescr = [str(l) for l in FL.ydescr]
        yop_all = np.mean([np.asarray(OP.Data[0].y).ravel() for OP in FL.OP_Data], axis=0)
        yop = lambda lab: yop_all[ydescr.index(lab)]
        T_lab = opts.get('T_op_label', 'ADRtAeroFxh_[N]')
        T_op = yop(T_lab) * (1e3 if 'kN' in T_lab else 1.0) if T_lab in ydescr else 0.0  # aero thrust (no rotor weight/inertia)
        # --- Structural block: kinematic rows, and accelerations (columns of the filter states only)
        for k in sub:
            A[iX[k], iX['d'+k]] = 1
        Aacc = np.column_stack([project(Am.loc[rowsF, l].values) for l in sStates])
        for i, k in enumerate(sub):
            A[iX['d'+k], iStates] = Aacc[i, :]
        if 'q_FA1' in sub and opts.get('aero_dq_scale', 1.0) != 1.0:   # damping of 1st FA mode (contains aero damping of rotor)
            A[iX['dq_FA1'], iX['dq_FA1']] *= opts['aero_dq_scale']
        iRows = [iX['d'+k] for k in sub]
        # --- Thrust
        if 'GenThrust' in iU:
            B[iRows, iU['GenThrust']] = project(colT(Bm, rowsF))
        bT = B[iRows, iU['GenThrust']] if 'GenThrust' in iU else 0
        # --- Static offset: f = A (x - xop) + B_T (T - T_op)
        static_mode = opts.get('static', 'op')
        if 'Static' in iU and static_mode == 'op':
            B[iRows, iU['Static']] = -A[np.ix_(iRows, iStates)] @ xop_s - bT * T_op
        elif 'Static' in iU and static_mode == 'gravity':
            # RNA weight (-M g at the RNA CoM) applied through the nacelle Fz/My input columns of the lin model
            B[iRows, iU['Static']] = -Mg * colB('NacFzN1_[N]') + rNG[0] * Mg * colB('NacMyN1_[Nm]')
        # --- Generalized force biases (closed-loop inverse-mass columns)
        for nmB, lab in (('Bx','PtfmFxN1_[N]'), ('Bphi','PtfmMyN1_[Nm]')):
            if nmB in iX:
                A[iRows, iX[nmB]] = colB(lab)
        if 'Bq' in iX:  # Tower-only bias (Case 2): use constant factor
            A[iX['dq_FA1'], iX['Bq']] = 1.0 / (KF.WTTN.twr.MM[6,6] + KF.WTTN.RNA.mass) if iX.get('dq_FA1') else 0

        # --- Hydro generalized force k_h * dq_h
        if 'q_h' in iX:
            k_h = KF._k_h()
            A[iRows, iX['dq_h']] = k_h[0]*colB('PtfmFxN1_[N]') + k_h[1]*colB('PtfmMyN1_[Nm]')
        # --- Outputs from Cm, Dm rows
        output_map = {'PtfmIMUAx':'QD2_Sg_[m/s^2]', 'PtfmIncly':'Q_P_[rad]'}
        if 'NcIMUAx' in iY:
            # Nacelle IMU acceleration written as combination of the generalized accelerations (QD2_Sg, QD2_P, QD2_TFA1):
            # weights from least-squares fit of the NcIMUTAxs row of (Cm,Dm) on the QD2_* rows (all DOFs of the lin. model).
            # Outputs are then consistent with the state equations (rows of A, B), including offsets.
            acc = {'x':'QD2_Sg_[m/s^2]', 'phi_y':'QD2_P_[rad/s^2]', 'q_FA1':'QD2_TFA1_[m/s^2]'}
            sF  = [lq[k] for k in full] + [lv[k] for k in full]
            cB  = ['HubFxN1_[N]', 'HubFzN1_[N]', 'PtfmFxN1_[N]', 'PtfmMyN1_[Nm]']
            def row(l): return np.concatenate([Cm.loc[l, sF].values, Dm.loc[l, cB].values])
            M_ = np.array([row(acc[k]) for k in full]).T; y = row('NcIMUTAxs_[m/s^2]'); sc = np.abs(M_).max(axis=1) + 1e-30
            w = np.linalg.lstsq(M_/sc[:, None], y/sc, rcond=None)[0]
            OF_lin['w_nc'] = w; OF_lin['full'] = full
            KF._set_nacelle_IMU(A, B, C, D, [dict(zip(full, w)).get(k, 0) for k in ('x', 'phi_y', 'q_FA1')])
        for on, lab in output_map.items():
            if on not in iY or lab not in ydescr: continue
            iy = iY[on]
            C[iy, :] = 0; D[iy, :] = 0
            C[iy, iStates] = Cm.loc[lab, sStates].values
            dT = colT(Dm, [lab])[0] if 'GenThrust' in iU else 0
            if 'GenThrust' in iU:
                D[iy, iU['GenThrust']] = dT
            if 'Static' in iU and static_mode == 'op':
                D[iy, iU['Static']] = yop(lab) - Cm.loc[lab, sStates].values @ xop_s - dT * T_op
            if 'Static' in iU and static_mode == 'gravity':
                D[iy, iU['Static']] = -Mg * Dm.loc[lab, 'NacFzN1_[N]'] + rNG[0] * Mg * Dm.loc[lab, 'NacMyN1_[Nm]']
            for nmB, labB in (('Bx','PtfmFxN1_[N]'), ('Bphi','PtfmMyN1_[Nm]')):
                if nmB in iX:
                    C[iy, iX[nmB]] = Dm.loc[lab, labB]
            if 'q_h' in iX:
                k_h = KF._k_h()
                C[iy, iX['dq_h']] = Dm.loc[lab, 'PtfmFxN1_[N]']*k_h[0] + Dm.loc[lab, 'PtfmMyN1_[Nm]']*k_h[1]
        OF_lin['xop_s'] = dict(zip(sStates, xop_s)); OF_lin['T_op'] = T_op

    def _set_hybrid_submatrix(KF, YS, OF, OF_lin):
        """ YAMS model in which the terms listed in setup_opts['hybrid'] are replaced by the OpenFAST ones.
        Hybrid model: best performance on Case 3 (full structure). Terms (columns/blocks of the acceleration rows):
        'KD': stiffness/damping (A, structural states), 'thrust': GenThrust column, 'static': Static column,
        'bias': Bx,Bphi columns, 'hydro': dq_h column.
        IMU outputs are recomputed from the hybrid rows (OpenFAST weights for the nacelle IMU if any term is replaced)."""
        iX, iU, iY = KF.iX, KF.iU, KF.iY
        HY = [m.copy() for m in YS]
        # Default hybrid terms: ['KD','static','hydro'] gives best performance on Case 3 across 16 cases
        # User can override via setup_opts['hybrid'] with any subset of: ['KD','thrust','static','bias','hydro']
        HYBRID_DEFAULTS = ['KD','static','hydro']
        terms = KF.setup_opts.get('hybrid', HYBRID_DEFAULTS)
        acc = [iX[k] for k in ('dx','dphi_y','dq_FA1') if k in iX]
        struct = [iX[k] for k in ('x','phi_y','q_FA1','dx','dphi_y','dq_FA1') if k in iX]
        def rep(im, rows, cols):  # Helper: copy OF block into HY[im]
            if len(cols)>0: HY[im][np.ix_(rows, cols)] = OF[im][np.ix_(rows, cols)]
        # Replace selected blocks with OpenFAST values
        if 'KD'     in terms: rep(0, acc, struct)
        if 'thrust' in terms and 'GenThrust' in iU: rep(1, acc, [iU['GenThrust']])
        if 'static' in terms and 'Static'    in iU: rep(1, acc, [iU['Static']])
        if 'bias'   in terms: rep(0, acc, [iX[k] for k in ('Bx','Bphi') if k in iX])
        if 'hydro'  in terms and 'q_h' in iX: rep(0, acc, [iX['dq_h']])
        if terms and 'w_nc' in OF_lin:
            A, B, C, D = HY
            if 'x' in iX:
                C[iY['PtfmIMUAx'], :] = A[iX['dx'],:]; D[iY['PtfmIMUAx'], :] = B[iX['dx'],:]
            wd = dict(zip(OF_lin['full'], OF_lin['w_nc']))
            KF._set_nacelle_IMU(A, B, C, D, [wd.get(k, 0) for k in ('x','phi_y','q_FA1')])
        return HY

    def compareMatrices(KF, tol=1e-12):
        """ Prints the non-zero elements of the YAMS and OpenFAST matrices side by side (with state/input/output names).
        Used for debugging with pdb or to see which YAMS terms differ from OpenFAST """
        names = dict(A=(KF.sX, KF.sX), B=(KF.sX, KF.sU), C=(KF.sY, KF.sX), D=(KF.sY, KF.sU))
        for im, nm in enumerate('ABCD'):
            Y, O = KF.mats['YS'][im], KF.mats['OF'][im]
            rows, cols = names[nm]
            print('--- {:s} matrix:  {:>14s} {:>14s} {:>14s} {:>10s}'.format(nm, 'YAMS', 'OpenFAST', 'diff', 'ratio'))
            for i in range(Y.shape[0]):
                for j in range(Y.shape[1]):
                    if abs(Y[i,j])>tol or abs(O[i,j])>tol:
                        ratio = O[i,j]/Y[i,j] if abs(Y[i,j])>tol else np.nan
                        print('{:>8s} {:>8s}: {:14.5e} {:14.5e} {:14.5e} {:10.3f}'.format(rows[i], cols[j], Y[i,j], O[i,j], O[i,j]-Y[i,j], ratio))

    def _load_openfast_lin(KF, linFiles):
        """ Reads the OpenFAST linear files (from pickle if available) and returns the averaged matrices """

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

        pklFile = linFiles[0].replace('.lin', f'_nLinFiles={len(linFiles)}.pkl')
        if os.path.exists(pklFile):
            INFO('Loading lin file pickle: ', pklFile)
            FL = FASTLin.from_pickle(pklFile)
        else:
            FL = FASTLin(linfiles=linFiles, prefix='', verbose=False, sX_sel=sX_sel, sU_sel=sU_sel, sY_sel=sY_sel, raiseNaNError=False)
            INFO('Writting pickle file:    ', pklFile)
            FL.save(pklFile)

#         for i, OP in enumerate(FL.OP_Data):
#             print('A OP', OP.linFiles[0])
#             printMat(OP.Data[0].A.values)
        Am, Bm, Cm, Dm = FL.average(return_dataframes=True)
        OF_lin = dict()
        OF_lin['A'] = Am
        OF_lin['B'] = Bm
        OF_lin['C'] = Cm
        OF_lin['D'] = Dm
        OF_lin['FL'] = FL
        return OF_lin

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
        if 'x' in KF.sX:
            fnd_x_q   = np.array([x[0],x[1]])
            fnd_xd_q  = np.array([x_dot[0], x_dot[1]])
            fnd_xdd_q = np.array([x_dot[4], x_dot[5]])
            q  ['Sg']   = x[0]
            q  ['P']    = x[1]
            j=1
            if 'q_FA1' in KF.sX:
                j+=1
                q  ['TFA1'] = x[j]
            if 'Psi' in KF.sX:
                j+=1
                q  ['Psi']  = x[j]

            qd ['Sg']   = x[4]
            qd ['P']    = x[5]
            j=5
            if 'q_FA1' in KF.sX:
                j+=1
                qd ['TFA1'] = x[j]
            if 'psi' in KF.sX:
                j+=1
                qd ['Psi']  = x[j]
            if x_dot is not None:
                qdd['Sg']   = x_dot[4]
                qdd['P']    = x_dot[5]
                j=5
                if 'q_FA1' in KF.sX:
                    j+=1
                    qdd['TFA1'] = x_dot[j]
                if 'psi' in KF.sX:
                    j+=1
                    qdd['Psi']  = x_dot[j]
        else:
            fnd_x_q  =None
            fnd_xd_q =None
            fnd_xdd_q=None
            q  ['TFA1'] = x[0]
            q  ['Psi']  = x[1]

            qd ['TFA1'] = x[0]
            qd ['Psi']  = x[1]
            if x_dot is not None:
                qdd['TFA1'] = x_dot[2]
                qdd['Psi']  = x_dot[3]

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
        """ 
        Inputs : Loads at rotor center (point R)

        """
        WT = KF.WT
        dInfo = KF.dInfo
        # ---
        p_ext      = np.zeros((3,len(KF.zDepth)))
        if p_hydro is not None:
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
        if 'x' in KF.sX:
            xdd_q *=0
            F_sec_h, M_sec_h, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext_for_h, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext = a_ext, PhiK=WT.fnd.PhiK)
            Fx_h_est = F_sec_h[0,0]
        else:
            Fx_h_est = 0
        #wet_nodes = KF.zDepth <= 0
        #Fx_h_est = np.trapezoid(p_hydro[wet_nodes], KF.zDepth[wet_nodes])

        if np.isnan(mnp_sec[4,0]):
            import pdb; pdb.set_trace()

        return rowOut, twr_sec, mnp_sec, Fx_h_est

    # --------------------------------------------------------------------------------}
    # --- Time loop
    # --------------------------------------------------------------------------------{
    def precomputeGains(KF):
        """ Kalman gains and covariances do not depend on the data (LTI model, constant Q and R): computed once,
        and frozen as soon as they converge (steady state)"""
        Ad, Bd, C, Q, R = KF.Xxd, KF.Xud, KF.Yx.values, KF.Q, KF.R
        nt, nX, nY = KF.nt, KF.nX, KF.nY
        Kt = np.zeros((nt, nX, nY)); Pt = np.zeros((nt, nX, nX))
        P = KF.P; I = np.eye(nX)
        Pt[0] = P
        for it in range(1, nt):
            P1m = Ad @ P @ Ad.T + Q
            K   = P1m @ C.T @ np.linalg.inv(C @ P1m @ C.T + R)
            Pn  = (I - K @ C) @ P1m
            Kt[it], Pt[it] = K, Pn
            if it>1 and np.max(np.abs(Pn-P)) <= 1e-12*np.max(np.abs(Pn)):  # steady state reached
                Kt[it+1:], Pt[it+1:] = K, Pn
                break
            P = Pn
        KF.Kt, KF.Pt, KF.P = Kt, Pt, Pt[-1]
        return Kt

    def timeLoop(KF):
        print(f'Time Loop, dt={KF.dt}, t=[{KF.time[0]} - {KF.time[-1]}]')
        nt, nX = KF.nt, KF.nX
        # --- Pre-computations (everything in the loop is numpy only, pandas is too slow)
        Kt     = KF.precomputeGains()
        Ad, Bd = KF.Xxd, KF.Xud
        A, B, C, D = KF.A.values, KF.B.values, KF.C.values, KF.D.values
        Ymeas  = KF.Y.values;       Uin = KF.U_clean.values
        X_hat  = np.zeros((nt, nX)); XD_hat = np.zeros((nt, nX)); Y_hat = np.zeros((nt, KF.nY)); U_hat = np.zeros((nt, KF.nU))
        S_hat  = KF.S_clean.values.copy()   # Row 0 = clean values (initial conditions)
        iS     = {s:i for i,s in enumerate(KF.sS)}
        Thrust_t = np.zeros(nt); Qaero_t = np.zeros(nt); GF_t = np.zeros(nt)
        iU, iY, iX = KF.iU, KF.iY, KF.iX
        hasWSE  = KF.wse is not None
        bFull   = 'x' in KF.sX or KF.method=='OpenFAST' # physical thrust is the input (generalized forces are built in B)
        bGenThrustU, bGenThrustY, bThrustY = 'GenThrust' in KF.sU, 'GenThrust' in KF.sY, 'Thrust' in KF.sY
        bHydro  = 'q_h' in KF.sX
        bPsi    = 'psi' in KF.iX
        bGenThrustX = 'GenThrust' in KF.sX
        # --- Initial conditions
        x = KF.X_clean.values[0].copy()
        X_hat[0] = x; Y_hat[0] = KF.Y_clean.values[0]; U_hat[0] = KF.U_clean.values[0]
        if 'Thrust' in KF.sS:
            Thrust = KF.S_clean['Thrust'].iloc[0]
        elif bGenThrustU:
            Thrust = KF.U_clean['GenThrust'].iloc[0]
        elif bGenThrustY:
            Thrust = KF.Y_clean['Thrust'].iloc[0]
        else:
            Thrust = 0
        GF = Thrust # Approximation at t=0
        if 'WS' in iS:
            WS_last = KF.S_clean['WS'].iloc[0]
        WS_t = np.zeros(nt)

        # --- Time loop: Kalman filter + wind speed estimator (these are the only sequential parts)
        for it in range(0, nt-1):
            y = Ymeas[it+1].copy()
            u = Uin[it+1].copy()
            if bGenThrustU: u[iU['GenThrust']] = GF # We use previous estimated generalized thrust as input.
            if bGenThrustY: y[iY['GenThrust']] = GF # We use previous estimated generalized thrust as measurement.
            if bThrustY   : y[iY['Thrust']] = Thrust # We use previous estimated thrust as measurement.
            # --- Kalman filter, predict and update
            xm = Ad @ x + Bd @ u
            x  = xm + Kt[it+1] @ (y - (C @ xm + D @ u))
            x_dot = A @ x + B @ u
            # --- Estimate Wind Speed and Thrust
            Thrust = 0; Qaero_hat = 0; WS_hat = 0
            if hasWSE:
                if KF.hacks['WSE'] == 'clean_inputs':
                    Qaero_hat = KF.X_clean['Qaero'].iloc[it+1]
                    omega     = KF.X_clean['dpsi'].iloc[it+1]
                else:
                    Qaero_hat = x[iX['Qaero']]
                    omega     = x[iX['dpsi']]
                pitch = u[iU['pitch']] * 180 / np.pi # deg
                WS_hat, _ = KF.wse.estimate(Qaero_hat, pitch=pitch, omega=omega, WS0=WS_last, relaxation=0, method='oper-crossing', t=KF.time[it+1])
                WS_last = float(WS_hat)
                if KF.hacks['thrust'] == 'clean':
                    Thrust = KF.U_clean['GenThrust'].iloc[it+1]
                else:
                    Thrust = KF.wse.Thrust(WS_hat, pitch=pitch, omega=omega)
            # Generalized Thrust Q_q_FA
            GF = Thrust if bFull else KF.WTTN.GF_lin(Thrust, x, bFull=True) # Physical thrust (full, or OpenFAST), otherwise generalized force of YAMS tower-only model
            # --- "Updated"/"hacked" states
            if bPsi:
                x[iX['psi']] = np.mod(x[iX['psi']], 2*np.pi)
            if bGenThrustX:
                x[iX['GenThrust']] = GF
            # --- Store
            X_hat[it+1] = x; XD_hat[it+1] = x_dot; U_hat[it+1] = u; Y_hat[it+1] = C @ x + D @ u
            Thrust_t[it+1] = Thrust; Qaero_t[it+1] = Qaero_hat; GF_t[it+1] = GF; WS_t[it+1] = WS_hat
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it, KF.time[it], WS_hat, Thrust))

        # --- Store results in dataframes
        KF.X_hat.iloc[:,:] = X_hat; KF.XD_hat.iloc[:,:] = XD_hat; KF.Y_hat.iloc[:,:] = Y_hat; KF.U_hat.iloc[:,:] = U_hat
        if 'WS' in iS       : S_hat[1:, iS['WS']]        = WS_t[1:]
        if bHydro           : S_hat[1:, iS['eta']]       = X_hat[1:, iX['q_h']]
        if 'Thrust' in iS   : S_hat[1:, iS['Thrust']]    = Thrust_t[1:]
        if 'GenThrust' in iS: S_hat[1:, iS['GenThrust']] = GF_t[1:]
        KF.S_hat.iloc[:,:] = S_hat
        # --- Section loads (post-processing, estimates are not fed back into the filter)
        KF.sectionLoads(X_hat, XD_hat, Thrust_t)

    # --------------------------------------------------------------------------------}
    # --- Section loads
    # --------------------------------------------------------------------------------{
    def _sl_eval(KF, x, x_dot, eta_dot, Thrust, it=None, full=False):
        """ Section loads for given states, state derivatives, eta_dot, and thrust """
        p_hydro = None
        if 'q_h' in KF.sX:
            p_hydro = KF.pHD['phi'] * eta_dot # p_h = k_h(z) q_h(t)
            p_hydro[KF.zDepth>0] = 0 # safety, shouldn't be necessary
        rowOut, twr_sec, mnp_sec, Fx_h_est = KF.computeSectionLoads(KF.time[1 if it is None else it+1], x, x_dot, p_hydro, Thrust, 0.0, it=it)
        if full:
            return rowOut, twr_sec, mnp_sec, Fx_h_est
        return np.array([mnp_sec[0,0], mnp_sec[4,0], twr_sec[0,0], twr_sec[4,0], Fx_h_est])

    def _sl_linear_map(KF, z0):
        """ Section loads Fx_sb, My_sb, Fx_i, My_i, Fx_h are (very close to) affine functions of z=[x, x_dot, eta_dot, Thrust].
        The affine map is obtained by central finite differences about the mean state z0, and then applied to the whole time series at once"""
        nX = KF.nX; n = 2*nX+2
        hs = np.r_[np.full(2*nX, 1e-1), 1e-1, 1e5]
        def f(z):
            return KF._sl_eval(z[:nX].copy(), z[nX:2*nX].copy(), z[2*nX], z[2*nX+1])
        c = f(z0); L = np.zeros((len(c), n))
        for j in range(n):
            e = np.zeros(n); e[j] = hs[j]
            L[:, j] = (f(z0+e) - f(z0-e)) / (2*hs[j])
        return c - L @ z0, L

    def sectionLoads(KF, X_hat, XD_hat, Thrust_t):
        """ Section loads along the structure, and storage of the main ones in S_hat
        mode (setup_opts['sl_mode']): 'linear' (default, fast: affine map), 'exact' (call to YAMS at each time step) """
        nt = KF.nt; iS = {s:i for i,s in enumerate(KF.sS)}
        mode = KF.setup_opts.get('sl_mode', 'linear')
        if any(KF.hacks[k] for k in ('SL_cleanQ', 'SL_cleanFtop', 'SL_cleanEtaDot', 'SL_cleanP')):
            mode = 'exact'
        KF.dfOut = KF.dInfo['dfOut']
        bHydro = 'q_h' in KF.sX
        eta_dot = X_hat[:, KF.iX['dq_h']] if bHydro else np.zeros(nt)
        S = np.zeros((nt, 5)) # Fx_sb, My_sb, Fx_i, My_i, Fx_h
        if mode == 'linear':
            Z = np.column_stack([X_hat, XD_hat, eta_dot, Thrust_t])[1:]
            c, L = KF._sl_linear_map(Z.mean(axis=0))
            S[1:] = c + Z @ L.T
            S[0] = S[1]
            # --- Quick verification against the exact calculation
            it = nt//2; ex = KF._sl_eval(X_hat[it], XD_hat[it], eta_dot[it], Thrust_t[it], it=it-1)
            err = np.max(np.abs(S[it]-ex)) / np.max(np.abs(ex))
            if err>1e-3: WARN('Linear section loads differ from exact ones by {:.1e} (relative). Use setup_opts["sl_mode"]="exact"'.format(err))
        else:
            cols = list(KF.dfOut.columns); out = np.zeros((nt, len(cols))); perm=None
            for it in range(0, nt-1):
                rowOut, twr_sec, mnp_sec, Fx_h = KF._sl_eval(X_hat[it+1], XD_hat[it+1], eta_dot[it+1], Thrust_t[it+1], it=it, full=True)
                if perm is None: perm = [rowOut.index.get_loc(c) for c in cols]
                out[it+1] = rowOut.values[perm]
                S[it+1] = [mnp_sec[0,0], mnp_sec[4,0], twr_sec[0,0], twr_sec[4,0], Fx_h]
            out[0] = out[1]; S[0] = S[1]
            KF.dfOut.iloc[:,:] = out
            KF.dfOut['Time_[s]'] = KF.time # TODO evaluate at t=0, for now we just replicate the value
        # --- Storage
        for s, j in (('Fx_sb',0), ('My_sb',1), ('Fx_i',2), ('My_i',3), ('Fx_h',4)):
            if s in iS:
                KF.S_hat[s] = S[:, j]
        KF.S_hat.iloc[0, :] = KF.S_clean.iloc[0, :]
        KF.sl_mode = mode

    def computeStats(KF, tRangeStats=None, stats='sigRatio,eps,R2'):
        """ Statistics (sigRatio, eps, R2) between clean values and estimates, for X, Y, S. Returns {channel: {stat: value}} """
        from welib.tools.stats import comparison_stats
        time = np.asarray(KF.time)
        IT = np.ones(len(time), dtype=bool) if tRangeStats is None else np.logical_and(time>tRangeStats[0], time<min(max(time), tRangeStats[1]))
        D = {}
        for clean, hat in ((KF.X_clean, KF.X_hat), (KF.Y_clean, KF.Y_hat), (KF.S_clean, KF.S_hat)):
            for s in hat.columns:
                D[s], _ = comparison_stats(time[IT], clean[s].values[IT], time[IT], hat[s].values[IT], stats=stats, method='1-2')
        return D

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
                if 'Thrust' in KF.sU:
                    Thrust = KF.U_clean['Thrust'].iloc[it]
                elif 'Thrust' in KF.sY:
                    Thrust = KF.Y_clean['Thrust'].iloc[it]
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

def main(fstFile, tmin=0, tmax=20, 
         show=False, export=True,
         compFile=None, hydroShapeFile=None, 
         aeroMapFile=None, operFile=None,
         linFiles=None, 
         method='YAMS', 
         hacks=None,
         Tp=None,
         nUnderSamp=10,
         setup_opts=None,
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
    if setup_opts['aero_est']:
        wse = TabulatedWSEstimator(fstFile=fstFile, operFile=operFile, aeroMapFile=aeroMapFile)
    else:
        wse=None

    eta_ref  = df_ref['Wave1Elev_[m]'].values
    time_ref = df_ref['Time_[s]'].values


    KF = KalmanFilterMTNS(WSE=wse, hacks=hacks, setup_opts=setup_opts)
    KF.setup_matrices(fstFile, 
                      compFile=compFile, hydroShapeFile=hydroShapeFile, Tp=Tp,
                      time_ref= time_ref, eta_ref = eta_ref,
                      method=method, linFiles=linFiles,
                      )

    # --- Loading "Measurements"
    # - Reference file is opened
    # - Measurements are extracted from it
    # - Other signals are extracted from the file, for comparison with estimates. These are referred as "clean" values
    # - Estimate sigmas from measurements (overriden in next section)
    KF.loadMeasurements(measFile=df_ref, tRange=[tmin,tmax], colMap=KF.colMap, timeCol='Time_[s]', raiseIfAbsent=True)

    if 'q_h' in KF.sX:
        KF.X_clean['dq_h'] = np.gradient(KF.X_clean['q_h'], KF.dt)
    # --- Storage for plot
    KF.prepareTimeStepping()
    # --- Process and measurement covariances
    # NOTE: sigs will be squared for P, Q, R
#     sigs = {'x':{}, 'y':{}, 'Q':{}}
#     sigs['y']['PtfmIMUAx'] = np.sqrt(1e-3)    
#     sigs['y']['PtfmIncly'] = np.sqrt(2.7e-7)
#     sigs['y']['dpsi']      = np.sqrt(1e-5)
#     sigs['y']['NcIMUAx']  = np.sqrt(1e-2)    
#     sigs['y']['Qgen']      = np.sqrt(1e-5)
# #     sQ  = ['x', 'phi_y', 'q_FA1', 'psi', 'dx', 'dphi_y', 'dq_FA1', 'dpsi'] # Mechanical states
# #     sQa = ['q_h', 'dq_h', 'Qaero']                                         # Augmented states
#     sigs['Q']['x']      = np.sqrt(KF.dt/dt_ref * 2e-6)
#     sigs['Q']['phi_y']  = np.sqrt(KF.dt/dt_ref * 2e-6)
#     sigs['Q']['q_FA1']  = np.sqrt(KF.dt/dt_ref * 2e-5)
#     sigs['Q']['psi']    = np.sqrt(KF.dt/dt_ref * 2e-6)
# 
#     sigs['Q']['dx']     = np.sqrt(KF.dt/dt_ref * 2e-3)
#     sigs['Q']['dphi_y'] = np.sqrt(KF.dt/dt_ref * 2e-6)
#     sigs['Q']['dq_FA1'] = np.sqrt(KF.dt/dt_ref * 2e-5)
#     sigs['Q']['dpsi']   = np.sqrt(KF.dt/dt_ref * 2e-5)
# 
#     sigs['Q']['q_h']    = np.sqrt(KF.dt/dt_ref * 1e-6)
#     sigs['Q']['dq_h']   = np.sqrt(KF.dt/dt_ref * 1e-4) #KF.Sw
#     sigs['Q']['Qaero']  = np.sqrt(KF.dt/dt_ref * 1e10)
#     sigs['x'] = sigs['Q'].copy()


    dt_ref = 0.02 * nUnderSamp
    sigs = {'x':{}, 'y':{}, 'Q':{}}
    # Measurements - For Matrix R
    if 'x' in KF.sX:
        sigs['y']['PtfmIMUAx'] = np.sqrt(1e-1)    
        sigs['y']['PtfmIncly'] = np.sqrt(2.7e-7)
    if 'dpsi' in KF.sY:
        sigs['y']['dpsi']      = np.sqrt(1e-5)
    sigs['y']['NcIMUAx']  = np.sqrt(1.0)    
    if 'Qgen' in KF.sY:
        sigs['y']['Qgen']      = np.sqrt(1e-5)
    
    # Clean variances of deviations are much smaller
    if 'x' in KF.sX:
        sigs['Q']['x']      = np.sqrt(1e-12)
        sigs['Q']['phi_y']  = np.sqrt(1e-12)
    if 'GenThrust' in KF.sX:
        sigs['Q']['GenThrust']  = np.sqrt(1e14)
    if 'q_FA1' in KF.sX:
        sigs['Q']['q_FA1']  = np.sqrt(1e-12)

    if 'x' in KF.sX:
        sigs['Q']['dx']     = np.sqrt(dt_ref * 1e-4)
        sigs['Q']['dphi_y'] = np.sqrt(dt_ref * 1e-7)
    sigs['Q']['dq_FA1'] = np.sqrt(dt_ref * 1e-5)
    if 'psi' in KF.sX:
        sigs['Q']['psi']    = np.sqrt(1e-12)
        sigs['Q']['dpsi']   = np.sqrt(dt_ref * 1e-5)

    if 'q_h' in KF.sX:
        sigs['Q']['q_h']    = np.sqrt(1e-12)
        sigs['Q']['dq_h']   = np.sqrt(dt_ref * 1e-3)
    if 'Qaero' in KF.sX:
        sigs['Q']['Qaero']  = np.sqrt(dt_ref * 1e10)
    if 'x' in KF.sX and 'q_FA1' in KF.sX: # Full structure with waves (Case 3), tuned variances
        sigs['Q']['Qaero']  = np.sqrt(6.4e13)
        sigs['Q']['dx']     = np.sqrt(5e-5)
        sigs['Q']['dphi_y'] = np.sqrt(6e-8)
        sigs['Q']['dq_FA1'] = np.sqrt(1e-5 if getattr(KF, 'method', 'YAMS') == 'OpenFAST' else 1e-3) # OpenFAST lin. dynamics are more accurate
        sigs['Q']['q_h']    = np.sqrt(2e-6)
        sigs['Q']['dq_h']   = np.sqrt(0.02)
        sigs['y']['NcIMUAx'] = np.sqrt(1e-2)
        if 'Bx' in KF.sX:
            sigs['Q']['Bx']   = np.sqrt(1e9)
            sigs['Q']['Bphi'] = np.sqrt(1e11)
    sigs['x'] = sigs['Q'].copy()


    if not setup_opts['monopileDOFs']:
        NOTE('Using sig value from KF_TNS')
        sigs = {'x':{}, 'y':{}, 'Q':{}}
        if 'q_FA1' in KF.sX:
            sigs['x']['q_FA1']    = 1.0
            sigs['x']['dq_FA1'] = 0.1
        if 'psi' in KF.sX:
            sigs['x']['psi']    = 0.1
            sigs['x']['dpsi']  = 0.1
    #         sigs['x']['Thrust'] = 1000000
        if 'Qaero' in KF.sX:
            sigs['x']['Qaero']  = 8*10**6*1.0
        if 'Bq' in KF.sX:
            sigs['x']['Bq']     = 0.1  # Tower mode bias
    #         sigs['x']['Qgen']   = 1.0*10**6
    #         sigs['x']['WS']     = 1.0
        sigs['Q'] = sigs['x'].copy()
        # Measurements - more or less half the std - For Matrix R
        sigs['y']['NcIMUAx'] = 0.08  # m/s^2
        if 'psi' in KF.sY:
            sigs['y']['dpsi'] = 0.05 # rad/s
        if 'Qgen' in KF.sY:
            sigs['y']['Qgen']  = 1*10**6
        if 'Thrust' in KF.sY:
            sigs['y']['Thrust'] = 1e11 / 1000
        if 'GenThrust' in KF.sY:
            sigs['y']['GenThrust'] = 1e11 / 1000
    #     sigs['y']['pitch'] = 2.00


    if (not setup_opts['q_FA1'] and not setup_opts['aero_est']):
        NOTE('Using sig value from KF_M')
        # TODO TODO TODO VALUES FROM KF_M
        dt_ref = 0.02 # NOTE: Q change with dt
        sigs = {'x':{}, 'y':{}, 'Q':{}}
        sigs['y']['PtfmIMUAx'] = np.sqrt(1e-3)
        sigs['Q']['x']         = np.sqrt(KF.dt/dt_ref * 2e-6)
        sigs['Q']['phi_y']     = np.sqrt(KF.dt/dt_ref * 2e-6)
        sigs['Q']['dx']        = np.sqrt(KF.dt/dt_ref * 2e-3)
        sigs['Q']['dphi_y']    = np.sqrt(KF.dt/dt_ref * 2e-6)
        sigs['Q']['q_h']       = np.sqrt(KF.dt/dt_ref * 2e-6)
        sigs['Q']['dq_h']      = np.sqrt(KF.dt/dt_ref * 2*KF.Sw)




    KF.setupCovariances(
            sigs=sigs,
            useDt=False, Pidentity=True, verbose=False)

    # Initial covariance of the (unknown, large) force biases and aero torque
    if 'Bx' in KF.sX:
        KF.P[KF.iX['Bx'], KF.iX['Bx']]     = KF.setup_opts['P0_Bx']
        KF.P[KF.iX['Bphi'], KF.iX['Bphi']] = KF.setup_opts['P0_Bphi']
    if 'Bq' in KF.sX:
        KF.P[KF.iX['Bq'], KF.iX['Bq']]     = KF.setup_opts['P0_Bq'] # Tower mode bias

    if (not setup_opts['q_FA1'] and not setup_opts['aero_est']):
        FAIL('Somehow the Q will end up different from the KF_M, debug that.')
        KF.Q[:] = np.diag([2e-6, 2e-6, 2e-3, 2e-6, 2e-6, 0.48])

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
    if KF.sl_mode == 'exact':
        file_sl = base + '_SectionLoads_KF_timeloop.outb'
        df_sl.to_outb(file_sl)
        print('Export:', file_sl)

    statsDict = {}    
#     try:

    base = os.path.splitext(fstFile)[0] + f'_KFMTNS_nX={KF.nX}_nU={KF.nU}_nY={KF.nY}_TMax{tmax}_method{method}'
    KF.statsDict = KF.computeStats(tRangeStats)
    statsDict    = KF.statsDict
    if show:
        fig = KF.plot_X( printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        plt.savefig(base + '_X.png')
        KF.plot_Y(printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        if (tmax-tmin)>400:
            plt.savefig(base + '_Y.png')
        fig = KF.plot_S(printStats=True, tRangeStats=tRangeStats, statsDict=statsDict)
        plt.savefig(base + '_S.png')
#         KF.plot_U()
#     except Exception as e:
#         FAIL('Plotting using KF plot functions failed:'+str(e))
    # KF.plot_P()
    # KF.plot_K()
    # KF.plot_innovation()
    

    #KFpkl = fstFile.replace('.fst', f'_KFMTNS_nX={KF.nX}_nU={KF.nU}_nY={KF.nY}_TMax{tRange[1]}.pkl')
    KFpkl = base + '.pkl'
    if export:
        if (tmax-tmin)<400:
            WARN('SKipping export')
        else:
            KF.wse=None # Safety
            print('>>> Writting:', KFpkl)
            KF.save(KFpkl)


    NOTE('Couple of issues to resolve:')
    NOTE(' - Matrices without aero states are quite different from KF_M.')
    NOTE(' - I might be missing the topload contribution to x and phi_y')
    NOTE(' - Time step dependency of the sigmas')


    return KF, df_ref, df_sl


def mainWrapper(fstFile, setup_opts, hacks, method, tRange=[0,600], show=False, linFiles=None):
    if linFiles is None:
        linFiles=[]
        linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1.1.lin')]
        linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim4mps.1.lin')]
        linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim7mps.1.lin')]
        linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim10mps.1.lin')]
        linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim12mps.1.lin')]

    compFile=None

    #compFile       = os.path.join(scriptDir, 'examples/_simulations/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv')
    aeroMapFile    = os.path.join(scriptDir, 'examples/_simulations/IEA-22-280-RWT/IEA-22-280-RWT_Cp_Ct_Cq.rpf')
    operFile       = os.path.join(scriptDir, 'examples/_simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv')

    hydroShapeFile = os.path.join(scriptDir, 'examples/_data/IEA22_HydroShapeFunction_Hs=8.1_Tp=12.7.csv')


    if hacks is None:
        pass
        # --- Super hack
        # hacks = {'thrust':'clean', 'WSE':'clean_inputs', 'SL_cleanQ':True, 
        #              'SL_cleanFtop':True, 'SL_cleanEtaDot':True, 'SL_cleanP':True}
        # --- Intermediate hack: States are exact - Hydro loads are exact
        #     hacks = {'SL_cleanQ':True, 'SL_cleanEtaDot':True, 'SL_cleanP':True} 

        # --- Intermediate hack
        #hacks = {'SL_cleanQ':True, 'SL_cleanEtaDot':True, 'SL_cleanFtop':True} # <<<< EXAMPLE

    KF, _, _ = main(fstFile=fstFile, linFiles=linFiles, 
         compFile=compFile,  hydroShapeFile=hydroShapeFile, Tp=12.7,
         aeroMapFile=aeroMapFile, operFile=operFile,
         hacks=hacks, show=show,
         nUnderSamp=1,
         tmin=tRange[0], tmax=tRange[1],
         method=method, setup_opts=setup_opts,
         export=True,
         )
    return KF, KF.statsDict


def MonopileOnly(method='YAMS', tRange=[0,600], hacks=None, show=False, sl_mode='linear'):

    linFiles=[]
    linFiles += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3_NoRNA_H1A0.1.lin')]
    fstFile    = os.path.join(scriptDir, 'examples/_simulations/06_Jonswap/OF_F3_NoRNA_H1A0_Hs=8.1_Tp=12.7.fst')

    # Explicit setup_opts for Case 1 (Monopile-only)
    setup_opts = {
        'hydro_states': True,           # Wave estimation enabled
        'monopileDOFs': True,           # Use monopile DOFs
        'aero_est': False,              # No aerodynamic estimation
        'q_FA1': False,                 # No tower mode
        'bias_states': True,            # Bias states for monopile (Bx, Bphi)
        'sl_mode': sl_mode,             # Section loads mode
        # Not used but documenting defaults:
        'kh_scale': [1.0, 1.0],         # Wave hydro scale
        'static_bias': [0, 0, 0],       # Static bias offset
        'P0_Bx': 1.0,                   # Covariance for platform surge bias
        'P0_Bphi': 1.0,                 # Covariance for platform pitch bias
        'hybrid': ['KD', 'static', 'hydro'],  # Hybrid model config
    }
    return mainWrapper(fstFile=fstFile, method=method, hacks=hacks, setup_opts=setup_opts, tRange=tRange, show=show, linFiles=linFiles)


def FullStructure_NoWave_NoMonopileDOFs(method='YAMS', tRange=[0,600], hacks=None, show=False, sl_mode='linear'):

    linFiles=[]
    linFiles+= [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H0A1.1.lin')]
    fstFile  =  os.path.join(scriptDir, 'examples/_simulations/06_Jonswap/OF_F3T1S1_H0A1.fst')                

    # Explicit setup_opts for Case 2 (Tower-only)
    setup_opts = {
        'hydro_states': False,          # NO wave estimation (no wave forcing)
        'monopileDOFs': False,          # NO monopile DOFs (tower-only)
        'aero_est': True,               # Aerodynamic force estimation
        'q_FA1': True,                  # Tower first bending mode
        'bias_states': True,            # Bias states: includes Bq for tower mode
        'sl_mode': sl_mode,             # Section loads mode
        # Explicit bias state covariances:
        'P0_Bx': 1.0,                   # Not used (no monopile DOFs)
        'P0_Bphi': 1.0,                 # Not used (no monopile DOFs)
        'P0_Bq': 0.1,                   # Covariance for tower mode bias (Bq)
        'hybrid': ['KD', 'static', 'hydro'],  # Hybrid model config
        'kh_scale': [1.0, 1.0],         # Not used (no wave states)
        'static_bias': [0, 0, 0],       # Static bias offset
    }
    return mainWrapper(fstFile=fstFile, method=method, hacks=hacks, setup_opts=setup_opts, tRange=tRange, show=show, linFiles=linFiles)


def FullStructure_WithWave(method='YAMS', tRange=[0,600], hacks=None, show=False, sl_mode='linear'):
    linFiles=[]
    linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1.1.lin')]
    linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim4mps.1.lin')]
    linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim7mps.1.lin')]
    linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim10mps.1.lin')]
    linFiles      += [os.path.join(scriptDir, 'examples/_simulations/00_EVA/OF_F3T1S1_H1A1_Trim12mps.1.lin')]

    fstFile  = os.path.join(scriptDir, 'examples/_simulations/06_Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst'); 

    # Explicit setup_opts for Case 3 (Full structure)
    setup_opts = {
        'hydro_states': True,           # Wave estimation enabled
        'monopileDOFs': True,           # Use monopile DOFs (platform surge, pitch)
        'aero_est': True,               # Aerodynamic force estimation
        'q_FA1': True,                  # Tower first bending mode
        'bias_states': True,            # Bias states: Bx, Bphi for monopile
        'sl_mode': sl_mode,             # Section loads mode
        # Explicit bias state covariances:
        'P0_Bx': 1.0,                   # Covariance for platform surge bias (Bx)
        'P0_Bphi': 1.0,                 # Covariance for platform pitch bias (Bphi)
        'P0_Bq': 0.1,                   # Not used (q_FA1 coupled to platform)
        'hybrid': ['KD', 'static', 'hydro'],  # Hybrid model config (blocks to replace from OpenFAST)
        'kh_scale': [1.0, 1.0],         # Wave hydro scale (x2 by default, no scale)
        'static_bias': [0, 0, 0],       # Static bias offset (platform forces)
    }
    return mainWrapper(fstFile=fstFile, method=method, hacks=hacks, setup_opts=setup_opts, tRange=tRange, show=show, linFiles=linFiles)


if __name__ == '__main__':

    show=False


#     tRange=[0,600]
    tRange=[210,250]
#     tRange=[210,310]
#     tRange=[5,600]
    
    # --- Case 1 - Monopile under wave - Works well
#     MonopileOnly(tRange=tRange, show=show)

    # --- Case 2 - Monopile and tower under wind: using only tower DOFs - Works well
#     FullStructure_NoWave_NoMonopileDOFs(tRange=tRange, show=show) #, method='OpenFAST')

    # --- Case 3 - Monopile and tower under wind and wave: using all DOFs and wave estimator'
    FullStructure_WithWave(tRange=tRange, show=show, method='OpenFAST')



    if show:
        plt.show()


# def FullStructure_WithWave_NoMonopileDOFs_NoWaveEstimator(method='YAMS', tRange=[0,600], hacks=None, show=False):
# 
#     fstFile  = os.path.join(scriptDir, 'examples/_simulations/06_Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst'); 
#     setup_opts = {'hydro_states':False, 'monopileDOFs':False, 'aero_est':True, 'q_FA1':True}
#     mainWrapper(fstFile=fstFile, method=method, hacks=hacks, setup_opts=setup_opts, show=show)

    # --- Case 2b - Monopile and tower under wind and wave: using only tower DOFs, with no wave estimator - works ok-ish as expected
#     FullStructure_WithWave_NoMonopileDOFs_NoWaveEstimator(tRange=tRange, show=show)
