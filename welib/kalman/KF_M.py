"""
Kalman filter model for "Monopile"

"""
import numpy as np
import pandas as pd
import os
# Welib
from welib.essentials import *
from welib.kalman.kalman import BuildSystem_Linear_MechOnly 
from welib.kalman.kalmanfilter import KalmanFilter

from welib.weio.fast_linearization_file import FASTLinearizationFile
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions
from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.yams.windturbine import monopileSetupFromOpenFAST

# --------------------------------------------------------------------------------}
# --- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that change from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterMonopile(KalmanFilter):

    def __init__(KF, debug=False):
        sQ  = ['x','phi_y','dx', 'dphi_y']
        sQa = ['q_h', 'dq_h']
        sY  = ['TTacc', 'phi_y']
        sU  = ['w'] # White noise
        sS  = ['My_sb','Fx_sb', 'eta', 'Fx_h']
        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)

    def setup_matrices(KF, 
                       fstFilename, 
                       hydro_shape_file=None, comp_file=None, Tp=None, zeta =0.12, qdhScale=1, # Hydro params
                       method='YAMS',
                       lin_file=None, # For method=='OpenFAST
              ):
        # --- Default arguments
        shapes_sub =[0,4] # TODO detemine this based on sQ
        
        # --- Windturbine model
        WT = FASTmodel2MTNSB(fstFilename, shapes_sub=shapes_sub, shapes_twr=[], shapes_bld=[],
                             DEBUG=False, bStiffening=True, main_axis='z', fixedShaft=True,
                             algo='OpenFAST').WT

		# --- ColMap
        KF.colMap={
                'x'      : 'Q_Sg_[m]' ,
                'phi_y'  : 'Q_P_[rad]' ,
                'dx'     : 'QD_Sg_[m/s]',
                'dphi_y' : 'QD_P_[rad/s]',
                'TTacc ' : 'NcIMUTAxs_[m/s^2]' ,
                'eta'    : 'Wave1Elev_[m]', 
                'q_h'    : '{Wave1Elev_[m]}',  # Hack to avoid deletion
                'Fx_h'   : 'HydroFxi_[N]',
                'Fx_sb'  : '-ReactFXss_[N]',
                'My_sb'  : '-ReactMYss_[N*m]',

            }



        if WT.pSS is not None:
            #NOTE('Setting Components', compFile)
            #WT.SS_setComponents(compFile)
            #NOTE('Setting Compute Eta')
            #WT.SS_computeEta(dfRef['Time_[s]'])
            if hydro_shape_file is not None:
                WT.HD_setShapeFunction(hydro_shape_file)
        # Ensure WT.MM contains the hydrodynamic mass if not already added
        GM_hydro = WT.pHD['GM_hydro']
        print('GM_hydro:\n', GM_hydro)
        for i, phi in enumerate(WT.fnd.PhiU):
            # NOTE: diagonal only?
            WT.MM[i,i]+=GM_hydro[i,i]

        pHD = WT.pHD
        # --- Method 2: using Legacy TNSB
        # pST, pSS, pHD, Sys, WT, ref = monopileSetupFromOpenFAST(fstFilename, shapes_sub=shapes_sub, TMIN=tRange[0], TMAX=tRange[1], nSubSample=nUnderSamp, hydroShape=hydroShape, reconHydro=True)
        # WT.fnd = WT.twr


        # --- Store turbine and hydro data in Kalman filter object
        KF.WT = WT
        KF.pHD = pHD

        KF.zDepth = WT.fnd.s_span - WT.WtrDpth

        # --- Setup state matrices, problem specific!
        if method=='YAMS':
            # State matrix A from MCK
            # Empty inputs/outputs B,C,D
            A,B,C,D = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK, nP=len(KF.sXa), nU=len(KF.sU), nY=len(KF.sY), Fp=None)
        elif method=='OpenFAST':
            raise NotImplementedError(f'method={method}')
            # print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> USING A LIN')
            # linFile = '../simulations/MT100/00_EVA/OF_NoHydro.2.lin'
            # linFileH = '../simulations/MT100/00_EVA/OF.3.lin'
            # SX = ['PtfmSurge_[m]', 'PtfmPitch_[rad]', 'd_PtfmSurge_[m/s]', 'd_PtfmPitch_[rad/s]']
            dfA = FASTLinearizationFile(lin_file).toDataFrame()['A']
            # Alin = dfA.loc[SX,SX]
            # dfAH= FASTLinearizationFile(linFileH).toDataFrame()['A']
            # AlinH = dfAH.loc[SX,SX]
            # A[:4,:4] = Alin.values[:4,:4]
            # A[:4,:4] = AlinH.values[:4,:4]
        else:
            raise NotImplementedError(f'method={method}')

        # --------------------------------------------------------------------------------}
        # ---  Code common to OpenFAST and YAMS
        # --------------------------------------------------------------------------------{
        KF.qdhScale = qdhScale

        Minv = np.linalg.inv(WT.MM)
        IQD   =[KF.iX['dx'], KF.iX['dphi_y']]
        A[IQD, KF.iX['dq_h']]  = Minv @ (pHD['k_h'][0], pHD['k_h'][1])*KF.qdhScale  # qd_h influence in mech DOF



        C[KF.iY['TTacc'], :] = A[KF.iX['dx'],:] # TTacc is assumed to be qdd_s
        D[KF.iY['TTacc'], :] = B[KF.iX['dx'],:] # TTacc is assumed to be qdd_s
        C[KF.iY['phi_y'], KF.iX['phi_y']] = 1


        # --- Shaping filter, Hydro state equation
        if Tp==12.7:
            KF.Sw= 2.3835e-01
        elif Tp==10.0:
            KF.Sw= 2.3835e-01/2
        else:
            raise NotImplementedError()
        KF.omega_p = 2*np.pi/Tp
        KF.zeta = zeta
        print('omega_p^2', KF.omega_p**2, '2 zeta omega_p', 2*KF.zeta*KF.omega_p)
        A[KF.iX['q_h'], KF.iX['dq_h']]  = 1
        A[KF.iX['dq_h'], KF.iX['q_h']]  = -KF.omega_p**2
        A[KF.iX['dq_h'], KF.iX['dq_h']] = -2 * KF.zeta * KF.omega_p
        B[KF.iX['dq_h'], KF.iU['w']] = 1 # White noise

        # --- Finally, we set the matrices
        KF.setMat(A, B, C, D)


    
    def timeLoop(KF):
        # --- Aliases to shorten notations
        WT = KF.WT

        # Prepare section output calculation
        dInfo = WT.calcOutputs_init(time=KF.time)
        
        # --- Initial conditions
        x = KF.initFromClean(var='x,y,u')
        

        # --- Time loop
        for it in range(0, KF.nt-1):    
            t = it * KF.dt
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values
            # --- Inputs
            u = KF.U_clean.iloc[it,:].values.copy()

            # --- Predictions of next time step based on current time step
            x, KF.P, _ = KF.estimateTimeStep(u, y, x, KF.P)

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

            Thrust = 0
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
            
            
            # --- Store extra info
            # Environment
            KF.S_hat.at[it+1, 'eta' ] = eta
            # Loads
            KF.S_hat.at[it+1, 'My_sb'] = M_sec[1,0]
            KF.S_hat.at[it+1, 'Fx_sb'] = F_sec[0,0]

            KF.S_hat.at[it+1, 'Fx_h' ] = F_sec_h[0,0]
            
            # --- Propagation to next time step
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f ' % (it,KF.time[it]))
