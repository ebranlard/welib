"""
Kalman filter model for "Monopile"

"""
import os
import numpy as np
import pandas as pd

# Welib
from welib.essentials import *
from welib.weio.fast_linearization_file import FASTLinearizationFile
# Kalman
from welib.kalman.kalman import BuildSystem_Linear_MechOnly 
from welib.kalman.kalmanfilter import KalmanFilter
# YAMS
from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions
from welib.yams.windturbine import monopileSetupFromOpenFAST




# --------------------------------------------------------------------------------}
# --- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterMonopile(KalmanFilter):

    def __init__(KF, debug=False, hacks=None):
        sQ  = ['x','phi_y','dx', 'dphi_y']
        sQa = ['q_h', 'dq_h']
        sU  = ['w'] # White noise
        sY  = ['PtfmIMUAx', 'phi_y']      
        sS  = ['My_sb','Fx_sb', 'eta', 'Fx_h']
        # --- Parent init
        KalmanFilter.__init__(KF, sX0=sQ, sXa=sQa, sU=sU, sY=sY, sS=sS)
        KF.debug = debug
        # Hacks
        hacks_def = {'SL_cleanQ':False, 'SL_cleanFtop':False, 'SL_cleanEtaDot':False, 'SL_cleanP':False}
        if hacks is None:
            KF.hacks = hacks_def
        else:
            hacks_def.update(hacks)
            KF.hacks = hacks_def
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
                       fstFile, 
                       hydroShapeFile=None, comp_file=None, Tp=None, zeta =0.12, qdhScale=1, # Hydro params
                       method='YAMS',
                       lin_file=None, # For method=='OpenFAST
                       ):
        # --- Default arguments
        shapes_sub =[0,4] # TODO detemine this based on sQ
        shapes_twr =[]
        # --- Windturbine model
        WT = FASTmodel2MTNSB(fstFile, shapes_sub=shapes_sub, shapes_twr=shapes_twr, shapes_bld=[],
                             DEBUG=False, bStiffening=True, main_axis='z', fixedShaft=True,
                             algo='OpenFAST').WT

		# --- ColMap
        KF.colMap={
                'x'      : 'Q_Sg_[m]' ,
                'phi_y'  : 'Q_P_[rad]' ,
                'dx'     : 'QD_Sg_[m/s]',
                'dphi_y' : 'QD_P_[rad/s]',
                #'MTacc ' : 'NcIMUTAxs_[m/s^2]' ,
                'PtfmIMUAx' : 'QD2_Sg_[m/s^2]' ,
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
            if hydroShapeFile is not None:
                WT.HD_setShapeFunction(hydroShapeFile)
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
        KF.WT  = WT
        KF.pHD = pHD
        KF.zDepth = WT.fnd.s_span - WT.WtrDpth

        # --- Setup state matrices, problem specific!
        MM_sub = WT.MM[:2,:2] # We only keep surge and pitch
        M_inv     = np.linalg.inv(WT.MM)
        M_inv_sub = np.linalg.inv(MM_sub)
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

        # --- Generalized hydro force 
        IQD   =[KF.iX['dx'], KF.iX['dphi_y']]
        if 'k_h' in pHD:
             A[IQD, KF.iX['dq_h']]  = M_inv_sub @ (pHD['k_h'][0], pHD['k_h'][1])*KF.qdhScale  # qd_h influence in mech DOF
        else:
            FAIL('k_h not present')

        # --- Outputs 
        # Monopile top acceleration
        C[KF.iY['PtfmIMUAx'], :] = A[KF.iX['dx'],:] # PtfmIMUAx is assumed to be qdd_s
        D[KF.iY['PtfmIMUAx'], :] = B[KF.iX['dx'],:] # PtfmIMUAx is assumed to be qdd_s
        # Inclination
        C[KF.iY['phi_y'], KF.iX['phi_y']] = 1   # We measure inclination        

        # Nacelle acceleration in x direction (including pitch coupling)




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

    # prepareTimeStepping 
    def prepareTimeStepping(KF, *args, **kwargs):
        KalmanFilter.prepareTimeStepping(KF, *args, **kwargs)
        KF.dInfo = KF.WT.calcOutputs_init(time=KF.time)

    

    def get_OF_DOFs(KF, x, x_dot=None):
        q   = KF.dInfo['q_default'].copy()
        qd  = KF.dInfo['q_default'].copy()
        qdd = KF.dInfo['q_default'].copy()
        # Model specific
        #        0      1         2      3  
        #sQ  = ['x', 'phi_y', 'dx', 'dphi_y'] # Mechanical states
        fnd_x_q   = np.array([x[0],x[1]])
        fnd_xd_q  = np.array([x_dot[0], x_dot[1]])
        fnd_xdd_q = np.array([x_dot[2], x_dot[3]])

        q  ['Sg']   = x[0]
        q  ['P']    = x[1]


        qd ['Sg']   = x[2]
        qd ['P']    = x[3]

        if x_dot is not None:
            qdd['Sg']   = x_dot[2]
            qdd['P']    = x_dot[3]

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
        qd['Sg']    = KF.df['dx'].iloc[it]
        qd['P']     = KF.df['dphi_y'].iloc[it]
        qdd['Sg']   = KF.df['ddx'].iloc[it]
        qdd['P']    = KF.df['ddphi_y'].iloc[it]
        return q, qd, qdd
    def timeLoop(KF):
        # --- Aliases to shorten notations
        WT = KF.WT
        dInfo = KF.dInfo

        # Prepare section output calculation
        KF.dfOut = dInfo['dfOut']
        
        # --- Initial conditions
        x = KF.initFromClean(var='x,y,u')
        

        # --- Time loop
        for it in range(0, KF.nt-1):    
            t = KF.time[it]
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values
            # --- Inputs
            u = KF.U_clean.iloc[it,:].values.copy()

            # --- Predictions of next time step based on current time step
            t = KF.time[it+1]
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
            
            # --- DOFs in the way expected by calcOutputs_step
            q, qd, qdd, x_q, xd_q, xdd_q = KF.get_OF_DOFs(x, x_dot=x_dot)

            # ---  Section Loads
            F_top = np.array((0.,0.,0.))
            M_top = np.array((0.,0.,0.))
            a_ext = np.array((0.,0.,-WT.gravity)) # external acceleration (gravity/earthquake)

            Thrust = 0
            Qaero = 0
            ser_Loads = pd.Series({'Fadd_R_xs':Thrust, 'Fadd_R_ys':0, 'Fadd_R_zs':0, 'Madd_R_xs':Qaero, 'Madd_R_ys':0, 'Madd_R_zs':0})
            # TEMPORARY for backward compatibility
#             ser_Loads['TwrBsFxt_[kN]'] = 0
#             ser_Loads['TwrBsFyt_[kN]'] = 0
#             ser_Loads['TwrBsFzt_[kN]'] = 0
#             ser_Loads['TwrBsMxt_[kN-m]'] = 0
#             ser_Loads['TwrBsMyt_[kN-m]'] = 0
#             ser_Loads['TwrBsMzt_[kN-m]'] = 0
#             dInfo['useInterfaceLoadsFromDF'] = True

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

            rowOut, twr_it, mnp_it = WT.calcOutputs_step(q, qd, qdd, dInfo, t=KF.time[it+1], ser_Loads=ser_Loads, mnp_p_ext=p_ext)
            KF.dfOut.loc[it+1] = rowOut
            F_sec = mnp_it[0:3,:]
            M_sec = mnp_it[3:6,:]

            #F_sec, M_sec, outD = beamSectionLoadsFromShapeFunctions    (x_q, xd_q, xdd_q, p_ext, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext=a_ext, corrections=0, PhiK=WT.fnd.PhiK)

            # No acceleration # TODO get it from calcOutputs 
            xdd_q *=0
            F_sec_h, M_sec_h, outD = beamSectionLoadsFromShapeFunctions(x_q, xd_q, xdd_q, p_ext_for_h, F_top, M_top, WT.fnd.s_span, WT.fnd.PhiU, WT.fnd.PhiV, WT.fnd.m, a_ext = a_ext, PhiK=WT.fnd.PhiK)
            Fx_h_est = F_sec_h[0,0]
            
            
            # --- Store extra info
            # Environment
            KF.S_hat.at[it+1, 'eta' ]  = eta
            # Loads
            KF.S_hat.at[it+1, 'My_sb']  = M_sec[1,0]
            KF.S_hat.at[it+1, 'Fx_sb']  = F_sec[0,0]
            KF.S_hat.at[it+1, 'Fx_h']   = Fx_h_est
            
            # --- Propagation to next time step
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f ' % (it,KF.time[it]))
        # TODO evaluate at t=0, for now we just replicate the value
        index = KF.dfOut.index
        cols = KF.dfOut.columns.difference(['Time_[s]'])
        KF.dfOut.loc[0, cols] = KF.dfOut.loc[1, cols]
        KF.dfOut.loc[index[-1], cols] = KF.dfOut.loc[index[-2], cols]
