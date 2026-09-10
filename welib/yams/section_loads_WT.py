""" 
REMEMBER:

- By Far it's dTwrMY/dq_FA1 (C Matrix) that is the most influencial factor for the section loads

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from welib.essentials import *
from welib.tools.strings import latexStrip
from welib.tools.stats import comparison_stats

from welib.yams.windturbine import FASTWindTurbine, rigidBodyKinematics
from welib.weio.dataframe import WEIODataFrame

from welib.fast.linmodel import FASTLinModelFTNSB

from welib.yams.rotations import R_x, R_y, R_z, rotMat
from welib.yams.kinematics import rigidBodyMotion2Points

class YAMSSectionLoadCalculator():
    def __init__(self, fstFile=None, WT=None, HD_compFile=None):
        self.fstFile = fstFile
        self.FAST = None
        self.dfRef = None

        if WT is None:
            self.WT = FASTWindTurbine(fstFile, algo='OpenFAST', HD_compFile=HD_compFile).WT
        else:
            self.WT = WT

    def emptyInputDF(self, nt, inputFrame='R_xs', units=True):
        if units:
            from welib.fast.dofs import COLMAP_QSHORT_TO_QOF
            q_qd_qd2_of = list(COLMAP_QSHORT_TO_QOF.keys())
            sLoads = ['Fadd_R_xs', 'Madd_R_xs']
            cols = ['Time_[s]']+q_qd_qd2_of +sLoads
        else:
            # TODO use OpenFAST Units
            DOFNames = ['Sg', 'Sw', 'Hv' ,'R', 'P', 'Y', 'TFA1', 'TSS1', 'Yaw']
            sq   = ['Q_'+s for s in DOFNames]
            sqd  = ['QD_'+s for s in DOFNames]
            sqdd = ['QD2_'+s for s in DOFNames]
            if inputFrame == 'R_xs':
                # Input at point R in coordinate system xs
                sLoads = ['Fadd_R_xs', 'Madd_R_xs']
            else:
                raise NotImplementedError()
    #             sqdd = ['FN', 'Madd_R_xs']
            cols = ['Time']+sq+sqd+sqdd +sLoads
        data = np.zeros( (nt, len(cols)))
        df   = pd.DataFrame(columns=cols, data=data)
        return df

    def fromDF(self, df, useTopLoadsFromDF=False, useInterfaceLoadsFromDF=False, noAcc=False, dt_resample=None, tRange=None, accMissing='raise'):
        """Compute OpenFAST-like section-load outputs in a single call.

        Typical usage:
            YSL = YAMSSectionLoadCalculator(fstFile=fstFile)
            dfOut = YSL.calc(dfIn)

        For monopile/foundation section outputs:
            dfOut, sections = YSL.fromSD(dfIn)
        """
        from scipy.interpolate import interp1d


        if tRange is not None:
            df=df[df['Time_[s]']>=tRange[0]]
            df=df[df['Time_[s]']<=tRange[1]]
            df.reset_index(inplace=True)

        if dt_resample is not None:
            INFO(f'Resampling at {dt_resample}')
            t_new = np.arange(df['Time_[s]'].min(), df['Time_[s]'].max(), dt_resample )
            f = interp1d(df['Time_[s]'], df.drop(columns=['Time_[s]']), axis=0, kind='linear', fill_value='extrapolate')
            df = pd.DataFrame(f(t_new), columns=df.columns.drop('Time_[s]'))
            df.insert(0, 'Time_[s]', t_new)

        df = WEIODataFrame(df) # Case incensitive
        self.dfRef = df


        dfOut, sections = self.WT.calcOutputsFromDF(df, useTopLoadsFromDF=useTopLoadsFromDF, useInterfaceLoadsFromDF=useInterfaceLoadsFromDF, noAcc=noAcc,
                                                    accMissing=accMissing)

        self.dfOut = dfOut
        self.sec   = sections

#         dfOut.export(outBase + '.YAMS.outb')

        zDepth   = None
        if self.WT.pHD is not None:
            zDepth   = self.WT.pHD['zDepth']

        return dfOut, sections


    # --------------------------------------------------------------------------------}
    # --- Plots 
    # --------------------------------------------------------------------------------{
    def plot_tower_section_loads(self, IsecTwr=None, component='MLyt', figFilename=None, printStats=False):
        dfRef = self.dfRef
        dfOut = self.dfOut

        if IsecTwr is None:
            IsecTwr = [0,4,8] # Section indices

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(12.4,12.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        for iiED,iED in enumerate(IsecTwr):
            sT = f'TwHt{iED+1}{component}_[kN-m]'
            t1 = dfRef['Time_[s]'].values;
            y1 = dfRef[sT].values
            t2 = dfOut['Time_[s]'].values
            y2 = dfOut[sT].values
            ax.plot(t1, y1, 'k-'                        , label='OpenFAST' if iiED==0 else None)
            ax.plot(t2, y2,  '--' , color=fColrs(iiED)  , label='Sim Ht{}'.format(iED+1))
            stats, sStats =  comparison_stats(t1,y1,t2,y2, stats='sigRatio,eps,R2', method='meanabs'); 
            addStats(ax, sT, sStats, printStats=printStats, factY=0.8)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('')
        ax.legend()
        if figFilename is not None:
            fig.savefig(figFilename)
        return fig

    def plot_tower_accelerations(self, IsecTwr=None, component='x', figFilename=None, printStats=False):
        dfRef = self.dfRef
        dfOut = self.dfOut

        if IsecTwr is None:
            IsecTwr = [0,4,8] # Section indices

        n = len(IsecTwr)
        if f'NcIMUTA{component}s_[m/s^2]' in dfRef:
            n+=1

        j=-1
        fig,axes = plt.subplots(n, 1, sharey=False, figsize=(6.4,5.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        if f'NcIMUTA{component}s_[m/s^2]' in dfRef:
            j=j+1;sig = 'NcIMUTAxs_[m/s^2]'; t1, y1, t2, y2 = dfRef['Time_[s]'].values ,dfRef[sig].values ,dfOut['Time_[s]'].values ,dfOut[sig].values; 
            axes[j].plot(t1, y1, 'k-'); axes[j].plot(t1, y2, '--'); axes[j].set_ylabel(sig); 
            stats, sStats =  comparison_stats(t1,y1,t2,y2, stats='sigRatio,eps,R2', method='meanabs'); 
            addStats(ax, sig, sStats, printStats=printStats, factY=0.8)
            #print(sig, stats) # TODO
        for iiED,iED in enumerate(IsecTwr):
            j=j+1; sig = f'TwHt{iED+1}AL{component}t_[m/s^2]'; t1, y1, t2, y2 = dfRef['Time_[s]'].values ,dfRef[sig].values ,dfOut['Time_[s]'].values ,dfOut[sig].values; 
            axes[j].plot(t1, y1, 'k-'); axes[j].plot(t1, y2, '--'); axes[j].set_ylabel(sig); 
            stats, sStats =  comparison_stats(t1,y1,t2,y2, stats='sigRatio,eps,R2', method='meanabs'); 
            addStats(ax, sig, sStats, printStats=printStats, factY=0.8)
            #print(sig, stats)
        if figFilename is not None:
            fig.savefig(figFilename)
        return fig



    def plot_monopile_section_loads_stats(self, IZ=None, component=0, figFilename=None, tRange=None, stat='meanabs'):
        F_sec    = self.sec['monopile']['F_sec']
        F_secRef = self.sec['monopile']['F_secRef']
        zDepth   = self.sec['monopile']['z']
        vTime    = self.dfRef['Time_[s]']
        IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
        fig = sec_plotFM(vTime, zDepth, F_sec, F_secRef, stat=stat, tRange=tRange)
        if figFilename is not None:
            fig.savefig(figFilename)

    
    def plot_monopile_section_loads(self, IZ=None, component=0, figFilename=None, tRange=None, printStats=False):
        """ 
         - component: 0 Fx, 4:My
        """

        F_sec    = self.sec['monopile']['F_sec']
        F_secRef = self.sec['monopile']['F_secRef']
        zDepth   = self.sec['monopile']['z']
        vTime    = self.dfRef['Time_[s]']


        fig, axes = plt.subplots(len(IZ), 1, sharey=True, sharex=True, figsize=(8.4,7.5))
        fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.05, hspace=0.07, wspace=0.20)
        for ii, iz in enumerate(IZ):
            time_plot (vTime, F_secRef[0, iz, :]/1e6,  F_sec[0, iz, :]/1e6, f'z={zDepth[iz]:.0f}m', tRange=tRange, ax=axes[ii])
            stats, sStats =  comparison_stats(vTime, F_secRef[0,iz,:]/1e6, vTime, F_sec[0,iz,:]/1e6, stats='sigRatio,eps,R2', method='meanabs')
            addStats(axes[ii], 'Fsec'+str(component), sStats, printStats=printStats, factY=0.8)
            #print(f'z {zDepth[iz]:5.0f}: ', stats)
        axes[-1].set_xlabel('Time [s]')
        fig.suptitle('FxSec [MN]')
        if figFilename is not None:
            fig.savefig(figFilename)

        return fig


    def plot_comp(self, sig, ylabel=None, figFilename=None, tRange=None, scale=1, ax=None, factY=0.8, printStats=True):
        if ylabel is None:
            ylabel=sig
        vTime = self.dfRef['Time_[s]']
        fig, ax = time_plot (vTime, self.dfRef[sig]/scale, self.dfOut[sig]/scale, ylabel, tRange=tRange, ax=ax)
        stats, sStats =  comparison_stats(vTime, self.dfRef[sig]/scale, vTime, self.dfOut[sig]/scale, stats='sigRatio,eps,R2', method='1-2')
        addStats(ax, sig, sStats, printStats=printStats, factY=0.8)

        if figFilename is not None:
            fig.savefig(figFilename)
        return fig, stats

    def plot_section_loads_stats(self, IZ=None, component=0, figFilename=None, tRange=None, stat='meanabs'):
        F_sec    = self.sec['combined']['F_sec']
        F_secRef = self.sec['combined']['F_secRef']
        zDepth   = self.sec['combined']['z']
        zRef     = self.sec['combined']['zRef']
        vTime    = self.dfRef['Time_[s]']
        fig = sec_plotFM(vTime, zDepth, F_sec, F_secRef, zRef=zRef, stat=stat, tRange=tRange)
        if figFilename is not None:
            fig.savefig(figFilename)
        return fig


# --------------------------------------------------------------------------------}
# --- Optimized functions or dedicated functions 
# --------------------------------------------------------------------------------{

# --- Optimized version see end of this script
class YAMSSectionLoadCalculatorOptimized_FTNS(YAMSSectionLoadCalculator):

    def fromDF(self, df, useTopLoadsFromDF=False, noAcc=False, **kwargs):
        dfSL = calcOutputsFromDFOptimized_FTNS(self.WT, df, useTopLoadsFromDF=useTopLoadsFromDF, noAcc=noAcc, **kwargs)
        return dfSL, None


def calcOutputsFromDFOptimized_FTNS(WT, df, noAcc=False, useTopLoadsFromDF=False, **kwargs):
    """ 
    NOTE: this is a copy of what is found in yams.windturbine.py
    """
    from welib.tools.tictoc import Timer
    from welib.fast.postpro import ED_TwrGag #, ED_TwrStations, getEDClass
    from welib.yams.section_loads import beamSectionLoads3D  # calls beamSectionLoads1D

    # Sanitization of input dataframe
    df = df.loc[:,~df.columns.duplicated()].copy()
    df.columns = [  v.split('_[')[0] for v in df.columns.values]
    df.reset_index(inplace=True)

    # --- States
    DOFNames = ['Sg', 'Sw', 'Hv' ,'R', 'P', 'Y', 'TFA1', 'TSS1', 'Yaw']
    sq   = ['Q_'+s for s in DOFNames]
    sqd  = ['QD_'+s for s in DOFNames]
    sqdd = ['QD2_'+s for s in DOFNames]
    sqall = sq+sqd+sqdd
    for s in sqall:
        if s not in df.keys():
            print('[WARN] Missing DOF from dataframe: {}'.format(s))
            df[s]=0
    Q   = df[sq]
    QD  = df[sqd]
    QDD = df[sqdd]
    Q.columns   = DOFNames
    QD.columns  = DOFNames
    QDD.columns = DOFNames
    if noAcc:
        QDD *=0

    # --------------------------------------------------------------------------------}
    # --- Constants 
    # --------------------------------------------------------------------------------{
    fnd = WT.fnd
    twr = WT.twr
    nac = WT.nac
    r_F0     = fnd.pos_global_init # np.array((0, 0, ED['PtfmRefzt']))
    r_T0     = twr.pos_global_init # np.array((0, 0, ED['TowerBsHt']))
    s_NGn0   = nac.masscenter # TODO
    tilt =WT.shaft_tilt
    rot_type = 'smallRot_OF'
    # --- DOFs
    DOF_f = ['Sg','Sw','Hv','R','P','Y']
    gravity = WT.gravity

    # --------------------------------------------------------------------------------}
    # --- ALLOC 
    # --------------------------------------------------------------------------------{
    nTwrSpan = len(twr.s_span)
    u_Ts_in_t = np.zeros((nTwrSpan,3))
    udd_Ts_in_t = np.zeros((nTwrSpan,3))
    theta_TTs_in_t = np.zeros((nTwrSpan,3))
    r_Ts = np.zeros((nTwrSpan,3))
    v_Ts = np.zeros((nTwrSpan,3))
    a_Ts = np.zeros((nTwrSpan,3))
    R_g2Ts = np.zeros((nTwrSpan,3,3)) 
    theta_TTs = np.zeros((nTwrSpan,3)) 
    theta_Ts  = np.zeros((nTwrSpan,3))
    omega_Ts  = np.zeros((nTwrSpan,3))
    omegad_Ts = np.zeros((nTwrSpan,3))

    nSpan = len(twr.s_span)
    p_ext = np.zeros(nSpan)
    a_struct_t = np.zeros((3,nSpan))



    # --- Outputs
    colOut = ['Time_[s]']
    colOut += sq + sqd + sqdd
    # ED Outputs
    HEDOut, I = ED_TwrGag(WT.ED, addBase=False)
    for iiSL,hED in enumerate(HEDOut):
        sT='TwHt{}'.format(iiSL+1)
        colOut+=[sT+'MLxt_[kN-m]', sT+'MLyt_[kN-m]', sT+'MLzt_[kN-m]']
    # TODO TODO TODO FIGURE OUT WHY THIS RETURN DTYPE OBJECT
    dfOut = pd.DataFrame(index=df.index, columns=colOut)

    # --- Calc Output per time step
    with Timer('Kinematics'):
        for it,t in enumerate(df['Time']):
            qDict   = Q.iloc[it,:].copy()
            qdDict  = QD.iloc[it,:].copy()
            qddDict = QDD.iloc[it,:].copy()
            #dd = kinematicsWT(WT, q, qd, qdd)

            # --------------------------------------------------------------------------------}
            # --- KINEMATICS 
            # --------------------------------------------------------------------------------{
            q_f   = np.array([qDict  [DOF] if DOF in qDict.keys()   else 0 for DOF in DOF_f])
            qd_f  = np.array([qdDict [DOF] if DOF in qdDict.keys()  else 0 for DOF in DOF_f])
            qdd_f = np.array([qddDict[DOF] if DOF in qddDict.keys() else 0 for DOF in DOF_f])
            DOF_t = np.array(['TFA1', 'TFA2', 'TSS1', 'TSS2'])[twr.shapes]
            q_t   = np.array([qDict[DOF] for DOF in DOF_t])
            qd_t  = np.array([qdDict[DOF] for DOF in DOF_t])
            qdd_t = np.array([qddDict[DOF] for DOF in DOF_t])

            # --- Ref point/fnd motion
            r_F      = r_F0 + q_f[:3]
            v_F      = qd_f[:3]
            a_F      = qdd_f[:3]
            theta_t  = q_f  [3:]
            omega_t  = qd_f [3:]
            omegad_t = qdd_f[3:]
            R_t2g    = rotMat(q_f[3:], rot=rot_type)
            R_g2t      = R_t2g.T

            s_FT0_in_f = r_T0-r_F0
            r_FT       = R_t2g.dot(s_FT0_in_f)
            r_T, v_T, a_T = rigidBodyMotion2Points(r_F, v_F, a_F, omega_t, omegad_t, r_FT) 

            # --- Tower section motions
            twr.updateFlexibleKinematics(q_t, qd_t, qdd_t) # yams.bodies.py
            for j in range(nTwrSpan):
                s_TTs0_in_t = twr.s_G0[:,j]  # undisplaced position
                u_Ts_in_t[j,:]   = twr.U[:,j]     # displacement field
                ud_Ts_in_t       = twr.UP[:,j]    # elastic velocity
                udd_Ts_in_t[j,:] = twr.UPP[:,j]    # elastic acceleration

                theta_TTs_in_t[j,:]  = np.array([-twr.V[1,j]  , twr.V[0,j] , 0])
                omega_TTs_in_t  = np.array([-twr.VP[1,j] , twr.VP[0,j], 0])
                omegad_TTs_in_t = np.array([-twr.VPP[1,j] , twr.VPP[0,j], 0])

                theta_TTs[j,:] =  R_t2g.dot(theta_TTs_in_t[j,:] )
                theta_Ts[j,:] =  theta_t + theta_TTs[j,:] # OK because small angle

                R_Ts2t = rotMat(theta_TTs_in_t[j,:], rot=rot_type)
                R_Ts2g = R_t2g.dot(R_Ts2t)
                R_g2Ts[j,:,:] = R_Ts2g.T

                omega_TTs = R_t2g.dot(omega_TTs_in_t)
                omegad_TTs = R_t2g.dot(omegad_TTs_in_t) 
                omega_Ts[j,:] = omega_t + omega_TTs
                omegad_Ts[j,:] = omegad_t + omegad_TTs + np.cross(omega_t, omega_TTs) # TODO double check extra contrib

                s_TTs_in_t  = s_TTs0_in_t + u_Ts_in_t[j,:] # displaced position
                r_TTs = R_t2g.dot(s_TTs_in_t)
                ud_Ts = R_t2g.dot(ud_Ts_in_t)
                udd_Ts = R_t2g.dot(udd_Ts_in_t[j,:])
                r_Ts[j,:] = r_T + r_TTs
                v_Ts[j,:] = v_T + np.cross(omega_t, r_TTs) + ud_Ts
                a_Ts[j,:] = a_T + np.cross(omega_t, np.cross(omega_t, r_TTs)) + np.cross(omegad_t, r_TTs) 
                a_Ts[j,:] += 2* np.cross(omega_t, ud_Ts) +  udd_Ts

            # --- Tower Top point (before Yaw)
            s_TTT0_in_t = twr.s_G0[:,-1] # undisplaced position
            r_TT0 =  r_T0 +  s_TTT0_in_t # undisplaced position of tower top 
            r_TT_undisp =  r_T +  R_t2g.dot(s_TTT0_in_t) # undisplaced, but rotated position of tower top 
            r_TT = r_Ts[-1,:]
            v_TT = v_Ts[-1,:]
            a_TT = a_Ts[-1,:]
            R_g2tt = R_g2Ts[-1,:,:] # To Tower Top
            omega_tt = omega_Ts[-1,:]
            omegad_tt = omegad_Ts[-1,:]
            R_g2p = R_g2tt

            # --- Nacelle Point/Body (last of tower)
            R_g2n  = R_g2tt
            r_N = r_TT
            v_N = v_TT
            a_N = a_TT
            omega_n = omega_tt  
            omegad_n = omegad_tt 

            # --- Shaft
            R_s2n = R_y(tilt)  # Rotation fromShaft to Nacelle
            R_g2s = (R_s2n.T).dot(R_g2n)

            # -- RNA (without Yaw Br) COG
            s_NGrna0_in_N = WT.RNA_noYawBr.masscenter
            dRNA = rigidBodyKinematics(s_NGrna0_in_N, r_N, R_g2n, v_N=v_N, omega_n=omega_n, a_N=a_N, omegad_n=omegad_n, point_name='Grna', source_name='N')
            # --------------------------------------------------------------------------------}
            # --- END KINEMATICS 
            # --------------------------------------------------------------------------------{
            dfOut.loc[it,'Time_[s]'] = t
            # --- Loads
            gravity_vec = np.array([0,0,-WT.gravity])
            # --- RNA (without Yaw Br) loads
            omd_n = omegad_n
            om_n = omega_n
            r_Grna = dRNA['r_Grna']
            a_Grna = dRNA['a_Grna']
            Mrna  = WT.RNA_noYawBr.mass
            JGrna = WT.RNA_noYawBr.masscenter_inertia
            JGrna_g = (R_g2n.T).dot(JGrna).dot(R_g2n)
            F_Grna_grav =  Mrna *gravity_vec
            r_NGrna = dRNA['r_NGrna']

            R_N   = Mrna * a_Grna - F_Grna_grav
            tau_N = np.cross(r_NGrna, R_N)
            tau_N += JGrna_g.dot(omd_n)
            tau_N += np.cross(om_n, JGrna_g.dot(om_n))
            # --- Force at N without YawBr Mass (such are "YawBr" sensors..) in global coordinates
            F_N = -R_N
            M_N = -tau_N   #np.cross(r_NGrna, F_Grna_grav)
            if not useTopLoadsFromDF:
                # Aero force
                # TODO gen?
                if 'Fadd_R_xs' in df.keys():
                    Fadd_R_in_g = R_g2s.T.dot((df.loc[it,'Fadd_R_xs'],0 ,0))
                    Madd_R_in_g = R_g2s.T.dot((df.loc[it,'Madd_R_xs'],0 ,0))
                    r_NR_in_n = WT.rot.pos_global # actually not pos_global but from N
                    r_NR_in_g = R_g2n.T.dot(r_NR_in_n)
                    Madd_R_N = np.cross(r_NR_in_g, Fadd_R_in_g)
                    Fadd_N = Fadd_R_in_g
                    Madd_N = Madd_R_in_g + Madd_R_N*0 # TODO experiment
                    F_N += Fadd_N
                    M_N += Madd_N
#                     else:
#                         raise Exception('Temporary safety')
            F_N_p = R_g2p.dot(F_N)
            M_N_p = R_g2p.dot(M_N)

            if useTopLoadsFromDF:
                F_N_p = np.array((df.loc[it,'YawBrFxp'], df.loc[it, 'YawBrFyp'], df.loc[it,'YawBrFzp']))*1000
                M_N_p = np.array((df.loc[it,'YawBrMxp'], df.loc[it, 'YawBrMyp'], df.loc[it,'YawBrMzp']))*1000
                F_N = (R_g2p.T).dot(F_N_p)
                M_N = (R_g2p.T).dot(M_N_p)
            
            # Yaw Brake contribution at N
            F_N_YawBr = WT.yawBr.mass * gravity_vec
            F_N += F_N_YawBr

            # --- Top Loads in tower coordinates
            F_N_t = R_g2t.dot(F_N)
            M_N_t = R_g2t.dot(M_N)

            # --- Section Loads
            for j in range(nSpan):
                a_struct_t[:,j] = R_g2t.dot(a_Ts[j,:])
            gravity_vec = np.array((0.,0.,-gravity)) # external acceleration (gravity/earthquake)
            a_ext = R_g2t.dot(gravity_vec)
            # NOTE: assumes that U,V, K have been computed using twr.updateFlexibleKinematics 
            F_sec, M_sec, outD =  beamSectionLoads3D(p_ext=p_ext, F_top=F_N_t, M_top=M_N_t, s_span=twr.s_span, m=twr.m, U=twr.U, V=twr.V, K=twr.K, a_struct=a_struct_t, a_ext=a_ext, corrections=1)

            for iiSL, hSL in enumerate(HEDOut):
                iSL = np.argmin(np.abs(hSL-WT.twr.s_span))
                hSL = WT.twr.s_span[iSL]

                sT='TwHt{}'.format(iiSL+1)
#                 dfOut[sT+'FLxt_[kN]'  ].loc[it] = F_sec[0, iSL]/1000
#                 dfOut[sT+'FLyt_[kN]'  ].loc[it] = F_sec[1, iSL]/1000
#                 dfOut[sT+'FLzt_[kN]'  ].loc[it] = F_sec[2, iSL]/1000
                dfOut.loc[it,sT+'MLxt_[kN-m]'] = M_sec[0, iSL]/1000
                dfOut.loc[it,sT+'MLyt_[kN-m]'] = M_sec[1, iSL]/1000
                dfOut.loc[it,sT+'MLzt_[kN-m]'] = M_sec[2, iSL]/1000
    return dfOut




# --------------------------------------------------------------------------------}
# --- Helper functions TODO
# --------------------------------------------------------------------------------{
from welib.tools.colors import MW_Orange, fColrs
import matplotlib
COLRS   = [fColrs(1), MW_Orange ,fColrs(2), fColrs(3)]
try:
    cmap = matplotlib.colormaps['viridis']
except:
    cmap = matplotlib.cm.get_cmap('viridis')
CMAP = [(cmap(v)[0],cmap(v)[1],cmap(v)[2],0.9) for v in np.linspace(0,1,5+1)]
COLRS=[CMAP[1], CMAP[4]]
colRef=COLRS[0]
colSim=COLRS[1]
LWRef=2.4
LWSim=1.5



def addStats(ax, sig, sStats, printStats=False, factY=0.8):
    Ylim = ax.get_ylim(); Xlim = ax.get_xlim()
    ax.text(Xlim[0]+(Xlim[1]-Xlim[0])/1000 ,Ylim[0]+(Ylim[1]-Ylim[0])*factY, sStats, fontsize=10)
    if printStats:
        print(f"{sig:10s} "+latexStrip(sStats))


def meanabs(x, **kwargs):
    return np.mean(np.abs(x), **kwargs)


def time_plot(t, ref=None, sim=None, label='', other=None, ax=None, tRange=None, refLab='OpenFAST', otherLab='Other'):
    figNotProvided = ax is None
    if ax is None:
        fig=plt.figure()
        fig.subplots_adjust(left=0.18, right=0.94, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    if ref is not None:
        ax.plot(t, ref,'-'  , c=colRef, lw=LWRef, label=refLab)
    if sim is not None:
        ax.plot(t, sim,'--' , c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(t, other,':' , label=otherLab)
    if tRange is not None:
        ax.axvline(x = tRange[0] , color = 'k', ls='--')
        ax.axvline(x = tRange[1], color = 'k', ls='--')
    if figNotProvided:
        ax.legend()
        ax.set_xlabel('Time [s]')
    ax.set_ylabel(label)
    ax.tick_params(direction='in')
    return fig, ax


def sec_plotFM(vTime, zDepth, F_sec, F_secRef=None, zRef=None, stat='mean', tRange=None, label='Section Force', other=None, otherLab='Other'):
    if tRange is None:
        IT = np.arange(0,F_sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    if zRef is None:
        zRef = zDepth

    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    if F_sec is not None:
        Fx_sec0 = fstat(F_sec   [0,:,IT],axis=0) # 3, nSpan nT
        My_sec0 = fstat(F_sec   [4,:,IT],axis=0)

    if F_secRef is not None:
        Fx_sec1 = fstat(F_secRef[0,:,IT],axis=0) # 3, nSpan nT
        My_sec1 = fstat(F_secRef[4,:,IT],axis=0)
        bNaN = np.isnan(Fx_sec1) # if sections are missing
        Fx_sec1 = Fx_sec1[~bNaN] 
        My_sec1 = My_sec1[~bNaN] 
        zRef    = zRef[~bNaN] 
        
    if other is not None:
        Fx_sec2 = fstat(other[0,:,IT],axis=0) # 3, nSpan nT
        My_sec2 = fstat(other[4,:,IT],axis=0)

    fig,axes = plt.subplots(1, 2, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    ax=axes[0]
    if F_secRef is not None:
        ax.plot(Fx_sec1/1e6, zRef,    '-' , c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(Fx_sec0/1e6, zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(Fx_sec2/1e6, zDepth,  ':' ,label=otherLab)

    ax.tick_params(direction='in')

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')

    ax=axes[1]
    if F_secRef is not None:
        ax.plot(My_sec1/1e6, zRef , '-',  c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(My_sec0/1e6, zDepth, '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(My_sec2/1e6, zDepth,  ':' ,label=otherLab)
    ax.set_xlabel('Section Moment My [MNm]')
    ax.tick_params(direction='in')
    if F_sec is not None:
        ax.legend()
    fig.suptitle(stat)
    return fig

def sec_plot(vTime, zDepth, sec, ref=None, other=None, stat='mean', tRange=None, label='Section Force', fig=None, refLab='OpenFAST', otherLab='Other'):
    if tRange is None:
        IT = np.arange(0, sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    sec0 = fstat(sec[:,IT], axis=1)
    if ref is not None:
        sec1 = fstat(ref[:,IT], axis=1)
    if other is not None:
        sec2 = fstat(other[:,IT], axis=1)

    if fig is None:
        fig, ax = plt.subplots(1, 1, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    else:
        ax=fig.axes[0]

    if ref is not None:
        ax.plot(sec1, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.plot(sec0,     zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(sec2, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.tick_params(direction='in')
    ax.legend()

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')
    fig.suptitle(stat)
    return fig




class OFLinSectionLoadCalculator():
    def __init__(self, fstLin, DOFs=None, Aero=True, Isec=None):

        if Isec is None:
            Isec = range(9)

        self.Aero=Aero

        model = FASTLinModelFTNSB(fstFilename=fstLin, usePickle=True)
        model.rename(verbose=False)
        sU=[]
        if Aero:
            sU += ['Thrust']
    #     sU += ['Qaero']
        sY =[]
        # sYF=['TwHt{}FLxt_[N]', 'TwHt{}FLyt_[N]', 'TwHt{}FLzt_[N]', 'TwHt{}MLxt_[Nm]', 'TwHt{}MLyt_[Nm]', 'TwHt{}MLzt_[Nm]']
        sYF=['TwHt{}MLyt_[Nm]'] # TODO select other components
        sY = [s.format(i+1) for i in Isec for s in sYF]
        model.extract(sU=sU, sY=sY, check=False)

        # Further reduce the model
        if DOFs is not None:
            model.extract(sX=DOFs.split(',') + ['d'+s for s in DOFs.split(',')], check=False)

        self.model = model

        self.sY = sY


    def fromOF(self, fstSim, qopMethod, uopMethod, yopMethod, tRange, loadMethod, dfWSE=None):
        """ """

        # --- NOTE NOTE NOTE Mean is cheating, use 'lin' for a more fair comparison
        time, dfOF = self.model.setupSimFromOF(fstFilename=fstSim, qopMethod=qopMethod, uopMethod=yopMethod, uMethod='DF', yopMethod=yopMethod, renameFS=True, tRange=tRange)
    #     time, dfOF = model.setupSimFromOF(fstFilename=fstSim, qopMethod='mean', uopMethod='mean', uMethod='DF', yopMethod=yopMethod, renameFS=True, tRange=[tMin, tMax])
    #     time, dfOF = model.setupSimFromOF(fstFilename=fstSim, qopMethod='lin', uopMethod='lin', uMethod='DF', yopMethod='lin', renameFS=True, tRange=[tMin, tMax])
        #A, B, C, D = model.toDataFrames()
        # C[C<1e-14]=0
        # D[D<1e-14]=0
        #print('-------------------- C MATRIX')
        #print(C)
        #print('-------------------- D MATRIX')
        #print(D)
        if self.Aero:
            if loadMethod=='WSE_est':
                dfIN = dfWSE.copy()
                dfIN['Thrust'] = dfWSE['Taero_est_[N]']
                dfIN['Qaero']  = dfWSE['Qaero_est_[N]']
                self.model.setupInputs(uopMethod=uopMethod, uMethod='DF', df=dfIN)
            elif loadMethod=='WSE_ref':
                dfIN = dfWSE.copy()
                dfIN['Thrust'] = dfWSE['Taero_ref_[N]']
                dfIN['Qaero'] = dfWSE['Qaero_ref_[N]']
                self.model.setupInputs(uopMethod=uopMethod, uMethod='DF', df=dfIN)
            elif loadMethod=='OFAero':
                pass

        # --- Option 1 simulate
        # model.simulate(out=True, calc='u,y')
        #dfLI = model.df 
        # --- Option 2, use states from OF and compute base on that
        dfLI = self.model.simulate_fake(calc='u,y')
        return dfLI


