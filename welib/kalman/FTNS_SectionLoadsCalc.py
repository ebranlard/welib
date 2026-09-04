""" 
REMEMBER:

- By Far it's dTwrMY/dq_FA1 (C Matrix) that is the most influencial factor for the section loads

"""
import numpy as np
# import os
import pandas as pd
# import matplotlib.pyplot as plt
# # Local 
from welib.yams.windturbine import FASTWindTurbine, rigidBodyKinematics
from welib.fast.linmodel import FASTLinModelFTNSB

from welib.yams.rotations import R_x, R_y, R_z, rotMat
from welib.yams.kinematics import rigidBodyMotion2Points



class YAMSSectionLoadCalculator():
    def __init__(self, fstFile=None, WT=None):
        if WT is None:
            self.WT = WT = FASTWindTurbine(fstFile, algo='OpenFAST').WT #, twrShapes=[0]) 
        else:
            self.WT = WT

    def emptyInputDF(self, nt, inputFrame='R_xs'):
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
        #df   = pd.DataFrame(columns=cols, data=data)
        df   = pd.DataFrame()
        return df
    
    def fromDF(self, df, useTopLoadsFromDF=False, noAcc=False):
        """ 
        Keys currently required in df:
          Q_Sg, Q_Hv
          QD_Sg
          QD2_Sg, etc
        """
        # TODO  lots of sanity check, better interfaces
        dfSL = self.WT.calcOutputsFromDF(df, useTopLoadsFromDF=useTopLoadsFromDF, noAcc=noAcc)
        return dfSL

# --- Optimized version see end of this script
class YAMSSectionLoadCalculatorOptimized(YAMSSectionLoadCalculator):

    def fromDF(self, df, useTopLoadsFromDF=False, noAcc=False):
        dfSL = calcOutputsFromDFOptimized(self.WT, df, useTopLoadsFromDF=useTopLoadsFromDF, noAcc=noAcc)
        return dfSL




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

# 
# 
# # --- Script parameters
# qopMethod ='lin' # mean, lin
# uopMethod ='zero' # zero
# yopMethod ='lin'
# # yopMethod ='lin'
# Isec = [0,4,7] # Section indices
# Isec = [0,1,2,3,4,5,6,7,8] # Section indices
# pklFilename = '../../article/code/data/SWT-3p6-130-ED-Stiff_CPCT.pkl'
# noAcc = False  # 
# loadMethod='WSE_est' # WSE_ref, 'WSE_est', 'OF'
# # loadMethod='WSE_ref' # WSE_ref, 'WSE_est', 'OF'
# # loadMethod='OFAero' # 
# # loadMethod='OF' # WSE_ref, 'WSE_est', 'OF'
# # noAcc = True  # 
# # tMin=300
# # tMax=320
# # tMax=20
# 
# tMin=0; tMax=700;
# DOFs=None # DOFs=['x','phi_y']
# 
# 
# # --- Derived parameters
# fstOut = fstSim.replace('.fst','.out')
# outFile = fstOut.replace('.out','').replace('.outb','')+'_t{}'.format(tMax)+'_'+loadMethod
# outFileLI = outFile + '_SL_OFLin'+'_yop={}_uop={}_qop={}'.format(yopMethod, uopMethod, qopMethod)
# outFileSL = outFile + '_SL_YAMS'
# if noAcc:
#     outFileSL += '_noAcc'
# outFileSL += '.csv'
# outFileLI += '.csv'
# outFileFig = outFile + '_SL.png'
# 
# # --- Properties
# WT = FASTWindTurbine(fstSim, algo='OpenFAST') #, twrShapes=[0])
# 
# 
# # --------------------------------------------------------------------------------}
# # --- WSE 
# # --------------------------------------------------------------------------------{
# # --- Aero Load WSE
# if Aero:
#     with Timer('Wind Speed estimation'):
#         wse = TabulatedWSEstimatorFloating(fstFile=fstSim, pickleFile=pklFilename)
#         dfWSE = wse.estimateTimeSeriesFromOF(fstSim, relaxation=0, tRange=[tMin,tMax], method='crossing-oper')
# 
# 
# # --------------------------------------------------------------------------------}
# # --- FASTLinModel 
# # --------------------------------------------------------------------------------{
# print('---------------------- LINEAR SIMULATION ----------------------------------')
# with Timer('Linear simulation'):
#     model = FASTLinModelFTNSB(fstFilename=fstLin, usePickle=True)
#     model.rename(verbose=False)
#     sU=[]
#     if Aero:
#         sU += ['Thrust']
# #     sU += ['Qaero']
#     sY =[]
#     # sYF=['TwHt{}FLxt_[N]', 'TwHt{}FLyt_[N]', 'TwHt{}FLzt_[N]', 'TwHt{}MLxt_[Nm]', 'TwHt{}MLyt_[Nm]', 'TwHt{}MLzt_[Nm]']
#     sYF=['TwHt{}MLyt_[Nm]']
#     sY = [s.format(i+1) for i in Isec for s in sYF]
#     model.extract(sU=sU, sY=sY, check=False)
# 
#     # Further reduce the model
#     if DOFs is not None:
#         model.extract(sX=DOFs.split(',') + ['d'+s for s in DOFs.split(',')], check=False)
# 
#     # --- NOTE NOTE NOTE Mean is cheating, use 'lin' for a more fair comparison
#     time, dfOF = model.setupSimFromOF(fstFilename=fstSim, qopMethod=qopMethod, uopMethod=yopMethod, uMethod='DF', yopMethod=yopMethod, renameFS=True, tRange=[tMin, tMax])
# #     time, dfOF = model.setupSimFromOF(fstFilename=fstSim, qopMethod='mean', uopMethod='mean', uMethod='DF', yopMethod=yopMethod, renameFS=True, tRange=[tMin, tMax])
# #     time, dfOF = model.setupSimFromOF(fstFilename=fstSim, qopMethod='lin', uopMethod='lin', uMethod='DF', yopMethod='lin', renameFS=True, tRange=[tMin, tMax])
#     A, B, C, D = model.toDataFrames()
#     # C[C<1e-14]=0
#     # D[D<1e-14]=0
#     print('-------------------- C MATRIX')
#     print(C)
#     print('-------------------- D MATRIX')
#     print(D)
#     if Aero:
#         if loadMethod=='WSE_est':
#             dfIN = dfWSE.copy()
#             dfIN['Thrust'] = dfWSE['Taero_est_[N]']
#             dfIN['Qaero']  = dfWSE['Qaero_est_[N]']
#             model.setupInputs(uopMethod=uopMethod, uMethod='DF', df=dfIN)
#         elif loadMethod=='WSE_ref':
#             dfIN = dfWSE.copy()
#             dfIN['Thrust'] = dfWSE['Taero_ref_[N]']
#             dfIN['Qaero'] = dfWSE['Qaero_ref_[N]']
#             model.setupInputs(uopMethod=uopMethod, uMethod='DF', df=dfIN)
#         elif loadMethod=='OFAero':
#             pass
# 
#     # --- Option 1 simulate
#     # model.simulate(out=True, calc='u,y')
#     #dfLI = model.df 
#     # --- Option 2, use states from OF and compute base on that
#     dfLI = model.simulate_fake(calc='u,y')
#     dfLI.to_csv(outFileLI, index=False)
# 
# #print(model._inputs_ts)
# 
# print(outFileLI, )
# fig = model.plotCompare(figSize=(6,4), useRenamedFS=True, columns=sY)
# fig.subplots_adjust(left=0.15, right=0.98, top=0.955, bottom=0.15, hspace=0.357, wspace=0.357)
# plt.show()
# # 
# raise Exception()
# 
# 
# 
# 
# # --------------------------------------------------------------------------------}
# # --- Section Loads using WT 
# # --------------------------------------------------------------------------------{
# print('--------------------- SECTION LOADS USING WT ------------------------------')
# 
# # --------------------------------------------------------------------------------}
# # --- PLOT 
# # --------------------------------------------------------------------------------{
# 
# # 
# fig,ax = plt.subplots(1, 1, sharey=False, figsize=(12.4,12.8)) # (6.4,4.8)
# fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
# for iiED,iED in enumerate(Isec):
#     sT = 'TwHt{}MLyt_[kN-m]'.format(iED+1)
#     sT2= 'TwHt{}MLyt_[Nm]'.format(iED+1)
# 
#     t1 = df['Time_[s]'].values
#     y1 = df[sT].values
#     t2 = dfSL['Time_[s]'].values
#     y2 = dfSL[sT].values
#     t3 = dfLI['Time_[s]'].values
#     y3 = dfLI[sT2]/1000
#     stats2, sStats2 =  comparison_stats(t1, y1, t2, y2, stats='sigRatio,eps,R2', method='mean')
#     stats3, sStats3 =  comparison_stats(t1, y1, t3, y3, stats='sigRatio,eps,R2', method='mean')
#     print('')
#     print(sStats2)
#     print(sStats3)
# 
#     ax.plot(t1, y1, 'k-'                        , label='OpenFAST' if iiED==0 else None)
#     ax.plot(t2, y2,  '--' , color=fColrs(iiED)  , label='Sim Ht{}'.format(iED+1))
#     ax.plot(t3, y3,  ':' , color=fColrs(iiED) )# , label='Lin')
# ax.set_xlabel('')
# ax.set_ylabel('')
# # ax.set_ylim([-60e6,  50e6])
# ax.legend()
# 
# fig.savefig(outFileFig)
# plt.show()
# 


# --------------------------------------------------------------------------------}
# --- COPY PASTE 
# --------------------------------------------------------------------------------{
def calcOutputsFromDFOptimized(WT, df, noAcc=False, useTopLoadsFromDF=False):
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



def kinematicsWT(WT, qDict, qdDict, qddDict):
    """ Update kinematics from fnd to blades """
    # see WT.kinematics


    return d

