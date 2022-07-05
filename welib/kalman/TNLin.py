import numpy as np
from .kalman import *
from .kalmanfilter import KalmanFilter
from .filters import moving_average
from welib.ws_estimator.tabulated import TabulatedWSEstimator
import yams
from welib.yams.TNSB_FAST import FASTmodel2TNSB
from welib.fast.linmodel import FASTLinModel

# --- External dependencies!
import welib.fast.fastlib as fastlib
import welib.weio as weio

#          'WS':'Wind1VelX', 'pitch':'BldPitch1','TTacc':'NcIMUTAxs'}
#          'Thrust':'RotThrust','Qaero':'RtAeroMxh','Qgen':'GenTq',
# NOTE: RotThrust contain gravity and inertia
DEFAULT_COL_MAP={
  ' ut1    ' : ' TTDspFA_[m]                   ' ,
  ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
  ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
  ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
  ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
  ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
#   ' Qgen   ' : ' {GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
  ' Qgen   ' : ' 97*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
  ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
  ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
  ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
}


class KalmanFilterTNLin(KalmanFilter):
    def __init__(KF, KM, FstFile, base, StateFile):
        """

        """
        super(KalmanFilterTNLin, KF).__init__(sX0=KM.sStates, sXa=KM.sAug, sU=KM.sInp, sY=KM.sMeas, sS=KM.sStor)
        iX = KF.iX
        iY = KF.iY
        iU = KF.iU

        # --- Mechanical system and turbine data
        WT2= FASTmodel2TNSB(FstFile , nShapes_twr=1,nShapes_bld=0, DEBUG=False, bStiffening=True, main_axis='z')    
        #WT2.DD      = WT2.DD*3.5 # increased damping to account for aero damping
        KF.WT2=WT2

        WT  = FASTLinModel(FstFile, StateFile=StateFile, DEBUG=False)
        KF.WT=WT
        A,B,C,D,M = WT.A, WT.B, WT.C, WT.D, WT.M # To Shorten notations

        # --- Creating a wind speed estimator (reads tabulated aerodynamic data)
        KF.wse = TabulatedWSEstimator(fst_file=FstFile)
        KF.wse.load_files(base=base,suffix='')
        # --- Build linear system
        nX = len(KM.sStates)+len(KM.sAug)
        nU = len(KM.sInp   )
        nY = len(KM.sMeas  )
        nq = len(KM.sStates)
        #
        nGear = WT.nGear
        Mqt      =  1/B.iloc[2,0]
        J_LSSnG  = -1/B.iloc[3,1]  # This is JLSS * nGear (scaled by nGea since it takes Torque at HSS and return influences at LSS). It is not J_HSS !!

        Mqt_ED   = M.iloc[0,0]
        J_LSS_ED = M.iloc[1,1]
        
        Xx, Xu, Yx, Yu = EmptyStateMat(nX, nU, nY)
        # --- Filling extended state matrices
        KF.KM = KM
        if KM.StateModel=='nt1_nx8' or KM.StateModel=='nt1_nx7': # sAug =  ['Thrust','Qaero','Qgen','WS']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values[:,2:]
            Yu[:  ,:   ] = D.values[:,2:]
            #----
            Xx[iX['omega'],iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            Xx[:nq,iX['Thrust']]  = B.values[:,0]
            Yx[:,  iX['Thrust']]  = D.values[:,0]
            Yx[iY['Qgen'],iX['Qgen']] = 1
            # --- Value Hack
            if KM.ThrustHack:
                Xx[iX['ut1dot'], iX['Thrust']] =  2.285e-06  # Thrust
#             Xx[iX['omega'],  iX['Qaero']]  =  2.345e-08  # Qa
#             Xx[2,0:4] =[ -6.132e+00,      0,   -5.730e-02,      0]
#             Xx[3,4  ] =  0
#             Xu[3,0  ] =  0
            # --- Consistency
            if KM.Qgen_LSS:
                Xx[iX['omega'],  iX['Qgen']]   =-Xx[iX['omega'],iX['Qaero']]
            else:
                Xx[iX['omega'],  iX['Qgen']]   =-Xx[iX['omega'],iX['Qaero']]*nGear
            Yx[0,0:] =Xx[2,0:]  # <<<< Important


        elif KM.StateModel=='nt1_nx6': # sAug = ['Qaero','Thrust']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values[:,1:]
            Yu[:  ,:   ] = D.values[:,1:]
            #----
            Xx[iX['omega'],  iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            Xx[:nq,iX['Thrust']] = B.values[:,0]
            Yx[:,  iX['Thrust']] = D.values[:,0]
            #  Value Hack
#             Xx[2,0:4] =[ -6.132e+00,      0,   -5.730e-02,      0]
            if KM.ThrustHack:
                Xx[2,iX['Thrust']] =  2.285e-06  # Thrust
#             Xx[3,iX['Qaero' ]] =  2.345e-08  # Torque
#             Xx[3,4  ] =  0
#             Xu[2,0  ] =  0
#             Yu[0,0  ] =  0
            # Consistency
            if KM.Qgen_LSS:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]
            else:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]*nGear
            Yx[0,0:6] =Xx[2,0:6]  # <<<< Important

        elif KM.StateModel=='nt1_nx5':  # sAug = ['Qaero']
            Xx[:nq,:nq ] = A.values
            Yx[:  ,:nq ] = C.values
            #----
            Xu[:nq,:nU ] = B.values
            Yu[:  ,:   ] = D.values
            #----
            Xx[iX['omega'],  iX['Qaero']] = 1/J_LSS_ED # ddpsi Qa # NOTE: LSS
            #  Value Hack
#             Xx[2,0:4] =[ -6.132e+00,      0,   -5.730e-02,      0]
            if KM.ThrustHack:
                Xu[2,0  ] =  2.285e-06  # Thrust
#             Xx[3,4]   =  2.345e-08  # Torque
            # Consistency
            if KM.Qgen_LSS:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]
            else:
                Xu[iX['omega'],iU['Qgen']]  =-Xx[iX['omega'],iX['Qaero']]*nGear
            Yx[0,0:4] = Xx[2,0:4]
            Yu[0,0]   = Xu[2,0]



        KF.setMat(Xx, Xu, Yx, Yu)



    def loadMeasurements(KF, MeasFile, nUnderSamp=1, tRange=None, ColMap=DEFAULT_COL_MAP):
        # --- Loading "Measurements"
        nGear  = KF.WT.ED['GBRatio']
        df=weio.read(MeasFile).toDataFrame()
        df=df.iloc[::nUnderSamp,:]                      # reducing sampling
        if tRange is not None:
            df=df[(df['Time_[s]']>= tRange[0]) & (df['Time_[s]']<= tRange[1])] # reducing time range
        time = df['Time_[s]'].values
        dt   = (time[-1] - time[0])/(len(time)-1)
        KF.df = fastlib.remap_df(df, ColMap, bColKeepNewOnly=False)
        # --- 
        KF.discretize(dt, method='exponential')
        KF.setTimeVec(time)
        KF.setCleanValues(KF.df)

        # --- Estimate sigmas from measurements
        sigX_c,sigY_c = KF.sigmasFromClean(factor=1)

    def prepareMeasurements(KF, NoiseRFactor=0, bFilterAcc=False, nFilt=15):
        # --- Creating noise measuremnts
        KF.setYFromClean(R=KF.R, NoiseRFactor=NoiseRFactor)
        if bFilterAcc:
            KF.set_vY('TTacc',  moving_average(KF.get_vY('TTacc'),n=nFilt) )

    def timeLoop(KF):
        # --- Initial conditions
        x = KF.initFromClean()
        WS_last     = KF.S_clean[KF.iS['WS'    ],0]
        KF.S_hat[KF.iS['WS'    ], 0]= WS_last

        if not KF.KM.bThrustInStates:
            Thrust_last = KF.S_clean[KF.iS['Thrust'],0]
            KF.S_hat[KF.iS['Thrust'], 0]= Thrust_last
        
        KF.X_hat[:,0]   = x
        P = KF.P

        WSavg      = np.zeros((50,1))
        WSavg[:]=WS_last

        for it in range(0,KF.nt-1):    
            # --- "Measurements"
            y  = KF.Y[:,it]

            # --- KF predictions
            u=KF.U_clean[:,it]
            if not KF.KM.bThrustInStates:
                u[0] = Thrust_last # (we don't know the thrust)
            x,P,_ = KF.estimateTimeStep(u,y,x,P,KF.Q,KF.R)

            # --- Estimate thrust and WS - Non generic code
            if KF.KM.bWSInStates:
                WS_last=x[KF.iX['WS']]
            pitch     = y[KF.iY['pitch']]*180/np.pi # deg
            Qaero_hat = x[KF.iX['Qaero']]
            omega     = x[KF.iX['omega']]
            WS_hat = KF.wse.estimate(Qaero_hat, pitch, omega, WS_last, relaxation = 0, WSavg=np.mean(WSavg))
            Qaero_hat = np.max(Qaero_hat,0)
            Thrust = KF.wse.Thrust(WS_hat, pitch, omega)

            GF = Thrust
            GF = KF.WT2.GF_lin(Thrust,x,bFull=True)

            # --- Store
            if KF.KM.bThrustInStates:
                x[KF.iX['Thrust']] = GF
            else:
                KF.S_hat[KF.iS['Thrust'], it+1]= GF
            if KF.KM.bWSInStates:
                x[KF.iX['WS']] = WS_hat
            KF.S_hat[KF.iS['WS'    ], it+1]= WS_hat
            x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)
            KF.X_hat[:,it+1]   = x
            KF.Y_hat[:,it+1]   = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)
            # --- Propagation to next time step
            Thrust_last = GF
            WS_last     = WS_hat
            WSavg[1:] = WSavg[0:-1]
            WSavg[0]  = WS_hat

            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it],WS_hat,Thrust))

        KF.P = P

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



def KalmanFilterTNLinSim(KM, FstFile, MeasFile, OutputFile, base, StateFile, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX=None, sigY=None, bExport=False, ColMap=DEFAULT_COL_MAP, debug=True):
    # ---
    KF=KalmanFilterTNLin(KM, FstFile, base, StateFile)
    if debug:
        print(KF.wse)
        print(KF.WT)
        print(KF)
    # --- Loading "Measurements"
    # Defining "clean" values 
    # Estimate sigmas from measurements
    KF.loadMeasurements(MeasFile, nUnderSamp=nUnderSamp, tRange=tRange, ColMap=ColMap)
    KF.sigX=sigX
    KF.sigY=sigY

    # --- Process and measurement covariances
    # --- Storage for plot
    KF.prepareTimeStepping()
    # --- Creating noise measuremnts
    KF.prepareMeasurements(NoiseRFactor=NoiseRFactor, bFilterAcc=bFilterAcc, nFilt=nFilt)

    # --- Time loop
    if debug:
        print(OutputFile)
    KF.timeLoop()
    KF.moments()

    if bExport:
        KF.export(OutputFile)
    return KF

