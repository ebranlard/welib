import numpy as np
from .kalman import *
from .kalmanfilter import KalmanFilter
from .filters import moving_average
from .TNSB_FAST import FASTmodel2TNSB

from welib.ws_estimator.tabulated import TabulatedWSEstimator
import welib.fast.fastlib as fastlib
import welib.weio as weio

DEFAULT_COL_MAP={
  ' ut1    ' : ' TTDspFA_[m]                   ' ,
  ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
  ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
  ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
  ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
  ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
  ' Qgen   ' : ' {GenTq_[kN-m]}  *1000 * nGear ' , # [kNm] -> [Nm] LSS
  ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
  ' pitch  ' : ' BldPitch1_[deg]'                , # [deg]->[rad]
  ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
}

class KalmanFilterTN(KalmanFilter):
    def __init__(KF, FstFile, base,  bThrustInStates=True, nShapes_twr=1):
        """

        """

        nShapes_bld   = 0 # Hard coded for TN
        nDOF_2nd      = nShapes_twr+1 # Mech DOFs     :  q = [u, psi]
        
        if nShapes_twr>1:
            raise NotImplementedError()
        if bThrustInStates:
            sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            sAug        = np.array(['Thrust' ,'Qaero'  ,'Qgen','WS'] )
            sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            sInp        = np.array(['pitch'])
        else:
            sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega','Qaero','Qgen'] )
            sAug        = np.array(['Qaero','Qgen','WS'] )
            sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            sInp        = np.array(['Thrust','pitch'])

        super(KalmanFilterTN, KF).__init__(sX0=sStates,sXa=sAug,sU=sInp,sY=sMeas)

        # --- Building state/outputs connection matrices
        M,C,K,Ya,Yv,Yq,Yp,Yu,Fp,Fu,Pp,Pq,Pv = EmptySystemMat (int(KF.nX0/2), KF.nY, KF.nP, KF.nU)

        # This below is problem specific
        if nShapes_twr==1 and bThrustInStates:
            Ya[0,0] = 1    # uddot                     = qddot[0]
            Yv[1,1] = 1    # psidot                    = qdot[1]
            Yp[2,2] = 1    # Direct feed-through of Mg
            Fp[0,0] = 1    # T                         = p[0]
            Fp[1,1] = 1    # dQ                        = p[1] -p[2]
            Fp[1,2] = -1   # dQ                        = p[1] -p[2]
            Yu[3,0] = 1    # pitch direct feedthrough
        else:
            raise NotImplementedError()


        # --- Mechanical system and turbine data
        WT = FASTmodel2TNSB(FstFile , nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, bStiffening=True, main_axis='z')
        if nShapes_twr==1:
            # TODO aerodamping
            WT.DD      = WT.DD*3.5 # increased damping to account for aero damping
        print(WT)
        KF.WT=WT

        # --- Creating a wind speed estimator (reads tabulated aerodynamic data)
        KF.wse = TabulatedWSEstimator(fst_file=FstFile)
        KF.wse.load_files(base=base,suffix='')
        print(KF.wse)
        # --- Building continuous and discrete state matrices
        M,C,K = WT.MM, WT.DD, WT.KK
        Xx,Xu,Yx,Yu = BuildSystem_Linear(M,C,K,Ya,Yv,Yq,Fp=Fp,Pp=Pp,Yp=Yp,Yu=Yu,Method='augmented_first_order')
        KF.setMat(Xx,Xu,Yx,Yu)



    def loadMeasurements(KF, MeasFile, nUnderSamp=1, tRange=None):
        # --- Loading "Measurements"
        ColMap = {'ut1':'TTDspFA', 'psi':'Azimuth','ut1dot':'NcIMUTVxs','omega':'RotSpeed',
                 'Thrust':'RtAeroFxh','Qaero':'RtAeroMxh','Qgen':'GenTq',
                 'WS':'RtVAvgxh', 'pitch':'BldPitch1','TTacc':'NcIMUTAxs'}
        #          'WS':'Wind1VelX', 'pitch':'BldPitch1','TTacc':'NcIMUTAxs'}
        #          'Thrust':'RotThrust','Qaero':'RtAeroMxh','Qgen':'GenTq',
        # NOTE: RotThrust contain gravity and inertia

        # TODO 
        nGear  = KF.WT.ED['GBRatio']
        df=weio.read(MeasFile).toDataFrame()
        df.columns = [  v.split('_[')[0] for v in df.columns.values] 
        if tRange is not None:
            df=df[(df['Time']>= tRange[0]) & (df['Time']<= tRange[1])] # reducing time range
        df=df.iloc[::nUnderSamp,:]                      # reducing sampling
        time = df['Time'].values
        dt   = (time[-1] - time[0])/(len(time)-1)

        df['GenSpeed'] *= 2*np.pi/60 # converted to rad/s
        df['RotSpeed'] *= 2*np.pi/60 # converted to rad/s
        df['GenTq']    *= 1000*nGear # Convert to rot torque
        df['Azimuth']  *= np.pi/180  # rad
        df['RotTorq']  *= 1000 # [kNm]->[Nm]
        df['RotThrust']*= 1000 # [kN]->[N]
        KF.df=df

        # --- 
        KF.discretize(dt, method='exponential')
        KF.setTimeVec(time)
        KF.setCleanValues(df,ColMap)

        # --- Estimate sigmas from measurements
        sigX_c,sigY_c = KF.sigmasFromClean(factor=1)

    def prepareTimeStepping(KF):
        # --- Process and measurement covariances
        KF.P, KF.Q, KF.R = KF.covariancesFromSig()
        # --- Storage for plot
        KF.initTimeStorage()

    def prepareMeasurements(KF, NoiseRFactor=0, bFilterAcc=False, nFilt=15):
        # --- Creating noise measuremnts
        KF.setYFromClean(R=KF.R, NoiseRFactor=NoiseRFactor)
        if bFilterAcc:
            KF.set_vY('TTacc',  moving_average(KF.get_vY('TTacc'),n=nFilt) )



    def timeLoop(KF):
        # --- Initial conditions
        x = KF.initFromClean()
        P = KF.P

        iY={lab: i   for i,lab in enumerate(KF.sY)}
        for it in range(0,KF.nt-1):    
            # --- "Measurements"
            y  = KF.Y[:,it]

            # --- KF predictions
            u=KF.U_clean[:,it]
            x,P,_ = KF.estimateTimeStep(u,y,x,P,KF.Q,KF.R)

            # --- Estimate thrust and WS - Non generic code
            WS0       = x[KF.iX['WS']]
            pitch     = y[KF.iY['pitch']]
            Qaero_hat = x[KF.iX['Qaero']]
            omega     = x[KF.iX['omega']]
            WS_hat = KF.wse.estimate(Qaero_hat, pitch, omega, WS0, relaxation = 0)
            Qaero_hat = np.max(Qaero_hat,0)
            Thrust = KF.wse.Thrust(WS_hat, pitch, omega)
            GF = KF.WT.GF_lin(Thrust,x,bFull=True)

            x[KF.iX['Thrust']] = GF
            x[KF.iX['WS']]     = WS_hat
            x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)

            # --- Store
            KF.X_hat[:,it+1] = x
            KF.Y_hat[:,it+1] = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)

            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it],x[7],x[4]))

        KF.P = P


    def moments(KF):
        WT=KF.WT
        z_test = fastlib.ED_TwrGag(WT.ED) - WT.ED['TowerBsHt']
        EI     = np.interp(z_test, WT.Twr.s_span, WT.Twr.EI[0,:])
        kappa  = np.interp(z_test, WT.Twr.s_span, WT.Twr.PhiK[0][0,:])
        qx    = KF.X_hat[KF.iX['ut1']]
        KF.M_sim = [qx*EI[i]*kappa[i]/1000 for i in range(len(z_test))]                 # in [kNm]
        KF.M_ref = [KF.df['TwHt{:d}MLyt'.format(i+1)].values for i in range(len(z_test)) ] # in [kNm]
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
        try:
            ax=fig.add_subplot(6,2,11)
            time_plot(ax,time,KF.M_ref[2], KF.M_sim[2])
            ax.set_ylabel('My [kNm]')
            ax=fig.add_subplot(6,2,12)
            spec_plot(ax,time,KF.M_ref[2], KF.M_sim[2])
        except:
            pass
        #                                         
    def plot_moments(KF,fig=None):
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
            ax.plot (KF.time, KF.M_ref[i], 'k-', color='k',       label='Reference' , lw=1)
            ax.plot (KF.time, KF.M_sim[i], '--', color=COLRS[i],label='Estimation', lw=0.8)
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



def KalmanFilterTNSim(FstFile, MeasFile, OutputFile, base, bThrustInStates, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX=None, sigY=None, bExport=False):
    # ---
    KF=KalmanFilterTN(FstFile, base,  bThrustInStates=bThrustInStates, nShapes_twr=1)
    print(KF)
    # --- Loading "Measurements"
    # Defining "clean" values 
    # Estimate sigmas from measurements
    KF.loadMeasurements(MeasFile, nUnderSamp=nUnderSamp, tRange=tRange)
    KF.sigX=sigX
    KF.sigY=sigY

    # --- Process and measurement covariances
    # --- Storage for plot
    KF.prepareTimeStepping()
    # --- Creating noise measuremnts
    KF.prepareMeasurements(NoiseRFactor=NoiseRFactor, bFilterAcc=bFilterAcc, nFilt=nFilt)

    # --- Time loop
    KF.timeLoop()
    KF.moments()

    if bExport:
        KF.export(OutputFile)
    return KF

