import os
import numpy as np

from welib.kalman.kalman import *
from welib.kalman.kalmanfilter import KalmanFilter
from welib.kalman.filters import moving_average

# --- External dependencies!
import welib.fast.fastlib as fastlib
import welib.weio as weio

#          'WS':'Wind1VelX', 'pitch':'BldPitch1','TTacc':'NcIMUTAxs'}
#          'Thrust':'RotThrust','Qaero':'RtAeroMxh','Qgen':'GenTq',
# NOTE: RotThrust contain gravity and inertia
DEFAULT_COL_MAP={
#   ' ut1    ' : ' TTDspFA_[m]                   ' ,
#   ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
#   ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
#   ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
#   ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
#   ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
#   ' Thrust ' : ' RtFldFxh_[N]                 ' ,
#   ' Qaero  ' : ' RtFldMxh_[N-m]               ' ,
# #   ' Qgen   ' : ' {GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
#   ' Qgen   ' : ' 97*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
#   ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
#   ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
#   ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
  ' x      ' : ' PtfmSurge_[m]                   '              ,
  ' y      ' : ' PtfmSway_[m]                   '               ,
  ' z      ' : ' {PtfmHeave_[m]}                '              ,
  ' phi_x  ' : ' {PtfmRoll_[deg]}   * np.pi/180               ' , # SI [deg] -> [rad]
  ' phi_y  ' : ' {PtfmPitch_[deg]}  * np.pi/180                ', # SI [deg] -> [rad]
  ' phi_z  ' : ' {PtfmYaw_[deg]}    * np.pi/180              '  , # SI [deg] -> [rad]
  #' q_FA1  ' : ' TTDspFA_[m]                   '                ,
  ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   '                , # SI [deg] -> [rad]
  ' q_FA1  ' : ' Q_TFA1_[m]                   '                ,
  ' q_SS1  ' : ' Q_TSS1_[m]                   '                ,
  ' dpsi  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 '                , # SI [rpm] -> [rad/s]
  ' dq_FA1 ' : ' QD_TFA1_[m/s]               '                ,
  ' dq_SS1 ' : ' QD_TSS1_[m/s]               '                ,
  ' dx     ' : ' QD_Sg_[m/s]               '              ,
  ' dy     ' : ' QD_Sw_[m/s]                '               ,
  ' dz     ' : ' QD_Hv_[m/s]               '              ,
  ' dphi_x ' : ' QD_R_[rad/s]                             ' ,
  ' dphi_y ' : ' QD_P_[rad/s]                             ',
  ' dphi_z ' : ' QD_Y_[rad/s]                             '  ,
  ' ddpsi  ' : 'QD2_GeAz_[rad/s^2]'  ,
  ' ddq_FA1' : 'QD2_TFA1_[m/s^2]               '                ,
  ' ddq_SS1' : 'QD2_TSS1_[m/s^2]               '                ,
  ' ddx    ' : 'QD2_Sg_[m/s^2]               '              ,
  ' ddy    ' : 'QD2_Sw_[m/s^2]                '               ,
  ' ddz    ' : 'QD2_Hv_[m/s^2]               '              ,
  ' ddphi_x' : 'QD2_R_[rad/s^2]                             ' ,
  ' ddphi_y' : 'QD2_P_[rad/s^2]                             ',
  ' ddphi_z' : 'QD2_Y_[rad/s^2]                             '  ,
  ' Thrust ' : ' RtFldFxh_[N]                 '                ,
  ' Qaero  ' : ' RtFldMxh_[N-m]               '                ,
#           ' Qgen   ' : ' {GenTq_[kN-m]}  *1000         '             , # [kNm] -> [Nm]
  ' Qgen   ' : ' {GenTq_[kN-m]}  *1000 '+'*{}'.format(1)             , # [kNm] -> [Nm] # <<<<<<< TODO TODO TODO nGear
  ' power  ' : ' {GenPwr_[kW]}  *1000 '            , # [kNm] -> [Nm] # <<<<<<< TODO TODO TODO nGear
  ' WS     ' : ' RtVAvgxh_[m/s]                '                ,
  ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 '                , # SI [deg]->[rad]
  ' NcIMUAx ' : ' NcIMUTAxs_[m/s^2]             ',
  ' NcIMUAy ' : ' NcIMUTAys_[m/s^2]             ',
  ' NcIMUAz ' : ' NcIMUTAzs_[m/s^2]             ',
  ' NcIMUVx ' : ' NcIMUTVxs_[m/s]             ',
  ' NcIMUVy ' : ' NcIMUTVys_[m/s]             ',
  ' NcIMUVz ' : ' NcIMUTVzs_[m/s]             ',
}


class KalmanFilterFTNSLin(KalmanFilter):
    def __init__(KF, KM, AE=None, Extrapolator=None):
        """

        """
        super(KalmanFilterFTNSLin, KF).__init__(sX0=KM.sQ, sXa=KM.sQa, sU=KM.sU, sY=KM.sY, sS=KM.sS)
        iX = KF.iX
        iY = KF.iY
        iU = KF.iU

        # --- Creating a wind speed estimator (reads tabulated aerodynamic data)
        if AE:
            # --- Wind Speed estimator
            if 'phi_y' not in KF.iX:
                print('[WARN] WSE: phi_y not in X, will assume phi_y=0')
            if 'pitch' not in KF.iU:
                print('[WARN] WSE: pitch not in U, will assume pitch=0')
            if not ('Qaero' in KF.iX or 'Qaero' in KF.iU):
                print('[WARN] WSE: not running WSE becasue Qaero is not in X or U')
                pickleFile=None
                KF.wse=None
            else:
                KF.wse=AE
        else:
            KF.wse = None

        KF.KM = KM
        KF.setMat(KM.Xx, KM.Xu, KM.Yx, KM.Yu)

    def loadMeasurements(KF, MeasFile, nUnderSamp=1, tRange=None, ColMap=DEFAULT_COL_MAP):
        # --- Loading "Measurements"
        ext = os.path.splitext(MeasFile)[1]
        if not os.path.exists(MeasFile):
            print('[WARN] Measurement file not found, trying with .out: {}'.format(MeasFile))
            MeasFile = MeasFile.replace(ext,'.out')
        #nGear  = KF.WT.ED['GBRatio']
        df = weio.read(MeasFile).toDataFrame()
        # Remapping/scaling columns to shortname variables
        df = fastlib.remap_df(df, ColMap, bColKeepNewOnly=False)

        nUnderSamp=max(nUnderSamp,1)
        df=df.iloc[::nUnderSamp,:]                      # reducing sampling
        if tRange is not None:
            df=df[(df['Time_[s]']>= tRange[0]) & (df['Time_[s]']<= tRange[1])] # reducing time range
        time = df['Time_[s]'].values
        dt   = (time[-1] - time[0])/(len(time)-1)
        KF.df = df
        # --- 
        KF.discretize(dt, method='exponential')
        KF.setTimeVec(time)
        KF.setCleanValues(KF.df)

        # --- Estimate sigmas from measurements
        sigX_c,sigY_c = KF.sigmasFromClean(factor=1)

    def prepareMeasurements(KF, NoiseRFactor=0, bFilterAcc=False, bFilterOm=False, nFilt=15, bFilterPhi=False):
        # --- Creating noise measuremnts
        KF.setYFromClean(R=KF.R, NoiseRFactor=NoiseRFactor)
        if bFilterAcc:
            if 'NcIMUAx' in KF.Y:
                KF.Y['NcIMUAx'] = moving_average(KF.Y['NcIMUAx'],n=nFilt) 
            if 'NcIMUAy' in KF.Y:
                KF.Y['NcIMUAy'] = moving_average(KF.Y['NcIMUAy'],n=nFilt) 
            if 'NcIMUAz' in KF.Y:
                KF.Y['NcIMUAz'] = moving_average(KF.Y['NcIMUAz'],n=nFilt) 
        if bFilterPhi:
            if 'phi_y' in KF.Y:
                print('>>>>> FILTERING PHI_y')
                KF.Y['phi_y'] = moving_average(KF.Y['phi_y'],n=nFilt) 
        if bFilterOm:
            if 'dpsi' in KF.Y:
                print('>>>>> FILTERING OMEGA')
                KF.Y['dpsi'] = moving_average(KF.Y['dpsi'],n=nFilt) 

    def timeLoop(KF):
        # --- Initial conditions
        x = KF.initFromClean()
        P = KF.P
        KF.U_hat.iloc[0,:] = KF.U_clean.iloc[0,:]

        # --- WSE
        if KF.wse:
            WS_last     = KF.S_clean['WS'][0]
            KF.S_hat['WS'][0]= WS_last
            WSavg      = np.zeros((50,1))
            WSavg[:]=WS_last


        if 'Thrust' in KF.sU:
            Thrust_last = KF.U_clean['Thrust'][0]
            iThrust = list(KF.sU).index('Thrust')

        for it in range(0,KF.nt-1):    
            t = it*KF.dt
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values

            # --- KF predictions
            u=KF.U_clean.iloc[it,:].values.copy()
            if 'Thrust' in KF.sU:
                # We use previous estimated thrust as input.
                u[iThrust] = Thrust_last  
            x,P,_ = KF.estimateTimeStep(u,y,x,P,KF.Q,KF.R)

            # --- Estimate thrust and WS - Non generic code
            if KF.wse:
                if 'WS' in KF.iX:
                    WS_last=x[KF.iX['WS']]
                if 'dpsi' in KF.iX:
                    omega     = x[KF.iX['dpsi']] # in rad/s for WSE
                if 'phi_y' in KF.iX:
                    phiy     = x[KF.iX['phi_y']] * 180/np.pi # in deg for WSE
                elif 'phi_y' in KF.iU:
                    phiy     = u[KF.iU['phi_y']] * 180/np.pi # in deg for WSE
                else:
                    raise Exception('Cannot run WSE if phi_y not in X') # Relax this later
                if 'pitch' in KF.iU:
                    pitch     = u[KF.iU['pitch']]*180/np.pi # in deg for WSE
                if 'Qaero' in KF.iX:
                    Qaero_hat = x[KF.iX['Qaero']]
                elif 'Qaero' in KF.iU:
                    Qaero_hat = u[KF.iU['Qaero']]
                else:
                    raise Exception('Cannot run WSE if Qaero not in X or U')

                #def estimate(self, Qa, omega, pitch,  phiy , WS0, relaxation=0, method='crossing', deltaWSMax=1, verbose=False, debug=False, t=0, WSref=np.nan): 
                WS_hat, _ = KF.wse.estimate(Qaero_hat, omega=omega, pitch=pitch, phiy=phiy, WS0=WS_last, relaxation=0, method='oper-crossing', t=t)
                Qaero_hat = np.max(Qaero_hat,0)
                Thrust = KF.wse.Thrust(WS_hat, omega=omega, pitch=pitch, phiy=phiy)
                GF = Thrust
                # GF = KF.WT2.GF_lin(Thrust,x,bFull=True) # TODO TODO
            else:
                WS_hat=0
                GF     = 0
                omega  = 0
                phiy   = 0
                omega  = 0
                pitch  = 0
                Thrust = 0

            # --- Store
            # TODO TODO WHY U IS NOT STORED?
            if 'Thrust' in KF.iX:
                x[KF.iX['Thrust']] = GF
            elif 'Thrust' in KF.iU:
                pass
                #u[KF.iU['Thrust']] = GF
            elif 'Thrust' in KF.iS:
                KF.S_hat.loc[it+1, 'Thrust']= GF

            if 'WS' in KF.iX:
                x[KF.iX['WS']] = WS_hat
            elif 'WS' in KF.iS:
                KF.S_hat.loc[it+1, 'WS'    ]= WS_hat
            if 'phi_y' in KF.iS:
                KF.S_hat.loc[it+1, 'phi_y'    ]= phiy*np.pi/180
            if 'Qaero' in KF.iS:
                KF.S_hat.loc[it+1, 'Qaero'    ]= Qaero_hat
            if 'dpsi' in KF.iS:
                KF.S_hat.loc[it+1, 'dpsi'    ]= omega
            if 'pitch' in KF.iS:
                KF.S_hat.loc[it+1, 'pitch'    ]= pitch*np.pi/180
            if 'WS0' in KF.iS:
                KF.S_hat.loc[it+1, 'WS0'    ]= WS_last

            if 'psi' in KF.iX:
                x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)

            KF.U_hat.iloc[it+1,:]   = u
            KF.X_hat.iloc[it+1,:]   = x
            KF.Y_hat.iloc[it+1,:]   = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)
            KF.XD_hat.iloc[it+1,:]  = np.dot(KF.Xx,x) + np.dot(KF.Xu,u) # Accelerations
            # --- Propagation to next time step
            if not np.isnan(GF):
                Thrust_last = GF
            WS_last     = WS_hat
            #WSavg[1:] = WSavg[0:-1]
            #WSavg[0]  = WS_hat
# 
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it],WS_hat,Thrust))
        KF.P = P

    # --------------------------------------------------------------------------------}
    # --- Extrapolation (calculation of moments  
    # --------------------------------------------------------------------------------{
    def virtualSensing(KF):
        """ """
        pass


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

