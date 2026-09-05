""" 
Kalman filter model for "Tower Nacelle Shaft" (based on yams TNSB)"

"""

import numpy as np
from .kalman import *
from .kalmanfilter import KalmanFilter
from .kalman_model import AugmentedLinModel
from .filters import moving_average
from welib.ws_estimator.tabulated import TabulatedWSEstimator
from welib.yams.models.TNSB_FAST import FASTmodel2TNSB
from welib.tools.stats import comparison_stats
import welib.fast.fastlib as fastlib
import welib.weio as weio

# --------------------------------------------------------------------------------}
# -- Augmented Linear Model 
# --------------------------------------------------------------------------------{
# This is the complicated step, setting up the state matrices based on various inputs 
class KalmanModelTN(AugmentedLinModel):
    def __init__(self, FstFile, bThrustInStates, debug=False):
        AugmentedLinModel.__init__(self)



        WT = FASTmodel2TNSB(FstFile , shapes_twr=[0],shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z').WT

        if bThrustInStates:
            self.sQ  = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sQa = np.array(['Thrust' ,'Qaero'  ,'Qgen','WS'] )
            self.sY  = np.array(['TTacc','omega','Qgen','pitch'])
            self.sU  = np.array(['pitch'])
        else:
            self.sQ  = np.array(['ut1'  ,'psi'  ,'ut1dot','omega','Qaero','Qgen'] )
            self.sQa = np.array(['Qaero','Qgen','WS'] )
            self.sY  = np.array(['TTacc','omega','Qgen','pitch'])
            self.sU  = np.array(['Thrust','pitch'])

        nGear = WT.ED['GBRatio']

        self.ColMap={
          ' ut1    ' : ' TTDspFA_[m]                   ' ,
          ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
          ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
          ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
          ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
          ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
          ' Qgen   ' : f'{nGear}'+'*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]  # NOTE: nGear
          ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
          ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
          ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
        }


        M,C,K,Ya,Yv,Yq,Yp,Yu,Fp,Fu,Pp,Pq,Pv = EmptySystemMat (len(self.sQ)//2, len(self.sY), len(self.sQa), len(self.sU))

        # This below is problem specific
        nShapes_twr   = 1 # Hard coded for TN
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
        if nShapes_twr==1:
            # TODO aerodamping
            WT.DD  = WT.DD*3.5 # increased damping to account for aero damping
        if debug:
            print(WT)

        self.WT=WT

        # --- Building continuous and discrete state matrices
        M,C,K = WT.MM, WT.DD, WT.KK
        A,B,C,D = BuildSystem_Linear(M,C,K,Ya,Yv,Yq,Fp=Fp,Pp=Pp,Yp=Yp,Yu=Yu,Method='augmented_first_order')

        self.A = A
        self.B = B
        self.C = C
        self.D = D


# --------------------------------------------------------------------------------}
# -- Kalman Filter 
# --------------------------------------------------------------------------------{
# The parts that changes from model to model are the time loop, potentially the measurement preps and postprocessing

class KalmanFilterTN(KalmanFilter):

    def __init__(KF, KM=None, WSE=None, debug=False):
        """

        """
        KalmanFilter.__init__(KF, KM=KM)
        KF.setMat(KM.A, KM.B, KM.C, KM.D)
        KF.WT = KM.WT
        KF.ColMap = KM.ColMap

        KF.wse = WSE # wind speed estimator 
        
    # --- Methods Common between TN and TNLin
    # loadMeasurements, prepareMeasurements    
    def loadMeasurements(KF, MeasFile, nUnderSamp=1, tRange=None, ColMap=None):
        # --- Loading "Measurements"
        nGear  = KF.WT.ED['GBRatio']
        df=weio.read(MeasFile).toDataFrame()
        df=df.iloc[::nUnderSamp,:]                      # reducing sampling
        if tRange is not None:
            df=df[(df['Time_[s]']>= tRange[0]) & (df['Time_[s]']<= tRange[1])] # reducing time range
        time = df['Time_[s]'].values
        dt   = (time[-1] - time[0])/(len(time)-1)
        if ColMap is None:
            ColMap = KF.ColMap
        KF.df = fastlib.remap_df(df, ColMap, bColKeepNewOnly=False)
        # --- 
        KF.discretize(dt, method='exponential')
        KF.setTimeVec(time)
        KF.setCleanValues(KF.df)

        # --- Estimate sigmas from measurements
        #KF.sigX_c, KF.sigY_c, KF.sigQ_c = KF.sigmasFromClean(factor=1)
        sigY, KF.R_c = KF.sigmasYFromClean(factor=1)


    def prepareMeasurements(KF, NoiseRFactor=0, bFilterAcc=False, nFilt=15):
        if KF.R_c is None:
            raise Exception('Cannot prepare measurements with noise with R_c if not set.')
        # --- Creating noise measuremnts
        KF.setYFromClean(R=KF.R_c, NoiseRFactor=NoiseRFactor)
        if bFilterAcc:
            KF.set_vY('TTacc',  moving_average(KF.get_vY('TTacc'),n=nFilt) )

    def timeLoop(KF):
        # --- Initial conditions
        x = KF.initFromClean()
        P = KF.P


        for it in range(0,KF.nt-1):    
            t = it*KF.dt
            # --- "Measurements"
            y  = KF.Y.iloc[it,:].values

            # --- KF predictions
            u=KF.U_clean.iloc[it,:].values.copy()
            x,P,_ = KF.estimateTimeStep(u,y,x,P,KF.Q,KF.R)

            # --- Estimate thrust and WS - Non generic code
            WS_last   = x[KF.iX['WS']]
            pitch     = y[KF.iY['pitch']]*180/np.pi # deg
            Qaero_hat = x[KF.iX['Qaero']]
            omega     = x[KF.iX['omega']]
            WS_hat, _ = KF.wse.estimate(Qaero_hat, pitch, omega, WS_last, relaxation = 0)
            Qaero_hat = np.max(Qaero_hat,0)
            Thrust = KF.wse.Thrust(WS_hat, pitch, omega)
            #GF = Thrust
            GF = KF.WT.GF_lin(Thrust,x,bFull=True)

            x[KF.iX['Thrust']] = GF
            x[KF.iX['WS']]     = WS_hat
            x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)

            # --- Store
            KF.X_hat.iloc[it+1,:] = x
            KF.Y_hat.iloc[it+1,:] = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)
            # --- Propagation to next time step

            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it],x[7],x[4]))

        KF.P = P

    # --- Methods Common between TN and TNLin
    # moments, export, plot_summary, plot_moments
    def moments(KF):
        WT=KF.WT
        z_test = fastlib.ED_TwrGag(WT.ED)[0] - WT.ED['TowerBsHt']
        EI     = np.interp(z_test, WT.twr.s_span, WT.twr.EI[0,:])
        kappa  = np.interp(z_test, WT.twr.s_span, WT.twr.PhiK[0][0,:])
        qx    = KF.X_hat['ut1']
        KF.M_sim = [qx*EI[i]*kappa[i]/1000 for i in range(len(z_test))]                 # in [kNm]
        KF.M_ref   = []
        KF.M_valid = [True]*len(z_test)
        for i in range(len(z_test)):
            try:
                val=KF.df['TwHt{:d}MLyt_[kN-m]'.format(i+1)].values
            except:
                try:
                    val=KF.df['TwHt{:d}MLyt'.format(i+1)].values
                except:
                    KF.M_valid[i] = False
                    val=KF.time*np.nan
            KF.M_ref.append(val)
        return KF.M_sim, KF.M_ref

    def export(KF,OutputFile):
        M=np.column_stack([KF.time]+[KF.X_clean[sj] for j,sj in enumerate(KF.sX)])
        M=np.column_stack([M]+[KF.X_hat  [sj] for j,sj in enumerate(KF.sX)])
        M=np.column_stack([M]+[KF.Y      [sj] for j,sj in enumerate(KF.sY)])
        M=np.column_stack([M]+[KF.Y_hat  [sj] for j,sj in enumerate(KF.sY)])
        if len(KF.sS)>0:
            M=np.column_stack([M]+[KF.S_clean[sj] for j,sj in enumerate(KF.sS)])
            M=np.column_stack([M]+[KF.S_hat  [sj] for j,sj in enumerate(KF.sS)])
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
        from welib.tools.colors import cmap_colors
        from welib.tools.spectral import fft_wrap

        COLRS = cmap_colors(4, 'viridis')

        STATS={}

        def spec_plot(ax,t,ref,sim):
            f1,S1,Info = fft_wrap(t,ref,output_type = 'PSD',averaging = 'Welch', nExp=10, detrend=True)
            f2,S2,Info = fft_wrap(t,sim,output_type = 'PSD',averaging = 'Welch', nExp=10, detrend=True)
            ax.plot(f1,S1,'-' , color=COLRS[0],label='Reference')
            ax.plot(f2,S2,'--', color=COLRS[1],label='simulation')
            ax.set_xlim([0,4])
            ax.set_xlabel('Frequency [Hz]')
            ax.set_yscale('log')

        def time_plot(ax,t,ref,sim, label=''):
            t=t[1:]
            ref=ref[0:-1]
            sim=sim[1:]

            ax.plot(t,ref,'-' , color=COLRS[0])
            ax.plot(t,sim,'--', color=COLRS[1])
            ax.set_ylabel(label)

            # Stats
            stats, sStatsL = comparison_stats(t, ref, t, sim, stats='sigRatio,eps,R2', method='1-2', latex=True)
            stats, sStats = comparison_stats (t, ref, t, sim, stats='sigRatio,eps,R2', method='1-2', latex=False)
            label = label.split('[')[0].strip()
            sStats = f'{label:10s} {sStats}'
            print(sStats)
            STATS[label] = stats

            Ylim=ax.get_ylim()
            Xlim=ax.get_xlim()
            ax.text(Xlim[0],Ylim[0]+(Ylim[1]-Ylim[0])*0.8, sStatsL, fontsize=11 )
            return stats, sStats

        # Aliases to shorten notations
        iX, iY, iS = KF.iX, KF.iY, KF.iS
        X_clean, X_hat = KF.X_clean, KF.X_hat
        S_clean, S_hat = KF.S_clean, KF.S_hat
        XS_clean = pd.concat([X_clean, S_clean], axis=1)
        XS_clean = XS_clean.loc[:, ~XS_clean.columns.duplicated()]
        XS_hat = pd.concat([X_hat, S_hat], axis=1)
        XS_hat = XS_hat.loc[:, ~XS_hat.columns.duplicated()]

        time = KF.time

        ##
        fig, axes = plt.subplots(8, 2, sharey=False, figsize=(13.8,8.8))
        #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        j=-1;
        j+=1; time_plot(axes[j,0], time, X_clean['Qaero']/ 1000, X_hat['Qaero']/ 1000, label='Qaero [kNm]'); 
        spec_plot(axes[j,1], time,X_clean['Qaero']/ 1000, X_hat['Qaero']/ 1000)

        j+=1; time_plot(axes[j,0], time, XS_clean['WS'], XS_hat['WS'], label='WS [m/s]'); 
        spec_plot(axes[j,1], time, XS_clean['WS'], XS_hat['WS'])
        j+=1; time_plot(axes[j,0], time, X_clean['omega'], X_hat['omega'], label='omega [rad/s]');
        spec_plot(axes[j,1], time, X_clean['omega'], X_hat['omega'])
        j+=1; time_plot(axes[j,0], time, XS_clean['Thrust']/1000, XS_hat['Thrust']/1000, label='Thrust [kN]'); 
        spec_plot(axes[j,1], time, XS_clean['Thrust']/1000, XS_hat['Thrust']/1000)
        j+=1; time_plot(axes[j,0], time, XS_clean['ut1'], XS_hat['ut1'], label='ut1 [m]'); 
        spec_plot(axes[j,1], time, XS_clean['ut1'], XS_hat['ut1'])
#         try:
        for i in range(len(KF.M_sim)):
            if KF.M_valid[i]:
                j+=1; time_plot(axes[j,0], time, KF.M_ref[i], KF.M_sim[i], label=f'M{i+1} [kNm]'); 
                spec_plot(      axes[j,1], time, KF.M_ref[i], KF.M_sim[i])
            if j>6:
                break


#         except:
#             pas
        return fig, STATS
        #                                         
    def plot_moments(KF,fig=None,scaleByMean=False):
        import matplotlib
        import matplotlib.pyplot as plt
        from welib.tools.colors import cmap_colors

        z_test = list(fastlib.ED_TwrGag(KF.WT.ED)[0] - KF.WT.ED['TowerBsHt'])
        print('z test:',z_test)
        n=len(z_test)
#         z_test.reverse()
        # --- Compare measurements
        COLRS = cmap_colors(n+1, 'viridis')

        if fig is None:
            fig=plt.figure()
        fig.set_size_inches(6.4,15.0,forward=True) # default is (6.4,4.8)
        for i,z in enumerate(z_test):
            ax = fig.add_subplot(n,1,i+1)
            M_sim =KF.M_sim[i]
            if scaleByMean:
                M_sim+=-np.mean(KF.M_sim[i])+np.mean(KF.M_ref[i])
            
            ax.plot (KF.time, KF.M_ref[i], '-' , color='k',       label='Reference' , lw=1)
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



# --------------------------------------------------------------------------------}
# --- Wrapper For Simulation 
# --------------------------------------------------------------------------------{
def KalmanFilterTNSim(FstFile, MeasFile, OutputFile, aeroMapFile, bThrustInStates, nUnderSamp, tRange, bFilterAcc, nFilt, NoiseRFactor, sigX=None, sigY=None, sigQ=None, bExport=False, ColMap=None, debug=False):
    # ---
    KM = KalmanModelTN(FstFile, bThrustInStates=bThrustInStates)

    # --- Creating a wind speed estimator (reads tabulated aerodynamic data)
    wse = TabulatedWSEstimator(fstFile=FstFile, aeroMapFile=aeroMapFile)

    # ---
    KF = KalmanFilterTN(KM, WSE=wse)
    if debug:
        print(KF.wse)    
        print(KF)
    # --- Loading "Measurements"
    # Defining "clean" values 
    # Estimate sigmas from measurements
    KF.loadMeasurements(MeasFile, nUnderSamp=nUnderSamp, tRange=tRange, ColMap=ColMap)
    KF.sigX=sigX
    KF.sigY=sigY
    KF.sigQ=sigQ
    if debug:
        KF.print_sigmas()

    # --- Process and measurement covariances
    KF.setupCovariances(useDt=False, Pidentity=True)
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


