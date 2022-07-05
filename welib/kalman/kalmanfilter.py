from .kalman import *
import numpy as np
import pandas as pd

class KalmanFilter(object):
    def __init__(self,sX0,sXa,sU,sY,sS=None):
        sS = [] if sS is None else sS
        self.sX0 = sX0
        self.sXa = sXa
        self.sU  = sU
        self.sY  = sY
        self.sS  = sS # Storage, "Misc" values

        #  State vector is States and Augmented states
        self.sX=np.concatenate((self.sX0,self.sXa))

        # --- Defining index map for convenience
        self.iX={lab: i   for i,lab in enumerate(self.sX)}
        self.iY={lab: i   for i,lab in enumerate(self.sY)}
        self.iU={lab: i   for i,lab in enumerate(self.sU)}
        self.iS={lab: i   for i,lab in enumerate(self.sS)}

        # Standard deviations and covariance matrix
        self.sigX_c = None
        self.sigY_c = None
        self.sigX   = None
        self.sigY   = None
        self.P = None
        self.Q = None
        self.R = None

    @property
    def nX(self):
        return len(self.sX)

    @property
    def nY(self):
        return len(self.sY)

    @property
    def nU(self):
        return len(self.sU)

    @property
    def nP(self):
        return len(self.sXa)

    @property
    def nX0(self):
        return len(self.sX0)

    @property
    def nS(self):
        return len(self.sS)

    def __repr__(self):
        def pretty_PrintMat(M,fmt='{:11.3e}',fmt_int='    {:4d}   ',sindent='   '):
            s=str(M)
#             s=sindent
#             for iline,line in enumerate(M):
#                 s+=''.join([(fmt.format(v) if int(v)!=v else fmt_int.format(int(v))) for v in line ])
#                 s+='\n'+sindent
            return s

        s=''
        s+='<kalman.KalmanFilter object> \n'
        s+='  sX  : {} \n'.format(self.sX)
        s+='  sX0 : {} \n'.format(self.sX0)
        s+='  sX1 : {} \n'.format(self.sXa)
        s+='  sU  : {} \n'.format(self.sU)
        s+='  sY  : {} \n'.format(self.sY)
        s+='  sS  : {} \n'.format(self.sS)
        if self.Xx is not None:
            s+=' Xx: State-State Matrix  \n'
            s+=pretty_PrintMat(self.Xx)+'\n'
        if self.Xu is not None:
            s+=' Xu: State-Input Matrix  \n'
            s+=pretty_PrintMat(self.Xu)+'\n'
        if self.Yx is not None:
            s+=' Yx: Output-State Matrix  \n'
            s+=pretty_PrintMat(self.Yx)+'\n'
        if self.Yu is not None:
            s+=' Yu: Output-Input Matrix  \n'
            s+=pretty_PrintMat(self.Yu)+'\n'
        return s

    @property
    def A(self):
        return self.Xx

    @property
    def B(self):
        return self.Xu

    @property
    def C(self):
        return self.Yx

    @property
    def D(self):
        return self.Yu


    def setMat(self, Xx, Xu, Yx, Yu):
        # --- 
        self.Xx, self.Xu, self.Yx, self.Yu= EmptyStateDF(self.nX,self.nU,self.nY, self.sX, self.sU, self.sY)

        if Xx.shape != self.Xx.shape:
            raise Exception('Shape of Xx ({}) not compatible with KF Xx shape ({}) '.format(Xx.shape, self.Xx.shape))
        if Xu.shape != self.Xu.shape:
            raise Exception('Shape of Xu ({}) not compatible with KF Xu shape ({}) '.format(Xu.shape, self.Xu.shape))
        if Yx.shape != self.Yx.shape:
            raise Exception('Shape of Yx ({}) not compatible with KF Yx shape ({}) '.format(Yx.shape, self.Yx.shape))
        if Yu.shape != self.Yu.shape:
            raise Exception('Shape of Yu ({}) not compatible with KF Yu shape ({}) '.format(Yu.shape, self.Yu.shape))
        self.Xx.iloc[:,:]=Xx
        self.Xu.iloc[:,:]=Xu
        self.Yx.iloc[:,:]=Yx
        self.Yu.iloc[:,:]=Yu

    def discretize(self,dt,method='exponential'):
        self.dt=dt
        self.Xxd,self.Xud = KFDiscretize(self.Xx, self.Xu, dt, method=method)

    def estimateTimeStep(self,u,y,x,P,Q,R):
        """
        OUTPUTS:
          z1: States at time n
          P1: Process covariance at time n
          Kk: Kalman gain
        """
        return EstimateKFTimeStep(u,y,x,self.Xxd,self.Xud,self.Yx.values,self.Yu.values,P,Q,R)

    def covariancesFromSig(self):
        if not hasattr(self,'sigX'):
            raise Exception('Set `sigX` before calling `covariancesFromSig` (e.g. `sigmasFromClean`)')
        if not hasattr(self,'sigY'):
            raise Exception('Set `sigY` before calling `covariancesFromSig` (e.g. `sigmasFromClean`)')

        for lab in self.sX:
            if self.sigX[lab]==0:
                print('[WARN] Sigma for x[{}] is zero, replaced by 1e-4'.format(lab))
                self.sigX[lab]=1e-4
        for lab in self.sY:
            if self.sigY[lab]==0:
                print('[WARN] Sigma for y[{}] is zero, replaced by 1e-4'.format(lab))
                self.sigY[lab]=1e-4

        P = np.eye(self.nX)
        Q = np.diag([self.sigX[lab]**2 for lab in self.sX])
        R = np.diag([self.sigY[lab]**2 for lab in self.sY])
        return P,Q,R



    # --------------------------------------------------------------------------------}
    # --- TIME, Optional convenient methods if a time vector is already available
    # --------------------------------------------------------------------------------{
    def setTimeVec(self,time):
        self.time=time

    @property
    def nt(self):
        return len(self.time)

    def setCleanValues(self,df,ColMap=None):
        if ColMap is None:
            ColMap=dict()
            for k in df.columns.values:
                ColMap[k]=k

        # --- Defining "clean" values 
        self.X_clean = pd.DataFrame(data=np.zeros((self.nt,self.nX)), columns=self.sX)
        self.Y_clean = pd.DataFrame(data=np.zeros((self.nt,self.nY)), columns=self.sY)
        self.U_clean = pd.DataFrame(data=np.zeros((self.nt,self.nU)), columns=self.sU)
        self.S_clean = pd.DataFrame(data=np.zeros((self.nt,self.nS)), columns=self.sS)
        for i,lab in enumerate(self.sX):
            try:
                self.X_clean[lab]=df[ColMap[lab]].values
            except:
                print('[WARN] Clean state not available      :', lab)

        for i,lab in enumerate(self.sY):
            try:
                self.Y_clean[lab]=df[ColMap[lab]].values
            except:
                print('[WARN] Clean measurement not available:', lab)
        for i,lab in enumerate(self.sU):
            try:
                self.U_clean[lab] =df[ColMap[lab]].values
            except:
                print('[WARN] Clean output not available     :', lab)
        for i,lab in enumerate(self.sS):
            try:
                self.S_clean[lab] =df[ColMap[lab]].values
            except:
                print('[WARN] Clean misc var not available   :', lab)

    def setY(self,df,ColMap=None):
        if ColMap is None:
            ColMap=dict()
            for k in df.columns.values:
                ColMap[k]=k

        for i,lab in enumerate(self.sY):
            self.Y[lab]=df[ColMap[lab]]

    def initTimeStorage(self):
        self.X_hat = pd.DataFrame(data=np.zeros((self.nt, self.nX)), columns=self.sX)
        self.Y_hat = pd.DataFrame(data=np.zeros((self.nt, self.nY)), columns=self.sY)
        self.Y     = pd.DataFrame(data=np.zeros((self.nt, self.nY)), columns=self.sY)
        self.S_hat = pd.DataFrame(data=np.zeros((self.nt, self.nS)), columns=self.sS)
    
#     # TODO use property or dict syntax
#     def get_vY(self,lab):
#         return self.Y[lab].values
#     def set_vY(self, lab, val ):
#         self.Y[lab]=val
# 
#     def get_vX_hat(self,lab):
#         return self.X_hat[lab].values
#     def set_vX_hat(self, lab, val ):
#         self.X_hat[lab]=val
# 
#     def get_Y(self,lab,it):
#         return self.Y[lab][it]
#     def set_Y(self, lab, val ):
#         self.Y[self.iY[lab],it]=val
# 
#     def get_X_hat(self,lab,it):
#         return self.X_hat[self.iX[lab],it]
#     def set_X_hat(self,lab, val ):
#         self.X_hat[self.iX[lab],it]=val


    def initFromClean(self):
        x = self.X_clean.iloc[0,:].values
        # x = np.zeros(nX)
        self.X_hat.iloc[0,:] = x
        self.Y_hat.iloc[0,:] = self.Y_clean.iloc[0,:]
        return x

    def initZero(self):
        return np.zeros(self.nX)


    def initFromSimulation(KF, measFile, nUnderSamp=1, tRange=None, colMap=None, timeCol='Time_[s]', dataDict=None):
        """" 
         - Open a simulation result file
         - Use dt to discretize the KF
         - Define clean values of measurements and states based on simulation
         - Define sigmas from the std of the clean signals

         dataDict: additional data provided for ColMap
        
        """
        import welib.fast.fastlib as fastlib 
        import welib.weio as weio

        # --- Loading "Measurements"
        df=weio.read(measFile).toDataFrame()
        df=df.iloc[::nUnderSamp,:]                      # reducing sampling
        if tRange is not None:
            df=df[(df[timeCol]>= tRange[0]) & (df[timeCol]<= tRange[1])] # reducing time range
        time = df[timeCol].values
        dt   = (time[-1] - time[0])/(len(time)-1)
        KF.df = fastlib.remap_df(df, colMap, bColKeepNewOnly=False, dataDict=dataDict)

        # --- 
        KF.discretize(dt, method='exponential')
        KF.setTimeVec(time)
        KF.setCleanValues(KF.df)

        # --- Estimate sigmas from measurements
        KF.sigX_c,KF.sigY_c = KF.sigmasFromClean(factor=1)

    def prepareTimeStepping(KF):
        # --- Process and measurement covariances
        KF.P, KF.Q, KF.R = KF.covariancesFromSig()
        # --- Storage for plot
        KF.initTimeStorage()


    def setYFromClean(self,NoiseRFactor=None,y_bias=None,R=None):
        """ 
        Create y vector from "clean" y values (when available with simulations for instance)
        Possibilty to add a constant bias

          y_bias: nY vector of bias for each measurements

          R : covariance matrix
          NoiseRFactor : factor for covariance matrix

        """
        if y_bias is None:
            y_bias = np.zeros(self.nY)

        if NoiseRFactor is not None:
            Ey = np.sqrt(R)*NoiseRFactor

        for it in range(0,self.nt):    
            self.Y.iloc[it,:] = self.Y_clean.iloc[it,:] + np.dot(Ey,np.random.randn(self.nY,1)).ravel() + y_bias

    def sigmasFromClean(self,factor=1):
        sigX   = dict()
        for iX,lab in enumerate(self.sX):
            std = np.std(self.X_clean[lab])
            if std==0:
                res=1
            else:
                res=10**(np.floor(np.log10(std))-1)
            sigX[lab]=np.floor(std/res)*res  * factor
        sigY   = dict()
        for iY,lab in enumerate(self.sY):
            std = np.std(self.Y_clean[lab])
            if std==0:
                res=1
            else:
                res=10**(np.floor(np.log10(std))-1)
            sigY[lab]=np.floor(std/res)*res * factor
        self.sigX=sigX.copy()
        self.sigY=sigY.copy()
        return sigX,sigY

    def print_sigmas(self):
        sigX_c=self.sigX_c
        sigY_c=self.sigY_c
        if sigX_c is not None:
            print('Sigma X            to be used     from inputs')
            for k,v in self.sigX.items():
                print('Sigma {:10s}: {:12.3f}  {:12.3f}'.format(k,v,sigX_c[k]))
        else:
            print('Sigma X            to be used')
            for k,v in self.sigX.items():
                print('Sigma {:10s}: {:12.3f}'.format(k,v))

        if sigY_c is not None:
            print('Sigma Y            to be used     from inputs')
            for k,v in self.sigY.items():
                print('Sigma {:10s}: {:12.3f}  {:12.3f}'.format(k,v,sigY_c[k]))
        else:
            print('Sigma Y            to be used')
            for k,v in self.sigY.items():
                print('Sigma {:10s}: {:12.3f}'.format(k,v))


    # --------------------------------------------------------------------------------}
    # --- Plot functions 
    # --------------------------------------------------------------------------------{

    def plot_X(KF, title='States X', **kwargs):
        return _plot(KF.time, KF.X_clean, KF.X_hat, KF.sX, title=title, **kwargs)

    def plot_Y(KF, title='Measurements Y', **kwargs):
        return _plot(KF.time, KF.Y_clean, KF.Y_hat, KF.sY, title=title, X_noisy=KF.Y, **kwargs)

    def plot_S(KF, title='Stored Values S',**kwargs):
        if KF.nS==0:
            return
        return _plot(KF.time, KF.S_clean, KF.S_hat, KF.sS, title=title, **kwargs)

    def save(KF, filename, fmt='pickle'):
        if fmt=='pickle':
            import pickle
            with open(filename,'wb') as f:
                pickle.dump(self,f)
        else:
            raise NotImplementedError()

    @staticmethod
    def load(filename):
        ext = os.path.splitext(filename)[1].lower()
        if ext=='.pkl':
            import pickle
            with open(filename,'rb') as f:
                dat=pickle.load(f)
        else:
            raise NotImplementedError()
        return dat


def _plot(time, X_clean, X_hat, sX, title='', X_noisy=None, fig=None, COLRS=None, channels=None, nPlotCols=1):
    import matplotlib
    import matplotlib.pyplot as plt
    # --- Compare States
    if COLRS is None:
        cmap = matplotlib.cm.get_cmap('viridis')
        COLRS = [(cmap(v)[0],cmap(v)[1],cmap(v)[2]) for v in np.linspace(0,1,3+1)]

    if channels is not None:
        sX=list(sX)
        I=[]
        for s in channels:
            try:
               I.append(sX.index(s))
            except:
                print('[FAIL] Signal {} not found '.format(s))
        if len(I)==0:
            I=np.arange(len(sX))
    else:
        I=np.arange(len(sX))

    if fig is None:

        if nPlotCols==2:
            fig,axes = plt.subplots(int(np.ceil(len(I)/2)), 2, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.07, right=0.98, top=0.955, bottom=0.05, hspace=0.20, wspace=0.20)
        else:
            fig,axes = plt.subplots(len(I), 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)



        if not hasattr(axes,'__len__'):
            axes=[axes]
    axes=(np.asarray(axes).T).ravel()
    
    for j,i in enumerate(I):
        s  = sX[i]
        ax = axes[j]
        ax.plot(time,X_clean[s],''  , color=COLRS[0],label='Reference')
        if X_noisy is not None:
            ax.plot(time,X_noisy[s],'-.',  color=COLRS[2] ,label='Noisy')
        ax.plot(time,X_hat[s],'--', color=COLRS[1],label='Estimate')
        ax.set_ylabel(s)
        ax.tick_params(direction='in')
    axes[0].set_title(title)
    axes[-1].set_xlabel('Time [s]')
    axes[-1].legend()
    return fig, axes


