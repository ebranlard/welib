# TODO Change notations
""" """
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
# Local
from .kalman import *
import welib.fast.fastlib as fastlib 
import welib.weio as weio
from welib.tools.strings import OK, INFO, WARN, FAIL, NOTE
from welib.tools.strings import prettyMat


def pretty_PrintMat(M,fmt='{:11.3e}',fmt_int='    {:4d}   ',sindent='   '):
    s = prettyMat(M, var=None, digits=2, nchar=None, sindent='   ', align='right', center0=True, newline=True, openChar='[',closeChar=']', sepChar=' ', xmin=1e-16)
    return s


class KalmanFilter(object):
    def __init__(self, sX0=None, sXa=None, sU=None, sY=None, sS=None, sXd=None, KM=None):
        # State names 
        self.sX0, self.sXa, self.sU, self.sY, self.sS, self.sXd = None, None, None, None, None, None

        # State matrices
        self.Xx, self.Xu, self.Yx, self.Yu = None, None, None, None

        # Standard deviations and covariance matrix
        self.sigX_c = None
        self.sigY_c = None
        self.sigX   = None
        self.sigY   = None
        self.P = None
        self.Q = None
        self.R = None
        self.R_c = None # From Measurements

        # Time
        self.time = None
        self.dt   = None
        self.it   = 0

        # Time storage (Dataframes)
        self.X_hat    = None # pd.DataFrame(data = np.zeros((self.nt, self.nX)), columns = self.sX) # Estimated state
        self.Y_hat    = None # pd.DataFrame(data = np.zeros((self.nt, self.nY)), columns = self.sY) # Estimate output / measurement
        self.Y        = None # pd.DataFrame(data = np.zeros((self.nt, self.nY)), columns = self.sY) # Actual measurement, with potential noise
        self.S_hat    = None # pd.DataFrame(data = np.zeros((self.nt, self.nS)), columns = self.sS)
        self.U_hat    = None # pd.DataFrame(data = np.zeros((self.nt, self.nU)), columns = self.sU)
        self.XD_hat   = None # pd.DataFrame(data = np.zeros((self.nt, self.nX)), columns = self.sXd)
        self.X_clean  = None # pd.DataFrame(data = np.zeros((self.nt,self.nX)), columns = self.sX)
        self.Y_clean  = None # pd.DataFrame(data = np.zeros((self.nt,self.nY)), columns = self.sY) # Measurement without noise
        self.U_clean  = None # pd.DataFrame(data = np.zeros((self.nt,self.nU)), columns = self.sU)
        self.S_clean  = None # pd.DataFrame(data = np.zeros((self.nt,self.nS)), columns = self.sS)
        self.XD_clean = None # pd.DataFrame(data = np.zeros((self.nt,self.nX)), columns = self.sXd)
        self.Pt       = None # np.zeros((self.nt, self.nX, self.nY))  # P is nx * nx
        self.Kt       = None # np.zeros((self.nt, self.nX, self.nY))  # K is nx * ny

        # --- Actually initialization
        if KM is not None:
            self.sX0 = KM.sQ
            self.sXa = KM.sQa
            self.sU  = KM.sU
            self.sY  = KM.sY
            self.sS  = KM.sS
            self.sXd = KM.sQd
            self.KM = KM
        
        else:
            self.sX0 = sX0
            self.sXa = sXa
            self.sU  = sU
            self.sY  = sY
            self.sS  = sS # Storage, "Misc" values
            self.sXd = sXd # Storage, "Misc" values
            self.KM = None

        #  State vector is States and Augmented states
        self.sX = np.concatenate((self.sX0, self.sXa))

        if self.sS is None :
            self.sS = []
        if self.sXd is None:
            self.sXd = ['d' + c for c in self.sX] # NOTE: might have duplication...

        # --- Defining index map for convenience
        self.iX = {lab: i   for i,lab in enumerate(self.sX)}
        self.iY = {lab: i   for i,lab in enumerate(self.sY)}
        self.iU = {lab: i   for i,lab in enumerate(self.sU)}
        self.iS = {lab: i   for i,lab in enumerate(self.sS)}
        # --- Define empty (nan) sigmas
        self._set_empty_sigs() 

        if KM is not None:
            self.setMat(KM.A, KM.B, KM.C, KM.D)


    @property
    def nX(self): return len(self.sX)
    @property
    def nY(self): return len(self.sY)
    @property
    def nU(self): return len(self.sU)
    @property
    def nP(self): return len(self.sXa)
    @property
    def nX0(self): return len(self.sX0)
    @property
    def nS(self): return len(self.sS)
    @property
    def A(self): return self.Xx
    @property
    def B(self): return self.Xu
    @property
    def C(self): return self.Yx
    @property
    def D(self): return self.Yu

    def __repr__(self):
        s='<{} object> with attributes:\n'.format(type(self).__name__)
        s+='  sX  : {} \n'.format(self.sX)
        s+='  sX0 : {} \n'.format(self.sX0)
        s+='  sX1 : {} \n'.format(self.sXa)
        s+='  sU  : {} \n'.format(self.sU)
        s+='  sY  : {} \n'.format(self.sY)
        s+='  sS  : {} \n'.format(self.sS)
        if self.Xx is not None:
            s+=' A: State-State Matrix (Xx) \n'
            s+=pretty_PrintMat(self.Xx)+'\n'
        if self.Xu is not None:
            s+=' B: State-Input Matrix (Xu) \n'
            s+=pretty_PrintMat(self.Xu)+'\n'
        if self.Yx is not None:
            s+=' C: Output-State Matrix (Yx) \n'
            s+=pretty_PrintMat(self.Yx)+'\n'
        if self.Yu is not None:
            s+=' D: Output-Input Matrix (Yu) \n'
            s+=pretty_PrintMat(self.Yu)+'\n'
        if self.P is not None:
            s+=' P: error covariance matrix\n'
            s+=pretty_PrintMat(self.P)+'\n'
        if self.Q is not None:
            s+=' Q: process noise\n'
            s+=pretty_PrintMat(self.Q)+'\n'
        if self.R is not None:
            s+=' R: measurement matrix\n'
            s+=pretty_PrintMat(self.R)+'\n'
        return s





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
        self.Xx.iloc[:,:] = Xx
        self.Xu.iloc[:,:] = Xu
        self.Yx.iloc[:,:] = Yx
        self.Yu.iloc[:,:] = Yu

        if np.any(np.isnan(Xx)): raise Exception('A matrix contains nan')
        if np.any(np.isnan(Xu)): raise Exception('B matrix contains nan')
        if np.any(np.isnan(Yx)): raise Exception('C matrix contains nan')
        if np.any(np.isnan(Yu)): raise Exception('D matrix contains nan')


    def _set_empty_sigs(self):
        self.sigX   = {s: np.nan for s in self.sX}
        self.sigY   = {s: np.nan for s in self.sY}
        self.sigQ   = {s: np.nan for s in self.sX}
        self.sigX_c = {s: np.nan for s in self.sX}
        self.sigY_c = {s: np.nan for s in self.sY}
        self.sigQ_c = {s: np.nan for s in self.sX}



    def checkObservability(self):
        import control
        O   = control.obsv(self.Xx, self.Yx)
        try:
            sys = control.StateSpace(self.Xx, self.Xu, self.Yx, self.Yu)
        except:
            FAIL('State space')
            pass
        try:
            Wc = control.gram(sys, 'c')
        except:
            FAIL('gramian Controlability')
            pass
        try:
            FAIL('gramian Observability')
            Wo  = control.gram(sys, 'o')
        except:
            pass




    def discretize(self, dt, method='exponential'):
        self.dt = dt
        self.Xxd, self.Xud = KFDiscretize(self.Xx, self.Xu, dt, method=method)

    def estimateTimeStep(self, u, y, x, P=None, Q=None, R=None, it=None):
        """
        OUTPUTS:
          z1: States at time n
          P1: Process covariance at time n
          Kk: Kalman gain
        """
        if it is None:
            it = self.it
        if Q is None:
            Q = self.Q
        if R is None:
            R = self.R
        if P is None:
            P = self.P

        x_new, P_new, Kk = EstimateKFTimeStep(u, y, x, self.Xxd, self.Xud, self.Yx.values,self.Yu.values, P, Q, R)

        self.P = P_new
        self.it = it+1
        
        # --- Store
        self.X_hat .iloc[it+1,:] = x_new
        self.XD_hat.iloc[it+1,:] = np.dot(self.A, x_new) + np.dot(self.B, u)
        self.Y_hat .iloc[it+1,:] = np.dot(self.C, x_new) + np.dot(self.D, u)

        self.Pt[it+1, :, :]   = P_new
        self.Kt[it+1, :, :]   = Kk

        return x_new, P_new, Kk



    # --------------------------------------------------------------------------------}
    # --- TIME, Optional convenient methods if a time vector is already available
    # --------------------------------------------------------------------------------{
    def setTimeVec(self, time):
        self.time = time

    @property
    def nt(self):
        return len(self.time)

    def setCleanValues(self, df, colMap=None, verbose=False):
        if colMap is None:
            colMap=dict()
            for k in df.columns.values:
                colMap[k]=k

        # --- Defining "clean" values 
        self.X_clean = pd.DataFrame(data=np.zeros((self.nt,self.nX)), columns=self.sX)
        self.Y_clean = pd.DataFrame(data=np.zeros((self.nt,self.nY)), columns=self.sY)
        self.U_clean = pd.DataFrame(data=np.zeros((self.nt,self.nU)), columns=self.sU)
        self.S_clean = pd.DataFrame(data=np.zeros((self.nt,self.nS)), columns=self.sS)
        self.XD_clean = pd.DataFrame(data=np.zeros((self.nt,self.nX)), columns=self.sXd)
        for i,lab in enumerate(self.sX):
            try:
                self.X_clean[lab]=df[colMap[lab]].values
            except:
                if verbose:
                    print('[WARN] Clean state not available      :', lab)
        for i,lab in enumerate(self.XD_clean.columns):
            try:
                self.XD_clean[lab]=df[colMap[lab]].values
            except:
                if verbose:
                    print('[WARN] Clean state not available      :', lab)

        for i,lab in enumerate(self.sY):
            try:
                self.Y_clean[lab]=df[colMap[lab]].values
            except:
                if verbose:
                    print('[WARN] Clean measurement not available:', lab)
        for i,lab in enumerate(self.sU):
            try:
                self.U_clean[lab] =df[colMap[lab]].values
            except:
                if verbose:
                    print('[WARN] Clean output not available     :', lab)
        for i,lab in enumerate(self.sS):
            try:
                self.S_clean[lab] =df[colMap[lab]].values
            except:
                if verbose:
                    print('[WARN] Clean misc var not available   :', lab)

    def setY(self,df,colMap=None):
        if colMap is None:
            colMap=dict()
            for k in df.columns.values:
                colMap[k]=k

        for i,lab in enumerate(self.sY):
            self.Y[lab]=df[colMap[lab]]

    def initTimeStorage(self):
        self.X_hat  = pd.DataFrame(data = np.zeros((self.nt, self.nX)), columns = self.sX)
        self.Y_hat  = pd.DataFrame(data = np.zeros((self.nt, self.nY)), columns = self.sY)
        self.Y      = pd.DataFrame(data = np.zeros((self.nt, self.nY)), columns = self.sY)
        self.S_hat  = pd.DataFrame(data = np.zeros((self.nt, self.nS)), columns = self.sS)
        self.U_hat  = pd.DataFrame(data = np.zeros((self.nt, self.nU)), columns = self.sU)
        self.XD_hat = pd.DataFrame(data = np.zeros((self.nt, self.nX)), columns = self.sXd)
        self.Pt     = np.zeros((self.nt, self.nX, self.nX))  # P is nx * nx
        self.Kt     = np.zeros((self.nt, self.nX, self.nY))  # K is nx * ny
    

    def initFromClean(self):
        """ Set initial conditions based on clean data"""
        x = self.X_clean.iloc[0,:].values.copy()
        # x = np.zeros(nX)
        self.X_hat.iloc[0,:] = x
        self.Y_hat.iloc[0,:] = self.Y_clean.iloc[0,:]
        return x

    def initZero(self):
        return np.zeros(self.nX)

    def loadMeasurements(KF, measFile, nUnderSamp=1, tRange=None, colMap=None, timeCol='Time_[s]', raiseIfAbsent=False):
        """" 
         - Open a simulation result file
         - Use dt to discretize the KF
         - Define clean values of measurements and states based on simulation
        """
        # --- Loading "Measurements"
        if isinstance(measFile, pd.DataFrame):
            df = measFile
        else:
            ext = os.path.splitext(measFile)[1]
            if not os.path.exists(measFile):
                if ext == '.outb':
                    WARN('Measurement file not found, trying with .out: {}'.format(measFile))
                    measFile = measFile.replace('.outb', '.out')
                elif ext == '.out':
                    WARN('Measurement file not found, trying with .outb: {}'.format(measFile))
                    measFile = measFile.replace('.out', '.outb')
            df=weio.read(measFile).toDataFrame()

        nUnderSamp=max(nUnderSamp,1)
        df=df.iloc[::nUnderSamp,:]                      # reducing sampling
        if tRange is not None:
            df=df[(df[timeCol]>= tRange[0]) & (df[timeCol]<= tRange[1])] # reducing time range
        time = df[timeCol].values
        dt   = (time[-1] - time[0])/(len(time)-1)
        # Remapping/scaling columns to shortname variables
        if colMap is not None:
            KF.df = fastlib.remap_df(df, colMap, bColKeepNewOnly=False, raiseIfAbsent=raiseIfAbsent)
        else:
            NOTE('Not performing any column mapping to loaded measurements')
            raise Exception() # Remove me, for checking for now
        # --- 
        KF.discretize(dt, method='exponential')
        KF.setTimeVec(time)
        KF.setCleanValues(KF.df)

        # --- Estimate sigmas from measurements
        sigY, KF.R_c = KF.sigmasYFromClean(factor=1)


    def sigmasYFromClean(self, factor=1):
        sigY   = dict()
        for iY,lab in enumerate(self.sY):
            y_data = np.asarray(self.Y_clean[lab])
            std = np.std(y_data)
            if std==0:
                res=1
            else:
                res=10**(np.floor(np.log10(std))-1)
            sigY[lab] = np.floor( std / res ) *res * factor
        # Setup clean noise
        self.R_c = np.diag([sigY[lab]**2 for lab in self.sY])
        return sigY, self.R_c

    def sigmasFromClean(self, factor=1, dt=None):
        if dt is None:
            WARN('sigmasFromClean should preferably use delta t')
        else:
            NOTE('Using dt for sigmas from Clean')
        sigX = dict()
        sigQ = dict()
        for lab in self.sX:
            x_data = np.asarray(self.X_clean[lab])
            # Handle sawtooth wrapping for angle states if present
            if lab == 'psi':
                x_data_unwrapped = np.unwrap(x_data)
                std_val = np.std(x_data)
            else:
                x_data_unwrapped = x_data
                std_val = np.std(x_data)
            
            res = 1 if std_val == 0 else 10**(np.floor(np.log10(std_val)) - 1)
            sigX[lab] = np.floor(std_val / res) * res * factor
            sigQ[lab] = sigX[lab]
            if dt is not None and len(x_data) > 1:
                # Calculate rate using unwrapped data for angle states
                dx_dt = np.diff(x_data_unwrapped) / dt
                std_rate = np.std(dx_dt)
                
                # Kinematic positions shouldn't have large process random walk
                if lab in ['x', 'y', 'z', 'phi_x', 'phi_y', 'phi_z', 'psi', 'q_FA1']:
                    pass
                    # Force process noise for kinematic position to be tiny/negligible
                    #sigQ[lab] = 1e-4
                    #sigQ[lab] = sigX[lab]
                else:
                    res_rate = 1 if std_rate == 0 else 10**(np.floor(np.log10(std_rate)) - 1)
                    sigQ[lab] = np.floor(std_rate / res_rate) * res_rate * factor
                sigQ[lab] = sigQ[lab] * dt # NOTE: we multiply by dt so that Q which is sigQ^2 is sigma^2 dt^2

        sigY   = dict()
        for iY,lab in enumerate(self.sY):
            y_data = np.asarray(self.Y_clean[lab])
            std = np.std(y_data)
            if std==0:
                res=1
            else:
                res=10**(np.floor(np.log10(std))-1)
            sigY[lab] = np.floor( std / res ) *res * factor
        # We store the clean
        self.sigX_c = sigX.copy()
        self.sigY_c = sigY.copy()
        self.sigQ_c = sigQ.copy()
        # And we store a ictionary for the user
        self.sigX   = sigX.copy()
        self.sigY   = sigY.copy()
        self.sigQ   = sigQ.copy()
        
        return sigX, sigY, sigQ


    def setupCovariances(KF, useDt=False, Pidentity=True):
        # TODO add tuning here
        dt = None
        if useDt:
            dt = KF.dt
        # --- Process and measurement covariances
        KF.P, KF.Q, KF.R = KF.covariancesFromSig(Pidentity=Pidentity, dt=dt)

    def covariancesFromSig(self, dt=None, Pidentity=True):
        # -- Safety
        has_nanX = any(np.isnan(v) for v in self.sigX.values())
        has_nanY = any(np.isnan(v) for v in self.sigY.values())
        has_nanQ = any(np.isnan(v) for v in self.sigQ.values())
        if has_nanX or has_nanY or has_nanQ:
            self.print_sigmas()
            raise Exception('covarianceFromSig, some sigs have NaN, call sigmasFromClean first for instance.')


        for lab in self.sX:
            if self.sigX[lab]==0:
                print('[WARN] Sigma for x[{}] is zero, replaced by 1e-4'.format(lab))
                self.sigX[lab]=1e-4
        for lab in self.sY:
            if self.sigY[lab]==0:
                print('[WARN] Sigma for y[{}] is zero, replaced by 1e-4'.format(lab))
                self.sigY[lab]=1e-4

        # Initial State Error Covariance P0
        if Pidentity:
            WARN('covarianceFromSig: P is set to identity')
            P = np.eye(self.nX)
        else:
            NOTE('covarianceFromSig: P is set using sigX')
            P = np.diag([self.sigX[lab]**2 for lab in self.sX])


        # Measurement Noise Covariance R
        R = np.diag([self.sigY[lab]**2 for lab in self.sY])

        # Process Noise Covariance Q
        Q = np.diag([self.sigQ[lab]**2 for lab in self.sX])
        # Consistent quadratic scaling: Var = (rate * dt)^2
        #Q = np.diag([self.sigQ[lab]**2 for lab in self.sX])
#         if dt is None:
#         else:
#             NOTE('covarianceFromSig: dt is used Q')
#             # sigQ represents rate std dev (units/s) -> step variance is (sigQ * dt)^2
#             Q = np.diag([(self.sigQ[lab] * dt)**2 for lab in self.sX])

        return P, Q, R



    def prepareTimeStepping(KF):
        # --- Storage for plot
        KF.initTimeStorage()


    def setYFromClean(self, NoiseRFactor=None, y_bias=None, R=None):
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
            if hasattr(NoiseRFactor, '__len__'):
                if len(NoiseRFactor)!=R.shape[0]:
                    raise Exception('NoiseRFactor has wrong length')
            else:
                NoiseRFactor=[NoiseRFactor] * R.shape[0]
            NoiseRFactor = np.asarray(NoiseRFactor)
            Ey = np.sqrt(R)*NoiseRFactor

        for it in range(0,self.nt):    
            self.Y.iloc[it,:] = self.Y_clean.iloc[it,:] + np.dot(Ey,np.random.randn(self.nY,1)).ravel() + y_bias


    def print_sigmas(self):
        sigX_c = getattr(self, 'sigX_c', None)
        sigY_c = getattr(self, 'sigY_c', None)
        sigQ_c = getattr(self, 'sigQ_c', None)
        sigQ   = getattr(self, 'sigQ', None)

        # Print State (X) and Process Noise (Q) Sigmas side by side
        header_x = 'Sigma X             from clean    to be used'
        header_q = ' |  Sigma Q             from clean    to be used'
        print(header_x + header_q)

        for k in self.sX:
            vx = self.sigX.get(k, np.nan)
            vxc = sigX_c.get(k, np.nan) if sigX_c is not None else np.nan
            
            line_x = 'Sigma {:10s}: {:12.3f}  {:12.3f}'.format(k, vxc, vx)
            
            vq  = sigQ.get  (k, np.nan) if sigQ is   not None else np.nan
            vqc = sigQ_c.get(k, np.nan) if sigQ_c is not None else np.nan
            
            line_q = ' |  Sigma Q {:8s}: {:12.3f}  {:12.3f}'.format(k, vqc, vq)
            print(line_x + line_q)


        # Print Measurement (Y) Sigmas side by side
        print('---')
        print('Sigma Y             from clean    to be used')
        for k in self.sY:
            vy = self.sigY.get(k, np.nan)
            vyc = sigY_c.get(k, np.nan) if sigY_c is not None else np.nan
            print('Sigma {:10s}: {:12.3f}  {:12.3f}'.format(k, vyc, vy))

        s = ''
        if getattr(self, 'Q', None) is not None:
            s += '\n Q: process noise matrix (diag):\n'
            s += pretty_PrintMat(np.diag(self.Q)) + '\n'
        if getattr(self, 'R', None) is not None:
            s += ' R: measurement matrix (diag)\n'
            s += pretty_PrintMat(np.diag(self.R)) + '\n'
        print(s)


    # --------------------------------------------------------------------------------}
    # --- Plot functions 
    # --------------------------------------------------------------------------------{

    def plot_X(KF, title='States X', **kwargs):
        """ plot states, return fig """
        return _plot(KF.time, KF.X_clean, KF.X_hat, KF.sX, title=title, **kwargs)


    def plot_U(KF, title='Inputs U', **kwargs):
        """ plot inputs, return fig """
        return _plot(KF.time, KF.U_clean, KF.U_hat, KF.sU, title=title, **kwargs)

    def plot_Y(KF, title='Measurements Y', **kwargs):
        """ plot measurements, return fig """
        return _plot(KF.time, KF.Y_clean, KF.Y_hat, KF.sY, title=title, X_noisy=KF.Y, **kwargs)

    def plot_S(KF, title='Stored Values S',**kwargs):
        """ plot stored values, return fig """
        if KF.nS==0:
            return
        return _plot(KF.time, KF.S_clean, KF.S_hat, KF.sS, title=title, **kwargs)

    def plot_K(self, title='Kalmab gain'):
        """ Plots the Kalman Gains K grouped by measurement. """
        # Kt shape: (nt, nX, nY)
        fig, axes = plt.subplots(self.nY, 1, sharex=True, figsize=(10, 4 * self.nY))
        if self.nY == 1: axes = [axes]

        for j, lab_y in enumerate(self.sY):
            ax = axes[j]
            for i, lab_x in enumerate(self.sX):
                ax.plot(self.time, self.Kt[:, i, j], label=f'Gain {lab_x}')
            
            ax.set_ylabel(f'Gain for {lab_y}')
            ax.set_title(f'Kalman Gains relative to Sensor: {lab_y}')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
            ax.grid(True, alpha=0.3)

        axes[-1].set_xlabel('Time [s]')
        plt.tight_layout()


    def plot_P(self):
        """ Plots the diagonal of the covariance matrix P over time. """
        fig, ax = plt.subplots(figsize=(10, 5))
        
        # Extract the diagonal elements for each state
        # Pt shape: (nt, nX, nX)
        for i, lab in enumerate(self.sX):
            # Taking the sqrt to plot standard deviation (more intuitive units)
            sigma = np.sqrt(self.Pt[:, i, i])
            ax.plot(self.time, sigma, label=fr'$\sigma$({lab})')
        
        ax.set_yscale('log')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Standard Deviation (log scale)')
        ax.set_title('State Uncertainty (Diagonal of P)')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
        ax.grid(True, which="both", alpha=0.3)
        plt.tight_layout()

    def plot_innovation(self):
        """
        Plots the residual (innovation) between measurements and predicted output.
        """
        # Calculate residuals
        residuals = self.Y - self.Y_hat
        
        fig, axes = plt.subplots(self.nY, 1, sharex=True, figsize=(10, 4 * self.nY))
        if self.nY == 1: axes = [axes]
        
        for j, lab_y in enumerate(self.sY):
            ax = axes[j]
            # Plotting the residual
            ax.plot(self.time, residuals[lab_y], label='Innovation (Y - Y_hat)', color='red', alpha=0.7)
            ax.axhline(0, color='black', linestyle='--', linewidth=1)
            
            ax.set_ylabel(f'Residual [{lab_y}]')
            ax.set_title(f'Diagnostic: Innovation for {lab_y}')
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3)

        axes[-1].set_xlabel('Time [s]')
        plt.tight_layout()



    def save(self, filename, fmt='pickle'):
        if fmt=='pickle':
            import pickle
            with open(filename,'wb') as f:
                pickle.dump(self,f)
        else:
            raise NotImplementedError()

    def saveOutputs(KF, filename, fmt='outb', df=None):

        if df is None:
            df = KF.toDataFrame()

        if fmt=='csv':
            pass
        elif fmt=='outb':
            from welib.weio.fast_output_file import writeDataFrame
            writeDataFrame(df, filename, binary=True)
        else:
            raise NotImplementedError()

        return df


    def toDataFrame(KF):
        """ Concatenante all info into a dataframe """

        def splitunit(s):
            iu=s.rfind('_[')
            if iu>1:
                return s[:iu], s[iu:]
            else:
                return s, ''

        def cleancol(l):
            return ['_clean'.join(splitunit(c)) for c in l]

        cols = []
        cols += list(KF.X_hat.columns)
        cols += cleancol(KF.X_clean.columns)
        cols += list(KF.U_hat.columns)
        cols += cleancol(KF.U_clean.columns)
        cols += ['OUT_'+c for c in list(KF.Y_hat.columns)]
        cols += ['OUT_'+c for c in cleancol(KF.Y_clean.columns)]
        cols += list(KF.S_hat.columns)
        cols += cleancol(KF.S_clean.columns)
        df = pd.concat((KF.X_hat, KF.X_clean, KF.U_hat, KF.U_clean, KF.Y_hat, KF.Y_clean, KF.S_hat, KF.S_clean), axis=1)
        df.columns = cols

        # "Accelerations" 
        dfAcc = pd.concat((KF.XD_hat, KF.XD_clean), axis=1)
        cols = list(KF.XD_hat.columns)
        cols += cleancol(KF.XD_clean.columns)
        dfAcc.columns=cols
        # We keep only the part that is not already in the state vector
        col_new = [c for c in cols if c not in df.columns]
        dfAcc=dfAcc[col_new] 

        # 
        df = pd.concat((df,dfAcc), axis=1)


        df.insert(0, 'Time_[s]', KF.time)
        return df


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


def _plot(time, X_clean, X_hat, sX, title='', X_noisy=None, fig=None, COLRS=None, channels=None, nPlotCols=1, figSize=(6.4,4.8), 
          stats='sigRatio,eps,R2', tRangeStats=None, printStats=False, statsDict=None,
          STY=None, refLast=False):
    from welib.tools.stats import comparison_stats
    from welib.tools.strings import latexStrip
    from welib.tools.colors import cmap_colors
    if statsDict is None:
        statsDict={}
    # --- Misc inits
    if COLRS is None:
        COLRS = cmap_colors(4, 'viridis')
    if STY is None:
        STY=['-','--','-','--']
    if printStats:
        import sys
        if sys.platform.startswith('win'):
            sys.stdout.reconfigure(encoding='utf-8')
            sys.stderr.reconfigure(encoding='utf-8')


    # --- Channel indices
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

    # --- Time indices for stats
    if tRangeStats is None:
        IT = np.arange(0, len(time))
    else:
        tRangeStats[1] = min( max(time), tRangeStats[1])
        IT = np.logical_and(time>tRangeStats[0], time<tRangeStats[1])
        if len(IT)==0:
            IT = np.arange(0, len(time))



    if fig is None:

        if nPlotCols==2:
            fig,axes = plt.subplots(int(np.ceil(len(I)/2)), 2, sharex=True, figsize=figSize) # (6.4,4.8)
            fig.subplots_adjust(left=0.07, right=0.98, top=0.955, bottom=0.05, hspace=0.20, wspace=0.20)
        else:
            fig,axes = plt.subplots(len(I), 1, sharex=True, figsize=figSize) # (6.4,4.8)
            fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)

        if not hasattr(axes,'__len__'):
            axes=[axes]
    else:
        axes = fig.axes
    axes=(np.asarray(axes).T).ravel()
    
    # --- Plot request signals
    for j,i in enumerate(I):
        s  = sX[i]
        ax = axes[j]
        if not refLast:
            if X_clean is not None:
                ax.plot(time,X_clean[s],STY[0]  , color=COLRS[0],label='Reference')
        if X_noisy is not None:
            ax.plot(time,X_noisy[s],STY[2],  color=COLRS[2] ,label='Noisy')
        ax.plot(time,X_hat[s],STY[1], color=COLRS[1],label='Estimate')
        if refLast:
            if X_clean is not None:
                ax.plot(time,X_clean[s],STY[0]  , color=COLRS[0],label='Reference')

        if stats:
            if X_clean is not None:
                statsDict[s], sStats = comparison_stats(time[IT], X_clean[s][IT], time[IT], X_hat[s][IT], stats=stats, method='1-2')
                Ylim = ax.get_ylim()
                Xlim = ax.get_xlim()
                ax.text(Xlim[0]+(Xlim[1]-Xlim[0])*0.01 ,Ylim[0]+(Ylim[1]-Ylim[0])*0.82, sStats, fontsize=10)
                if printStats:
                    print(f"{s:10s} "+latexStrip(sStats))

        if tRangeStats is not None:
            ax.axvline(x = tRangeStats[0], color = 'k', ls='--')
            ax.axvline(x = tRangeStats[1], color = 'k', ls='--')

        ax.set_ylabel(s)
        ax.tick_params(direction='in')


    axes[0].set_title(title)
    axes[-1].set_xlabel('Time [s]')
    axes[len(I)-1].legend()
    # Remove unnecessary axes
    if len(axes)>len(I):
        for j in range(len(I), len(axes)):
            axes[j].axis('off')
    return fig


