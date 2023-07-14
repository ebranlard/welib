import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.optimize import minimize_scalar

from welib.weio.fast_input_deck import FASTInputDeck
from welib.tools.signal_analysis import zero_crossings

# ---
def interp2d_pairs(*args,**kwargs):
    """ Same interface as interp2d but the returned interpolant will evaluate its inputs as pairs of values.
    Inputs can therefore be arrays

    example:
       f = interp2d_pairs(vx, vy, M, kind='cubic')

    vx: array of length nx
    vy: array of length ny
    M : array of shape nx x ny
    f : interpolant function
          v = f(x,y) : if x,y are array of length n, v is of length n
                       with  v_i = f(x_i, y_i)
    author: E. Branlard
    """
    # Internal function, that evaluates pairs of values, output has the same shape as input
    def interpolant(x,y,f):
        x,y = np.asarray(x), np.asarray(y)
        return (si.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x.ravel(), y.ravel())[0]).reshape(x.shape)
    # Wrapping the scipy interp2 function to call out interpolant instead
    return lambda x,y: interpolant(x,y,si.interp2d(*args,**kwargs))


def Paero(WS, Pitch, Omega, R, rho, fCP):
    """ Taero returns the aerodynamic power
         - Pitch  [deg]
         - Omega [rad/s]
         - R : the blade radius [m]
         - fCP : an interpolant for CP(Pitch,lambda) as returned by interp2d_paris
         - rho : the air density [kg/m^3]
    """
    Lambda = Omega * R / WS
    CP     = fCP(Pitch, Lambda)
    P      = 1/2*rho*np.pi*R**2*WS**3*CP
    return P


def Qaero(WS, Pitch, Omega, R, rho, fCP):
    """ Qaero returns the aerodynamic torque
         - Pitch [deg]
         - Omega [rad/s]
         - R : the blade radius
         - fCP : an interpolant for CP(Pitch,lambda)
         - rho : the air density
    """
    Pitch = np.asarray(Pitch)
    WS    = np.asarray(WS)
    Omega = np.asarray(Omega)
    Lambda = Omega * R / WS
    CP = fCP(Pitch, Lambda)
    Q = 1/2*rho*np.pi*R**2*WS**3/Omega*CP
    return Q

def Taero(WS, Pitch, Omega, R, rho, fCT):
    """ Taero returns the aerodynamic thrust of a given turbine
         - Pitch [deg]
         - Omega [rad/s]
         - R : the blade radius
         - fCP : an interpolant for CP(Pitch,lambda)
         - rho : the air density
    """
    Lambda = Omega * R / WS
    CT = fCT(Pitch,Lambda)
    T = 1/2*rho*np.pi*R**2*WS**2*CT
    return T


class TabulatedWSEstimatorBase():


    def __init__(self, R=None, rho=1.225, fstFile=None):
        """ 
        INPUTS:
          either:
          - fstFile: FAST Input file, used to obtain rotor radius and airdensity
          or
          - rho : air density [kg/m^3]
          - R       : rotor radius [m]
        """
        # --- Data
        self.WSmax = 35
        self.CP    = None
        self.CT    = None
        self.OP    = None

        # ---
        if fstFile:
            fst = FASTInputDeck(fstFile)
            R       = fst.ED['TipRad']
            if fst.AD is None:
                raise Exception('AeroDyn file not read but needed for wind speed estimator, while reading {}'.format(fstFile))
            rho_AD = fst.AD['AirDens']
            try:
                rho_main = fst.fst['AirDens']
            except:
                rho_main = rho_AD
            if isinstance(rho_AD, str):
                rho_AD = rho_main

        self.fstFile  = fstFile
        self.R       = R
        self.rho = rho

    def _sanitizeOP(self, OP, expectedCols=None, onlyExpected=True):
        if expectedCols is None:
            expectedCols = ['WS_[m/s]', 'Pitch_[deg]', 'RotSpeed_[rpm]']

        # --- Trying to be nice about column names
        OP.columns = [c.lower().replace(' ','_').replace('(','[').replace(')',']') for c in OP.columns]

        d =dict([(k, 'WS_[m/s]') for k in ['ws_[m/s]', 'ws']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'Pitch_[deg]') for k in ['bldpitch1_[deg]', 'bldpitch_[deg]','pitch_[deg]', 'pitch']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'RotSpeed_[rpm]') for k in ['rotspeed_[rpm]', 'rpm', 'rpm_[rpm]', 'omega_[rpm]']])
        OP.rename(columns = d, inplace=True)

        if 'rtaeromxh_[kn-m]' in OP.keys(): # TODO standardize Units WE
            OP['rtaeromxh_[n-m]'] =  OP['rtaeromxh_[kn-m]'].values*1000
        d =dict([(k, 'Qaero_[Nm]') for k in ['rtaeromxh_[n-m]']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'PhiY_[deg]') for k in ['phiy_[deg]']])
        OP.rename(columns = d, inplace=True)

        for c in expectedCols:
            if c not in OP:
                print('>>> Columns', OP.keys())
                raise Exception('OP is missing: {}'.format(c))

        if onlyExpected:
            OP = OP[expectedCols]

        return OP


class TabulatedWSEstimator(TabulatedWSEstimatorBase):

    def __init__(self, R=None, rho=1.225, fstFile=None, basename=None, operFile=None, aeroMapFile=None, OmegaLow=0, OmegaRated=10):
        """ 
        INPUTS:
          either:
          - fstFile: FAST Input file, used to obtain rotor radius and airdensity
          or
          - rho : air density [kg/m^3]
          - R       : rotor radius [m]
        """
        # Initialize parent class
        TabulatedWSEstimatorBase.__init__(self, R, rho, fstFile)

        if basename is not None:
            aeroMapFile = basename+'_CPCTCQ.txt'
            operFile    = basename+'_Oper.csv'

        # --- DATA
        self.Pitch    = None
        self.Lambda   = None
        # Operating condition
        self.OmegaLow   = OmegaLow
        self.OmegaRated = OmegaRated
        self.OP     = None
        # Files
        self.operFile = operFile
        self.aeroMapFile = aeroMapFile

        if operFile is not None:
            self.loadOper(operFile)

        if aeroMapFile is not None:
            self.loadAeroMap(aeroMapFile)


    def loadOper(self, operFile):
        if not os.path.exists(operFile):
            print('[WARN] Operating point file not found: ',operFile)
        else:
            #print('>>> Loading oper file: ',operFile)
            import welib.weio as weio
            OP = weio.read(operFile).toDataFrame()

            OP = self._sanitizeOP(OP, ['WS_[m/s]', 'Pitch_[deg]', 'RotSpeed_[rpm]', 'Qaero_[Nm]'], onlyExpected=False)

            self.WS   =OP['WS_[m/s]'].values
            self.Omega=OP['RotSpeed_[rpm]'].values*2*np.pi/60
            self.OmegaRated=np.max(self.Omega)
            self.OmegaLow  =0.4*self.OmegaRated
            self.WSRated=np.interp(self.OmegaRated*0.98,self.Omega,self.WS)
            self.WSCutOff=28
            self.OP=OP
            self.operFile = operFile

    def loadAeroMap(self, aeroMapFile):
        """ Load file containing aeromap: CP(Lambda, Pitch) """
        #print('>>>> Loading aeroMapFile', aeroMapFile)
        # TODO more file formats
        from welib.weio.rosco_performance_file import ROSCOPerformanceFile
        rs = ROSCOPerformanceFile(aeroMapFile)
        self.Pitch  = rs['pitch']
        self.Lambda = rs['TSR']
        self.CP     = rs['CP']
        self.CT     = rs['CT']
        # Trigger
        self.aeroMapFile = aeroMapFile
        self.computeWeights()

    def loadFromBasename(self, basename=None, suffix=''):
        """  """
        aeroMapFile = basename+'_CPCTCQ'+suffix+'.txt'
        operFile      = basename+'_Oper'+suffix+'.csv'

        if os.path.exists(aeroMapFile):
            self.loadAeroMap(aeroMapFile)
        else:
            # Old format
            LambdaFile  = basename + '_Lambda'+suffix+'.csv'
            PitchFile   = basename + '_Pitch'+suffix+'.csv'
            CPFile      = basename + '_CP'+suffix+'.csv'
            CTFile      = basename + '_CT'+suffix+'.csv'
            self.Pitch  = pd.read_csv(PitchFile ,header = None).values.ravel()
            self.Lambda = pd.read_csv(LambdaFile,header = None).values.ravel()
            self.CP     = pd.read_csv(CPFile,header     = None).values
            self.CP[self.CP<=0]=0
            self.CT     = pd.read_csv(CTFile,header     = None).values
            self.CT[self.CT<=0]=0
            # Trigger
            self.computeWeights()

        if os.path.exists(operFile):
            self.loadOper(operFile)

    def computeWeights(self):
        # Compute interpolants
        self.CP[self.CP<=0]=0
        self.fCP = interp2d_pairs(self.Pitch,self.Lambda,self.CP,kind='cubic')
        if self.CT is not None:
            self.CT[self.CT<=0]=0
            self.fCT = interp2d_pairs(self.Pitch,self.Lambda,self.CT,kind='cubic')
        else:
            self.fCT = None

        # --- Lambda at higher res
        LambdaMid = self.Lambda[:-1] + np.diff(self.Lambda)/2
        Lambda = np.sort(np.concatenate((self.Lambda, LambdaMid)))
        self.LambdaHR = Lambda

    def Power(self,WS,Pitch,Omega):
        return Paero(WS, Pitch, Omega, self.R, self.rho, self.fCP)

    def Thrust(self,WS,Pitch,Omega):
        return Taero(WS, Pitch, Omega, self.R, self.rho, self.fCT)

    def Torque(self,WS,Pitch,Omega):
        return Qaero(WS, Pitch, Omega, self.R, self.rho, self.fCP)

    def TorqueAt(self, Pitch, Omega):
        """ 
        Return Torque(WS) curve for a given pitch and rotational speed
        Pitch,Omega: scalar
        """
        WS     = Omega * self.R / self.LambdaHR[-1::-1] # NOTE using LambdaHR to benefit from cubic interpolation
        WS = WS[WS<self.WSmax]
        vPitch = np.array([Pitch]*len(WS))
        vOmega = np.array([Omega]*len(WS))
        return WS, self.Torque(WS, vPitch, vOmega)

    def estimate(self, Qa, pitch, omega, WS0, relaxation=0, WSavg=None, debug=False, method='min', deltaWSMax=1): 
        """
        INPUTS:
         - Qa: aerodynamic torque [Nm]
         - omega: rotational speed [rad/s]
         - pitch: pitch angle [deg]
         - WS0:  wind speed guess/previous estimate [m/s]
         - method: method
         # TODO compute rolling average on the fly

        NOTE: 
          - 'min'/ fCP : uses cubic interpolation
          - 'crossing': uses linear interpolation (but at higher res thanks)
        """
        info=None
        if debug:
            # Store state
            info={}
            info['Qa']=Qa; info['pitch']=pitch; info['omega']=omega; info['WS0']=WS0
            info['relaxation']=relaxation; info['method']=method; info['deltaWSMax']=deltaWSMax;


        def estim(WS0, delta, maxiter=50, tol=0.0001):
            vWS, vQ = self.TorqueAt(Pitch=pitch, Omega=omega)
            try:
                fQ = si.interp1d(vWS, vQ, kind='cubic')
            except:
                raise Exception()

            #fun = lambda WS : abs(Qa - fQ(WS))
            def fun(WS):
                res = abs(Qa - fQ(WS))
                #print('WS:{:.3f} res:{:15.5f} Qa:{:10.1f} {:10.1f}'.format(*(WS, res, Qa, fQ(WS))))
                return res
            WSmin = max(max(0, WS0-delta), vWS[0])
            WSmax = min(min(WS0+delta, self.WSmax), vWS[-1])
            #print('Bounds:', WSmin, WSmax, 'Mid:', (WSmax+WSmin)/2, vWS[0], vWS[-1])
            #fun = lambda WS : abs(Qa - Qaero(WS, pitch, omega, self.R, self.rho, self.fCP)) # OLD
            res = minimize_scalar(fun, bounds=[WSmin, WSmax], method='bounded', options={'xatol': tol, 'maxiter': maxiter})
            residual = Qa - Qaero(res.x, pitch, omega, self.R, self.rho, self.fCP)
            return res.x, residual

#         if WSavg is not None:
#             WS0=(WS0+WSavg)/2
        WS_est = WS0
        if method.find('min')>=0:
            if omega<=0.1:
                WS_est = WS0
            elif omega<self.OmegaLow:
                if self.OP is not None:
                    ws_guess=np.interp(omega, self.Omega, self.WS)
                else:
                    ws_guess=WS0 # TODO
                WS1,residual = estim(WS0, delta=2)
                WS_est=(4*WS1+ws_guess)/5
            else:
                WS_est, residual = estim(WS0, delta=deltaWSMax)
                #print('residual',residual)

        WSoper = None
        # --- Use operating conditions
        if method.find('oper')>=0:
            if self.OP is None:
                raise Exception('Cannot use method `oper`, operFile was not provided')
            WSOP    = self.OP['WS_[m/s]'].values   #[m/s]
            QaeroOP = self.OP['Qaero_[Nm]'].values #[Nm]
            PitchOP = self.OP['Pitch_[deg]'].values #[deg]
            # np.interp(Qa, self.OP['RtAeroMxh_[N-m]'], self.OP['WS_[m/s]'])
            WScrossOP, _, _ = zero_crossings(QaeroOP-Qa, x=WSOP)
            if len(WScrossOP)==0:
                # Can happen if torque below minimum or above maximum torque
                #print('>>> No crossing OP', pitch, omega, Qa)
                #import matplotlib.pyplot as plt
                #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                #ax.plot(WSOP, QaeroOP, label='Omega = {:3.1f}'.format(omega))
                #ax.plot(WSOP, Qa+WSOP*0)
                #ax.set_xlabel('WS [m/s]')
                #ax.set_ylabel('Q [N]')
                #ax.set_xlim([0,30])
                #ax.set_ylim([0,2e7])
                #ax.legend()
                #plt.show()
                WSoper = WS0
            elif len(WScrossOP)==1:
                WSoper = WScrossOP[0]
            elif len(WScrossOP)>=1:
                #print('{} OPcross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossOP), pitch, omega, Qa, WS0), WScrossOP)
                WSoper=WScrossOP[0]
                #import matplotlib.pyplot as plt
                #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                #ax.plot(WSOP, QaeroOP, label='Omega = {:3.1f}'.format(omega))
                #ax.plot(WSOP, Qa+WSOP*0)
                #ax.set_xlabel('WS [m/s]')
                #ax.set_ylabel('Q [N]')
                #ax.set_xlim([0,30])
                #ax.set_ylim([0,2e7])
                #ax.legend()
                #plt.show()
                #raise Exception('TODO')
            WSest = WSoper
            if debug:
                info['WS_oper'] = WSoper
                info['WS_crossOP'] = WScrossOP

        # --- Find torque for given pitch and omega
        if method.find('crossing')>=0:
            if omega>0:
                iNear = None
                vWS, vQ = self.TorqueAt(Pitch=pitch, Omega=omega)
                WScross, iBef, sign = zero_crossings(vQ-Qa, x=vWS)
                if len(WScross)==0:
                    #print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScross), pitch, omega, Qa, WS0), WScross)
                    #import matplotlib.pyplot as plt
                    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                    #ax.plot(vWS, vQ, label='Omega = {:3.1f}'.format(omega))
                    #ax.plot(vWS, Qa+vWS*0)
                    #ax.set_xlabel('WS [m/s]')
                    #ax.set_ylabel('Q [N]')
                    #ax.set_xlim([0,30])
                    #ax.set_ylim([0,2e7])
                    #ax.legend()
                    #plt.show()
                    WS_est = WS_est # WS0
                elif len(WScross)==1:
                    #print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScross), pitch, omega, Qa, WS0), WScross)
                    WS_est = WScross[0]
                elif len(WScross)==2:
                    #print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScross), pitch, omega, Qa, WS0), WScross)
                    if WSoper is not None:
                        iNear = np.argmin(abs(WScross-WSoper))
                    else:
                        iNear = np.argmin(abs(WScross-WS0))
                    WS_est = WScross[iNear] # TODO
                else:
                    #print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScross), pitch, omega, Qa, WS0), WScross)
                    # Find WS based on OP
                    if WSoper is not None:
                        iNear = np.argmin(abs(WScross-WSoper))
                    else:
                        iNear = np.argmin(abs(WScross-WS0))
                    WS_est = WScross[iNear] # TODO
        #             import matplotlib.pyplot as plt
        #             fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #             fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #             ax.plot(vWS, vQ, label='Omega = {:3.1f}'.format(omega))
        #             ax.plot(vWS, Qa+vWS*0)
        #             ax.set_xlabel('WS [m/s]')
        #             ax.set_ylabel('Q [N]')
        #             ax.set_xlim([0,30])
        #             ax.set_ylim([0,2e7])
        #             ax.legend()
        #             plt.show()

                # --- Use closest point
                if WSoper is not None:
                    WS_guess = (WSoper+WS0)/2
                else:
                    WS_guess = WS0
                if abs(WS_est - WS_guess) > deltaWSMax:
                    #print('>>>', WS_est)
                    dist = np.sqrt( ((vWS-WS_guess)/WS_guess)**2 + ((vQ-Qa)/Qa)**2 )
                    i = np.argmin(dist)
                    WS_est = vWS[i]

                if debug:
                    info['WS_cross'] = WScross
                    info['iNear'] = iNear
                    #info['cross_vQ'] = vQ
                    #info['cross_vWS'] = vWS

            



#         if omega<self.OmegaRated*0.95:
#             #  Below rated, we have the quasi steady ws as functin of omega as a best guess
#             ws_qs=np.interp(omega,self.Omega,self.WS)
#             if omega<self.OmegaLow:
#                 WS1,residual = estim(WS0,3)
#                 WS=(4*WS1+ws_qs)/5
#             else:
#                 WS,residual = estim(WS0,3)
# 
#             if np.abs(WS-ws_qs)>4:
#                 WS,residual = estim(ws_qs,3, maxiter=1000)
# 
#             if np.abs(residual)/Qa>0.1:
#                 WS,residual = estim(ws_qs, 17, maxiter=1000)
# 
#             if np.abs(residual)/Qa>0.1:
#                 WS,residual = estim(ws_qs+10, 10, maxiter=1000)
# 
#             if np.abs(residual)/Qa>0.1:
#                 if WSavg is not None:
#                     print('NOT GOOD 1 - WS={:.1f} WSqs={:.1f} WSavg={:.1f} - om={:.2f} pitch={:.2f}'.format(WS,ws_qs,WSavg,omega,pitch))
#                 else:
#                     print('NOT GOOD 1 - WS={:.1f} WSqs={:.1f} - om={:.2f} pitch={:.2f}'.format(WS,ws_qs,omega,pitch))
#         else:
#             # above omega rated, we are between WSrated-3 and WSCutoff
#             WS,residual = estim(WS0,3)
#             if WS<self.WSRated:
#                 WSmid=(self.WSCutOff+self.WSRated)/2
#                 WS,residual = estim(WSmid, 16, maxiter=1000)
# 
# #                 if WSavg is not None:
# #                     if np.abs(WS-WSavg)>4:
# #                         WS,residual = estim(WSavg, 6, maxiter=1000)
# 
#             if np.abs(residual)/Qa>0.1:
#                 print('NOT GOOD 2 - WS={:.1f} WS0={:.1f} - om={:.2f} pitch={:.2f}'.format(WS,WS0,omega,pitch))

        WS = WS0*relaxation + (1-relaxation)*WS_est

        if debug:
            info['WS_est1'] = WS_est
            info['WS_est2'] = WS
            self._debug_info=info

        return WS, info


    def estimateTimeSeries(self, Qaero, Pitch, Omega, WS_prev=None, WS_ref=None, debug=False, **kwargs):
        """ 
        Perform wind speed estimation given a time series of aerodynamic torque, pitch and rotational speed
        """
        print('Estimating WS on time series...')
        WS_est = np.zeros(Omega.shape)
        if WS_prev is None:
            WS_prev = 1
        ts_info = None
        if debug:
            # Storage for debug
            pass
        for i,(Qa, pitch, omega) in enumerate(zip(Qaero, Pitch, Omega)):
            ws_hat, info    = self.estimate(Qa, pitch, omega, WS_prev, debug=debug, **kwargs)
            WS_est[i] = ws_hat
            WS_prev   = ws_hat
            if debug:
                if WS_ref is not None:
                    if abs(WS_est[i] - WS_ref[i])>6:
                        print('[FAIL] Error too big:\n')
                        info['WS_ref'] = WS_ref[i]
                        self.debugPlot(info=info, HR=False)

                        return WS_est, ts_info
        return WS_est, ts_info

    def estimateTimeSeriesFromOF(self, outFilename, tRange=None, **kwargs):
        """" 
         - tRange: tuple (tmin, tmax) to limit the time used
        """
        import welib.weio as weio
        df = weio.read(outFilename).toDataFrame()
        if tRange is not None:
            df = df[np.logical_and(df['Time_[s]']>=tRange[0],df['Time_[s]']<=tRange[1])]
        time       = df['Time_[s]'].values
        WS_ref     = df['RtVAvgxh_[m/s]'].values # Rotor avg
        Pitch      = df['BldPitch1_[deg]'].values
        Qaero_ref  = df['RtFldMxh_[N-m]'].values
        Omega      = df['RotSpeed_[rpm]'].values*2*np.pi/60 # rad/s
        lambda_ref = Omega*self.R/WS_ref
        # Estimating wind speed on time series
        WS_est, ts_info = self.estimateTimeSeries(Qaero_ref, Pitch, Omega, WS_prev=WS_ref[0]*0.9, WS_ref=WS_ref, **kwargs)
        # Evaluating torque
        Qaero_eval = self.Torque(WS_ref, Pitch, Omega)
        Qaero_est  = self.Torque(WS_est, Pitch, Omega)
        # Storing data into a dataframe
        M    = np.column_stack((time, WS_ref, WS_est, Qaero_ref, Qaero_eval, Qaero_est, Omega, Pitch))
        cols = ['Time_[s]','WS_ref_[m/s]','WS_est_[m/s]','Qaero_ref_[N]','Qaero_eval_[N]','Qaero_est_[N]','Omega_[rad/s]','Pitch_[deg]']
        dfOut = pd.DataFrame(data=M, columns=cols)
        return dfOut

    def debugPlot(self, info=None, HR=False):
        from welib.tools.colors import python_colors


        if info is None:
            info = self._debug_info
        for k,v in info.items():
            print('{:15s}={}'.format(k,v))


        pitch = info['pitch']
        omega = info['omega']
        Qa    = info['Qa']
        WS_guess= info['WS0']
        try:
            WS_est= info['WS_est1']
            Qeval1 = self.Torque(info['WS_est1'], pitch, omega)
        except:
            WS_est = None
            Qeval1 = None
        try:
            WS_ref= info['WS_ref']
        except:
            WS_ref = None

        WS0, Q0 = self.TorqueAt(Pitch=pitch, Omega=omega)
        WS1, Q1 = self.TorqueAt(Pitch=pitch, Omega=omega*0.95)
        WS2, Q2 = self.TorqueAt(Pitch=pitch, Omega=omega*1.05)


        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(WS0   , Q0      , '-+', c=python_colors(0),       label='Interpolated data (at LambdaHR)')
        ax.plot(WS1   , Q1      , '--', c=python_colors(0), alpha=0.5)
        ax.plot(WS2   , Q2      , '--', c=python_colors(0), alpha=0.5)
        if WS_ref is not None:
            ax.plot(WS_ref, Qa      , 'ko',       label='Ref')
        if WS_est is not None:
            ax.plot(WS_est, Qeval1  , 'd', c=python_colors(1) , ms=5, label='Est.')
        ax.plot(WS0   , WS0*0+Qa, 'k--',      label='Qa')
        ax.plot([WS_guess, WS_guess], [np.min(Q0), np.max(Q0)], 'k:', label='WS guess/prev')

        if 'WS_cross' in info.keys():
            for iWS,WS in enumerate(info['WS_cross']):
                ax.plot(WS, Qa  , '*' , c=python_colors(3), ms=4, label='Cross (n={})'.format(len(info['WS_cross'])) if iWS==0 else None)

        if 'WS_crossOP' in info.keys():
            for iWS,WS in enumerate(info['WS_crossOP']):
                ax.plot(WS, Qa  , '^' , c=python_colors(4), ms=7, label='CrossOP (n={})'.format(len(info['WS_crossOP'])) if iWS==0 else None)
            ax.plot(info['WS_oper'], Qa, 'o',  label='Oper', ms=4)

            WSOP    = self.OP['WS_[m/s]'].values   #[m/s]
            QaeroOP = self.OP['Qaero_[Nm]'].values #[Nm]
            ax.plot(WSOP, QaeroOP, 'k--', label='OP')


        if HR:
            nHR=1000
            WS_HR = np.linspace(WS0[0],WS0[-1],nHR)
            Qeval_HR = self.Torque(WS_HR, [pitch]*nHR, [omega]*nHR)
            ax.plot(WS_HR, Qeval_HR, 'k-', lw=0.5, label='High Res eval')
        ax.set_xlabel('Wind Speed [m/s]')
        ax.set_ylabel('Torque [N]')
        ax.legend()

    def __repr__(self):
        s=''
        s+='<ws_estimator.TabulatedWSEstimator object> \n'
        s+=' - Lambda : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.Lambda),np.max(self.Lambda),self.Lambda[1]-self.Lambda[0], len(self.Lambda))
        s+=' - Pitch  : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.Pitch) ,np.max(self.Pitch) ,self.Pitch[1]-self.Pitch[0]  , len(self.Pitch))
        s+=' - CP     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CP),np.max(self.CP),self.CP.shape[0],self.CP.shape[1])
        if self.CT is not None:
            s+=' - CT     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CT),np.max(self.CT),self.CT.shape[0],self.CT.shape[1])
        s+=' - R      : {}  \n'.format(self.R)
        s+=' - rho    : {}  \n'.format(self.rho)
        s+=' - OmegaLow:{}  \n'.format(self.OmegaLow)
        # files
        s+=' - fstFile:     {}  \n'.format(self.fstFile)
        s+=' - aeroMapFile: {}  \n'.format(self.aeroMapFile)
        s+=' - operFile:    {}  \n'.format(self.operFile)
        return s


if __name__=='__main__':
    import pandas as pd
    import matplotlib.pyplot as plt
    from spectral import fft_wrap
    import welib.weio as weio
    # --- Parameters
    # InputFile = 'GeneratorDynamics.outb'
    InputFile = 'DLC120_ws13_ye000_s1_r1.outb'

    # --- Turbine data
    turbine = dict()
    turbine['R']         = 63
    turbine['rho']   = 1.225
    g=9.81

    # --- Reading aerodynamic data for the turbine
    Pitch  = pd.read_csv('Pitch_data.csv',header  = -1).values
    Lambda = pd.read_csv('Lambda_data.csv',header = -1).values
    CP     = pd.read_csv('CP_data.csv',header     = -1).values
    CT     = pd.read_csv('CT_data.csv',header     = -1).values
    # Create the interpolant for CP and CT, CP(pitch,lambda) (same interface interp2d) 
    turbine['fCP'] = interp2d_pairs(Pitch,Lambda,CP,kind='cubic')
    turbine['fCT'] = interp2d_pairs(Pitch,Lambda,CT,kind='cubic')



    # --- Reading in some "measuremensts" or "simulation"
    # --- Removing units from columns
    df=weio.read(InputFile).toDataFrame()
    df.columns = [  v.split('_[')[0] for v in df.columns.values] 
    time      = df['Time'].values
    genspeed  = df['GenSpeed'].values * 2*np.pi/60 # converted to rad/s
    rotspeed  = df['RotSpeed'].values * 2*np.pi/60 # converted to rad/s
    thrust    = df['RotThrust']*1000
    gentq     = df['GenTq']*1000*97                # Convert to rot torque
    azimuth   = df['Azimuth']                      # deg
    windspeed = df['Wind1VelX']
    pitch     = df['BldPitch3']
    rottorq   = df['RotTorq']*1000
    rottorq2  = df['RtFldMxh']
    thrust2   = df['RtFldFxh']




    # --- Evaluate the interpolant on each pairs of x and y values
    F = Taero(windspeed, pitch, rotspeed, turbine['R'], turbine['rho'], turbine['fCT'])
    
    Q = Qaero(windspeed, pitch, rotspeed, turbine['R'], turbine['rho'], turbine['fCP'])
    
    # --- normalize F, Q, rottorq2 and thrust2 to imrpove spectra analysis
    Fnorm = F - np.average(F)
    Qnorm = Q - np.average(Q)
    thrust2norm = thrust2 - np.average(thrust2)
    rottorq2norm = rottorq2 - np.average(rottorq2)
    # --- Import functio for calculating spectra from Q and F data
#    f1, S1, Info  = fft_wrap(time,Fnorm,output_type='amplitude',averaging='Welch')
#    f2, S2, Info  = fft_wrap(time,rottorq2norm,output_type='amplitude',averaging='Welch')
    f1, S1, Info  = fft_wrap(time,Fnorm,output_type='PSD',averaging='Welch')
    f2, S2, Info  = fft_wrap(time,thrust2norm,output_type='PSD',averaging='Welch')
    
    
    
    # --- Figures
    
    # -- Figure 1: Qaero
    plt.figure(num=1, figsize=(11, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(time,Q/1000,'r--',label='Qaero calc')
    plt.plot(time,rottorq2/1000,'k',label = 'Qaero fast')
    plt.gca().legend()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('time (s)')
    plt.ylabel('Qaeo (kN.m)')
    
    # --- Figure 2: Faero
    plt.figure(num=2, figsize=(11, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(time,F      ,'c--',label='Faero calc')
    plt.plot(time,thrust2,'m'  ,label='Faero (FAST)')
    plt.plot(time,thrust ,'b:' ,label='Faero (FAST Fstruct - Fweight)')
    plt.gca().legend()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('time (s)')
    plt.ylabel('Faero (kN)')

    # --- Figure 3: Spectra between Faero calc and Faero FAST
    plt.figure(num=3, figsize=(11, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(f1,S1,'r:',label='Faero calc')
    plt.plot(f2,S2,'k--',label='Faero (FAST)')
    plt.gca().legend()
    plt.axvline(x=.201)
    plt.axvline(x=3*.201)
    plt.axvline(x=9*.201)
    
    plt.xlim([-.0001, 5])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density (Welch Avg.)') 
    plt.yscale('log')
