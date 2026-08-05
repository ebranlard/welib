import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.optimize import minimize_scalar

from welib.weio.fast_input_deck import FASTInputDeck
from welib.tools.signal_analysis import zero_crossings
from welib.tools.strings import WARN, OK, INFO, FAIL
from welib.tools.clean_exceptions import *
from welib.tools.colors import python_colors, fColrs, lighten_color
from welib.tools.stats import rsquare, mean_rel_err, comparison_stats

# ---
def interp2d_pairs(X, Y, Z, kind='cubic', **kwargs):
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
    # --- OLD
    #     def interpolant(x,y,f):
    #         x,y = np.asarray(x), np.asarray(y)
    #         return (si.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x.ravel(), y.ravel())[0]).reshape(x.shape)
    #return lambda x,y: interpolant(x, y, si.interp2d(X, Y, Z, kind=kind))
    # --- NEW
    xmin, xmax = np.min(X.flatten()), np.max(X.flatten())
    ymin, ymax = np.min(Y.flatten()), np.max(Y.flatten())
    Finterp = si.RegularGridInterpolator((X,Y), Z.T, method=kind)
    #r = si.RectBivariateSpline(X, Y, Z.T)
    def interpolant(x,y):
        x,y = np.asarray(x), np.asarray(y)
        x = np.clip(x, xmin, xmax)
        y = np.clip(y, ymin, ymax)
        return Finterp((x,y))
    return interpolant


def Paero(WS, pitch, omega, R, rho, fCP):
    """ Taero returns the aerodynamic power
         - pitch  [deg]
         - omega [rad/s]
         - R : the blade radius [m]
         - fCP : an interpolant for CP(pitch,lambda) as returned by interp2d_paris
         - rho : the air density [kg/m^3]
    """
    Lambda = omega * R / WS
    CP     = fCP(pitch, Lambda)
    P      = 1/2*rho*np.pi*R**2*WS**3*CP
    return P


def Qaero(WS, pitch, omega, R, rho, fCP):
    """ Qaero returns the aerodynamic torque
         - pitch [deg]
         - omega [rad/s]
         - R : the blade radius
         - fCP : an interpolant for CP(pitch,lambda)
         - rho : the air density
    """
    pitch = np.asarray(pitch)
    WS    = np.asarray(WS)
    omega = np.asarray(omega)
    Lambda = omega * R / WS
    CP = fCP(pitch, Lambda)
    Q = 1/2*rho*np.pi*R**2*WS**3/omega*CP
    return Q

def Taero(WS, pitch, omega, R, rho, fCT):
    """ Taero returns the aerodynamic thrust of a given turbine
         - pitch [deg]
         - omega [rad/s]
         - R : the blade radius
         - fCP : an interpolant for CP(pitch,lambda)
         - rho : the air density
    """
    Lambda = omega * R / WS
    CT = fCT(pitch,Lambda)
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
        self.CQ    = None
        self.OP    = None

        self.fCP    = None
        self.fCT    = None
        self.fCQ    = None

        # ---
        if fstFile:
            fst = FASTInputDeck(fstFile, readlist=['Fst', 'ED', 'AD'])
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
            outAD = [c.strip('"').strip().lower() for c in fst.AD['OutList']]
            outAD = [c for c in outAD if len(c)>0]
            if 'rtaeromxh' not in outAD:
                WARN('RtAeroMxh not found in AD out list')

        self.fstFile = fstFile
        self.R       = R
        self.rho     = rho

    def _sanitizeOP(self, OP, expectedCols=None, onlyExpected=True):
        if expectedCols is None:
            expectedCols = ['WS_[m/s]', 'Pitch_[deg]', 'RotSpeed_[rpm]']

        # --- Trying to be nice about column names
        OP.columns = [c.lower().replace(' ','_').replace('(','[').replace(')',']') for c in OP.columns]

        d =dict([(k, 'WS_[m/s]') for k in ['ws_[m/s]', 'ws', 'wind_[m/s]', 'wind_speed_[m/s]']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'Pitch_[deg]') for k in ['bldpitch1_[deg]', 'bldpitch_[deg]','pitch_[deg]', 'pitch']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'RotSpeed_[rpm]') for k in ['rotspeed_[rpm]', 'rpm', 'rpm_[rpm]', 'omega_[rpm]', 'rotor_speed_[rpm]']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'AeroPower_[kW]') for k in ['mech_power_[kw]', 'rtaeropwr_[kw]']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'TSR_[-]') for k in ['tsr', 'rtaerotsr_[-]']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'CP_[-]') for k in ['cp', 'cp_[-]']])
        OP.rename(columns = d, inplace=True)

        # TODO standardize Units WE
        if 'rtaeromxh_[kn-m]' in OP.keys(): 
            OP['rtaeromxh_[n-m]'] =  OP['rtaeromxh_[kn-m]'].values*1000
        if 'torque_[knm]' in OP.keys(): 
            OP['torque_[nm]'] =  OP['torque_[knm]'].values*1000
        d =dict([(k, 'Qaero_[Nm]') for k in ['rtaeromxh_[n-m]', 'torque_[nm]']])
        OP.rename(columns = d, inplace=True)

        d =dict([(k, 'PhiY_[deg]') for k in ['phiy_[deg]']])
        OP.rename(columns = d, inplace=True)

        missing = [c for c in expectedCols if c not in OP]
        if len(missing)>0:
            print('>>> Columns', OP.keys())
            raise Exception('Missing columns : {}'.format(missing))


        if 'CP_[-]' not in OP and 'AeroPower_[kW]' in OP:
            OP['CP_[-]'] = OP['AeroPower_[kW]' ]*1000 / (1/2* self.rho * OP['WS_[m/s]']**3 * self.R**2 * np.pi )

        if onlyExpected:
            OP = OP[expectedCols]

        return OP


class TabulatedWSEstimator(TabulatedWSEstimatorBase):

    def __init__(self, R=None, rho=1.225, fstFile=None, basename=None, operFile=None, aeroMapFile=None, omegaLow=0, omegaRated=10):
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
        self.pitch    = None
        self.Lambda   = None
        # Operating condition
        self.omegaLow   = omegaLow
        self.omegaRated = omegaRated
        self.OP     = None
        # Files
        self.operFile = operFile
        self.aeroMapFile = aeroMapFile

        if aeroMapFile is not None:
            self.loadAeroMap(aeroMapFile)
        if operFile is not None:
            self.loadOper(operFile)


    def loadOper(self, operFile):

        if not os.path.exists(operFile):
            print('[WARN] Operating point file not found: ',operFile)
            operFile+'    [NOT FOUND]'
        else:
            #print('>>> Loading oper file: ',operFile)
            import welib.weio as weio
            OP = weio.read(operFile).toDataFrame()
            self._setOP(OP)
        self.operFile = operFile
        self._interpOP() # (Needs weights to be computed first)

    def setDB(self, WS, pitch, rpm, phiy, CP, CT):
        print('TODO taken from Floating, need care')
        self.WS     = np.asarray(WS)
        self.pitch  = np.asarray(pitch)
        self.omega  = np.asarray(rpm)*np.pi/30
        self.CP     = CP
        self.CT     = CT
        MWS     = np.full((len(WS), len(rpm), len(pitch)), np.nan)
        Momega  = np.full((len(WS), len(rpm), len(pitch)), np.nan)
        assert(MWS.shape == CP.shape)
        assert(CP.shape == CT.shape)
        # TODO vectorize
        for i,ws in enumerate(WS): 
            for j,om in enumerate(self.omega): 
                for k,pit in enumerate(pitch): 
                    MWS   [i,j,k,l] = ws
                    Momega[i,j,k,l] = om
        #         MWS     = MWS             # TODO generated from WS...
        #         Momega  = MRPM*np.pi/30
        self.CP[np.isnan(self.CP)]=0
        self.CP[self.CP<0]=0
        self.CT[np.isnan(self.CT)]=0
        self.CT[self.CT<0]=0

        # --- Computing weights
        if self.R is None:
            raise Exception('R should be set')
        if self.rho is None:
            raise Exception('rho should be set')
        P = self.CP * 1/2 * self.rho * np.pi * self.R**2 * MWS**3
        self.P = P
        if self.CT is not None:
            T = self.CT * 1/2 * self.rho * np.pi * self.R**2 * MWS**2
        Q = P/Momega
        self.Q = Q
        self.computeWeights(P, Q, T) # TODO

        if self.OP is not None:
            # Trigger
            self._interpOP() # (Needs weights to be computed first)


    def setFromTimeSeries(self, df, nWS=6, nRPM=6, nPitch=5, nPhi=4):
        """ """
        print('TODO taken from Floating, need care')
        # --- Time series
        Q     = df['Qaero'].values
        WS    = df['WS'].values
        omega = df['dpsi'].values          # rad/s
        rpm   = omega * 30/ np.pi
        pitch = df['pitch'].values*180/np.pi # deg
        P     = df['power'].values
        T     = df['Thrust'].values # ...Hoping that's correct
        R   = self.R
        rho = self.rho
        CP = P / (1/2 * rho * np.pi * R**2 * WS**3)
        CT = T / (1/2 * rho * np.pi * R**2 * WS**2)

        # ---- Bins
        WSb    = np.linspace(np.min(WS)  , np.max(WS), nWS)
        RPMb   = np.linspace(np.min(rpm) , np.max(rpm), nRPM)
        Pitchb = np.linspace(np.min(pitch), np.max(pitch), nPitch)
        dws = np.diff(WSb)[0]
        dom = np.diff(RPMb)[0]
        dpi = np.diff(Pitchb)[0]
        dph = np.diff(Phiyb)[0]

        # --- Filling up CP/CT
        MCP = np.full((len(WSb), len(RPMb), len(Pitchb)), np.nan)
        MCT = np.full((len(WSb), len(RPMb), len(Pitchb)), np.nan)
        for i,ws in enumerate(WSb): 
            for j,om in enumerate(RPMb): 
                for k,pit in enumerate(Pitchb): 
                    bWS = np.logical_and(WS    >= ws-dws , WS    <= ws+dws)
                    bOM = np.logical_and(rpm   >= om-dom , rpm   <= om+dom)
                    bPI = np.logical_and(pitch >= pit-dpi, pitch <= pit+dpi)
                    bPH = np.logical_and(phiy  >= ph-dph , phiy  <= ph+dph)
                    if len(bWS)==0:
                        print('>>> NO WS SELECTION')
                    if len(bOM)==0:
                        print('>>> NO OMEGA SELECTION')
                    if len(bPI)==0:
                        print('>>> NO PITCH SELECTION')
                    bAll = np.logical_and.reduce((bWS,bOM,bPI,bPH))
                    IAll = np.where(bAll)[0]
                    #if len(IAll)==0:
                    #    print('>>> NO Total SELECTION')
                    MCP[i,j,k] = np.mean(CP[bAll])
                    MCT[i,j,k] = np.mean(CT[bAll])
        # --- Reset DB
        self.setDB(WSb, Pitchb, RPMb, MCP, MCT)

    def _setOP(self, OP):
        # --- Operating conditions
        OP = self._sanitizeOP(OP, expectedCols=['WS_[m/s]', 'Pitch_[deg]', 'RotSpeed_[rpm]', 'Qaero_[Nm]'], onlyExpected=False)
 
        self.OP=OP
        self.WS_op      = OP['WS_[m/s]'].values
        self.omega_op   = OP['RotSpeed_[rpm]'].values*2*np.pi/60 # [rad/s]
        self.omegaRated = np.max(self.omega_op)
        self.omegaLow   = 0.4 * self.omegaRated                                      # TODO
        self.WSRated    = np.interp(self.omegaRated*0.98, self.omega_op, self.WS_op) # TODO
        self.WSCutOff   = 28


    def _interpOP(self):
        # Needs Weights to be computed
        # ---  Compute interpolated values at Operating points to be consistent
        WS    = self.OP['WS_[m/s]'].values[:]
        omega = self.OP['RotSpeed_[rpm]'].values[:] *np.pi/30
        pitch = self.OP['Pitch_[deg]'].values[:]
        self.OP['Paero_i_[W]'] = self.Power(WS, omega, pitch) # Interpolated power
        self.OP['Qaero_i_[Nm]'] = self.OP['Paero_i_[W]']/omega
#         Q2                     = self.Torque(WS, omega, pitch, phiy) # Interpolated torque
#         Q1                      = self.OP['Paero_i_[W]']/omega
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(WS, Q1 ,'-'   , label='')
#         ax.plot(WS, Q2 ,'--'   , label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         plt.show()

    def loadAeroMap(self, aeroMapFile):
        """ Load file containing aeromap: CP(Lambda, Pitch) """
        #print('>>>> Loading aeroMapFile', aeroMapFile)
        # TODO more file formats
        from welib.weio.rosco_performance_file import ROSCOPerformanceFile
        rs = ROSCOPerformanceFile(aeroMapFile)
        self.pitch  = rs['pitch']
        self.Lambda = rs['TSR']
        self.CP     = rs['CP']
        self.CT     = rs['CT']
        self.CQ     = rs['CQ']
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
            self.pitch  = pd.read_csv(PitchFile ,header = None).values.ravel()
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
        self.fCP = interp2d_pairs(self.pitch, self.Lambda, self.CP, kind='cubic')

        if self.CQ is not None:
            self.CQ[self.CQ<=0]=0
            self.fCQ = interp2d_pairs(self.pitch, self.Lambda, self.CQ, kind='cubic')

        if self.CT is not None:
            self.CT[self.CT<=0]=0
            self.fCT = interp2d_pairs(self.pitch, self.Lambda, self.CT, kind='cubic')


        # --- Lambda at higher res
        LambdaMid = self.Lambda[:-1] + np.diff(self.Lambda)/2
        Lambda = np.sort(np.concatenate((self.Lambda, LambdaMid)))
        self.LambdaHR = Lambda


    def CP_eval(self, WS, pitch, omega):
        #P = Paero(WS, pitch, omega, self.R, self.rho, self.fCP)
        Lambda = omega * self.R / WS
        CP     = self.fCP(pitch, Lambda)
        return CP

    def Power(self, WS, pitch, omega):
        """
        Return power from fCP
         - WS: wind speed [m/s]
         - omega: rotational speed [rad/s]
         - pitch: pitch angle [deg]
         - phiy: platform pitch angle [deg]
         """
        return Paero(WS, pitch, omega, self.R, self.rho, self.fCP)

    def Thrust(self, WS, pitch, omega):
        """
        Return Thrust from fCP
         - WS: wind speed [m/s]
         - omega: rotational speed [rad/s]
         - pitch: pitch angle [deg]
         - phiy: platform pitch angle [deg]
         """
        return Taero(WS, pitch, omega, self.R, self.rho, self.fCT)

    def Torque(self, WS, pitch, omega):
        """
        Return Torque from fCP
         - WS: wind speed [m/s]
         - omega: rotational speed [rad/s]
         - pitch: pitch angle [deg]
         - phiy: platform pitch angle [deg]
        """
        return Qaero(WS, pitch, omega, self.R, self.rho, self.fCP)

    def TorqueFromCQ(self, WS, pitch, omega):
        pitch = np.asarray(pitch)
        WS    = np.asarray(WS)
        omega = np.asarray(omega)
        Lambda = omega * self.R / WS
        CQ = self.fCQ(pitch, Lambda)
        Q = 1/2 * self.rho * np.pi * self.R**3 * WS**2 * CQ
        return Q

    def TorqueAt(self, pitch, omega):
        """ 
        Return Torque(WS) curve for a given pitch and rotational speed
        Pitch,omega: scalar
        """
        WS     = omega * self.R / self.LambdaHR[-1::-1] # NOTE using LambdaHR to benefit from cubic interpolation
        WS = WS[WS<self.WSmax]
        vPitch = np.array([pitch]*len(WS))
        vomega = np.array([omega]*len(WS))
        return WS, self.Torque(WS, vPitch, vomega)

    def estimate(self, Qa, pitch, omega, WS0, relaxation=0, WSavg=None, method='min', deltaWSMax=1, verbose=False, debug=False, t=0): 
        """
        INPUTS:
         - Qa: aerodynamic torque [Nm]
         - omega: rotational speed [rad/s]
         - pitch: pitch angle [deg]
         - WS0:  wind speed guess/previous estimate [m/s]
         - method: method in 
              'min'     : use minimize_scalar optimization
              'oper'    : use operating conditions only
              'crossing': use crossings with Cp curve
         # TODO compute rolling average on the fly

        NOTE: 
          - 'min'/ fCP : uses cubic interpolation
          - 'crossing': uses linear interpolation (but at higher res thanks)
        """
        info=None

        def saveState(info=None):
            if info is None:
                info={}
            # Store state
            info['Qa']=Qa; info['pitch']=pitch; info['omega']=omega; info['WS0']=WS0
            info['relaxation']=relaxation; info['method']=method; info['deltaWSMax']=deltaWSMax;
            info['t']=t
            return info
            
        if debug:
            info = saveState()


        def estim(WS0, delta, maxiter=50, tol=0.0001):
            vWS, vQ = self.TorqueAt(pitch=pitch, omega=omega)
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
            elif omega<self.omegaLow:
                if self.OP is not None:
                    ws_guess=np.interp(omega, self.omega_op, self.WS_op)
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
                if verbose:
                    print('>>> No crossing OP', pitch, omega, Qa)
                    print('{} OPcross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossOP), pitch, omega, Qa, WS0), WScrossOP)
                #import matplotlib.pyplot as plt
                #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                #ax.plot(WSOP, QaeroOP, label='omega = {:3.1f}'.format(omega))
                #ax.plot(WSOP, Qa+WSOP*0)
                #ax.set_xlabel('WS [m/s]')
                #ax.set_ylabel('Q [N]')
                #ax.set_xlim([0,30])
                #ax.set_ylim([0,2e7])
                #ax.legend()
                #plt.show()
                WSoper = WS0
            elif len(WScrossOP)==1:
                if verbose:
                    print('{} OPcross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossOP), pitch, omega, Qa, WS0), WScrossOP)
                WSoper = WScrossOP[0]
            elif len(WScrossOP)>=1:
                if verbose:
                    print('{} OPcross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossOP), pitch, omega, Qa, WS0), WScrossOP)
                WSoper=WScrossOP[0]
                #import matplotlib.pyplot as plt
                #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                #ax.plot(WSOP, QaeroOP, label='omega = {:3.1f}'.format(omega))
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
                vWS, vQ = self.TorqueAt(pitch=pitch, omega=omega)
                WScross, iBef, sign = zero_crossings(vQ-Qa, x=vWS)
                if len(WScross)==0:
                    #print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScross), pitch, omega, Qa, WS0), WScross)
                    #import matplotlib.pyplot as plt
                    #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                    #ax.plot(vWS, vQ, label='omega = {:3.1f}'.format(omega))
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
        #             ax.plot(vWS, vQ, label='omega = {:3.1f}'.format(omega))
        #             ax.plot(vWS, Qa+vWS*0)
        #             ax.set_xlabel('WS [m/s]')
        #             ax.set_ylabel('Q [N]')
        #             ax.set_xlim([0,30])
        #             ax.set_ylim([0,2e7])
        #             ax.legend()
        #             plt.show()

                # --- If estimate is far away, change it
                # --- Use closest point
                if WSoper is not None:
                    WS_guess = (WSoper+WS0)/2
                else:
                    WS_guess = WS0
                if abs(WS_est - WS_guess) > deltaWSMax:
                    if verbose:
                        print('Icross, jump too big')
                    dist = np.sqrt( ((vWS-WS_guess)/WS_guess)**2 + ((vQ-Qa)/Qa)**2 )
                    i = np.nanargmin(dist)
                    WS_est = vWS[i]

                if debug:
                    info['WS_cross'] = WScross
                    info['iNear'] = iNear
                    #info['cross_vQ'] = vQ
                    #info['cross_vWS'] = vWS

            

        WS = WS0*relaxation + (1-relaxation)*WS_est

        if debug:
            info['WS_est1'] = WS_est
            info['WS_est2'] = WS
            self._debug_info=info

        return WS, info


    def estimateTimeSeries(self, Qaero, pitch, omega, WS_prev=None, WS_ref=None, debug=False, **kwargs):
        """ 
        Perform wind speed estimation given a time series of aerodynamic torque, pitch and rotational speed
        """
        WS_est = np.zeros(omega.shape)
        if WS_prev is None:
            WS_prev = 1
        ts_info = None
        if debug:
            # Storage for debug
            pass
        for i,(Qa, pitch, omega) in enumerate(zip(Qaero, pitch, omega)):
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
        pitch      = df['BldPitch1_[deg]'].values
        try:
            Qaero_ref  = df['RtAeroMxh_[N-m]'].values
        except:
            Qaero_ref  = df['RtFldMxh_[N-m]'].values
        omega      = df['RotSpeed_[rpm]'].values*2*np.pi/60 # rad/s
        lambda_ref = omega*self.R/WS_ref
        # Estimating wind speed on time series
        WS_est, ts_info = self.estimateTimeSeries(Qaero_ref, pitch, omega, WS_prev=WS_ref[0]*0.9, WS_ref=WS_ref, **kwargs)
        # Evaluating torque
        Qaero_eval = self.Torque(WS_ref, pitch, omega)
        Qaero_est  = self.Torque(WS_est, pitch, omega)
        # Storing data into a dataframe
        M    = np.column_stack((time, WS_ref, WS_est, Qaero_ref, Qaero_eval, Qaero_est, omega, pitch))
        cols = ['Time_[s]','WS_ref_[m/s]','WS_est_[m/s]','Qaero_ref_[N]','Qaero_eval_[N]','Qaero_est_[N]','omega_[rad/s]','Pitch_[deg]']
        dfOut = pd.DataFrame(data=M, columns=cols)

        self.df = dfOut
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

        WS0, Q0 = self.TorqueAt(pitch=pitch, omega=omega)
        WS1, Q1 = self.TorqueAt(pitch=pitch, omega=omega*0.95)
        WS2, Q2 = self.TorqueAt(pitch=pitch, omega=omega*1.05)


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

    def operPlot(self):
        if self.OP is None:
            raise Exception()


    
        WS     = self.OP['WS_[m/s]'].values
        omega  = self.OP['RotSpeed_[rpm]'].values * 2*np.pi/60
        pitch  = self.OP['Pitch_[deg]'].values
        Qa_ref = self.OP['Qaero_[Nm]'].values



        # --- Method 1 just call raw method fCP
        Lambda = omega * self.R / WS
        CP = self.fCP(pitch, Lambda)
        Qa = self.Torque(WS, pitch, omega)
        Pa = self.Power (WS, pitch, omega)
        Qa2 = self.TorqueFromCQ(WS, pitch, omega)

        print(self.OP.keys())

        # --- Torque
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(WS, Qa_ref, 'k-', label='Reference from Oper')
        ax.plot(WS, Qa , '--', label='From CP')
        ax.plot(WS, Qa2, ':', label='From CQ')
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('Torque [Nm]')
        ax.legend()
# 
#         if 'TSR_[-]' in self.OP:
#             fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
#             fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#             ax.plot(WS, self.OP['TSR_[-]'], 'k-', label='Reference from Oper')
#             ax.plot(WS, Lambda, '--')
#             ax.set_xlabel('Wind speed [m/s]')
#             ax.set_ylabel('TSR [-]')
#             ax.legend()
# 
#         if 'AeroPower_[kW]' in self.OP:
#             fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
#             fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#             ax.plot(WS, self.OP['AeroPower_[kW]']*1000, 'k-', label='Reference from Oper')
#             ax.plot(WS, Pa, '--')
#             ax.set_xlabel('Wind speed [m/s]')
#             ax.set_ylabel('Power [kW]')
#             ax.legend()

        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        if 'CP_[-]' in self.OP:
            ax.plot(WS, self.OP['CP_[-]'], 'k-')
        ax.plot(WS, CP, '--')
#         ax.plot(WS, Qa, '--', label='From CPCTCQ')
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('CP [-]')
        ax.legend()


        return fig


    def plotTimeSeriesEstimation(self, axes=None):
        df = self.df
        # --- Plot
        fig,axes = plt.subplots(3, 1, sharex=False, figsize=(13.4,7.0)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)


        # TODO get those using  max of WS and Q and min max of oper
#         Ylim1 = [0,20]    # WS
#         Ylim2 = [0,4.5e6] # Q
#         Ylim3 = [-11,25] # Oper


        # --- Where the data is invalid
        #bInv = np.logical_or.reduce((df['omega_[rad/s]'] <self.omega[0],   df['omega_[rad/s]']  >self.omega[-1]))
        #bInv = np.logical_or.reduce((df['Pitch_[deg]']   <self.pitch[0],   df['Pitch_[deg]']    >self.pitch[-1], bInv))
        #bInv = np.logical_or.reduce((df['WS_ref_[m/s]']   <self.WS[0],     df['WS_ref_[m/s]']   >self.WS[-1]   , bInv))
        #bInv2 = np.isnan(df['Qaero_est_[N]'])
        ## #with Timer('ValidValues'):
        ## #    bInv2 =  ~self.validValues(df['WS_ref_[m/s]'], df['omega_[rad/s]'], df['Pitch_[deg]']  , df['PtfmPitch_[deg]'])
        ## #import pdb; pdb.set_trace()
        #bInv3 = np.logical_or(bInv, bInv2)
        #bInv = bInv3
        #b    = ~bInv3
        # 
        # --- WS plot
        ax=axes[0]

        # #ax.fill_between(t, Ylim1[0], Ylim1[1], where=bInv2, alpha=0.1, color=python_colors(1))
        # # ax.fill_between(t, Ylim1[0], Ylim1[1], where=bInv, alpha=0.1, color=(0.5,0.5,0.5))
        # 
        # # stats, sStats = comparison_stats(t[b], df['WS_ref_[m/s]'].values[b], t[b], df['WS_est_[m/s]'].values[b])
        # # ax.text(2,Ylim1[0]+(Ylim1[1]-Ylim1[0])*0.89, sStats, fontsize=11 )
        # 
        # # ax.axhline(y = self.WS[0 ], color=python_colors(0), linestyle = '--', lw=0.5)
        # # ax.axhline(y = self.WS[-1], color=python_colors(0), linestyle = '--', lw=0.5)
        # # ax.set_ylim(Ylim1)


        ax.plot(df['Time_[s]'], df['WS_ref_[m/s]'],     color=fColrs(1), label='OpenFAST')
        ax.plot(df['Time_[s]'], df['WS_est_[m/s]'], ':',color=fColrs(4), label='Estimated')
        ax.set_ylabel('Wind speed [m/s]')

        #  --- Qplot
        ax=axes[1]
        # #ax.fill_between(t, Ylim2[0], Ylim2[1], where=bInv2, alpha=0.1, color=python_colors(1))
        # # ax.fill_between(t, Ylim2[0], Ylim2[1], where=bInv, alpha=0.1, color=(0.5,0.5,0.5))
        # 
        # stats, sStats = comparison_stats(t[b], df['Qaero_ref_[N]'].values[b], t[b], df['Qaero_est_[N]'].values[b])
        # ax.text(2,Ylim2[0]+(Ylim2[1]-Ylim2[0])*0.89, sStats, fontsize=11 )
        # 
        ax.plot(df['Time_[s]'], df['Qaero_ref_[N]' ]    , color=fColrs(1),   label='OpenFAST')
        ax.plot(df['Time_[s]'], df['Qaero_est_[N]' ],':', color=fColrs(4),  label='From WS Estimated')
        # #ax.plot(df['Time_[s]'], df['Qaero_eval_[N]'], '--', label='Evaluated')
        ax.set_ylabel('Qaero [N]')
        # ax.set_ylim(Ylim2)
        ax.legend(loc='center left')
        # 
        # --- Oper
        ax=axes[2]
        # #ax.fill_between(t, Ylim3[0], Ylim3[1], where=bInv2, alpha=0.1, color=python_colors(1))
        # # ax.fill_between(t, Ylim3[0], Ylim3[1], where=bInv, alpha=0.1, color=(0.5,0.5,0.5))
        # 
        colrs=[python_colors(0), python_colors(1), python_colors(2)]
        colrs=[fColrs(1), lighten_color(fColrs(1),0.3), lighten_color(fColrs(1),0.6)]

        ax.plot(df['Time_[s]'], df['omega_[rad/s]']*30/np.pi, '-',c=colrs[0], label='omega [rpm]')
        ax.plot(df['Time_[s]'], df['Pitch_[deg]']           , '--',c=colrs[1], label='Pitch [deg]')
        # ax.plot(df['Time_[s]'], df['PtfmPitch_[deg]']       , '-.',c=colrs[2], label='PhiY [deg]')
        # # 
        # #             ax.axhline(y = wse.omega[0 ]*30/np.pi, color=colrs[0], linestyle = '--', lw=0.5)
        # #             ax.axhline(y = wse.omega[-1]*30/np.pi, color=colrs[0], linestyle = '--', lw=0.5)
        # #             ax.axhline(y = wse.pitch[0 ],          color=colrs[1], linestyle = '--', lw=0.5)
        # #             ax.axhline(y = wse.pitch[-1],          color=colrs[1], linestyle = '--', lw=0.5)
        # 
        # # ax.set_ylim(Ylim3)
        ax.set_xlabel('Time [s]')
        ax.legend()

        for ax in axes.flatten():
            #ax.set_xlim([0,600])
            ax.tick_params(direction='in')





    def __repr__(self):
        s=''
        s+='<ws_estimator.TabulatedWSEstimator object> \n'
        if self.Lambda is not None:
            s+=' - Lambda : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.Lambda),np.max(self.Lambda),self.Lambda[1]-self.Lambda[0], len(self.Lambda))
        else:
            s+=' - Lambda : {}\n'.format(self.Lambda)
        if self.pitch is not None:
            s+=' - pitch  : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.pitch) ,np.max(self.pitch) ,self.pitch[1]-self.pitch[0]  , len(self.pitch))
        else:
            s+=' - pitch  : {}\n'.format(self.pitch)
        if self.CP is not None:
            s+=' - CP     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CP),np.max(self.CP),self.CP.shape[0],self.CP.shape[1])
        else:
            s+=' - CP     : {} \n'.format(self.CP)
        if self.CT is not None:
            s+=' - CT     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CT),np.max(self.CT),self.CT.shape[0],self.CT.shape[1])
        s+=' - R      : {}  \n'.format(self.R)
        s+=' - rho    : {}  \n'.format(self.rho)
        s+=' - omegaLow:{}  \n'.format(self.omegaLow)
        # files
        s+=' - fstFile:     {}  \n'.format(self.fstFile)
        s+=' - aeroMapFile: {}  \n'.format(self.aeroMapFile)
        s+=' - operFile:    {}  \n'.format(self.operFile)
        s+=' - OP      :\n  {}  \n'.format(self.OP)
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
