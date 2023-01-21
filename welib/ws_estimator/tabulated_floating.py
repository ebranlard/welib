import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.optimize import minimize_scalar


from welib.ws_estimator.tabulated import TabulatedWSEstimatorBase

from welib.weio.fast_input_deck import FASTInputDeck
from welib.weio.pickle_file import PickleFile

from welib.tools.signal_analysis import zero_crossings


class TabulatedWSEstimatorFloating(TabulatedWSEstimatorBase):
    """ 
    WS, RPM, pitch, phiy

    """
    def __init__(self, R=None, rho=1.225, fstFile=None, pickleFile=None):
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

        # Data
        self.pickleFile = None

        if pickleFile:
            self.loadPickle(pickleFile)

    def loadPickle(self, pickleFile):
        self.pickleFile = pickleFile
        pkl = PickleFile(pickleFile)

        # --- Operating conditions
        try:
            OP = pkl['Oper']
        except:
            OP = pkl['OP']
        print(OP)
        #from welib.fast.tools.lin import matToSIunits
        #OP = matToSIunits(OP, name='OP', verbose=False, row=False, col=True)
        #print(OP)
        OP = self._sanitizeOP(OP, expectedCols=['WS_[m/s]', 'Pitch_[deg]', 'RotSpeed_[rpm]', 'PhiY_[deg]'], onlyExpected=False)

        print(OP)
        self.OP=OP
        self.WS_op   =OP['WS_[m/s]'].values
        self.omega_op=OP['RotSpeed_[rpm]'].values*2*np.pi/60 # [rad/s]
        self.omegaRated=np.max(self.omega_op)
#         self.omegaLow  =0.4*self.omegaRated
        self.WSRated=np.interp(self.omegaRated*0.98, self.omega_op, self.WS_op)
#             self.WSCutOff=28

        # ---
        self.WS     = pkl['WS']
        self.pitch  = pkl['Pitch']
        self.omega  = pkl['RPM' ]*np.pi/30
        self.phiy   = pkl['PhiY']
        self.CP     = pkl['CP']
        self.CT     = pkl['CT']
        self.CP[np.isnan(self.CP)]=0
        self.CP[self.CP<0]=0
        self.CT[np.isnan(self.CT)]=0
        self.CT[self.CT<0]=0

        # --- Computing weights
        MWS     = pkl['MWS']
        Momega  = pkl['MOmega']*np.pi/30
        print(pkl.keys())
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


        # ---  Compute interpolated values at Operating points to be consistent
        WS    = self.OP['WS_[m/s]'].values[:]
        Omega = self.OP['RotSpeed_[rpm]'].values[:] *np.pi/30
        pitch = self.OP['Pitch_[deg]'].values[:]
        phiy  = self.OP['PhiY_[deg]'].values[:]
        self.OP['Paero_i_[W]'] = self.Power(WS, Omega, pitch, phiy) # Interpolated power
        Q2                     = self.Torque(WS, Omega, pitch, phiy) # Interpolated torque
        Q1                      = self.OP['Paero_i_[W]']/Omega
        self.OP['Qaero_i_[Nm]'] = self.OP['Paero_i_[W]']/Omega
        print(Q1)
        print(Q2)
        print(Q1-Q2)
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(WS, Q1 ,'-'   , label='')
        ax.plot(WS, Q2 ,'--'   , label='')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.legend()
#         plt.show()
#         import pdb; pdb.set_trace()

    def computeWeights(self, P, Q, T):
        # Compute interpolants
        self.fP=si.RegularGridInterpolator( (self.WS, self.omega, self.pitch, self.phiy), P, fill_value=np.nan, bounds_error=False)
        self.fQ=si.RegularGridInterpolator( (self.WS, self.omega, self.pitch, self.phiy), Q, fill_value=np.nan, bounds_error=False)
        if self.CT is not None:
            self.fT=si.RegularGridInterpolator( (self.WS, self.omega, self.pitch, self.phiy), T, fill_value=np.nan, bounds_error=False)

        # --- Lambda at higher res
        #LambdaMid = self.Lambda[:-1] + np.diff(self.Lambda)/2
        #Lambda = np.sort(np.concatenate((self.Lambda, LambdaMid)))
        #self.LambdaHR = Lambda
        WS = self.OP['WS_[m/s]'].values
        #WS = self.WS
        self.WSHR = np.linspace(WS[0], WS[-1], 2*len(WS))

    def Power(self, WS, omega, pitch, phiy):
        return self.fP( (WS, omega, pitch, phiy) )

    def Thrust(self, WS, omega, pitch, phiy):
        return self.fT( (WS, omega, pitch, phiy) )

    def Torque(self, WS, omega, pitch, phiy):
        return self.fQ( (WS, omega, pitch, phiy) )

    def TorqueAt(self, omega, pitch, phiy):
        """ 
        Return Torque(WS) curve for a given pitch and rotational speed
        pitch,omega: scalar
        """
        WS     = self.WSHR
        WS = WS[WS<self.WSmax]
        vpitch = np.array([pitch]*len(WS))
        vomega = np.array([omega]*len(WS))
        vphiy  = np.array([phiy]*len(WS))
        return WS, self.Torque(WS, omega=vomega, pitch=vpitch, phiy=vphiy)

    def PowerAt(self, omega, pitch, phiy):
        """ 
        Return Power(WS) curve for a given pitch and rotational speed
        pitch,omega: scalar
        """
        WS     = self.WSHR
        WS = WS[WS<self.WSmax]
        vpitch = np.array([pitch]*len(WS))
        vomega = np.array([omega]*len(WS))
        vphiy  = np.array([phiy]*len(WS))
        return WS, self.Power(WS, omega=vomega, pitch=vpitch, phiy=vphiy)


    def estimate(self, Qa, omega, pitch,  phiy , WS0, relaxation=0, WSavg=None, debug=False, method='min', deltaWSMax=1, verbose=False): 
        """
        INPUTS:
         - Qa: aerodynamic power [W]
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
            info['Qa']=Qa; info['pitch']=pitch; info['omega']=omega; info['phiy']=phiy; info['WS0']=WS0
            info['relaxation']=relaxation; info['method']=method; info['deltaWSMax']=deltaWSMax;


        # --- Use operating conditions
        if method.find('oper')>=0:
            if self.OP is None:
                raise Exception('Cannot use method `oper`, operFile was not provided')
            WSOP    = self.OP['WS_[m/s]'].values   #[m/s]
            QaeroOP = self.OP['Qaero_i_[Nm]'].values 
            PitchOP = self.OP['Pitch_[deg]'].values #[deg]
            PhiYOP  = self.OP['PhiY_[deg]'].values #[deg]
            # np.interp(Qa, self.OP['RtAeroMxh_[N-m]'], self.OP['WS_[m/s]'])
            WScrossOP, _, _ = zero_crossings(QaeroOP-Qa, x=WSOP) #, bouncingZero=True)
            if len(WScrossOP)==0:
                # Can happen if torque below minimum or above maximum torque
                if verbose:
                    print('>>> No crossing OP', pitch, omega, Qa)
                    print('{} OPcross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossOP), pitch, omega, Qa, WS0), WScrossOP)
                #import pdb; pdb.set_trace()
                #import matplotlib.pyplot as plt
                #fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
                #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
                #ax.plot(WSOP, QaeroOP, label='Omega = {:3.1f}'.format(omega))
                #ax.plot(WSOP, Qa+WSOP*0)
                #ax.set_xlabel('WS [m/s]')
                #ax.set_ylabel('Q [N]')
                ##ax.set_xlim([0,30])
                ##ax.set_ylim([0,2e7])
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
                #ax.plot(WSOP, QaeroOP, label='Omega = {:3.1f}'.format(omega))
                #ax.plot(WSOP, Qa+WSOP*0)
                #ax.set_xlabel('WS [m/s]')
                #ax.set_ylabel('Q [N]')
                #ax.set_xlim([0,30])
                #ax.set_ylim([0,2e7])
                #ax.legend()
                #plt.show()
                #raise Exception('TODO')
            WS_est = WSoper
            if debug:
                info['WS_oper'] = WSoper
                info['WS_crossOP'] = WScrossOP

        # --- Find torque for given pitch and omega
        if method.find('crossing')>=0:
            if omega>0:
                iNear = None
                vWS, vQ = self.TorqueAt(omega=omega, pitch=pitch, phiy=phiy)
                WScrossI, iBef, sign = zero_crossings(vQ-Qa, x=vWS)
                if len(WScrossI)==0:
                    if verbose:
                        print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
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
                    WS_I = WS_est # WS0
                elif len(WScrossI)==1:
                    if verbose:
                        print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
                    WS_I = WScrossI[0]
                elif len(WScrossI)==2:
                    if verbose:
                        print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
                    if WSoper is not None:
                        iNear = np.argmin(abs(WScrossI-WSoper))
                    else:
                        iNear = np.argmin(abs(WScrossI-WS0))
                    WS_I = WScrossI[iNear] # TODO
                else:
                    if verbose:
                        print('{} cross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
                    # Find WS based on OP
                    if WSoper is not None:
                        iNear = np.argmin(abs(WScrossI-WSoper))
                    else:
                        iNear = np.argmin(abs(WScrossI-WS0))
                    WS_I = WScrossI[iNear] # TODO
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

                # --- If estimate is far away, change it
                # We use the oper and previous point as a guess
                if WSoper is not None:
                    WS_guess = (WSoper+WS0)/2
                else:
                    WS_guess = WS0
                if abs(WS_I - WS_guess) > deltaWSMax:
                    #print('>>>', WS_I)
                    p0 = (WS0, Qa) # Previous point, at current torque
                    dist = np.sqrt( ((vWS-WS_guess)/WS_guess)**2 + ((vQ-Qa)/Qa)**2 )
                    i = np.argmin(dist)
                    WS_est = vWS[i]
                else:
                    WS_est = WS_I

                if debug:
                    info['WS_crossI'] = WScrossI
                    info['WS_I']      = WS_I
                    info['iNear'] = iNear
                    #info['cross_vQ'] = vQ
                    #info['cross_vWS'] = vWS


        WS = WS0*relaxation + (1-relaxation)*WS_est
        if debug:
            info['WS_est'] = WS_est
            info['WS'] = WS

        self._debug_info=info

        return WS, info

    def debugPlot(self, info=None, HR=False):
        from welib.tools.colors import python_colors

        if info is None:
            info = self._debug_info
        for k,v in info.items():
            print('{:15s}={}'.format(k,v))


        pitch = info['pitch']
        omega = info['omega']
        phiy  = info['phiy']
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

        WS0, Q0 = self.TorqueAt(pitch=pitch     , omega=omega     , phiy=phiy)
        WS1, Q1 = self.TorqueAt(pitch=pitch     , omega=omega*0.95, phiy=phiy)
        WS2, Q2 = self.TorqueAt(pitch=pitch     , omega=omega*1.05, phiy=phiy)
        WS3, Q3 = self.TorqueAt(pitch=pitch     , omega=omega     , phiy=phiy*0.95)
        WS4, Q4 = self.TorqueAt(pitch=pitch     , omega=omega     , phiy=phiy*1.05)
        WS5, Q5 = self.TorqueAt(pitch=pitch*0.95, omega=omega     , phiy=phiy)
        WS6, Q6 = self.TorqueAt(pitch=pitch*1.05, omega=omega     , phiy=phiy)


        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(WS0   , Q0      , '-+', c=python_colors(0),       label='Interp. data')
        ax.plot(WS1   , Q1      , '--', c=python_colors(0), alpha=0.5)
        ax.plot(WS2   , Q2      , '--', c=python_colors(0), alpha=0.5)
        ax.plot(WS2   , Q3      , '--', c=python_colors(0), alpha=0.5)
        ax.plot(WS2   , Q4      , '--', c=python_colors(0), alpha=0.5)
        ax.plot(WS2   , Q5      , '--', c=python_colors(0), alpha=0.5)
        ax.plot(WS2   , Q6      , '--', c=python_colors(0), alpha=0.5)
        if WS_ref is not None:
            ax.plot(WS_ref, Qa      , 'ko',       label='Ref')
        if WS_est is not None:
            ax.plot(WS_est, Qeval1  , 'd', c=python_colors(1) , ms=5, label='Est.')
        ax.plot(WS0   , WS0*0+Qa, 'k--',      label='Qa')
        ax.plot([WS_guess, WS_guess], [np.min(Q0), np.max(Q0)], 'k:', label='WS guess/prev')

        if 'WS_crossI' in info.keys():
            col=python_colors(3)
            n=len(info['WS_crossI'])
            for iWS,WS in enumerate(info['WS_crossI']):
                ax.plot(WS      , Qa, 'o' , c=col, ms=6, label='Interp. Cross (n={})'.format(n) if iWS==0 else None)
            if n>1:
                ax.plot(info['WS_I'], Qa, 'o',  c=col, ms=11, label='Interp. Cross selected', markerfacecolor='none')

        if 'WS_crossOP' in info.keys():
            col = python_colors(4)
            WSOP    = self.OP['WS_[m/s]'].values   #[m/s]
            QaeroOP = self.OP['Qaero_i_[Nm]'].values #[Nm]
            ax.plot(WSOP, QaeroOP, ':', label='OP', c=col)
            n=len(info['WS_crossOP'])
            for iWS,WS in enumerate(info['WS_crossOP']):
                ax.plot(WS,          Qa, '^' , c=col, ms=6 , label='OP Cross (n={})'.format(n) if iWS==0 else None)
            if n>1:
                ax.plot(info['WS_oper'], Qa, '^' , c=col, ms=11 , label='OP Cross selected', markerfacecolor='none')

        col = 'k'
        ax.plot(info['WS_est'], Qa, 's',  c=col, ms=3, label='WS est')
        ax.plot(info['WS_est'], Qa, 's',  c=col, ms=10, label='WS est (with relaxation)', markerfacecolor='none')



        if HR:
            nHR=1000
            WS_HR = np.linspace(WS0[0],WS0[-1],nHR)
            Qeval_HR = self.Torque(WS_HR, [pitch]*nHR, [omega]*nHR)
            ax.plot(WS_HR, Qeval_HR, 'k-', lw=0.5, label='High Res eval')
        ax.set_xlabel('Wind Speed [m/s]')
        ax.set_ylabel('Torque [N]')
        ax.legend(fontsize=9)
        return fig

    def __repr__(self):
        s=''
        s+='<ws_estimator.TabulatedWSEstimator object> \n'
        s+=' - WS : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.WS),np.max(self.WS),self.WS[1]-self.WS[0], len(self.WS))
        s+=' - omega : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.omega),np.max(self.omega),self.omega[1]-self.omega[0], len(self.omega))
        s+=' - pitch  : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.pitch) ,np.max(self.pitch) ,self.pitch[1]-self.pitch[0]  , len(self.pitch))
        s+=' - CP     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CP),np.max(self.CP),self.CP.shape[0],self.CP.shape[1])
        if self.CT is not None:
            s+=' - CT     : [min={:8.3f}, max={:8.3f}, n={}x{}]  \n'.format(np.min(self.CT),np.max(self.CT),self.CT.shape[0],self.CT.shape[1])
        s+=' - R      : {}  \n'.format(self.R)
        s+=' - rho    : {}  \n'.format(self.rho)
#         s+=' - omegaLow:{}  \n'.format(self.omegaLow)
#         # files
#         s+=' - fstFile:     {}  \n'.format(self.fstFile)
#         s+=' - aeroMapFile: {}  \n'.format(self.aeroMapFile)
        s+=' - pickleFile:    {}  \n'.format(self.pickleFile)
        return s
