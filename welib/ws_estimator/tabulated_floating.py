import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.optimize import minimize_scalar


from welib.ws_estimator.tabulated import TabulatedWSEstimatorBase

import welib.weio as weio
from welib.weio.fast_input_deck import FASTInputDeck
from welib.weio.pickle_file import PickleFile
from welib.tools.signal_analysis import zero_crossings
from welib.tools.dictlib import renameDictKey
from welib.yams.models.simulator import _loadOFOut


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
        # --- Sanitize pickle
        pkl = renameDictKey(pkl, ['Oper' ,'OP'] , 'OP')
        pkl = renameDictKey(pkl, ['Operi','OPi'], 'OPi')
        self.pkl=pkl
        # --- Operating conditions
        OP = pkl['OP']
        #from welib.fast.tools.lin import matToSIunits
        #OP = matToSIunits(OP, name='OP', verbose=False, row=False, col=True)
        OP = self._sanitizeOP(OP, expectedCols=['WS_[m/s]', 'Pitch_[deg]', 'RotSpeed_[rpm]', 'PhiY_[deg]'], onlyExpected=False)

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
        #print(pkl.keys())
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
#         fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#         ax.plot(WS, Q1 ,'-'   , label='')
#         ax.plot(WS, Q2 ,'--'   , label='')
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.legend()
#         plt.show()

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


    def clip(self, omega, pitch, phiy, verbose=False):
        """ Bound values so that they remain within sampled data bounds """
        omega_ = omega
        pitch_ = pitch
        phiy_  = phiy
        clipped = False
        # --- Bounding
        if omega_<self.omega[0]:
            omega_ = self.omega[0]
            clipped = True
            if verbose:
                print('omegaLow for p={:8.3f} om={:8.3} ph={:8.3f}'.format(pitch, omega, phiy))
        if pitch_<self.pitch[0]:
            pitch_ = self.pitch[0]
            clipped = True
            if verbose:
                print('pitchLow for p={:8.3f} om={:8.3} ph={:8.3f}'.format(pitch, omega, phiy))
        if phiy_<self.phiy[0]:
            phiy_ = self.phiy[0]
            clipped = True
            if verbose:
                print('phiyLow  for p={:8.3f} om={:8.3} ph={:8.3f}'.format(pitch, omega, phiy))
        if omega_>self.omega[-1]:
            omega_ = self.omega[-1]
            clipped = True
            if verbose:
                print('omegaHig for p={:8.3f} om={:8.3} ph={:8.3f}'.format(pitch, omega, phiy))
        if pitch_>self.pitch[-1]:
            pitch_ = self.pitch[-1]
            clipped = True
            if verbose:
                print('pitchHig for p={:8.3f} om={:8.3} ph={:8.3f}'.format(pitch, omega, phiy))
        if phiy_>self.phiy[-1]:
            phiy_ = self.phiy[-1]
            clipped = True
            if verbose:
                print('phiyHig  for p={:8.3f} om={:8.3} ph={:8.3f}'.format(pitch, omega, phiy))
        return omega_, pitch_, phiy_, clipped


    def estimate(self, Qa, omega, pitch,  phiy , WS0, relaxation=0, method='crossing', deltaWSMax=1, verbose=False, debug=False, t=0, WSref=np.nan): 
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

        def saveState(info=None):
            if info is None:
                info={}
            # Store state
            info['Qa']=Qa; info['pitch']=pitch; info['omega']=omega; info['phiy']=phiy; info['WS0']=WS0
            info['relaxation']=relaxation; info['method']=method; info['deltaWSMax']=deltaWSMax;
            info['t']=t
            return info

        if debug:
            info = saveState()

        # --- Use operating conditions
        if method.find('oper')>=0:
            if self.OP is None:
                raise Exception('Cannot use method `oper`, operFile was not provided')
            WSOP    = self.OP['WS_[m/s]'].values   #[m/s]
            QaeroOP = self.OP['Qaero_i_[Nm]'].values 
            PitchOP = self.OP['Pitch_[deg]'].values #[deg]
            PhiYOP  = self.OP['PhiY_[deg]'].values #[deg]
            # np.interp(Qa, self.OP['RtFldMxh_[N-m]'], self.OP['WS_[m/s]'])
            WScrossOP, _, _ = zero_crossings(QaeroOP-Qa, x=WSOP) #, bouncingZero=True)
            if len(WScrossOP)==0:
                # Can happen if torque below minimum or above maximum torque
                if verbose:
                    print('>>> No crossing OP', pitch, omega, Qa)
                    print('{} OPcross p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossOP), pitch, omega, Qa, WS0), WScrossOP)
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
                # Bound values
                omega_, pitch_, phiy_, clipped = self.clip(omega, pitch, phiy, verbose=verbose)
                # 
                vWS, vQ = self.TorqueAt(omega=omega_, pitch=pitch_, phiy=phiy_)
                if np.all(np.isnan(vQ)):
                    print('ICross - All NaN for p={:8.3f} om={:8.3} ph={:8.3f} (WSref={})'.format(pitch, omega, phiy, WSref))
                    print('>>>TODO')
                WScrossI, iBef, sign = zero_crossings(vQ-Qa, x=vWS)
                if len(WScrossI)==0:
                    if verbose:
                        print('{} Icross  p={:8.3f} om={:8.3f} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
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
                        print('{} Icross  p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
                    WS_I = WScrossI[0]
                elif len(WScrossI)==2:
                    if verbose:
                        print('{} Icross  p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
                    if WSoper is not None:
                        iNear = np.argmin(abs(WScrossI-WSoper))
                    else:
                        iNear = np.argmin(abs(WScrossI-WS0))
                    WS_I = WScrossI[iNear] # TODO
                else:
                    if verbose:
                        print('{} Icross  p={:8.3f} om={:8.3} Qa={:10.2f} WS0={:7.3f}'.format(len(WScrossI), pitch, omega, Qa, WS0), WScrossI)
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
                    if verbose:
                        print('Icross, jump too big')
                    dist = np.sqrt( ((vWS-WS_guess)/WS_guess)**2 + ((vQ-Qa)/Qa)**2 )
                    try:
                        i = np.nanargmin(dist)
                        WS_est = vWS[i]
                    except ValueError:
                        if verbose:
                            print('Giving up')
                            WS_est=WS_guess
                else:
                    WS_est = WS_I

                if debug:
                    info['WS_crossI'] = WScrossI
                    info['WS_I']      = WS_I
                    info['iNear'] = iNear
                    #info['cross_vQ'] = vQ
                    #info['cross_vWS'] = vWS
        try: 
            WS_est
        except NameError:
            raise Exception('Invalid method {}'.format(method))
#         if (WS0-WS_est)>5:
#             if debug is False:
#                 info = saveState()
#                 self.debugRerun(info)



        WS = WS0*relaxation + (1-relaxation)*WS_est
        if debug:
            info['WS_est'] = WS_est
            info['WS'] = WS

        self._debug_info=info

        return WS, info



    def estimateTimeSeries(self, Qaero, omega, pitch, phiy, WS_prev=None, WS_ref=None, debug=False, time=None, **kwargs):
        """ 
        Perform wind speed estimation given a time series of aerodynamic torque, pitch and rotational speed
        """
        print('Estimating WS on time series...')
        WS_est = np.zeros(omega.shape)
        if WS_prev is None:
            WS_prev = 1
        ts_info = None
        if debug:
            # Storage for debug
            pass
        if time is None:
            time = np.arange(len(Qaero))
        if WS_ref is None:
            WSref = time*np.nan
        else:
            WSref = np.array(WS_ref)
        for i,(Qa, om, pit, ph, t) in enumerate(zip(Qaero, omega, pitch, phiy, time)):
            ws_hat, info    = self.estimate(Qa, omega=om, pitch=pit, phiy=ph, WS0=WS_prev, debug=debug, t=t, WSref=WS_ref[i], **kwargs)
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


    def validValues(self, windspeed, omega, pitch, phiy):
        """ Return boolean array of valid values"""
        def findSurroundingIndices(x,v):
            i=np.argmin(np.abs(v-x)) 
            if v[i]>x and i!=0: 
                i=i-1
            if v[i]==x and i==len(v): 
                j=i
                i=i-1
            if i+1<len(v):
                j = i+1
            else:
                j=i
            return i,j
        bValidMat = self.pkl['bValid']
        bValid = np.array([True]*len(windspeed))
        for i, (ws, om, pi, ph) in enumerate(zip(windspeed, omega, pitch,phiy)):
            # Find closest surrounding indices
            IWS = findSurroundingIndices(ws, self.WS)
            IOM = findSurroundingIndices(om, self.omega)
            IPI = findSurroundingIndices(pi, self.pitch)
            IPH = findSurroundingIndices(ph, self.phiy)
            b=True
            for iws in IWS:
                for iom in IOM:
                    for ipi in IPI:
                        for iph in IPH:
                            b=b and bValidMat[iws, iom, ipi, iph]
            bValid[i] = b
        return bValid

    def estimateTimeSeriesFromOF(self, fstFilename, tRange=None, **kwargs):
        """" 
         - tRange: tuple (tmin, tmax) to limit the time used
        """
        df, time = _loadOFOut(fstFilename, tRange=tRange)
        WS_ref     = df['RtVAvgxh_[m/s]'].values # Rotor avg
        Pitch      = df['BldPitch1_[deg]'].values
        Qaero_ref  = df['RtFldMxh_[N-m]'].values
        Omega      = df['RotSpeed_[rpm]'].values*2*np.pi/60 # rad/s
        PhiY       = df['PtfmPitch_[deg]'].values
        # Estimating wind speed on time series
        WS_est, ts_info = self.estimateTimeSeries(Qaero_ref, omega=Omega, pitch=Pitch, phiy=PhiY, WS_prev=WS_ref[0]*0.9, WS_ref=WS_ref, time=time, **kwargs)
        # Evaluating torque
        Qaero_eval = self.Torque(WS_ref, omega=Omega, pitch=Pitch, phiy=PhiY)
        Qaero_est  = self.Torque(WS_est, omega=Omega, pitch=Pitch, phiy=PhiY)
        # Monitoring if nan occured
        bNaN = np.isnan(Qaero_eval)
        if sum(bNaN)>0:
            fig = self.plotValidSamples()
            axes=fig.axes
            axes[0].plot(WS_ref[bNaN], Pitch[bNaN], 'k.', label='NaN')
            axes[1].plot(WS_ref[bNaN], Omega[bNaN], 'k.', label='NaN')
            axes[2].plot(WS_ref[bNaN], PhiY[bNaN] , 'k.', label='NaN')
            axes[0].legend()
        # Storing data into a dataframe
        M    = np.column_stack((time, WS_ref, WS_est, Qaero_ref, Qaero_eval, Qaero_est, Omega, Pitch, PhiY))
        cols = ['Time_[s]','WS_ref_[m/s]','WS_est_[m/s]','Qaero_ref_[N]','Qaero_eval_[N]','Qaero_est_[N]','Omega_[rad/s]','Pitch_[deg]', 'PtfmPitch_[deg]']
        dfOut = pd.DataFrame(data=M, columns=cols)
        return dfOut

    def debugOP(self, omega, pitch, phiy, ws=None):
        """
        Look at raw data and interpolated around a given an operating point
        """

        WS_bef, Q_bef = self.TorqueAt(omega=omega, pitch=pitch, phiy=phiy)
        if np.all(np.isnan(Q_bef)):
            print('>>> All Q are NaN')

        omega, pitch, phiy, clipped = self.clip(omega, pitch, phiy, verbose=True)

        if clipped:
            print('>>> Clipping was needed!')

        def findSurroundingIndices(x,v):
            i=np.argmin(np.abs(v-x)) 
            if v[i]>x and i!=0: 
                i=i-1
            if v[i]==x and i==len(v): 
                j=i
                i=i-1
            if i+1<len(v):
                j = i+1
            else:
                j=i
            return i,j

        # Find closest surrounding indices
        iom,jom = findSurroundingIndices(omega, self.omega)
        ipi,jpi = findSurroundingIndices(pitch, self.pitch)
        ipy,jpy = findSurroundingIndices(phiy, self.phiy)
        if ws is not None:
            iws,jws = findSurroundingIndices(ws, self.WS)
        print('omega: {:8.3f} <= {:8.3f} <= {:8.3f} - index: {} {}'.format(self.omega[iom], omega, self.omega[jom], iom, jom ))
        print('RPM  : {:8.3f} <= {:8.3f} <= {:8.3f} - index: {} {}'.format(self.omega[iom]*30/np.pi, omega*30/np.pi, self.omega[jom]*30/np.pi, iom, jom))
        print('pitch: {:8.3f} <= {:8.3f} <= {:8.3f} - index: {} {}'.format(self.pitch[ipi], pitch, self.pitch[jpi], ipi, jpi ))
        print('phiy : {:8.3f} <= {:8.3f} <= {:8.3f} - index: {} {}'.format(self.phiy [ipy], phiy , self.phiy [jpy], ipy, jpy))
        if ws is not None:
            print('ws   : {:8.3f} <= {:8.3f} <= {:8.3f} - index: {} {}'.format(self.WS [iws], ws , self.WS [jws], iws, jws))

        WS1, P1 = self.PowerAt (omega=omega, pitch=pitch, phiy=phiy)
        WS1, Q1 = self.TorqueAt(omega=omega, pitch=pitch, phiy=phiy)

        # --- Plot nearest raw data, OP, and interpolated
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(self.WS, self.Q[:, iom , ipi , ipy]  , 'k-', label='raw')
        if ws is not None:
            ax.plot(self.WS[iw], self.Q[iw, iom , ipi , ipy]  , 'ko')
        ax.plot(self.WS[:], self.Q[:, jom, ipi, ipy], 'k--', alpha=0.1)
        ax.plot(self.WS[:], self.Q[:, iom, jpi, ipy], 'k--', alpha=0.1)
        ax.plot(self.WS[:], self.Q[:, iom, ipi, jpy], 'k--', alpha=0.1)
        ax.plot(self.WS[:], self.Q[:, jom, jpi, jpy], 'k--', alpha=0.1)
        ax.plot(WS1, Q1 , label='Interpolated')
        ax.plot(self.OP['WS_[m/s]'], self.OP['Qaero_i_[Nm]'], ':', label='OP')
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('')
        ax.legend()
        return fig


    def plotValidSamples(self):
        """ """
        from welib.tools.colors import python_colors
        pkl = self.pkl
        Operi  = pkl['OPi']
        MWS    = pkl['MWS']
        MOmega = pkl['MOmega']
        MPitch = pkl['MPitch']
        MPhi   = pkl['MPhi']
        bValid = pkl['bValid']
        alpha=0.02
        fig,axes = plt.subplots(1, 3, sharey=False, figsize=(12.4,3.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.30)
        ax = axes[0]
        ax.plot(Operi['WS_[m/s]'], Operi['Pitch_[deg]'],     'k-'    , label='Operating conditions')
        ax.plot(MWS[bValid].flatten(), MPitch[bValid].flatten(), 'o'   , alpha=0.02)
        ax.plot(10,0,'o', alpha=0.5, label='Sampled sim.', c=python_colors(0))
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('Pitch [deg]')
        ax.legend()
        ax = axes[1]
        ax.plot(Operi['WS_[m/s]'], Operi['RPM_[rpm]'],     'k-'    , label='Reference')
        ax.plot(MWS[bValid].flatten(), MOmega[bValid].flatten(), 'o'    , label='Points', alpha=0.02)
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('Rotor speed [rpm]')
        ax = axes[2]
        ax.plot(Operi['WS_[m/s]'], Operi['PhiY_[deg]'],     'k-'    , label='Reference')
        ax.plot(MWS[bValid].flatten(), MPhi[bValid].flatten(), 'o'    , label='Points', alpha=0.02)
        ax.set_xlabel('Wind speed [m/s]')
        ax.set_ylabel('Platform Pitch [deg]')
        return fig





    def debugPrint(self, info):
        print('-------------------- DEBUG DUMP -----------------')
        for k,v in info.items():
            print('{:15s}={}'.format(k,v))

    def debugRerun(self, info, plot=True):
        Qa         = info['Qa']
        omega      = info['omega']
        pitch      = info['pitch']
        phiy       = info['phiy']
        WS0        = info['WS0']
        relaxation = info['relaxation']
        method     = info['method']
        deltaWSMax = info['deltaWSMax']
        self.debugPrint(info)
        ws_est, info2 = self.estimate(Qa, omega, pitch,  phiy , WS0, relaxation=relaxation, method=method, deltaWSMax=deltaWSMax, debug=True, verbose=True) 
        if plot:
            self.debugPlot(info2)


    def debugPlot(self, info=None, HR=False):
        from welib.tools.colors import python_colors

        if info is None:
            info = self._debug_info
        self.debugPrint(info)


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
        s+=' - WS     : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.WS),np.max(self.WS),self.WS[1]-self.WS[0], len(self.WS))
        s+=' - omega  : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.omega),np.max(self.omega),self.omega[1]-self.omega[0], len(self.omega))
        s+=' - pitch  : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.pitch) ,np.max(self.pitch) ,self.pitch[1]-self.pitch[0]  , len(self.pitch))
        s+=' - phiy   : [min={:8.3f}, max={:8.3f}, delta={:8.4f}, n={}]  \n'.format(np.min(self.phiy) ,np.max(self.phiy) ,self.phiy[1]-self.phiy[0]  , len(self.phiy))
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
