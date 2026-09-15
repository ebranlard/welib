"""
Tools for ServoDyn

- Plot the simple generator model
- Fit a simple generator model


"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from welib.tools.strings import WARN
from welib.weio.fast_input_file import FASTInputFile
from welib.weio.rosco_discon_file import ROSCODISCONFile
# from pydatview.tools.curve_fitting import gentorque, GeneratorTorqueFitter, model_fit
from welib.tools.curve_fitting import gentorque, GeneratorTorqueFitter

RPM2RADS=   2*np.pi/60;
RADS2RPM=1/(2*np.pi/60);

class ServoDyn:

    def __init__(self, svdFilename_or_data=None, TP=None, load_discon=True):
        """ 
        Initialize a ServoDyn object either with:
          - svdFilename: a servody input file name
          - svdData: an instance of FASTInputFile
        """

        # --- Data
        self.File=None
        self.DISCON=None
        self._data_rpm    = None # User Data / measurements
        self._data_torque = None # User Data / measurements
        self._fit_rpm     = None # Fit from user data
        self._fit_torque  = None # Fit from user data

        # Read SubDyn file
        if svdFilename_or_data is not None:
            if hasattr(svdFilename_or_data,'startswith'): # if string
                self.File = FASTInputFile(svdFilename_or_data)
            else:
                self.File = svdFilename_or_data

            if load_discon:
                discon_in_file = os.path.normpath(self.File['DLL_InFile'].strip('"'))
                if not os.path.isabs(discon_in_file):
                    discon_in_path = os.path.join( os.path.dirname(self.File.filename) , discon_in_file)
                else:
                    discon_in_path = discon_in_file
                if os.path.exists(discon_in_path):
                    self.DISCON = ROSCODISCONFile(discon_in_path)


    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|properties:\n'
        s+='|* VS_RtGnSp: {} [rpm]\n'.format(self.VS_RtGnSp)
        s+='|* VS_RtTq:   {} [Nm]\n'.format(self.VS_RtTq)
        s+='|* VS_Rgn2K:  {} [N-m/rpm^2]\n'.format(self.VS_Rgn2K)
        s+='|* VS_SlPc:   {} [%]\n'.format(self.VS_SlPc)
        s+='|- File: (input file data)\n'
        s+=f'|- DISCON: is None? {self.DISCON is None}\n'
        s+='|methods:\n'
        s+='|- VS_print()\n'
        s+='|- VS_plot\n'
        s+='|- VS_dataframe\n'
        return s

    # --------------------------------------------------------------------------------}
    # --- VS properties 
    # --------------------------------------------------------------------------------{
#         if discon: 
#             if self.DISCON is None:
#                 raise Exception('Cannot get DISCON parameters, file was not read')
#             Rgn2K   /= RADS2RPM**2              # [N-m/rpm^2] 
#             SlPc=0 
# # !------- VS TORQUE CONTROL ------------------------------------------------
# # 94.40000            ! VS_GenEff			- Generator efficiency mechanical power -> electrical power, [should match the efficiency defined in the generator properties!], [%]
# # 4.30935e+04         ! VS_ArSatTq		- Above rated generator torque PI control saturation, [Nm]
# # 4.00000e+04         ! VS_MaxRat			- Maximum torque rate (in absolute value) in torque controller, [Nm/s].
# # 4.74029e+04         ! VS_MaxTq			- Maximum generator torque in Region 3 (HSS side), [Nm].
# # 0.00000e+00         ! VS_MinTq			- Minimum generator torque (HSS side), [Nm].
# # 34.64286            ! VS_MinOMSpd		- Minimum generator speed [rad/s]
# # 5.00000e+06         ! VS_RtPwr			- Wind turbine rated power [W]
# # 1                   ! VS_n				- Number of generator PI torque controller gains
# # -6.97771e+02        ! VS_KP				- Proportional gain for generator PI torque controller [-]. (Only used in the transitional 2.5 region if VS_ControlMode =/ 2)
# # -1.04507e+02        ! VS_KI				- Integral gain for generator PI torque controller [s]. (Only used in the transitional 2.5 region if VS_ControlMode =/ 2)
# # 7.50000             ! VS_TSRopt		    - Power-maximizing region 2 tip-speed-ratio. Only used in VS_ControlMode = 2.
# 
#         else:
#             GenModel = self.File['GenModel']  # Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]
#             GenEff   = self.File['GenEff']    # Generator efficiency [ignored by the Thevenin and user-defined generator models] (%)
#             GenTiStr = self.File['GenTiStr']  # Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)
#             GenTiStp = self.File['GenTiStp']  # Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)
#             SpdGenOn = self.File['SpdGenOn']  # Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]
#             TimGenOn = self.File['TimGenOn']  # Time to turn on the generator for a startup (s) [used only when GenTiStr=True]
#             TimGenOf = self.File['TimGenOf']  # Time to turn off the generator (s) [used only when GenTiStp=True]
#             SlPc     = self.File['VS_SlPc']   # Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) [used only when VSContrl=1]
    @property
    def VS_RtGnSp(self):
        """Rated generator speed [rpm]"""
        if self.DISCON is not None:
            return self.DISCON['VS_RefSpd'] * RADS2RPM # NOTE: discon is in rad/s. Rated generator speed [rad/s]
        return self.File['VS_RtGnSp']                  # [RPM] Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]

    @property
    def VS_RtTq(self):
        """Rated generator torque [Nm]"""
        if self.DISCON is not None:
            return self.DISCON['VS_RtTq'] # Rated torque, [Nm].
        return self.File['VS_RtTq']       # Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]

    @property
    def VS_Rgn2K(self):
        """Generator torque constant in Region 2 [N-m/rpm^2]"""
        if self.DISCON is not None:
            return self.DISCON['VS_Rgn2K'] / (RADS2RPM**2) # [N-m/(rad/s)^2] NOTE: different than OpenFAST. Generator torque constant in Region 2 (HSS side). Only used in VS_ControlMode = 1,3,4
        return self.File['VS_Rgn2K'] # [N-m/rpm^2] Generator torque constant in Region 2 for simple variable-speed generator control (HSS side)  [used only when VSContrl=1]



    @property
    def VS_SlPc(self):
        """Rated generator slip percentage [%]"""
        if self.DISCON is not None:
            return 0.1
        return self.File['VS_SlPc']

    def VS_print(self):
        """Print variable speed control parameters"""

        RtGnSp = self.VS_RtGnSp
        RtTq   = self.VS_RtTq
        Rgn2K  = self.VS_Rgn2K
        SlPc   = self.VS_SlPc

        if self.DISCON is None:
            SpdGenOn = self.File['SpdGenOn']
            print('SpdGenOn : {:20.3f} [rpm]'.format(SpdGenOn))

        print('VS_RtGnSp: {:20.3f} [rpm],       {:15.3f} [rad/s]'.format(RtGnSp, RtGnSp * RPM2RADS))
        print('VS_RtTq  : {:20.3f} [Nm]'.format(RtTq))
        print('VS_Rgn2K : {:20.3f} [N-m/rpm^2], {:15.3f} [N-m/(rad/s)^2]'.format(Rgn2K, Rgn2K * (RADS2RPM**2)))
        print('VS_SlPc  : {:20.3f} [%]'.format(SlPc))
        RtTqCheck = Rgn2K*RtGnSp**2
        if RtTqCheck>RtTq:
            print('TqCheck  : {:20.3f} [Nm] <? {:20.3f}'.format(RtTqCheck, RtTq))
            WARN(' Rgn2K*RtGnSp**2 ({})  > RtTq ({})'.format(RtTqCheck, RtTq))


    def VS_DataFrame(self, rpm=None, rpm_start=None, addCornerRPM=True, nRPM=100, fact_start=0.5, fact_max=1.3):
        """ Return a dataframe with the data from the simple variable speed model from the File

        INPUTS:
         - addCornerRPM: when true, SpdGenOn and RtGenSp are added to the rpm vector

        """
        SpdGenOn = self.File['SpdGenOn']   # Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]
        RtGnSp   = self.VS_RtGnSp  # Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]
        RtTq     = self.VS_RtTq    # Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]
        Rgn2K    = self.VS_Rgn2K   # Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]
        SlPc     = self.VS_SlPc    # Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) [used only when VSContrl=1]
        
        # Default values
        if rpm_start is None:
            rpm_start = SpdGenOn
            if RtGnSp < rpm_start:
                rpm_start=0
        if rpm is None:
            rpm = np.linspace(rpm_start*fact_start, RtGnSp*fact_max, nRPM) # RPM, HSS


        # Add known RPMs values
        if addCornerRPM:
            rpm = np.unique(np.append(rpm, RtGnSp))           # Ensure RtGnSp is part of it
            if SpdGenOn > rpm[0] and SpdGenOn<rpm[-1]:
                rpm = np.unique(np.append(rpm, SpdGenOn))     # Ensure SpdGenOn is part of it

        GenTq    = gentorque(rpm, (RtGnSp, RtTq  , Rgn2K , SlPc , rpm_start))

        M = np.column_stack((rpm, GenTq))
        cols = ['Generator_Speed_[rpm]', 'Generator_Torque_[Nm]']
        df = pd.DataFrame(data=M, columns=cols)
        return df


    def VS_plot(self, rpm=None, rpm_start=None,
                ax=None, label='ServoDyn VS', ls='-', color='k' # Plot options
                ):

        df = self.VS_DataFrame(rpm=rpm, rpm_start=rpm_start)


        if ax is None:
            fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)

        ax.plot(df['Generator_Speed_[rpm]'], df['Generator_Torque_[Nm]']/1000, label=label, ls=ls, color=color)


        if self._data_rpm is not None:
            ax.plot(self._data_rpm, self._data_torque, 'o', label='Data', ms=3)
        if self._fit_rpm is not None:
            ax.plot(self._fit_rpm, self._fit_torque, '--', label='Fit')

        ax.legend()

        ax.set_xlabel('Generator speed [rpm]')
        ax.set_ylabel('Generator Torque [kNm]')
        return ax


    def VS_fit(self, rpm, torque, verbose=False):
        # --- Option 1a - pyDatView, Dedicated/auto, calls GeneratorTorqueFitter

        #>>> Model fit sFunc : fitter: gentorque
        #>>> Model fit p0    : None
        #>>> Model fit bounds: None
        #>>> Model fit kwargs: {}
        #y_fit, pfit, fitter = model_fit('fitter: gentorque', rpm, torque, bounds=None, p0=None, verbose=True)

        #>>> Model fit sFunc : predef: gentorque
        #>>> Model fit p0    : RtGnSp=100 , RtTq=1000  , Rgn2K=0.01 ,SlPc=5 , SpdGenOn=0
        #>>> Model fit bounds: RtGnSp=(0.1,inf) , RtTq=(1,inf), Rgn2K=(0.0,0.1) ,SlPc=(0,20) , SpdGenOn=(0,inf)
        #>>> Model fit kwargs: OrderedDict()

        # bounds='RtGnSp=(10.4,11.8) , RtTq=(4.9e6,5.2e6), Rgn2K=(10000.0,30000) ,SlPc=(0,15) , SpdGenOn=(1,6)'
        # p0 = [RtGnSp, RtTq, Rgn2K, SlPc, SpdGenOn]
        # y_fit, pfit, fitter = model_fit('predef: gentorque', rpm, torque) #, bounds=bounds, p0=p0)

        # --- Option 2
        p0     = None
        bounds = None
        #RtGnSp=7 , RtTq=32000  , Rgn2K=1000 ,SlPc=5 , SpdGenOn=2
        #RtGnSp=(0.1,inf) , RtTq=(1,inf), Rgn2K=(0.0,inf) ,SlPc=(0,20) , SpdGenOn=(0,inf)
        fitter = GeneratorTorqueFitter(rpm, torque, p0=p0, bounds=bounds)

        # Store data
        self._data_rpm = rpm.copy()
        self._data_torque = torque.copy()

        self._fitter = fitter
        self._fit_rpm = rpm
        self._fit_torque = fitter.model['fitted_function'](rpm)


        if verbose:
            coeffs = fitter.model['coeffs']
            def toStringVLD(lab,val,descr=''):
                val='{}'.format(val)
                lab='{}'.format(lab)
                # Trying to reproduce WISDEM format
                if len(val)<22:
                    val='{:22s}'.format(val)
                if len(lab)<11:
                    lab='{:11s}'.format(lab)
                return val+' '+lab+' - '+descr.strip().lstrip('-').lstrip()


            print('Fitted parameters:')
            print(toStringVLD('VS_RtGnSp', coeffs['RtGnSp']))
            print(toStringVLD('VS_RtTq  ', coeffs['RtTq']  ))
            print(toStringVLD('VS_Rgn2K ', coeffs['Rgn2K'] ))
            print(toStringVLD('VS_SlPc  ', coeffs['SlPc']  ))
            print('')
            print('Comparison with File:')
            print("{:22s}{:11s}".format(str(np.around(coeffs['RtGnSp'],6)), "RtGnSp") + "File:" + "{}".format(self.File['VS_RtGnSp'])) # TODO, RelError)
            print("{:22s}{:11s}".format(str(np.around(coeffs['RtTq'],6)),   "RtTq")   + "File:" + "{}".format(self.File['VS_RtTq']))
            print("{:22s}{:11s}".format(str(np.around(coeffs['Rgn2K'],6)),  "Rgn2K")  + "File:" + "{}".format(self.File['VS_Rgn2K']))
            print("{:22s}{:11s}".format(str(np.around(coeffs['SlPc'],6)),   "SlPc")   + "File:" + "{}".format(self.File['VS_SlPc']))




        return self._fit_torque, fitter

