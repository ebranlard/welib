"""
Tools for ServoDyn

- Plot the simple generator model
- Fit a simple generator model


"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from welib.weio.fast_input_file import FASTInputFile
# from pydatview.tools.curve_fitting import gentorque, GeneratorTorqueFitter, model_fit
from welib.tools.curve_fitting import gentorque, GeneratorTorqueFitter

RPM2RADS=   2*np.pi/60;
RADS2RPM=1/(2*np.pi/60);

class ServoDyn:

    def __init__(self, svdFilename_or_data=None, TP=None):
        """ 
        Initialize a ServoDyn object either with:
          - svdFilename: a servody input file name
          - svdData: an instance of FASTInputFile
        """

        # --- Data
        self.File=None
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


    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        s+='|properties:\n'
        s+='|- File: (input file data)\n'
#         s+='|- TP  : {} \n'.format(self._TP)
#         s+='|* graph: (Nodes/Elements/Members)\n'
#         s+='|* pointsMJ, pointsMN, pointsMNout\n'
        s+='|methods:\n'
        s+='|- VS_print\n'
        s+='|- VS_plot\n'
        s+='|- VS_dataframe\n'
#         s+='|- setTopMass\n'
#         s+='|- beamDataFrame, beamFEM, beamModes\n'
#         s+='|- toYAMSData\n'
        return s


    def VS_print(self):
        VSContrl = self.File['VSContrl']  # Variable-speed control mode {0: none, 1: simple VS, 3: user-defined from routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)
        GenModel = self.File['GenModel']  # Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]
        GenEff   = self.File['GenEff']    # Generator efficiency [ignored by the Thevenin and user-defined generator models] (%)
        GenTiStr = self.File['GenTiStr']  # Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)
        GenTiStp = self.File['GenTiStp']  # Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)
        SpdGenOn = self.File['SpdGenOn']  # Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]
        TimGenOn = self.File['TimGenOn']  # Time to turn on the generator for a startup (s) [used only when GenTiStr=True]
        TimGenOf = self.File['TimGenOf']  # Time to turn off the generator (s) [used only when GenTiStp=True]
        RtGnSp   = self.File['VS_RtGnSp'] # Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]
        RtTq     = self.File['VS_RtTq']   # Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]
        Rgn2K    = self.File['VS_Rgn2K']  # Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]
        SlPc     = self.File['VS_SlPc']   # Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) [used only when VSContrl=1]
        print('SpdGenOn : {:20.3f} [rpm]'      .format(SpdGenOn))
        print('VS_RtGnSp: {:20.3f} [rpm]'      .format(RtGnSp),', ',np.around(RtGnSp*RPM2RADS,4),'rad/s')
        print('VS_RtTq  : {:20.3f} [Nm]'       .format(RtTq))
        print('VS_Rgn2K : {:20.3f} [N-m/rpm^2]'.format(Rgn2K))
        print('VS_SlPc  : {:20.3f} [%]'        .format(SlPc))
        print('')
        print('VS_Rgn2K : {:20.3f} [N-m/(rad/s)^2]'.format(Rgn2K*RADS2RPM**2))
        print('')
        RtTqCheck = Rgn2K*RtGnSp**2
        print('TqCheck  : {:20.3f} [Nm]'.format(RtTqCheck))
        if RtTqCheck>RtTq:
            WARN(' Rgn2K*RtGnSp**2 ({})  > RtTq ({})'.format(RtTqCheck, RtTq))


    def VS_DataFrame(self, rpm=None, rpm_start=None, addCornerRPM=True, nRPM=100, fact_start=0.5, fact_max=1.3):
        """ Return a dataframe with the data from the simple variable speed model from the File

        INPUTS:
         - addCornerRPM: when true, SpdGenOn and RtGenSp are added to the rpm vector

        """
        SpdGenOn = self.File['SpdGenOn']   # Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]
        RtGnSp   = self.File['VS_RtGnSp']  # Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]
        RtTq     = self.File['VS_RtTq']    # Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]
        Rgn2K    = self.File['VS_Rgn2K']   # Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]
        SlPc     = self.File['VS_SlPc']    # Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) [used only when VSContrl=1]
        
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

        ax.plot(df['Generator_Speed_[rpm]'], df['Generator_Torque_[Nm]'], label=label, ls=ls, color=color)


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

