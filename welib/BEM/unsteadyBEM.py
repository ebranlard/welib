""" 
Python implementation of an unsteady BEM code

Reference:
   [1]:  E. Branlard (2017) "Wind Turbine Aerodynamics and Vorticity Based Method", Chapter 10, Springer

"""
import numpy as np
import os
from numpy import cos, sin, arctan2, pi, arccos, exp, sqrt
from numpy.linalg import norm
from scipy.interpolate import interp1d
import copy
import pandas as pd
import matplotlib.pyplot as plt
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid





# Load more models
# try:
from welib.BEM.highthrust import a_Ct
from welib.BEM.highthrust import a_k
from welib.yams.rotations import BodyXYZ_fromDCM, BodyXYZ_A, BodyXYZ_DCM # EulerExtract, EulerConstruct
# except: 
#     pass



class BEMStates:
    def __init__(self, nB, nr):
        pass
        # Induction
        # Dynamic wake
        # Dynamic stall

class BEMDiscreteStates:
    def __init__(self, nB, nr):
        self.t = None
        self.it=-1
        # Induction
        self.Vind_g  = np.zeros((nB,nr,3)) # Dynamic induced velocity with skew and dyn wake, global coordinates
        self.Vind_p  = np.zeros((nB,nr,3)) # Dynamic induced velocity with skew and dyn wake, polar coordinates
        self.a_qs    = np.zeros((nB,nr)) # axial induction
        # Dynamic wake
        self.Vind_qs_p  = np.zeros((nB,nr,3)) # Quasi-steady velocity, polar coordinates
        self.Vind_qs_g  = np.zeros((nB,nr,3)) # Quasi-steady velocity
        self.Vind_int_p = np.zeros((nB,nr,3)) # Intermediate velocity, polar coordinates
        self.Vind_dyn_p = np.zeros((nB,nr,3)) # Dynamic induced velocity (before skew/yaw), polar coordinates
        self.Vind_dyn_g = np.zeros((nB,nr,3)) # Dynamic induced velocity (before skew/yaw), global coordinates
        # Dynamic stall
        self.fs = np.zeros((nB,nr)) # Separation 

# --------------------------------------------------------------------------------}
# --- Main Class UnsteadyBEM 
# --------------------------------------------------------------------------------{
class UnsteadyBEM():
    """ 
    Perform unsteady BEM calculations
    """
    def __init__(self, filename=None):
        """ 
        filename: an OpenFAST (.fst) or AeroDyn driver (.dvr) input file.
                  will be used to initialize the parameters and operating conditions
        """
        # Aero Data
        self.chord  = None
        self.polars = None
        # Environment
        self.rho     = None
        self.kinVisc = None # Kinematic viscosity [kg/m^3]

        # Structural "Inputs"
        # position, velocity, and orientation of all blade station
        #self.cone   = None
        #self.twist  = None
        self.nB     = None
        self.r      = None # radial stations
        self.meanLineAC = None #(nB, nr, 3)) x,y,z position of the AC (prebend, presweep, "r")

        # Algorithm
        self.setDefaultOptions()

        # Init from given input file
        if filename is not None:
            self.init_from_FAST(filename) # Init from a FAST input file 

    def setDefaultOptions(self):
        #self.projMod = 'noSweepPitchTwist' # 
        #self.projMod = 'liftingLine' # 
        self.projMod = 'polar' # 
        self.algorithm = 'legacy' # Main switch for algorithm
        self.nIt = 200  # maximum number of iterations in BEM
        self.aTol = 10 ** -6 # tolerance for axial induction factor convergence
        self.relaxation = 0.5  # relaxation factor in axial induction factor
        self.CTcorrection = 'AeroDyn15'  #  type of CT correction more model implementated in the future like 'spera'
        self.swirlMethod  = 'AeroDyn' # type of swirl model
        self.bUseCm = True  # Use Moment 
        self.bSwirl = True  # swirl flow model enabled / disabled
        self.bTipLoss = True # enable / disable tip loss model
        self.bHubLoss = False # enable / disable hub loss model
        self.bTipLossCl = False # enable / disable Cl loss model
        self.TipLossMethod = 'Glauert'  # type of tip loss model
        self.bDynaStall = True # dynamic stall model
        self.bDynaWake = True # dynamic stall model
#           1   DBEMT_Mod          - Type of dynamic BEMT (DBEMT) model {1=constant tau1, 2=time-dependent tau1} (-) [used only when WakeMod=2]
#           4   tau1_const         - Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1]

        self.bYawModel = True # Yaw correction
        self.bAIDrag = True # influence on drag coefficient on normal force coefficient
        self.bTIDrag = True # influence on drag coefficient on tangential force coefficient
        self.WakeMod=1 # 0: no inductions, 1: BEM inductions

    def init_from_FAST(self, FASTFileName):
        from welib.weio.fast_input_deck import FASTInputDeck
        F = FASTInputDeck(FASTFileName,readlist=['AD','ED','ADbld','AF','IW'])
        driver =  F.version=='AD_driver'
        # --- Safety checkes
        if F.AD is None:
            raise Exception('steadyBEM: Cannot open AD file referenced in:'.format(FASTFileName))
        if driver:
            dvr = F.fst
            if dvr['NumTurbines']>1:
                raise NotImplementedError('steadyBEM: Number of turbines should be 1')
            if not dvr['BasicHAWTFormat(1)']:
                raise NotImplementedError('steadyBEM: BasicHAWTFormat should be true for now')
        if not driver:
            if F.ED is None:
                raise Exception('steadyBEM: Cannot open ED file referenced in:'.format(FASTFileName))

        # --- Environment
        try:
            self.rho     = float(F.fst['FldDens'])  # New OF > 3.0
        except:
            try:
                self.rho     = float(F.fst['AirDens'])  # New OF > 3.0
            except:
                self.rho     = float(F.AD['AirDens'])   # Old OF <=3.0
        try:
            self.kinVisc = float(F.fst['KinVisc'])  # New OF > 3.0
        except:
            self.kinVisc = float(F.AD['KinVisc'])   # Old OF <= 3.0

        # --- Geometry
        if driver:
            self.nB       = dvr['NumBlades(1)']
            r_hub         = dvr['HubRad(1)']
            # Input geometry (may be overriden by motion/simulations)
            self.cone0    = dvr['PreCone(1)']
            self.tilt0    = dvr['ShftTilt(1)']  # OpenFAST convention, negative about yn [deg]
            self.OverHang = dvr['OverHang(1)']
            self.Twr2Shft = dvr['Twr2Shft(1)']

            HubHt = dvr['HubHt(1)']
            self.TowerHt  = HubHt - self.Twr2Shft + self.OverHang*sin(self.tilt0*np.pi/180) # TODO double check

        else:
            self.nB    = F.ED['NumBl']
            r_hub      = F.ED['HubRad']
            # Input geometry (may be overriden by motion/simulations)
            self.cone0     = F.ED['PreCone(1)']
            self.tilt0     = F.ED['ShftTilt']  # OpenFAST convention, negative about yn [deg]
            self.TowerHt   = F.ED['TowerHt']
            self.OverHang  = F.ED['OverHang']
            self.Twr2Shft  = F.ED['Twr2Shft']


        # --- Aerodynamics
        self.r        = F.AD.Bld1['BldAeroNodes'][:,0] + r_hub
        self.chord    = np.stack([F.AD.Bld1['BldAeroNodes'][:,5]]*self.nB)
        nr = len(self.r)
        self.meanLineAC = np.zeros((self.nB, nr,3))
        for iB in range(self.nB):
            # TODO Use AD.Bld123
            self.meanLineAC[iB,:,0] = F.AD.Bld1['BldAeroNodes'][:,1] # along xb
            self.meanLineAC[iB,:,1] = F.AD.Bld1['BldAeroNodes'][:,2] # along yb
            self.meanLineAC[iB,:,2] = self.r
        self.twist    = F.AD.Bld1['BldAeroNodes'][:,4]*np.pi/180
        polars=[]
        ProfileID=F.AD.Bld1['BldAeroNodes'][:,-1].astype(int)
        for ipolar in  ProfileID:
            nTabs = F.AD.AF[ipolar-1]['NumTabs']
            if nTabs==1:
                polars.append(F.AD.AF[ipolar-1]['AFCoeff'])
            else:
                print('[WARN] unsteadyBEM multiple polar present')
                vRe   = [F.AD.AF[ipolar-1]['re_{}'.format(i+1)] for i in range(nTabs)]
                vCtrl = [F.AD.AF[ipolar-1]['Ctrl_{}'.format(i+1)] for i in range(nTabs)]
                print(vRe)
                print(vCtrl)
                # Taking last polar...
                polars.append(F.AD.AF[ipolar-1]['AFCoeff_{}'.format(nTabs)])
        self.polars = polars
        self.bAIDrag  = F.AD['AIDrag']
        self.bTIDrag  = F.AD['TIDrag']
        self.bHubLoss = F.AD['HubLoss']
        self.bTipLoss = F.AD['TipLoss']
        self.bSwirl   = F.AD['TanInd']
        self.bUseCm   = F.AD['UseBlCm']

        # --- Operating conditions
        # Maybe...

	# trigger once data is set
        self._init()
     
     
    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)
        # Aero Data
        s+='Aerodynamic data:\n'
        s+=' - r         : size: {}, min: {}, max: {}\n'.format(len(self.r), np.min(self.r), np.max(self.r))
        s+=' - chord     : size: {}, min: {}, max: {}\n'.format(len(self.chord), np.min(self.chord), np.max(self.chord))
        s+=' - twist     : size: {}, min: {}, max: {} [deg]\n'.format(len(self.twist), np.min(self.twist), np.max(self.twist))
        s+=' - polars    : size: {}\n'.format(len(self.polars))
        s+=' - nB        : {}\n'.format(self.nB)
        s+=' - cone0     : {} [deg] (OpenFAST convention, about yh)\n'.format(self.cone0 )
        s+=' - tilt0     : {} [deg] (OpenFAST convention, negative along yn)\n'.format(self.tilt0 )
        # Environment
        s+='Environmental conditions:\n'
        s+=' - rho       : {} [kg/m^3]\n'.format(self.rho    )
        s+=' - kinVisc   : {} [kg/m^3]\n'.format(self.kinVisc)
        # Latest operating conditions
        #s+='Latest operating conditions:\n'
        #s+=' - Omega0    : {} [rpm]\n'.format(self.Omega0)
        #s+=' - pitch0    : {} [deg]\n'.format(self.pitch0)
        #s+=' - V0        : {} [m/s]\n'.format(self.V0)
        s+='Algorithm options:\n'
        s+=' - projMod   : {}\n'.format(self.projMod )
        s+=' - algorithm : {}\n'.format(self.algorithm )
        s+=' - nIt       : {}\n'.format(self.nIt      )
        s+=' - aTol      : {}\n'.format(self.aTol      )
        s+=' - relaxation: {}\n'.format(self.relaxation)
        s+=' - bSwirl    : {}\n'.format(self.bSwirl    )
        s+=' - bTipLoss  : {}\n'.format(self.bTipLoss  )
        s+=' - bHubLoss  : {}\n'.format(self.bHubLoss  )
        s+=' - bAIDrag   : {}\n'.format(self.bAIDrag   )
        s+=' - bTIDrag   : {}\n'.format(self.bTIDrag  )
        s+=' - WakeMod   : {}\n'.format(self.WakeMod  )
        return s

    def _init(self):
        # Creating interpolation functions for each polar, now in rad!
        # TODO RE/Ctrl
        self.fPolars = [interp1d(p[:,0]*np.pi/180,p[:,1:],axis=0) for p in self.polars]

    def getInitStates(self):
        return BEMDiscreteStates(self.nB, len(self.r))

    def timeStepInit(self, t0, tmax, dt):
        """ Allocate storage for tiem step values"""
        self.time=np.arange(t0,tmax+dt/2,dt)
        nt = len(self.time)
        nB = self.nB
        nr = len(self.r)
        # --- Spanwise data
        # Coeffients
        self.TSR    = np.zeros(nt)
        self.Cl_qs  = np.zeros((nt,nB,nr))
        self.Cd_qs  = np.zeros((nt,nB,nr))
        self.Cl     = np.zeros((nt,nB,nr))
        self.Cd     = np.zeros((nt,nB,nr))
        self.Cm     = np.zeros((nt,nB,nr))
        self.cn     = np.zeros((nt,nB,nr))
        self.ct     = np.zeros((nt,nB,nr))
        self.Cx_a   = np.zeros((nt,nB,nr))
        self.Cy_a   = np.zeros((nt,nB,nr))
        self.Ct     = np.zeros((nt,nB,nr))
        self.Cq     = np.zeros((nt,nB,nr))
        # Velocities
        self.Vrel_p = np.zeros((nt,nB,nr,3)) # Un,Ut,Ur
        self.Vrel_xa = np.zeros((nt,nB,nr)) # 
        self.Vrel_ya = np.zeros((nt,nB,nr)) # 
        self.Vrel_za = np.zeros((nt,nB,nr)) # 
        self.Vind_g = np.zeros((nt,nB,nr,3))
        self.Vind_h = np.zeros((nt,nB,nr,3))
        self.Vind_p = np.zeros((nt,nB,nr,3))
        self.Vind_s = np.zeros((nt,nB,nr,3))
        self.Vind_qs_p = np.zeros((nt,nB,nr,3))
        self.Vflw_p = np.zeros((nt,nB,nr,3)) # Vwnd-Vstr
        self.Vwnd_g = np.zeros((nt,nB,nr,3))
        self.Vwnd_h = np.zeros((nt,nB,nr,3))
        self.Vwnd_p = np.zeros((nt,nB,nr,3))
        self.Vwnd_s = np.zeros((nt,nB,nr,3))
        self.Vwnd_a = np.zeros((nt,nB,nr,3))
        self.Vstr_g = np.zeros((nt,nB,nr,3))
        self.Vstr_h = np.zeros((nt,nB,nr,3))
        self.Vstr_p = np.zeros((nt,nB,nr,3))
        self.Vstr_s = np.zeros((nt,nB,nr,3))
        self.Vstr_xa = np.zeros((nt,nB,nr))
        self.Vstr_ya = np.zeros((nt,nB,nr))
        self.Vrel   = np.zeros((nt,nB,nr))
        # BEM var
        self.AxInd  = np.zeros((nt,nB,nr))
        self.TnInd  = np.zeros((nt,nB,nr))
        self.AxInd_qs = np.zeros((nt,nB,nr))
        self.TnInd_qs = np.zeros((nt,nB,nr))
        self.BEM_F    = np.zeros((nt,nB,nr))
        self.BEM_k    = np.zeros((nt,nB,nr))
        self.BEM_kp   = np.zeros((nt,nB,nr))
        self.BEM_CT_qs= np.zeros((nt,nB,nr))
        # Loads per span
        self.L      = np.zeros((nt,nB,nr))
        self.D      = np.zeros((nt,nB,nr))
        self.Mm     = np.zeros((nt,nB,nr))
        self.Fn     = np.zeros((nt,nB,nr))
        self.Ft     = np.zeros((nt,nB,nr))
        self.F_a   = np.zeros((nt,nB,nr,3)) # Force in a
        self.F_i   = np.zeros((nt,nB,nr,3)) # Force  in inertial/global
        self.M_i   = np.zeros((nt,nB,nr,3)) # Moment in inertial/global
        self.F_s    = np.zeros((nt,nB,nr,3))
        self.Gamma  = np.zeros((nt,nB,nr))
        self.alpha  = np.zeros((nt,nB,nr))
        self.phi    = np.zeros((nt,nB,nr))
        self.Re     = np.zeros((nt,nB,nr))
        # AeroDyn outputs in ill-defined / algo-dependent systems
        self.AD_F_o = np.zeros((nt,nB,nr,3)) # Outputs of AeroDyn, algo specific
        self.AD_M_o = np.zeros((nt,nB,nr,3)) # Outputs of AeroDyn, algo specific
        self.AD_Fn = np.zeros((nt,nB,nr))
        self.AD_Ft = np.zeros((nt,nB,nr))
        # Integrated values
        self.BAeroF_i  = np.zeros((nt,nB,3))
        self.BAeroM_i  = np.zeros((nt,nB,3))
        self.RtAeroF_i = np.zeros((nt,3))
        self.RtAeroM_i = np.zeros((nt,3))
        self.RtAeroF_r = np.zeros((nt,3))
        self.RtAeroM_r = np.zeros((nt,3))
        self.Thrust   = np.zeros(nt)
        self.Torque   = np.zeros(nt)
        self.Power    = np.zeros(nt)
        self.chi      = np.zeros(nt)
        self.chi0     = np.zeros(nt)
        self.RtVAvg   = np.zeros((nt,3))
        self.psi      = np.zeros(nt)
        self.Omega    = np.zeros(nt)  # [rpm]
        self.RtArea   = np.zeros(nt)
        self.SkewAzimuth  = np.zeros((nt,nB))
        # Blade blades
        self.BladeTorque = np.zeros((nt,nB))
        self.BladeThrust = np.zeros((nt,nB))
        self.BladeEdge   = np.zeros((nt,nB))
        self.BladeFlap   = np.zeros((nt,nB))



    def calcOutput(self):
        """ 
        Compute output for current state variables and inputs.
        """
        R = np.sqrt(self.RtArea/pi)
        q = 0.5*self.rho*self.RtArea*self.RtVAvg[:,0]**2
        self.CT = self.Thrust/(q)
        self.CQ = self.Torque/(q*R)
        self.CP = self.Power /(q*self.RtVAvg[:,0])


    # --------------------------------------------------------------------------------}
    # --- IO 
    # --------------------------------------------------------------------------------{
    def toDataFrame(self, BldNd_BladesOut=None, BldNd_BlOutNd=None, BldNd_OutList=None):
        """ Export time series to a pandas dataframe
        Column names are set to match OpenFAST outputs

        BldNd_BladesOut: if None, all blades nodal output is output. Otherwise set it to an integer < self.nB
        BldNd_BlOutNd:   if None, all radial position are output.    Otherwise set it to a list of radial indices [0,2,4, etc] 
        BldNd_OutList:   if None, all possible channels are output.  Otherwise, set it to a list of channel names with units. 
                         for instance:  ['Fx_[N/m]', 'Vx_[m/s]', 'AxInd_[-]', Phi_[deg]'] 
        """
        p=self # Alias


        columns=['Time_[s]']
        columns+=['Thrust_[N]']
        columns+=['Torque_[N/m]']

        #if not hasattr(self,'CP'):
        self.calcOutput()

        df = pd.DataFrame()
        df['Time_[s]']        = self.time
        df['Azimuth_[deg]']   = np.mod(self.psi,360)
        df['RtSpeed_[rpm]']   = self.Omega
        df['RtAeroFxh_[N]']   = self.RtAeroF_r[:,0]
        df['RtAeroFyh_[N]']   = self.RtAeroF_r[:,1]
        df['RtAeroFzh_[N]']   = self.RtAeroF_r[:,2]
        df['RtAeroMxh_[N-m]'] = self.RtAeroM_r[:,0]
        df['RtAeroMyh_[N-m]'] = self.RtAeroM_r[:,1]
        df['RtAeroMzh_[N-m]'] = self.RtAeroM_r[:,2]
        df['RtAeroPwr_[W]']   = self.Power
        df['RtAeroCt_[-]']    = self.CT
        df['RtAeroCq_[-]']    = self.CQ
        df['RtAeroCp_[-]']    = self.CP
        df['RtTSR_[-]']        = self.TSR
        df['RtVAvgxh_[m/s]']  = self.RtVAvg[:,0]
        df['RtVAvgyh_[m/s]']  = self.RtVAvg[:,1]
        df['RtVAvgzh_[m/s]']  = self.RtVAvg[:,2]
        df['RtArea_[m^2]']    = self.RtArea
        df['RtSkew_[deg]']    = self.chi0

        # TODO These should be renamed to "i"
        df['B1AeroFxg_[N]']   = self.BAeroF_i[:,0,0]
        df['B1AeroFyg_[N]']   = self.BAeroF_i[:,0,1]
        df['B1AeroFzg_[N]']   = self.BAeroF_i[:,0,2]
        df['B1AeroMxg_[N-m]'] = self.BAeroM_i[:,0,0]
        df['B1AeroMyg_[N-m]'] = self.BAeroM_i[:,0,1]
        df['B1AeroMzg_[N-m]'] = self.BAeroM_i[:,0,2]
        df['RtAeroFxg_[N]']   = self.RtAeroF_i[:,0]
        df['RtAeroFyg_[N]']   = self.RtAeroF_i[:,1]
        df['RtAeroFzg_[N]']   = self.RtAeroF_i[:,2]
        df['RtAeroMxg_[N-m]'] = self.RtAeroM_i[:,0]
        df['RtAeroMyg_[N-m]'] = self.RtAeroM_i[:,1]
        df['RtAeroMzg_[N-m]'] = self.RtAeroM_i[:,2]

        for iB in np.arange(self.nB):
            df['B'+str(iB+1)+'Azimuth_[deg]']  = np.mod(self.psi+self.SkewAzimuth[:,iB],360)

        if self.projMod=='polar': # TODO replace with "outProj"
            Vflw_o = self.Vwnd_p-self.Vstr_p
            Vwnd_o = self.Vwnd_p.copy()
        elif self.projMod=='noSweepPitchTwist':
            Vflw_o = self.Vwnd_s-self.Vstr_s
            Vwnd_o = self.Vwnd_s.copy()
        else:
            raise NotImplementedError()

        # --- Creating all the possible column names for all blades and radial position, matching OpenFAST convention
        if BldNd_BladesOut is None: 
            BldNd_BladesOut=np.arange(self.nB)
        else:
            BldNd_BladesOut=np.arange(min(BldNd_BladesOut, self.nB))
        if BldNd_BlOutNd is None: 
            BldNd_BlOutNd = np.arange(len(self.r))
        if BldNd_OutList is None:

            BldNd_OutList =[]
            BldNd_OutList+=['VDisx', 'VDisy' ,'VDisz']
            BldNd_OutList+=['VDisxi','VDisyi','VDiszi']
            BldNd_OutList+=['VDisxp','VDisyp','VDiszp']
            BldNd_OutList+=['VDisxh','VDisyh','VDiszh']
            BldNd_OutList+=['STVx','STVy','STVz']
            BldNd_OutList+=['STVxi' ,'STVyi', 'STVzi']
            BldNd_OutList+=['STVxp' ,'STVyp', 'STVzp']
            BldNd_OutList+=['STVxh' ,'STVyh', 'STVzh']
            BldNd_OutList+=['Vx','Vy']
            BldNd_OutList+=['Vindx','Vindy']
            BldNd_OutList+=['Vindxi','Vindyi','Vindzi']
            BldNd_OutList+=['Vindxp','Vindyp','Vindzp']
            BldNd_OutList+=['Vindxh','Vindyh','Vindzh']
            BldNd_OutList+=['AxInd_qs','TnInd_qs']
            BldNd_OutList+=['BEM_F','BEM_k','BEM_kp','BEM_CT_qs']
            BldNd_OutList+=['AxInd'   ,'TnInd'  ]
            BldNd_OutList+=['Vrel','Phi' ,'Alpha']
            BldNd_OutList+=['Cl','Cd']
            BldNd_OutList+=['Re', 'Gam']
            BldNd_OutList+=['Fx','Fy']
            BldNd_OutList+=['Fn','Ft','Mm']
            BldNd_OutList+=['Fxi','Fyi','Fzi']
            BldNd_OutList+=['Mxi','Myi','Mzi']

        RefColUnits={
            'fx'        : 'Fx'       +'_[N/m]'   ,
            'fy'        : 'Fy'       +'_[N/m]'   ,
            'fxi'       : 'Fxi'      +'_[N/m]'   ,
            'fyi'       : 'Fyi'      +'_[N/m]'   ,
            'fzi'       : 'Fzi'      +'_[N/m]'   ,
            'mxi'       : 'Mxi'      +'_[N-m/m]'   ,
            'myi'       : 'Myi'      +'_[N-m/m]'   ,
            'mzi'       : 'Mzi'      +'_[N-m/m]'   ,
            'vx'        : 'Vx'       +'_[m/s]'   ,
            'vy'        : 'Vy'       +'_[m/s]'   ,
            'vdisx'     : 'VDisx'    +'_[m/s]'   ,
            'vdisy'     : 'VDisy'    +'_[m/s]'   ,
            'vdisz'     : 'VDisz'    +'_[m/s]'   ,
            'vdisxi'    : 'VDisxi'   +'_[m/s]'   ,
            'vdisyi'    : 'VDisyi'   +'_[m/s]'   ,
            'vdiszi'    : 'VDiszi'   +'_[m/s]'   ,
            'vdisxh'    : 'VDisxh'   +'_[m/s]'   ,
            'vdisyh'    : 'VDisyh'   +'_[m/s]'   ,
            'vdiszh'    : 'VDiszh'   +'_[m/s]'   ,
            'vdisxp'    : 'VDisxp'   +'_[m/s]'   ,
            'vdisyp'    : 'VDisyp'   +'_[m/s]'   ,
            'vdiszp'    : 'VDiszp'   +'_[m/s]'   ,
            'stvx'      : 'STVx'     +'_[m/s]'   ,
            'stvy'      : 'STVy'     +'_[m/s]'   ,
            'stvz'      : 'STVz'     +'_[m/s]'   ,
            'stvxi'     : 'STVxi'    +'_[m/s]'   ,
            'stvyi'     : 'STVyi'    +'_[m/s]'   ,
            'stvzi'     : 'STVzi'    +'_[m/s]'   ,
            'stvxh'     : 'STVxh'    +'_[m/s]'   ,
            'stvyh'     : 'STVyh'    +'_[m/s]'   ,
            'stvzh'     : 'STVzh'    +'_[m/s]'   ,
            'stvxp'     : 'STVxp'    +'_[m/s]'   ,
            'stvyp'     : 'STVyp'    +'_[m/s]'   ,
            'stvzp'     : 'STVzp'    +'_[m/s]'   ,
            'vrel'      : 'Vrel'     +'_[m/s]'   ,
            'axind'     : 'AxInd'    +'_[-]'     ,
            'tnind'     : 'TnInd'    +'_[-]'     ,
            'axind_qs'  : 'AxInd_qs' +'_[-]'     ,
            'tnind_qs'  : 'TnInd_qs' +'_[-]'     ,
            'phi'       : 'Phi'      +'_[deg]'   ,
            'vindx'     : 'Vindx'    +'_[m/s]'   ,
            'vindy'     : 'Vindy'    +'_[m/s]'   ,
            'vindxi'    : 'Vindxi'   +'_[m/s]'   ,
            'vindyi'    : 'Vindyi'   +'_[m/s]'   ,
            'vindzi'    : 'Vindzi'   +'_[m/s]'   ,
            'vindxh'    : 'Vindxh'   +'_[m/s]'   ,
            'vindyh'    : 'Vindyh'   +'_[m/s]'   ,
            'vindzh'    : 'Vindzh'   +'_[m/s]'   ,
            'vindxp'    : 'Vindxp'   +'_[m/s]'   ,
            'vindyp'    : 'Vindyp'   +'_[m/s]'   ,
            'vindzp'    : 'Vindzp'   +'_[m/s]'   ,
            'alpha'     : 'Alpha'    +'_[deg]'   ,
            'fn'        : 'Fn'       +'_[N/m]'   ,
            'ft'        : 'Ft'       +'_[N/m]'   ,
            'mm'        : 'Mm'       +'_[N-m/m]' ,
            'cl'        : 'Cl'       +'_[-]'     ,
            'cd'        : 'Cd'       +'_[-]'     ,
            're'        : 'Re'       +'_[-]'     ,
            'gam'       : 'Gam'      +'_[m^2/s]' ,
            'bem_f'     : 'BEM_F'    +'_[-]'     ,
            'bem_k'     : 'BEM_k'    +'_[-]'     ,
            'bem_kp'    : 'BEM_kp'   +'_[-]'     ,
            'bem_ct_qs' : 'BEM_CT_qs'+'_[-]'     }

        # All columns
        BldNd_columns = ['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+ RefColUnits[c.lower()] for iB in BldNd_BladesOut  for c in BldNd_OutList for ir in BldNd_BlOutNd]
        # Dataframe
        df_B_r = pd.DataFrame(np.zeros((len(self.time), len(BldNd_columns))), columns=BldNd_columns)

        df= pd.concat((df,df_B_r),axis=1)

        def set3dVar(ll,  var):
            lab = RefColUnits[ll]
            #print('>>> Setting {:15s}  - min:{:13.3e} max:{:13.3e}'.format(lab, np.min(var.flatten()), np.max(var.flatten())))
            for iB in BldNd_BladesOut:
                for ir in np.arange(len(self.r)): 
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+lab] = var[:,iB,ir]


        # AeroDyn x-y is "section coord" s
        # AeroDyn n-t is "airfoil coord" a but y is switched
        # AeroDyn doesn't have polar coord..
        # TODO make it case/unit insensitive and based on AeroDyn OutList!
        for l in BldNd_OutList:
            ll = l.lower()
            if ll == 'vx'       : set3dVar(ll, Vflw_o        [: ,:, :,0] )
            elif ll == 'vy'       : set3dVar(ll, Vflw_o        [: ,:, :,1] )
            elif ll == 'vdisx'    : set3dVar(ll, Vwnd_o        [: ,:, :,0] )
            elif ll == 'vdisy'    : set3dVar(ll, Vwnd_o        [: ,:, :,1] )
            elif ll == 'vdisz'    : set3dVar(ll, Vwnd_o        [: ,:, :,2] )
            elif ll == 'vdisxi'   : set3dVar(ll, self.Vwnd_g   [: ,:, :,0] )
            elif ll == 'vdisyi'   : set3dVar(ll, self.Vwnd_g   [: ,:, :,1] )
            elif ll == 'vdiszi'   : set3dVar(ll, self.Vwnd_g   [: ,:, :,2] )
            elif ll == 'vdisxh'   : set3dVar(ll, self.Vwnd_h   [: ,:, :,0] )
            elif ll == 'vdisyh'   : set3dVar(ll, self.Vwnd_h   [: ,:, :,1] )
            elif ll == 'vdiszh'   : set3dVar(ll, self.Vwnd_h   [: ,:, :,2] )
            elif ll == 'vdisxp'   : set3dVar(ll, self.Vwnd_p   [: ,:, :,0] )
            elif ll == 'vdisyp'   : set3dVar(ll, self.Vwnd_p   [: ,:, :,1] )
            elif ll == 'vdiszp'   : set3dVar(ll, self.Vwnd_p   [: ,:, :,2] )
            elif ll == 'stvx'     : set3dVar(ll, self.Vstr_s   [: ,:, :,0] )
            elif ll == 'stvy'     : set3dVar(ll, self.Vstr_s   [: ,:, :,1] )
            elif ll == 'stvz'     : set3dVar(ll, self.Vstr_s   [: ,:, :,2] )
            elif ll == 'stvxi'    : set3dVar(ll, self.Vstr_g   [: ,:, :,0] )
            elif ll == 'stvyi'    : set3dVar(ll, self.Vstr_g   [: ,:, :,1] )
            elif ll == 'stvzi'    : set3dVar(ll, self.Vstr_g   [: ,:, :,2] )
            elif ll == 'stvxh'    : set3dVar(ll, self.Vstr_h   [: ,:, :,0] )
            elif ll == 'stvyh'    : set3dVar(ll, self.Vstr_h   [: ,:, :,1] )
            elif ll == 'stvzh'    : set3dVar(ll, self.Vstr_h   [: ,:, :,2] )
            elif ll == 'stvxp'    : set3dVar(ll, self.Vstr_p   [: ,:, :,0] )
            elif ll == 'stvyp'    : set3dVar(ll, self.Vstr_p   [: ,:, :,1] )
            elif ll == 'stvzp'    : set3dVar(ll, self.Vstr_p   [: ,:, :,2] )
            elif ll == 'vrel'     : set3dVar(ll, self.Vrel     [: ,:, :]   )
            elif ll == 'axind'    : set3dVar(ll, self.AxInd    [: ,:, :]   )
            elif ll == 'tnind'    : set3dVar(ll, self.TnInd    [: ,:, :]   )
            elif ll == 'axind_qs' : set3dVar(ll, self.AxInd_qs [: ,:, :]   )
            elif ll == 'tnind_qs' : set3dVar(ll, self.TnInd_qs [: ,:, :]   )
            elif ll == 'phi'      : set3dVar(ll, self.phi      [: ,:, :]   )
            elif ll == 'vindx'    : set3dVar(ll, self.Vind_s   [: ,:, :,0] )
            elif ll == 'vindy'    : set3dVar(ll, self.Vind_s   [: ,:, :,1] )
            elif ll == 'vindxi'   : set3dVar(ll, self.Vind_g   [: ,:, :,0] )
            elif ll == 'vindyi'   : set3dVar(ll, self.Vind_g   [: ,:, :,1] )
            elif ll == 'vindzi'   : set3dVar(ll, self.Vind_g   [: ,:, :,2] )
            elif ll == 'vindxh'   : set3dVar(ll, self.Vind_h   [: ,:, :,0] )
            elif ll == 'vindyh'   : set3dVar(ll, self.Vind_h   [: ,:, :,1] )
            elif ll == 'vindzh'   : set3dVar(ll, self.Vind_h   [: ,:, :,2] )
            elif ll == 'vindxp'   : set3dVar(ll, self.Vind_p   [: ,:, :,0] )
            elif ll == 'vindyp'   : set3dVar(ll, self.Vind_p   [: ,:, :,1] )
            elif ll == 'vindzp'   : set3dVar(ll, self.Vind_p   [: ,:, :,2] )
            elif ll == 'alpha'    : set3dVar(ll, self.alpha    [: ,:, :]   )
            elif ll == 'fxi'      : set3dVar(ll, self.F_i      [: ,:, :,0] )
            elif ll == 'fyi'      : set3dVar(ll, self.F_i      [: ,:, :,1] )
            elif ll == 'fzi'      : set3dVar(ll, self.F_i      [: ,:, :,2] )
            elif ll == 'mxi'      : set3dVar(ll, self.M_i      [: ,:, :,0] )
            elif ll == 'myi'      : set3dVar(ll, self.M_i      [: ,:, :,1] )
            elif ll == 'mzi'      : set3dVar(ll, self.M_i      [: ,:, :,2] )
            elif ll == 'mm'       : set3dVar(ll,self.Mm       [: ,:, :]    ) 
            elif ll == 'cl'       : set3dVar(ll, self.Cl       [: ,:, :]   )
            elif ll == 'cd'       : set3dVar(ll, self.Cd       [: ,:, :]   )
            elif ll == 're'       : set3dVar(ll, self.Re       [: ,:, :]   )
            elif ll == 'gam'      : set3dVar(ll, self.Gamma    [: ,:, :]   )
            elif ll == 'bem_f'    : set3dVar(ll, self.BEM_F    [: ,:, :]   )
            elif ll == 'bem_k'    : set3dVar(ll, self.BEM_k    [: ,:, :]   )
            elif ll == 'bem_kp'   : set3dVar(ll, self.BEM_kp   [: ,:, :]   )
            elif ll == 'bem_ct_qs': set3dVar(ll, self.BEM_CT_qs[: ,:, :]   )
            # --- AeroDyn algo-dependent outputs
            #elif ll == 'fn'       : set3dVar(ll, self.F_s      [: ,:, :,0] ) # TODO should be F_a
            #elif ll == 'ft'       : set3dVar(ll,-self.F_s      [: ,:, :,1] ) # TODO
#             elif ll == 'fx'       : set3dVar(ll, self.F_s      [: ,:, :,0] ) # TODO should be F_a
#             elif ll == 'fy'       : set3dVar(ll,-self.F_s      [: ,:, :,1] ) # TODO
            elif ll == 'fn'       : set3dVar(ll, self.AD_Fn       [: ,: ,:]   )    #  m%X(IdxNode,IdxBlade)*ct - m%Y(IdxNode,IdxBlade)*st
            elif ll == 'ft'       : set3dVar(ll, self.AD_Ft       [: ,: ,:]   )    # -m%X(IdxNode,IdxBlade)*st - m%Y(IdxNode,IdxBlade)*ct
            elif ll == 'fx'       : set3dVar(ll, self.AD_F_o      [: ,:, :,0] ) #  m%X(IdxNode,IdxBlade)
            elif ll == 'fy'       : set3dVar(ll,-self.AD_F_o      [: ,:, :,1] ) # -m%Y(IdxNode,IdxBlade)
            else:
                raise Exception('Unavaible channel',l)
        return df




    def toDataFrameRadial(self, it=-1):
        df = pd.DataFrame()
        df['r/R_[-]']   = self.r/max(self.r)
        #df['r_[m]']     = self.r
        for iB in np.arange(self.nB):
            B='B'+str(iB+1)
            df[B+'AxInd_[-]'] = self.AxInd[it,iB,:]
            df[B+'TnInd_[-]'] = self.AxInd[it,iB,:]
        return df


    def timeStep(self, t, dt, xd0, psi, psiB0,
            origin_pos_gl, omega_gl, R_r2g,  # Kinematics of rotor origin
            R_bld2SB, # from Blade to Shaft-Blade for each blade
            pos_gl, Vstr_gl, R_a2g,            # Kinematics of nodes
            Vwnd_gl, # Wind at each positions in global
            firstCallEquilibrium=False,
            kappa = None # cant angle, TODO TODO TODO get rid of me, temporary. User EulerExtract of R_ntr2g
            ):
        """ 
        xBEM0: BEM states at t-1

        INPUTS:
         - t: current time step [s]
         - dt: time interval [s]
         - xd0: discreate states at current time step. Instance of BEMDiscreteStates
         - psi: current azimuth [rad]
         - psiB0: azimuthal offsets for each blades compared to psi 
         - origin_pos_gl: position of rotor origin in global coordinates
         - omega_gl:  rotational speed of rotor in global coordinates
         - R_r2g   :  transformation from rotor coordinates to global
         - R_bld2SB :  transformation from blade to Shaft-Blade (for each blade)        (nB x 3 x 3)
                      typically consist of azimuth, cone and pitch
         - pos_gl: positions of all blade nodes in global              (nB x nr x 3)
         - Vstr_gl: structural velocity of a ll blade nodes in global   (nB x nr x 3)
         - Vwnd_gl: wind velocity of a ll blade nodes in global        (nB x nr x 3)
         - R_a2g  : transformation matrix from airfoil to global       (nB x nr x 3 x 3)
        """
        xd1   = copy.deepcopy(xd0)
        xd1.t = t
        xd1.it = xd0.it+1 # Increase time step 
        # Safety
        if xd0.t is not None:
            if np.abs((xd0.t+dt-t))>dt/10:
                raise Exception('timeStep method expects to be called on regular intervals')

        # Time step storage for vectorization
        nB, nr, _ = pos_gl.shape
        p = self # alias
        ### Loop on blades
        # --------------------------------------------------------------------------------
        # --- Step 0a: geometry 
        # --------------------------------------------------------------------------------
        x_hat_disk = R_r2g[:,0]

        # --- "n-t-r" and rotating blade hub system
        # "n-t-r"  is the same but with t=-yh
        # Precalculate the M_ph matrix -- no reason to recalculate for each output
        R_r2h = np.zeros((3,3))  # Rotor to blade rotating hub
        R_g2h = np.zeros((self.nB,3,3))  # Rotor to blade rotating hub
        R_ntr2g = np.zeros((self.nB,3,3)) # "n-t-r" to global
        for  iB in range(self.nB):
             psi_hub = 2*pi*iB/self.nB
             R_r2h[0,:] = [ 1, 0           , 0            ]
             R_r2h[1,:] = [ 0, cos(psi_hub), sin(psi_hub) ]
             R_r2h[2,:] = [ 0,-sin(psi_hub), cos(psi_hub) ]
             R_g2h[iB, :, :] = R_r2h.dot( R_r2g.T ) # inertial to Blade rotating hub
             # "n-t-r"  is the same but with t=-yh
             RR=(R_g2h[iB].T).copy()
             RR[:,1]*=-1
             R_ntr2g[iB] = RR

        # --- Polar grid coordinates and radial-coordinate 
        R_g2p = np.zeros((self.nB, nr, 3, 3)) # transformation from global to staggered-polar, "orientationAnnulus"
        r_p   = np.zeros((self.nB, nr))       # radius in polar grid, "rLocal"
        drdz  = np.zeros((self.nB, nr))       # projection of element dz
        r_center2nodes = np.zeros((nB, nr, 3)) # Lever arm vector from rotor center to blade nodes
        Rs_p=[0]*nB
        P_rotor_i = origin_pos_gl
        x_hat_disk = R_r2g[:,0] 
        for iB in range(self.nB):
            r_RA      = np.zeros((nr,3))  # Vector from rotor center to Airfoil
            r_RA_orth = np.zeros((nr,3))  # Vector from rotor center to airfoil projected on polar plane
            for ie in range(nr):
                R_g2p[iB,ie,:,:], r_RA[ie,:], r_RA_orth[ie,:] = polarCoord(P_rotor_i, x_hat_disk, pos_gl[iB,ie])  
                r_center2nodes[iB, ie, : ] = r_RA[ie,:] # Used for moment calculation
                r_p[iB, ie] = norm(r_RA_orth[ie,:])  # "r_local"
            dr = norm(r_RA_orth[1:]-r_RA_orth[0:-1], axis=1)
            dz = norm(r_RA     [1:]-r_RA     [0:-1], axis=1)
            dr = np.concatenate(([dr[0]], dr))
            dz = np.concatenate(([dz[0]], dz))
            drdz[iB, dz!=0] = dr/dz
        rhub_p  = np.mean(r_p[iB,0])
        R_p     = np.max (r_p[iB,-1]) # Rotor radius projected onto polar grid
        if self.algorithm=='polarProj':
            pass
        else:
            drdz = 1
        R_s2g = np.zeros((nB, nr, 3, 3))
        if self.algorithm=='polarProj':
            # NOTE: section is ill-defined
            # AeroDyn uses "orientaion annulus everywhere..
            for iB in range(self.nB):
                for ie in range(nr):
                    R_s2g[iB,ie,:,:] = R_g2p[iB,ie,:,:].T
        else:
            pass
            #raise Exception()

        # --- Curvilinear coordinates
        # See Init_BEMTmodule, "zRoot, zTip, zLocal"
        # NOTE: in OpenFAST this does is done at Init only (no aeroelastic deformation)
        sHub = np.zeros(self.nB)              # "zRoot"
        sBldNodes = np.zeros((self.nB, nr))   # "zLocal"
        for iB in range(self.nB):
            P_root_i = pos_gl[iB, 0] # NOTE: that's an approximation
            sHub[iB], sBldNodes[iB] = curvilinearCoord(P_rotor_i, P_root_i, pos_gl[iB,:])

        # --- From polar to airfoil (includes aeroelastic deformation)
        R_p2a = np.zeros((self.nB, nr, 3, 3)) 
        toe   = np.zeros((self.nB, nr))
        cant  = np.zeros((self.nB, nr))
        beta   = np.zeros((self.nB, nr)) # including pitch and aeroelastic torsion
        for iB in range(self.nB):
            for ie in range(nr):
                R_p2a[iB,ie] = (R_a2g[iB,ie].T).dot(R_g2p[iB,ie].T)
                phis =  BodyXYZ_fromDCM(R_p2a[iB,ie]) # "EulerExtract"
                toe [iB, ie] =  phis[0]
                cant[iB, ie] =  phis[1]
                beta  [iB, ie]= -phis[2] # NOTE: - for wind energy convention
        #iB=0
        #print('')
        #print('toe' ,toe[iB]*180/np.pi)
        #print('cant',cant[iB]*180/np.pi)
        #print('r_p',r_p[iB])
        #print('rHub_p',rhub_p)
        #print('R_p',   R_p)
        #print('>>> sHub',sHub[iB])
        #print('>>> sBldNodes',sBldNodes[iB])
        # --- TipLoss HubLoss - Constants
        tipLossConst = np.zeros((self.nB, nr))
        hubLossConst = np.zeros((self.nB, nr))
        for iB in range(self.nB):
            tipLossConst[iB] = self.nB*(sBldNodes[iB,-1] - sBldNodes[iB,:])/(2.0*sBldNodes[iB,:])
            # NOTE different conventions are possible for hub losses
            hubLossConst[iB] = self.nB*(sBldNodes[iB, :] - sHub[iB]       )/(2.0*sHub[iB])

        bFixedInduction = np.zeros((self.nB, nr), dtype=bool)
        if p.bTipLoss:
            bFixedInduction = np.logical_or(bFixedInduction, tipLossConst==0)
        if p.bHubLoss:
            bFixedInduction = np.logical_or(bFixedInduction, hubLossConst==0)
        tipLossConst[bFixedInduction] = 1
        hubLossConst[bFixedInduction] = 1
        bFValid = np.logical_not(bFixedInduction)



        # --- Rotor speed for power
        omega_r = R_r2g.T.dot(omega_gl) # rotational speed in rotor coordinate system
        Omega = omega_r[0] # rotation speed of shaft (along x)


        # --------------------------------------------------------------------------------}
        # --- Step 0b: Wind and structural Velocities - Skew coordinate system
        # --------------------------------------------------------------------------------{
        # --- Flow
        Vstr_p  = np.zeros((nB,nr,3))
        Vwnd_p  = np.zeros((nB,nr,3))
        for iB in np.arange(nB):
            for ie in np.arange(nr):
                R_p2g = R_g2p[iB,ie].T
                # Polar coordinates
                Vstr_p[iB,ie] = (R_p2g.T).dot(Vstr_gl[iB,ie]) # Structural velocity in polar coordinates
                Vwnd_p[iB,ie] = (R_p2g.T).dot(Vwnd_gl[iB,ie]) # Wind Velocity in polar coordinates
        Vflw_p  = Vwnd_p- Vstr_p  # Relative flow velocity, including wind and structural motion
        Vflw_g  = Vwnd_gl-Vstr_gl # Relative flow velocity, including wind and structural motion
        #print('Vwndxp,',Vwnd_p[0,:,0])
        #print('Vwndyp,',Vwnd_p[0,:,1])
        #print('Vwndzp,',Vwnd_p[0,:,2])
        #print('Vstrxp,',Vstr_p[0,:,0])
        #print('Vstryp,',Vstr_p[0,:,1])
        #print('Vstrzp,',Vstr_p[0,:,2])
        #print('Vflw_xp,',Vflw_p[0,:,0])
        #print('Vflw_yp,',Vflw_p[0,:,1])
        #print('Vflw_zp,',Vflw_p[0,:,2])

        # --- Average wind
        # Average wind in global, and rotor coord
        Vwnd_avg_g = np.mean(np.mean(Vwnd_gl[:,:,:],axis=0),axis=0)
        Vwnd_avg_r = (R_r2g.T).dot(Vwnd_avg_g)
        # Average relative wind (Wnd-Str)
        Vflw_avg_g = np.mean(np.mean(Vflw_g[:,:,:],axis=0),axis=0) # "V_diskAvg"
        # Coordinate system with "y" in the cross wind direction for skew model
        V_dot_x  = np.dot(Vflw_avg_g, x_hat_disk)
        V_ytmp   = V_dot_x * x_hat_disk - Vflw_avg_g
        V_ynorm  = sqrt(V_ytmp[0]**2+V_ytmp[1]**2+V_ytmp[2]**2)
        if abs(V_ynorm)<1e-8:
            y_hat_disk = R_r2g[:,1]
            z_hat_disk = R_r2g[:,2]
        else:
            y_hat_disk = V_ytmp / V_ynorm
            z_hat_disk = np.cross(Vflw_avg_g, x_hat_disk ) / V_ynorm
        # --- Skew
        # Fake "Azimuth angle" used for skew model
        SkewAzimuth=np.zeros(nB)
        for iB in np.arange(nB):
            z_hat = R_ntr2g[iB][:,2] # TODO TODO
            tmp_sz_y = -1.0*np.dot(z_hat,y_hat_disk)
            tmp_sz   =      np.dot(z_hat,z_hat_disk)
            if np.abs(tmp_sz_y)<1e-8 and np.abs(tmp_sz)<1e-8:
                SkewAzimuth[iB]=0
            else:
                SkewAzimuth[iB] = arctan2( tmp_sz_y, tmp_sz )
        # Skew angle without induction
        Vw_r = (R_r2g.T).dot(Vflw_avg_g)
        Vw_rn     = Vw_r[0] # normal to disk
        Vw_r_norm = sqrt(Vw_r[0]**2+Vw_r[1]**2+Vw_r[2]**2)
        chi0      = np.arccos(Vw_rn / Vw_r_norm)

        # --- TSR
        if V_dot_x==0:
            TSR = 0
        else:
            TSR =  Omega * R_p / V_dot_x


        # --- START ITERATION FOR QUASI-STEADY INDUCTIONS
        # NOTE: Typically no iteration is done. 
        #       We assume that the induced velocity is known (in xd0.Vind_p, xd0.Vind_g)
        #       and from that we compute a non converged quasi-steady induction.
        # NOTE: to be rigorous we would need to use inputs at t and t+dt
        if firstCallEquilibrium:
            nItMax=50
        else:
            nItMax=1
        for iteration in np.arange(nItMax):
            #print('')
            #print('>>>>>>>>>>>> TIME ', t, '>>>>>>>>>>>>>>>>>>>>> ITERATION ', iteration)
            # --------------------------------------------------------------------------------
            # --- Step 1: velocity components based on induced velocities in xd0
            # --------------------------------------------------------------------------------
            Vrel_a  = np.zeros((nB,nr,3))
            Vrel_p  = np.zeros((nB,nr,3))
            Vrel_k  = np.zeros((nB,nr,3))
            Vind_g  = np.zeros((nB,nr,3))
            for iB in np.arange(nB):
                for ie in np.arange(nr):
                    #R_p2g = R_ntr2g[iB] 
                    #R_p2g = R_g2h[iB].T 
                    R_p2g = R_g2p[iB,ie].T
                    # NOTE: inductions from previous time step, in polar grid (more realistic than global)
                    # Global
                    #Vind_g[iB,ie] = xd0.Vind_g[iB,ie] # dynamic inductions at previous time step
                    Vind_g[iB,ie] = (R_p2g).dot(xd0.Vind_p[iB,ie]) # dynamic inductions at previous time step
                    Vrel_g = Vwnd_gl[iB,ie]+Vind_g[iB,ie]-Vstr_gl[iB,ie]
                    # Polar coordinates
                    Vrel_p[iB,ie] = (R_p2g.T).dot(Vrel_g)
                    # Kappa coordinates
                    Vrel_k[iB,ie,0] = Vrel_p[iB,ie,0]*np.cos(cant[iB,ie]) # n 
                    Vrel_k[iB,ie,1] = Vrel_p[iB,ie,1] # t
                    # Airfoil coordinates
                    Vrel_a[iB,ie] = (R_a2g[iB,ie].T).dot(Vrel_g) # TODO use R_p2a instead, and remove zp component
            #print('xd0Vind_xp,',xd0.Vind_p[0,:,0])
            #print('xd0Vind_yp,',xd0.Vind_p[0,:,1])
            #print('xd0Vind_zp,',xd0.Vind_p[0,:,2])

            # Velocity norm and Reynolds
            Vrel_norm_k = sqrt(Vrel_k[:,:,0]**2 + Vrel_k[:,:,1]**2)
            Vrel_norm_a = sqrt(Vrel_a[:,:,0]**2 + Vrel_a[:,:,1]**2)
            # TODO Change AeroDyn for Re and Vrel calculations
            a_dyn      =-xd0.Vind_p[:,:,0]/Vflw_p[:,:,0] 
            aprime_dyn = xd0.Vind_p[:,:,1]/Vflw_p[:,:,1]
            VrelUA, v_ac_x, v_ac_y = RelativeVelocityUA( a_dyn, aprime_dyn, Vflw_p[:,:,0], Vflw_p[:,:,1], cant, xVelCorr=0)
            Re        = VrelUA*p.chord/p.kinVisc/10**6 # Reynolds in million
            #Re        = Vrel_norm_a*p.chord/p.kinVisc/10**6 # Reynolds in million
            #ReUA = ReynoldsNumberUA(iteration, self.algorithm, a_dyn, aprime_dyn, Vflw_p[:,:,0], Vflw_p[:,:,1], Vflw_p[:,:,2], p.chord, p.kinVisc, toe, cant, twist)

            # --------------------------------------------------------------------------------
            # --- Step 2: Flow Angle and tip loss
            # --------------------------------------------------------------------------------
            phi_k = np.arctan2(Vrel_k[:,:,0], Vrel_k[:,:,1]) # flow angle [rad] in kappa system
            phi_p = np.arctan2(Vrel_p[:,:,0], Vrel_p[:,:,1]) # NOTE: using polar grid for phi
            if self.algorithm=='legacy':
                phi_tl = phi_p
                phi    = phi_p
            elif self.algorithm=='polarProj':
                #phi_tl = phi_kappa(phi_k, cant)
                phi_tl = phi_p                    # TODO Check equivalence with above 
                phi    = phi_k
            else:
                raise Exception()
            # --- Tip and hub losses
            F    = np.ones((nB,nr))
            Ftip = np.ones((nB,nr))
            Fhub = np.ones((nB,nr))
            bFValid = np.logical_and(bFValid, phi_tl>0)
            sphi   = sin(phi_tl[bFValid])
            if (p.bTipLoss): #Glauert tip correction
                # --- Keep me
                #b=sin(phi_tl)>0.01
                #Ftip[b] = 2./pi*arccos(exp(-(nB *(R_p-r_p[b]))/(2*r_p[b]*sin(phi_tl[b]))))
                #b2=abs(r_p-R_p)<1e-3
                #Ftip[b2]=0.001
                factorTip     = np.clip(tipLossConst[bFValid] / sphi , -1    , np.inf)
                expFactor     = np.clip(np.exp(-factorTip)           , 0     ,  1)
                Ftip[bFValid] = np.clip(2./pi*arccos(expFactor)      , 0.0001,  1)

            # --- Hub loss
            if (p.bHubLoss): #Glauert hub loss correction
                # TODO TODO TODO
                factorHub     = np.clip(hubLossConst[bFValid] / sphi , -1    , np.inf)
                expFactor     = np.clip(np.exp(-factorHub)           , 0     ,  1)
                Fhub[bFValid] = np.clip(2./pi*arccos(expFactor)      , 0.0001,  1)
                #Fhub[bFValid] = 2./pi*arccos(exp(-nB/2. *(r_p-rhub_p)/ (rhub_p*np.sin(phi_tl))))
            F = Ftip * Fhub

            # --------------------------------------------------------------------------------
            # --- Step 3: Angle of attack
            # --------------------------------------------------------------------------------
            alpha = np.arctan2(Vrel_a[:,:,0],Vrel_a[:,:,1])        # angle of attack [rad]
            #print('ALPHA\n', alpha*180/np.pi)
            # --------------------------------------------------------------------------------
            # --- Step 4: Aerodynamic Coefficients
            # --------------------------------------------------------------------------------
            ClCdCm = np.array([p.fPolars[ie](alpha[iB,ie]) for iB in np.arange(nB) for ie in np.arange(nr)]).reshape((nB,nr,3))
            Cl=ClCdCm[:,:,0]
            Cd=ClCdCm[:,:,1]
            Cm=ClCdCm[:,:,2]
            # Project to airfoil coordinates
            C_xa,        C_ya        =  ProjectClCd_to_a(Cl, Cd, alpha)
            C_xa_noDrag, C_ya_noDrag =  ProjectClCd_to_a(Cl, Cd*0, alpha)
            # Project to polar coordinates
            C_p        = np.zeros((nB,nr,3))
            C_g        = np.zeros((nB,nr,3))
            C_p_noDrag = np.zeros((nB,nr,3))
            for iB in np.arange(nB):
                for ie in np.arange(nr):
                    R_p2g = R_g2p[iB,ie].T
                    C_g        [iB,ie]=R_a2g[iB,ie].dot(np.array([C_xa[iB,ie], C_ya[iB,ie], 0]))
                    C_p        [iB,ie]=(R_p2g.T).dot(C_g[iB,ie])
                    C_p_noDrag [iB,ie]=(R_p2g.T).dot(R_a2g[iB,ie]).dot(np.array([C_xa_noDrag[iB,ie], C_ya_noDrag[iB,ie], 0]))
            # Project elementary radial element ds vs dr
            # 
            # Cn and Ct 
            if (p.bAIDrag):
                cnForAI = C_p[:,:,0]
            else:
                cnForAI = C_p_noDrag[:,:,0]
            if (p.bTIDrag):
                ctForTI = C_p[:,:,1]
            else:
                ctForTI = C_p_noDrag[:,:,1]
            # L = 0.5 * p.rho * Vrel_norm**2 * p.chord[ie]*Cl
            # --------------------------------------------------------------------------------
            # --- Step 5: Quasi-steady induction
            # --------------------------------------------------------------------------------
            # NOTE: all is done in polar grid
            #lambda_r = Vstr_p[:,:,1]/Vwnd_p[:,:,0] # "omega r/ U0n" defined in polar grid # TODO TODO TODO
            lambda_r = -Vflw_p[:,:,1]/Vflw_p[:,:,0] # "omega r/ U0n" defined in polar grid # TODO TODO TODO
            #lambda_r = Vflw_p[:,:,1]/Vflw_p[:,:,0] # "omega r/ U0n" defined in polar grid # TODO TODO TODO
            #V0       = np.sqrt(Vwnd_p[:,:,0]**2 + Vwnd_p[:,:,1]**2) # TODO think about that # TODO TODO TOD
            V0       = np.sqrt(Vflw_p[:,:,0]**2 + Vwnd_p[:,:,1]**2) # TODO think about that
            sigma    = p.chord*p.nB/(2*pi*r_p) # NOTE: using radius in polar grid
            #a,aprime,CT = fInductionCoefficients(a_last,Vrel_in4,Un,Ut,V0_in3,V0_in4,nnW_in4,omega,chord(e),F,Ftip,CnForAI,CtForTI,lambda_r,sigma(e),phi,Algo)
            if p.WakeMod==0:
                a_qs   = V0*0
                aprime_qs = V0*0
            else:
                if p.algorithm =='legacy':
                    a_qs,aprime_qs,CT = _fInductionCoefficients(Vrel_norm_k, V0, F, cnForAI, ctForTI, lambda_r, sigma, phi, 
                            bSwirl=p.bSwirl, CTcorrection=p.CTcorrection, swirlMethod=p.swirlMethod,
                            relaxation=p.relaxation, a_last=xd0.a_qs,
                            algorithm=self.algorithm, drdz=drdz
                    )

                    #k =  sigma*cnForAI/(4*F*np.sin(phi)**2)*drdz # NOTE: OpenFAST formulation
                    #kp =-sigma*ctForTI/(4*F*np.sin(phi)*np.cos(phi))
                    k  = 1/(1-a_qs)
                    kp = aprime_qs/(1+aprime_qs)
                else:

                    Vrel_xp = clip_zeros(Vrel_p[:,:,0], 1e-5) # To avoid division by zero
                    Vrel_yp = clip_zeros(Vrel_p[:,:,1], 1e-5)

                    # --- k factor
                    k =  sigma*cnForAI/(4*F*np.sin(phi)**2)*drdz # NOTE: OpenFAST formulation
                    k2 = sigma*cnForAI/(4*F*np.sin(phi)**2)*cos(cant)**2/drdz   
                    k1 = sigma*cnForAI/(4*F)*Vrel_norm_a**2/(Vrel_xp**2)/drdz
                    # TODO Solve for a_qs
                    # TODO Solve for a_qs
                    # TODO Solve for a_qs
                    if abs(chi0)>0:
                        a_qs = V0*0
                        for iB in np.arange(nB):
                            for ie in np.arange(nr):
                                if bFValid[iB,ie]:
                                    a_qs[iB,ie] = a_k(k[iB, ie],  chi=chi0, phi=phi[iB, ie], F=F[iB,ie], method='MomentumGlauertSkewRootsWithHT')
                    else:
                        a_qs = k/(k+1.0)
                    effectiveYaw = chi0 # TODO TODO
                    #print('chi0',chi0)
                    #print('>>xd0a',xd0.a_qs[0])
                    #print('>>xd1a',xd1.a_qs[0])
                    a_qs_last = xd0.a_qs # TODO Temporary because right now a is solved above
                    a0_an = np.sqrt(1+(np.tan(effectiveYaw)/(1.0-a_qs_last))**2)
                    a0_an = 1
                    kpCorrectionFactor = 1
                    kp =-sigma*ctForTI/(4*F*np.sin(phi)*np.cos(phi))*(kpCorrectionFactor)/a0_an            
                    kp1=-sigma*ctForTI/(4*F)*Vrel_norm_a**2/(Vrel_xp*Vrel_yp) /drdz
#                     if iteration==3:
#                         import pdb; pdb.set_trace()

                    # TODO Solve for aprime_qs
                    # TODO Solve for aprime_qs
                    # TODO Solve for aprime_qs
                    aprime_qs = kp/(1-kp)
                    CT=0 # TODO TODO TODO

            # --- Wherever the tiploss factor is undefined 
            k        [bFixedInduction] = 0
            kp       [bFixedInduction] = 0
            a_qs     [bFixedInduction] = 1.00
            aprime_qs[bFixedInduction] = 0.00
            # --- HACK
            #a_qs[:,:]      = 0.3
            #aprime_qs[:,:] = 0.05
            #a_qs     [bFixedInduction] = 0.20
            #aprime_qs[bFixedInduction] = 0.002

            if np.any(np.isnan(a_qs)):
                print('>> BEM crashing')

            # Storing last values, for relaxation
            xd1.a_qs = a_qs.copy()
            # Quasi steady inductions, polar and global coordinates
            xd1.Vind_qs_p = np.zeros((nB,nr,3))
            for iB in np.arange(nB):
                for ie in np.arange(nr):
                    R_p2g = R_g2p[iB,ie].T
                    xd1.Vind_qs_p[iB,ie] = np.array([-a_qs[iB,ie]*Vflw_p[iB,ie,0],  aprime_qs[iB,ie]*Vflw_p[iB,ie,1], 0])
                    xd1.Vind_qs_g[iB,ie] = R_p2g.dot(xd1.Vind_qs_p[iB,ie]) # global

            if firstCallEquilibrium:
                # We update the previous states induction
                xd0.a_qs   = a_qs.copy()
                xd0.Vind_g = xd1.Vind_qs_g.copy()
                xd0.Vind_p = xd1.Vind_qs_p.copy()

            # --- END OF ITERATION LOOP

        # --- ITERATION LOOP OVER


        if firstCallEquilibrium:
            # Initialize dynamic wake variables
            xd0.Vind_qs_p  = xd1.Vind_qs_p.copy()
            xd1.Vind_qs_p  = xd1.Vind_qs_p.copy()
            xd0.Vind_int_p = xd1.Vind_qs_p.copy()
            xd0.Vind_dyn_p = xd1.Vind_qs_p.copy()
        # --------------------------------------------------------------------------------
        # --- Dynamic wake model, in polar coordinates (for "constant" structural velocity)
        # --------------------------------------------------------------------------------
        if (p.bDynaWake):
            #raise Exception()
            a_avg = min([np.mean(a_qs),0.5])
            V_avg = max([np.mean(V0),0.001])
            tau1 = 1.1 / (1 - 1.3 *a_avg)*R_p/V_avg
            tau2 = (0.39 - 0.26 * (r_p/R_p)**2) * tau1
            tau2 = np.tile(tau2[:,:,None],3)
            # Oye's dynamic inflow model, discrete time integration
            H              = xd1.Vind_qs_p + 0.6 * tau1 * (xd1.Vind_qs_p - xd0.Vind_qs_p) /dt
            xd1.Vind_int_p = H + (xd0.Vind_int_p - H) * exp(-dt/tau1) # intermediate velocity
            xd1.Vind_dyn_p = xd1.Vind_int_p + (xd0.Vind_dyn_p - xd1.Vind_int_p) * exp(-dt/tau2)
            # In global
            for iB in np.arange(nB):
                for ie in np.arange(nr):
                    #R_p2g = R_ntr2g[iB] # TODO TODO TODO
                    R_p2g = R_g2p[iB,ie].T
                    xd1.Vind_dyn_g[iB,ie] = R_p2g.dot(xd1.Vind_dyn_p[iB,ie]) # global
        else:
            xd1.Vind_dyn_g = xd1.Vind_qs_g.copy()
            xd1.Vind_dyn_p = xd1.Vind_qs_p.copy()


        # --------------------------------------------------------------------------------
        # ---  Yaw model, repartition of the induced velocity
        # --------------------------------------------------------------------------------
        if p.bYawModel:
            #raise Exception('YawModel')
            #xd1.Vind_g = xd1.Vind_dyn_g.copy()
            #xd1.Vind_p = xd1.Vind_dyn_p.copy()
            #psi0 = np.arctan( Vwnd_avg_g[2]/Vwnd_avg_r[1])  # TODO
            # Sections that are about 0.7%R
            Ir= np.logical_and(r_p[0]>=0.5*R_p, r_p[0] <=0.8*R_p)
            if len(Ir)==0:
                Ir=r_p[0]>0
            Vind_avg_g = np.mean(np.mean(xd1.Vind_dyn_g[:,Ir,:],axis=0),axis=0)
            Vind_avg_r = (R_r2g.T).dot(Vind_avg_g)
            # Skew angle with induction
            V_r      = Vwnd_avg_r + Vind_avg_r
            V_rn     = V_r[0] # normal to disk
            V_r_norm = sqrt(V_r[0]**2+V_r[1]**2+V_r[2]**2)
            chi = np.arccos(V_rn/V_r_norm)
            #print('chi0',chi0*180/pi,'chi',chi*180/pi,'psi0',psi0*180/np.pi)
            if np.abs(chi)>pi/2:
                print('>>> chi too large')
            yawCorrFactor = 15*np.pi/32 # close to 3/2
            for iB in np.arange(nB):
                for ie in np.arange(nr):
                    #R_p2g = R_ntr2g[iB]
                     R_p2g = R_g2p[iB,ie].T
                     xd1.Vind_p[iB,ie]   = xd1.Vind_dyn_p[iB,ie]
                     xd1.Vind_p[iB,ie,0] = xd1.Vind_dyn_p[iB,ie,0] * (1 + yawCorrFactor*r_p[iB,ie]/R_p * np.tan(chi/2)*np.sin(SkewAzimuth[iB])) #* np.cos(psiB0[iB]+psi - psi0))
                     xd1.Vind_g[iB,ie]   = R_p2g.dot(xd1.Vind_p[iB,ie]) # global
                     # AeroDyn:
                    #chi = (0.6_ReKi*a + 1.0_ReKi)*chi0
                    #a = a * (1.0 +  yawCorrFactor * yawCorr_tan * (tipRatio) * sin(azimuth))
        else:
           xd1.Vind_g = xd1.Vind_dyn_g.copy()
           xd1.Vind_p = xd1.Vind_dyn_p.copy()
        # Now the induced velocities are known


        # --------------------------------------------------------------------------------
        # --- Unstead airfoil aerodynamics/ dynamics stall
        # --------------------------------------------------------------------------------
        def setInputsForUA():
            # INPUTS for UA
            # TODO Change AeroDyn for Re and Vrel calculations
            a_dyn      =-xd1.Vind_p[:,:,0]/Vflw_p[:,:,0]  # NOTE: choose xd1 or xd0
            aprime_dyn = xd1.Vind_p[:,:,1]/Vflw_p[:,:,1]  # NOTE: choose xd1 or xd0
            VrelUA, v_ac_x, v_ac_y = RelativeVelocityUA( a_dyn, aprime_dyn, Vflw_p[:,:,0], Vflw_p[:,:,1], cant, xVelCorr=0)
            Re        = VrelUA*p.chord/p.kinVisc/10**6 # Reynolds in million
            #Re        = Vrel_norm_a*p.chord/p.kinVisc/10**6 # Reynolds in million
            #ReUA = ReynoldsNumberUA(iteration, self.algorithm, a_dyn, aprime_dyn, Vflw_p[:,:,0], Vflw_p[:,:,1], Vflw_p[:,:,2], p.chord, p.kinVisc, toe, cant, twist)
            return  VrelUA, v_ac_x, v_ac_y
        VrelUA, v_ac_x, v_ac_y = setInputsForUA()

        # Would likely occur here
        Cl_dyn = Cl.copy()
        Cd_dyn = Cd.copy()
        Cm_dyn = Cm.copy()
        C_xa_dyn, C_ya_dyn =  ProjectClCd_to_a(Cl, Cd, alpha)

        # --------------------------------------------------------------------------------
        # --- Step 6: Output Calculations 
        # --------------------------------------------------------------------------------
        it = xd1.it # time step

        # --- Dynamic induced velocity
        Vind_g     = xd1.Vind_g
        a_dyn      =-xd1.Vind_p[:,:,0]/Vflw_p[:,:,0] 
        aprime_dyn = xd1.Vind_p[:,:,1]/Vflw_p[:,:,1]

        # --- Loads
        # NOTE: Using VrelUA and not Vrel_norm_a for compatibility with AeroDyn
        q_dyn_c = 0.5 * p.rho * VrelUA**2 * p.chord # dynamic pressure with chord
        self.L[it]    = q_dyn_c * Cl_dyn
        self.D[it]    = q_dyn_c * Cd_dyn
        self.Mm[it]   = q_dyn_c * Cm_dyn * p.chord  # z component of "momentAirfoil", stored in "m%M"
        self.F_a[it,:,:,0] = q_dyn_c * C_xa_dyn     # "forceAirfoil"
        self.F_a[it,:,:,1] = q_dyn_c * C_ya_dyn
#         self.Fn[it]   = q_dyn * C_p[:,:,0]  # NOTE: this is what would make sense but not what AeroDyn does
#         self.Ft[it]   = q_dyn * C_p[:,:,1]
        #for  iB in range(self.nB):
        #    for ie in np.arange(nr):
        #        F_g = q_dyn[iB,ie] * C_g[iB,ie]
        #        self.F_s[it,iB,ie] = (R_s2g[iB,ie].T).dot(F_g)

        # --- Loads the AeroDyn way
        # NOTE NOTE NOTE: THIS IS AeroDyn and could use some rethinking
        # It has different meaning depending on the algorithm
        if p.algorithm=='legacy':
            Cx, Cy = Transform_ClCd_to_CxCy(Cl_dyn, Cd_dyn, phi)
            Cz  = 0
            Cmx = Cm_dyn*0
            Cmy = Cm_dyn*0
            Cmz = Cm_dyn
        else:
            Cx, Cy, Cz, Cmx, Cmy, Cmz = Transform_ClCdCm_to_CxCyCzCmxCmyCmz(Cl_dyn, Cd_dyn, Cm_dyn, toe, cant, beta, alpha)
        self.AD_F_o[it,:,:,0]  =  Cx       * q_dyn_c           #  "m%X"  For polarProj, this is Cxp
        self.AD_F_o[it,:,:,1]  = -Cy       * q_dyn_c           #  "m%Y"  For polarProj, this is Cyp
        self.AD_F_o[it,:,:,2]  =  Cz       * q_dyn_c           #  "m%Z"  For polarProj, this is Czp
        self.AD_M_o[it,:,:,0]  =  Cmx      * q_dyn_c * p.chord #  "m%MX" 
        self.AD_M_o[it,:,:,1]  =  Cmy      * q_dyn_c * p.chord #  "m%MY" 
        self.AD_M_o[it,:,:,2]  =  Cmz      * q_dyn_c * p.chord #  "m%MZ"  
        ct = cos(beta)
        st = sin(beta)
        self.AD_Fn[it,:,:] =  self.AD_F_o[it,:,:,0]*ct - self.AD_F_o[it,:,:,1]*st   #  m%X(IdxNode,IdxBlade)*ct - m%Y(IdxNode,IdxBlade)*st
        self.AD_Ft[it,:,:] = -self.AD_F_o[it,:,:,0]*st - self.AD_F_o[it,:,:,1]*ct   # -m%X(IdxNode,IdxBlade)*st - m%Y(IdxNode,IdxBlade)*ct
          
          



        # --------------------------------------------------------------------------------
        # --- Step 7: Storage and projection
        # --------------------------------------------------------------------------------
        # --- Velocities
        self.AxInd[it]    = a_dyn
        self.TnInd[it]    = aprime_dyn
        self.AxInd_qs[it] = a_qs
        self.TnInd_qs[it] = aprime_qs
        #self.Vrel[it]    = Vrel_norm_a
        self.Vrel[it]     = VrelUA # TODO Change AeroDyn
        self.Re[it]       = Re     # TODO Change AeroDyn
        # Global system
        self.Vwnd_g[it] = Vwnd_gl
        self.Vstr_g[it] = Vstr_gl
        self.Vind_g[it] = Vind_g
        # H system
        for  iB in range(self.nB):
            for ie in np.arange(nr):
                self.Vwnd_h[it,iB,ie] = R_g2h[iB].dot( Vwnd_gl[iB,ie,:])
                self.Vstr_h[it,iB,ie] = R_g2h[iB].dot( Vstr_gl[iB,ie,:])
                self.Vind_h[it,iB,ie] = R_g2h[iB].dot( Vind_g [iB,ie,:])
        for iB in np.arange(nB):
            for ie in np.arange(nr):
                R_p2g = R_g2p[iB,ie].T
                Vind_s = (R_s2g[iB,ie].T).dot(Vind_g[iB,ie]) # Induced velocity in section coordinates
                #Vind_a = (R_a2g[iB,ie].T).dot(Vind_g) # Induced velocity in airfoil coordinates
                Vind_p = (R_p2g       .T).dot(Vind_g[iB,ie]) # Induced velocity in polar coordinates
                self.Vind_p[it,iB,ie,:] =Vind_p
                self.Vind_s[it,iB,ie,:] =Vind_s
                ## Wind
                Vwnd_s = (R_s2g[iB,ie].T).dot(Vwnd_gl[iB,ie]) # Wind Velocity in section coordinates
                #Vwnd_a = (R_a2g[iB,ie].T).dot(Vwnd_gl) # Wind Velocity in airfoil coordinates
                self.Vwnd_s[it,iB,ie,:]=Vwnd_s
                ## Structural velocity
                #Vstr_a = (R_a2g[iB,ie].T).dot(Vstr_g) # Structural velocity in airfoil coordinates
                self.Vstr_s[it,iB,ie,:]=(R_s2g[iB,ie].T).dot(Vstr_gl[iB,ie])

        # --- Coefficients
        self.TSR[it] = TSR
        self.Cl[it]   = Cl
        self.Cd[it]   = Cd
        self.Cm[it]   = Cm
        self.cn[it]   = C_p[:,:,0]
        self.ct[it]   = C_p[:,:,1]
        # C_g also available


        # polar system (missing Vind)
        self.Vrel_p[it]  = Vrel_p[:,:,:] # NOTE: Vrel is using previous inductions..
        self.Vstr_p[it]  = Vstr_p[:,:,:]
        self.Vwnd_p[it]  = Vwnd_p[:,:,:]
        self.Vflw_p[it]  = Vflw_p[:,:,:]
        self.Vind_qs_p[it] = xd1.Vind_qs_p[:,:,:]
        self.RtVAvg[it]  = (R_r2g.T).dot(Vflw_avg_g) # in Hub/rotor coordinate
        self.SkewAzimuth[it]  = SkewAzimuth*180/pi
        self.chi0[it]  = chi0*180/pi
        # airfoil system
        self.Vrel_xa[it] = Vrel_a[:,:,0]
        self.Vrel_ya[it] = Vrel_a[:,:,1]
        self.Vrel_za[it] = Vrel_a[:,:,2]
        # --- BEM Variables 
        self.BEM_F[it]  = F
        self.BEM_k[it]  = k
        self.BEM_kp[it] = kp
        self.BEM_CT_qs[it] = 4*k*F*(1-a_qs)**2

        #self.alpha[it] = alpha*180./pi
        self.alpha[it] = (phi - beta)*180/np.pi # TODO: Change AeroDyn - y%WriteOutput( OutIdx )  = Rad2M180to180Deg( m%BEMT_y%phi(IdxNode,IdxBlade) - m%BEMT_u(Indx)%theta(IdxNode,IdxBlade) )
        self.phi[it]   = phi*180./pi
        # --- Misc
        self.Gamma[it]  = 0.5*Re*Cl*p.kinVisc*10**6 # Circulation [m^2/s]
        # Gamma y%WriteOutput( OutIdx ) = 0.5 * p%BEMT%chord(IdxNode,IdxBlade) * m%BEMT_y%Vrel(IdxNode,IdxBlade) * m%BEMT_y%Cl(IdxNode,IdxBlade) ! "Gam" [m^2/s]

        self.psi[it]  = psi*180/pi
        self.Omega[it]  = Omega*60/(2*np.pi) # [rpm]
        self.RtArea[it]  = pi*R_p**2

        # --- Loads
        if p.algorithm=='legacy':
            self.F_i[it,:,:,:] = 0
            self.M_i[it,:,:,:] = 0
            #y%BladeLoad(k)%Force(:,j)  = matmul( force,  m%orientationAnnulus(:,:,j,k) )  ! force per unit length of the jth node in the kth blade
            #y%BladeLoad(k)%Moment(:,j) = matmul( moment, m%orientationAnnulus(:,:,j,k) )  ! moment per unit length of the jth node in the kth blade
        else:
            for iB in np.arange(nB):
                for ie in np.arange(nr):
                    self.F_i[it,iB,ie] = R_a2g[iB,ie].dot(self.F_a[it,iB,ie]          )    #   y%BladeLoad(k)%Force(:,j) 
                    self.M_i[it,iB,ie] = R_a2g[iB,ie].dot( [0, 0, self.Mm[it,iB,ie]]  )    #   y%BladeLoad(k)%Moment(:,j)


        # --- Integral quantities for rotor
        F_i = self.F_i[it,:nB,:,:3] # For all blades
        M_i = self.M_i[it,:nB,:,:3] # For all blades
        M_rxF_i = np.zeros((nB, nr, 3)) # Force lever arm contribution to moments
        M_rxF_i = np.cross(r_center2nodes, F_i) # Force lever arm contribution to moments
        if self.bUseCm:
            # Important contribution in cone/prebend
            # QMz = trapezoid(self.Mm[it,:]*np.sin(cant), r_p)
            M_i += M_rxF_i 
        else:
            M_i = M_rxF_i 

        # Integrate forces
        self.BAeroF_i[it, :nB , 0] = trapezoid(F_i[:,:,0], sBldNodes)
        self.BAeroF_i[it, :nB , 1] = trapezoid(F_i[:,:,1], sBldNodes)
        self.BAeroF_i[it, :nB , 2] = trapezoid(F_i[:,:,2], sBldNodes)
        # Integrate moments and add force contrib
        self.BAeroM_i[it, :nB , 0] = trapezoid(M_i[:,:,0], sBldNodes)
        self.BAeroM_i[it, :nB , 1] = trapezoid(M_i[:,:,1], sBldNodes)
        self.BAeroM_i[it, :nB , 2] = trapezoid(M_i[:,:,2], sBldNodes)
        # Sum over blade
        self.RtAeroF_i[it, :3] =  np.sum(self.BAeroF_i[it, : , :],axis=0)
        self.RtAeroM_i[it, :3] =  np.sum(self.BAeroM_i[it, : , :],axis=0)

        # --- In hub coordinate system
        for iB in range(nB):
            self.RtAeroF_r[it] =  R_r2g.T .dot( self.RtAeroF_i[it, :3] )
            self.RtAeroM_r[it] =  R_r2g.T .dot( self.RtAeroM_i[it, :3] )

        # --- Integral quantities for rotor
        #        # TODO Might need rethinking
        #        if self.bUseCm:
        #        else:
        #            QMz = np.zeros(nB)
        #
# #       m%AllOuts( RtAeroPwr ) = omega * m%AllOuts( RtAeroMxh )
        #self.BladeThrust[it,:] = trapezoid(self.Fn[it,:]    , r_p)       # Normal to rotor plane
        #self.BladeTorque[it,:] = trapezoid(self.Ft[it,:]*r_p, r_p) + QMz[:]  # About shaft 
        #self.Thrust[it] = np.sum(trapezoiddeThrust[it,:])            # Normal to rotor plane
        #self.Torque[it] = np.sum(self.BladeTorque[it,:])
        self.Thrust[it] = self.RtAeroF_r[it,0]
        self.Torque[it] = self.RtAeroM_r[it,0] # NOTE: coord sys
        self.Power[it]  = Omega*self.RtAeroM_r[it, 0]
            # TODO TODO
            #self.BladeEdge   = np.zeros((nt,nB))
            #self.BladeFlap   = np.zeros((nt,nB))
                #         if (WT.Sources.Format=='wtperf'):
                #             RES.BladeThrust[idB] = sum(np.multiply(Rotor.dr,(Pn * np.cos(np.pi/180*cone))))
                #             RES.BladeTorque[idB] = sum(np.multiply(np.multiply(Rotor.dr,Pt),(r * np.cos(np.pi/180*cone))))
                #         else:
                #             RES.BladeTorque[idB] = getTorqueFromBlade(r,Pt * np.cos(np.pi/180*cone),R)
                #             RES.BladeThrust[idB] = getThrustFromBlade(r,Pn * np.cos(np.pi/180*cone),R)
                #         RES.BladeFlap[idB] = sum(np.multiply(np.multiply(Rotor.dr,(np.transpose(Pn) * np.cos(np.pi/180*cone))),(r - rhub)))
                #         RES.BladeEdge[idB] = sum(np.multiply(np.multiply(np.multiply(Rotor.dr,np.transpose(Pt)),(r * np.cos(np.pi/180*cone))),(r - rhub)))
                # ### Torque momentum at hub
                # #RES.BladeFlap(idB)=trapezoidotor.r;R],[Rotor.r.*(Pz*0+Pn);0]);
                # #RES.BladeEdge(idB)=trapezoidotor.r;R],[Rotor.r.*(Py*0+Pt);0]);
                #  #### Returning AerotrapezoidForces
                #  RES.Flap = sum(RES.BladeFlap)
                #  RES.Edge = sum(RES.BladeEdge)
        return xd1

    def simulationConstantRPM(self, time, RPM, windSpeed=None, windExponent=None, windRefH=None, windFunction=None, cone=None, tilt=None, yaw=None, hubHeight=None, firstCallEquilibrium=True,
            BldNd_BladesOut=None, BldNd_BlOutNd=None  # Blade outputs
            ):
        """ 
        wrapper function to perform a simple simulation at constant RPM
       
        INPUTS:
        Different ways to specify wind:
          - windSpeed: scalar, wind speed (at hub height) along x
          - windExponent power law exponent speed for wind speed (None=uniform wind)
          - windRefH reference height for power law
        OR
          - windFunction: function with interface: f(x,y,z,t)=u,v,w
                       with x,y,z,u arrays of arbitrary shapes

        - cone: override cone values, in degrees, OpenFAST convention [deg]
        - tilt: override tilt values, in degrees, OpenFAST convention (negative about yn !) [deg]
        - hubHeight: if provided, override the hubheight that is computed based on OverHang, Twr2Shaft and Tilt
        - firstCallEquilibrium: if true, the inductions are set to the equilibrium values at t=0 (otherwise,0)

        """
        # --- Define motion
        motion = PrescribedRotorMotion()
        motion.init_from_BEM(self, tilt=tilt, cone=cone, psi0=0, yaw=yaw)
        motion.setType('constantRPM', RPM=RPM)
        if hubHeight is not None:
            motion.origin_pos_gl0=np.array([0,0,hubHeight])

        # --- Define wind function
        if windFunction is None:
            if windExponent is None:
                windFunction = lambda x,y,z,t : (np.ones(x.shape)*windSpeed, np.zeros(x.shape), np.zeros(y.shape))
            else:
                if windRefH is None:
                    raise Exception('Hub height needs to be provided')
                windFunction = lambda x,y,z,t : (np.ones(x.shape)*windSpeed*(z/windRefH)**windExponent, np.zeros(x.shape), np.zeros(x.shape))

        # --- Perform time loop
        dt=time[1]-time[0]
        xdBEM = self.getInitStates()
        self.timeStepInit(time[0],time[-1],dt) 
        for it,t in enumerate(self.time):
            motion.update(t)
            u,v,w = windFunction(motion.pos_gl[:,:,0], motion.pos_gl[:,:,1], motion.pos_gl[:,:,2], t)  
            Vwnd_g = np.moveaxis(np.array([u,v,w]),0,-1) # nB x nr x 3
            xdBEM = self.timeStep(t, dt, xdBEM, motion.psi, motion.psi_B0,
                    motion.origin_pos_gl, motion.omega_gl, motion.R_SB2g, 
                    motion.R_bld2SB, # From blades 2 Shaft-Blade
                    motion.pos_gl, motion.vel_gl, motion.R_a2g,
                    Vwnd_g,
                    firstCallEquilibrium= it==0 and firstCallEquilibrium,
                    )
            #if np.mod(t,1)<dt/2:
            #    print(t)
            # --- Aditional storage
        df = self.toDataFrame(BldNd_BladesOut=BldNd_BladesOut, BldNd_BlOutNd=BldNd_BlOutNd)
        return df

# --------------------------------------------------------------------------------}
# --- Utils common between steady and unsteady BEM
# --------------------------------------------------------------------------------{

def ProjectClCd_to_a(Cl, Cd, alpha):
    C_xa =  Cl*cos(alpha)+ Cd*sin(alpha) 
    C_ya = -Cl*sin(alpha)+ Cd*cos(alpha)
    return C_xa, C_ya

def Transform_ClCd_to_CxCy(Cl, Cd,  phi):
    """ Transform the aerodynamic coefficients (Cl,Cd,Cm) (directed based on Vrel_xy_a )
    from the airfoil coordinate system (a) to the without sweep pitch coordinate system (w)
    NOTE: "Cy" is currently "-Cyw" """
    cphi = cos(phi)
    sphi = sin(phi)
    # resolve into normal (x) and tangential (y) forces
    Cx = Cl*cphi + Cd*sphi # Cx = Cxw
    Cy = Cl*sphi - Cd*cphi  # Cy = -Cyw
    return Cx, Cy

def Transform_ClCdCm_to_CxCyCzCmxCmyCmz(Cl, Cd, Cm, toe, cant, theta, alpha):
    """
    Transform the aerodynamic coefficients (Cl,Cd,Cm) (directed based on Vrel_xy_a )
    from the airfoil coordinate system (a) to the polar coordinate system (p)
    NOTE: "Cy" is currently "-Cyp" """
    xa_p, ya_p, za_p =  airfoilAxesInPolar(toe, cant, theta)
    #R_p2a = rotPolar2Airfoil(tau, kappa, theta)
    # transform force coefficients into airfoil frame
    Cxa, Cya = ProjectClCd_to_a(Cl, Cd, alpha)
    # Put force coefficients back into rotor plane reference frame
    Cxp  =   Cxa * xa_p[0]  + Cya*ya_p[0]   #  Cxp  and  cn
    mCyp =- (Cxa * xa_p[1]  + Cya*ya_p[1])  # -Cyp      ct
    Czp  =   Cxa * xa_p[2]  + Cya*ya_p[2]   #  Czp
    # Put moment coefficients into the rotor reference frame
    Cmx = Cm * za_p[0]
    Cmy = Cm * za_p[1]
    Cmz = Cm * za_p[2]
    return Cxp, mCyp, Czp, Cmx, Cmy, Cmz
 

def _fInductionCoefficients(Vrel_norm, V0, F, cnForAI, ctForTI,
        lambda_r, sigma, phi, relaxation=0.4, a_last=None, bSwirl=True, 
        drdz=1, algorithm='legacy', 
        CTcorrection='AeroDyn15', swirlMethod='AeroDyn'):
    """Compute the induction coefficients

        Inputs
        ----------
        a_last    : last iteration axial induction factor
        Vrel_norm : normed relative velocity
        V0        : free stream velocity
        F         : total loss
        cnForAI   : normal force coefficient
        ctForTI   : tangential force coefficient
        lambda_r  : speed ratio distribution
        sigma     : blade solidity
        phi       : flow angle [deg]
        relaxation: relaxation factor in axial induction factor
        bSwirl    : swirl flow model enabled / disabled

        Outputs
        ----------
        a: axial induction factor
        aprime: tangential induction factor
        Ct: local thrust coefficient
    """
    # --- Default a and CT
    if algorithm=='legacy':
        a = 1. / ((4.*F*sin(phi)**2)/(sigma*(cnForAI+10**-8))+1) # NOTE singularity avoided
    elif algorithm=='polarProj':
        a = 1. / ((4.*F*sin(phi)**2)/(drdz*sigma*(cnForAI+10**-8))+1) # NOTE singularity avoided
    else:
        raise NotImplementedError()
    # CT=(1-a_last).^2.*sigma.*CnForAI./((sind(phi)).^2)
    Ct = Vrel_norm**2 * sigma * cnForAI/(V0**2)  # that's a CT loc
    # AeroDyn
    #k = sigma*cn/4.0_ReKi/F/sphi/sphi
    #if (k <= 2/3) then  ! momentum state for a < 0.4
    # a = k/(1+k)

    # --- Hight thrust correction
    if CTcorrection=='GlauertCT':
        # Glauert correction as default
        #>>> NOTE this is:  a = a_Ct(Ct, a, method='Glauert') from highthrust
        ac = 0.3 
        bHigh = a > ac
        fg = 0.25*(5.-3.*a[bHigh])
        a[bHigh] = Ct[bHigh]/(4.*F[bHigh]*(1.-fg*a[bHigh]))
    else:
        a = a_Ct(Ct=Ct, a=a, F=F, method=CTcorrection)

    a[F<0.01]=1 # HACK to match aerodyn # TODO make that an option

    # --- Relaxation for high Ct
    if a_last is not None:
        bHigh = a>0.3
        a[bHigh] = a[bHigh]*relaxation + (1.-relaxation)*a_last[bHigh]

    # --- Swirl
    if bSwirl is True:
        if swirlMethod=='AeroDynOld':
            aprime=0.5*(sqrt(1+4*a*F*(1-a)/lambda_r**2)-1);

        elif swirlMethod=='AeroDyn':
            # NOTE: AeroDyn has more tests (e.g. if cos(phi)=0)
            aprime=np.zeros(a.shape)
            b0 = np.logical_or(np.abs(a-1)<1e-5, np.abs(phi)<1e-5)
            b1 = np.logical_not(b0)
            aprime[b0] = 0 
            kp         = sigma[b1]*ctForTI[b1]/(4*F[b1]*sin(phi[b1])*cos(phi[b1]))
            aprime[b1] = kp/(1-kp)
        elif swirlMethod=='HAWC2':
            aprime = (Vrel_norm**2*ctForTI*sigma)/(4.*(1.-a)*V0**2*lambda_r)
        elif swirlMethod=='Default': # Need a better name
            aprime=1/((4*F*sin(phi)*cos(phi)) /(sigma*ctForTI)  -1 );
        else:
            raise NotImplementedError()
    else:
        aprime = a * 0.

    # Bounding values for safety
    a     [np.isnan(a)]      = 0
    aprime[np.isnan(aprime)] = 0
    aprime = np.clip(aprime,-1,1.0) 
    a      = np.clip(a     ,-1,1.5)
    Ct     = np.clip(Ct    ,-1,3)
    return a, aprime, Ct



def clip_zeros(x, tol=1e-5):
    b = np.abs(x)<tol
    #x[b] = tol*np.sign(x[b]) # NOTE: sign return 0 if x==0
    x[b] = tol* (2*(x[b] >= 0) - 1)
    return x

# --------------------------------------------------------------------------------}
# --- Helper class to prescribe a motion
# --------------------------------------------------------------------------------{
from welib.yams.utils import R_x, R_y, R_z
class PrescribedRotorMotion():
    """ 
    Class to return:
        position, velocity, and orientation of all blade station

     - R_a2g : from airfoil to global (this is well defined, also called "n-t" system in AeroDyn)
     - R_s2g : from section to global (this is ill-defined), this coordinate is used to define the "axial" and "tangential" inductions

    """
    def __init__(self):
        # body kinematics, Body is "rotor"
        self.R_SB2g = np.eye(3) # orientation from body to global
        self.R_SB2g_0 = np.eye(3) # orientation from body to global at t=0 (typically: tilt)
        self.origin_pos_gl0 = np.array([0,0,0]) # body origin at t=0 (position of rotor center)
        self.origin_pos_gl = np.array([0,0,0]) # body origin   (position of rotor center)
        self.origin_vel_gl = np.array([0,0,0])
        self.omega_gl      = np.array([0,0,0])

        # Blades 
        self.R_bld2SB=None # rotation matrices from blades to shaft-blade (i.e. rotor), contains azimuth and cone

        # --- DATA
        self.tilt = None # About yn, negative of OpenFAST convention [rad]
        self.yaw = None 
        self.cone = None 


    def init_from_inputs(self,  nB, r, twist, rotorOrigin, tilt, cone, yaw=0, psi0=0, meanLineAC=None):
        """ 
        INPUTS:
         - tilt: tilt angle [deg], with OpenFAST convention, negative about yn
         - cone: cone angle [deg], with OpenFAST convention
 x_hat_disk.dot(elemPosRelToHub)

         - meanLineAC: position of AC for each blade, in each blade coordinate system (nB x nr x 3)
                        meanLineAC [iB, :, 0]: prebend
                        meanLineAC [iB, :, 1]: presweep
                        meanLineAC [iB, :, 2]: radial 
        """
        # TODO TODO Pitch!
        self.nB   =  nB
        if yaw is None:
            yaw=0
        self.yaw  =  yaw*np.pi/180
        self.cone =  cone*np.pi/180
        self.tilt = -tilt*np.pi/180  # input has OpenFAST convention negative about yn, add negative signe
        self.r     = r               # spanwise position, from hub center (typically r= HubRad->TipRad)
        self.twist = twist
        self.allocate()
        self.origin_pos_gl0 = rotorOrigin

        # Basic geometries for nacelle
        #self.R_SB2g_0 = R_y(self.tilt)  # Rotation fromShaft to Nacelle
        self.R_SB2g_0 = R_y(self.tilt).dot(R_z(self.yaw))  # Rotation fromShaft to Nacelle

        # Orientation of blades with respect to body
        self.R_bld2SB = [np.eye(3)]*self.nB
        self.psi_B0  = np.zeros(self.nB)
        for iB in np.arange(self.nB):
            self.psi_B0[iB]= psi0 + iB*2*np.pi/self.nB
            R_SB = R_x(self.psi_B0[iB]) 
            self.R_bld2SB[iB] = R_SB.dot(R_y(self.cone)) # blade2shaft / blade2rotor
            #self.R_ntr2b[iB] = R_SB.dot(np.array([[1,0,0],[0,-1,0],[0,0,1]])) # "n-t-r" to shaft

        # Set initial positions and orientations in body coordinates
        nr=len(self.r)
        if meanLineAC is None:
            meanLineAC = np.zeros((self.nB, nr, 3))
            meanLineAC[:,:,2]=self.r
        for iB in np.arange(self.nB):
            CrvAC = meanLineAC[iB,:,0]
            Cant = calcCantAngle(CrvAC, self.r-self.r[0])
            #print('r'    , self.r-self.r[0])
            #print('CrvAC', CrvAC )
            #print('Cant',  Cant )
            Cant *= np.pi/180

            for ir in np.arange(nr):
                self.pos0[iB,ir,:] = self.R_bld2SB[iB].dot(meanLineAC[iB,ir,:])
                R_a2bld = BodyXYZ_A(0, Cant[ir],  -self.twist[ir])
                #R_a2bld = BodyXYZ_A(0, 0,  -self.twist[ir])  # airfoil to blade, twist only
                #R_a2bld = R_z(-self.twist[ir])              # airfoil to blade, twist only
                self.R_a2SB_0[iB,ir,:,:] = self.R_bld2SB[iB].dot(R_a2bld)
                #R_g2a_0 =  (self.R_a2SB_0[iB,ir,:,:]).T .dot(self.R_SB2g_0.T) # "BladeMotion%Orientation"
                #print('orientationlL',R_a2bld[:,0])
                #print('orientationlL',R_a2bld[:,1])
                #print('orientationlL',R_a2bld[:,2])
                #print('BLADE MOTION REF ORI',R_g2a_0[0,:])
                #print('BLADE MOTION REF ORI',R_g2a_0[1,:])
                #print('BLADE MOTION REF ORI',R_g2a_0[2,:])
                #self.R_s2SB_0[iB,ir,:,:] = self.R_a2SB_0[iB,ir,:,:] #self.R_bld2SB[iB]  # TODO without pitch twist

    def init_from_BEM(self, BEM, tilt=None, cone=None, psi0=0, yaw=None):
        """ 
        Initializes motion from a BEM class
        Possibility to override the tilt and cone geometry:
        INPUTS:
         - tilt: tilt angle [deg], with OpenFAST convention, negative about yn
         - cone: cone angle [deg], with OpenFAST convention
        
        """
        if tilt is None:
            tilt=BEM.tilt0
        if cone is None:
            cone=BEM.cone0

        zTT = BEM.TowerHt+BEM.Twr2Shft
        rotorOrigin =[BEM.OverHang*np.cos(-tilt*np.pi/180), 0, -BEM.OverHang*np.sin(-tilt*np.pi/180)]
        #print('Rotor Origin', rotorOrigin)
        if yaw is not None:
            rotorOrigin = R_z(yaw).dot(rotorOrigin)
        #print('Rotor Origin', rotorOrigin)
        rotorOrigin[2] += zTT


        meanLineAC = BEM.meanLineAC
        self.init_from_inputs(BEM.nB, BEM.r, BEM.twist, rotorOrigin, tilt=tilt, cone=cone, yaw=yaw, psi0=0, meanLineAC=meanLineAC)

    def allocate(self):
        nr = len(self.r)
        # Section nodes
        self.pos0    = np.zeros((self.nB,nr,3))   # position of nodes at t= 0 in body coordinates 
        #self.R_s2SB_0  = np.zeros((self.nB,nr,3,3)) # Orientation section to (Shaft-Blade body) at t=0
        self.R_a2SB_0  = np.zeros((self.nB,nr,3,3)) # Orientation airfoil to (Shaft-Blade body) at t=0
        self.pos_gl = np.zeros((self.nB,nr,3))   # position of all nodes
        self.vel_gl = np.zeros((self.nB,nr,3))   # linear velocities
        self.R_a2g  = np.zeros((self.nB,nr,3,3)) # Orientation airfoil to global

    def setType(self, sType, **kwargs):
        self.sType=sType
        self.opts=kwargs

    def rigidbodyKinUpdate(self, P_gl, vel_gl, omega_gl, R_SB2g):
        """
        Update position, velocities, orientations of nodes assuming a rigid body motion
        of the origin given as input
        """
        self.origin_pos_gl = P_gl
        self.origin_vel_gl = vel_gl
        self.omega_gl      = omega_gl
        self.R_SB2g         = R_SB2g

        # Update of positions
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                s_OP =  R_SB2g.dot(self.pos0[iB,ir,:])
                self.pos_gl[iB,ir,:] = P_gl   + s_OP
                self.vel_gl[iB,ir,:] = vel_gl + np.cross(omega_gl, s_OP)
                self.R_a2g[iB,ir,:,:] = R_SB2g.dot(self.R_a2SB_0[iB,ir,:])

    def update(self, t):
        if self.sType=='constantRPM':
            omega = self.opts['RPM']*2.*np.pi/(60.)
            psi = t*omega
            pos = self.origin_pos_gl0
            vel = np.array([0,0,0])
            ome = self.R_SB2g_0.dot(np.array([omega,0,0]))
            R_SB2g = self.R_SB2g_0.dot(R_x(psi))
            self.rigidbodyKinUpdate(pos, vel, ome, R_SB2g)
            self.psi=psi # hack

        elif self.sType=='x-oscillation':
            omega = self.opts['frequency']*2.*np.pi
            A = self.opts['amplitude']
            x    = A*np.sin(omega*t)
            xdot = A*omega*np.cos(omega*t)
            pos = self.origin_pos_gl0+np.array([x,0,0])
            vel = np.array([xdot,0,0])
            ome = np.array([0,0,0])
            R_SB2g =self.R_SB2g_0
            self.rigidbodyKinUpdate(pos, vel, ome, R_SB2g)

        elif self.sType=='constantRPM x-oscillation':
            omega = self.opts['frequency']*2.*np.pi
            A = self.opts['amplitude']
            x    = A*np.sin(omega*t)
            xdot = A*omega*np.cos(omega*t)
            omegapsi = self.opts['RPM']*2.*np.pi/(60.)
            psi = t*omegapsi
            pos = self.origin_pos_gl0+np.array([x,0,0])
            vel = np.array([xdot,0,0])
            ome = self.R_SB2g_0.dot(np.array([omegapsi,0,0]))
            R_SB2g = self.R_SB2g_0.dot(R_x(psi))
            self.rigidbodyKinUpdate(pos, vel, ome, R_SB2g)
            self.psi=psi # hack
        else:
            raise NotImplementedError(sType)

    def plotCurrent(self, ax=None, fig=None, lines=None):
        """ Plot current blade positions """
        if fig is None and ax is None:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax  = fig .add_subplot(111, projection='3d', proj_type='ortho')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.tick_params(direction='in')
            ax.axis('equal')
            ax.set_aspect('equal')
        if lines is None:
            lines=[]
            for iB in np.arange(self.nB):
                lnp, = ax.plot(self.pos_gl[iB,:,0], self.pos_gl[iB,:,1], self.pos_gl[iB,:,2])
                lines.append(lnp)
        else:
            for iB, l in enumerate(lines):
                l.set_data(motion.pos_gl[iB,:,0:2].T)
                l.set_3d_properties(motion.pos_gl[iB,:,2])
        return fig, ax, lines
    
    def plotAnim(self, t0, tmax, dt, interval=10, xlim=None,ylim=None,zlim=None):
        """ Animate motion """
        # 3d animation
        fig, ax, lnp = self.plotCurrent()
        def init():
            if xlim is not None:
                ax.set_xlim3d(xlim)
            if ylim is not None:
                ax.set_ylim3d(ylim)
            if zlim is not None:
                ax.set_zlim3d(zlim)
            ax.set_aspect('equal')
            ax.axis('equal')
            return lnp
        def update(i):
            t=t0+i*dt
            self.update(t)
            self.plotCurrent(ax, lines=lnp)
            return lnp
        ani = FuncAnimation(fig, update, frames=np.arange(0,int((tmax-t0)/dt)), init_func=init, blit=False, interval=interval)
        return ani

# --------------------------------------------------------------------------------}
# --- Geometrical Utils 
# --------------------------------------------------------------------------------{
def rotPolar2Airfoil(tau, kappa, beta):
    """ Return R_p2a based on three angles. Watch out for Beta convention!
    INPUTS:
     - tau:   toe angle [rad]
     - kappa: cant angle [rad]
     - beta:  twist angle with wind turbine convention (negative about zb) [rad]
    """
    M33= np.array([
        [ cos(kappa)*cos(beta),   -cos(tau)*sin(beta)+sin(tau)*sin(kappa)*cos(beta),  -sin(tau)*sin(beta) - cos(tau)*sin(kappa)*cos(beta)],
        [  cos(kappa)*sin(beta),    cos(tau)*cos(beta)+sin(tau)*sin(kappa)*sin(beta),   sin(tau)*cos(beta) - cos(tau)*sin(kappa)*sin(beta)],
        [  sin(kappa)          ,   -sin(tau)*cos(kappa)                             ,        cos(tau)*cos(kappa)                        ]]
        , dtype='object'
        )
    return M33.astype(float)
    #return EulerConstruct(tau, kappa, -beta):
    #return BodyXYZ_DCM(tau, kappa, -beta)

def airfoilAxesInPolar(tau, kappa, beta):
    """ 
    See OpenFAST routine getAirfoilOrientation
    INPUTS:
     - tau:   toe angle [rad]
     - kappa: cant angle [rad]
     - beta:  twist angle with wind turbine convention (negative about zb) [rad]
    OUTPUTS:
     - xa_p : 3 x nB x nr  if  toe is nB x nr
     - ya_p :
     - za_p :
    """
    R_p2a = rotPolar2Airfoil(tau, kappa, beta)
    xa_p = R_p2a[0,:] # "afNormalVec"
    ya_p = R_p2a[1,:] # "afAxialVec"
    za_p = R_p2a[2,:] # "afRadialVec"
    return xa_p, ya_p, za_p

def polarCoord(P_rotor, x_hat_disk, P_bld):
    """ 
    Compute transformation matrix from global to polar for each blade node

    See OPENFAST subroutine Calculate_MeshOrientation_Rel2Hub
    INPUTS:
      - P_rotor:   3-vector, position of rotor center in global
      - x_hat_disk: 3-vector, normal vector, hub or shaft x axis
      - P_bld: 3-vector, position of blade node in global
    OUTPUTS:
      - R_g2p: from global to polar, 3x3
    """
    R_g2p = np.zeros((3, 3))
    # Project element position onto the rotor plane
    elemPosRelToRotor = P_bld - P_rotor 
    elemPosRotorProj = elemPosRelToRotor - x_hat_disk * ( x_hat_disk.dot(elemPosRelToRotor) )

    # Get unit vectors of the local annulus reference frame
    z_hat_annulus = elemPosRotorProj/ np.linalg.norm( elemPosRotorProj )
    x_hat_annulus = x_hat_disk
    y_hat_annulus = np.cross( z_hat_annulus, x_hat_annulus )

    # Form a orientation matrix for the annulus reference frame
    R_g2p[0, :]  = x_hat_annulus
    R_g2p[1, :]  = y_hat_annulus
    R_g2p[2, :]  = z_hat_annulus
    return R_g2p, elemPosRelToRotor, elemPosRotorProj


def curvilinearCoord(P_r, P_root, P_b):
    """
    Approximate curvilinear coordinates along the blade, 0 at the rotor center (not blade root)
    # See Init_BEMTmodule, "zRoot, zTip, zLocal"
    INPUTS:
      - P_r:   3-vector, position of rotor center in global coordinates
      - P_root:   3-vector, position of blade root in global coordinates
                  Note: typically P_root=P_b[0] but it might not necesarily be the case
      - P_b: (nrx3)-array, positions of blade nodes in global coordinates
    OUTPUTS:
      -
    """
    assert(P_b.shape[1]==3)
    nr  = len(P_b)
    s = np.zeros(nr)
    s_h = norm(P_root - P_r)  # hub
    s[0] = s_h + norm(P_b[0] - P_root)  # NOTE: second term is usually zero if blade starts at root
    for ie in range(nr-1):
        s[ie+1] = s[ie] + norm( P_b[ie+1] - P_b[ie] )

    return s_h, s


def calcCantAngle(f, xi, stencilSize=3):
    """ 
    Compute Cant angle from prebend using second order accurate gradient

    INPUTS:
     - f : (prebend) function values , array of size n
     - xi: (radius)  x value where the function is evaluated, array of size n
     - stencilSize: stencil for ! For 
    OUTPUTS:
     - cantAngel in [deg]
    """
    from welib.mesh.gradient import differ_stencil
    #dfdx = np.gradient(f,xi)
    #dxdx = np.gradient(xi,xi)
    #crvAng = np.degrees(np.arctan2(dfdx,dxdx))
    #print('dfdx',dfdx)
    #print('dxdx',dxdx)
    #print('dfdx',dfdx*180/np.pi)
    #print('crvAng',crvAng)
    #dx = np.gradient(xi)
    #df = np.gradient(f)
    #crvAng = np.degrees(np.arctan2(df,dx))
    #print('crvAng',crvAng)

    nx = len(xi)
    ns = stencilSize

    cPrime   = np.zeros(nx)
    fPrime   = np.zeros(nx)
    cantAngle = np.zeros(nx)
    cx   = np.zeros(ns)
    cf   = np.zeros(ns)
    xiIn = np.zeros(ns)
    fIn  = np.zeros(ns)

    for i in range(nx):
        if i==0:
            fIn  = f [0:ns]
            xiIn = xi[0:ns]
        elif i==nx-1:
            fIn  = f [nx-ns:nx]
            xiIn = xi[nx-ns:nx]
        else:
            fIn  = f [i-1:i+2]
            xiIn = xi[i-1:i+2]
        cx = differ_stencil( xi[i], 1, 2, xiIn)
        for j in range(ns):
            cPrime[i] += cx[j]*xiIn[j]
            fPrime[i] += cx[j]*fIn [j]            
    cantAngle = np.arctan2(fPrime, cPrime)*180/np.pi
    return cantAngle


def phi_kappa(phi, kappa):
    y = sin(phi)
    x = cos(phi)*cos(kappa)
    b = np.logical_not(np.logical_and(y== 0, y == 0))
    phi_k = np.zeros_like(phi)
    #phi_k[b] = 0
    phi_k[b] = np.arctan2(y[b], x[b])
    return phi_k




def ReynoldsNumberUA(BEMMod, axInduction, tanInduction, Vx, Vy, Vz, chord, nu, toe, cant, theta):
    """ Return Reynolds number in million
    NOTE: ugly  funciton, temporary for now to match AeroDyn
    """
    Vx_p = Vx*(1-axInduction) 
    Vy_p = Vy*(1+tanInduction)
    Vz_p = Vz
    ## Project inflow vector onto airfoil plane
    if BEMMod=='polarProj':
        xa_p, ya_p, za_p = airfoilAxesInPolar(toe, cant, theta)
        #xa_p = xa_p.astype(float)
        #ya_p = ya_p.astype(float)
        #za_p = za_p.astype(float)
        dotProduct = Vx_p*za_p[0,:,:] + Vy_p*za_p[1,:,:] + Vz_p*za_p[2,:,:]
        Vrel_x = Vx_p - dotProduct*za_p[0]  
        Vrel_y = Vy_p - dotProduct*za_p[1]  
        Vrel_z = Vz_p - dotProduct*za_p[2]  
        #inflowVecInAirfoilPlane = inflowVec - dot_product( inflowVec, za_p ) * za_p 
    else:
        raise Exception()
#          # TODO TODO TODO EB CHECK THAT THE SAME MIGHT BE OBTAINED IF cant=0, toe=0
#          inflowVecInAirfoilPlane(1) = inflowVec(1)
#          inflowVecInAirfoilPlane(2) = inflowVec(2)
#          inflowVecInAirfoilPlane(3) = 0.0_ReKi
    # Wxy: resultant velocity in the xy airfoil plane.
    Wxy = np.sqrt(Vrel_x**2 + Vrel_y**2)
    Re =  Wxy * chord / nu
    Re[Re<0.001] = 0.001 
    Re = Re/10**6
    return Re

def RelativeVelocityUA( axInduction, tanInduction, Vx, Vy, cantAngle, xVelCorr=0):
      v_ac_x = (Vx*cos(cantAngle)+xVelCorr)*(1 - axInduction)
      v_ac_y =                           Vy*(1 + tanInduction)
      Vrel = np.sqrt(v_ac_x**2 + v_ac_y**2)
      return Vrel, v_ac_x, v_ac_y




if __name__=="__main__":
    """ See examples/ for more examples """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.animation import FuncAnimation
    import welib.essentials

    # --- CONE DEBUG
#     inputfile = 'C:/W0/Work/2018-NREL/BAR-Cone-BEM/simulations_parametric/cone/np2_DEBUG/ad_driver_p{:02d}.dvr'.format(cone)
#     cone = 0; 
#     BEM = UnsteadyBEM(inputfile)
#     BEM.algorithm='polarProj'
#     BEM.bYawModel = yaw!=0 # Yaw correction
#     BEM.bDynaStall = False # dynamic stall model
#     BEM.bDynaWake = False # dynamic stall model
#     print(BEM)
#     time=np.arange(0,10,0.1)
#     windSpeed = 10
#     RPM       = 10
#     dfTS = BEM.simulationConstantRPM(time, RPM, windSpeed=windSpeed, cone=cone, firstCallEquilibrium=True, BldNd_BladesOut=1)
#     dfTS.to_csv(inputfile.replace('.dvr','_welib.csv'), index=False)

    # --- YAW DEBUG
    yaw = np.array([0,50])
    yaw = 50

    # --- UNSTEADY BEM
    inputfile = 'C:/W0/Work/2018-NREL/BAR-Cone-BEM/simulations_parametric/yaw/250_DEBUG/ad_driver_p{:02d}.dvr'.format(yaw)
    BEM = UnsteadyBEM(inputfile)
    BEM.algorithm='polarProj'
    BEM.bYawModel = True # Yaw correction
    BEM.bDynaStall = False # dynamic stall model
    BEM.bDynaWake = False # dynamic stall model
    print(BEM)
    time=np.arange(0,10,0.1)
    windSpeed = 10
    RPM       = 10
    print('>>> YAW ',yaw)
    dfTS = BEM.simulationConstantRPM(time, RPM, windSpeed=windSpeed, yaw=yaw, firstCallEquilibrium=True, BldNd_BladesOut=1)
    dfTS.to_csv(inputfile.replace('.dvr','_welib.csv'), index=False)






















# 
# 
# 
# 
# 
# 
#     # --- Read a FAST model to get Aerodynamic parameters
#     BEM = UnsteadyBEM('./Main_Onshore.fst')
#     BEM.CTcorrection='AeroDyn' # High Thrust correction
#     BEM.swirlMethod ='AeroDyn' # type of swirl model
# #    BEM.swirlMethod ='HAWC2' # type of swirl model
# #     BEM.bSwirl = True  # swirl flow model enabled / disabled
#     BEM.WakeMod=1 # 0: no inductions, 1: BEM inductions
# #     BEM.bTipLoss = True # enable / disable tip loss model
# #     BEM.bHubLoss = False # enable / disable hub loss model
# #     BEM.bTipLossCl = False # enable / disable Cl loss model
# #     BEM.TipLossMethod = 'Glauert'  # type of tip loss model
# #     BEM.bDynaStall = True # dynamic stall model
#     BEM.bDynaWake = True # dynamic inflow model
# #    BEM.bDynaWake = True # dynamic inflow model
#     BEM.bYawModel = True # Yaw correction
# #     BEM.bYawModel = False # Yaw correction
# #     BEM.bAIDrag = True # influence on drag coefficient on normal force coefficient
# #     BEM.bTIDrag = True # influence on drag coefficient on tangential force coefficient
#     BEM.relaxation = 0.3
# 
#     time=np.arange(0,10,0.05)
#     RPM=10
#     df= BEM.simulationConstantRPM(time, RPM, windSpeed=10, windExponent=0.0, windRefH=90, windFunction=None, cone=None, tilt=None, firstCallEquilibrium=True)
#     df.to_csv('_BEM.csv', index=False)

#     # --- Read a FAST model to get structural parameters for blade motion
#     motion = PrescribedRotorMotion()
#     motion.init_from_FAST('./Main_Onshore.fst', tilt=0, cone=0)
#     motion.setType('constantRPM', RPM=10.0)
#     #motion.setType('constantRPM x-oscillation', RPM=12.1, frequency=1.1, amplitude=20)
# 
#     dt=0.05
#     dtRadOut=1.0
#     tmax=10
#     
#     xdBEM = BEM.getInitStates()
# 
#     # Allocate
#     BEM.timeStepInit(0,tmax,dt) 
# 
#     for it,t in enumerate(BEM.time):
#         motion.update(t)
#         xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
#                 motion.origin_pos_gl, motion.omega_gl, motion.R_SB2g, 
#                 motion.R_ntr2g, motion.R_bld2SB,
#                 motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g,
#                 firstCallEquilibrium=it==0
# #                 firstCallEquilibrium=it==0
#                 )
#         if np.mod(t,dtRadOut)<dt/2:
#             #print(t)
#             #dfRad = BEM.toDataFrameRadial(it)
#             #dfRad.to_csv('_BEMRad_t{:04d}.csv'.format(it), index=False)
#             pass
#     df = BEM.toDataFrame()
#     df.to_csv('_BEM.csv', index=False)
# 
#     #ani = motion.plotAnim(0,10,0.01, ylim=[-80,80], zlim=[-80,80])
#     plt.show()

