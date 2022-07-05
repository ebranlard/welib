""" 
Python implemenation of an unsteady BEM code

Default models are provided in this file so it can be shipped as "standalone".
More models are provided based on other "welib" packages.

Reference:
   [1]: Branlard, 2017, Wind Turbines Aerodynamics and Vorticity Based Methods: Fundamentals and recent applications, Springer

"""
import numpy as np
from numpy import cos, sin, arctan2, pi, arccos, exp, abs, min, sqrt
from scipy.interpolate import interp1d
import copy
import pandas as pd

# Load more models
# try:
from welib.BEM.highthrust import a_Ct
# except: 
#     pass


def _fInductionCoefficients(Vrel_norm, V0, F, cnForAI, ctForTI,
        lambda_r, sigma, phi, relaxation=0.4, a_last=None, bSwirl=True, CTcorrection='AeroDyn', swirlMethod='AeroDyn'):
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
    a = 1. / ((4.*F*sin(phi)**2)/(sigma*(cnForAI+10**-8))+1) # NOTE singularity avoided
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
        a = a_Ct(Ct, a, F, method=CTcorrection)

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
        self.a       = np.zeros((nB,nr)) # axial induction
        # Dynamic wake
        self.Vind_qs_p  = np.zeros((nB,nr,3)) # Quasi-steady velocity, polar coordinates
        self.Vind_qs_g  = np.zeros((nB,nr,3)) # Quasi-steady velocity
        self.Vind_int_p = np.zeros((nB,nr,3)) # Intermediate velocity, polar coordinates
        self.Vind_dyn_p = np.zeros((nB,nr,3)) # Dynamic induced velocity (before skew/yaw), polar coordinates
        self.Vind_dyn_g = np.zeros((nB,nr,3)) # Dynamic induced velocity (before skew/yaw), global coordinates
        # Dynamic stall
        self.fs = np.zeros((nB,nr)) # Separation 

class AeroBEM:
    """ 
    Perform unsteady BEM calculations
    """
    def __init__(self):

        # Aero Data
        self.chord  = None
        self.polars = None
        # Environment
        self.rho     = None
        self.kinVisc = None

        # Structural "Inputs"
        # position, velocity, and orientation of all blade station
        #self.cone   = None
        #self.twist  = None
        self.nB     = None
        self.r      = None # radial stations

        self.setDefaultOptions()

    def setDefaultOptions(self):
        self.nbIt = 200  # maximum number of iterations in BEM
        self.aTol = 10 ** -6 # tolerance for axial induction factor convergence
        self.relaxation = 0.5  # relaxation factor in axial induction factor
        self.CTcorrection = 'AeroDyn'  #  type of CT correction more model implementated in the future like 'spera'
        self.swirlMethod  = 'AeroDyn' # type of swirl model
        self.Ngrid = 1.0
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
        self.bReInterp = False # interpolate the input tabulated airfoil data for Reynolds variation
        self.bThicknessInterp = True # interpolate the input tabulated airfoil data for thickness variation
        self.WakeMod=1 # 0: no inductions, 1: BEM inductions
        self.bRoughProfiles = False # use rough profiles for input airfoil data

    def init_from_FAST(self, FASTFileName):
        from welib.weio.fast_input_deck import FASTInputDeck
        F = FASTInputDeck(FASTFileName,readlist=['AD','ED','ADbld','AF'])

        # Environment
        try:
            self.rho     = float(F.fst['AirDens'])  # New OF > 3.0
        except:
            self.rho     = float(F.AD['AirDens'])   # Old OF <=3.0
        try:
            self.kinVisc = float(F.fst['KinVisc'])  # New OF > 3.0
        except:
            self.kinVisc = float(F.AD['KinVisc'])   # Old OF <= 3.0


        # Aerodynamics
        self.nB    = F.ED['NumBl']
        self.r     = F.AD.Bld1['BldAeroNodes'][:,0] + F.ED['HubRad']
        chord      = F.AD.Bld1['BldAeroNodes'][:,-2]
        self.chord = np.stack([chord]*self.nB)
        self.twist = F.AD.Bld1['BldAeroNodes'][:,-3]*np.pi/180
        polars=[]
        ProfileID=F.AD.Bld1['BldAeroNodes'][:,-1].astype(int)
        for ipolar in  ProfileID:
            polars.append(F.AD.AF[ipolar-1]['AFCoeff'])
        self.polars = polars

        # Input geometry (may be overriden by motion/simulations)
        self.cone0     = F.ED['PreCone(1)']
        self.tilt0     = F.ED['ShftTilt']
        self.TowerHt   = F.ED['TowerHt']
        self.OverHang  = F.ED['OverHang']
        self.Twr2Shft  = F.ED['Twr2Shft']

        self._init()

    def _init(self):
        # Creating interpolation functions for each polar, now in rad!
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
        self.Cl_qs  = np.zeros((nt,nB,nr))
        self.Cd_qs  = np.zeros((nt,nB,nr))
        self.Cl     = np.zeros((nt,nB,nr))
        self.Cd     = np.zeros((nt,nB,nr))
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
        self.Vind_p = np.zeros((nt,nB,nr,3))
        self.Vind_s = np.zeros((nt,nB,nr,3))
        self.Vind_qs_p = np.zeros((nt,nB,nr,3))
        self.Vflw_p = np.zeros((nt,nB,nr,3)) # Vwnd-Vstr
        self.Vflw_s = np.zeros((nt,nB,nr,3)) # Vwnd-Vstr
        self.Vwnd_p = np.zeros((nt,nB,nr,3))
        self.Vwnd_s = np.zeros((nt,nB,nr,3))
        self.Vwnd_a = np.zeros((nt,nB,nr,3))
        self.Vstr_p = np.zeros((nt,nB,nr,3))
        self.Vstr_s = np.zeros((nt,nB,nr,3))
        self.Vstr_xa = np.zeros((nt,nB,nr))
        self.Vstr_ya = np.zeros((nt,nB,nr))
        self.Vrel   = np.zeros((nt,nB,nr))
        self.AxInd  = np.zeros((nt,nB,nr))
        self.TnInd  = np.zeros((nt,nB,nr))
        # Loads per span
        self.L      = np.zeros((nt,nB,nr))
        self.D      = np.zeros((nt,nB,nr))
        self.Fn     = np.zeros((nt,nB,nr))
        self.Ft     = np.zeros((nt,nB,nr))
        self.F_a   = np.zeros((nt,nB,nr,3))
        self.F_s    = np.zeros((nt,nB,nr,3))
        self.Gamma  = np.zeros((nt,nB,nr))
        self.alpha  = np.zeros((nt,nB,nr))
        self.phi    = np.zeros((nt,nB,nr))
        self.Re     = np.zeros((nt,nB,nr))
        # Integrated values
        self.Thrust   = np.zeros(nt)
        self.Torque   = np.zeros(nt)
        self.Power    = np.zeros(nt)
        self.chi      = np.zeros(nt)
        self.chi0     = np.zeros(nt)
        self.RtVAvg   = np.zeros((nt,3))
        self.psi      = np.zeros(nt)
        self.RtArea   = np.zeros(nt)
        self.SkewAzimuth  = np.zeros((nt,nB))
        # Blade blades
        self.BladeTorque = np.zeros((nt,nB))
        self.BladeThrust = np.zeros((nt,nB))
        self.BladeEdge   = np.zeros((nt,nB))
        self.BladeFlap   = np.zeros((nt,nB))



    def calcOutput(self):
        R = np.sqrt(self.RtArea/pi)
        q = 0.5*self.rho*self.RtArea*self.RtVAvg[:,0]**2
        self.CT=self.Thrust/(q)
        self.CQ=self.Torque/(q*R)
        self.CP=self.Power /(q*self.RtVAvg[:,0])


    def toDataFrame(self, BldNd_BladesOut=None, BldNd_BlOutNd=None, BldNd_OutList=None):
        """ Export time series to a pandas dataframe
        Column names are set to match OpenFAST outputs

        BldNd_BladesOut: if None, all blades nodal output is output. Otherwise set it to an integer < self.nB
        BldNd_BlOutNd:   if None, all radial position are output.    Otherwise set it to a list of radial indices [0,2,4, etc] 
        BldNd_OutList:   if None, all possible channels are output.  Otherwise, set it to a list of channel names with units. 
                         for instance:  ['Fx_[N/m]', 'Vx_[m/s]', 'AxInd_[-]', Phi_[deg]'] 
        """
        columns=['Time_[s]']
        columns+=['Thrust_[N]']
        columns+=['Torque_[N/m]']

        #if not hasattr(self,'CP'):
        self.calcOutput()

        df = pd.DataFrame()
        df['Time_[s]']        = self.time
        df['Azimuth_[deg]']   = np.mod(self.psi,360)
        df['RtAeroFxh_[N]']   = self.Thrust
        df['RtAeroMxh_[N-m]'] = self.Torque
        df['RtAeroPwr_[W]']   = self.Power
        df['RtAeroCt_[-]']    = self.CT
        df['RtAeroCq_[-]']    = self.CQ
        df['RtAeroCp_[-]']    = self.CP

        df['RtVAvgxh_[m/s]']  = self.RtVAvg[:,0]
        df['RtVAvgyh_[m/s]']  = self.RtVAvg[:,1]
        df['RtVAvgzh_[m/s]']  = self.RtVAvg[:,2]
        df['RtArea_[m^2]']    = self.RtArea
        df['RtSkew_[deg]']    = self.chi0



        for iB in np.arange(self.nB):
            df['B'+str(iB+1)+'Azimuth_[deg]']  = np.mod(self.SkewAzimuth[:,iB],360)

        Vflw_s = self.Vwnd_s-self.Vstr_s


        # --- Creating all the possible column names for all blades and radial position, matching OpenFAST convention
        if BldNd_BladesOut is None: 
            BldNd_BladesOut=np.arange(self.nB)
        else:
            BldNd_BladesOut=np.arange(min(BldNd_BladesOut, self.nB))
        if BldNd_BlOutNd is None: 
            BldNd_BlOutNd = np.arange(len(self.r))
        if BldNd_OutList is None:
            BldNd_OutList=['Fx_[N/m]','Fy_[N/m]','Vx_[m/s]','Vy_[m/s]','VDisx_[m/s]','VDisy_[m/s]','STVx_[m/s]','STVy_[m/s]','STVz_[m/s]','Vrel_[m/s]','TnInd_[-]','AxInd_[-]','Phi_[deg]','Vindx_[m/s]','Vindy_[m/s]','Alpha_[deg]','Fn_[N/m]','Ft_[N/m]','Cl_[-]','Cd_[-]']
        # All columns
        BldNd_columns = ['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+ c for iB in BldNd_BladesOut  for c in BldNd_OutList for ir in BldNd_BlOutNd]
        # Dataframe
        df_B_r = pd.DataFrame(np.zeros((len(self.time), len(BldNd_columns))), columns=BldNd_columns)

        df= pd.concat((df,df_B_r),axis=1)

        # AeroDyn x-y is "section coord" s
        # AeroDyn n-t is "airfoil coord" a
        # AeroDyn doesn't have polar coord..
        for iB in np.arange(self.nB):
            if 'Fx_[N/m]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Fx_[N/m]'] = self.F_s[:,iB,ir,0]
            if 'Fy_[N/m]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Fy_[N/m]'] =-self.F_s[:,iB,ir,1] # NOTE: weird sign
            if 'Vx_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vx_[m/s]'] =      Vflw_s[:,iB,ir,0]
            if 'Vy_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vy_[m/s]'] =      Vflw_s[:,iB,ir,1]
            if 'VDisx_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'VDisx_[m/s]'] = self.Vwnd_s[:,iB,ir,0]
            if 'VDisy_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'VDisy_[m/s]'] = self.Vwnd_s[:,iB,ir,1]
            if 'STVx_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'STVx_[m/s]'] = self.Vstr_s[:,iB,ir,0]
            if 'STVy_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'STVy_[m/s]'] = self.Vstr_s[:,iB,ir,1]
            if 'STVz_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'STVz_[m/s]'] = self.Vstr_s[:,iB,ir,2]
            if 'Vrel_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vrel_[m/s]'] = self.Vrel[:,iB,ir]
            if 'TnInd_[-]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'TnInd_[-]'] = self.TnInd[:,iB,ir]
            if 'AxInd_[-]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'AxInd_[-]'] = self.AxInd[:,iB,ir]
            if 'Phi_[deg]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Phi_[deg]'] = self.phi[:,iB,ir]
            if 'Vindx_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vindx_[m/s]'] = self.Vind_s[:,iB,ir,0]
            if 'Vindy_[m/s]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vindy_[m/s]'] = self.Vind_s[:,iB,ir,1]
            if 'Alpha_[deg]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Alpha_[deg]'] = self.alpha[:,iB,ir]
            #AeroDyn "n-t", is almost like xa but y is switched
            if 'Fn_[N/m]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Fn_[N/m]'] = self.F_a[:,iB,ir,0]
            if 'Ft_[N/m]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Ft_[N/m]'] =-self.F_a[:,iB,ir,1]
            if 'Cl_[-]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Cl_[-]'] = self.Cl[:,iB,ir]
            if 'Cd_[-]' in BldNd_OutList:
                for ir in np.arange(len(self.r)):
                    df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Cd_[-]'] = self.Cd[:,iB,ir]
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
            R_ntr2g, R_bld2r, # "polar grid to global" for each blade
            pos_gl, Vstr_gl, R_s2g, R_a2g,            # Kinematics of nodes
            Vwnd_gl, # Wind at each positions in global
            firstCallEquilibrium=False
            ):
        """ 
        xBEM0: BEM states at t-1
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
        # --- Step 0: geometry 
        # --------------------------------------------------------------------------------
        # --- Compute rotor radius, hub radius, and section radii
        r = np.zeros((nB,nr))
        Rs=[0]*nB
        for iB in np.arange(nB):
            R_p2g = R_ntr2g[iB]
            # radial position (in polar grid) of first and last node taken
            Rs[iB] = (R_p2g.T).dot(pos_gl[iB,-1,:]-origin_pos_gl)[2]
            rhub   = (R_p2g.T).dot(pos_gl[iB,0,:] -origin_pos_gl)[2]
            # loop on elements
            for ie in np.arange(nr):
                r[iB,ie] = (R_p2g.T).dot(pos_gl[iB,ie,:]-origin_pos_gl)[2] # radial position in polar grid
        R = np.max(Rs)
        # --- Rotor speed for power
        omega_r = R_r2g.T.dot(omega_gl) # rotational speed in rotor coordinate system
        Omega = omega_r[0] # rotation speed of shaft (along x)

        if firstCallEquilibrium:
            nit=50
        else:
            nit=1
        for iterations in np.arange(nit):
            # --------------------------------------------------------------------------------
            # --- Step 1: velocity components
            # --------------------------------------------------------------------------------
            Vrel_a  = np.zeros((nB,nr,3))
            Vrel_p  = np.zeros((nB,nr,3))
            Vstr_p  = np.zeros((nB,nr,3))
            Vwnd_p  = np.zeros((nB,nr,3))
            for iB in np.arange(nB):
                R_p2g = R_ntr2g[iB]
                for ie in np.arange(nr):
                    # Velocity in global
                    Vwnd_g = Vwnd_gl[iB,ie]
                    # NOTE: inductions from previous time step, in polar grid (more realistic than global)
                    #Vind_g = xd0.Vind_g[iB,ie] # dynamic inductions at previous time step
                    Vind_g = (R_p2g).dot(xd0.Vind_p[iB,ie]) # dynamic inductions at previous time step
                    Vstr_g = Vstr_gl[iB,ie]
                    Vrel_g = Vwnd_g+Vind_g-Vstr_g
                    # Airfoil coordinates
                    Vrel_a[iB,ie] = (R_a2g[iB,ie].T).dot(Vrel_g)
                    # Polar coordinates
                    Vstr_p[iB,ie] = (R_p2g.T).dot(Vstr_g) # Structural velocity in polar coordinates
                    Vrel_p[iB,ie] = (R_p2g.T).dot(Vrel_g)
                    Vwnd_p[iB,ie] = (R_p2g.T).dot(Vwnd_g) # Wind Velocity in polar coordinates
            Vflw_p  = Vwnd_p-Vstr_p # Relative flow velocity, including wind and structural motion
            Vflw_g  = Vwnd_gl-Vstr_gl # Relative flow velocity, including wind and structural motion

            # Velocity norm and Reynolds
            Vrel_norm = sqrt(Vrel_a[:,:,0]**2 + Vrel_a[:,:,1]**2)
            Re        = Vrel_norm*p.chord/p.kinVisc/10**6 # Reynolds in million
            # --------------------------------------------------------------------------------
            # --- Step 2: Flow Angle and tip loss
            # --------------------------------------------------------------------------------
            phi_p = np.arctan2(Vrel_p[:,:,0],-Vrel_p[:,:,1])  # NOTE: using polar grid for phi
            # --- Tip loss
            F = np.ones((nB,nr))
            if (p.bTipLoss): #Glauert tip correction
                b=sin(phi_p)>0.01
                F[b] = 2./pi*arccos(exp(-(nB *(R-r[b]))/(2*r[b]*sin(phi_p[b]))))
                b2=abs(r-R)<1e-3
                F[b2]=0.001
            # --- Hub loss
            if (p.bHubLoss): #Glauert hub loss correction
                F = F* 2./pi*arccos(exp(-nB/2. *(r-rhub)/ (rhub*np.sin(phi_p))))
            #F[F<=1e-3]=0.5
            # --------------------------------------------------------------------------------
            # --- Step 3: Angle of attack
            # --------------------------------------------------------------------------------
            alpha = np.arctan2(Vrel_a[:,:,0],Vrel_a[:,:,1])        # angle of attack [rad]
            # --------------------------------------------------------------------------------
            # --- Step 4: Aerodynamic Coefficients
            # --------------------------------------------------------------------------------
            ClCdCm = np.array([p.fPolars[ie](alpha[iB,ie]) for iB in np.arange(nB) for ie in np.arange(nr)]).reshape((nB,nr,3))
            Cl=ClCdCm[:,:,0]
            Cd=ClCdCm[:,:,1]
            # Project to airfoil coordinates
            C_xa       ,C_ya        = Cl*cos(alpha)+ Cd*sin(alpha  )   ,  -Cl*sin(alpha)+ Cd*cos(alpha)
            C_xa_noDrag,C_ya_noDrag = Cl*cos(alpha)                    ,  -Cl*sin(alpha)
            # Project to polar coordinates
            C_p        = np.zeros((nB,nr,3))
            C_g        = np.zeros((nB,nr,3))
            C_p_noDrag = np.zeros((nB,nr,3))
            for iB in np.arange(nB):
                R_p2g = R_ntr2g[iB]
                for ie in np.arange(nr):
                    C_g        [iB,ie]=R_a2g[iB,ie].dot(np.array([C_xa       [iB,ie], C_ya       [iB,ie], 0]))
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
            sigma    = p.chord*p.nB/(2*pi*r)
            #a,aprime,CT = fInductionCoefficients(a_last,Vrel_in4,Un,Ut,V0_in3,V0_in4,nnW_in4,omega,chord(e),F,Ftip,CnForAI,CtForTI,lambda_r,sigma(e),phi,Algo)
            if p.WakeMod==0:
                a      = V0*0
                aprime = V0*0
            else:
                a,aprime,CT = _fInductionCoefficients(Vrel_norm, V0, F, cnForAI, ctForTI, lambda_r, sigma, phi_p, 
                        bSwirl=p.bSwirl, CTcorrection=p.CTcorrection, swirlMethod=p.swirlMethod,
                        relaxation=p.relaxation, a_last=xd0.a
                )

            if np.any(np.isnan(a)):
                print('>> BEM crashing')

            # Storing last values, for relaxation
            xd1.a=a.copy()
            # Quasi steady inductions, polar and global coordinates
            xd1.Vind_qs_p = np.zeros((nB,nr,3))
            for iB in np.arange(nB):
                R_p2g = R_ntr2g[iB]
                for ie in np.arange(nr):
                    # NOTE: Vind is negative along n and t!
                    xd1.Vind_qs_p[iB,ie] = np.array([-a[iB,ie]*Vflw_p[iB,ie,0],  aprime[iB,ie]*Vflw_p[iB,ie,1], 0])
                    xd1.Vind_qs_g[iB,ie] = R_p2g.dot(xd1.Vind_qs_p[iB,ie]) # global

            if firstCallEquilibrium:
                # We update the previous states induction
                xd0.a      = a.copy()
                xd0.Vind_g = xd1.Vind_qs_g.copy()
                xd0.Vind_p = xd1.Vind_qs_p.copy()
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
            a_avg = min([np.mean(a),0.5])
            V_avg = max([np.mean(V0),0.001])
            tau1 = 1.1 / (1 - 1.3 *a_avg)*R/V_avg
            tau2 = (0.39 - 0.26 * (r/R)**2) * tau1
            tau2 = np.tile(tau2[:,:,None],3)
            # Oye's dynamic inflow model, discrete time integration
            H              = xd1.Vind_qs_p + 0.6 * tau1 * (xd1.Vind_qs_p - xd0.Vind_qs_p) /dt
            xd1.Vind_int_p = H + (xd0.Vind_int_p - H) * exp(-dt/tau1) # intermediate velocity
            xd1.Vind_dyn_p = xd1.Vind_int_p + (xd0.Vind_dyn_p - xd1.Vind_int_p) * exp(-dt/tau2)
            # In global
            for iB in np.arange(nB):
                R_p2g = R_ntr2g[iB]
                for ie in np.arange(nr):
                    xd1.Vind_dyn_g[iB,ie] = R_p2g.dot(xd1.Vind_dyn_p[iB,ie]) # global
        else:
            xd1.Vind_dyn_g = xd1.Vind_qs_g.copy()
            xd1.Vind_dyn_p = xd1.Vind_qs_p.copy()

        # --------------------------------------------------------------------------------}
        # --- Disk averaged quantities
        # --------------------------------------------------------------------------------{
        # Average wind in global, and rotor coord
        Vwnd_avg_g = np.mean(np.mean(Vwnd_gl[:,:,:],axis=0),axis=0)
        Vwnd_avg_r = (R_r2g.T).dot(Vwnd_avg_g)
        # Average relative wind (Wnd-Str)
        Vflw_avg_g = np.mean(np.mean(Vflw_g[:,:,:],axis=0),axis=0)
        x_hat_disk = R_r2g[:,0]
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
        # Fake "Azimuth angle" used for skew model
        SkewAzimuth=np.zeros(nB)
        for iB in np.arange(nB):
            z_hat = R_ntr2g[iB][:,2]
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

        # --------------------------------------------------------------------------------
        # ---  Yaw model, repartition of the induced velocity
        # --------------------------------------------------------------------------------
        if p.bYawModel:
           #xd1.Vind_g = xd1.Vind_dyn_g.copy()
           #xd1.Vind_p = xd1.Vind_dyn_p.copy()
           #psi0 = np.arctan( Vwnd_avg_g[2]/Vwnd_avg_r[1])  # TODO
           # Sections that are about 0.7%R
           Ir= np.logical_and(r[0]>=0.5*R, r[0] <=0.8*R)
           if len(Ir)==0:
               Ir=r[0]>0
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
               R_p2g = R_ntr2g[iB]
               for ie in np.arange(nr):
                   xd1.Vind_p[iB,ie] = xd1.Vind_dyn_p[iB,ie]
                   xd1.Vind_p[iB,ie,0] = xd1.Vind_dyn_p[iB,ie,0] * (1 + yawCorrFactor*r[iB,ie]/R * np.tan(chi/2)*np.sin(SkewAzimuth[iB])) #* np.cos(psiB0[iB]+psi - psi0))
                   xd1.Vind_g[iB,ie] = R_p2g.dot(xd1.Vind_p[iB,ie]) # global
                   # AeroDyn:
                   #chi = (0.6_ReKi*a + 1.0_ReKi)*chi0
                   #a = a * (1.0 +  yawCorrFactor * yawCorr_tan * (tipRatio) * sin(azimuth))


        else:
           xd1.Vind_g = xd1.Vind_dyn_g.copy()
           xd1.Vind_p = xd1.Vind_dyn_p.copy()
        # --------------------------------------------------------------------------------
        # --- Step 6: Outputs
        # --------------------------------------------------------------------------------
        it = xd1.it # time step
        # --- Coefficients
        self.Cl[it]   = Cl
        self.Cd[it]   = Cd
        self.cn[it]   = C_p[:,:,0]
        self.ct[it]   = C_p[:,:,1]
        # C_g also available
        # --- Loads
        q_dyn = 0.5 * p.rho * Vrel_norm**2 * p.chord # dynamic pressure
        self.L[it]    = q_dyn * Cl
        self.D[it]    = q_dyn * Cd
        self.Fn[it]   = q_dyn * C_p[:,:,0]
        self.Ft[it]   = q_dyn * C_p[:,:,1]
        self.F_a[it,:,:,0] = q_dyn * C_xa
        self.F_a[it,:,:,1] = q_dyn * C_ya
        # --- Velocities
        a_dyn      =-xd1.Vind_p[:,:,0]/Vflw_p[:,:,0] 
        aprime_dyn = xd1.Vind_p[:,:,1]/Vflw_p[:,:,1]
        self.AxInd[it] = a_dyn      
        self.TnInd[it] = aprime_dyn 
        self.Vrel[it]  = Vrel_norm
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
        # --- Misc
        self.alpha[it] = alpha*180./pi
        self.phi[it]   = phi_p*180./pi
        self.Gamma[it]  = 0.5*Re*Cl*p.kinVisc*10**6 # Circulation [m^2/s]
        self.psi[it]  = psi*180/pi
        self.RtArea[it]  = pi*R**2

        for iB in np.arange(nB):
            R_p2g = R_ntr2g[iB]
            for ie in np.arange(nr):
                Vind_g = xd1.Vind_g[iB,ie] # dynamic inductions at current time step
                Vind_s = (R_s2g[iB,ie].T).dot(Vind_g) # Induced velocity in section coordinates
                #Vind_a = (R_a2g[iB,ie].T).dot(Vind_g) # Induced velocity in airfoil coordinates
                Vind_p = (R_p2g       .T).dot(Vind_g) # Induced velocity in polar coordinates
                self.Vind_p[it,iB,ie,:] =Vind_p
                self.Vind_s[it,iB,ie,:] =Vind_s
                ## Wind
                Vwnd_g = Vwnd_gl[iB,ie]
                Vwnd_s = (R_s2g[iB,ie].T).dot(Vwnd_g) # Wind Velocity in section coordinates
                #Vwnd_a = (R_a2g[iB,ie].T).dot(Vwnd_gl) # Wind Velocity in airfoil coordinates
                self.Vwnd_s[it,iB,ie,:]=Vwnd_s
                ## Structural velocity
                Vstr_g = Vstr_gl[iB,ie]
                #Vstr_a = (R_a2g[iB,ie].T).dot(Vstr_g) # Structural velocity in airfoil coordinates
                self.Vstr_s[it,iB,ie,:]=(R_s2g[iB,ie].T).dot(Vstr_g)
                # --- Loads
                F_g = q_dyn[iB,ie] * C_g[iB,ie]
                self.F_s[it,iB,ie] = (R_s2g[iB,ie].T).dot(F_g)


        # Blade integrated loads
        self.BladeThrust[it,:] = np.trapz(self.Fn[it,:]  , r) # Normal to rotor plane
        self.BladeTorque[it,:] = np.trapz(self.Ft[it,:]*r, r) # About shaft 
        self.Thrust[it] = np.sum(self.BladeThrust[it,:])    # Normal to rotor plane
        self.Torque[it] = np.sum(self.BladeTorque[it,:])
        self.Power[it]  = Omega*self.Torque[it]
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
                # #RES.BladeFlap(idB)=trapz([Rotor.r;R],[Rotor.r.*(Pz*0+Pn);0]);
                # #RES.BladeEdge(idB)=trapz([Rotor.r;R],[Rotor.r.*(Py*0+Pt);0]);
                #  #### Returning Aerodynamic Forces
                #  RES.Flap = sum(RES.BladeFlap)
                #  RES.Edge = sum(RES.BladeEdge)
        return xd1

    def simulationConstantRPM(self, time, RPM, windSpeed=None, windExponent=None, windRefH=None, windFunction=None, cone=0, tilt=0, hubHeight=None, firstCallEquilibrium=True):
        """ 
        wrapper function to perform a simple simulation at constant RPM
       
        IINPUTS:
        Different ways to specify wind:
          - windSpeed: scalar, wind speed (at hub height) along x
          - windExponent power law exponent speed for wind speed (None=uniform wind)
          - windRefH reference height for power law
        OR
          - windFunction: function with interface: f(x,y,z,t)=u,v,w
                       with x,y,z,u arrays of arbitrary shapes

        - cone: override cone values, in degrees, OpenFAST convention
        - tilt: override tilt values, in degrees, OpenFAST convention
        - hubHeight: if provided, override the hubheight that is computed based on OverHang, Twr2Shaft and Tilt
        - firstCallEquilibrium: if true, the inductions are set to the equilibrium values at t=0 (otherwise,0)

        """
        # --- Define motion
        motion = PrescribedRotorMotion()
        motion.init_from_BEM(self, tilt=tilt, cone=cone, psi0=0)
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
                    motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                    motion.R_ntr2g, motion.R_bld2b,
                    motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g,
                    Vwnd_g,
                    firstCallEquilibrium= it==0 and firstCallEquilibrium
                    )
            #if np.mod(t,1)<dt/2:
            #    print(t)
        df = self.toDataFrame()
        return df

# --------------------------------------------------------------------------------}
# --- Helper class to prescribe a motion
# --------------------------------------------------------------------------------{
from welib.yams.utils import R_x, R_y, R_z
class PrescribedRotorMotion():
    """ 
    Class to return:
        position, velocity, and orientation of all blade station

     - R_a2g : from airfoil to global (this is well defined, also called "n-t" system in AeroDyn)
     - R_s2g : from section to blobal (this is ill-defined), this coordinate is used to define the "axial" and "tangential" inductions

    """
    def __init__(self):
        # body kinematics
        self.R_b2g = np.eye(3) # orientation from body to global
        self.R_b2g0 = np.eye(3) # orientation from body to global at t=0
        self.origin_pos_gl0 = np.array([0,0,0]) # body origin at t=0
        self.origin_pos_gl = np.array([0,0,0]) # body origin
        self.origin_vel_gl = np.array([0,0,0])
        self.omega_gl      = np.array([0,0,0])

        # Blades 
        self.R_bld2b=None # rotation matrix from blade to body

    def init_from_inputs(self,  nB, r, twist, rotorOrigin, tilt, cone, psi0=0):
        self.nB   =  nB
        self.cone = cone*np.pi/180
        self.tilt = -tilt*np.pi/180 
        self.r     = r
        self.twist = twist
        self.allocate()
        self.origin_pos_gl0 = rotorOrigin

        # Basic geometries for nacelle
        self.R_b2g0 = R_y(self.tilt)  # Rotation fromShaft to Nacelle

        # Orientation of blades with respect to body
        self.R_bld2b = [np.eye(3)]*self.nB
        self.R_ntr2b = [np.eye(3)]*self.nB
        self.psi_B0  = np.zeros(self.nB)
        for iB in np.arange(self.nB):
            self.psi_B0[iB]= psi0 + iB*2*np.pi/self.nB
            R_SB = R_x(self.psi_B0[iB]) 
            self.R_bld2b[iB] = R_SB.dot(R_y(self.cone)) # blade2shaft
            self.R_ntr2b[iB] = R_SB.dot(np.array([[1,0,0],[0,-1,0],[0,0,1]])) # "n-t-r" to shaft

        # Set initial positions and orientations in body coordinates
        nr=len(self.r)
        meanLine=np.zeros((nr,3))
        meanLine[:,0]=0 # TODO
        meanLine[:,1]=0 # TODO
        meanLine[:,2]=self.r
        for iB in np.arange(self.nB):
            for ir in np.arange(nr):
                self.pos0[iB,ir,:] = self.R_bld2b[iB].dot(meanLine[ir,:])
                R_a2bld = R_z(-self.twist[ir]) # section to blade
                self.R_a02b[iB,ir,:,:] = self.R_bld2b[iB].dot(R_a2bld)   # TODO curvature
                self.R_s02b[iB,ir,:,:] = self.R_bld2b[iB]                # TODO curvature

    def init_from_BEM(self, BEM, tilt=None, cone=None, psi0=0):
        """ 
        Initializes motion from a BEM class
        Possibility to override the tilt and cone geometry:
        tilt: tilt angle in deg, with OpenFAST convention
        cone: cone angle in deg, with OpenFAST convention
        
        """
        if tilt is None:
            tilt=BEM.tilt0
        if cone is None:
            cone=BEM.cone0

        rotorOrigin =[BEM.OverHang*np.cos(-tilt*np.pi/180), 0, BEM.TowerHt+BEM.Twr2Shft-BEM.OverHang*np.sin(-tilt*np.pi/180)]
        self.init_from_inputs(BEM.nB, BEM.r, BEM.twist, rotorOrigin, tilt=tilt, cone=cone, psi0=0)

    def allocate(self):
        nr = len(self.r)
        # Section nodes
        self.pos0    = np.zeros((self.nB,nr,3))   # position of nodes at t= 0 in body coordinates 
        self.R_s02b  = np.zeros((self.nB,nr,3,3)) # Orientation section to body at t=0
        self.R_a02b  = np.zeros((self.nB,nr,3,3)) # Orientation airfoil to body at t=0
        self.pos_gl = np.zeros((self.nB,nr,3))   # position of all nodes
        self.vel_gl = np.zeros((self.nB,nr,3))   # linear velocities
        self.R_s2g  = np.zeros((self.nB,nr,3,3)) # Orientation section to global
        self.R_a2g  = np.zeros((self.nB,nr,3,3)) # Orientation airfoil to global
        self.R_ntr2g = [np.eye(3)]*self.nB

    def setType(self, sType, **kwargs):
        self.sType=sType
        self.opts=kwargs

    def rigidbodyKinUpdate(self, P_gl, vel_gl, omega_gl, R_b2g):
        """
        Update position, velocities, orientations of nodes assuming a rigid body motion
        of the origin given as input
        """
        self.origin_pos_gl = P_gl
        self.origin_vel_gl = vel_gl
        self.omega_gl      = omega_gl
        self.R_b2g         = R_b2g

        # Update of positions
        for iB in np.arange(self.nB):
            self.R_ntr2g[iB] = R_b2g.dot(self.R_ntr2b[iB])
            for ir in np.arange(len(self.r)):
                s_OP =  R_b2g.dot(self.pos0[iB,ir,:])
                self.pos_gl[iB,ir,:] = P_gl   + s_OP
                self.vel_gl[iB,ir,:] = vel_gl + np.cross(omega_gl, s_OP)
                self.R_s2g[iB,ir,:,:] = R_b2g.dot(self.R_s02b[iB,ir,:])
                self.R_a2g[iB,ir,:,:] = R_b2g.dot(self.R_a02b[iB,ir,:])

    def update(self, t):
        if self.sType=='constantRPM':
            omega = self.opts['RPM']*2.*np.pi/(60.)
            psi = t*omega
            pos = self.origin_pos_gl0
            vel = np.array([0,0,0])
            ome = self.R_b2g0.dot(np.array([omega,0,0]))
            R_b2g = self.R_b2g0.dot(R_x(psi))
            self.rigidbodyKinUpdate(pos, vel, ome, R_b2g)
            self.psi=psi # hack

        elif self.sType=='x-oscillation':
            omega = self.opts['frequency']*2.*np.pi
            A = self.opts['amplitude']
            x    = A*np.sin(omega*t)
            xdot = A*omega*np.cos(omega*t)
            pos = self.origin_pos_gl0+np.array([x,0,0])
            vel = np.array([xdot,0,0])
            ome = np.array([0,0,0])
            R_b2g =self.R_b2g0
            self.rigidbodyKinUpdate(pos, vel, ome, R_b2g)

        elif self.sType=='constantRPM x-oscillation':
            omega = self.opts['frequency']*2.*np.pi
            A = self.opts['amplitude']
            x    = A*np.sin(omega*t)
            xdot = A*omega*np.cos(omega*t)
            omegapsi = self.opts['RPM']*2.*np.pi/(60.)
            psi = t*omegapsi
            pos = self.origin_pos_gl0+np.array([x,0,0])
            vel = np.array([xdot,0,0])
            ome = self.R_b2g0.dot(np.array([omegapsi,0,0]))
            R_b2g = self.R_b2g0.dot(R_x(psi))
            self.rigidbodyKinUpdate(pos, vel, ome, R_b2g)
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

if __name__=="__main__":
    """ See examples/ for more examples """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.animation import FuncAnimation

    # --- Read a FAST model to get Aerodynamic parameters
    BEM = AeroBEM()
    #BEM.init_from_FAST('../../data/NREL5MW/Main_Onshore.fst')
    BEM.init_from_FAST('./Main_Onshore.fst')
    BEM.CTcorrection='AeroDyn' # High Thrust correction
    BEM.swirlMethod ='AeroDyn' # type of swirl model
#    BEM.swirlMethod ='HAWC2' # type of swirl model
#     BEM.bSwirl = True  # swirl flow model enabled / disabled
    BEM.WakeMod=1 # 0: no inductions, 1: BEM inductions
#     BEM.bTipLoss = True # enable / disable tip loss model
#     BEM.bHubLoss = False # enable / disable hub loss model
#     BEM.bTipLossCl = False # enable / disable Cl loss model
#     BEM.TipLossMethod = 'Glauert'  # type of tip loss model
#     BEM.bDynaStall = True # dynamic stall model
    BEM.bDynaWake = True # dynamic inflow model
#    BEM.bDynaWake = True # dynamic inflow model
    BEM.bYawModel = True # Yaw correction
#     BEM.bYawModel = False # Yaw correction
#     BEM.bAIDrag = True # influence on drag coefficient on normal force coefficient
#     BEM.bTIDrag = True # influence on drag coefficient on tangential force coefficient
    BEM.relaxation = 0.3

    time=np.arange(0,10,0.05)
    RPM=10
    df= BEM.simulationConstantRPM(time, RPM, windSpeed=10, windExponent=0.0, windRefH=90, windFunction=None, cone=None, tilt=None, firstCallEquilibrium=True)
    df.to_csv('_BEM.csv', index=False)

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
#                 motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
#                 motion.R_ntr2g, motion.R_bld2b,
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

