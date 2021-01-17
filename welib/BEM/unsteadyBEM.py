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
        self.bYawModel = True # Yaw correction
        self.bAIDrag = True # influence on drag coefficient on normal force coefficient
        self.bTIDrag = True # influence on drag coefficient on tangential force coefficient
        self.bReInterp = False # interpolate the input tabulated airfoil data for Reynolds variation
        self.bThicknessInterp = True # interpolate the input tabulated airfoil data for thickness variation
        self.WakeMod=1 # 0: no inductions, 1: BEM inductions
        self.bRoughProfiles = False # use rough profiles for input airfoil data

    def init_from_FAST(self, FASTFileName):
        import welib.weio as weio
        F=weio.FASTInputDeck(FASTFileName,readlist=['AD','ED'])

        # Environment
        self.rho     = F.AD['AirDens']
        self.kinVisc = F.AD['KinVisc']

        self.nB   =  F.ED['NumBl']
        self.r     = F.AD.Bld1['BldAeroNodes'][:,0] + F.ED['HubRad']
        chord = F.AD.Bld1['BldAeroNodes'][:,-2] 
        self.chord=np.stack([chord]*self.nB)
        #self.twist = F.AD.Bld1['BldAeroNodes'][:,-3]*nself.pi/180
        #self.cone = -F.ED['PreCone(1)']*nself.pi/180
        polars=[]
        ProfileID=F.AD.Bld1['BldAeroNodes'][:,-1].astype(int)
        for ipolar in  ProfileID:
            polars.append(F.AD.AF[ipolar-1]['AFCoeff'])
        self.polars = polars

        self._init()

    def _init(self):
        # Creating interpolation functions for each polar, now in rad!
        self.fPolars = [interp1d(p[:,0]*np.pi/180,p[:,1:],axis=0) for p in self.polars]

    def getInitStates(self):
        return BEMDiscreteStates(self.nB, len(self.r))

    def timeStepInit(self, t0, tmax, dt):
        """ Allocate storage for tiem step values"""
        self.time=np.arange(t0,tmax,dt)
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
        self.Vrel_n = np.zeros((nt,nB,nr)) # Un
        self.Vrel_t = np.zeros((nt,nB,nr)) # Ut
        self.Vrel_r = np.zeros((nt,nB,nr)) # Ur
        self.Vrel_xa = np.zeros((nt,nB,nr)) # 
        self.Vrel_ya = np.zeros((nt,nB,nr)) # 
        self.Vrel_za = np.zeros((nt,nB,nr)) # 
        self.Vind_n = np.zeros((nt,nB,nr))
        self.Vind_t = np.zeros((nt,nB,nr))
        self.Vind_qs_n = np.zeros((nt,nB,nr))
        self.Vind_qs_t = np.zeros((nt,nB,nr))
        self.Vwnd_n = np.zeros((nt,nB,nr))
        self.Vwnd_t = np.zeros((nt,nB,nr))
        self.Vwnd_r = np.zeros((nt,nB,nr))
        self.Vwnd_xs = np.zeros((nt,nB,nr))
        self.Vwnd_ys = np.zeros((nt,nB,nr))
        self.Vwnd_xa = np.zeros((nt,nB,nr))
        self.Vwnd_ya = np.zeros((nt,nB,nr))
        self.Vstr_n = np.zeros((nt,nB,nr))
        self.Vstr_t = np.zeros((nt,nB,nr))
        self.Vstr_r = np.zeros((nt,nB,nr))
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
        self.Fx_a   = np.zeros((nt,nB,nr))
        self.Fy_a   = np.zeros((nt,nB,nr))
        self.Gamma  = np.zeros((nt,nB,nr))
        self.alpha  = np.zeros((nt,nB,nr))
        self.phi    = np.zeros((nt,nB,nr))
        self.Re     = np.zeros((nt,nB,nr))
        # Advanced
        self.Cx_s   = np.zeros((nt,nB,nr))
        self.Cy_s   = np.zeros((nt,nB,nr))
        self.Vstr_xs = np.zeros((nt,nB,nr))
        self.Vstr_ys = np.zeros((nt,nB,nr))
        self.Px_s   = np.zeros((nt,nB,nr))
        self.Py_s   = np.zeros((nt,nB,nr))
        # Integrated values
        self.Thrust   = np.zeros(nt)
        self.Torque   = np.zeros(nt)
        self.Power    = np.zeros(nt)
        self.chi      = np.zeros(nt)
        self.RtVAvgxh   = np.zeros(nt)
        self.RtVAvgyh   = np.zeros(nt)
        self.RtVAvgzh   = np.zeros(nt)
        self.psi      = np.zeros(nt)
        self.RtArea   = np.zeros(nt)
        # Blade blades
        self.BladeTorque = np.zeros((nt,nB))
        self.BladeThrust = np.zeros((nt,nB))
        self.BladeEdge   = np.zeros((nt,nB))
        self.BladeFlap   = np.zeros((nt,nB))

    def toDataFrame(self):
        """ Export time series to a pandas dataframe
        Column names are set to match OpenFAST outputs
        """
        columns=['Time_[s]']
        columns+=['Thrust_[N]']
        columns+=['Torque_[N/m]']

        df = pd.DataFrame()
        df['Time_[s]']        = self.time
        df['Azimuth_[deg]']   = np.mod(self.psi,360)
        df['RtAeroFxh_[N]']   = self.Thrust
        df['RtAeroMxh_[N-m]'] = self.Torque
        df['RtVAvgxh_[m/s]']  = self.RtVAvgxh
        df['RtArea_[m^2]']    = self.RtArea
        # AeroDyn Fx/ Fy is "section coord" s, but we use polar here
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Fx_[N/m]'] = self.Fn[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Fy_[N/m]'] = self.Ft[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'STVx_[m/s]'] = self.Vstr_n[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'STVy_[m/s]'] = - self.Vstr_t[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'STVz_[m/s]'] = self.Vstr_r[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vrel_[m/s]'] = self.Vrel[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'VDisx_[m/s]'] = self.Vwnd_n[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'VDisy_[m/s]'] =-self.Vwnd_t[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'TnInd_[-]'] = self.TnInd[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'AxInd_[-]'] = self.AxInd[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Phi_[deg]'] = self.phi[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vindx_[m/s]'] = self.Vind_n[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Vindy_[m/s]'] =-self.Vind_t[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Alpha_[deg]'] = self.alpha[:,iB,ir]
        # AeroDyn "n-t", is almost like xa but y is switched
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Fn_[N/m]'] = self.Fx_a[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Ft_[N/m]'] =-self.Fy_a[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Cl_[-]'] = self.Cl[:,iB,ir]
        for iB in np.arange(self.nB):
            for ir in np.arange(len(self.r)):
                df['AB'+str(iB+1)+'N{:03d}'.format(ir+1)+'Cd_[-]'] = self.Cd[:,iB,ir]
        return df

    def toDataFrameRadial(self, it):
        df = pd.DataFrame()
        df['r/R_[-]']   = self.r/max(self.r)
        #df['r_[m]']     = self.r
        for iB in np.arange(self.nB):
            B='B'+str(iB+1)
            df[B+'AxInd_[-]'] = self.AxInd[it,iB,:]
            df[B+'TnInd_[-]'] = self.AxInd[it,iB,:]
        return df


    def timeStep(self, t, dt, xd0, psi,
            origin_pos_gl, omega_gl, R_b2g,  # Kinematics of rotor origin
            R_ntr2g, R_bld2b, # "polar grid to global" for each blade
            pos_gl, vel_gl, R_s2g, R_a2g,            # Kinematics of nodes
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

# #     dr = WT.Rotor.dr
# # # normal vectors
# #     n_thrust_in3 = np.array([[0],[0],[- 1]])
# #     #n_rotor_in4=[0; 0; 1];         #normal vector to the rotor plan in 4
# #     n_rotor_in2 = np.array([[0],[0],[1]])
        # Time step storage for vectorization
        nB, nr, _ = pos_gl.shape
# #     # yaw model initialization
# #     khi = WT.Aero.last.chi
# #     psi0 = 0
# #     ## Elasticity dependence
# # # Angles between blades updated depending on Dynamic degree of freedom 2
# #     Vpsi0 = np.mod(np.arange(0,(360 / nB) * (nB - 1)+(360 / nB),(360 / nB)),360)
# #     Vpsi = np.mod(Vpsi0 + x(2) * 180 / pi,360)
# #     omega = v(2)
# #     # Shaft vector coordinates
# #     rs_in1 = np.transpose(a12) * WT.Shaft.rs_in2
        p = self # alias
        Vwnd_gl = np.array([10,0,0]) # TODo
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
            Rs[iB] = (R_p2g).dot(pos_gl[iB,-1,:]-origin_pos_gl)[2]
            rhub = (R_p2g).dot(pos_gl[iB,0,:]-origin_pos_gl)[2]
            # loop on elements
            for ie in np.arange(nr):
                r[iB,ie] = (R_p2g).dot(pos_gl[iB,ie,:]-origin_pos_gl)[2] # radial position in polar grid
        R = np.max(Rs)

        if firstCallEquilibrium:
            print('>>>> EQUILIBRIUM')
            nit=200
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
                    Vwnd_g = Vwnd_gl # TODO
                    # NOTE: inductions from previous time step, in polar grid (more realistic than global)
                    #Vind_g = xd0.Vind_g[iB,ie] # dynamic inductions at previous time step
                    Vind_g = (R_p2g).dot(xd0.Vind_p[iB,ie]) # dynamic inductions at previous time step
                    Vstr_g = vel_gl[iB,ie]
                    Vrel_g = Vwnd_g+Vind_g-Vstr_g
                    # Airfoil coordinates
                    Vrel_a[iB,ie] = (R_a2g[iB,ie].T).dot(Vrel_g)
                    # Polar coordinates
                    Vstr_p[iB,ie] = (R_p2g.T).dot(Vstr_g) # Structural velocity in polar coordinates
                    Vrel_p[iB,ie] = (R_p2g.T).dot(Vrel_g)
                    Vwnd_p[iB,ie] = (R_p2g.T).dot(Vwnd_g) # Wind Velocity in polar coordinates

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
            # --- Hub loss
            if (p.bHubLoss): #Glauert hub loss correction
                F = F* 2./pi*arccos(exp(-nB/2. *(r-rhub)/ (rhub*np.sin(phi_p))))
            F[F<=1e-3]=0.5
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
            lambda_r = Vstr_p[:,:,1]/Vwnd_p[:,:,0] # "omega r/ U0n" defined in polar grid
            V0       = np.sqrt(Vwnd_p[:,:,0]**2 + Vwnd_p[:,:,1]**2) # TODO think about that
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
                    xd1.Vind_qs_p[iB,ie] = np.array([-a[iB,ie]*Vwnd_p[iB,ie,0], -aprime[iB,ie]*Vstr_p[iB,ie,1], 0])
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
            #tau1=4
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

        # --------------------------------------------------------------------------------
        # ---  Yaw model, repartition of the induced velocity
        # --------------------------------------------------------------------------------
        if p.bYawModel:
           xd1.Vind_g = xd1.Vind_dyn_g.copy() # TODO
           xd1.Vind_p = xd1.Vind_dyn_p.copy() # TODO
#                     ### Yaw model, Skew angle and psi0
#                    # if (e == Rotor.e_ref_for_khi and idB == 1):
#                    #     ### Determination of psi0
#                    #     r_hub = WT.Tower.rt_in1 + rs_in1
#                    #     # Incoming wind at hub
#                    #     V0_in1 = getPointIncomingWindLegacy02(r_hub,psi,WT,Wind,Algo)
#                    #     V0_in2 = a12 * V0_in1
#                    #     # psi0
#                    #     psi0 = atan2(V0_in2(2),V0_in2(1)) * 180 / pi
#                    #     ### Determination of skew angle
#                    #     # Averaging Wn on each blade
#                    #     meanWn_in4 = np.array([[0],[0],[mean(W0(3,Rotor.e_ref_for_khi,:))]])
#                    #     meanWn_in2 = np.transpose(a34) * meanWn_in4
#                    #     V_prime_for_khi_in2 = V0_in2 + meanWn_in2
#                    #     khi = acosd(np.dot(n_rotor_in2,V_prime_for_khi_in2) / norm(V_prime_for_khi_in2))
#                    #     W[:,e,idB] = W0(:,e,idB) * (1 + r(e) / R * tand(khi / 2) * np.cos(np.pi/180*Vpsi(idB) - psi0))
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
        self.Fx_a[it] = q_dyn * C_xa
        self.Fy_a[it] = q_dyn * C_ya
        # --- Velocities
        #Vrel_a  = np.zeros((nB,nr,3))
        self.AxInd[it] = a
        self.TnInd[it] = aprime
        self.Vrel[it]  = Vrel_norm
        # polar system (missing Vind)
        self.Vrel_n[it]  = Vrel_p[:,:,0] # NOTE: Vrel is using previous inductions..
        self.Vrel_t[it]  = Vrel_p[:,:,1]
        self.Vrel_r[it]  = Vrel_p[:,:,2]
        self.Vstr_n[it]  = Vstr_p[:,:,0]
        self.Vstr_t[it]  = Vstr_p[:,:,1]
        self.Vstr_r[it]  = Vstr_p[:,:,2]
        self.Vwnd_n[it]  = Vwnd_p[:,:,0]
        self.Vwnd_t[it]  = Vwnd_p[:,:,1]
        self.Vwnd_r[it]  = Vwnd_p[:,:,2]
        self.Vind_qs_n[it] = xd1.Vind_qs_p[:,:,0]
        self.Vind_qs_t[it] = xd1.Vind_qs_p[:,:,1]
        self.RtVAvgxh[it]  = np.mean(Vwnd_p[:,:,0])
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
                #Vind_s = (R_s2g[iB,ie].T).dot(Vind_g) # Induced velocity in section coordinates
                #Vind_a = (R_a2g[iB,ie].T).dot(Vind_g) # Induced velocity in airfoil coordinates
                Vind_p = (R_p2g       .T).dot(Vind_g) # Induced velocity in polar coordinates
                self.Vind_n[it,iB,ie] =Vind_p[0]
                self.Vind_t[it,iB,ie] =Vind_p[1]
                ## Wind
                #Vwnd_g = Vwnd_gl
                #Vwnd_s = (R_s2g[iB,ie].T).dot(Vwnd_gl) # Wind Velocity in section coordinates
                #Vwnd_a = (R_a2g[iB,ie].T).dot(Vwnd_gl) # Wind Velocity in airfoil coordinates
                ## Structural velocity
                #Vstr_g = vel_gl[iB,ie]
                #Vstr_s = (R_s2g[iB,ie].T).dot(Vstr_g) # Structural velocity in section coordinates
                #Vstr_a = (R_a2g[iB,ie].T).dot(Vstr_g) # Structural velocity in airfoil coordinates
            # Blade integrated loads
            self.BladeThrust[it] = np.trapz(self.Fn[it]  , r) # Normal to rotor plane
            self.BladeTorque[it] = np.trapz(self.Ft[it]*r, r) # About shaft 
            self.Thrust[it] = np.sum(self.BladeThrust[it])    # Normal to rotor plane
            self.Torque[it] = np.sum(self.BladeTorque[it])
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
                #  RES.Torque = sum(RES.BladeTorque)
                #  RES.Thrust = sum(RES.BladeThrust)
                #  RES.Flap = sum(RES.BladeFlap)
                #  RES.Edge = sum(RES.BladeEdge)
                #  RES.Power = omega * RES.Torque
        return xd1

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
        self.origin_pos_gl0 = np.array([0,0,0]) # body origin at t=0
        self.origin_pos_gl = np.array([0,0,0]) # body origin
        self.origin_vel_gl = np.array([0,0,0])
        self.omega_gl      = np.array([0,0,0])

        # Blades 
        self.R_bld2b=None # rotation matrix from blade to body

    def init_from_FAST(self, FASTFileName):
        import welib.weio as weio
        F=weio.FASTInputDeck(FASTFileName,readlist=['AD','ED'])

        self.nB   =  F.ED['NumBl']
        self.cone = -F.ED['PreCone(1)']*np.pi/180*0 # <<<<<<<<<<<<<<<<<<<<< HATODO TODO TODO TOO TODOTODO TODO TODO TOO TODOCK
        self.r     = F.AD.Bld1['BldAeroNodes'][:,0] + F.ED['HubRad']
        self.chord = F.AD.Bld1['BldAeroNodes'][:,-2] 
        self.twist = F.AD.Bld1['BldAeroNodes'][:,-3]*np.pi/180
        self.allocate()

        # Geometry
        ED=F.ED
        r_ET_inE    = np.array([0,0,ED['TowerBsHt']               ]) 
        r_TN_inT    = np.array([0,0,ED['TowerHt']-ED['TowerBsHt'] ])
        # Basic geometries for nacelle
        theta_tilt_y = -ED['ShftTilt']*np.pi/180*0 # <<<<<<<<<<<<<<HACK HACK TODO TODO TODO TOO TODO TODO   # NOTE: tilt has wrong orientation in FAST
        R_NS = R_y(theta_tilt_y)  # Rotation fromShaft to Nacelle
        r_NS_inN    = np.array([0             , 0, ED['Twr2Shft']]) # Shaft start in N
        r_SR_inS    = np.array([ED['OverHang'], 0, 0             ]) # Rotor center in S
        r_SGhub_inS = np.array([ED['HubCM']   , 0, 0             ]) + r_SR_inS # Hub G in S
        r_NR_inN    = r_NS_inN + R_NS.dot(r_SR_inS)                 # Rotor center in N
        r_NGnac_inN = np.array([ED['NacCMxn'],0,ED['NacCMzn']    ]) # Nacelle G in N
        r_RGhub_inS = - r_SR_inS + r_SGhub_inS

        # Orientation of blades with respect to body
        self.R_bld2b = [np.eye(3)]*self.nB
        self.R_ntr2b = [np.eye(3)]*self.nB
        for iB in np.arange(self.nB):
            psi_B= -iB*2*np.pi/self.nB
            R_SB = R_x(psi_B) 
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
            ome = np.array([omega,0,0])
            R_b2g = R_x(psi)
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
            R_b2g = np.eye(3)
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
            ome = np.array([omegapsi,0,0])
            R_b2g = R_x(psi)
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
    BEM.init_from_FAST('../../data/NREL5MW/Main_Onshore_OF2.fst')
    BEM.CTcorrection='AeroDyn' # High Thrust correction
    BEM.swirlMethod ='AeroDyn' # type of swirl model
#     BEM.bSwirl = True  # swirl flow model enabled / disabled
#     BEM.bTipLoss = True # enable / disable tip loss model
#     BEM.bHubLoss = False # enable / disable hub loss model
#     BEM.bTipLossCl = False # enable / disable Cl loss model
#     BEM.TipLossMethod = 'Glauert'  # type of tip loss model
#     BEM.bDynaStall = True # dynamic stall model
    BEM.bDynaWake = True # dynamic inflow model
#     BEM.bYawModel = True # Yaw correction
#     BEM.bAIDrag = True # influence on drag coefficient on normal force coefficient
#     BEM.bTIDrag = True # influence on drag coefficient on tangential force coefficient
    #BEM.relaxation = 1.0

    # --- Read a FAST model to get structural parameters for blade motion
    motion = PrescribedRotorMotion()
    motion.init_from_FAST('../../data/NREL5MW/Main_Onshore_OF2.fst')
    motion.setType('constantRPM', RPM=10.0)
    #motion.setType('constantRPM x-oscillation', RPM=12.1, frequency=1.1, amplitude=20)

    dt=0.1
    dtRadOut=1.0
    tmax=70
    
    xdBEM = BEM.getInitStates()

    # Allocate
    BEM.timeStepInit(0,tmax,dt) 

    for it,t in enumerate(BEM.time):
        motion.update(t)
        xdBEM = BEM.timeStep(t, dt, xdBEM, motion.psi,
                motion.origin_pos_gl, motion.omega_gl, motion.R_b2g, 
                motion.R_ntr2g, motion.R_bld2b,
                motion.pos_gl, motion.vel_gl, motion.R_s2g, motion.R_a2g,
                firstCallEquilibrium=it==0
                )
        if np.mod(t,dtRadOut)<dt/2:
            #print(t)
            #dfRad = BEM.toDataFrameRadial(it)
            #dfRad.to_csv('_BEMRad_t{:04d}.csv'.format(it), index=False)
            pass
    df = BEM.toDataFrame()
    df.to_csv('_BEM.csv', index=False)

    #ani = motion.plotAnim(0,10,0.01, ylim=[-80,80], zlim=[-80,80])
    plt.show()

