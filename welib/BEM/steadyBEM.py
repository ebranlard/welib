# -*- coding: utf-8 -*-
""" 
Reference:
     E.Branlard: Wind turbine Aerodynamics and Vorticity Based Method - Chapter 10, Springer, 2017
"""
import numpy as np
import os
from numpy import pi, cos, exp, sqrt, sin, arctan2, arccos
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------}
# --- Main Class SteadyBEM 
# --------------------------------------------------------------------------------{
class SteadyBEM():
    def __init__(self, FASTFileName=None):
        # Aero Data
        self.chord  = None
        self.polars = None
        # Environment
        self.rho     = None
        self.kinVisc = None # Kinmematic viscosity [kg/m^3]

        # Structural "Inputs"
        self.twist  = None # [deg]
        self.nB     = None # Number of blades
        self.cone0  = None # cone angle [deg]
        self.r      = None # radial stations
        # Algorithm
        self.setDefaultOptions()

        # Data from runs
        self.lastRun = None

        # Init from a FAST input file 
        if FASTFileName is not None:
            self.init_from_FAST(FASTFileName)

    def setDefaultOptions(self):
        self.nbIt = 200  # maximum number of iterations in BEM
        self.aTol = 10 ** -6 # tolerance for axial induction factor convergence
        self.relaxation = 0.5  # relaxation factor in axial induction factor
#         self.CTcorrection = 'AeroDyn'  #  type of CT correction more model implementated in the future like 'spera'
#         self.swirlMethod  = 'AeroDyn' # type of swirl model
#         self.Ngrid = 1.0
        self.bSwirl = True  # swirl flow model enabled / disabled
        self.bTipLoss = True # enable / disable tip loss model
        self.bHubLoss = False # enable / disable hub loss model
#         self.bTipLossCl = False # enable / disable Cl loss model
#         self.TipLossMethod = 'Glauert'  # type of tip loss model
        self.bAIDrag = True # influence on drag coefficient on normal force coefficient
        self.bTIDrag = True # influence on drag coefficient on tangential force coefficient
#         self.bReInterp = False # interpolate the input tabulated airfoil data for Reynolds variation
#         self.bThicknessInterp = True # interpolate the input tabulated airfoil data for thickness variation
        self.WakeMod=1 # 0: no inductions, 1: BEM inductions
#         self.bRoughProfiles = False # use rough profiles for input airfoil data

    def init_from_FAST(self, FASTFileName):
        # TODO unify with FASTFile2SteadyBEM(FASTFileName)
        from welib.weio.fast_input_deck import FASTInputDeck
        F = FASTInputDeck(FASTFileName,readlist=['AD','ED','ADbld','AF'])
        if F.AD is None:
            raise Exception('Cannot open AD file referenced in:'.format(FASTFileName))
        if F.ED is None:
            raise Exception('Cannot open ED file referenced in:'.format(FASTFileName))

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
        chord = F.AD.Bld1['BldAeroNodes'][:,-2] 
        self.chord = chord
        chord      = F.AD.Bld1['BldAeroNodes'][:,-2]
        self.chord = chord # np.stack([chord]*self.nB)
        self.twist = F.AD.Bld1['BldAeroNodes'][:,-3] # TODO unsteady is in rad! *np.pi/180
        polars=[]
        ProfileID=F.AD.Bld1['BldAeroNodes'][:,-1].astype(int)
        for ipolar in  ProfileID:
            polars.append(F.AD.AF[ipolar-1]['AFCoeff'])
        self.polars = polars

        # Input geometry (may be overriden by motion/simulations)
        self.cone0     =-F.ED['PreCone(1)'] # TODO decide on sine convention. Unsteady use opposite TODO rad or deg
        #self.tilt0     = F.ED['ShftTilt']   # TODO rad or deg
        #self.TowerHt   = F.ED['TowerHt']
        #self.OverHang  = F.ED['OverHang']
        #self.Twr2Shft  = F.ED['Twr2Shft']
        # TODO hub height maybe?

    def __repr__(self):
        s='<{} object>:\n'.format(type(self).__name__)

        # Aero Data
        s+='Aerodynamic data:\n'
        s+=' - r         : size: {}, min: {}, max: {}\n'.format(len(self.r), np.min(self.r), np.max(self.r))
        s+=' - chord     : size: {}, min: {}, max: {}\n'.format(len(self.chord), np.min(self.chord), np.max(self.chord))
        s+=' - twist     : size: {}, min: {}, max: {} [deg]\n'.format(len(self.twist), np.min(self.twist), np.max(self.twist))
        s+=' - polars    : size: {}\n'.format(len(self.polars))
        s+=' - nB        : {}\n'.format(self.nB)
        s+=' - cone0     : {} [deg]\n'.format(self.cone0 )
        # Environment
        s+='Environmental conditions:\n'
        s+=' - rho       : {} [kg/m^3]\n'.format(self.rho    )
        s+=' - kinVisc   : {} [kg/m^3]\n'.format(self.kinVisc)
        s+='Algorithm options:\n'
        s+=' - nbIt      : {}\n'.format(self.nbIt      )
        s+=' - aTol      : {}\n'.format(self.aTol      )
        s+=' - relaxation: {}\n'.format(self.relaxation)
        s+=' - bSwirl    : {}\n'.format(self.bSwirl    )
        s+=' - bTipLoss  : {}\n'.format(self.bTipLoss  )
        s+=' - bHubLoss  : {}\n'.format(self.bHubLoss  )
        s+=' - bAIDrag   : {}\n'.format(self.bAIDrag   )
        s+=' - bTIDrag   : {}\n'.format(self.bTIDrag  )
        s+=' - WakeMod   : {}\n'.format(self.WakeMod  )
        return s


    def calcOutput(self, Omega, pitch, V0, 
            xdot=0, u_turb=0,
            a_init=None, ap_init=None):
        out = calcSteadyBEM(Omega, pitch, V0, xdot, u_turb,
            nB=self.nB, cone=self.cone0, r=self.r, chord=self.chord, twist=self.twist, polars=self.polars, # Rotor
            rho=self.rho, KinVisc=self.kinVisc,    # Environment
            nItMax=self.nbIt, aTol=self.aTol, bTipLoss=self.bTipLoss, bHubLoss=self.bHubLoss, 
            bAIDrag=self.bAIDrag, bTIDrag=self.bTIDrag, bSwirl=self.bSwirl, relaxation=self.relaxation, 
            a_init=a_init, ap_init=ap_init)
        return out

    def parametric(self, rotSpeed, pitch, windSpeed, outputFilename=None, radialBasename=None, plot=False, **kwargs):
        """ Run the BEM on multiple operating conditions
            Inputs:
            -------
            rotSpeed   [rpm]: rotational speed, array-like, size n
            pitch      [deg]: pitch angle, array-like, size n
            windSpeed  [m/s]: wind speed, array-like, size n
        """
        # --- Sanity
        createParentDir(outputFilename)
        createParentDir(radialBasename)
        nSim = len(rotSpeed)

        # --- Running BEM simulations for each operating conditions
        a0 , ap0 = None,None # Initial guess for inductions, to speed up BEM
        dfOut=None
        for i,(u0,rpm,theta), in enumerate(zip(windSpeed,rotSpeed,pitch)):
            xdot   = 0        # structrual velocity [m/s] 
            u_turb = 0        # turbulence fluctuation [m/s] 
            BEM = self.calcOutput(rpm, theta, u0, xdot, u_turb, a_init=a0, ap_init=ap0)
            # Store previous values to speed up convergence
            a0, ap0 = BEM.a, BEM.aprime
            # Export radial data to file
            if radialBasename:
                #filenameRadial = radialBasename+'_i{:04d}_ws{:02.0f}_radial.csv'.format(i,u0)
                filenameRadial = radialBasename+'ws{:02.0f}_radial.csv'.format(u0)
                BEM.WriteRadialFile(filenameRadial)
            # Cumulative storage of integrated values
            dfOut = BEM.StoreIntegratedValues(dfOut)

            self.lastRun=BEM

        if outputFilename:
            dfOut.to_csv(outputFilename, index=False, sep='\t')

        if plot:
            label = kwargs['label'] if 'label' in kwargs.keys() else None
            import matplotlib.pyplot as plt
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            ax.plot(dfOut['WS_[m/s]'].values, dfOut['AeroPower_[kW]'].values, label=label)
            ax.set_xlabel('Wind speed [m/s]')
            ax.set_ylabel('[-]')
            ax.legend()
            ax.set_title('BEM Steady - Performance curve')
            ax.tick_params(direction='in')
            return dfOut, fig

        return dfOut, None

    def findPitchForPowerCurve(self, windSpeed,  omega, powerRated=None, powerTarget=None, pitch0=0):
        """ 
        INPUTS:
         - windSpeed: wind speed (WS), array-like, size n    [m/s]
         - omega    :  rotor speed as function of WS, array like of size n [rpm]
         - powerRated: rated power, scalar   [W]
                  if provided, pitch is optimized only above rated power
             OR
           powerTarget: target power as function of WS, array like of size n [W]
                  if provided, pitch is optimized over all range of power

         - pitch0: initial pitch/pitch to apply in Region 2
        """
        import scipy.optimize as sco
        # See fWTFindPitch.m
        # --- Sanity
        if (powerTarget is None and powerRated is None) or (powerTarget is not None and powerRated is not None):
            raise Exception('Provide either powerTarget or powerRated')
        
        optimPitchForAllWS= powerRated is None

        pitch = np.zeros_like(windSpeed)
        power = np.zeros_like(windSpeed)
        # --- SubFunction to optimize pitch
        def _findPitch(Ptarget, lastPitch):
            def dP(pitch):
                return self.calcOutput(rpm, pitch, ws, a_init=a0, ap_init=ap0).Power-Ptarget
            fun = dP
            lastPitch = sco.fsolve(fun,lastPitch)
            # Call again for verif
            out = self.calcOutput(rpm, lastPitch, ws)
            dp = out.Power-Ptarget
            return lastPitch, out, dp

        # --- Find pitch 
        a0, ap0 = None, None
        lastPitch = pitch0
        if optimPitchForAllWS:
            for i,(ws, rpm, pTarget) in enumerate(zip(windSpeed, omega, powerTarget)):
                lastPitch+=1
                lastPitch, out, dp = _findPitch(pTarget, lastPitch)
                a0, ap0 = out.a, out.aprime
                pitch[i] = lastPitch
                power[i] = out.Power
                print('U={:5.1f}, P={:7.1f}kW, pitch={:5.2f}deg, RPM={:4.1f} (DP={:7.1f}kW)'.format(ws, power[i]/1000, pitch[i], rpm, dp/1000))
        else:
            for i,(ws,rpm) in enumerate(zip(windSpeed, omega)):
                # First run to check whether optimization is needed
                out = self.calcOutput(rpm, lastPitch, ws)
                if out.Power<powerRated:
                    lastPitch = pitch0
                    dp = 0
                else:
                    lastPitch = lastPitch+1
                    lastPitch, out, dp = _findPitch(powerRated, lastPitch)
                pitch[i] = lastPitch
                power[i] = out.Power
                print('U={:5.1f}, P={:7.1f}kW, pitch={:5.2f}deg, RPM={:4.1f} (DP={:7.1f}kW)'.format(ws, power[i]/1000, pitch[i], rpm, dp/1000))

        return pitch, power

    # --------------------------------------------------------------------------------}
    # --- IO 
    # --------------------------------------------------------------------------------{
    def radialDataFrame(self):
        return self.lastRun.radialDataFrame() # not very pretty


# --------------------------------------------------------------------------------}
# --- Utils 
# --------------------------------------------------------------------------------{
def createParentDir(basename):
    if basename is None:
        return
    parentDir = os.path.dirname(basename)
    if not os.path.exists(parentDir):
        os.makedirs(parentDir)

def _fAeroCoeffWrap(fPolars, alpha, phi, bAIDrag=True, bTIDrag=True):
    """Tabulated airfoil data interpolation
        Inputs
        ----------
        Polars: interpolant function for each alpha
        alpha: Angle Of Attack [rad]
        phi  : flow angle  [rad]

        Outputs
        ----------
        Cl,Cd  : lift and drag coefficients
        cnForAI: normal and tangential coefficient  (for induction computation)
    """
    alpha[alpha<-pi] += 2*pi
    alpha[alpha> pi] -= 2*pi
    Cl = np.zeros(alpha.shape)
    Cd = np.zeros(alpha.shape)
    for i,(fPolar,alph) in enumerate(zip(fPolars,alpha)):
        ClCdCm = fPolar(alph)
        Cl[i], Cd[i] = ClCdCm[0], ClCdCm[1]
    # --- Normal and tangential
    cn = Cl * cos(phi) + Cd * sin(phi)
    ct = Cl * sin(phi) - Cd * cos(phi)
    if (bAIDrag):
        cnForAI = cn
    else:
        cnForAI = Cl*cos(phi) # cnNoDrag
    if (bTIDrag):
        ctForTI = ct
    else:
        ctForTI =  Cl * sin(phi) # ctNoDrag
    return Cl, Cd, cnForAI, ctForTI

def _fInductionCoefficients(a_last, Vrel_norm, V0, F, cnForAI, ctForTI,
        lambda_r, sigma, phi, relaxation=0.4,bSwirl=True):
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
    a = 1. / ((4.*F*sin(phi)**2)/(sigma*(cnForAI+10**-8))+1) # NOTE simgularity avoided
    # CT=(1-a_last).^2.*sigma.*CnForAI./((sind(phi)).^2)
    Ct = Vrel_norm**2 * sigma * cnForAI/(V0**2)  # that's a CT loc
    # --- Hight thrust correction
    # Glauert correction
    #>>> NOTE this is:  a = a_Ct_a(Ct, a, method='Glauert') from HighThrust
    ac = 0.3
    bHigh = a > ac
    fg = 0.25*(5.-3.*a[bHigh])
    a[bHigh] = Ct[bHigh]/(4.*F[bHigh]*(1.-fg*a[bHigh]))
    #a_high=0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2)^2+4*(K*ac^2-1)));
    # --- Relaxation
    a = a*relaxation + (1.-relaxation)*a_last

    # --- Swirl
    if bSwirl is True:
        #aprime = 1/((4*F*sin(phi)*cos(phi))/(sigma*ctForTI+10**-8)-1)
        # HAWC2 method:
        aprime = (Vrel_norm**2*ctForTI*sigma)/(4.*(1.-a)*V0**2*lambda_r)
    else:
        aprime = a * 0.
    # Bounding values for safety
    aprime = np.clip(aprime,-1,1.5) 
    a      = np.clip(a     ,-1,1.5)
    Ct     = np.clip(Ct    ,-1,3)
    return a, aprime, Ct


# --------------------------------------------------------------------------------}
# --- Class to store BEM outputs 
# --------------------------------------------------------------------------------{
class SteadyBEM_Outputs:

    def radialDataFrame(BEM):
        header='r_[m] a_[-] a_prime_[-] Ct_[-] Cq_[-] Cp_[-] cn_[-] ct_[-] phi_[deg] alpha_[deg] Cl_[-] Cd_[-] Pn_[N/m] Pt_[N/m] Vrel_[m/s] Un_[m/s] Ut_[m/s] F_[-] Re_[-] Gamma_[m^2/s] uia_[m/s] uit_[m/s] u_turb_[m/s]'
        header=header.split()
        #header=' '.join(['{:14s}'.format(s) for s in header.split()])
        M=np.column_stack((BEM.r,BEM.a,BEM.aprime,BEM.Ct,BEM.Cq,BEM.Cp,BEM.cn,BEM.ct,BEM.phi,BEM.alpha,BEM.Cl,BEM.Cd,BEM.Pn,BEM.Pt,BEM.Vrel,BEM.Un,BEM.Ut,BEM.F,BEM.Re,BEM.Gamma,BEM.uia,BEM.uit,BEM.u_turb))
        return pd.DataFrame(data=M, columns=header)

    def WriteRadialFile(BEM, filename):
        df = BEM.radialDataFrame()
        np.savetxt(filename, df.values, header=' '.join(['{:14s}'.format(s) for s in df.columns]), fmt='%14.7e',comments='#')
        #df.to_csv(filename, index=False, sep='\t')

    def StoreIntegratedValues(BEM, df=None):
        S=pd.Series(dtype='float64')
        S['WS_[m/s]']         = BEM.V0
        S['RotSpeed_[rpm]']   = BEM.Omega *60/(2*np.pi)
        S['Pitch_[deg]']      = BEM.Pitch
        S['AeroThurst_[kN]']  = BEM.Thrust/1000
        S['AeroTorque_[kNm]'] = BEM.Torque/1000
        S['AeroPower_[kW]']   = BEM.Power/1000
        S['AeroFlap_[kNm]']   = BEM.Flap/1000
        S['AeroEdge_[kNm]']   = BEM.Edge/1000
        S['AeroCT_[-]']       = BEM.CT
        S['AeroCP_[-]']       = BEM.CP
        S['AeroCQ_[-]']       = BEM.CQ

        if df is None:
            df = pd.DataFrame([S])
        else:
            df = pd.concat((df, pd.DataFrame([S])))
        return df

# --------------------------------------------------------------------------------}
# --- Low level function  
# --------------------------------------------------------------------------------{
def calcSteadyBEM(Omega,pitch,V0,xdot,u_turb,
        nB, cone, r, chord, twist, polars, # Rotor
        rho=1.225,KinVisc=15.68*10**-6,    # Environment
        nItMax=100, aTol=10**-6, bTipLoss=True, bHubLoss=False, bAIDrag=True, bTIDrag=True, bSwirl=True, relaxation=0.4, a_init=None, ap_init=None):
    """ Run the BEM main loop
        Inputs:
        -------
        Omega [rpm]:
        pitch [deg]:
        twist [deg]:
        cone  [deg]:
        r     [m]  : from rhub to R
        chord [m]  :
        polars     : nSpan matrices  

        Outputs
        ----------
        BEM : class with attributes, such as BEM.r, BEM.a, BEM.Power
    """
    VHubHeight = V0
    # --- Converting units
    fulltwist = (twist+pitch) *pi/180    # [rad]
    Omega    = Omega*2*pi/60 # [rad/s]
    # --- Derived params
    rhub, R  = r[0], r[-1]
    # Computing a dr, such that sum(dr)=R-rhub
    dr    = np.diff(r) 
    MidPointAfter = np.concatenate((  r[0:-1]+dr/2 , [R] ))
    MidPointBefore= np.concatenate(( [r[0]] ,  r[1:]-dr/2))
    dr    = MidPointAfter-MidPointBefore
    cCone    = cos(cone*pi/180.)
    sigma    = chord * nB / (2.0 * pi * r * cCone)
    lambda_r = Omega * r * cCone/ V0
    # Creating interpolation functions for each polar, now in rad!
    fPolars = [interp1d(p[:,0]*pi/180,p[:,1:],axis=0) for p in polars]
    # Initializing outputs
    if a_init is None:
        a_init = np.ones((len(r)))*0.2
    if ap_init is None:
        ap_init = np.ones(len(r))*0.01
    # --- Vectorized BEM algorithm
    # Radial inputs: a_init,ap_init,r,chord,fulltwist,sigma,lambda_r,fPolars
    a, aprime  = a_init, ap_init
    for i in np.arange(nItMax):
        # --------------------------------------------------------------------------------
        # --- Step 0: Relative wind
        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------
        # --- Step 1: Wind Components
        # --------------------------------------------------------------------------------
        Ut = Omega * r * (1. + aprime)
        Un = V0 * (1. - a) - xdot + u_turb
        Vrel_norm = np.sqrt(Un** 2 + Ut** 2)
        # --------------------------------------------------------------------------------
        # --- Step 2: Flow Angle
        # --------------------------------------------------------------------------------
        phi = arctan2(Un, Ut) # flow angle [rad]
        # --------------------------------------------------------------------------------
        # --- Tip loss
        # --------------------------------------------------------------------------------
        Ftip = np.ones((len(r)))
        Fhub = np.ones((len(r)))
        IOK=sin(phi)>0.01
        try:
            if bTipLoss:
                # Glauert tip correction
                Ftip[IOK] = 2/pi*arccos(exp(-nB/2*(R-r[IOK])/(r[IOK]*sin(phi[IOK]))))
            if bHubLoss:
                # Prandtl hub loss correction
                Fhub[IOK] = 2/pi*arccos(exp(-nB/2*(r[IOK]-rhub)/(rhub*sin(phi[IOK]))));
        except:
            raise
        F=Ftip*Fhub;
        F[F<=0]=0.5 # To avoid singularities
        # --------------------------------------------------------------------------------
        # --- Step 3: Angle of attack
        # --------------------------------------------------------------------------------
        alpha = phi - fulltwist # [rad], contains pitch
        # --------------------------------------------------------------------------------
        # --- Step 4: Profile Data
        # --------------------------------------------------------------------------------
        Cl, Cd, cnForAI, ctForTI = _fAeroCoeffWrap(fPolars, alpha, phi, bAIDrag, bTIDrag)
        # --------------------------------------------------------------------------------
        # --- Step 5: Induction Coefficients
        # --------------------------------------------------------------------------------
        # Storing last values
        a_last      = a
        aprime_last = aprime
        a, aprime, CT_loc = _fInductionCoefficients(a_last,Vrel_norm,V0, F, cnForAI, ctForTI,
                                               lambda_r, sigma, phi, relaxation, bSwirl)

        if (i > 3 and (np.mean(np.abs(a-a_last)) + np.mean(np.abs(aprime - aprime_last))) < aTol):  # used to be on alpha
            break
    nIt = i + 1
    if i == nItMax-1:
        print('Maximum iterations reached : Omega=%.2f V0=%.2f' % (Omega,V0))
    #print('Converged: V0=%2.f om=%5.2f pit=%3.1f nIt=%d' % (V0,Omega,pitch, nIt))
    # --------------------------------------------------------------------------------
    # --- Step 6: Outputs
    # --------------------------------------------------------------------------------
    BEM=SteadyBEM_Outputs();
    BEM.a,BEM.aprime,BEM.phi,BEM.Cl,BEM.Cd,BEM.Un,BEM.Ut,BEM.Vrel,BEM.F,BEM.nIt = a,aprime,phi,Cl,Cd,Un,Ut,Vrel_norm,F,nIt
    # L = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Cl
    # D = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Cd
    # Radial quantities (recomputed since thought as derived outputs)
    BEM.cn = BEM.Cl * cos(BEM.phi) + BEM.Cd * sin(BEM.phi)
    BEM.ct = BEM.Cl * sin(BEM.phi) - BEM.Cd * cos(BEM.phi)
    BEM.Pn    = 0.5 * rho * BEM.Vrel**2 * chord * BEM.cn   # [N/m]
    BEM.Pt    = 0.5 * rho * BEM.Vrel**2 * chord * BEM.ct   # [N/m] 
    BEM.alpha = (BEM.phi - fulltwist)*180/pi               # [deg]
    BEM.phi   = BEM.phi*180/pi                             # [deg]
    BEM.Re    = BEM.Vrel * chord / KinVisc / 10**6  # Reynolds number in Millions
    BEM.Gamma = 0.5 * BEM.Vrel * chord * BEM.Cl   # Circulation [m^2/s]
    # Radial quantities, "dr" formulation
    BEM.ThrLoc   = dr * BEM.Pn * cCone
    BEM.ThrLocLn = BEM.Pn * cCone
    BEM.Ct       = nB * BEM.ThrLoc / (0.5 * rho * VHubHeight** 2 * (2*pi * r * cCone * dr))
    BEM.TqLoc    = dr * r * BEM.Pt * cCone
    BEM.TqLocLn  = r * BEM.Pt * cCone
    BEM.Cq      = nB * BEM.TqLoc / (0.5 * rho * VHubHeight** 2 * (2*pi * r * cCone)) * dr * r * cCone
    BEM.Cp      = BEM.Cq*lambda_r
    # --- Integral quantities
    BEM.Torque = nB * np.trapz(r * (BEM.Pt * cCone), r)  # Rotor shaft torque [N]
    BEM.Thrust = nB * np.trapz(     BEM.Pn * cCone, r)   # Rotor shaft thrust [N]
    BEM.Flap   = np.trapz( (BEM.Pn * cCone) * (r - rhub), r)      # Flap moment at blade root [Nm]
    BEM.Edge   = np.trapz(  BEM.Pt * (r * cCone) * (r - rhub), r) # Edge moment at blade root [Nm]
    BEM.Power = Omega * BEM.Torque
    BEM.CP = BEM.Power  / (0.5 * rho * V0**3 * pi * R**2)
    BEM.CT = BEM.Thrust / (0.5 * rho * V0**2 * pi * R**2)
    BEM.CQ = BEM.Torque / (0.5 * rho * V0**2 * pi * R**3)
    BEM.r=r
    BEM.R=R
    BEM.uia    = V0 * BEM.a
    BEM.uit    = Omega * r * BEM.aprime
    BEM.u_turb = np.ones(r.shape)*u_turb
    BEM.Omega = Omega
    BEM.Pitch = pitch
    BEM.V0 = V0
    return BEM


# --------------------------------------------------------------------------------}
# --- Tools to extract necessary inputs for a BEM simulation
# --------------------------------------------------------------------------------{
def FASTFile2SteadyBEM(FASTFileName):
    from welib.weio.fast_input_deck import FASTInputDeck
    F = FASTInputDeck(FASTFileName,readlist=['AD','ED','ADbld','AF'])
    if F.AD is None:
        raise Exception('Cannot open AD file referenced in:'.format(FASTFileName))
    if F.ED is None:
        raise Exception('Cannot open ED file referenced in:'.format(FASTFileName))

    try:
        rho     = float(F.fst['AirDens'])  # New OF > 3.0
    except:
        rho     = float(F.AD['AirDens'])   # Old OF <=3.0
    try:
        KinVisc = float(F.fst['KinVisc'])  # New OF > 3.0
    except:
        KinVisc = float(F.AD['KinVisc'])   # Old OF <= 3.0

    nB   =  F.ED['NumBl']
    cone = -F.ED['PreCone(1)']
    r     = F.AD.Bld1['BldAeroNodes'][:,0] + F.ED['HubRad']
    chord = F.AD.Bld1['BldAeroNodes'][:,-2] 
    twist = F.AD.Bld1['BldAeroNodes'][:,-3]
    polars=[]
    ProfileID=F.AD.Bld1['BldAeroNodes'][:,-1].astype(int)
    for ipolar in  ProfileID:
        polars.append(F.AD.AF[ipolar-1]['AFCoeff'])
    return nB,cone,r,chord,twist,polars,rho,KinVisc

def FLEX2SteadyBEM(BladeFile,ProfileFile, cone, tilt, r_hub, rho=1.225, nB=3):
    import welib.weio as weio
    pro=weio.read(ProfileFile)
    bld=weio.read(BladeFile).toDataFrame()
    bld.columns=[v.split('_[')[0] for v in bld.columns.values]
    # TODO
    #rho=1.225;
    #r_hub=2.253;
    #nB   = 3
    #cone = 5
    #tilt = 6; # TODO
    r     = bld['r'].values+r_hub
    chord = bld['Chord'].values
    twist = bld['AeroTwist'].values
    thick = bld['RelThickness'].values
    polars=[]
    ProfileSet = bld['ProfileSet'].astype(int)
    for iset,thick in zip(ProfileSet,thick):
        PSet=pro.ProfileSets[iset-1]
        AOA=PSet.polars[0][:,0]
        CL = np.column_stack([v[:,1] for v in PSet.polars])
        CD = np.column_stack([v[:,2] for v in PSet.polars])
        CM = np.column_stack([v[:,3] for v in PSet.polars])
        fCL = interp1d(PSet.thickness,CL,axis=1)
        fCD = interp1d(PSet.thickness,CD,axis=1)
        fCM = interp1d(PSet.thickness,CM,axis=1)
        myCL=fCL(thick)
        myCD=fCD(thick)
        myCM=fCM(thick)
        polar=np.zeros(PSet.polars[0].shape)
        polar[:,0]=AOA
        polar[:,1]=myCL
        polar[:,2]=myCD
        polar[:,3]=myCM
        polars.append(polar)
    return nB,cone,r,chord,twist,polars,rho




if __name__=="__main__":
    """ See examples/ for more examples """

    # --- Read a FAST model to get the main parameters needed
    nB,cone,r,chord,twist,polars,rho,KinVisc = FASTFile2SteadyBEM('../../data/NREL5MW/Main_Onshore.fst')

    # --- Run BEM on a set of operating points
    WS =[5,10]
    RPM=[7,12]
    a0, ap0  = None,None  # inductions, used to speed up parametric study
    for i,ws in enumerate(WS):
        V0        = WS[i]
        Omega     = RPM[i]
        pitch=2     #[deg]
        xdot=0      #[m/s]
        u_turb=0    #[m/s]
        BEM=calcSteadyBEM(Omega,pitch,V0,xdot,u_turb,
                    nB,cone,r,chord,twist,polars,
                    rho=rho,KinVisc=KinVisc,bTIDrag=False,bAIDrag=True,
                    a_init =a0,
                    ap_init=ap0
                    )
        a0, ap0 = BEM.a, BEM.aprime
        print('WS ',V0, 'Power',BEM.Power,'Thrust',BEM.Thrust)

        # --- Save radial distribution to a csv file
        filename='_BEM_ws{}_radial.csv'.format(V0)
        BEM.WriteRadialFile(filename)

