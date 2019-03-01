# -*- coding: utf-8 -*-
""" Multiple classes definition from the Blade Element Momentum
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
This class structure follows the Matlab implementation of Dr Emmnanuel Branlard
<emmanuel.branlard@gmail.com>. The present steady BEM implementation is a very limited part
of the Matlab BEM, therefore some classes are limited in size as only steady state related variables
are used. Future development will include unsteady BEM implementation
with non symmetric inflow conditions (shear and veer) as well as blade structural properties.
"""

from math import gamma, pi
import numpy as np


class InitRotor:
    """ This class holds rotor related parameters
    """
    def __init__(self):
        # Required Params
        self.nB = -1.0  # number of blades
        self.R = -1.0  # Rotor radius [m]
        self.BladeLength = -1.0  # Rotor radius [m]
        self.rhub = 0.0  # Hub  length [m]
        self.cone = 0.0  # cone angle [deg]
        self.ne = -1.0 # number of blade elements in BEL
        self.HubHeight = 60.0 # hub height im m
        self.Omega = -1.  # Nominal Rotational speed [rad/s]

    def Set(self, **parRotor):
        """Parsing data ."""
        for key, value in parRotor.iteritems():
            setattr(self, key, value)


class InitEnv:
    """ This class holds ambient conditions parameters
    """
    def __init__(self):
        self.rho = 1.225  # air density [kg/m^3]
        self.g = 9.81  # gravity [m/s^2]
        self.KinVisc = 15.68 * 10 ** -6  # Knematic viscosity nu [m^2/s]
    def Set(self, **parEnv):
        """Parsing data ."""
        for key, value in parEnv.iteritems():
            setattr(self, key, value)

class InitSpec:
    """ This class holds turbine operating specifications
    """
    def __init__(self, Env, Rotor):
        # Required Params
        self.Cp_rated = -1.0  # estimation of Cp_max at mean wind speed from 31784
        self.P_rated = -1.0  # 500kW
        # self.V_rated_thumbs = gamma(1. + 1. / Env.k) * Env.A + 6.  # rule of thumb mean WS + 6m/s
        self.V_rated = (self.P_rated / (self.Cp_rated * 0.5 * Env.rho * Rotor.R ** 2. * pi)) ** (1. / 3.)
        self.Omega_rated = 0.
        self.TSR_rated = 0. # tip speed ratio at rated
        self.RPM = [] # rotational speed in RPM
        self.Pitch = [] # pitch angle in deg
        self.WS = [] # mean hub height wind speed

    def Set(self, **parSpec):
        """Parsing data ."""
        for key, value in parSpec.iteritems():
            setattr(self, key, value)

class InitMisc:
    """ This class holds turbine operating initial states
    """
    def __init__(self):
        # Required Params
        self.Profiles = -1. # profile sets  in AE HAWC2 file
        self.Format = 'hawc' # format of turbine files
        self.WTname = 'NREL5MW' # turbine name for input files

    def Set(self, **parMisc):
        """Parsing data ."""
        for key, value in parMisc.iteritems():
            setattr(self, key, value)


class InitAlgo:
    """ This class holds BEM algorithm flags and parameters
    """
    def __init__(self):
        self.nbIt = 200  # maximum number of iterations in BEM
        self.aTol = 10 ** -6 # tolerance for axial induction factor convergence
        self.relaxation = 0.5  # relaxation factor in axial induction factor
        self.CTcorrection = 'GlauertCT'  #  type of CT correction more model implementated in the future like 'spera'
        self.Ngrid = 1.0
        self.bSwirl = True  # swirl flow model enabled / disabled
        self.SwirlMethod = 'HAWC' # type of swirl model
        self.bTipLoss = True # enable / disable tip loss model
        self.bHubLoss = False # enable / disable hub loss model
        self.bTipLossCl = False # enable / disable Cl loss model
        self.TipLossMethod = 'Glauert'  # type of tip loss model
        self.bDynaStall = 0 # dynamic stall model
        self.bAIDrag = True # influence on drag coefficient on normal force coefficient
        self.bTIDrag = True # influence on drag coefficient on tangential force coefficient
        self.bReInterp = False # interpolate the input tabulated airfoil data for Reynolds variation
        self.bThicknessInterp = True # interpolate the input tabulated airfoil data for thickness variation
        self.bRoughProfiles = False # use rough profiles for input airfoil data

    def Set(self, **parAlgo):
        """Parsing data ."""
        for key, value in parAlgo.iteritems():
            setattr(self, key, value)


class InitSim:
    """ This class holds further simulation parameters
    """
    def __init__(self):
        self.WT = '' # turbine name
        self.Name = '' # turbine name, redundant to Name, adaptation from Matlab library
        self.rho = 1.225 # air density
        self.KinVisc = 15.68 * 10 ** -6 # kinematic viscosity
        self.WS = 0. # wind speed m/s
        self.RPM = 0. # rotational speed RPM
        self.PITCH = 0. # pitch angle [deg]
        self.YAW = 0. # yaw angle [deg]
        self.Omega = 0. # rotational speed rad/s

    def Set(self, **parSim):
        """Parsing data ."""
        for key, value in parSim.iteritems():
            setattr(self, key, value)


class InitWind:
    """ This class holds inflow wind model. Further implementation will include non symmetric inflow
    (shear and veer) as well as inflow turbulence
    """
    def __init__(self):
        self.V0 = 0. # free stream velocity at hub height

    def Set(self, **parWind):
        """Parsing data ."""
        for key, value in parWind.iteritems():
            setattr(self, key, value)

class InitBEM:
    """ This class holds outputs from the BEM simulations
    """
    def __init__(self, ne):
        """ Initializing the data ."""
        self.phi = np.zeros((ne)) # flow angle in deg
        self.alpha = np.zeros((ne)) # angle of attack in deg
        self.a = np.zeros((ne)) # axial induction
        self.a_last = np.zeros((ne)) # last iteration axial induction for iterative loop
        self.aprime = np.zeros((ne)) # tangential induction
        self.aprime_last = np.zeros((ne)) # last iteration tangential induction
        self.Cd = np.zeros((ne)) # drag coefficient
        self.Cl = np.zeros((ne)) # lift coefficient
        self.Cn = np.zeros((ne)) # normal force coefficient
        self.Ct = np.zeros((ne)) # tangential force coefficient
        self.CT = np.zeros((ne)) # thrust coefficient
        self.Vrel = np.zeros((ne)) # relative velocity in m/s
        self.Re = np.zeros((ne)) # Reynolds number based on chord
        self.F = np.zeros((ne)) # loss corrections
        self.Fperf = np.zeros((ne))
        self.Fshen = [] # Shen's model loss corrections
        self.Un = np.zeros((ne)) # normal velocity
        self.Ut = np.zeros((ne)) # tangential velocity
        self.nIt = np.zeros((ne)) # number of iterations
        self.omega = [] # rotational speed in rad/s
        self.r= np.zeros((ne))
    def Set(self, **parBEM):
        """Parsing data ."""
        for key, value in parBEM.iteritems():
            setattr(self, key, value)

