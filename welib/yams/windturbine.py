"""
Generic wind turbine class for structural models
based on Generic bodies classes

Intended to be used by other models such as TNSB, FNTSB, for yams_rec, yams_sympy
Unification not done yet.

Example:

    WT = WindTurbineStructure.fromFAST('Main.fst')
    print(WT.nac)
    print(WT.twr)
    print(WT.RNA)


"""

import os
import numpy as np
import pandas as pd
import copy
import sys
import welib.weio as weio
from collections import OrderedDict
from welib.essentials import *
from welib.yams.bodies import RigidBody, FlexibleBody, FASTBeamBody, BeamBody
from welib.yams.yams_rec import YAMSRecRigidBody
from welib.yams.yams_rec import YAMSRecFASTBeamBody
from welib.yams.rotations import R_x, R_y, R_z, rotMat
from welib.yams.kinematics import rigidBodyMotion2Points
from welib.yams.utils import translateInertiaMatrixToCOG, translateInertiaMatrix
from welib.weio.dataframe import WEIODataFrame

from welib.tools.pandalib import remap_df, remove_duplicated_col_df
from welib.tools.decorators import require_attrs

from welib.tools.compare import compare

from welib.yams.section_loads import beamSectionLoads3D, beamSectionLoadsFromShapeFunctions
from welib.hydro.morison import monopileHydroLoads1D







class WindTurbineStructure():
    def __init__(self):
        self.bld = None # list of blades
        self.hub = None
        self.gen = None
        self.rot = None         # origin at R, rigid body bld+hub
        self.rotgen = None      # origin at R, rigid body bld+hub+genLSS
        self.RNA    = None      # origin at N, rigid body bld+hub+gen+nac+yawBr
        self.RNA_noYawBr = None # origin at N, rigid body bld+hub+gen+nac
        self.hubgen = None         # Hub + Gen
        self.nac = None
        self.yawBr = None        # origin at N
        self.yaw = None   # TODO yaw or yawbr
        self.twr = None
        self.fnd = None  # mnp or floater
        self.grd = None 
        # Geometry
        self.r_ET_inE=None
        self.r_TN_inT=None
        self.r_NS_inN = None
        self.r_NR_inN = None
        self.r_SR_inS = None

        self.shaft_tilt = None # [rad] # TODO introduce nac_titl
        self.blade_cone = None # [rad]
        self.nac_yaw    = None # TODO decide what this means (offset, DOF)

        # --- From TNSB
        self.MM   = None
        self.KK   = None
        self.DD   = None
        self._q   = None
        self.bTiltBeforeNac = None

        # --- Environment
        self._gravity = None
        self.pHD = None # Hydrodynamics
        self.pAD = None # Aerodynamics
        self.pSS = None

        # --- Simulation data
        self.DOF    = None # List of {'name':'', 'active':None, 'q0':None, 'qd0':None,'q_channel':None, 'qd_channel':None, 'qdd_channel':None}

        # --- Any other data
        self.additional_properties=[] # for user output in __repr__, so I remember what we have in the object

        # Derived properties
        #self.DOFname
        #self.q0     
        #self.qd0    

    def copy(self):
        return copy.deepcopy(self)

    def __repr__(B):
        s='<Generic {} object>:\n'.format(type(B).__name__)
        try:
            s+=' - DOF: list of dict with keys {} \n'.format(B.DOF[0].keys())
        except:
            pass
        s+=' * DOFname    : {} \n'.format(B.DOFname)
        s+=' * q0         : {}\n'.format(B.q0)
        s+=' * qd0        : {}\n'.format(B.qd0)
        s+=' * q_channels : {}\n'.format(B.q_channels)
        s+=' * qd_channels: {}\n'.format(B.qd_channels)
        s+=' * qdd_channels: {}\n'.format(B.qdd_channels)
        s+=' Bodies: bld, hub, gen, rot, nac, twr, fnd \n'

        s+=' useful getters: activeDOFs, channels:\n'
        s+=' useful functions:  \n'
        s+='    simulate_py\n'
        s+='    simulate_py_lin\n'
        return s

    @staticmethod
    def fromFAST(fstFilename):
        return FASTWindTurbine(fstFilename).WT

    def reshape_3array_to_atleast_2d(self):
        for name, value in vars(self).items():
            if isinstance(value, np.ndarray) and value.shape == (3,):
                setattr(self, name, value.reshape((3,1)))

    @property
    def q(self):
        if self._q is None:
            return None
        return np.asarray(self._q).reshape((-1, 1))

    @q.setter
    def q(self, q):
        self._q = np.asarray(q).reshape((-1, 1))

    @property
    def q0(self):
        return pd.Series(np.array([dof['q0'] for dof in self.DOF if dof['active']]), index=self.DOFname)
    @property
    def qd0(self):
        return pd.Series(np.array([dof['qd0'] for dof in self.DOF if dof['active']]), index=self.dDOFname)
    @property
    def z0(self):
        z0 = np.concatenate((np.asarray(self.q0), np.asarray(self.qd0)))
        return pd.Series(z0, index=self.zname)

    @property
    def DOFname(self):
        return np.array([dof['name'] for dof in self.DOF if dof['active']])
    @property
    def dDOFname(self):
        return np.array(['d'+dof['name'] for dof in self.DOF if dof['active']])
    @property
    def zname(self):
        qnames = self.DOFname
        qdnames = ['d' + s for s in qnames]
        return list(qnames) + list(qdnames)

    @property
    def q_channels(self):
        return np.array([dof['q_channel'] for dof in self.DOF if dof['active']])
    @property
    def qd_channels(self):
        return np.array([dof['qd_channel'] for dof in self.DOF if dof['active']])
    @property
    def qdd_channels(self):
        return np.array([dof['qdd_channel'] for dof in self.DOF if dof['active']])

    @property
    def FASTDOFScales(self):
        """ scales to go from raw DOF to the same units as OpenFAST """
        DOFscales=[]
        for s in self.channels:
            if s.find('[deg]')>0:
                DOFscales.append(180/np.pi)
            elif s.find('TSS1')>0:
                DOFscales.append(-1)
            elif s.find('[rpm]')>0:
                DOFscales.append(60/(2*np.pi))
            else:
                DOFscales.append(1)
        return DOFscales


    @property
    def activeDOFs(self):
        return [dof for dof in self.DOF if dof['active']]

    @property
    def channels(self):
        """ """
        chan = self.pos_channels
        chan+= self.vel_channels
        return chan

    @property
    def pos_channels(self):
        """ """
        return [dof['q_channel']  if dof['q_channel']  is not None else dof['name']     for dof in self.DOF if dof['active']]
        
    @property
    def vel_channels(self):
        return [dof['qd_channel'] if dof['qd_channel'] is not None else 'd'+dof['name'] for dof in self.DOF if dof['active']]

    @property
    def acc_channels(self):
        return [dof['qdd_channel'] if dof['qdd_channel'] is not None else 'dd'+dof['name'] for dof in self.DOF if dof['active']]

    @property
    def mass(self):
        return self.WT_rigid.mass # useful for hydro Fz

    @property
    def gravity(self):
        if self._gravity is not None:
            return self._gravity
        else:
            raise Exception('Gravity not set')

    @gravity.setter
    def gravity(self, value):
        self._gravity = value


    def yams_parameters(self, flavor='allbodies', 
            J_at_Origin=True,
            MoorAtRef=True,
            ):
        """ 
        Export parameters needed by yams based on WT bodies
        INPUTS:
        - flavor: ['rigidbody', 'allbodies']
        """
        WT=self;
        # --- Dict needed by structural script 
        p = dict()
        p['g'] = self.gravity

        p['M_B']      = WT.WT_rigid.mass # useful for hydro Fz
        p['tilt']     =-WT.ED['ShftTilt']*np.pi/180 # in rad
        if flavor=='onebody':
            # One body for all turbine
            p['x_BG']     = WT.WT_rigid.masscenter[0]  # From body origin to body COG
            p['y_BG']     = WT.WT_rigid.masscenter[1]
            p['z_BG']     = WT.WT_rigid.masscenter[2]
            p['z_B0']     = - WT.WT_rigid.pos_global[2] # From body origin to refHeight taken as (0,0,0)
            p2=dict()
            p2['J_xx_BG']  = WT.WT_rigid.masscenter_inertia[0,0]
            p2['J_yy_BG']  = WT.WT_rigid.masscenter_inertia[1,1]
            p2['J_zz_BG']  = WT.WT_rigid.masscenter_inertia[2,2]
            p2['J_yx_BG']  = WT.WT_rigid.masscenter_inertia[0,1]
            p2['J_zx_BG']  = WT.WT_rigid.masscenter_inertia[0,2] # Likely not needed
            p2['J_zy_BG']  = WT.WT_rigid.masscenter_inertia[1,2] # Likely not needed
            p2['J_xx_BO']  = WT.WT_rigid.inertia[0,0]
            p2['J_yy_BO']  = WT.WT_rigid.inertia[1,1]
            p2['J_zz_BO']  = WT.WT_rigid.inertia[2,2]
            p2['J_yx_BO']  = WT.WT_rigid.inertia[0,1]
            p2['J_zx_BO']  = WT.WT_rigid.inertia[0,2]
            p2['J_zy_BO']  = WT.WT_rigid.inertia[1,2]
            if J_at_Origin:
                p['J_xx_B']   = p2['J_xx_BO']
                p['J_yy_B']   = p2['J_yy_BO']
                p['J_zz_B']   = p2['J_zz_BO']
                p['J_yx_B']   = p2['J_yx_BO']
                p['J_zx_B']   = p2['J_zx_BO']
                p['J_zy_B']   = p2['J_zy_BO']
            else:
                p['J_xx_B']   = p2['J_xx_BG']
                p['J_yy_B']   = p2['J_yy_BG']
                p['J_zz_B']   = p2['J_zz_BG']
                p['J_yx_B']   = p2['J_yx_BG']
                p['J_zx_B']   = p2['J_zx_BG']
                p['J_zy_B']   = p2['J_zy_BG']
            p['z_OT']     = WT.twr.pos_global[2]         # distance from "Origin" (MSL) to tower base
            p['z_EF']     = WT.fnd.pos_global[2]         # distance from "Origin" (MSL) to PtfmRefz

            # IMU
            p['x_BI'] = WT.ED['NcIMUxn']
            p['y_BI'] = WT.ED['NcIMUyn']
            p['z_BI'] = WT.ED['NcIMUzn'] + WT.ED['TowerHt'] - WT.ED['PtfmRefzt'] # TODO or PtfmRefzt?

        if flavor in ['allbodies']:
            p['z_FG']     = WT.fnd.masscenter[2]
            p['M_F']      = WT.fnd.mass
            p['J_xx_F']   = WT.fnd.masscenter_inertia[0,0]
            p['J_yy_F']   = WT.fnd.masscenter_inertia[1,1]
            p['J_zz_F']   = WT.fnd.masscenter_inertia[2,2]
            p['x_NR']     = WT.r_NR_inN[0]                    # x-coord from N to R in nac-coord
            p['z_NR']     = WT.r_NR_inN[2]                    # z-coord from N to R in nac-coord
            p['x_RNAG']   = WT.RNA.masscenter[0]            # x-coord from N to RNA_G in nac-coord
            p['z_RNAG']   = WT.RNA.masscenter[2]            # z-coord from N to RNA_G in nac-coord
            p['M_RNA']    = WT.RNA.mass                   # Total mass of RNA
            p['J_xx_RNA'] = WT.RNA.masscenter_inertia[0,0]           # Inertia of RNA at RNA_G in nac-coord
            p['J_yy_RNA'] = WT.RNA.masscenter_inertia[1,1]           # Inertia of RNA at RNA_G in nac-coord
            p['J_zz_RNA'] = WT.RNA.masscenter_inertia[2,2]           # Inertia of RNA at RNA_G in nac-coord
            p['J_zx_RNA'] = WT.RNA.masscenter_inertia[0,2]           # Inertia of RNA at RNA_G in nac-coord
            p['L_T']      = WT.twr.length
            p['z_OT']     = WT.twr.pos_global[2]         # distance from "Origin" (MSL) to tower base
            p['z_EF']     = WT.fnd.pos_global[2]         # distance from "Origin" (MSL) to PtfmRefz
            p['M_T']      = WT.twr.MM[0,0]
            p['z_TG']     = WT.twr.masscenter[2]
            p['J_xx_T']   = WT.twr.masscenter_inertia[0,0]
            p['J_yy_T']   = WT.twr.masscenter_inertia[1,1]
            p['J_zz_T']   = WT.twr.masscenter_inertia[2,2]
            p['Oe_T']     = WT.twr.Oe6
            p['Gr_T']     = WT.twr.Gr
            p['Ge_T']     = WT.twr.Ge
            p['MM_T']     = WT.twr.MM
            p['u_xT1c']   = 1 # TODO remove in equations
            # TODO TODO TODO THIS IS NOT GENERIC
            p['v_yT1c']   = WT.twr.Bhat_t_bc[1,0]  # Mode 1  3 x nShapes
            if WT.twr.nf>1:
                try:
                    p['v_xT2c']   = WT.twr.Bhat_t_bc[0,1]  # Mode 2
                except:
                    # TODO use directions
                    WARN('YAMS: WindTurbine, which mode shape is which for tower top rotations factors, TODO.')
            p['DD_T']     = WT.twr.DD
            p['KK_T']     = WT.twr.KK
            # Rotor (blades + hub)
            p['Jxx_R']    = WT.rotgen.masscenter_inertia[0,0]
            p['JO_R']     = WT.rotgen.masscenter_inertia[1,1]
            p['M_R']      = WT.rot.mass
            # Nacelle 
            p['J_xx_N']   = WT.nac.masscenter_inertia[0,0]
            p['J_yy_N']   = WT.nac.masscenter_inertia[1,1]
            p['J_zz_N']   = WT.nac.masscenter_inertia[2,2]
            p['J_zx_N']   = WT.nac.masscenter_inertia[0,2]   
            p['M_N']      = WT.nac.mass
            p['x_NG']     = WT.nac.masscenter[0]            # x-coord from N to nac G in nac-coord
            p['z_NG']     = WT.nac.masscenter[2]            # z-coord from N to nac G in nac-coord

            # Blades
            p['r_h']      = WT.ED['HubRad']
            p['theta_c']  = WT.ED['PreCone(1)']
            p['theta_p']  = WT.ED['BlPitch(1)']
            p['psi_0']    = WT.ED['Azimuth']
            p['z_BG']     = WT.bld[0].masscenter[2]
            p['J_xx_B']   = WT.bld[0].masscenter_inertia[0,0]
            p['J_yy_B']   = WT.bld[0].masscenter_inertia[1,1]
            p['J_zz_B']   = WT.bld[0].masscenter_inertia[2,2]
            p['Oe_B']     = WT.bld[0].Oe6
            p['Gr_B']     = WT.bld[0].Gr
            p['Ge_B']     = WT.bld[0].Ge
            p['MM_B']     = WT.bld[0].MM
            p['DD_B']     = WT.bld[0].DD
            p['KK_B']     = WT.bld[0].KK
            p['MM1_B']    = np.array(WT.bld[0].MM1)

            # IMU
            p['x_TI'] = WT.ED['NcIMUxn']
            p['y_TI'] = WT.ED['NcIMUyn']
            p['z_TI'] = WT.ED['NcIMUzn'] 

            # Flap 1
    #         p['MM_B'][6,6] =  9.249926E+02
    #         p['KK_B'][6,6] =  1.669818E+04
    #         p['DD_B'][6,6] =  3.752971E+01

            # Flap 2
    #         p['MM_B'][6,6] =  5.614598E+02 
    #         p['KK_B'][6,6] =  8.495747E+04
    #         p['DD_B'][6,6] =  6.595256E+01
            # Edge 1
#             p['MM_B'][6,6] =  1.430779E+03
#             p['KK_B'][6,6] =  6.682207E+04
#             p['DD_B'][6,6] =  9.337225E+01


# Blade generalized mass matrix, Blade 1:
#   9.249926E+02  0.000000E+00  0.000000E+00
#   0.000000E+00  5.614598E+02  0.000000E+00
#   0.000000E+00  0.000000E+00  1.430779E+03
# Blade generalized stiffness matrix, Blade 1:
#   1.669818E+04 -1.453662E+03  0.000000E+00
#  -1.453662E+03  8.495747E+04  0.000000E+00
#   0.000000E+00  0.000000E+00  6.682207E+04
# Blade generalized damping matrix, Blade 1:
#   3.752971E+01 -1.128479E+00  0.000000E+00
#  -3.267153E+00  6.595256E+01  0.000000E+00
#   0.000000E+00  0.000000E+00  9.337225E+01


# Blade generalized mass matrix, Blade 1:
#   9.249926E+02  0.000000E+00  0.000000E+00
#   0.000000E+00  5.614598E+02  0.000000E+00
#   0.000000E+00  0.000000E+00  1.430779E+03
# Blade generalized stiffness matrix, Blade 1:
#   1.669818E+04 -1.453662E+03  0.000000E+00
#  -1.453662E+03  8.495747E+04  0.000000E+00
#   0.000000E+00  0.000000E+00  6.682207E+04
# Blade generalized damping matrix, Blade 1:
#   3.752971E+01 -1.128479E+00  0.000000E+00
#  -3.267153E+00  6.595256E+01  0.000000E+00

#   0.000000E+00  0.000000E+00  9.337225E+01

        # Mooring restoring
        if flavor=='onebody':
            p['z_BM']      = 0
        elif flavor in ['allbodies']:
            p['z_BM']      = 0
            p['z_TM']      = 0
            p['K_x_M']     = 0
            
            p['K_z_M']     = 0
            p['K_phi_x_M'] = 0
            p['K_phi_y_M'] = 0
            p['K_phi_z_M'] = 0
        if WT.MAP is not None:
            if WT.MAP._K_lin is None:
                if MoorAtRef:
                    #MoorP = (0,0,p['z_OT'])
                    MoorP = (0,0,p['z_EF'])
                else:
                    MoorP = (0,0,0)
                K_Moor,_ = WT.MAP.stiffness_matrix(epsilon=1e-2, point=MoorP)
            else:
                K_Moor = WT.MAP._K_lin
            for i in range(6):
                for j in range(6):
                    if j>=i:
                        p['KM_{}{}'.format(i,j)] = K_Moor[i,j] # TODO the best is to add it to the stiffness matrix..
        else:
            for i in range(6):
                for j in range(6):
                    if j>=i:
                        p['KM_{}{}'.format(i,j)] = 0


        p['z_T0'] = -p['z_OT'] # TODO get rid of it
        #p['z_T0'] = -p['z_OT'] # TODO get rid of it

        ## Buoyancy
        #if flavor=='onebody':
        #    print('>>> TODO Buoyancy point')
        #    p['z_BB']      = 0
        #else:
        #    p['z_TB']      = 0

        return p

    def checkPackage(self, pkg, verbose=False):
        # --- Checking package info
        info = pkg.info()
        nDOFExpected= info['nq']
        if verbose:
            print('DOFs:', self.DOFname, 'Model:',info['name'], 'nDOF:',nDOFExpected )
        if len(self.DOFname)!=nDOFExpected:
            raise Exception('Inconsistency in number of DOFs')

    def simulate_py(self, pkg, p, time, u=None, acc=False, forcing=False, verbose=False):
        """ Perform non-linear simulation based on a `model` (python package) generated by yams_sympy """
        from welib.yams.models.packman import simulate

        # --- Checking package info
        self.checkPackage(pkg)

        if verbose:
            print('------------------ WINDTURBINE: NON LINEAR SIMULATION -----------------------')
        q0  = np.asarray(self.q0)
        qd0 = np.asarray(self.qd0)
        resNL, sysNL, dfNL = simulate(pkg, time, q0, qd0=qd0, p=p, u=u, acc=acc, forcing=forcing, 
                sStates=self.channels, Factors=self.FASTDOFScales, sAcc=self.acc_channels)
        if verbose:
            print('-----------------------------------------------------------------------------')

        return resNL, sysNL, dfNL


    def py_lin(self, pkg, p, time, uop=None, qop=None, qdop=None, du=None, MCKextra=None, MCKu=None, noBlin=False, verbose=False):
        """ Perform linear simulation based on a model (python package) generated by yams_sympy """
        # TODO TODO TODO MOVE ME TO packman
        from welib.yams.models.packman import linearModel
        if verbose:
            print('------------------ WINDTURBINE: LINEAR SIMULATION----------------------------')
        # --- Checking package info
        self.checkPackage(pkg)
        if qop is None:
            qop = self.q0*0 
        qop=np.asarray(qop)
        if qdop is None:
            qdop= self.qd0*0 
        qdop=np.asarray(qdop)
        # --- Initial conditions (of pertubations)
        dq0  = np.asarray(self.q0 ) - qop
        dqd0 = np.asarray(self.qd0) - qdop

        sysLI = linearModel(pkg, p, dq0=dq0, dqd0=dqd0, time=time, uop=uop, qop=qop, qdop=qdop, du=du, MCKextra=MCKextra, MCKu=MCKu, noBlin=noBlin,
                sX=self.q_channels, sXd=self.qd_channels)
        print('-----------------------------------------------------------------------------')
        return sysLI


    def simulate_py_lin(self, pkg, p, time, uop=None, qop=None, qdop=None, du=None, MCKextra=None, MCKu=None, acc=False, forcing=False, noBlin=False, verbose=False):
        """ Perform linear simulation based on a model (python package) generated by yams_sympy """
        if verbose:
            print('-------------------------- LINEAR SIMULATION --------------------------------')
        # --- Checking package info
        self.checkPackage(pkg)
        info = pkg.info()
        # --- Getting linear model
        #sysLI = MechSystem(M=M_lin, K=K_lin, C=C_lin, F=fF, x0=dq0, xdot0=dqd0)
        sysLI = self.py_lin(pkg, p, time, uop=uop, qop=qop, qdop=qdop, du=du, MCKextra=MCKextra, MCKu=MCKu, noBlin=noBlin)

        # --- Setup Mech system (for time integration)
        resLI, _ = sysLI.integrate(time, method='RK45') # **options):

        # --- Convert to dataframe
        calc=''
        if acc:
            calc+='xdd,'
        if forcing:
            calc+='f,'
        dfLI = sysLI.res2DataFrame(resLI, sStates = self.channels, Factors=self.FASTDOFScales, x0=qop, xd0=qdop, calc=calc, sAcc=self.acc_channels)

        if verbose:
            print('-----------------------------------------------------------------------------')
        return resLI, sysLI, dfLI

    def picklable(self):
        """ Make the object picklable..."""
        if self.MAP:
            self.MAP=None # Library is ctype, not picklable...



    def kinematics(self, qDict, qdDict, qddDict=None):
        """ Update kinematics from fnd to blades """

        fnd = self.fnd
        twr = self.twr
        nac = self.nac
        r_F0     = fnd.pos_global_init # np.array((0, 0, ED['PtfmRefzt']))
        r_T0     = twr.pos_global_init # np.array((0, 0, ED['TowerBsHt']))
        s_NGn0   = nac.masscenter # TODO
        d = kinematics(qDict, qdDict, qddDict, r_F0=r_F0, r_T0=r_T0, twr=twr, 
                s_NGn0=s_NGn0,
                tilt=self.shaft_tilt,
                algo = self.algo)

        # -- RNA (without Yaw Br) COG
        s_NGrna0_in_N = self.RNA_noYawBr.masscenter
        dRNA = rigidBodyKinematics(s_NGrna0_in_N, d['r_N'], d['R_g2n'], v_N=d['v_N'], omega_n=d['omega_n'], a_N=d['a_N'], omegad_n=d['omegad_n'], point_name='Grna', source_name='N')
        d.update(dRNA)

        # -- IMU Kinematics
        s_NIMU_in_N = np.array([self.ED['NcIMUxn'], self.ED['NcIMUyn'], self.ED['NcIMUzn']])
        dIMU = rigidBodyKinematics(s_NIMU_in_N, d['r_N'], d['R_g2n'], v_N=d['v_N'], omega_n=d['omega_n'], a_N=d['a_N'], omegad_n=d['omegad_n'], point_name='IMU', source_name='N')
        d.update(dIMU)

        # Store in bodies
        fnd.pos_global = d['r_F']
        fnd.R_b2g      = d['R_g2f'].T

        return d


    # --------------------------------------------------------------------------------}
    # --- Sea state related (might need an object in the future)
    # --------------------------------------------------------------------------------{
    @require_attrs(['pSS', 'gravity'])
    def SS_setComponents(self, compFile=None, ap=None, fp=None, epsp=None):
        from welib.hydro.wavekin import wavenumber
        pSS = self.pSS
        if compFile is not None:
            dfComp = weio.read(compFile).toDataFrame()
            pSS['ap']   = dfComp['Amplitude_[m]']
            pSS['fp']   = dfComp['Frequency_[Hz]']
            pSS['epsp'] = dfComp['Phase_[rad]']
        else:
            print('>>> Using Wave of amplitude 3 and period 12 for now')
            pSS['ap']   = np.array([3])     # Amplitudes
            pSS['fp']   = np.array([1/12])  # frequencies [Hz]
            pSS['epsp'] = np.array([np.pi]) # Deterministic phas
        # --- Wave Kinematics
        pSS['kp'] = wavenumber(pSS['fp'], pSS['WaterDepth'], self.gravity) # Wave numbers
        pSS['compFile'] = compFile
        return pSS

    @require_attrs(['pSS'])
    def SS_computeEta(self, time):
        from welib.hydro.wavekin import elevation2d
        pSS = self.pSS
        if 'ap' not in pSS:
            raise Exception('windturbine: Cannot compute Eta, components not set yet (call SS_setComponents)')

        time = np.asarray(time, dtype=float)
        eta_sim = elevation2d(pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], time, x=0)
        if len(time)>0:
            dt = time[1]-time[0]
            eta_dot = np.gradient(eta_sim, dt) # Velocity state
        else:
            print('[WARN] time is empty, eta dot is zero')
            eta_dot = eta_sim*0
        pSS['eta_time'] = time
        pSS['eta']      = eta_sim
        pSS['eta_dot']  = eta_dot
        return pSS


    # --------------------------------------------------------------------------------}
    # --- Hydro 
    # --------------------------------------------------------------------------------{
    @require_attrs(['pHD', 'fnd'])
    def HD_setShapeFunction(self, hydroShape):
        if hydroShape is None:
            raise Exception('hydroShape is None')
        pHD = self.pHD
        dfH = weio.read(hydroShape).toDataFrame()
        if not np.array_equal(dfH['z_[m]'],pHD['zDepth']):
            raise Exception('z depth different with shape function and hydrodyn, TODO interpolate')
        pHD['phi']  = dfH['phi_[Ns/m^2]'].values # Used to be called k_h_z
        pHD['phit'] = dfH['phit_[-]'].values
        bWet = pHD['zDepth']<=0
        k_h = np.zeros(len(self.WT.fnd.PhiU))
        for i, phi in enumerate(self.WT.fnd.PhiU):
            phi_x = phi[0,:]
            k_h[i] = np.trapezoid(pHD['phi'][bWet] * phi_x[bWet], pHD['zDepth'][bWet])
        print('k_h     : ', k_h)
        pHD['k_h'] = k_h



    def _insertOFDOFsInDF(self, df, fill_value = 0, verbose=False):
        """ 
        Insert missing DOF time series in dataframe.
            If a DOF is turned off, we insert zero
            Otherwise, we warn the user
        """ 
        from welib.fast.dofs import COLMAP_OFout_TO_QOF, QOF

        #return np.array([dof['q_channel'] for dof in self.DOF if dof['active']])
        _sq   = self.q_channels
        _sqd  = self.qd_channels
        _sqdd = self.qdd_channels

        keys = list(df.keys())
        df = remap_df(df, COLMAP_OFout_TO_QOF, bColKeepNewOnly=False, inPlace=False, verbose=verbose, raiseIfAbsent=False)
        df = remove_duplicated_col_df(df)

        active_OF_DOF = set(self.q_channels + self.qd_channels + self.qdd_channels)
        df_columns = set(df.columns)
        #missing_dofs = set(QOF) - set(df.columns)
        for OF_dof in QOF:
            if OF_dof not in df_columns:
                if OF_dof in active_OF_DOF:
                    FAIL(f"DOF '{OF_dof}' is missing from DataFrame but it is Active. Filling it with {fill_value}.")
                df[OF_dof] = np.ones(len(df))*fill_value

        return df


    def calcOutputsFromDF(WT, df, noAcc=False, useTopLoadsFromDF=False):
        """ 
        Given a dataFrame containing time series of DOF
        Compute outputs using OpenFAST Naming Convention
        
        INPUTS: 
         - df: dataframe with time series of degrees of fredom
               For instance df= weio.read('main.outb').toDataFrame()
               Columns Names: 'Time_[s]'
               Q_, QD_, QDD_ ['Sg', 'Sw', 'Hv' ,'R', 'P', 'Y', 'TFA1', 'TSS1', 'Yaw']
         - noAcc: set accelerations to zero
        """
        from welib.tools.tictoc import Timer
        from welib.fast.postpro import ED_TwrGag #, ED_TwrStations, getEDClass
        if len(df)==0:
            raise Exception('No Data in dataframe, make sure you selected a proper time range')

        def towerSectionLoads(twr, F_top_t, M_top_t, kin, gravity):
            nSpan = len(twr.s_span)
            p_ext = np.zeros(nSpan)
            a_struct_t = np.zeros((3,nSpan))
            R_g2t = kin['R_g2t']
            for j in range(nSpan):
                a_struct_t[:,j] = R_g2t.dot(kin['a_Ts'][j,:])
            gravity_vec = np.array((0.,0.,-gravity)) # external acceleration (gravity/earthquake)
            a_ext = R_g2t.dot(gravity_vec)
            # NOTE: assumes that U,V, K have been computed using twr.updateFlexibleKinematics 
            F_sec, M_sec, outDBG =  beamSectionLoads3D(p_ext=p_ext, F_top=F_top_t, M_top=M_top_t, s_span=twr.s_span, m=twr.m, U=twr.U, V=twr.V, K=twr.K, a_struct=a_struct_t, 
                     a_ext=a_ext, corrections=1)
            return F_sec, M_sec

        df = WT._insertOFDOFsInDF(df)

        # --- States
        sq   = [ "Q_Sg_[m]"       , "Q_Sw_[m]"       , "Q_Hv_[m]"       , "Q_R_[rad]"       , "Q_P_[rad]"       , "Q_Y_[rad]"       , "Q_TFA1_[m]"      , "Q_TFA2_[m]"       , "Q_TSS1_[m]"       , "Q_TSS2_[m]"       , "Q_Yaw_[rad]"       , ]
        sqd  = [ "QD_Sg_[m/s]"    , "QD_Sw_[m/s]"    , "QD_Hv_[m/s]"    , "QD_R_[rad/s]"    , "QD_P_[rad/s]"    , "QD_Y_[rad/s]"    , "QD_TFA1_[m/s]"   , "QD_TFA2_[m/s]"    , "QD_TSS1_[m/s]"    , "QD_TSS2_[m/s]"    , "QD_Yaw_[rad/s]"    , ]
        sqdd = [ "QD2_Sg_[m/s^2]" , "QD2_Sw_[m/s^2]" , "QD2_Hv_[m/s^2]" , "QD2_R_[rad/s^2]" , "QD2_P_[rad/s^2]" , "QD2_Y_[rad/s^2]" , "QD2_TFA1_[m/s^2]", "QD2_TFA2_[m/s^2]" , "QD2_TSS1_[m/s^2]" , "QD2_TSS2_[m/s^2]" , "QD2_Yaw_[rad/s^2]" , ]
        missing_dofs = set(sq+sqd+sqdd) - set(df.columns)
        if len(missing_dofs)>0:
            raise Exception(f'Some DOFS are missing from dataframe, implementation error {missing_dofs}')

        # --- DOFs
        Q   = df[sq]
        QD  = df[sqd]
        QDD = df[sqdd]
        # TODO TODO Sort out issue of convention in OpenFAST
        Q['Q_TSS1_[m]']      *= -1
        QD['QD_TSS1_[m/s]']  *= -1
        QDD['QD2_TSS1_[m/s^2]'] *= -1
        DOFNames_Short = ['Sg','Sw','Hv','R','P','Y','TFA1','TFA2','TSS1','TSS2','Yaw']
        Q.columns   = DOFNames_Short
        QD.columns  = DOFNames_Short
        QDD.columns = DOFNames_Short
        if noAcc:
            QDD *=0

        # --- Outputs
        colOut = ['Time_[s]']
        colOut += sq + sqd + sqdd
        # IMU
        colOut += ['NcIMUTVxs','NcIMUTVys','NcIMUTVzs']
        colOut += ['NcIMUTAxs','NcIMUTAys','NcIMUTAzs']
        colOut += ['NcIMURVxs','NcIMURVys','NcIMURVzs']
        colOut += ['NcIMURAxs','NcIMURAys','NcIMURAzs']
        # Tower Top
        colOut+= ['TwrTpTDxi','TwrTpTDyi','TwrTpTDzi']
        # Yaw Brake
        colOut+= ['YawBrTDxp','YawBrTDyp','YawBrTDzp']
        colOut+= ['YawBrTDxt','YawBrTDyt','YawBrTDzt']
        colOut+= ['YawBrTVxp','YawBrTVyp','YawBrTVzp']
        colOut+= ['YawBrTAxp','YawBrTAyp','YawBrTAzp']
        colOut+= ['YawBrRVxp','YawBrRVyp','YawBrRVzp']
        colOut+= ['YawBrRAxp','YawBrRAyp','YawBrRAzp']
        colOut+= ['YawBrFxp','YawBrFyp','YawBrFzp']
        colOut+= ['YawBrMxp','YawBrMyp','YawBrMzp']
        # ED Outputs
        HEDOut, I = ED_TwrGag(WT.ED, addBase=False)
#         hSL = WT.twr.s_span[iSL] # TODO TODO
#         twr_IOut = [
#         iSL = np.argmin(np.abs(hSL-WT.twr.s_span))
#         sT='TwHt{}'.format(iiSL+1)


        for iiSL,hED in enumerate(HEDOut):
            sT='TwHt{}'.format(iiSL+1)
            colOut+=[sT+'TPxi_[m]'  , sT+'TPyi_[m]'  , sT+'TPzi_[m]']
            colOut+=[sT+'TDxt_[m]'  , sT+'TDyt_[m]'  , sT+'TDzt_[m]']
            colOut+=[sT+'RDxt_[deg]', sT+'RDyt_[deg]', sT+'RDzt_[deg]']
            colOut+=[sT+'ALxt_[m/s^2]', sT+'ALyt_[m/s^2]', sT+'ALzt_[m/s^2]']
            colOut+=[sT+'FLxt_[kN]', sT+'FLyt_[kN]', sT+'FLzt_[kN]']
            colOut+=[sT+'MLxt_[kN-m]', sT+'MLyt_[kN-m]', sT+'MLzt_[kN-m]']
        # TODO TODO TODO FIGURE OUT WHY THIS RETURN DTYPE OBJECT
        #dfOut = pd.DataFrame(index=df.index, columns=colOut, dtype=float)
        dfOut = WEIODataFrame(index=df.index, columns=colOut, dtype=float)
        gravity_vec = np.array([0,0,-WT.gravity])

        # --- Calc Output per time step

        with Timer('Time Loop'):
            for it,t in enumerate(df['Time_[s]']):
                if np.mod(it,1000)==0:
                    print(f'Time Loop {it}/{len(df)}')
                # --- Main DOFs
                q   = Q.iloc[it,:].copy()
                qd  = QD.iloc[it,:].copy()
                qdd = QDD.iloc[it,:].copy()
                dfOut.loc[it, sq]   = q.values
                dfOut.loc[it, sqd]  = qd.values
                dfOut.loc[it, sqdd] = qdd.values

                # --------------------------------------------------------------------------------}
                # --- Kinematics 
                # --------------------------------------------------------------------------------{
                # --- Kinematics
                dd = WT.kinematics(q, qd, qdd)
                dfOut.loc[it, 'Time_[s]'] = t
                # TDi includes all platform motions
                dfOut.loc[it, 'TwrTpTDxi'] = dd['u_N_tot'][0] 
                dfOut.loc[it, 'TwrTpTDyi'] = dd['u_N_tot'][1]
                dfOut.loc[it, 'TwrTpTDzi'] = dd['u_N_tot'][2]

                # Alias
                u_N     = dd['u_N']
                v_N     = dd['v_N']
                a_N     = dd['a_N']
                om_N    = dd['omega_n']
                omd_N   = dd['omegad_n']
                u_N_p   = dd['R_g2p'].dot(u_N)
                u_N_t   = dd['R_g2t'].dot(u_N)
                v_N_p   = dd['R_g2p'].dot(v_N)
                a_N_p   = dd['R_g2p'].dot(a_N)
                om_N_p  = dd['R_g2p'].dot(om_N)
                omd_N_p = dd['R_g2p'].dot(omd_N)

                dfOut.loc[it, 'YawBrTDxt'] = u_N_t[0]
                dfOut.loc[it, 'YawBrTDyt'] = u_N_t[1]
                dfOut.loc[it, 'YawBrTDzt'] = u_N_t[2]
                dfOut.loc[it, 'YawBrTDxp'] = u_N_p[0]
                dfOut.loc[it, 'YawBrTDyp'] = u_N_p[1]
                dfOut.loc[it, 'YawBrTDzp'] = u_N_p[2]

                dfOut.loc[it, 'YawBrTVxp'] = v_N_p[0]
                dfOut.loc[it, 'YawBrTVyp'] = v_N_p[1]
                dfOut.loc[it, 'YawBrTVzp'] = v_N_p[2]
                dfOut.loc[it, 'YawBrTAxp'] = a_N_p[0]
                dfOut.loc[it, 'YawBrTAyp'] = a_N_p[1]
                dfOut.loc[it, 'YawBrTAzp'] = a_N_p[2]
                dfOut.loc[it, 'YawBrRVxp'] = om_N_p[0] * 180/np.pi
                dfOut.loc[it, 'YawBrRVyp'] = om_N_p[1] * 180/np.pi
                dfOut.loc[it, 'YawBrRVzp'] = om_N_p[2] * 180/np.pi
                dfOut.loc[it, 'YawBrRAxp'] = omd_N_p[0] * 180/np.pi
                dfOut.loc[it, 'YawBrRAyp'] = omd_N_p[1] * 180/np.pi
                dfOut.loc[it, 'YawBrRAzp'] = omd_N_p[2] * 180/np.pi

                # Alias
                a_IMU     = dd['a_IMU']
                v_IMU     = dd['v_IMU']
                om_IMU    = dd['omega_n']
                omd_IMU   = dd['omegad_n']
                a_IMU_s   = dd['R_g2s'].dot(a_IMU)
                v_IMU_s   = dd['R_g2s'].dot(v_IMU)
                om_IMU_s  = dd['R_g2s'].dot(om_IMU)
                omd_IMU_s = dd['R_g2s'].dot(omd_IMU)

                dfOut.loc[it, 'NcIMUTVxs'] = v_IMU_s[0]
                dfOut.loc[it, 'NcIMUTVys'] = v_IMU_s[1]
                dfOut.loc[it, 'NcIMUTVzs'] = v_IMU_s[2]
                dfOut.loc[it, 'NcIMUTAxs'] = a_IMU_s[0]
                dfOut.loc[it, 'NcIMUTAys'] = a_IMU_s[1]
                dfOut.loc[it, 'NcIMUTAzs'] = a_IMU_s[2]
                dfOut.loc[it, 'NcIMURVxs'] = om_IMU_s[0] * 180/np.pi
                dfOut.loc[it, 'NcIMURVys'] = om_IMU_s[1] * 180/np.pi
                dfOut.loc[it, 'NcIMURVzs'] = om_IMU_s[2] * 180/np.pi
                dfOut.loc[it, 'NcIMURAxs'] = omd_IMU_s[0] * 180/np.pi
                dfOut.loc[it, 'NcIMURAys'] = omd_IMU_s[1] * 180/np.pi
                dfOut.loc[it, 'NcIMURAzs'] = omd_IMU_s[2] * 180/np.pi

                # --- RNA (without Yaw Br) loads
                omd_n       = dd['omegad_n']
                om_n        = dd['omega_n']
                R_g2p       = dd['R_g2p']
                R_g2n       = dd['R_g2n']
                r_Grna      = dd['r_Grna']
                a_Grna      = dd['a_Grna']

                # --------------------------------------------------------------------------------}
                # --- Loads
                # --------------------------------------------------------------------------------{
                rowDF_in = df.iloc[it] # Prescribed loads from input for Hacking only
                Mrna        = WT.RNA_noYawBr.mass
                JGrna       = WT.RNA_noYawBr.masscenter_inertia
                JGrna_g     = (R_g2n.T).dot(JGrna).dot(R_g2n)
                F_Grna_grav = Mrna *gravity_vec
                r_NGrna     = dd['r_NGrna']

                R_N   = Mrna * a_Grna - F_Grna_grav
                tau_N = np.cross(r_NGrna, R_N)
                tau_N += JGrna_g.dot(omd_n)
                tau_N += np.cross(om_n, JGrna_g.dot(om_n))

                # --- Force at N without YawBr Mass (such are "YawBr" sensors..) in global coordinates
                F_N = -R_N
                M_N = -tau_N   #np.cross(r_NGrna, F_Grna_grav)
                if not useTopLoadsFromDF:
                    # Aero force
                    # TODO gen?
                    R_g2s = dd['R_g2s']
                    if 'Fadd_R_xs' in df.keys():
                        Fadd_R_in_g = R_g2s.T.dot((df['Fadd_R_xs'].loc[it],0 ,0))
                        Madd_R_in_g = R_g2s.T.dot((df['Madd_R_xs'].loc[it],0 ,0))
                        r_NR_in_n = WT.rot.pos_global # actually not pos_global but from N
                        r_NR_in_g = R_g2n.T.dot(r_NR_in_n)
                        Madd_R_N = np.cross(r_NR_in_g, Fadd_R_in_g)
                        Fadd_N = Fadd_R_in_g
                        Madd_N = Madd_R_in_g + Madd_R_N*0 # TODO experiment
                        F_N += Fadd_N
                        M_N += Madd_N
#                     else:
#                         raise Exception('Temporary safety')
                F_N_p = R_g2p.dot(F_N)
                M_N_p = R_g2p.dot(M_N)

                dfOut.loc[it, 'YawBrFxp'] = F_N_p[0]/1000
                dfOut.loc[it, 'YawBrFyp'] = F_N_p[1]/1000
                dfOut.loc[it, 'YawBrFzp'] = F_N_p[2]/1000
                dfOut.loc[it, 'YawBrMxp'] = M_N_p[0]/1000
                dfOut.loc[it, 'YawBrMyp'] = M_N_p[1]/1000
                dfOut.loc[it, 'YawBrMzp'] = M_N_p[2]/1000


                if useTopLoadsFromDF:
                    F_N_p = np.array((df['YawBrFxp_[kN]'].loc[it], df['YawBrFyp_[kN]'].loc[it], df['YawBrFzp_[kN]'].loc[it]))*1000
                    M_N_p = np.array((df['YawBrMxp_[kN-m]'].loc[it], df['YawBrMyp_[kN-m]'].loc[it], df['YawBrMzp_[kN-m]'].loc[it]))*1000
                    F_N = (R_g2p.T).dot(F_N_p)
                    M_N = (R_g2p.T).dot(M_N_p)
                
                # Yaw Brake contribution at N
                F_N_YawBr = WT.yawBr.mass * gravity_vec
                F_N += F_N_YawBr

                # --- Top Loads in tower coordinates
                R_g2t = dd['R_g2t']
                F_N_t = R_g2t.dot(F_N)
                M_N_t = R_g2t.dot(M_N)
                TopLoad_t = np.concatenate((F_N_t, M_N_t))

                # --------------------------------------------------------------------------------}
                # ---  Tower Section Loads and Kinematics
                # --------------------------------------------------------------------------------{
                F_sec, M_sec = towerSectionLoads(WT.twr, F_N_t, M_N_t, kin=dd, gravity=WT.gravity)

                dfOut.loc[it, 'TwrBsFxt_[kN]']   = F_sec[0, 0] /1000
                dfOut.loc[it, 'TwrBsFyt_[kN]']   = F_sec[1, 0] /1000
                dfOut.loc[it, 'TwrBsFzt_[kN]']   = F_sec[2, 0] /1000
                dfOut.loc[it, 'TwrBsMxt_[kN-m]'] = M_sec[0, 0] /1000
                dfOut.loc[it, 'TwrBsMyt_[kN-m]'] = M_sec[1, 0] /1000
                dfOut.loc[it, 'TwrBsMzt_[kN-m]'] = M_sec[2, 0] /1000

                for iiSL, hSL in enumerate(HEDOut): # At ED output sections only
                    iSL = np.argmin(np.abs(hSL-WT.twr.s_span)) # TODO precompute
                    hSL = WT.twr.s_span[iSL]
                    sT='TwHt{}'.format(iiSL+1)
                    dfOut.loc[it, sT+'FLxt_[kN]']   = F_sec[0, iSL]/1000
                    dfOut.loc[it, sT+'FLyt_[kN]']   = F_sec[1, iSL]/1000
                    dfOut.loc[it, sT+'FLzt_[kN]']   = F_sec[2, iSL]/1000
                    dfOut.loc[it, sT+'MLxt_[kN-m]'] = M_sec[0, iSL]/1000
                    dfOut.loc[it, sT+'MLyt_[kN-m]'] = M_sec[1, iSL]/1000
                    dfOut.loc[it, sT+'MLzt_[kN-m]'] = M_sec[2, iSL]/1000

                    dfOut.loc[it, sT+'TDxt_[m]'] = dd['u_Ts_in_t'][iSL,0]
                    dfOut.loc[it, sT+'TDyt_[m]'] = dd['u_Ts_in_t'][iSL,1]
                    dfOut.loc[it, sT+'TDzt_[m]'] = dd['u_Ts_in_t'][iSL,2]
                    dfOut.loc[it, sT+'RDxt_[deg]'] = dd['theta_TTs_in_t'][iSL,0]*180/np.pi
                    dfOut.loc[it, sT+'RDyt_[deg]'] = dd['theta_TTs_in_t'][iSL,1]*180/np.pi
                    dfOut.loc[it, sT+'RDzt_[deg]'] = dd['theta_TTs_in_t'][iSL,2]*180/np.pi
                    a_Ts = R_g2t.dot(dd['a_Ts'][iSL])
                    dfOut.loc[it, sT+'ALxt_[m/s^2]'] = a_Ts[0]
                    dfOut.loc[it, sT+'ALyt_[m/s^2]'] = a_Ts[1]
                    dfOut.loc[it, sT+'ALzt_[m/s^2]'] = a_Ts[2]

                    dfOut.loc[it, sT+'TPxi_[m]'] = dd['r_Ts'][iSL,0]
                    dfOut.loc[it, sT+'TPyi_[m]'] = dd['r_Ts'][iSL,1]
                    dfOut.loc[it, sT+'TPzi_[m]'] = dd['r_Ts'][iSL,2]
                    
        # --- Combine Section loads into a dicitonary                    
        sections=None
        return dfOut, sections



# --------------------------------------------------------------------------------}
# --- Helpers 
# --------------------------------------------------------------------------------{
def rigidBlades(blds, hub=None, r_O=[0,0,0]):
    """ return a rigid body for the three blades
    All bodies should be in a similar frame
    """
    blades = blds[0].toRigidBody()
    for B in blds[1:]:
        B_rigid = B.toRigidBody()
        blades = blades.combine(B_rigid, r_O=r_O)
    blades.name='blades'
    return blades


# --------------------------------------------------------------------------------}
# --- Converters 
# --------------------------------------------------------------------------------{
class FASTWindTurbine():
    """ 
    Constructor for Wind tubine class
    """
    def __init__(self, fstFilename=None, main_axis='z', 
                    nSpanTwr=None, twrShapes=None, 
                    nSpanBld=None, bldShapes=None,
                    nSpanSub=None, subShapes=None,
                    fixedShaft=False,
                    algo='', bldStartAtRotorCenter=True,
                    gravity=None,
                    WT=None,
                    verbose=False,
                    # HD options
                    HD_compFile=None, # Nasty HD component files for wave elevation
                    # SD options
                    SD_bOverride = False ,
                    SD_FEM_method= 'full', # full (SubDYn) or cbeam, will affect shape functions
                 ):
        """
        INPUTS:
         - twrShapes: Select shapes to use for tower. If None, twrShapes=[0,1,2,3]

        """
        # --- Main Data
        self.ED          = None
        self.FST         = None
        self.twrFile     = None        # TODO
        self.bldFile     = None
        self.main_axis   = main_axis
        self.fstFilename = fstFilename
        self.pBld        = None
        self.pTwr        = None
        self.algo        = None

        # --- Sanity checks
        if WT is not None and fstFilename is not None:
            raise Exception('Cannot provide both a WT and a fstFilename')

        # --- Store direct inputs
        if WT is None:
            self.WT = WindTurbineStructure()
        else:
            self.WT = WT
        self.WT.algo = algo #<<<<< Setting algo early on
        self.verbose = verbose
        self.twrShapes = twrShapes # will be overriden later
        self.bldShapes = bldShapes # will be overriden later
        self.subShapes = subShapes # will be overriden later

        if fstFilename is None:
            return
        # --- Read fast input files
        self.loadFST(fstFilename) # self.FST, self.ED, self.gravity
        self.setGravity(gravity) # self.FST, self.ED, self.gravity
        self.setupEDGeom()
        self.setupEDHub()
        self.setupEDGen()
        self.setupEDNac()
        self.setupEDBld(shapes=bldShapes, nSpan=nSpanBld, bldStartAtRotorCenter=bldStartAtRotorCenter)
        self.setupEDRot()
        self.setupEDYaw()
        self.setupEDRNA()
        if self.FST['CompSub']>0:
            FAIL('windturbine.py: SubDyn `fnd` not implemented, only ED rigid body platform included.')
            self.setupEDRigidFloat()
        else:
            self.setupEDRigidFloat()
        self.setupEDTwr(shapes=twrShapes, nSpan=nSpanTwr)
        self.setupWTRigid()
        self.setupMAP()
        self.setupEDDOFs()
        
        self.WT.ED = self.ED # TODO 


    def loadFST(self, fstFilename, readlist=None):
        from welib.weio.fast_input_deck import FASTInputDeck

        if readlist is None:
            readlist = ['Fst', 'ED', 'EDtwr', 'EDbld', 'SD', 'HD', 'SS']

        # --- Reading main OpenFAST files
        ext=os.path.splitext(fstFilename)[1]
        if ext.lower()!='.fst':
            raise Exception('FNSB requires a fst file as input')
        self.DCK     = FASTInputDeck(fstFilename, readlist = readlist)
        self.FST     = self.DCK.fst_vt['Fst']
        self.ED      = self.DCK.fst_vt['ElastoDyn']
        self.bldFile = self.DCK.fst_vt['ElastoDynBlade']
        self.twrFile = self.DCK.fst_vt['ElastoDynTower']
        # TODO, MoorDyn, BeamDyn 
        self.SDFile = self.DCK.fst_vt['SubDyn']
        self.HDFile = self.DCK.fst_vt['HydroDyn']
        self.SSFile = self.DCK.fst_vt['SeaState']


    def setGravity(self, gravity=None):
        if gravity is not None:
            self.WT.gravity = gravity
        else:
            try:
                self.WT.gravity = self.FST['gravity']
            except:
                try:
                    self.WT.gravity = self.ED['gravity']
                except:
                    raise Exception('Variable gravity not found in FST file or ED file.')


    def _defaultNSpanTwr(self, nSpan=None, verbose=False, fallback=101):
        nSpanED = self.ED['TwrNodes']
        if nSpan is None:
            if self.WT.algo=='OpenFAST':
                nSpan = nSpanED
                if verbose:
                    print('[INFO] windturbine: Using number of tower nodes ({}) from OpenFAST Input file.'.format(nSpan))
            else:
                nSpan=fallback
                if verbose:
                    print('[INFO] windturbine: Using default of tower nodes ({}).'.format(nSpan))
        else:
            if self.WT.algo=='OpenFAST':
                if nSpan!=nSpanED:
                    INFO('windturbine: Using user-specified number of tower nodes ({}) instead of ED({}).'.format(nSpan, nSpanED))
        return nSpan

    def _defaultNSpanBld(self, nSpan=None, verbose=False, fallback=61):
        nSpanED = self.ED['BldNodes']
        if nSpan is None:
            if self.WT.algo=='OpenFAST':
                nSpan = nSpanED
                if verbose:
                    print('[INFO] windturbine: Using number of blade nodes ({}) from OpenFAST Input file.'.format(nSpan))
            else:
                nSpan=fallback
                if verbose:
                    print('[INFO] windturbine: Using default number of blade nodes ({}).'.format(nSpan))
        else:
            if self.WT.algo=='OpenFAST':
                if nSpan!=nSpanED:
                    WARN('windturbine: Using user-specified number of blade nodes ({}) instead of ED({}).'.format(nSpan, nSpanED))
        return nSpan


    def setupSDInit(self):
        from welib.fast.subdyn import SubDyn
        if self.FST['CompSub']==0:
            raise Exception('Windturbine: SubDyn cannot be initialized, CompSub=0.')
            self.SD=None
            return

        # Mostly to get zBot..
        self.SD = SubDyn(self.SDFile) # TODO TODO
        # We store everything in SD for convenience
        self.SD.graph__ = self.SDFile.toGraph() 
        self.SD.graph__.divideElements(self.SDFile['NDiv'])
        self.SD.graph__.sortNodesBy('z')
        df = self.SD.graph__.nodalDataFrame()
        self.SD.zBot = np.min(df['z'])
        self.SD.zTop = np.max(df['z'])
        self.SD.RayleighCoeff = None
        self.SD.DampMat       = None
        if self.SDFile['GuyanDampMod']==1:
            # Rayleigh Damping
            self.SD.RayleighCoeff=self.SDFile['RayleighDamp']
            #if RayleighCoeff[0]==0:
            #    damp_zeta=omega*RayleighCoeff[1]/2. 
        elif self.SDFile['GuyanDampMod']==2:
            # Full matrix
            self.SD.DampMat = self.SDFile['GuyanDampMatrix']
            self.SD.DampMat = self.SD.DampMat[np.ix_(shapes,shapes)]


    def setupEDGeom(self, zBot=0, bTiltBeforeNac=False, flavor=''):

        ED = self.ED
        WT = self.WT

        if self.main_axis=='x':
           WT.r_EF_inE = np.array([zBot                         ,0,0]) 
           WT.r_ET_inE = np.array([ED['TowerBsHt']              ,0,0]) 
           WT.r_FT_inF = np.array([ED['TowerBsHt']-zBot         ,0,0]) 
           WT.r_TN_inT = np.array([ED['TowerHt']-ED['TowerBsHt'],0,0])

           WT.shaft_tilt   = ED['ShftTilt']*np.pi/180    # NOTE: tilt has wrong orientation in FAST
           WT.blade_cone = -ED['Precone(1)']*np.pi/180

           if bTiltBeforeNac:
               raise NotImplementedError()
               WT.R_NS0 = np.eye(3)
               WT.R_TN0 = R_y(WT.shaft_tilt)
           else:
               WT.R_NS0 = R_y(WT.shaft_tilt)
               WT.R_TN0 = np.eye(3)
               WT.R_NS  = R_y(WT.shaft_tilt)
               WT.r_NGnac_inN = np.array([ED['NacCMzn'],ED['NacCMyn'],ED['NacCMxn']] )
               WT.r_NS_inN    = np.array([ED['Twr2Shft'] ,0,0]) # S on tower axis
           WT.r_SR_inS    = np.array([0,0,ED['OverHang']] ) # S and R 
           WT.r_SGhub_inS = np.array([0,0,ED['OverHang']+ED['HubCM']]   ) # 

        elif self.main_axis=='z':
            # 
            WT.r_EF_inE    = np.array([0,0,zBot  ]) 
            WT.r_EPtfm_inE = np.array([0,0,ED['PtfmRefzt']      ])  # TODO TODO TODO
            WT.r_FT_inF    = np.array([0,0,ED['TowerBsHt']-zBot ]) 
            WT.r_ET_inE    = np.array([0,0,ED['TowerBsHt']      ])  # TODO TODO TODO
            WT.r_TN_inT    = np.array([0,0,ED['TowerHt']-ED['TowerBsHt'] ])

            # Basic geometries for nacelle
            WT.shaft_tilt = -ED['ShftTilt']*np.pi/180  # NOTE: tilt has wrong orientation in FAST
            WT.blade_cone = ED['Precone(1)']*np.pi/180

            if bTiltBeforeNac:
                raise NotImplementedError()
                R_NS0 = np.eye(3)
                R_TN0 = R_y(WT.shaft_tilt)
            else:
                WT.R_NS0 = R_y(WT.shaft_tilt)  # Rotation fromShaft to Nacelle
                WT.R_TN0 = np.eye(3)
                WT.R_NS  = R_y(WT.shaft_tilt)  # Rotation fromShaft to Nacelle
                WT.r_NGnac_inN = np.array([ED['NacCMxn'],ED['NacCMyn'], ED['NacCMzn']    ])                  # Nacelle G in N
                WT.r_NS_inN    = np.array([0             , 0, ED['Twr2Shft']]) # Shaft start in N
            WT.r_SR_inS    = np.array([ED['OverHang'], 0, 0             ]) # Rotor center in S
            WT.r_SGhub_inS = np.array([ED['HubCM']   , 0, 0             ]) + WT.r_SR_inS # Hub G in S



        # --- Common
        WT.r_NR_inN    = WT.r_NS_inN + WT.R_NS.dot(WT.r_SR_inS)       # Rotor center in N
        WT.r_RGhub_inS = - WT.r_SR_inS + WT.r_SGhub_inS


        # --- Monopile
        if flavor=='monopile_is_tower':
            WT.r_ET_inE = WT.r_EF_inE
            WT.r_TN_inT = WT.r_FT_inF+WT.r_TN_inT # assume that F and T are in system E here


        # --- OpenFAST compatibility
        if WT.algo.lower()=='openfast':
            #print('[INFO] YAMS Wind Turbine - Using Algo OpenFAST')
            from welib.fast.elastodyn import rotorParameters, bladeParameters, towerParameters, bladeDerivedParameters, towerDerivedParameters
            pBld = bladeParameters(self.ED.filename)
            self.pBld = bladeDerivedParameters(pBld, inertiaAtBladeRoot=False)
            pRot, pbld, phub = rotorParameters(self.ED.filename, identicalBlades=False)
            RotMass = pRot['RotMass']
            pTwr = towerParameters(self.ED.filename, RotMass=RotMass, gravity=WT.gravity)
            self.pTwr = towerDerivedParameters(pTwr)

    def setupEDHub(self, bHubMass=1, flavor= ''):
        # --- Hub  (defined using point N and nacelle coord as ref)
        ED = self.ED
        WT = self.WT
        M_hub      = ED['HubMass']*bHubMass
        JxxHub_atR = ED['HubIner']*bHubMass
        hub = RigidBody('Hub', M_hub, (JxxHub_atR,0,0), s_OG=WT.r_SGhub_inS, R_b2g=WT.R_NS, s_OP=WT.r_SR_inS, r_O=WT.r_NS_inN)
        WT.hub = hub
#         if flavor=='yams_rec':
#             raise NotImplementedError()

    def setupEDHubGen(self, bHubMass=1, flavor=''):
        """ sft = hub + gen"""
        ED = self.ED
        WT = self.WT
        M_hub      = ED['HubMass']*bHubMass
        JxxHub_atR = (ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2) * bHubMass
        if flavor=='yams_rec':
            # --- Hub
            IR_hub = np.zeros((3,3))
            if self.main_axis=='x':
                IR_hub[2,2] = JxxHub_atR
            elif self.main_axis=='z':
                IR_hub[0,0] = JxxHub_atR
            IG_hub = translateInertiaMatrix(I_A=IR_hub, Mass=M_hub, r_BG=np.array([0,0,0]), r_AG = WT.r_RGhub_inS)
            hubgen = YAMSRecRigidBody('ShaftHubGen', M_hub, IG_hub, WT.r_SGhub_inS)
        else:
            raise NotImplementedError()
        WT.hubgen = hubgen


    def setupEDGen(self, flavor=''):
        # --- Generator (Low speed shaft) (defined using point N and nacelle coord as ref)
        ED = self.ED
        WT = self.WT
        Jp = ED['GenIner']*ED['GBRatio']**2
        if flavor=='yams_rec':
            IR_gen = np.zeros((3,3))
            if self.main_axis=='x':
                IR_gen[2,2] = Jp
            elif self.main_axis=='z':
                IR_gen[0,0] = Jp
            # Generator has no mass, no need to translate inertia
            #IG_gen = translateInertiaMatrix(I_A=IR_gen, Mass=0, r_BG=np.array([0,0,0]), r_AG = WT.r_RGhub_inS)
            gen = YAMSRecRigidBody('Gen', 0, IR_gen, WT.r_SGhub_inS)
        else:
            gen = RigidBody('Gen', 0, (Jp,0,0), s_OG=[0,0,0], R_b2g=WT.R_NS, r_O=WT.r_NS_inN) 
        WT.gen = gen

    def setupEDNac(self, bNacMass=1, flavor=''):
        # --- Nacelle (defined using point N and nacelle coord as ref)
        ED = self.ED
        WT = self.WT
        M_nac       = ED['NacMass'] * bNacMass
        JyyNac_atN  = ED['NacYIner'] *bNacMass # Inertia of nacelle at N in N

        if flavor=='yams_rec':
            I0_nac = np.zeros((3,3)) 
            # ElastoDyn NacYIner is the inertia about nacelle local y-axis.
            I0_nac[1,1] = ED['NacYIner']
            I0_nac = I0_nac * bNacMass
            IG_nac = translateInertiaMatrixToCOG(I0_nac, M_nac, WT.r_NGnac_inN)
            # Nacelle Body
            nac = YAMSRecRigidBody('Nacelle', M_nac, IG_nac, WT.r_NGnac_inN)

        else:
            nac = RigidBody('Nac', M_nac, (0,JyyNac_atN,0), WT.r_NGnac_inN, s_OP = [0,0,0])
        WT.nac = nac

    def setupEDBld(self, shapes=None, nSpan=None, 
                   spanFrom0=False, bBldMass=1, bldStartAtRotorCenter=True,
                   massExpected=None,
                   flavor=''
                   ):
        # TODO TODO TODO  Harmonize with fast.elastodyn when algo is OpenFAST
        ED = self.ED
        WT = self.WT
        # --- Default arguments 
        if shapes is None:  # if not provided we respect the ED file
            shapes=[]
            if ED['FlapDOF1']:
                shapes+=[0]
            if ED['FlapDOF2']:
                shapes+=[1]
            if ED['EdgeDOF']:
                shapes+=[2]
        self.bldShapes = shapes
        nSpan = self._defaultNSpanBld(nSpan, verbose=self.verbose)

        m    = self.bldFile['BldProp'][:,3]
        jxxG=0*m
        #if algo=='OpenFAST':
        #    jxxG=0*m
        #else:
        #    jxxG = 0*m     # NOTE: unknown
        #    print('>>> windturbine.py: TODO: using unknown jxxG')
        nB = ED['NumBl']
        bld=np.zeros(nB,dtype=object)

        if flavor=='yams_rec':
            # ----------- YAMS REC------------------------------------------------------------
            bld[0] = YAMSRecFASTBeamBody('blade', ED, self.bldFile,
                                         Mtop=0, 
                                         shapes=shapes, nSpan=nSpan, 
                                         main_axis=self.main_axis, 
                                         spanFrom0=spanFrom0, algo=WT.algo,
                                         massExpected=massExpected, 
                                         gravity=WT.gravity) # NOTE: legacy spanfrom0

        else:
            # ----------- GENERIC BODY -------------------------------------------------------
            bld[0] = FASTBeamBody(ED, self.bldFile, 
                                  Mtop=0, 
                                  jxxG=jxxG, 
                                  shapes=shapes, nSpan=nSpan, 
                                  main_axis=self.main_axis, 
                                  spanFrom0=spanFrom0, bldStartAtRotorCenter=bldStartAtRotorCenter, algo=WT.algo,
                                  gravity=WT.gravity
                                  ) 

        if WT.algo.lower()=='openfast':
            # Overwrite blade generalized matrices with OpenFAST-compatible values.
            M = self.pBld['BldMass']
            mdCM = np.asarray(self.pBld['mdCM']).ravel()
            J = np.asarray(self.pBld['J'])

            MM_of = np.zeros_like(bld[0].MM)
            MM_of[0,0] = M
            MM_of[1,1] = M
            MM_of[2,2] = M
            MM_of[0:3,3:6] = -np.array([ [0, -mdCM[2], mdCM[1]], [mdCM[2], 0, -mdCM[0]], [-mdCM[1], mdCM[0], 0], ])
            MM_of[3:6,0:3] = MM_of[0:3,3:6].T
            MM_of[3:6,3:6] = J
            if len(shapes)>0:
                I = np.asarray(shapes, dtype=int)
                MM_of[0:3,6:] = self.pBld['Ct'][I,:].T
                MM_of[3:6,6:] = self.pBld['Cr'][I,:].T
                MM_of[6:,0:3] = MM_of[0:3,6:].T
                MM_of[6:,3:6] = MM_of[3:6,6:].T
                MM_of[6:,6:] = self.pBld['Me'][np.ix_(I, I)]

                bld[0].MM[:,:] = MM_of
                bld[0].KK0[6:,6:] = self.pBld['Ke0'][np.ix_(I, I)] # NOTE: OpenFAST uses Ke0
                bld[0].KK [6:,6:] = self.pBld['Ke'] [np.ix_(I, I)]
                bld[0].DD [6:,6:] = self.pBld['De'] [np.ix_(I, I)]

        # --- Common code
        bld[0].MM *=bBldMass
        # Copy blades
        for iB in range(nB-1):
            bld[iB+1] = copy.deepcopy(bld[0])
        # Set Rotation matrices from shaft to each blade
        for iB,B in enumerate(bld):
            B.name='bld'+str(iB+1)
            psi_B= -iB*2*np.pi/len(bld) 
            if self.main_axis=='x':
                R_SB = R_z(0*np.pi + psi_B) # TODO psi offset and psi0
            elif self.main_axis=='z':
                R_SB = R_x(0*np.pi + psi_B) # TODO psi0
            # IMPORTANT FOR RNA set R_b2g
            if bldStartAtRotorCenter:
                R_SB = np.dot(R_SB, R_y(ED['PreCone({})'.format(iB+1)]*np.pi/180)) # blade2shaft
                B.R_b2g= R_SB
            else:
                print('>>>> TODO TODO TODO Wind Turbine R_SB')
                R_SB = np.dot(R_SB, R_y(ED['PreCone({})'.format(iB+1)]*np.pi/180)) # blade2shaft
                B.R_b2g= R_SB
        WT.bld = bld

        # --- Blades (with origin R, using N as "global" ref)
        blds_rigid = rigidBlades(bld, r_O = [0,0,0])
        blds_rigid.pos_global = WT.r_NR_inN
        blds_rigid.R_b2g      = WT.R_NS
        WT.blds_rigid = blds_rigid


    def setupEDRot(self):
        ED = self.ED
        WT = self.WT
        if WT.hub is None or WT.gen is None or WT.blds_rigid is None:
            raise Exception('Setup hub and gen first')
        # --- Rotor = Hub + Blades (with origin R, using N as global ref)
        rot = WT.blds_rigid.combine(WT.hub, R_b2g=WT.R_NS, r_O=WT.blds_rigid.pos_global)
        rot.name='rotor'
        rotgen = rot.combine(WT.gen, R_b2g=WT.R_NS, r_O=WT.blds_rigid.pos_global)
        #print(rotgen)
        WT.rot        = rot        # origin at R, rigid body bld+hub
        WT.rotgen     = rotgen     # origin at R, rigid body bld+hub+genLSS

    def setupEDYaw(self, flavor=''):
        ED = self.ED
        WT = self.WT
        #--- Yaw bearing, at tower top
        M_yawBr = ED['YawBrMass']
        if flavor=='yams_rec':
            WT.yawBr = YAMSRecRigidBody('YawBearing',M_yawBr,(0,0,0),(0,0,0));
            if M_yawBr>0:
                print('[WARN] TODO YAW BEARING MASS NOT FULLY IMPLEMENTED IN TNSB')
        else:
            WT.yawBr = RigidBody('YawBr', M_yawBr, J=(0,0,0), s_OG=(0,0,0))

    def setupEDRNA(self):
        ED = self.ED
        WT = self.WT
        if WT.hub is None or WT.gen is None or WT.rot is None:
            raise Exception('Setup hub and gen first')
        # --- RNA 
        RNA = WT.rot.combine(WT.gen).combine(WT.nac,r_O=[0,0,0]).combine(WT.yawBr, r_O=[0,0,0])
        RNA.name='RNA'

        # --- RNA without YawBr Mass
        RNA_noYawBr = WT.rot.combine(WT.gen).combine(WT.nac,r_O=[0,0,0])

        WT.RNA = RNA
        WT.RNA_noYawBr = RNA_noYawBr

    def setupEDRigidFloat(self, DOFs=None, flavor=''):
        """ 
        - DOFs : e.g. [0,4] for surge and pitch
        """
        ED = self.ED
        WT = self.WT
        dofs = ['PtfmSgDOF', 'PtfmSwDOF', 'PtfmHvDOF', 'PtfmRDOF', 'PtfmPDOF', 'PtfmYDOF']
        if DOFs is None:
            self.subShapes = [i for i, dof in enumerate(dofs) if ED[dof]]
        else:
            self.subShapes = DOFs
        # --- Fnd (defined wrt ground/MSL "E")
        # print(FST.keys())
        M_fnd = ED['PtfmMass']
        r_EGfnd_inF = np.array([ED['PtfmCMxt'],ED['PtfmCMyt'],ED['PtfmCMzt']])
        r_EPtfm_inF    = np.array([0             ,0             ,ED['PtfmRefzt']]) # TODO, this is wrong
        r_PtfmGfnd_inF = -r_EPtfm_inF + r_EGfnd_inF
        if flavor=='yams_rec':
            raise NotImplementedError()
        else:
            WT.fnd = RigidBody('fnd', M_fnd, (ED['PtfmRIner'], ED['PtfmPIner'], ED['PtfmYIner']), s_OG=r_PtfmGfnd_inF, r_O=r_EPtfm_inF) 
        
    def setupEDTwr(self, shapes=None, nSpan=None, 
                   flavor='', 
                   bAxialCorr=False, bStiffening=True
                   ):
        # TODO TODO TODO  Harmonize with fast.elastodyn when algo is OpenFAST
        ED = self.ED
        WT = self.WT
        # --- Twr
        if shapes is None:  # If not provided we respect the ED setttings
            shapes=[] 
            # TODO WATCH OUT ORDER
            if ED['TwFADOF1']:
                shapes+=[0]
            if ED['TwFADOF2']:
                shapes+=[1]
            if ED['TwSSDOF1']:
                shapes+=[2]
            if ED['TwSSDOF2']:
                shapes+=[3]
        self.twrShapes = shapes
        nSpan = self._defaultNSpanTwr(nSpan, verbose=self.verbose)


        if flavor=='yams_rec':
            # Tower Body
            twr = YAMSRecFASTBeamBody('tower', ED, self.twrFile, 
                                      Mtop=WT.RNA.mass, 
                                      shapes=shapes, nSpan=nSpan, 
                                      main_axis=self.main_axis, 
                                      bStiffening=bStiffening, algo=WT.algo,
                                      gravity=WT.gravity
                                      )
            if WT.algo.lower()=='openfast':
                WARN('Windturbine: applying targeted OpenFAST tower override (Kg_SW only) for yams_rec parity.')

                M = self.pTwr['TwrMass']
                mdCM = np.asarray(self.pTwr['mdCM']).ravel()
                J = np.asarray(self.pTwr['J'])
                MM_OF = np.zeros_like(twr.MM)
                MM_OF[0,0] = M
                MM_OF[1,1] = M
                MM_OF[2,2] = M
                MM_OF[0:3,3:6] = -np.array([ [0, -mdCM[2], mdCM[1]], [mdCM[2], 0, -mdCM[0]], [-mdCM[1], mdCM[0], 0], ])
                MM_OF[3:6,0:3] = MM_OF[0:3,3:6].T
                MM_OF[3:6,3:6] = J
                if len(shapes)>0:
                    I = np.asarray(shapes, dtype=int)
                    MM_OF[0:3,6:] = self.pTwr['Ct'][I,:].T
                    MM_OF[3:6,6:] = self.pTwr['Cr'][I,:].T
                    MM_OF[6:,0:3] = MM_OF[0:3,6:].T
                    MM_OF[6:,3:6] = MM_OF[3:6,6:].T
                    MM_OF[6:,6:] = self.pBld['Me'][np.ix_(I, I)]
#                     twr.MM[:,:] = MM_OF

                #twr.MM[0,0]         = self.pTwr['TwrMass']
                #twr.MM[1,1]         = self.pTwr['TwrMass']
                #twr.MM[2,2]         = self.pTwr['TwrMass']
#                 #twr.MM      [6:,6:] = self.pTwr['Me']   [np.ix_(shapes,shapes)] # TODO TODO a bit too strong for second mode
                twr.KK      [6:,6:] = self.pTwr['Ke']   [np.ix_(shapes,shapes)]
                twr.KK0     [6:,6:] = self.pTwr['Ke0']  [np.ix_(shapes,shapes)]
                twr.KKg_self[6:,6:] = self.pTwr['Kg_SW'][np.ix_(shapes,shapes)]
                twr.KKg_Mtop[6:,6:] = self.pTwr['Kg_TM'][np.ix_(shapes,shapes)]
                #twr.KK  = twr.KK0 + twr.KKg
                twr.DD      [6:,6:] = self.pTwr['De']   [np.ix_(shapes,shapes)]
        else:
            twr = FASTBeamBody(ED, self.twrFile, 
                               Mtop=WT.RNA.mass, 
                               shapes=shapes, nSpan=nSpan, 
                               main_axis=self.main_axis, 
                               bAxialCorr=bAxialCorr, bStiffening=bStiffening, algo=WT.algo, 
                               gravity=WT.gravity
                               ) 

        if flavor!='yams_rec':
            # TODO impose this always?
            if WT.algo.lower()=='openfast':
                twr.MM[0,0]         = self.pTwr['TwrMass']
                twr.MM[1,1]         = self.pTwr['TwrMass']
                twr.MM[2,2]         = self.pTwr['TwrMass']
                twr.MM      [6:,6:] = self.pTwr['Me']   [np.ix_(shapes,shapes)]
                twr.KK      [6:,6:] = self.pTwr['Ke']   [np.ix_(shapes,shapes)]
                twr.KK0     [6:,6:] = self.pTwr['Ke0']  [np.ix_(shapes,shapes)]
                twr.KKg_self[6:,6:] = self.pTwr['Kg_SW'][np.ix_(shapes,shapes)]
                twr.KKg_Mtop[6:,6:] = self.pTwr['Kg_TM'][np.ix_(shapes,shapes)]
                twr.DD      [6:,6:] = self.pTwr['De']   [np.ix_(shapes,shapes)]

        twr_rigid  = twr.toRigidBody()
        twr_rigid.pos_global = WT.r_ET_inE

        WT.twr = twr
        WT.twr_rigid = twr_rigid



    @require_attrs(['SD'], 'call setupSDInit first')
    def setupSD(self, Mtop=0, shapes=None, nSpan=None, 
                bStiffening=True, bCI=True, bOverride=True, # Algo options
                FEM_method='cbeam',
                ):
        """ 
        INPUTS:
           -bOverride :  Override some MM, KK values based on internal FEM CB values, closer to SubDyn
           -bCI       :  Turn on or off the concentrated masses
           -FEM_method:  'cbeam' simplified continuous beam FEM
                         'full'  similar to SubDyn, uses graph, recommended
        """
        # --- Default arguments
        ED = self.ED
        if shapes is None:
            shapes=[] 
            dofs = ['PtfmSgDOF', 'PtfmSwDOF', 'PtfmHvDOF', 'PtfmRDOF', 'PtfmPDOF', 'PtfmYDOF']
            shapes = [i for i, dof in enumerate(dofs) if ED[dof]]
            # TODO CB
            self.subShapes = shapes

        CI = None
        if bCI:
            CI = self.SD.concentrated_masses
        else:
            WARN('Concentrated inertia for SubDyn turned off!')

        fnd = YAMSRecFASTBeamBody('substructure', self.ED, self.SD, Mtop=Mtop, shapes=shapes, nSpan=nSpan, 
                                  main_axis=self.main_axis, bStiffening=bStiffening, gravity=self.WT.gravity,
                                  FEM_method=FEM_method,
                                  concentrated_inertias=CI) # TODO, we could remove that to avoid double counting

        # Optional exact SubDyn reduced-matrix matching for selected Guyan coordinates.
        # This bypasses GMBeam-integrated modal MM/KK for the foundation flexible block.

#         printMat('KKg', fnd.KKg[6:,6:])
#         printMat('KKg_self', fnd.KKg_self[6:,6:])
#         printMat('KKg_Mtop', fnd.KKg_Mtop[6:,6:])
#         printMat('KKg_rot', fnd.KKg_rot[6:,6:])


        #if self.WT.algo=='OpenFAST' and bOverride:
        if bOverride:
            if self.SD._FEM is None or self.SD._FEM.MM_CB is None or self.SD._FEM.KK_CB is None:
                raise Exception('SubDyn reduced matrices not available. Ensure SD.init/applyCB was run before override.')
            if shapes is None:
                raise Exception('For override, `shapes` should be be provided')

            WARN('Windturbine: OVERRIDDING SubDyn values with FEM M_CB and K_CB computation')
            if len(shapes)>0:
                I = [int(i) for i in shapes]
                MM_CB = self.SD._FEM.MM_CB
                KK_CB = self.SD._FEM.KK_CB
                if np.max(I) >= MM_CB.shape[0] or np.max(I) >= KK_CB.shape[0]:
                    raise Exception('Shape index outside SubDyn reduced matrix size')

                MM_sel = MM_CB[np.ix_(I, I)].copy()
                KK_sel = KK_CB[np.ix_(I, I)].copy()

                fnd.MM[6:,6:] = MM_sel
                fnd.KK0[6:,6:] = KK_sel

#             fnd.KKg_self[6:,6:] = 0
#             fnd.KKg_Mtop[6:,6:] = 0
#             fnd.KKg_rot[6:,6:] = 0

                fnd.KKg = fnd.KKg_self + fnd.KKg_Mtop + fnd.KKg_rot  
                fnd.KK  = fnd.KK0 + fnd.KKg

                if self.SD.RayleighCoeff is not None:
                    fnd.DD[6:,6:] = MM_sel*self.SD.RayleighCoeff[0] + KK_sel*self.SD.RayleighCoeff[1]

        #print(Fnd)
        #print('Fnd MM\n',Fnd.MM[6:,6:])
        #print('Fnd KK\n',Fnd.KK[6:,6:])
        # HACK here because doesn't handle this for now
        if self.SDFile['GuyanDampMod']==1:
            fnd.DD[6:,6:] = fnd.MM[6:,6:]*self.SD.RayleighCoeff[0] + fnd.KK[6:,6:]*self.SD.RayleighCoeff[1] 
        self.WT.fnd = fnd


    def setupWTRigid(self):
        WT = self.WT
        if WT.fnd is None or WT.twr_rigid is None or WT.RNA is None:
            raise Exception('Setup fnd and twr_rigid and RNA first')
        # --- Full WT rigid body equivalent, with point T as ref
        RNAb = copy.deepcopy(WT.RNA)
        RNAb.pos_global = WT.r_TN_inT+WT.r_ET_inE
        RNAb.R_b2g      = np.eye(3)
        WT_rigid = RNAb.combine(WT.twr_rigid, r_O=WT.r_ET_inE).combine(WT.fnd, r_O=WT.r_ET_inE) # TODO TODO TODO T or Ptfm
        #WT_rigid = RNAb.combine(twr_rigid, r_O=r_EPtfm_inE).combine(fnd, r_O=r_EPtfm_inE) # TODO T or Ptfm
        WT.WT_rigid = WT_rigid

    def setupMAP(self):
        WT = self.WT
        if WT.fnd is None:
            raise Exception('Setup fnd and twr_rigid first')
        # --- Moorings
        MAP = None
        K_Moor=np.zeros((6,6))
        if self.FST['CompMooring']==1:
            from welib.moor.mappp import Map
    #         try:
            MAP = Map(self.fstFilename)
            #zRef  = twr.pos_global[2]  
            zRef  = WT.fnd.pos_global[2]  
            K_Moor,_ = MAP.stiffness_matrix(epsilon=1e-2, point=(0,0,zRef))
            #K_Moor2,_ = MAP.stiffness_matrix(epsilon=1e-2, point=(0,0,0))
    #         except:
    #             print('YAMS Wind Turbine: problem loading MAP model (only supported on windows)')
        elif self.FST['CompMooring']==2:
            print('YAMS Wind Turbine: TODO MoorDyn')
        WT.MAP    = MAP
        WT.K_Moor = K_Moor     # HACK..

    def setupEDDOFs(self):
        # NOTE: may be overriden by setActiveDOFs to adapt to a given model
        ED = self.ED
        # --- Degrees of freedom
        DOFs=[]
        if hasattr(self, 'SD'):
            if self.SD is not None:
                if self.SD.File['Nmodes']>0:
                    #print('TODO, Need to figure out DOF order with CB')
                    # TODO SubDyn, let's figure out order...
                    for iCB in range(int(self.SD.File['Nmodes'])):
                        DOFs+=[{'name': f'CB{iCB+1}', 'active':True, 'q0':0, 'qd0':0, 'q_channel': f'QCB{iCB+1}_[-]' , 'qd_channel':f'QDCB{iCB+1}_[-]','qdd_channel':f'QD2_CB{iCB+1}_[-]'}]

        # TODO TODO TODO handle alias Q_Sg
        DOFs+=[{'name':'x'      , 'active':ED['PtfmSgDOF'], 'q0': ED['PtfmSurge']  , 'qd0':0 , 'q_channel':'Q_Sg_[m]' , 'qd_channel':'QD_Sg_[m/s]','qdd_channel':'QD2_Sg_[m/s^2]'}]
        DOFs+=[{'name':'y'      , 'active':ED['PtfmSwDOF'], 'q0': ED['PtfmSway']   , 'qd0':0 , 'q_channel':'Q_Sw_[m]' , 'qd_channel':'QD_Sw_[m/s]','qdd_channel':'QD2_Sw_[m/s^2]'}]
        DOFs+=[{'name':'z'      , 'active':ED['PtfmHvDOF'], 'q0': ED['PtfmHeave']  , 'qd0':0 , 'q_channel':'Q_Hv_[m]' , 'qd_channel':'QD_Hv_[m/s]','qdd_channel':'QD2_Hv_[m/s^2]'}]

        # TODO TODO TODO issue here with deg and rad
        DOFs+=[{'name':'phi_x' , 'active':ED['PtfmRDOF'] , 'q0': ED['PtfmRoll']*np.pi/180  , 'qd0':0 , 'q_channel':'Q_R_[rad]'  , 'qd_channel':'QD_R_[rad/s]', 'qdd_channel':'QD2_R_[rad/s^2]'}]
        DOFs+=[{'name':'phi_y' , 'active':ED['PtfmPDOF'] , 'q0': ED['PtfmPitch']*np.pi/180 , 'qd0':0 , 'q_channel':'Q_P_[rad]' , 'qd_channel':'QD_P_[rad/s]', 'qdd_channel':'QD2_P_[rad/s^2]'}]
        DOFs+=[{'name':'phi_z' , 'active':ED['PtfmYDOF'] , 'q0': ED['PtfmYaw']*np.pi/180   , 'qd0':0 , 'q_channel':'Q_Y_[rad]'   , 'qd_channel':'QD_Y_[rad/s]', 'qdd_channel':'QD2_Y_[rad/s^2]'}]

        DOFs+=[{'name':'q_FA1'  , 'active':ED['TwFADOF1'] , 'q0': ED['TTDspFA']  , 'qd0':0 , 'q_channel':'Q_TFA1_[m]', 'qd_channel':'QD_TFA1_[m/s]', 'qdd_channel':'QD2_TFA1_[m/s^2]'}]
        DOFs+=[{'name':'q_SS1'  , 'active':ED['TwSSDOF1'] , 'q0': ED['TTDspSS']  , 'qd0':0 , 'q_channel':'Q_TSS1_[m]', 'qd_channel':'QD_TSS1_[m/s]', 'qdd_channel':'QD2_TSS1_[m/s^2]'}]
        DOFs+=[{'name':'q_FA2'  , 'active':ED['TwFADOF2'] , 'q0': ED['TTDspFA']  , 'qd0':0 , 'q_channel':'Q_TFA2_[m]', 'qd_channel':'QD_TFA2_[m/s]', 'qdd_channel':'QD2_TFA2_[m/s^2]'}]
        DOFs+=[{'name':'q_SS2'  , 'active':ED['TwSSDOF2'] , 'q0': ED['TTDspSS']  , 'qd0':0 , 'q_channel':'Q_TSS2_[m]', 'qd_channel':'QD_TSS2_[m/s]', 'qdd_channel':'QD2_TSS2_[m/s^2]'}]

        DOFs+=[{'name':'theta_y','active':ED['YawDOF']  , 'q0': ED['NacYaw']*np.pi/180   , 'qd0':0 ,          'q_channel':'Q_Yaw_[rad]' , 'qd_channel':'QD_Yaw_[rad/s]', 'qdd_channel':'QD2_Yaw_[rad/s^2]'}]
        DOFs+=[{'name':'psi'    ,'active':ED['GenDOF']  , 'q0': ED['Azimuth']*np.pi/180  , 'qd0':ED['RotSpeed']*2*np.pi/60 , 'q_channel':'Q_GeAz_[rad]', 'qd_channel':'QD_GeAz_[rad/s]', 'qdd_channel': 'QD2_GeAz_[rad/s^2]'}]

        DOFs+=[{'name':'nu'     ,'active':ED['DrTrDOF'] , 'q0': 0  , 'qd0':0 , 'q_channel':'Q_DrTr_[rad]', 'qd_channel':'QD_DrTr_[rad/s]', 'qdd_channel':'QD2_DrTr_[rad/s^2]'}]

        # 
        for ib in np.arange(ED['NumBl']):
            B=str(ib+1)
            DOFs+=[{'name':'q_B{}F1'.format(B), 'active':ED['FlapDOF1'] , 'q0': ED['OOPDefl'], 'qd0':0, 'q_channel':'Q_B{}F1_[m]'.format(B), 'qd_channel':'QD_B{}F1_[m/s]'.format(B), 'qdd_channel':'QD2_B{}F1_[m/s^2]'.format(B)}]
            DOFs+=[{'name':'q_B{}E1'.format(B), 'active':ED['FlapDOF2'] , 'q0': ED['OOPDefl'], 'qd0':0, 'q_channel':'Q_B{}F2_[m]'.format(B), 'qd_channel':'QD_B{}E1_[m/s]'.format(B), 'qdd_channel':'QD2_B{}E1_[m/s^2]'.format(B)}]
            DOFs+=[{'name':'q_B{}E1'.format(B), 'active':ED['EdgeDOF']  , 'q0': ED['IPDefl'] , 'qd0':0, 'q_channel':'Q_B{}E1_[m]'.format(B), 'qd_channel':'QD_B{}E1_[m/s]'.format(B), 'qdd_channel':'QD2_B{}E1_[m/s^2]'.format(B)}]
        self.WT.DOF = DOFs

    def setActiveDOFs(self, fixedShaft=False, shapes_sub=None, shapes_twr=None, shapes_bld=None, verbose=False):
        if shapes_sub is None:
            shapes_sub =[]
        # Override based on model
        SUB_NAMES =  ['x', 'y', 'z', 'phi_x', 'phi_y', 'phi_z', 'CB1', 'CB2', 'CB3', 'CB4', 'CB5', 'CB6', 'CB7', 'CB8', 'CB9', 'CB10'] # Ptfm
        TWR_NAMES =  ['q_FA1', 'q_FA2', 'q_SS1', 'q_SS2'] # Twr
        BLD_NAMES =  ['q_B{}F1', 'q_B{}E1', 'q_B{}E2'] # Twr

        NAMEOFF=[]
        NAMEOFF += [SUB_NAMES[i] for i in range(len(SUB_NAMES)) if i not in shapes_sub]
        NAMEOFF += [TWR_NAMES[i] for i in range(4) if i not in shapes_twr]
        NAMEOFF += ['theta_y'] # Yaw
        NAMEOFF += ['nu'] # Shaft torsion
        if fixedShaft is None:
            pass
        elif fixedShaft:
            NAMEOFF += ['psi'] # Shaft torsion
        for iB in range(3):
            NAMEOFF += [BLD_NAMES[i].format(iB+1) for i in range(3) if i not in shapes_bld]

        for dof in self.WT.DOF:
            if dof['name'] in NAMEOFF:
                if dof['active']:
                    if verbose:
                        print('Deactivating {:10s} ({:20s}) eventhough it was active in ED'.format(dof['name'], dof['q_channel']))
                    dof['active']=False
            else:
                if not dof['active']:
                    if verbose:
                        print('Activating   {:10s} ({:20s}) eventhough it was inactive in ED'.format(dof['name'], dof['q_channel']))
                    dof['active']=True



    def setupSeaState(self, compFile=None):
        if self.verbose:
            INFO('Setting up SeaState')
        pSS = None

        # --- Return None if not active
        if 'CompSeaSt' not in self.FST:
            #WARN('Not supporting sea state for old OpenFAST input files')
            self.WT.pSS = None # Store data in WT class
            return pSS
        if self.FST['CompSeaSt']==0:
            #WARN('Sea State is None')
            self.WT.pSS = None # Store data in WT class
            return pSS

        # --- Setting main data from input file
        pSS = dict() # TODO might deserve an object
        self.WT.pSS = pSS # Store data in WT class

        # --- TODO get rid of these
        try:
            self.WT.WtrDens   = self.FST['WtrDens']
            self.WT.WtrDpth   = self.FST['WtrDpth']
        except:
            # Legacy
            self.WT.WtrDens   = self.HDFile['WtrDens']
            self.WT.WtrDpth   = self.HDFile['WtrDpth']

        pSS['rho']        = self.WT.WtrDens      # [kg/m3]
        pSS['WaterDepth'] = self.WT.WtrDpth      # [m]


        if compFile is not None:
            WARN('Not recommended to set component at this stage')
            self.WT.SS_setComponents(compFile)

    def setupHydro(self, hydroShape=None):
        if self.verbose:
            INFO('Setting up Hydro')

        pHD = None
        if 'CompSeaSt' not in self.FST:
            WARN('Not supporting hydro for old OpenFAST input files')
            self.WT.pHD = None # Store data in WT class
            return 
        # --- Return None if not active
        if self.FST['CompHydro']==0:
            #WARN('Hydro is None')
            self.pHD = pHD # Store data in class
            return

        # --- Safety checks
        if self.WT.pSS is None:
            raise Exception('windturbine: setupHydro: Cannot setup Hydro, call setupSeaState first.')
        if self.WT.fnd is None:
            raise Exception('windturbine: setupHydro: Cannot setup Hydro, need fnd setup.')

        # --- Setting main data from input file
        pHD = dict() # TODO might deserve an object
        self.WT.pHD = pHD # Store data in class

        # --- TODO get rid of these
        self.WT.Hydro     = self.FST['CompHydro']>0 # FST['CompSeaSt']>0 and 
        self.WT.HD        = self.HDFile
        # TODO get rid of these
        try:
            self.WT.WtrDens   = self.FST['WtrDens']
            self.WT.WtrDpth   = self.FST['WtrDpth']
        except:
            # Legacy
            self.WT.WtrDens   = self.HDFile['WtrDens']
            self.WT.WtrDpth   = self.HDFile['WtrDpth']

        # --- Hydro floaters
        if self.WT.fnd.SD is None:
            FAIL('Setting up hydro with out SD not supported yet')
            return pHD

        # --- Hydro Monopile specific
        if not hasattr(self.WT.fnd, 'PhiU'):
            raise Exception('For hydro monopile setup we need a shape function')


        zBeam = np.array(sorted(self.WT.fnd.SD.pointsMN['z'].unique()))
        zDepth = self.WT.fnd.s_span - self.WT.pSS['WaterDepth'] # For check
        zDepth = zBeam
        HD = self.HDFile
        try:
            cprop = HD.getTab('SectionPropCyl')
        except:
            cprop = HD.getTab('SectionProp')
        D_ = cprop['PropD'].values[0]
        sprop = HD.getTab('SmplPropCyl')
        Cp_ = sprop['SimplCp'].values[0]
        Cd_ = sprop['SimplCd'].values[0]
        Ca_ = sprop['SimplCa'].values[0]
        print(f'HD props: Cp={Cp_} Cd={Cd_} Ca={Ca_} D={D_} ')
        PlaceHolderOnes = np.ones([len(zDepth),1])
        pHD['zDepth']   = zDepth
        pHD['D']        = D_  * PlaceHolderOnes
        pHD['Cd']       = Cd_ * PlaceHolderOnes
        pHD['CM']       = (Ca_+Cp_) *PlaceHolderOnes
        pHD['Ca']       = Ca_ * PlaceHolderOnes
        pHD['Cp']       = Cp_ * PlaceHolderOnes
        pHD['m_hydro']  = self.WT.pSS['rho'] * np.pi / 4 * pHD['D'].ravel()**2 * pHD['Ca'].ravel()
            
        # Generalized hydro mass matrix 
        GM_hydro = np.zeros((len(self.WT.fnd.PhiU),len(self.WT.fnd.PhiU)))
        bWet = zDepth<=0
        for i, phi in enumerate(self.WT.fnd.PhiU):
            phi_x = phi[0,:]
            GM_hydro[i,i] = np.trapezoid(pHD['m_hydro'][bWet] * phi_x[bWet]**2, zDepth[bWet])

        pHD['GM_hydro']= GM_hydro
        INFO('Remember to add GM_hydro to mass matrix for models with full assembly')
        # TODO add GM_hydro after assembly
#         print('GM_hydro:\n', GM_hydro)
#         for i, phi in enumerate(WT.fnd.PhiU):
#             # NOTE: diagonal only?
#             WT.MM[i,i]+=GM_hydro[i,i]
#     if tuneM:
#         raise Exception('Removed ?')
#         # TODO remove this, it's application specific
#         # Tuning
#         factM1=1.009
#         factM2=1.003
#         WT2.MM[0,0]*=factM1
#         WT2.MM[1,1]*=factM2
#         WT2.KK[0,0]*=factM1
#         WT2.KK[1,1]*=factM2

        # --- Load hydroShape
        if hydroShape is not None:
            WARN('Not recommended to set hydroShape at this stage')
            self.WT.HD_setShapeFunction(hydroShape)

# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def kinematics(qDict, qdDict, qddDict=None, r_F0=None, r_T0=None, twr=None, fnd=None, s_NGn0=None, 
        tilt=0,
        algo='OpenFAST'):
    """ 
    INPUTS:
     - qDict: dictionary for degrees of freedom with optional keys:
            DOF_f=['Sg','Sw','Hv','R','P','Y']
            DOF_t=['TFA1', 'TFA2', 'TSS1','TSS2']
            DOF_n=['Yaw']
     - r_F0:     undisplaced position of platform ref point, in global coord: (0, 0, PtfmRefzt)
     - r_FT0:    undisplaced position of tower base        , in global coord  (0, 0, TowerBsHt)
     - s_NGn0:   undisplaced position of nacelle COG      , in nacelle coord (NacCMxn, NacCMyn, NacCMzn)
     - twr: BeamBody 
     - fnd: BeamBody 

    """

    # --- Dealing with optional arguments
    if qddDict is None:
        qddDict = dict([(k,0) for k in qDict.keys()])
    if r_F0 is None:
        r_F0=np.array([0,0,0])
    if r_T0 is None:
        r_T0=np.array([0,0,0])
    if algo=='OpenFAST':
        rot_type = 'smallRot_OF'
    else:
        raise NotImplementedError()

    # --- DOFs
    DOF_f = ['Sg','Sw','Hv','R','P','Y']
    q_f   = np.array([qDict  [DOF] if DOF in qDict.keys()   else 0 for DOF in DOF_f])
    qd_f  = np.array([qdDict [DOF] if DOF in qdDict.keys()  else 0 for DOF in DOF_f])
    qdd_f = np.array([qddDict[DOF] if DOF in qddDict.keys() else 0 for DOF in DOF_f])
    DOF_t = np.array(['TFA1', 'TFA2', 'TSS1', 'TSS2'])[twr.shapes]
    q_t   = np.array([qDict[DOF] for DOF in DOF_t])
    qd_t  = np.array([qdDict[DOF] for DOF in DOF_t])
    qdd_t = np.array([qddDict[DOF] for DOF in DOF_t])
    qYaw   = qDict['Yaw']
    qdYaw  = qdDict['Yaw']
    qddYaw = qddDict['Yaw']

    d = dict() # Outputs

    # --- Monopile flexible motion
    if fnd is not None:
        fnd.updateFlexibleKinematics(q_f, qd_f, qdd_f)




    # --- Ref point/fnd motion
    r_F      = r_F0 + q_f[:3]
    v_F      = qd_f[:3]
    a_F      = qdd_f[:3]
    theta_f  = q_f  [3:]
    omega_f  = qd_f [3:]
    omegad_f = qdd_f[3:]
    R_f2g    = rotMat(q_f[3:], rot=rot_type)
    R_g2f    = R_f2g.T

    # Store in dict
    d['r_F'] = r_F
    d['v_F'] = v_F
    d['a_F'] = a_F
    d['R_g2f'] = R_g2f
    d['theta_f'] = theta_f
    d['omega_f'] = omega_f
    d['omegad_f'] = omegad_f

    # --- Tower base motion
    R_t2g      = R_f2g.copy()
    R_g2t      = R_t2g.T
    s_FT0_in_f = r_T0-r_F0
    r_FT       = R_f2g.dot(s_FT0_in_f)
    r_T, v_T, a_T = rigidBodyMotion2Points(r_F, v_F, a_F, omega_f, omegad_f, r_FT) 
    theta_t = theta_f.copy()
    omega_t = omega_f.copy()
    omegad_t = omegad_f.copy()

    d['R_g2t'] = R_g2t
    d['r_T'] = r_T
    d['v_T'] = v_T
    d['a_T'] = a_T
    d['theta_ft'] = np.array((0,0,0))
    d['theta_t']  = theta_t
    d['omega_t']  = omega_t
    d['omegad_t']  = omegad_t

    # --- Tower section motions
    nTwrSpan = len(twr.s_span)
    u_Ts_in_t = np.zeros((nTwrSpan,3))
    udd_Ts_in_t = np.zeros((nTwrSpan,3))
    theta_TTs_in_t = np.zeros((nTwrSpan,3))
    r_Ts = np.zeros((nTwrSpan,3))
    v_Ts = np.zeros((nTwrSpan,3))
    a_Ts = np.zeros((nTwrSpan,3))
    R_g2Ts = np.zeros((nTwrSpan,3,3)) 
    theta_TTs = np.zeros((nTwrSpan,3)) 
    theta_Ts  = np.zeros((nTwrSpan,3))
    omega_Ts  = np.zeros((nTwrSpan,3))
    omegad_Ts = np.zeros((nTwrSpan,3))
    twr.updateFlexibleKinematics(q_t, qd_t, qdd_t) # yams.bodies.py <<<<<<<<<<<<<<<<<<<<<<<<<
    for j in range(nTwrSpan):
        # TODO TODO TODO
        # Missing dipsplacement, velocity, and acceleration due to shoterning of beam
        # TODO TODO TODO Accelerations need debugging
        s_TTs0_in_t = twr.s_G0[:,j]  # undisplaced position
        u_Ts_in_t[j,:]   = twr.U[:,j]     # displacement field
        ud_Ts_in_t       = twr.UP[:,j]    # elastic velocity
        udd_Ts_in_t[j,:] = twr.UPP[:,j]    # elastic acceleration

        if twr.main_axis=='z':
            theta_TTs_in_t[j,:]  = np.array([-twr.V[1,j]  , twr.V[0,j] , 0])
            omega_TTs_in_t  = np.array([-twr.VP[1,j] , twr.VP[0,j], 0])
            omegad_TTs_in_t = np.array([-twr.VPP[1,j] , twr.VPP[0,j], 0])
        else:
            raise NotImplementedError()

        theta_TTs[j,:] =  R_t2g.dot(theta_TTs_in_t[j,:] )
        theta_Ts[j,:] =  theta_t + theta_TTs[j,:] # OK because small angle

        R_Ts2t = rotMat(theta_TTs_in_t[j,:], rot=rot_type)
        R_Ts2g = R_t2g.dot(R_Ts2t)
        R_g2Ts[j,:,:] = R_Ts2g.T

        omega_TTs = R_t2g.dot(omega_TTs_in_t)
        omegad_TTs = R_t2g.dot(omegad_TTs_in_t) 
        omega_Ts[j,:] = omega_t + omega_TTs
        omegad_Ts[j,:] = omegad_t + omegad_TTs + np.cross(omega_t, omega_TTs) # TODO double check extra contrib

        s_TTs_in_t  = s_TTs0_in_t + u_Ts_in_t[j,:] # displaced position
        r_TTs = R_t2g.dot(s_TTs_in_t)
        ud_Ts = R_t2g.dot(ud_Ts_in_t)
        udd_Ts = R_t2g.dot(udd_Ts_in_t[j,:])
        r_Ts[j,:] = r_T + r_TTs
        v_Ts[j,:] = v_T + np.cross(omega_t, r_TTs) + ud_Ts
        a_Ts[j,:] = a_T + np.cross(omega_t, np.cross(omega_t, r_TTs)) + np.cross(omegad_t, r_TTs) 
        a_Ts[j,:] += 2* np.cross(omega_t, ud_Ts) +  udd_Ts

    d['theta_TTs_in_t']  = theta_TTs_in_t
    d['u_Ts_in_t']  = u_Ts_in_t
#     d['udd_Ts_in_t']  = udd_Ts_in_t
    d['r_Ts']   = r_Ts
    d['v_Ts']   = v_Ts
    d['a_Ts']   = a_Ts
    d['R_g2Ts'] = R_g2Ts
    d['theta_Ts']  = theta_Ts
    d['theta_fTs'] = theta_TTs
    d['omega_Ts']  = omega_Ts
    d['omegad_Ts'] = omegad_Ts

    # --- Tower Top point (before Yaw)
    s_TTT0_in_t = twr.s_G0[:,-1] # undisplaced position
    r_TT0 =  r_T0 +  s_TTT0_in_t # undisplaced position of tower top 
    r_TT_undisp =  r_T +  R_t2g.dot(s_TTT0_in_t) # undisplaced, but rotated position of tower top 
    r_TT = r_Ts[-1,:]
    v_TT = v_Ts[-1,:]
    a_TT = a_Ts[-1,:]
    R_g2tt = R_g2Ts[-1,:,:] # To Tower Top
    omega_tt = omega_Ts[-1,:]
    omegad_tt = omegad_Ts[-1,:]
    d['R_g2p'] = R_g2tt

    # --- Nacelle Point/Body (last of tower)
    R_tt2n = R_z(-qYaw)
    R_g2n  = R_tt2n.dot(R_g2tt)
    R_n2g  = R_g2n.T
    omega_tt2n_in_t  = np.array((0,0,qdYaw))
    omegad_tt2n_in_t = np.array((0,0,qddYaw))
    omega_tt2n      = R_n2g.dot(omega_tt2n_in_t)
    omegad_tt2n     = R_n2g.dot(omegad_tt2n_in_t)
    r_N = r_TT
    v_N = v_TT
    a_N = a_TT
    omega_n = omega_tt   + omega_tt2n
    omegad_n = omegad_tt + omegad_tt2n + np.cross(omega_tt, omega_tt2n)
    d['u_N_tot'] = r_N-r_TT0
    d['u_N']    = r_N-r_TT_undisp
    d['r_N']    = r_N
    d['v_N']    = v_N
    d['a_N']    = a_N
    d['R_tt2n'] = R_tt2n
    d['R_g2n']  = R_g2n
    d['omega_n']  = omega_n
    d['omegad_n'] = omegad_n

    # --- Nacelle COG
    dGn = rigidBodyKinematics(s_NGn0, r_N, R_g2n, v_N, omega_n, a_N=a_N, omegad_n=omegad_n, point_name='Gn', source_name='N')
    d.update(dGn)

    # --- Shaft
    R_s2n = R_y(tilt)  # Rotation fromShaft to Nacelle
    R_g2s = (R_s2n.T).dot(R_g2n)
    d['R_g2s'] = R_g2s


    return d

def rigidBodyKinematics(s_NP0_in_n, r_N, R_g2n, v_N, omega_n, a_N=None, omegad_n=None, point_name='P', source_name='O'):
    """ 
    Simple rigid body motion kinematics for a point in the nacelle
    """
    if a_N is None:
        a_N = np.array([0,0,0])
    if omegad_n is None:
        omegad_n = np.array([0,0,0])

    r_NP = (R_g2n.T).dot(s_NP0_in_n)
    r_P, v_P, a_P = rigidBodyMotion2Points(r_N, v_N, a_N, omega_n, omegad_n, r_NP) 
    s = point_name
    d = {'r_'+source_name+s:r_NP, 'r_'+s:r_P, 'v_'+s:v_P, 'a_'+s:a_P}
    return d





if __name__ == '__main__':
    np.set_printoptions(linewidth=300, precision=2)
    WT = FASTWindTurbine('../../data/NREL5MW/Main_Onshore.fst')
    print(WT)

