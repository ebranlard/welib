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
import welib.weio as weio
from collections import OrderedDict
from welib.yams.bodies import RigidBody, FlexibleBody, FASTBeamBody
from welib.yams.rotations import R_x, R_y, R_z, rotMat
from welib.yams.kinematics import rigidBodyMotion2Points

class WindTurbineStructure():
    def __init__(self):
        self.bld = None # list of blades
        self.hub = None
        self.gen = None
        self.rot = None # hub+blds
        self.nac = None
        self.twr = None
        self.fnd = None
        # Geometry
        self.r_NS_inN = None
        self.r_NR_inN = None
        self.r_SR_inS = None
        self.shaft_tilt = None # [rad]

        # Information relevant for simulation
        self.DOF    ={'name':'', 'active':None, 'q0':None, 'qd0':None,'q_channel':None, 'qd_channel':None, 'qdd_channel':None}

        self.algo=None # 'OpenFAST'
        # Derived properties
        #self.DOFname
        #self.q0     
        #self.qd0    


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
        return FASTWindTurbine(fstFilename)

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
        try:
            return self.FST['Gravity'] # NEW
        except:
            return self.ED['Gravity'] # OLD OpenFAST


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
        try:
            p['g']        = WT.FST['Gravity'] # NEW
        except:
            p['g']        = WT.ED['Gravity']

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
            try:
                p['v_xT2c']   = WT.twr.Bhat_t_bc[0,1]  # Mode 2
            except:
                print('[FAIL] YAMS: WindTurbine, TODO which mode shape is which for tower top rotations factors')
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

    def checkPackage(self, pkg):
        # --- Checking package info
        info = pkg.info()
        nDOFExpected= info['nq']
        print('DOFs:', self.DOFname, 'Model:',info['name'], 'nDOF:',nDOFExpected )
        if len(self.DOFname)!=nDOFExpected:
            raise Exception('Inconsistency in number of DOFs')

    def simulate_py(self, pkg, p, time, u=None, acc=False, forcing=False):
        """ Perform non-linear simulation based on a `model` (python package) generated by yams_sympy """
        from welib.yams.models.packman import simulate

        # --- Checking package info
        self.checkPackage(pkg)

        print('------------------ WINDTURBINE: NON LINEAR SIMULATION -----------------------')
        q0  = np.asarray(self.q0)
        qd0 = np.asarray(self.qd0)
        resNL, sysNL, dfNL = simulate(pkg, time, q0, qd0=qd0, p=p, u=u, acc=acc, forcing=forcing, 
                sStates=self.channels, Factors=self.FASTDOFScales, sAcc=self.acc_channels)
        print('-----------------------------------------------------------------------------')

        return resNL, sysNL, dfNL


    def py_lin(self, pkg, p, time, uop=None, qop=None, qdop=None, du=None, MCKextra=None, MCKu=None, noBlin=False):
        """ Perform linear simulation based on a model (python package) generated by yams_sympy """
        # TODO TODO TODO MOVE ME TO packman
        from welib.yams.models.packman import linearModel
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


    def simulate_py_lin(self, pkg, p, time, uop=None, qop=None, qdop=None, du=None, MCKextra=None, MCKu=None, acc=False, forcing=False, noBlin=False):
        """ Perform linear simulation based on a model (python package) generated by yams_sympy """
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
        #from welib.beams.theory import UniformBeamGuyanModes, UniformBeamBendingModes
        #from welib.yams.yams import UniformBeamBody
        from welib.fast.postpro import ED_TwrGag #, ED_TwrStations, getEDClass
        #from welib.fast.elastodyn import rotorParameters, towerParameters, ED_Parameters, ED_CoordSys, ED_Positions
        #from welib.fast.elastodyn import ED_AngPosVelPAcc, ED_LinVelPAcc, ED_qDict2q
        #from welib.weio.fast_input_deck import FASTInputDeck
        #from welib.yams.windturbine import FASTWindTurbine
        from welib.yams.flexibility import beamSectionLoads3D  # calls beamSectionLoads1D
        #from welib.yams.rotations import rotMat, R_y, R_z, R_x
        #from welib.yams.kinematics import *

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
            F_sec, M_sec =  beamSectionLoads3D(p_ext=p_ext, F_top=F_top_t, M_top=M_top_t, s_span=twr.s_span, m=twr.m, U=twr.U, V=twr.V, K=twr.K, a_struct=a_struct_t, 
                     a_ext=a_ext, corrections=1)
            return F_sec, M_sec


        # Sanitization of input dataframe
        df = df.loc[:,~df.columns.duplicated()].copy()
        df.columns = [  v.split('_[')[0] for v in df.columns.values]
        df.reset_index(inplace=True)

        # --- States
        DOFNames = ['Sg', 'Sw', 'Hv' ,'R', 'P', 'Y', 'TFA1', 'TSS1', 'Yaw']
        sq   = ['Q_'+s for s in DOFNames]
        sqd  = ['QD_'+s for s in DOFNames]
        sqdd = ['QD2_'+s for s in DOFNames]
        sqall = sq+sqd+sqdd
        for s in sqall:
            if s not in df.keys():
                print('[WARN] Missing DOF from dataframe: {}'.format(s))
                df[s]=0
        Q   = df[sq]
        QD  = df[sqd]
        QDD = df[sqdd]
        Q.columns   = DOFNames
        QD.columns  = DOFNames
        QDD.columns = DOFNames
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
        for iiSL,hED in enumerate(HEDOut):
            sT='TwHt{}'.format(iiSL+1)
            colOut+=[sT+'TPxi_[m]'  , sT+'TPyi_[m]'  , sT+'TPzi_[m]']
            colOut+=[sT+'TDxt_[m]'  , sT+'TDyt_[m]'  , sT+'TDzt_[m]']
            colOut+=[sT+'RDxt_[deg]', sT+'RDyt_[deg]', sT+'RDzt_[deg]']
            colOut+=[sT+'ALxt_[m/s^2]', sT+'ALyt_[m/s^2]', sT+'ALzt_[m/s^2]']
            colOut+=[sT+'FLxt_[kN]', sT+'FLyt_[kN]', sT+'FLzt_[kN]']
            colOut+=[sT+'MLxt_[kN-m]', sT+'MLyt_[kN-m]', sT+'MLzt_[kN-m]']
        # TODO TODO TODO FIGURE OUT WHY THIS RETURN DTYPE OBJECT
        dfOut = pd.DataFrame(index=df.index, columns=colOut)

        # --- Calc Output per time step
        with Timer('Kinematics'):
            for it,t in enumerate(df['Time']):
                q   = Q.iloc[it,:].copy()
                qd  = QD.iloc[it,:].copy()
                qdd = QDD.iloc[it,:].copy()
                # TODO TODO Sort out issue of convention in OpenFAST
                q['TSS1'] = -q['TSS1']
                qd['TSS1'] = -qd['TSS1']
                qdd['TSS1'] = -qdd['TSS1']
    #             q['TSS2'] = -q['TSS2']

                dd = WT.kinematics(q, qd, qdd)
                dfOut['Time_[s]'].loc[it] = t
                # TDi includes all platform motions
                dfOut['TwrTpTDxi'].loc[it] = dd['u_N_tot'][0] 
                dfOut['TwrTpTDyi'].loc[it] = dd['u_N_tot'][1]
                dfOut['TwrTpTDzi'].loc[it] = dd['u_N_tot'][2]

                u_N = dd['u_N']
                v_N = dd['v_N']
                a_N = dd['a_N']
                om_N = dd['omega_n']
                omd_N = dd['omegad_n']
                u_N_p = dd['R_g2p'].dot(u_N)
                u_N_t = dd['R_g2t'].dot(u_N)
                v_N_p = dd['R_g2p'].dot(v_N)
                a_N_p = dd['R_g2p'].dot(a_N)
                om_N_p = dd['R_g2p'].dot(om_N)
                omd_N_p = dd['R_g2p'].dot(omd_N)

                dfOut['YawBrTDxt'].loc[it] = u_N_t[0]
                dfOut['YawBrTDyt'].loc[it] = u_N_t[1]
                dfOut['YawBrTDzt'].loc[it] = u_N_t[2]
                dfOut['YawBrTDxp'].loc[it] = u_N_p[0]
                dfOut['YawBrTDyp'].loc[it] = u_N_p[1]
                dfOut['YawBrTDzp'].loc[it] = u_N_p[2]

                dfOut['YawBrTVxp'].loc[it] = v_N_p[0]
                dfOut['YawBrTVyp'].loc[it] = v_N_p[1]
                dfOut['YawBrTVzp'].loc[it] = v_N_p[2]
                dfOut['YawBrTAxp'].loc[it] = a_N_p[0]
                dfOut['YawBrTAyp'].loc[it] = a_N_p[1]
                dfOut['YawBrTAzp'].loc[it] = a_N_p[2]
                dfOut['YawBrRVxp'].loc[it] = om_N_p[0] * 180/np.pi
                dfOut['YawBrRVyp'].loc[it] = om_N_p[1] * 180/np.pi
                dfOut['YawBrRVzp'].loc[it] = om_N_p[2] * 180/np.pi
                dfOut['YawBrRAxp'].loc[it] = omd_N_p[0] * 180/np.pi
                dfOut['YawBrRAyp'].loc[it] = omd_N_p[1] * 180/np.pi
                dfOut['YawBrRAzp'].loc[it] = omd_N_p[2] * 180/np.pi

                a_IMU = dd['a_IMU']
                v_IMU = dd['v_IMU']
                om_IMU = dd['omega_n']
                omd_IMU = dd['omegad_n']

                a_IMU_s = dd['R_g2s'].dot(a_IMU)
                v_IMU_s = dd['R_g2s'].dot(v_IMU)
                om_IMU_s = dd['R_g2s'].dot(om_IMU)
                omd_IMU_s = dd['R_g2s'].dot(omd_IMU)

                dfOut['NcIMUTVxs'].loc[it] = v_IMU_s[0]
                dfOut['NcIMUTVys'].loc[it] = v_IMU_s[1]
                dfOut['NcIMUTVzs'].loc[it] = v_IMU_s[2]
                dfOut['NcIMUTAxs'].loc[it] = a_IMU_s[0]
                dfOut['NcIMUTAys'].loc[it] = a_IMU_s[1]
                dfOut['NcIMUTAzs'].loc[it] = a_IMU_s[2]
                dfOut['NcIMURVxs'].loc[it] = om_IMU_s[0] * 180/np.pi
                dfOut['NcIMURVys'].loc[it] = om_IMU_s[1] * 180/np.pi
                dfOut['NcIMURVzs'].loc[it] = om_IMU_s[2] * 180/np.pi
                dfOut['NcIMURAxs'].loc[it] = omd_IMU_s[0] * 180/np.pi
                dfOut['NcIMURAys'].loc[it] = omd_IMU_s[1] * 180/np.pi
                dfOut['NcIMURAzs'].loc[it] = omd_IMU_s[2] * 180/np.pi
        

                # --- Loads
                gravity_vec = np.array([0,0,-WT.gravity])
                # --- RNA (without Yaw Br) loads
                omd_n = dd['omegad_n']
                om_n = dd['omega_n']
                R_g2p = dd['R_g2p']
                R_g2n = dd['R_g2n']
                r_Grna = dd['r_Grna']
                a_Grna = dd['a_Grna']
                Mrna  = WT.RNA_noYawBr.mass
                JGrna = WT.RNA_noYawBr.masscenter_inertia
                JGrna_g = (R_g2n.T).dot(JGrna).dot(R_g2n)
                F_Grna_grav =  Mrna *gravity_vec
                r_NGrna = dd['r_NGrna']

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

                dfOut['YawBrFxp'].loc[it] = F_N_p[0]/1000
                dfOut['YawBrFyp'].loc[it] = F_N_p[1]/1000
                dfOut['YawBrFzp'].loc[it] = F_N_p[2]/1000
                dfOut['YawBrMxp'].loc[it] = M_N_p[0]/1000
                dfOut['YawBrMyp'].loc[it] = M_N_p[1]/1000
                dfOut['YawBrMzp'].loc[it] = M_N_p[2]/1000


                if useTopLoadsFromDF:
                    F_N_p = np.array((df['YawBrFxp'].loc[it], df['YawBrFyp'].loc[it], df['YawBrFzp'].loc[it]))*1000
                    M_N_p = np.array((df['YawBrMxp'].loc[it], df['YawBrMyp'].loc[it], df['YawBrMzp'].loc[it]))*1000
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

                # --- Section Loads
                F_sec, M_sec = towerSectionLoads(WT.twr, F_N_t, M_N_t, kin=dd, gravity=WT.gravity)

                for iiSL, hSL in enumerate(HEDOut):
                    iSL = np.argmin(np.abs(hSL-WT.twr.s_span))
                    hSL = WT.twr.s_span[iSL]

                    sT='TwHt{}'.format(iiSL+1)
                    dfOut[sT+'FLxt_[kN]'  ].loc[it] = F_sec[0, iSL]/1000
                    dfOut[sT+'FLyt_[kN]'  ].loc[it] = F_sec[1, iSL]/1000
                    dfOut[sT+'FLzt_[kN]'  ].loc[it] = F_sec[2, iSL]/1000
                    dfOut[sT+'MLxt_[kN-m]'].loc[it] = M_sec[0, iSL]/1000
                    dfOut[sT+'MLyt_[kN-m]'].loc[it] = M_sec[1, iSL]/1000
                    dfOut[sT+'MLzt_[kN-m]'].loc[it] = M_sec[2, iSL]/1000

                    dfOut[sT+'TDxt_[m]'].loc[it] = dd['u_Ts_in_t'][iSL,0]
                    dfOut[sT+'TDyt_[m]'].loc[it] = dd['u_Ts_in_t'][iSL,1]
                    dfOut[sT+'TDzt_[m]'].loc[it] = dd['u_Ts_in_t'][iSL,2]
                    dfOut[sT+'RDxt_[deg]'].loc[it] = dd['theta_TTs_in_t'][iSL,0]*180/np.pi
                    dfOut[sT+'RDyt_[deg]'].loc[it] = dd['theta_TTs_in_t'][iSL,1]*180/np.pi
                    dfOut[sT+'RDzt_[deg]'].loc[it] = dd['theta_TTs_in_t'][iSL,2]*180/np.pi
                    a_Ts = R_g2t.dot(dd['a_Ts'][iSL])
                    dfOut[sT+'ALxt_[m/s^2]'].loc[it] = a_Ts[0]
                    dfOut[sT+'ALyt_[m/s^2]'].loc[it] = a_Ts[1]
                    dfOut[sT+'ALzt_[m/s^2]'].loc[it] = a_Ts[2]

                    dfOut[sT+'TPxi_[m]'].loc[it] = dd['r_Ts'][iSL,0]
                    dfOut[sT+'TPyi_[m]'].loc[it] = dd['r_Ts'][iSL,1]
                    dfOut[sT+'TPzi_[m]'].loc[it] = dd['r_Ts'][iSL,2]

        return dfOut



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
def FASTWindTurbine(fstFilename, main_axis='z', nSpanTwr=None, twrShapes=None, nSpanBld=None, algo='',
        bldStartAtRotorCenter=True):
    """
    INPUTS:
     - twrShapes: Select shapes to use for tower. If None, twrShapes=[0,1,2,3]

    """
    # TODO TODO TODO  Harmonize with TNSB.py
    # TODO TODO TODO  Harmonize with fast.elastodyn when algo is OpenFAST
    # --- Reading main OpenFAST files
    ext     = os.path.splitext(fstFilename)[1]
    FST     = weio.read(fstFilename)
    rootdir = os.path.dirname(fstFilename)
    EDfile  = os.path.join(rootdir,FST['EDFile'].strip('"')).replace('\\','/')
    ED      = weio.read(EDfile)
    rootdir = os.path.dirname(EDfile)
    try:
        bldfile = os.path.join(rootdir,ED['BldFile(1)'].strip('"')).replace('\\','/')
    except:
        bldfile = os.path.join(rootdir,ED['BldFile1'].strip('"')).replace('\\','/')
    twrfile = os.path.join(rootdir,ED['TwrFile'].strip('"')).replace('\\','/')
    # TODO SubDyn, MoorDyn, BeamDyn 
    try:
        gravity = FST['gravity']
    except:
        try:
            gravity = ED['gravity']
        except:
            raise Exception('Variable gravity not found in FST file or ED file.')

    r_EPtfm_inE = np.array([0,0,ED['PtfmRefzt']               ])  # TODO TODO TODO
    r_ET_inE    = np.array([0,0,ED['TowerBsHt']               ])  # TODO TODO TODO
    r_TN_inT    = np.array([0,0,ED['TowerHt']-ED['TowerBsHt'] ])
    # Basic geometries for nacelle
    theta_tilt_y = -ED['ShftTilt']*np.pi/180  # NOTE: tilt has wrong orientation in FAST
    R_NS = R_y(theta_tilt_y)  # Rotation fromShaft to Nacelle
    r_NS_inN    = np.array([0             , 0, ED['Twr2Shft']]) # Shaft start in N
    r_SR_inS    = np.array([ED['OverHang'], 0, 0             ]) # Rotor center in S
    r_SGhub_inS = np.array([ED['HubCM']   , 0, 0             ]) + r_SR_inS # Hub G in S
    r_NR_inN    = r_NS_inN + R_NS.dot(r_SR_inS)                 # Rotor center in N
    r_NGnac_inN = np.array([ED['NacCMxn'],0,ED['NacCMzn']    ]) # Nacelle G in N
    r_RGhub_inS = - r_SR_inS + r_SGhub_inS
    if main_axis=='x':
        raise NotImplementedError()

    # --- OpenFAST compatibility
    if algo.lower()=='openfast':
        #print('[INFO] YAMS Wind Turbine - Using Algo OpenFAST')
        from welib.fast.elastodyn import rotorParameters, bladeParameters, towerParameters, bladeDerivedParameters, towerDerivedParameters
        pBld = bladeParameters(EDfile)
        pBld = bladeDerivedParameters(pBld, inertiaAtBladeRoot=False)
        pRot, pbld, phub = rotorParameters(EDfile, identicalBlades=False)
        RotMass = pRot['RotMass']
        pTwr = towerParameters(EDfile, RotMass=RotMass, gravity=gravity)
        pTwr = towerDerivedParameters(pTwr)

    # --- Hub  (defined using point N and nacelle coord as ref)
    M_hub  = ED['HubMass']
    JxxHub_atR = ED['HubIner']
    hub = RigidBody('Hub', M_hub, (JxxHub_atR,0,0), s_OG=r_SGhub_inS, R_b2g=R_NS, s_OP=r_SR_inS, r_O=r_NS_inN) 

    # --- Generator (Low speed shaft) (defined using point N and nacelle coord as ref)
    gen = RigidBody('Gen', 0, (ED['GenIner']*ED['GBRatio']**2,0,0), s_OG=[0,0,0], R_b2g=R_NS,  r_O=r_NS_inN) 

    # --- Nacelle (defined using point N and nacelle coord as ref)
    M_nac = ED['NacMass']
    JyyNac_atN = ED['NacYIner'] # Inertia of nacelle at N in N
    nac = RigidBody('Nac', M_nac, (0,JyyNac_atN,0), r_NGnac_inN, s_OP=[0,0,0])

    # --- Blades 
    bldFile = weio.read(bldfile)
    m    = bldFile['BldProp'][:,3]
    jxxG=0*m
    #if algo=='OpenFAST':
    #    jxxG=0*m
    #else:
    #    jxxG = 0*m     # NOTE: unknown
    #    print('>>> windturbine.py: TODO: using unknown jxxG')
    nB = ED['NumBl']
    bld=np.zeros(nB,dtype=object)
    bld[0] = FASTBeamBody(ED, bldFile, Mtop=0, main_axis=main_axis, jxxG=jxxG, spanFrom0=False, bldStartAtRotorCenter=bldStartAtRotorCenter, nSpan=nSpanBld, gravity=gravity, algo=algo) 
    if algo.lower()=='openfast':
        # Overwrite blade props until full compatibility implemented
        bld[0].MM[0,0]    = pBld['BldMass']
        bld[0].MM[1,1]    = pBld['BldMass']
        bld[0].MM[2,2]    = pBld['BldMass']
        # TODO TODO bldShapes
        bld[0].MM [6:,6:] = pBld['Me']
        bld[0].KK0[6:,6:] = pBld['Ke0'] # NOTE: Ke has no stiffening
        bld[0].DD [6:,6:] = pBld['De']

    for iB in range(nB-1):
        bld[iB+1]=copy.deepcopy(bld[0])
        #bld[iB+1].R_b2g
    for iB,B in enumerate(bld):
        B.name='bld'+str(iB+1)
        psi_B= -iB*2*np.pi/len(bld) 
        if main_axis=='x':
            R_SB = R_z(0*np.pi + psi_B) # TODO psi offset and psi0
        elif main_axis=='z':
            R_SB = R_x(0*np.pi + psi_B) # TODO psi0
        if bldStartAtRotorCenter:
            R_SB = np.dot(R_SB, R_y(ED['PreCone({})'.format(iB+1)]*np.pi/180)) # blade2shaft
            B.R_b2g= R_SB
        else:
            print('>>>> TODO TODO TODO Wind Turbine R_SB')
            R_SB = np.dot(R_SB, R_y(ED['PreCone({})'.format(iB+1)]*np.pi/180)) # blade2shaft
            B.R_b2g= R_SB

    # --- Blades (with origin R, using N as "global" ref)
    blds_rigid = rigidBlades(bld, r_O = [0,0,0])
    blds_rigid.pos_global = r_NR_inN
    blds_rigid.R_b2g      = R_NS

    # --- Rotor = Hub + Blades (with origin R, using N as global ref)
    rot = blds_rigid.combine(hub, R_b2g=R_NS, r_O=blds_rigid.pos_global)
    rot.name='rotor'
    rotgen = rot.combine(gen, R_b2g=R_NS, r_O=blds_rigid.pos_global)
    #print(rotgen)

    #--- Yaw bearing, at tower top
    M_yawBr = ED['YawBrMass']
    yawBr = RigidBody('YawBr', M_yawBr, J=(0,0,0), s_OG=(0,0,0))


    # --- RNA 
    RNA = rot.combine(gen).combine(nac,r_O=[0,0,0]).combine(yawBr, r_O=[0,0,0])
    RNA.name='RNA'
    #print(RNA)
    M_RNA = RNA.mass

    # --- RNA without YawBr Mass
    RNA_noYawBr = rot.combine(gen).combine(nac,r_O=[0,0,0])

    # --- Fnd (defined wrt ground/MSL "E")
#     print(FST.keys())
    M_fnd = ED['PtfmMass']
    r_EGfnd_inF = np.array([ED['PtfmCMxt'],ED['PtfmCMyt'],ED['PtfmCMzt']])
    r_EPtfm_inF    = np.array([0             ,0             ,ED['PtfmRefzt']]) # TODO, this is wrong
    r_PtfmGfnd_inF = -r_EPtfm_inF + r_EGfnd_inF
    fnd = RigidBody('fnd', M_fnd, (ED['PtfmRIner'], ED['PtfmPIner'], ED['PtfmYIner']), s_OG=r_PtfmGfnd_inF, r_O=r_EPtfm_inF) 

    # --- Twr
    if twrShapes is None: 
        twrShapes=[]
        if ED['TwFADOF1']:
            twrShapes+=[0]
        if ED['TwFADOF2']:
            twrShapes+=[1]
        if ED['TwSSDOF1']:
            twrShapes+=[2]
        if ED['TwSSDOF2']:
            twrShapes+=[3]
    twrFile = weio.read(twrfile)
    twr = FASTBeamBody(ED, twrFile, Mtop=M_RNA, main_axis='z', bAxialCorr=False, bStiffening=True, shapes=twrShapes, nSpan=nSpanTwr, algo=algo, gravity=gravity) # TODO options
    twr_rigid  = twr.toRigidBody()
    twr_rigid.pos_global = r_ET_inE
    if algo=='OpenFAST':
        twr.MM[0,0]   = pTwr['TwrMass']
        twr.MM[1,1]   = pTwr['TwrMass']
        twr.MM[2,2]   = pTwr['TwrMass']
        twr.MM      [6:,6:] = pTwr['Me']   [np.ix_(twrShapes,twrShapes)]
        twr.KK      [6:,6:] = pTwr['Ke']   [np.ix_(twrShapes,twrShapes)]
        twr.KK0     [6:,6:] = pTwr['Ke0']  [np.ix_(twrShapes,twrShapes)]
        twr.KKg_self[6:,6:] = pTwr['Kg_SW'][np.ix_(twrShapes,twrShapes)]
        twr.KKg_Mtop[6:,6:] = pTwr['Kg_TM'][np.ix_(twrShapes,twrShapes)]
        twr.DD      [6:,6:] = pTwr['De']   [np.ix_(twrShapes,twrShapes)]


    # --- Full WT rigid body equivalent, with point T as ref
    RNAb = copy.deepcopy(RNA)
    RNAb.pos_global = r_TN_inT+r_ET_inE
    RNAb.R_b2g      = np.eye(3)
    WT_rigid = RNAb.combine(twr_rigid, r_O=r_ET_inE).combine(fnd, r_O=r_ET_inE) # TODO TODO TODO T or Ptfm
    #WT_rigid = RNAb.combine(twr_rigid, r_O=r_EPtfm_inE).combine(fnd, r_O=r_EPtfm_inE) # TODO T or Ptfm

    # --- Moorings
    MAP = None
    K_Moor=np.zeros((6,6))
    if FST['CompMooring']==1:
        from welib.moor.mappp import Map
#         try:
        MAP = Map(fstFilename)
        #zRef  = twr.pos_global[2]  
        zRef  = fnd.pos_global[2]  
        K_Moor,_ = MAP.stiffness_matrix(epsilon=1e-2, point=(0,0,zRef))
        #K_Moor2,_ = MAP.stiffness_matrix(epsilon=1e-2, point=(0,0,0))
#         except:
#             print('YAMS Wind Turbine: problem loading MAP model (only supported on windows)')

    elif FST['CompMooring']==2:
        print('YAMS Wind Turbine: TODO MoorDyn')


    # --- Degrees of freedom
    DOFs=[]
    DOFs+=[{'name':'x'      , 'active':ED['PtfmSgDOF'], 'q0': ED['PtfmSurge']  , 'qd0':0 , 'q_channel':'PtfmSurge_[m]' , 'qd_channel':'QD_Sg_[m/s]','qdd_channel':'QD2_Sg_[m/s^2]'}]
    DOFs+=[{'name':'y'      , 'active':ED['PtfmSwDOF'], 'q0': ED['PtfmSway']   , 'qd0':0 , 'q_channel':'PtfmSway_[m]'  , 'qd_channel':'QD_Sw_[m/s]','qdd_channel':'QD2_Sw_[m/s^2]'}]
    DOFs+=[{'name':'z'      , 'active':ED['PtfmHvDOF'], 'q0': ED['PtfmHeave']  , 'qd0':0 , 'q_channel':'PtfmHeave_[m]' , 'qd_channel':'QD_Hv_[m/s]','qdd_channel':'QD2_Hv_[m/s^2]'}]

    DOFs+=[{'name':'phi_x' , 'active':ED['PtfmRDOF'] , 'q0': ED['PtfmRoll']*np.pi/180  , 'qd0':0 , 'q_channel':'PtfmRoll_[deg]'  , 'qd_channel':'QD_R_[rad/s]', 'qdd_channel':'QD2_R_[rad/s^2]'}]
    DOFs+=[{'name':'phi_y' , 'active':ED['PtfmPDOF'] , 'q0': ED['PtfmPitch']*np.pi/180 , 'qd0':0 , 'q_channel':'PtfmPitch_[deg]' , 'qd_channel':'QD_P_[rad/s]', 'qdd_channel':'QD2_P_[rad/s^2]'}]
    DOFs+=[{'name':'phi_z' , 'active':ED['PtfmYDOF'] , 'q0': ED['PtfmYaw']*np.pi/180   , 'qd0':0 , 'q_channel':'PtfmYaw_[deg]'   , 'qd_channel':'QD_Y_[rad/s]', 'qdd_channel':'QD2_Y_[rad/s^2]'}]

    DOFs+=[{'name':'q_FA1'  , 'active':ED['TwFADOF1'] , 'q0': ED['TTDspFA']  , 'qd0':0 , 'q_channel':'Q_TFA1_[m]', 'qd_channel':'QD_TFA1_[m/s]', 'qdd_channel':'QD2_TFA1_[m/s^2]'}]
    DOFs+=[{'name':'q_SS1'  , 'active':ED['TwSSDOF1'] , 'q0': ED['TTDspSS']  , 'qd0':0 , 'q_channel':'Q_TSS1_[m]', 'qd_channel':'QD_TSS1_[m/s]', 'qdd_channel':'QD2_TSS1_[m/s^2]'}]
    DOFs+=[{'name':'q_FA2'  , 'active':ED['TwFADOF2'] , 'q0': ED['TTDspFA']  , 'qd0':0 , 'q_channel':'Q_TFA2_[m]', 'qd_channel':'QD_TFA2_[m/s]', 'qdd_channel':'QD2_TFA2_[m/s^2]'}]
    DOFs+=[{'name':'q_SS2'  , 'active':ED['TwSSDOF2'] , 'q0': ED['TTDspSS']  , 'qd0':0 , 'q_channel':'Q_TSS2_[m]', 'qd_channel':'QD_TSS2_[m/s]', 'qdd_channel':'QD2_TSS2_[m/s^2]'}]

    DOFs+=[{'name':'theta_y','active':ED['YawDOF']  , 'q0': ED['NacYaw']*np.pi/180   , 'qd0':0 ,          'q_channel':'NacYaw_[deg]' , 'qd_channel':'QD_Yaw_[rad/s]', 'qdd_channel':'QD2_Yaw_[rad/s^2]'}]
    DOFs+=[{'name':'psi'    ,'active':ED['GenDOF']  , 'q0': ED['Azimuth']*np.pi/180  , 'qd0':ED['RotSpeed']*2*np.pi/60 , 'q_channel':'Azimuth_[deg]', 'qd_channel':'RotSpeed_[rpm]', 'qdd_channel': 'QD2_GeAz_[rad/s^2]'}]

    DOFs+=[{'name':'nu'     ,'active':ED['DrTrDOF'] , 'q0': 0  , 'qd0':0 , 'q_channel':'Q_DrTr_[rad]', 'qd_channel':'QD_DrTr_[rad/s]', 'qdd_channel':'QD2_DrTr_[rad/s^2]'}]

    # 
    for ib in np.arange(ED['NumBl']):
        B=str(ib+1)
        DOFs+=[{'name':'q_B{}Fl1'.format(B), 'active':ED['FlapDOF1'] , 'q0': ED['OOPDefl'], 'qd0':0, 'q_channel':'Q_B{}F1_[m]'.format(B), 'qd_channel':'QD_B{}F1_[m/s]'.format(B), 'qdd_channel':'QD2_B{}F1_[m/s^2]'.format(B)}]
        DOFs+=[{'name':'q_B{}Ed1'.format(B), 'active':ED['FlapDOF2'] , 'q0': ED['OOPDefl'], 'qd0':0, 'q_channel':'Q_B{}F2_[m]'.format(B), 'qd_channel':'QD_B{}E1_[m/s]'.format(B), 'qdd_channel':'QD2_B{}E1_[m/s^2]'.format(B)}]
        DOFs+=[{'name':'q_B{}Ed1'.format(B), 'active':ED['EdgeDOF']  , 'q0': ED['IPDefl'] , 'qd0':0, 'q_channel':'Q_B{}E1_[m]'.format(B), 'qd_channel':'QD_B{}E1_[m/s]'.format(B), 'qdd_channel':'QD2_B{}E1_[m/s^2]'.format(B)}]

# ---------------------- DEGREES OF FREEDOM --------------------------------------
# False          FlapDOF1    - First flapwise blade mode DOF (flag)
# False          FlapDOF2    - Second flapwise blade mode DOF (flag)
# False          EdgeDOF     - First edgewise blade mode DOF (flag)
# ---------------------- INITIAL CONDITIONS --------------------------------------
#           0   OoPDefl     - Initial out-of-plane blade-tip displacement (meters)
#           0   IPDefl      - Initial in-plane blade-tip deflection (meters)


    # --- Return
    WT = WindTurbineStructure()
    WT.bld        = bld        # origin at R
    WT.blds_rigid = blds_rigid # origin at R
    WT.hub        = hub        # origin at S
    WT.rot        = rot        # origin at R, rigid body bld+hub
    WT.rotgen     = rotgen     # origin at R, rigid body bld+hub+genLSS
    WT.gen        = gen        # origin at S
    WT.nac        = nac        # origin at N
    WT.yawBr      = yawBr      # origin at N
    WT.twr        = twr        # origin at T
    WT.fnd        = fnd        # origin at T
    WT.MAP        = MAP        # typically at (0,0,0) NOTE: not a body
    WT.K_Moor     = K_Moor     # HACK..
    WT.RNA        = RNA          # origin at N, rigid body bld+hub+gen+nac+yawBr
    WT.RNA_noYawBr = RNA_noYawBr # origin at N, rigid body bld+hub+gen+nac
    WT.WT_rigid   = WT_rigid   # rigid wind turbine body, origin at MSL

    WT.DOF= DOFs

    #WT.r_ET_inE = 
    #WT.r_TN_inT
    # --- Geometry
    WT.r_NS_inN = r_NS_inN
    WT.r_NR_inN = r_NR_inN
    WT.r_SR_inS = r_SR_inS
    WT.shaft_tilt = -ED['ShftTilt']*np.pi/180  # NOTE: tilt has wrong orientation in FAST
    #R_NS = R_y(WT.shaft_tilt)  # Rotation fromShaft to Nacelle


    WT.ED=ED
    WT.FST=FST
    WT.algo=algo

    return WT



# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def kinematics(qDict, qdDict, qddDict=None, r_F0=None, r_T0=None, twr=None, s_NGn0=None, 
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
    twr.updateFlexibleKinematics(q_t, qd_t, qdd_t) # yams.bodies.py
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

