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
import copy
import welib.weio as weio
from collections import OrderedDict
from welib.yams.utils import R_x, R_y, R_z
from welib.yams.bodies import RigidBody, FlexibleBody, FASTBeamBody

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

        # Information relevant for simulation
        self.DOF    ={'name':'', 'active':None, 'q0':None, 'qd0':None,'q_channel':None, 'qd_channel':None, 'qdd_channel':None}

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
        return np.array([dof['q0'] for dof in self.DOF if dof['active']])
    @property
    def qd0(self):
        return np.array([dof['qd0'] for dof in self.DOF if dof['active']])
    @property
    def DOFname(self):
        return np.array([dof['name'] for dof in self.DOF if dof['active']])
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
            p['v_yT1c']   = WT.twr.Bhat_t_bc[1,0]  # Mode 1  3 x nShapes
            p['v_xT2c']   = WT.twr.Bhat_t_bc[0,1]  # Mode 2
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
                    K_Moor,_ = WT.MAP.stiffness_matrix(epsilon=1e-2, point=(0,0,p['z_OT']))
                else:
                    K_Moor,_ = WT.MAP.stiffness_matrix(epsilon=1e-2, point=(0,0,0))
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

        print('---------------------- NON LINEAR SIMULATION --------------------------------')
        q0  = self.q0
        qd0 = self.qd0
        resNL, sysNL, dfNL = simulate(pkg, time, q0, qd0=qd0, p=p, u=u, acc=acc, forcing=forcing, 
                DOFs=self.channels, Factors=self.FASTDOFScales, sAcc=self.acc_channels)
        print('-----------------------------------------------------------------------------')

        return resNL, sysNL, dfNL


    def py_lin(self, pkg, p, time, uop=None, qop=None, qdop=None, du=None, MCKextra=None, MCKu=None, noBlin=False):
        """ Perform linear simulation based on a model (python package) generated by yams_sympy """
        # TODO TODO TODO MOVE ME TO packman
        from welib.yams.models.packman import linearModel
        print('-------------------------- LINEAR SIMULATION --------------------------------')
        # --- Checking package info
        self.checkPackage(pkg)
        if qop is None:
            qop = self.q0*0 
        if qdop is None:
            qdop= self.qd0*0 
        # --- Initial conditions (of pertubations)
        dq0  = self.q0  - qop
        dqd0 = self.qd0 - qdop

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
        resLI=sysLI.integrate(time, method='RK45') # **options):

        # --- Convert to dataframe
        dfLI = sysLI.toDataFrame(self.channels, self.FASTDOFScales, x0=qop, xd0=qdop, acc=acc, forcing=forcing, sAcc=self.acc_channels)

        print('-----------------------------------------------------------------------------')
        return resLI, sysLI, dfLI


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

    """
    # TODO TODO TODO  Harmonize with TNSB.py

    # --- Reading main OpenFAST files
    ext     = os.path.splitext(fstFilename)[1]
    FST     = weio.read(fstFilename)
    rootdir = os.path.dirname(fstFilename)
    EDfile  = os.path.join(rootdir,FST['EDFile'].strip('"')).replace('\\','/')
    ED      = weio.read(EDfile)
    rootdir = os.path.dirname(EDfile)
    bldfile = os.path.join(rootdir,ED['BldFile(1)'].strip('"')).replace('\\','/')
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
    if algo=='OpenFAST':
        jxxG=0*m
    else:
        jxxG = m     # NOTE: unknown
        print('>>> windturbine.py: TODO: using unknown jxxG')
    nB = ED['NumBl']
    bld=np.zeros(nB,dtype=object)
    bld[0] = FASTBeamBody(ED, bldFile, Mtop=0, main_axis=main_axis, jxxG=jxxG, spanFrom0=False, bldStartAtRotorCenter=bldStartAtRotorCenter, nSpan=nSpanBld, gravity=gravity, algo=algo) 
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

    # --- Fnd (defined wrt ground/MSL "E")
#     print(FST.keys())
    M_fnd = ED['PtfmMass']
    r_OGfnd_inF = np.array([ED['PtfmCMxt'],ED['PtfmCMyt'],ED['PtfmCMzt']])
    r_OT_inF    = np.array([0             ,0             ,ED['PtfmRefzt']])
    r_TGfnd_inF = -r_OT_inF + r_OGfnd_inF
    fnd = RigidBody('fnd', M_fnd, (ED['PtfmRIner'], ED['PtfmPIner'], ED['PtfmYIner']), s_OG=r_TGfnd_inF, r_O=r_OT_inF) 

    # --- Twr
    twrFile = weio.read(twrfile)
    twr = FASTBeamBody(ED, twrFile, Mtop=M_RNA, main_axis='z', bAxialCorr=False, bStiffening=True, shapes=twrShapes, nSpan=nSpanTwr, algo=algo, gravity=gravity) # TODO options
    twr_rigid  = twr.toRigidBody()
    twr_rigid.pos_global = r_ET_inE


    # --- Full WT rigid body equivalent, with point T as ref
    RNAb = copy.deepcopy(RNA)
    RNAb.pos_global = r_TN_inT+r_ET_inE
    RNAb.R_b2g      = np.eye(3)
    WT_rigid = RNAb.combine(twr_rigid, r_O=r_ET_inE).combine(fnd, r_O=r_ET_inE) # TODO T or Ptfm
    #WT_rigid = RNAb.combine(twr_rigid, r_O=r_EPtfm_inE).combine(fnd, r_O=r_EPtfm_inE) # TODO T or Ptfm

    # --- Moorings
    MAP = None
    if FST['CompMooring']==1:
        from welib.moor.mappp import Map
#         try:
        MAP = Map(fstFilename)
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
    DOFs+=[{'name':'q_SS2'  , 'active':ED['TwSSDOF2'] , 'q0': ED['TTDspSS']  , 'qd0':0 , 'q_channel':'Q_TSS1_[m]', 'qd_channel':'QD_TSS1_[m/s]', 'qdd_channel':'QD2_TSS1_[m/s^2]'}]

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
    WT.MAP        = MAP       # typically at (0,0,0) NOTE: not a body
    WT.RNA        = RNA        # origin at N, rigid body bld+hub+gen+nac
    WT.WT_rigid   = WT_rigid   # rigid wind turbine body, origin at MSL

    WT.DOF= DOFs

    #WT.r_ET_inE = 
    #WT.r_TN_inT
    WT.r_NS_inN = r_NS_inN
    WT.r_NR_inN = r_NR_inN
    WT.r_SR_inS = r_SR_inS
    WT.ED=ED
    WT.FST=FST

    return WT


if __name__ == '__main__':
    np.set_printoptions(linewidth=300, precision=2)
    WT = FASTWindTurbine('../../data/NREL5MW/Main_Onshore.fst')
    print(WT)

