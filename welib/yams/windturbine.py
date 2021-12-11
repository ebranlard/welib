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
        self.DOF    ={'name':'', 'active':None, 'q0':None, 'qd0':None,'q_channel':None, 'qd_channel':None}

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
        s+=' * DOFname: {} \n'.format(B.DOFname)
        s+=' * q0     : {}\n'.format(B.q0)
        s+=' * qd0    : {}\n'.format(B.qd0)
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
#         for dof in self.DOF:
#             if dof['q_channel'] is None:
#                 chan.append(dof['name'])
#             elif hasattr(dof['q_channel'] , __len__):
# 
#             if 
        chan = [dof['q_channel'] if dof['q_channel'] is not None else dof['name']  for dof in self.DOF if dof['active']]
        chan+= [dof['qd_channel'] if dof['qd_channel'] is not None else 'd'+dof['name'] for dof in self.DOF if dof['active']]
        return chan

    def yams_parameters(self):
        WT=self;
        # --- Dict needed by structural script 
        p = dict()
        p['z_FG']     = WT.fnd.masscenter[2]
        p['M_F']      = WT.fnd.mass
        p['J_xx_F']   = WT.fnd.masscenter_inertia[0,0]
        p['J_yy_F']   = WT.fnd.masscenter_inertia[1,1]
        p['J_zz_F']   = WT.fnd.masscenter_inertia[2,2]
        p['g']        = WT.ED['Gravity']
        p['tilt']     =-WT.ED['ShftTilt']*np.pi/180 # in rad
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
        # One body for all turbine
        p['M_B']      = WT.WT_rigid.mass
        p['x_BG']     = WT.WT_rigid.masscenter[0]
        p['y_BG']     = WT.WT_rigid.masscenter[1]
        p['z_BG']     = WT.WT_rigid.masscenter[2]
        p['J_xx_B']   = WT.WT_rigid.masscenter_inertia[0,0]
        p['J_yy_B']   = WT.WT_rigid.masscenter_inertia[1,1]
        p['J_zz_B']   = WT.WT_rigid.masscenter_inertia[2,2]
        p['J_zx_B']   = WT.WT_rigid.masscenter_inertia[0,2]

        # Mooring restoring
        p['z_TM']      = 0
        p['K_x_M']     = 0
        p['K_y_M']     = 0
        p['K_z_M']     = 0
        p['K_phi_x_M'] = 0
        p['K_phi_y_M'] = 0
        p['K_phi_z_M'] = 0
        # Buoyancy
        p['z_TB']      = 0
        return p

    def simulate_py(self, model, p, time=None, refOut=None):
        """ Perform non-linear simulation based on a model (python package) generated by yams_sympy """
        from welib.system.mech_system import MechSystem

        # --- Reference simulation
        if refOut is not None:
            df=weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
            #time = np.linspace(0,50,5000)
            time = df['Time_[s]'].values
        elif time is None:
            raise Exception('Time vector must be provided if no reference simulation is given')

        # --- Checking package info
        info = model.info()
        nDOFExpected= info['nq']
        print('DOFs:', self.DOFname, 'Model:',info['name'], 'nDOF:',nDOFExpected )
        if len(self.DOFname)!=nDOFExpected:
            raise Exception('Inconsistency in number of DOFs')

        # --- Initial conditions
        q0  = self.q0
        qd0 = self.qd0
        print('q0 :',q0)
        print('qd0:',qd0)


        # --- Non linear
        u=dict()
        for key in info['su']:
            # TODO get it from ref simulation
            u[key]= lambda t: 0
        t=0
        MM      = model.mass_matrix(q0,p)
        forcing = model.forcing(t,q0,qd0,p,u)

        # --- integrate non-linear system
        fM = lambda x: model.mass_matrix(x, p)
        fF = lambda t,x,xd: model.forcing(t, x, xd, p=p, u=u)
        sysNL = MechSystem(fM, F=fF, x0=q0, xdot0=qd0 )
        resNL = sysNL.integrate(time, method='RK45')
        return resNL, sysNL

    def simulate_py_lin(self, model, p, time=None, refOut=None, qop=None, qdop=None):
        """ Perform linear simulation based on a model (python package) generated by yams_sympy """
        from welib.system.mech_system import MechSystem
        # --- Reference simulation
        if refOut is not None:
            df=weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
            #time = np.linspace(0,50,5000)
            time = df['Time_[s]'].values
        elif time is None:
            raise Exception('Time vector must be provided if no reference simulation is given')

        # --- Checking package info
        info = model.info()
        nDOFExpected= info['nq']
        print('DOFs:', self.DOFname, 'Model:',info['name'], 'nDOF:',nDOFExpected )
        if len(self.DOFname)!=nDOFExpected:
            raise Exception('Inconsistency in number of DOFs')

        # --- Operating point
        if qop is None:
            qop = self.q0*0   # TODO
        if qdop is None:
            qdop= self.qd0*0  # TODO
        uop=dict() # Inputs at operating points
        for key in info['su']:
            # TODO get it from ref simulation
            uop[key]= 0 

        # --- Initial conditions (of pertubations)
        q0  = self.q0  - qop
        qd0 = self.qd0 - qdop

        # --- Evaluate linear structural model at operating point
        M_lin   = model.M_lin(qop,p)
        C_lin   = model.C_lin(qop,qdop,p,uop)
        K_lin   = model.K_lin(qop,qdop,p,uop) 
        B_lin   = model.B_lin(qop,qdop,p,uop)

        # --- Integrate linear system
        fF = lambda t,x,xd: np.array([0]*len(qop)) # TODO
        sysLI = MechSystem(M=M_lin, K=K_lin, C=C_lin, F=fF, x0=q0, xdot0=qd0)
        resLI=sysLI.integrate(time, method='RK45') # **options):

        # --- Print linearized mass damping 
        #     print('--------------------')
        #     print('Linear Mass Matrix: ')
        #     print(M_lin)
        #     print('--------------------')
        #     print('Linear Damping Matrix: ')
        #     print(C_lin)
        #     print('--------------------')
        #     print('Linear Stifness Matrix: ')
        #     print(K_lin)
        #     print('--------------------')
        #     print('Linear RHS: ')
        #     print(B_lin)
        return resLI, sysLI


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
def FASTWindTurbine(fstFilename, main_axis='z', nSpanTwr=None, twrShapes=None, nSpanBld=None, algo=''):
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

    r_ET_inE    = np.array([0,0,ED['TowerBsHt']               ]) 
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
    jxxG = m     # NOTE: unknown
    print('>>> windturbine.py: TODO: using unknown jxxG')
    nB = ED['NumBl']
    bld=np.zeros(nB,dtype=object)
    bld[0] = FASTBeamBody(ED, bldFile, Mtop=0, main_axis=main_axis, jxxG=jxxG, spanFrom0=False, nSpan=nSpanBld) 
    for iB in range(nB-1):
        bld[iB+1]=copy.deepcopy(bld[0])
        bld[iB+1].R_b2g
    for iB,B in enumerate(bld):
        B.name='bld'+str(iB+1)
        psi_B= -iB*2*np.pi/len(bld) 
        if main_axis=='x':
            R_SB = R_z(0*np.pi + psi_B) # TODO psi offset and psi0
        elif main_axis=='z':
            R_SB = R_x(0*np.pi + psi_B) # TODO psi0
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
    twr = FASTBeamBody(ED, twrFile, Mtop=M_RNA, main_axis='z', bAxialCorr=False, bStiffening=True, shapes=twrShapes, nSpan=nSpanTwr, algo=algo) # TODO options
    twr_rigid  = twr.toRigidBody()
    twr_rigid.pos_global = r_ET_inE


    # --- Full WT rigid body equivalent, with point T as ref
    RNAb = copy.deepcopy(RNA)
    RNAb.pos_global = r_TN_inT+r_ET_inE
    RNAb.R_b2g      = np.eye(3)
    WT_rigid = RNAb.combine(twr_rigid, r_O=r_ET_inE).combine(fnd, r_O=r_ET_inE)

    # --- Degrees of freedom
    DOFs=[]
    DOFs+=[{'name':'x'      , 'active':ED['PtfmSgDOF'], 'q0': ED['PtfmSurge']  , 'qd0':0 , 'q_channel':'PtfmSurge_[m]' , 'qd_channel':'QD_Sg_[m/s]'}]
    DOFs+=[{'name':'y'      , 'active':ED['PtfmSwDOF'], 'q0': ED['PtfmSway']   , 'qd0':0 , 'q_channel':'PtfmSway_[m]'  , 'qd_channel':'QD_Sw_[m/s]'}]
    DOFs+=[{'name':'z'      , 'active':ED['PtfmHvDOF'], 'q0': ED['PtfmHeave']  , 'qd0':0 , 'q_channel':'PtfmHeave_[m]' , 'qd_channel':'QD_Hv_[m/s]'}]

    DOFs+=[{'name':'\phi_x' , 'active':ED['PtfmRDOF'] , 'q0': ED['PtfmRoll']*np.pi/180  , 'qd0':0 , 'q_channel':'PtfmRoll_[deg]'  , 'qd_channel':'QD_R_[rad/s]'}]
    DOFs+=[{'name':'\phi_y' , 'active':ED['PtfmPDOF'] , 'q0': ED['PtfmPitch']*np.pi/180 , 'qd0':0 , 'q_channel':'PtfmPitch_[deg]' , 'qd_channel':'QD_P_[rad/s]'}]
    DOFs+=[{'name':'\phi_z' , 'active':ED['PtfmYDOF'] , 'q0': ED['PtfmYaw']*np.pi/180   , 'qd0':0 , 'q_channel':'PtfmYaw_[deg]'   , 'qd_channel':'QD_Y_[rad/s]'}]

    DOFs+=[{'name':'q_FA1'  , 'active':ED['TwFADOF1'] , 'q0': ED['TTDspFA']  , 'qd0':0 , 'q_channel':'Q_TFA1_[m]', 'qd_channel':'QD_TFA1_[m/s]'}]
    DOFs+=[{'name':'q_SS1'  , 'active':ED['TwSSDOF1'] , 'q0': ED['TTDspSS']  , 'qd0':0 , 'q_channel':'Q_TSS1_[m]', 'qd_channel':'QD_TSS1_[m/s]'}]
    DOFs+=[{'name':'q_FA2'  , 'active':ED['TwFADOF2'] , 'q0': ED['TTDspFA']  , 'qd0':0 , 'q_channel':'Q_TFA2_[m]', 'qd_channel':'QD_TFA2_[m/s]'}]
    DOFs+=[{'name':'q_SS2'  , 'active':ED['TwSSDOF2'] , 'q0': ED['TTDspSS']  , 'qd0':0 , 'q_channel':'Q_TSS1_[m]', 'qd_channel':'QD_TSS1_[m/s]'}]

    DOFs+=[{'name':'\\theta_y','active':ED['YawDOF']  , 'q0': ED['NacYaw']*np.pi/180   , 'qd0':0 ,          'q_channel':'NacYaw_[deg]' , 'qd_channel':'QD_Yaw_[rad/s]'}]
    DOFs+=[{'name':'\\psi'    ,'active':ED['GenDOF']  , 'q0': ED['Azimuth']*np.pi/180  , 'qd0':ED['RotSpeed']*2*np.pi/60 , 'q_channel':'Azimuth_[deg]', 'qd_channel':'RotSpeed_[rpm]'}]

    DOFs+=[{'name':'\\nu'     ,'active':ED['DrTrDOF'] , 'q0': 0  , 'qd0':0 , 'q_channel':'Q_DrTr_[rad]', 'qd_channel':'QD_DrTr_[rad/s]'}]

    DOFs+=[{'name':'q_Fl1'  , 'active':ED['FlapDOF1'] , 'q0': ED['OOPDefl']  , 'qd0':0 , 'q_channel':'Q_B1F1_[m]', 'qd_channel':'QD_B1F1_[m/s]'}]
    DOFs+=[{'name':'q_Ed1'  , 'active':ED['EdgeDOF']  , 'q0': ED['IPDefl']   , 'qd0':0 , 'q_channel':'Q_B1E1_[m]', 'qd_channel':'QD_B1E1_[m/s]'}]

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
    WT.RNA        = RNA        # origin at N, rigid body bld+hub+gen+nac
    WT.WT_rigid   = WT_rigid   # rigid wind turbine body, origin at MSL

    WT.DOF= DOFs

    #WT.r_ET_inE = 
    #WT.r_TN_inT
    WT.r_NS_inN = r_NS_inN
    WT.r_NR_inN = r_NR_inN
    WT.r_SR_inS = r_SR_inS
    WT.ED=ED

    return WT


if __name__ == '__main__':
    np.set_printoptions(linewidth=300, precision=2)
    WT = FASTWindTurbine('../../data/NREL5MW/Main_Onshore.fst')
    print(WT)

