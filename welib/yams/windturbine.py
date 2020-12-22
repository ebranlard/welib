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
from welib.yams.utils import R_x, R_y, R_z
from welib.yams.bodies import RigidBody, FlexibleBody, FASTBeamBody

class WindTurbineStructure():
    def __init__(self):
        self.hub = None
        self.gen = None
        self.nac = None
        self.twr = None

    @staticmethod
    def fromFAST(fstFilename):
        return FASTWindTurbine(fstFilename)


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
def FASTWindTurbine(fstFilename, main_axis='z'):
    """

    """
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
    nB = ED['NumBl']
    bld=np.zeros(nB,dtype=object)
    bld[0] = FASTBeamBody(ED, bldFile, Mtop=0, main_axis=main_axis, jxxG=jxxG, spanFrom0=False) 
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
    blades = rigidBlades(bld, r_O = [0,0,0])
    blades.pos_global = r_NR_inN
    blades.R_b2g      = R_NS

    # --- Rotor = Hub + Blades (with origin R, using N as global ref)
    rot = blades.combine(hub, R_b2g=R_NS, r_O=blades.pos_global)
    rot.name='rotor'

    # --- RNA
    RNA = rot.combine(gen).combine(nac,r_O=[0,0,0])
    RNA.name='RNA'
    print(RNA)

    # --- RNA
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
    twr = FASTBeamBody(ED, twrFile, Mtop=M_RNA, main_axis='z', bAxialCorr=False, bStiffening=True) # TODO options

    # --- Return
    WT = WindTurbineStructure()
    WT.hub = hub # origin at S
    WT.gen = gen # origin at S
    WT.nac = nac # origin at N
    WT.twr = twr # origin at T
    WT.fnd = fnd # origin at T
    WT.bld = bld # origin at R
    WT.rot = rot # origin at R, rigid body bld+hub
    WT.RNA = RNA # origin at N, rigid body bld+hub+gen+nac
    #WT.r_ET_inE
    #WT.r_TN_inT
    #WT.r_NS_inN
    #WT.r_SR_inS

    return WT


if __name__ == '__main__':
    np.set_printoptions(linewidth=300, precision=2)
    FASTWindTurbine('../../data/NREL5MW/Main_Onshore_OF2.fst')
