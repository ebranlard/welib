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
# --- Converters 
# --------------------------------------------------------------------------------{
def FASTWindTurbine(fstFilename):

    # --- Reading main OpenFAST files
    ext=os.path.splitext(fstFilename)[1]
    FST=weio.read(fstFilename)
    rootdir = os.path.dirname(fstFilename)
    EDfile = os.path.join(rootdir,FST['EDFile'].strip('"')).replace('\\','/')
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

    # Blades 
    bldFile = weio.read(bldfile)

    m= bldFile['BldProp'][:,3]  
    jxxG= m # NOTE: unknown
    bld=FASTBeamBody(ED, bldFile, Mtop=0, main_axis='z', jxxG=jxxG)
    print(bld.MM)
    # TODO reconsider hubrad there
    psi=0
    bld1 = RigidBody('Bld1', 1.684475202e+04, (1.225603507e+07, 1.225603507e+07, 1.684475202e+04), s_OG=[0,0,22.00726], s_OP=[0,0,0], R_b2g=R_x(psi+ 0        ).dot(R_y(ED['PreCone(1)']*np.pi/180))  )
    bld2 = RigidBody('Bld2', 1.684475202e+04, (1.225603507e+07, 1.225603507e+07, 1.684475202e+04), s_OG=[0,0,22.00726], s_OP=[0,0,0], R_b2g=R_x(psi+ 2*np.pi/3).dot(R_y(ED['PreCone(2)']*np.pi/180))  )
    bld3 = RigidBody('Bld3', 1.684475202e+04, (1.225603507e+07, 1.225603507e+07, 1.684475202e+04), s_OG=[0,0,22.00726], s_OP=[0,0,0], R_b2g=R_x(psi+ 4*np.pi/3).dot(R_y(ED['PreCone(3)']*np.pi/180))  )
    print(bld1)
    print('Blade mass matrix\n',np.around(bld1.mass_matrix,5))
    blades = bld1.combine(bld2).combine(bld3, name='Blades')
    print(blades)
    # 
    blades.pos_global = r_NR_inN # Setting origin w.r.t. N
    blades.R_b2g      = R_NS     # Setting frame w.r.t. N, as rotated
    print(blades)

    # --- Rotor = Hub + Blades
    hubgen= hub.combine(gen, R_b2g=R_NS, r_O=hub.pos_global)
    rotor = blades.combine(hubgen, R_b2g=R_NS, r_O=hub.pos_global)
    print(rotor)

    rotor = blades.combine(hub).combine(gen, r_O=r_NS_inN, R_b2g=R_NS)
    print(rotor)
    rotor = blades.combine(hub, r_O=r_NS_inN, R_b2g=R_NS)
    print(rotor)

    print(hub)
    print(gen)
    print(nac)
    # --- Fnd (defined wrt ground/MSL "E")
    print(FST.keys())
    M_fnd = ED['PtfmMass']
    r_OGfnd_inF = [ED['PtfmCMxt'],ED['PtfmCMyt'],ED['PtfmCMzt']]
    r_OT_inF    = [0             ,0             ,ED['PtfmRefzt']]
    r_TGfnd_inF = r_OT_inF + r_OGfnd_inF
    fnd = RigidBody('fnd', M_fnd, (ED['PtfmRIner'], ED['PtfmPIner'], ED['PtfmYIner']), s_OG=r_TGfnd_inF, r_O=r_OT_inF) 

    M_RNA= rotor.mass + M_nac

    # --- Twr
    twrFile = weio.read(twrfile)
    twr = FASTBeamBody(ED, twrFile, Mtop=M_RNA, main_axis='z', bAxialCorr=False, bStiffening=True) # TODO options

#     z   = twr['HtFract_[-]']*(ED['TowerHt']-ED['TowerBsHt'])
#     m   = twr['TMassDen_[kg/m]']  
#     nSpan = len(z)
#     PhiU         = np.zeros((4,3,nSpan))     # Shape functions
#     PhiU[0][0,:] = twr['ShapeForeAft1_[-]']  # along x
#     PhiU[1][0,:] = twr['ShapeForeAft2_[-]']  # along x
#     PhiU[2][1,:] = twr['ShapeSideSide1_[-]'] # along y
#     PhiU[3][1,:] = twr['ShapeSideSide2_[-]'] # along y
#     s_G      = np.zeros((3,nSpan))       # COG location
#     s_G[2,:] = z
#     jxxG= z*0 + m # NOTE: unknown
#     MM, Gr, Ge, Oe, Oe6 = GMBeam(s_G, z, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis='z', split_outputs=False, rot_terms=True)

    # --- Return
    WT = WindTurbineStructure()
    WT.hub = hub
    WT.gen = gen
    WT.nac = nac
    WT.twr = twr
    #WT.r_ET_inE
    #WT.r_TN_inT
    #WT.r_NS_inN
    #WT.r_SR_inS

    return WT


if __name__ == '__main__':
    np.set_printoptions(linewidth=300, precision=2)
    FASTWindTurbine('../../data/NREL5MW/Main_Onshore_OF2.fst')
