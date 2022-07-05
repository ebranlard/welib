import unittest

import os
import numpy as np
from welib.yams.bodies import *
from welib.yams.utils import *
import welib.weio as weio

MyDir=os.path.dirname(__file__)

# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestUtils(unittest.TestCase):

    def test_inertia(self):
        np.set_printoptions(linewidth=300, precision=2)
        # read ElastoDyn file
        edFile=os.path.join(MyDir,'./../../../data/NREL5MW/offshore/NREL5MW_ED_Offshore.dat')
        ED = weio.FASTInputFile(edFile)

        #def __init__(self, name, mass, J, s_OG, R_b2g=np.eye(3), s_OP=None, r_O=[0,0,0]):
        theta_tilt_y = -ED['ShftTilt']*np.pi/180  # NOTE: tilt has wrong orientation in FAST
        R_NS = R_y(theta_tilt_y)
        r_NS_inN    = np.array([0             , 0, ED['Twr2Shft']]) # Shaft start in N
        r_SR_inS    = np.array([ED['OverHang'], 0, 0             ]) # Rotor center in S
        r_SGhub_inS = np.array([ED['HubCM']   , 0, 0             ]) + r_SR_inS # Hub G in S
        r_NR_inN    = r_NS_inN + R_NS.dot(r_SR_inS)                 # Rotor center in N
        r_NGnac_inN = np.array([ED['NacCMxn'],0,ED['NacCMzn']    ]) # Nacelle G in N
        r_RGhub_inS = - r_SR_inS + r_SGhub_inS

        # --- Hub 
        M_hub  = ED['HubMass']
        JxxHub_atR = ED['HubIner'] #+ ED['GenIner']*ED['GBRatio']**2 
        hub = RigidBody('Hub', M_hub, (JxxHub_atR,0,0), s_OG=r_SGhub_inS, R_b2g=R_NS, s_OP=r_SR_inS, r_O=r_NS_inN) 
        #print(hub)

        # --- Generator 
        gen = RigidBody('Gen', 0, (ED['GenIner']*ED['GBRatio']**2,0,0), s_OG=[0,0,0], R_b2g=R_NS,  r_O=r_NS_inN) 
        #print(gen)

        # --- Nacelle 
#         M_nac = ED['NacMass']
#         JyyNac_atN = ED['NacYIner'] # Inertia of nacelle at N in N
#         nac = RigidBody('Nac', M_nac, (0,JyyNac_atN,0), r_NGnac_inN, s_OP=[0,0,0])
#         print(nac)

        # Blades 
        # TODO reconsider hubrad there
        psi=0
        bld1 = RigidBody('Bld1', 1.684475202e+04, (1.225603507e+07, 1.225603507e+07, 1.684475202e+04), s_OG=[0,0,22.00726], s_OP=[0,0,0], R_b2g=R_x(psi+ 0        ).dot(R_y(ED['PreCone(1)']*np.pi/180))  )
        bld2 = RigidBody('Bld2', 1.684475202e+04, (1.225603507e+07, 1.225603507e+07, 1.684475202e+04), s_OG=[0,0,22.00726], s_OP=[0,0,0], R_b2g=R_x(psi+ 2*np.pi/3).dot(R_y(ED['PreCone(2)']*np.pi/180))  )
        bld3 = RigidBody('Bld3', 1.684475202e+04, (1.225603507e+07, 1.225603507e+07, 1.684475202e+04), s_OG=[0,0,22.00726], s_OP=[0,0,0], R_b2g=R_x(psi+ 4*np.pi/3).dot(R_y(ED['PreCone(3)']*np.pi/180))  )
        #print(bld1)
        #print('Blade mass matrix\n',np.around(bld1.mass_matrix,5))
        blades = bld1.combine(bld2).combine(bld3, name='Blades')
        #print(blades)
        # 
        blades.pos_global = r_NR_inN # Setting origin w.r.t. N
        blades.R_b2g      = R_NS     # Setting frame w.r.t. N, as rotated
        #print(blades)

        # --- Rotor = Hub + Blades
        hubgen= hub.combine(gen, R_b2g=R_NS, r_O=hub.pos_global)
        rotor = blades.combine(hubgen, R_b2g=R_NS, r_O=hub.pos_global)
        #print(rotor)

        rotor = blades.combine(hub).combine(gen, r_O=r_NS_inN, R_b2g=R_NS)
        #print(rotor)
        rotor = blades.combine(hub, r_O=r_NS_inN, R_b2g=R_NS)
        #print(rotor)


if __name__=='__main__':
    unittest.main()
