"""
MTNSB refer to : Foundation (monopile or Jacket) Tower Nacelle, Shaft, Blades

This scripts provides some helper functions to simulates such an assembly of bodies using the Rayleigh-Rizt approximation and joint coordinates. 

The theory is provided in the reference below. The article also contains an example in the its section, which is reproduced in the file test_TNSB.py.


Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

##
import numpy as np
import os

from welib.yams.yams_rec import YAMSRecGroundBody
from welib.yams.rotations import R_x, R_y, R_z

from welib.yams.windturbine import rigidBlades
from welib.yams.windturbine import WindTurbineStructure, FASTWindTurbine


class MTNSBStructure(WindTurbineStructure):
    def __init__(self, 
                 main_axis='x',
                 bTiltBeforeNac=False,
                 ):
        #, main_axis='x', theta_tilt=0, theta_yaw=0, theta_cone=0, bTiltBeforeNac=False):
        WindTurbineStructure.__init__(self)	

        # --- Additional properties not in Parent class
        self.main_axis      = main_axis
        self.nac_yaw        = 0 # TODO
        self.bTiltBeforeNac = bTiltBeforeNac
        self.additional_properties=[] # for user output in __repr__, so I remember what we have in the object


    def __repr__(self):
        s='<TNSB {} object> with fields:\n'.format(type(self).__name__)
        s+=f' - grd twr yaw nac hubgen bld: RigidBody or FASTBeamBody\n'
        s+=f' - MM KK DD: matrices\n'
        if self.q is not None:
            s+=f' - q : {self.q.flatten()}\n'
        s+=f' - r_EF_inE: {self.r_EF_inE.flatten()}\n'
        s+=f' - r_FT_inF: {self.r_FT_inF.flatten()}\n'
        s+=f' - r_ET_inE: {self.r_ET_inE.flatten()}\n'
        s+=f' - r_TN_inT: {self.r_TN_inT.flatten()}\n'
        s+=f' - r_NS_inN: {self.r_NS_inN.flatten()}\n'
        s+=f' - r_SR_inS: {self.r_SR_inS.flatten()}\n'
        s+=f' - main_axis     : {self.main_axis}\n'
        s+=f' - shaft_tilt    : {self.shaft_tilt*180/np.pi} [deg] (but stored in rad)\n'
        s+=f' - blade_cone    : {self.blade_cone*180/np.pi} [deg] (but stored in rad)\n'
        s+=f' - nac_yaw       : {self.nac_yaw  *180/np.pi}  [deg] (but stored in rad)\n'
        s+=f' - bTiltBeforeNac: {self.bTiltBeforeNac}\n'
        s+='----------------------------------------------------------------\n'
        s+=f' * Origin F  : {self.fnd.pos_global.T}\n'
        s+=f' * Origin T  : {self.twr.pos_global.T}\n'
        s+=f' * Origin N  : {self.nac.pos_global.T}\n'
        s+=f' * Origin R  : {self.bld[0].pos_global.T}\n'
        s+=f' * Origin S  : {self.hubgen.pos_global.T}\n'
        s+='----------------- RNA ---------------------------------------\n'
        # TODO
        s+=f'M_RNA       {self.RNA.mass:.4f}\n'
#         s+=f'r_NGrna_inN {np.asarray(self.r_NGrna_inN).flatten()}\n'
#         s+=f'     r_NGnac_inN {np.asarray(self.r_NGnac_inN).flatten().round(4)} M_nac {self.nac.mass:.4f}\n'
#         s+=f'     r_NGhub_inN {np.asarray(self.r_NGhub_inN).flatten().round(4)} M_hub {self.hubgen.mass:.4f}\n'
#         s+=f'     r_NGrot_inN {np.asarray(self.r_NGrot_inN).flatten().round(4)} M_rot {self.M_rot:.4f}\n'
        s+='---------------Couplings ---------------------------------------\n'
        s+='Constant: (Bhat_t)\n'
        s+=str(self.twr.Bhat_t_bc)+'\n'
        s+='Time varying:\n'
        s+=str(self.twr.alpha_couplings)+'\n' # Time varying function of twr.gzf
        s+='---------------More --------------------------------------------\n'
        s+= ' - Additional Props: {}\n'.format(self.additional_properties)
        return s

    # --------------------------------------------------------------------------------}
    # --- Creating a TNSB model automatically 
    # --------------------------------------------------------------------------------{
    def auto_assembly(self, q=None, DEBUG=False, fixedShaft=False):
        if q is None:
            raise NotImplementedError('Figure out DOF dimension, zero q vector in windturbine')

        # TODO gen
        # --- Aliases
        sft = self.hubgen
        gen = self.gen
        nac = self.nac
        yaw = self.yawBr
        bld = self.bld
        twr = self.twr
        fnd = self.fnd # Monopile
        theta_yaw    =  self.nac_yaw
        r_EF_inE     =  self.r_EF_inE
        r_FT_inF     =  self.r_FT_inF
        r_ET_inE     =  self.r_ET_inE
        r_TN_inT     =  self.r_TN_inT    
        r_NS_inN     =  self.r_NS_inN   
        r_SR_inS     =  self.r_SR_inS    

        # --- Connection transformation matrices 
        if self.main_axis=='x':
            #R_NS     = np.dot(R_y(-tilt_up),R_z(q_psi + np.pi)) # << tilt 
            if self.bTiltBeforeNac:
                R_cn0 = np.dot(R_x (self.nac_yaw) , R_y (self.shaft_tilt))
                R_cs0 = R_z (np.pi)
            else:
                R_cn0 = R_x (self.nac_yaw) 
                R_cs0 = np.dot( R_y(self.shaft_tilt) , R_z (np.pi)) # Note: OrientBefore
            Shaft_axis='z'
        if self.main_axis=='z':
            if self.bTiltBeforeNac:
                R_cn0 = np.dot(R_z (self.nac_yaw) , R_y(self.shaft_tilt))
                R_cs0 = R_x (np.pi)
            else:
                R_cn0 = R_z (self.nac_yaw) 
                R_cs0 = np.dot(R_y(self.shaft_tilt) , R_x (np.pi) )# Note: OrientBefore
            Shaft_axis='x'
        nB=len(bld)


        # Creating reference frame
        grd = YAMSRecGroundBody()

        # Connections between bodies
        #grd.connectTo(twr, Point=r_ET_inE, Type='Rigid')
        grd.connectTo(fnd, Point=r_EF_inE, Type='Rigid')
        fnd.connectTo(twr, Point=r_FT_inF, Type='Rigid')
        twr.connectTo(nac, Point=r_TN_inT, Type='Rigid', RelOrientation = R_cn0 , OrientAfter=True)
        twr.connectTo(yaw, Point=r_TN_inT, Type='Rigid', RelOrientation = R_cn0 , OrientAfter=True)
        if fixedShaft:
            nac.connectTo (sft   , Point=r_NS_inN, Type='Rigid', RelOrientation = R_cs0, OrientAfter=False)
        else:
            nac.connectTo (sft   , Point=r_NS_inN, Type='SphericalJoint',JointRotations=[Shaft_axis],RelOrientation = R_cs0, OrientAfter=False)
        for i,B in enumerate(bld):
            psi_B= -i*2*np.pi/nB # 0 -2pi/2 2pi/3  or 0 pi
            if self.main_axis=='x':
                R_SB = R_z(0*np.pi + psi_B)
            elif self.main_axis=='z':
                R_SB = R_x(0*np.pi + psi_B)
            R_SB = np.dot(R_SB, R_y(self.blade_cone))
            sft.connectTo(B, Point=r_SR_inS, Type='Rigid', RelOrientation = R_SB)

        # Setting DOF index for all bodies and connections 
        nq = grd.setupDOFIndex();
        if nq!=len(q):
           print('>>> ',nq,len(q))
           raise Exception('Wrong number of dof')

        # TODO TODO
        self.grd  = grd # TODO?
        self.fixedShaft = fixedShaft # TODO?

        self.updateKinematics(q)
# 

    def updateKinematics(self, q):
        q = np.asarray(q).reshape((len(q),1))
        print(q)
        self.grd.updateChildrenKinematicsNonRecursive(q)
        self.fnd.updateChildrenKinematicsNonRecursive(q)
        self.twr.updateChildrenKinematicsNonRecursive(q)
        self.yawBr.updateChildrenKinematicsNonRecursive(q)
        self.nac.updateChildrenKinematicsNonRecursive(q)
        self.hubgen.updateChildrenKinematicsNonRecursive(q)

        # --- Full system
        nq = len(q)
        MM = self.grd.M
        KK = self.grd.K
        DD = self.grd.D

        MM[np.abs(MM)< 1e-09] = 0

        # --- returning everthin in a structure class
        # TODO
        self.MM   = MM
        self.KK   = KK
        self.DD   = DD
        self.q    = q
        #self.init_trigger()



# --------------------------------------------------------------------------------}
# --- Creating a FNSB model from a FAST model
# --------------------------------------------------------------------------------{
class FASTmodel2MTNSB(FASTWindTurbine):

    def __init__(self, FST_file, 
                 shapes_sub=[0,4], nSpan_sub=None,
                 shapes_twr=None, nSpan_twr=None,
                 shapes_bld=None, nSpan_bld=None, 
                 bHubMass=1, bNacMass=1, bBldMass=1, 
                 bStiffening=True, bTiltBeforeNac=False, spanFrom0=True, # TODO for legacy, we keep this for now..
                 DEBUG=False, verbose=False,
                 main_axis ='x',
                 q=None, 
                 fixedShaft=False,
                 bladeMassExpected=None,
                 gravity=None,
                 algo='', # TODO replace with OpenFAST
                 ):
        """ 
        Returns the following structure
          WT.mnp :  BeamBody
          WT.twr :  BeamBody
          WT.shft:  RigiBody
          WT.nac :  RigidBody
          WT.blds:  List of BeamBodies

          MM, KK, DD : mass, stiffness and damping matrix of full system

          shapes_sub: 0:ux, 1:uy, 2:uz, 3:vx, 4:vy, 5:vz 
                  E.g. for surge and pitch (ux, vy): [0,4]
        """
        if shapes_sub is None:
            shapes_sub=[]
        if shapes_twr is None:
            shapes_twr=[0,1]
        if shapes_bld is None:
            shapes_bld=[]
        self.shapes_sub = shapes_sub # we store fo convenience
        self.shapes_twr = shapes_twr # we store fo convenience
        self.shapes_bld = shapes_bld # we store fo convenience

        # --- Defining a default MNSB structure
        WT = MTNSBStructure(
                         main_axis=main_axis, 
                         bTiltBeforeNac=bTiltBeforeNac,
                         )
        # --- Calling Parent with that structure 
        FASTWindTurbine.__init__(self, WT=WT,
                                 main_axis=main_axis, 
                                 algo=algo)

        # --- Read fast input files
        readlist = ['Fst', 'ED', 'EDtwr', 'EDbld', 'SD']
        self.loadFST(FST_file, readlist=readlist)
        self.setGravity(gravity)

        # --- Reading SubDyn file
        self.setupSDInit()
        zBot = self.SD.zBot

        # --- Default arguments (needs ED loaded)
        nSpan_twr = self._defaultNSpanTwr(nSpan_twr, verbose=verbose)
        nSpan_bld = self._defaultNSpanBld(nSpan_bld, verbose=verbose)

        # --------------------------------------------------------------------------------}
        ## --- Creating bodies
        # --------------------------------------------------------------------------------{
        ## --- Strucural and geometrical Inputs
        self.setupEDGeom(zBot=zBot, bTiltBeforeNac=bTiltBeforeNac, flavor='')
        # --- Sft = Hub + Gen
        self.setupEDHubGen(bHubMass=bHubMass, flavor='yams_rec') 
        # --- Gen only
        self.setupEDGen(flavor='yams_rec')
        # --- Nac
        self.setupEDNac(bNacMass=bNacMass, flavor='yams_rec')
        # --- Yaw
        self.setupEDYaw(flavor='yams_rec')
        # --- Hub
        self.setupEDHub() #flavor='yams_rec')
        # WT.bld
        self.setupEDBld(shapes=shapes_bld, nSpan=nSpan_bld,
                        spanFrom0=spanFrom0, massExpected=bladeMassExpected,
                        flavor='yams_rec')
        self.setupEDRot()   # WT.rot and rotgen  (Generic Rigid Body)
        self.setupEDRNA()   # WT.RNA             (Generic Rigid Body)
        #--------------------------- HUB NAC YAW RNA COMMON WITH FTNSB 
        # WT.twr
        self.setupEDTwr(shapes=shapes_twr, nSpan=nSpan_twr, 
                        bStiffening=bStiffening,
                        flavor='yams_rec')
        # --- FND body
        Mtop = self.WT.RNA.mass
        self.setupSD(shapes=shapes_sub, nSpan=nSpan_sub,
                     Mtop = Mtop,
                     bStiffening=bStiffening,
                     flavor='yams_rec')

        # --------------------------------------------------------------------------------}
        # --- Initial conditions and DOFs
        # --------------------------------------------------------------------------------{
        # --- Initial conditions
        self.setupEDDOFs()
        self.setActiveDOFs(shapes_sub=self.shapes_sub, shapes_twr=self.shapes_twr, shapes_bld=shapes_bld, fixedShaft=fixedShaft, verbose=verbose)
        nDOF = len(self.WT.q0) # 1 + len(shapes_twr) + len(shapes_bld) * nB # +1 for Shaft
        if DEBUG:
            print('Initial conditions:')
            print(self.WT.q0)
            print(self.WT.qd0)
            print(self.WT.z0)
        nDOF = len(self.WT.q0) # 1 + len(shapes_twr) + len(shapes_bld) * nB # +1 for Shaft
        if q is None:
            q = np.zeros((nDOF,1)) # Only pos, not vel here.
        # --------------------------------------------------------------------------------}
        # --- Assembly 
        # --------------------------------------------------------------------------------{
        self.WT.auto_assembly(q=q, DEBUG=DEBUG, fixedShaft=fixedShaft)

        # --- Useful data
        WT=self.WT
        WT.DCK = self.DCK
        WT.FST = self.FST
        WT.ED  = self.ED
        WT.SD  = self.SD
        WT.Hydro     = self.FST['CompHydro']>0 # FST['CompSeaSt']>0 and 
        WT.HD        = self.DCK.fst_vt['HydroDyn']
        try:
            WT.WtrDens   = self.FST['WtrDens']
            WT.WtrDpth   = self.FST['WtrDpth']
        except:
            WT.WtrDens   = self.HD['WtrDens']
            WT.WtrDpth   = self.HD['WtrDpth']

        WT.DampMat       = self.SD.DampMat
        WT.RayleighCoeff = self.SD.RayleighCoeff
        WT.additional_properties +=['DCK', 'FST', 'ED', 'DampMat', 'RayleighCoeff', 'WaterDepth','Hydro']






if __name__=='__main__':
    from welib.yams.models.TNSB_FAST import FASTmodel2TNSB
    fstFile = 'C:/Users/ebranlard/Documents/Work/2024-10-OESI-Digitwin/DigiTwinMonopile/code5_wt/simulations_wt/06_Jonswap/OF_F3T1S1_H1A0_Hs=2.5_Tp=10.fst'

    print('=============================================================================')
    print('--------------- REF ----------------')
    q = np.array([[10],[0]])
    WT = FASTmodel2TNSB(fstFile, q=q, shapes_twr=[0], shapes_bld=[], main_axis='z', assembly='auto').WT
    print(WT)

    print('=============================================================================')
    print('--------------- REF ----------------')
    q = [0, 0.0, 0, 0.000]
    WT = FASTmodel2MTNSB(fstFile, q=q, shapes_sub=[0,4], shapes_twr=[0], main_axis='z').WT
    print(WT)
    print('=============================================================================')
    print('--------------- SURGE---------------')
    q = [10, 0.0, 0, 0.000]
    WT.updateKinematics(q)
    print(WT)
    print('=============================================================================')
    print('--------------- PITCH---------------')
    q = [0, 0.5, 0, 0.000]
    WT.updateKinematics(q)
    print(WT)
    print('=============================================================================')
    print('--------------- TOWER---------------')
    q = [0, 0.0, 10, 0.000]
    WT.updateKinematics(q)
    print(WT)

    
