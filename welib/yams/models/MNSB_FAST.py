"""
MTNSB refer to : Foundation (monopile or Jacket) Tower Nacelle, Shaft, Blades

This scripts provides some helper functions to simulates such an assembly of bodies using the Rayleigh-Rizt approximation and joint coordinates. 

The theory is provided in the reference below. The article also contains an example in the its section, which is reproduced in the file test_TNSB.py.


Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

##
import numpy as np
import copy
import os


from welib.yams.windturbine import FASTWindTurbine
from welib.yams.yams_rec import YAMSRecFASTBeamBody
from welib.yams.models.TNSB import TNSBStructure

from welib.weio.fast_input_deck import FASTInputDeck



# --------------------------------------------------------------------------------}
# --- Creating a FNSB model from a FAST model
# --------------------------------------------------------------------------------{
class FASTmodel2MNSB(FASTWindTurbine):

    def __init__(self, FST_file, 
                 shapes_sub=[0,4], nSpan_sub=None,
                 shapes_bld=[], nSpan_bld=None,
                 bHubMass=1, bNacMass=1, bBldMass=1, 
                 DEBUG=False, 
                 main_axis ='x', bStiffening=True, assembly='manual', q=None, bTiltBeforeNac=False,
                 fixedShaft=False,
                 spanFrom0=True, # TODO for legacy, we keep this for now..
                 bladeMassExpected=None,
                 gravity=None,
                 algo='', # TODO replace with OpenFAST
                 FEM_method=None,                 
                 SD_bOverride=True,
                 SD_bCI = False,
                 verbose=False
	    ):
        """ 
        Returns the following structure
          WT.fnd :  BeamBody
          WT.shft:  RigiBody
          WT.nac :  RigidBody
          WT.bld :  List of BeamBodies

          MM, KK, DD : mass, stiffness and damping matrix of full system

          shapes_sub: 0:ux, 1:uy, 2:uz, 3:vx, 4:vy, 5:vz 
                  E.g. for surge and pitch (ux, vy): [0,4]


          NOTE/TODO: compare this with "windturbine.py"
        """
        shapes_twr=[]

        self.shapes_sub = shapes_sub # we store fo convenience
        self.shapes_twr = shapes_twr # we store fo convenience
        self.shapes_bld = shapes_bld # we store fo convenience

        # --- Defining a default MNSB structure
        WT = TNSBStructure(
                         main_axis=main_axis, 
                         bTiltBeforeNac=bTiltBeforeNac,
                         )
        # --- Calling Parent with that structure 
        FASTWindTurbine.__init__(self, WT=WT,
                                 main_axis=main_axis, 
                                 algo=algo)

        # --- Read fast input files
        readlist = ['Fst', 'ED', 'EDtwr', 'EDbld', 'SD','HD']
        self.loadFST(FST_file, readlist=readlist)
        self.setGravity(gravity)
        #if SD is None:
        #    raise Exception('Couldnt read SubDyn file')
        #if bld is None:
        #    raise Exception('Couldnt read blade file')

        # --- Reading SubDyn file
        self.setupSDInit()
        zBot = self.SD.zBot

        # --- Default arguments (needs ED loaded)
        #nSpan_twr = self._defaultNSpanTwr(nSpan_twr, verbose=verbose)
        nSpan_bld = self._defaultNSpanBld(nSpan_bld, verbose=verbose)

        # --------------------------------------------------------------------------------}
        ## --- Creating bodies
        # --------------------------------------------------------------------------------{
        ## --- Strucural and geometrical Inputs
        self.setupEDGeom(zBot=zBot, bTiltBeforeNac=bTiltBeforeNac, flavor='monopile_is_tower')
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
        # self.setupEDTwr(shapes=shapes_twr, nSpan=nSpan_twr, 
        #                 bStiffening=bStiffening,
        #                 flavor='yams_rec')
        # Tower Body
        #   None for now

        # --- FND body
        Mtop = self.WT.RNA.mass
        self.setupSD(shapes=shapes_sub, nSpan=nSpan_sub,
                     Mtop = Mtop,
                     bStiffening=bStiffening,
                     bOverride = SD_bOverride,
                     bCI       = SD_bCI,
                     FEM_method=FEM_method,
                     flavor='yams_rec'
                     )

        self.WT.twr = self.WT.fnd # TNSBStructure will ignore fnd
        self.WT.fnd = None

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

        if assembly=='manual':
            self.WT.manual_assembly(q=q, DEBUG=DEBUG, fixedShaft=fixedShaft)
        else:
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
    FASTmodel2FNSB()
    pass
