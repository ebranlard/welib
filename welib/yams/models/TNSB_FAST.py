##
import numpy as np
import copy
import matplotlib.pyplot as plt
import os

from welib.yams.windturbine import FASTWindTurbine
from welib.yams.yams_rec import YAMSRecFASTBeamBody, YAMSRecRigidBody
from welib.yams.utils import *
from welib.yams.models.TNSB import TNSBStructure

import welib.weio as weio
from welib.weio.fast_input_file import FASTInputFile
from welib.weio.fast_input_deck import FASTInputDeck

# --------------------------------------------------------------------------------}
# --- Creating a TNSB model from a FAST model
# --------------------------------------------------------------------------------{
# TODO TODO TODO
# TODO TODO TODO HARMONIZE WITH WINDTURBINE.PY AND TNSB..
# TODO TODO TODO
class FASTmodel2TNSB(FASTWindTurbine):
    """ 
    Constructor for a TNSB Wind turbine Structure
    """
    
    def __init__(self, FST_file,nB=3, shapes_twr=None, shapes_bld=None, 
                   nSpan_twr=None, nSpan_bld=None, 
                   bHubMass=1, bNacMass=1, bBldMass=1, 
                   DEBUG=False, main_axis ='x', bStiffening=True, assembly='manual', q=None, bTiltBeforeNac=False,
                   spanFrom0=True, # TODO for legacy, we keep this for now..
                   bladeMassExpected=None,
                   gravity=None,
                   algo='', # TODO replace with OpenFAST
                   verbose=False
        ):
        """ 
        Returns the following structure
          WT.twr :  BeamBody
          WT.shft:  RigiBody
          WT.nac :  RigidBody
          WT.bld:  List of BeamBodies

          MM, KK, DD : mass, stiffness and damping matrix of full system


          NOTE/TODO: compare this with "windturbine.py"
        """
        if shapes_twr is None:
            shapes_twr=[0,1]
        if shapes_bld is None:
            shapes_bld=[]

        # --- Defining a default TNSB structure
        WT = TNSBStructure(
                         main_axis=main_axis, 
                         bTiltBeforeNac=bTiltBeforeNac,
                         )

        # --- Calling Parent with that structure 
        FASTWindTurbine.__init__(self, WT=WT,
                                 main_axis=main_axis, 
                                 algo=algo)


        # --- Read fast input files
        readlist = ['Fst', 'ED', 'EDtwr', 'EDbld']
        self.loadFST(FST_file, readlist=readlist)

        # --- Default arguments (needs ED loaded)
        nSpan_twr = self._defaultNSpanTwr(nSpan_twr, verbose=verbose)
        nSpan_bld = self._defaultNSpanBld(nSpan_bld, verbose=verbose)

        # --------------------------------------------------------------------------------}
        ## --- Creating bodies
        # --------------------------------------------------------------------------------{
        ## --- Strucural and geometrical Inputs
        self.setupEDGeom(zBot=0, bTiltBeforeNac=bTiltBeforeNac)
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
        #print('Stiffnening', bStiffening)
        #print('Ttw.KKg   \n', Twr.KKg[6:,6:])
        if DEBUG:
            self.setupDebug()
        # --------------------------------------------------------------------------------}
        # --- Assembly 
        # --------------------------------------------------------------------------------{
        nDOF = 1 + len(shapes_twr) + len(shapes_bld) * nB # +1 for Shaft
        if q is None:
            q = np.zeros((nDOF,1)) # TODO, full account of q not done

        if assembly=='manual':
            self.WT.manual_assembly(q=q, DEBUG=DEBUG)
        else:
            self.WT.auto_assembly(q=q, DEBUG=DEBUG)

        # --- Initial conditions
        ED = self.ED
        omega_init = ED['RotSpeed']*2*np.pi/60 # rad/s
        psi_init   = ED['Azimuth']*np.pi/180   # rad
        FA_init    = ED['TTDspFA']
        iPsi     = self.WT.iPsi
        nDOFMech = len(self.WT.MM)
        q_init   = np.zeros(2*nDOFMech) # x2, state space

        if len(shapes_twr)>0:
            q_init[0] = FA_init

        q_init[iPsi]          = psi_init
        q_init[nDOFMech+iPsi] = omega_init

        self.WT.q_init = q_init
        if DEBUG:
            print('Initial conditions:')
            print(q_init)

        # --- Useful data
        WT=self.WT
        self.WT.ED=ED


# --------------------------------------------------------------------------------}
# --- Read Relevant fields from an outb file 
# --------------------------------------------------------------------------------{
def readFASTOut():
    pass




if __name__=='__main__':
    bStiffening=True
    shapes_twr=[0]
    shapes_bld=[]
    nDOF = 1 + nShapes_twr + nShapes_bld * 3
    q = np.zeros((nDOF,1)) # TODO, full account of q not done
    q[[0]]= 0          # Twr 1
#     q[[1]]=0.1        # Twr 2
#     q[[2]]=0*np.pi/4. # psi

    np.set_printoptions(linewidth=500)
    assembly='auto'
    main_axis='z'
    #StructA= FASTmodel2TNSB('../data/NREL5MW_ED.dat', shapes_twr=shapes_twr,shapes_bld=shapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    StructA= FASTmodel2TNSB('examples/_F0T2RNA/Spar_ED_ForED.dat', shapes_twr=shapes_twr,shapes_bld=shapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    assembly='manual'
#     assembly='auto'
#     main_axis='x'
#     #StructM= FASTmodel2TNSB('../data/NREL5MW_ED.dat', shapes_twr=shapes_twr,shapes_bld=shapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    StructM= FASTmodel2TNSB('examples/_F0T2RNA/Spar_ED_ForED.dat', shapes_twr=shapes_twr,shapes_bld=shapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
#     print('------------------')
    from scipy.linalg import block_diag
#     print('RR')
    RR = np.eye(3)
#     RR = np.zeros((3,3))
#     RR[0,2]=1 # send z to x
#     RR[1,1]=-1 # send y to -y
#     RR[2,0]=1  # send x to z
    RR=block_diag(RR,RR)
#     print(RR)

#     print('Twr: B_T:')
#     print(StructA.Twr.B_inB)
#     print(StructM.Twr.B_inB)
#     print(StructA.Twr.r_O)
#     print(StructM.Twr.r_O)
# 
#     print('Twr.alpha_y:')
#     print(StructA.alpha)
#     print(StructM.alpha)
# 
# 
#     print('Nac: B_N:')
#     print(StructM.Nac.B_inB)
#     print(np.dot(RR, StructA.Nac.B_inB))
#     print('Nac R_B:')
#     print(StructA.Nac.R_b2g)
#     print(StructM.Nac.R_b2g)
# 
#     print('sft: R_S:')
#     print(StructA.sft.R_b2g)
#     print(StructM.sft.R_b2g)
#     print('sft: B_S:')
#     print(StructA.sft.B_inB)
#     print(np.dot(RR,StructM.sft.B_inB))
#     print(np.dot(RR,StructM.sft.BB_inB)-StructA.sft.BB_inB)

#     print('Bld1 R_B:')
#     print(StructA.Blds[0].R_b2g)
#     print(StructM.Blds[0].R_b2g)
#     print('Bld1: B_S:')
#     print(StructA.Blds[0].B_inB)
#     print(np.dot(RR,StructM.Blds[0].B_inB))
#     print(np.dot(RR,StructM.Blds[0].BB_inB)-StructA.Blds[0].BB_inB)
#     print('Bld2: B_S:')
#     print(StructM.Blds[1].B_inB)
#     print(StructA.Blds[1].B_inB)
# #     print(StructM.Blds[1].BB_inB-StructA.Blds[1].BB_inB)
#     print('Bld3: B_S:')
#     print(StructM.Blds[2].B_inB)
#     print(np.dot(RR,StructA.Blds[2].B_inB ))
#     print(np.dot(RR,StructM.Blds[2].BB_inB)-StructA.Blds[2].BB_inB)


#     print('Fields available in `Struct`:')
#     print(Struct.__dict__.keys())
#     print('Twr Damp matrix:')
#     print(StructA.Twr.DD)
#     print(StructM.Twr.DD)
#     print('Twr KK matrix:')
#     print(StructA.Twr.KK)
#     print(StructM.Twr.KK)
    print('Twr Mass matrix:')
    print(StructA.Twr.MM[6:,6:])
#     print(StructM.Twr.MM[6:,6:])
#     print(StructA.Twr.MM[6:,6:]-StructM.Twr.MM[6:,6:])

#     print('Bld Mass matrix:')
#     print(StructM.Blds[0].MM[3:,3:])
#     print(StructA.Blds[0].MM[3:,3:])
#     print('Bld Mass matrix:')
#     print(np.dot(RR.T,StructM.Blds[0].MM).dot(RR))
#     print(StructA.Blds[0].MM-np.dot(RR.T,StructM.Blds[0].MM).dot(RR))
#     print(np.dot(StructA.Blds[0].BB_inB.T,StructA.Blds[0].MM).dot(StructA.Blds[0].BB_inB))
#     print(np.dot(StructM.Blds[0].BB_inB.T,StructM.Blds[0].MM).dot(StructM.Blds[0].BB_inB))
#     print(StructA.Blds[0].MM-StructM.Blds[0].MM)

    print('Damp matrix:')
#     print(StructA.DD)
#     print(StructM.DD)
    print(StructM.DD-StructA.DD)

    print('Stiff matrix:')
#     print(StructA.KK)
#     print(StructM.KK)
    print(StructM.KK-StructA.KK)

    print('Mass matrix:')
    print(StructA.MM)
#     print(StructM.MM)
#     print(StructM.MM-StructA.MM)
#     print(StructA)
#     print(StructM)
#     print('Origin E :',StructM.Grd.r_O.T)
