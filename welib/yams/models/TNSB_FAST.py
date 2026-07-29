##
import numpy as np
import copy
import matplotlib.pyplot as plt
import os

from welib.yams.windturbine import FASTWindTurbine
from welib.yams.yams import FASTBeamBody, YAMSRecRigidBody
from welib.yams.utils import *
from welib.yams.models.TNSB import manual_assembly, auto_assembly, TNSBStructure

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
    
    def __init__(self, FST_file,nB=3,nShapes_twr=2, nShapes_bld=0, 
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

        WT = TNSBStructure()
        FASTWindTurbine.__init__(self, WT=WT,
                                 main_axis=main_axis, 
                                 algo=algo)


        # --- Read fast input files
        readlist = ['Fst', 'ED', 'EDtwr', 'EDbld']
        self.loadFST(FST_file, readlist=readlist)

        
        nDOF = 1 + nShapes_twr + nShapes_bld * nB # +1 for Shaft
        if q is None:
            q = np.zeros((nDOF,1)) # TODO, full account of q not done

        # --- LEGACY
        ED = self.ED

        # --- Default arguments
        if nSpan_twr is None:
            if algo=='OpenFAST':
                nSpan_twr = ED['TwrNodes']
                if verbose:
                    print('[INFO] TNSB_FAST: Using number of tower nodes ({}) from OpenFAST Input file.'.format(nSpan_twr))
            else:
                nSpan_twr=101
                if verbose:
                    print('[INFO] TNSB_FAST: Using default of tower nodes ({}).'.format(nSpan_twr))
        else:
            if algo=='OpenFAST':
                if verbose:
                    print('[INFO] TNSB_FAST: Using user-specified number of tower nodes ({}).'.format(nSpan_twr))
        if nSpan_bld is None:
            if algo=='OpenFAST':
                nSpan_bld = ED['BldNodes']
                if verbose:
                    print('[INFO] TNSB_FAST: Using number of blade nodes ({}) from OpenFAST Input file.'.format(nSpan_bld))
            else:
                nSpan_bld=61
                if verbose:
                    print('[INFO] TNSB_FAST: Using default number of blade nodes ({}).'.format(nSpan_bld))
        else:
            if algo=='OpenFAST':
                if verbose:
                    print('[INFO] TNSB_FAST: Using user-specified number of blade nodes ({}).'.format(nSpan_bld))


        ## --- Strucural and geometrical Inputs
        self.setupEDGeom(zBot=0, bTiltBeforeNac=bTiltBeforeNac)
        self.reshape_3array_to_atleast_2d()
        # --- Legacy
        theta_tilt_y =  self.WT.shaft_tilt
        theta_cone_y =  self.WT.theta_cone_y
        r_ET_inE     =  self.WT.r_ET_inE
        r_TN_inT     =  self.WT.r_TN_inT    
        R_NS0        =  self.WT.R_NS0 
        R_TN0        =  self.WT.R_TN0 
        r_NGnac_inN  =  self.WT.r_NGnac_inN
        r_NS_inN     =  self.WT.r_NS_inN   
        r_SR_inS     =  self.WT.r_SR_inS    
        r_SGhub_inS  =  self.WT.r_SGhub_inS 
        r_RGhub_inS  =  self.WT.r_RGhub_inS 

        # --- Hub
        # TODO, here hub and Gen put together...
        M_hub   = ED['HubMass']*bHubMass
        IR_hub = np.zeros((3,3))
        if main_axis=='x':
            IR_hub[2,2] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2
        elif main_axis=='z':
            IR_hub[0,0] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2
        IR_hub = IR_hub * bHubMass
        IG_hub = translateInertiaMatrix(I_A=IR_hub, Mass=M_hub, r_BG=np.array([0,0,0]), r_AG=r_RGhub_inS)

        # --- Nac
        self.setupEDNac(bNacMass=bNacMass, flavor='yams_rec')
        nac = self.WT.nac

        # --- Yaw
        M_yaw   = ED['YawBrMass']
        # Yaw Bearing # TODO TODO TODO
        Yaw=YAMSRecRigidBody('YawBearing',M_yaw,(0,0,0),(0,0,0));
        if M_yaw>0:
            print('[WARN] TODO YAW BEARING MASS NOT FULLY IMPLEMENTED IN TNSB')



        # --------------------------------------------------------------------------------}
        ## --- Creating bodies
        # --------------------------------------------------------------------------------{
        # Bld
        Blds=[]
        Blds.append(FASTBeamBody('blade',ED,self.bldFile,Mtop=0,nShapes=nShapes_bld, nSpan=nSpan_bld, main_axis=main_axis, spanFrom0=spanFrom0, massExpected=bladeMassExpected, gravity=gravity, algo=algo)) # NOTE: legacy spanfrom0
        Blds[0].MM *=bBldMass
        for iB in range(nB-1):
            Blds.append(copy.deepcopy(Blds[0]))
        # IMPORTANT FOR RNA set R_b2g
        for iB,B in enumerate(Blds):
            B.name='bld'+str(iB+1)
            psi_B= -iB*2*np.pi/len(Blds) 
            if main_axis=='x':
                R_SB = R_z(0*np.pi + psi_B) # TODO psi offset and psi0
            elif main_axis=='z':
                R_SB = R_x(0*np.pi + psi_B) # TODO psi0
            R_SB = np.dot(R_SB, R_y(ED['PreCone({})'.format(iB+1)]*np.pi/180)) # blade2shaft
            B.R_b2g= R_SB

        # ShaftHubGen Body  NOTE: generator!!! This is ugly
        Sft=YAMSRecRigidBody('ShaftHubGen',M_hub,IG_hub,r_SGhub_inS)
        
        # Gen only
        Gen=YAMSRecRigidBody('Gen', 0, IG_hub, r_SGhub_inS)

        #print('>>> IG_hub',IG_hub, r_SGhub_inS)

        M_rot= sum([B.mass for B in Blds])
        M_RNA= M_rot + Sft.mass + self.WT.nac.mass + Yaw.mass
        # Tower Body
        Twr = FASTBeamBody('tower',ED,self.twrFile,Mtop=M_RNA,nShapes=nShapes_twr, nSpan=nSpan_twr, main_axis=main_axis,bStiffening=bStiffening, gravity=gravity, algo=algo)
        #print('Stiffnening', bStiffening)
        #print('Ttw.KKg   \n', Twr.KKg[6:,6:])
        if DEBUG:
            print('HubMass',Sft.mass)
            print('NacMass',nac.mass)
            print('RotMass',M_rot)
            print('RNAMass',M_RNA)
            print('IG_hub')
            print(IG_hub)
            print('IG_nac')
            print(IG_nac)
            print('I_gen_LSS', ED['GenIner']*ED['GBRatio']**2)
            print('I_hub_LSS', ED['hubIner'])
            print('I_rot_LSS', nB*Blds[0].MM[5,5])
            print('I_tot_LSS', nB*Blds[0].MM[5,5]+ED['hubIner']+ED['GenIner']*ED['GBRatio']**2) 
            print('r_NGnac_inN',r_NGnac_inN.T)
            print('r_SGhub_inS',r_SGhub_inS.T)
        # --------------------------------------------------------------------------------}
        # --- Assembly 
        # --------------------------------------------------------------------------------{
        if assembly=='manual':
             manual_assembly(Twr,Yaw,nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac, WT=self.WT)
        else:
            auto_assembly(Twr,Yaw,nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac, WT=self.WT)

        # --- Initial conditions
        omega_init = ED['RotSpeed']*2*np.pi/60 # rad/s
        psi_init   = ED['Azimuth']*np.pi/180   # rad
        FA_init    = ED['TTDspFA']
        iPsi     = self.WT.iPsi
        nDOFMech = len(self.WT.MM)
        q_init   = np.zeros(2*nDOFMech) # x2, state space

        if nShapes_twr>0:
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
    nShapes_twr=1
    nShapes_bld=0
    nDOF = 1 + nShapes_twr + nShapes_bld * 3
    q = np.zeros((nDOF,1)) # TODO, full account of q not done
    q[[0]]= 0          # Twr 1
#     q[[1]]=0.1        # Twr 2
#     q[[2]]=0*np.pi/4. # psi

    np.set_printoptions(linewidth=500)
    assembly='auto'
    main_axis='z'
    #StructA= FASTmodel2TNSB('../data/NREL5MW_ED.dat', nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    StructA= FASTmodel2TNSB('examples/_F0T2RNA/Spar_ED_ForED.dat', nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    assembly='manual'
#     assembly='auto'
#     main_axis='x'
#     #StructM= FASTmodel2TNSB('../data/NREL5MW_ED.dat', nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    StructM= FASTmodel2TNSB('examples/_F0T2RNA/Spar_ED_ForED.dat', nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
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
#     print(StructA.Nac.R_0b)
#     print(StructM.Nac.R_0b)
# 
#     print('Sft: R_S:')
#     print(StructA.Sft.R_0b)
#     print(StructM.Sft.R_0b)
#     print('Sft: B_S:')
#     print(StructA.Sft.B_inB)
#     print(np.dot(RR,StructM.Sft.B_inB))
#     print(np.dot(RR,StructM.Sft.BB_inB)-StructA.Sft.BB_inB)

#     print('Bld1 R_B:')
#     print(StructA.Blds[0].R_0b)
#     print(StructM.Blds[0].R_0b)
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
