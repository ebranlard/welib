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


from welib.yams.windturbine import rigidBlades
from welib.yams.windturbine import FASTWindTurbine
from welib.yams.yams_rec import YAMSRecFASTBeamBody, YAMSRecRigidBody
from welib.yams.utils import *
from welib.yams.models.TNSB import manual_assembly, auto_assembly, TNSBStructure

import welib.weio as weio
from welib.weio.fast_input_file import FASTInputFile
from welib.weio.fast_input_deck import FASTInputDeck

# --------------------------------------------------------------------------------}
# --- Creating a FNSB model from a FAST model
# --------------------------------------------------------------------------------{
# TODO TODO TODO
# TODO TODO TODO HARMONIZE WITH WINDTURBINE.PY AND TNSB..
# TODO TODO TODO
class FASTmodel2FNSB(FASTWindTurbine):

    def __init__(self, FST_file, 
                 shapes_sub=[0,4], nSpan_sub=None,
                 nShapes_bld=0, nSpan_bld=None,
                 bHubMass=1, bNacMass=1, bBldMass=1, 
                 DEBUG=False, 
                 main_axis ='x', bStiffening=True, assembly='manual', q=None, bTiltBeforeNac=False,
                 fixedShaft=False,
                 spanFrom0=True, # TODO for legacy, we keep this for now..
                 bladeMassExpected=None,
                 gravity=None,
                 algo='', # TODO replace with OpenFAST
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

        WT = TNSBStructure()
        FASTWindTurbine.__init__(self, WT=WT)
        # Override the WT

        # --- Read fst file
        ext=os.path.splitext(FST_file)[1]
        if ext.lower()!='.fst':
            raise Exception('FNSB requires a fst file as input')

        DCK = FASTInputDeck(FST_file)
        FST = DCK.fst_vt['Fst']
        ED  = DCK.fst_vt['ElastoDyn']
        SD  = DCK.fst_vt['SubDyn']
        bld  = DCK.fst_vt['ElastoDynBlade']
        if gravity is None:
            try:
                gravity = FST['gravity']
            except:
                gravity = ED['gravity'] # Old interface

        # --- Reading SubDyn file
        nShapes_sub = len(shapes_sub)
        graph = SD.toGraph() # NOTE: this is repeated in bodies.py...
        graph.divideElements(SD['NDiv'])
        graph.sortNodesBy('z')
        df = graph.nodalDataFrame()
        zBot = np.min(df['z'])
        zTop = np.max(df['z'])
        RayleighCoeff=None
        DampMat=None
        if SD['GuyanDampMod']==1:
            # Rayleigh Damping
            RayleighCoeff=SD['RayleighDamp']
            #if RayleighCoeff[0]==0:
            #    damp_zeta=omega*RayleighCoeff[1]/2. 
        elif SD['GuyanDampMod']==2:
            # Full matrix
            DampMat = SD['GuyanDampMatrix']
            DampMat=DampMat[np.ix_(shapes,shapes)]

        # --- Default arguments
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

        nB = ED['NumBl']
        if fixedShaft:
            nDOF = nShapes_sub + nShapes_bld * nB # 
        else:
            nDOF = 1 + nShapes_sub + nShapes_bld * nB # +1 for Shaft
        if q is None:
            q = np.zeros((nDOF,1)) # TODO, full account of q not done

        ## --- Strucural and geometrical Inputs
        if main_axis=='x':
            theta_tilt_y= ED['ShftTilt']*np.pi/180 # NOTE: tilt has wrong orientation in FAST
            theta_cone_y=-ED['Precone(1)']*np.pi/180
            r_EF_inE    = np.array([[zBot]                         ,[0],[0]]) 
            r_ET_inE    = np.array([[ED['TowerBsHt']]              ,[0],[0]]) 
            r_FT_inF    = np.array([[ED['TowerBsHt']-zBot]         ,[0],[0]]) 
            r_TN_inT    = np.array([[ED['TowerHt']-ED['TowerBsHt']],[0],[0]])
            if bTiltBeforeNac:
                raise NotImplementedError()
                R_NS0 = np.eye(3)
                R_TN0 = R_y(theta_tilt_y)
            else:
                R_NS0 = R_y(theta_tilt_y)
                R_TN0 = np.eye(3)
                r_NGnac_inN = np.array([[ED['NacCMzn']]                ,[0],[ED['NacCMxn']]] )
                r_NS_inN    = np.array([[ED['Twr2Shft']]               ,[0],[0]]) # S on tower axis
            r_SR_inS    = np.array([[0]                            ,[0],[ED['OverHang']]] ) # S and R 
            r_SGhub_inS = np.array([[0]                            ,[0],[ED['OverHang']+ED['HubCM']]]   ) # 
        elif main_axis=='z':
            theta_tilt_y=-ED['ShftTilt']*np.pi/180 # NOTE: tilt has wrong orientation in FAST
            theta_cone_y= ED['Precone(1)']*np.pi/180

            r_EF_inE    = np.array([[0], [0] , [zBot]  ]) 
            r_ET_inE    = np.array([[0], [0], [ED['TowerBsHt']]        ]) 
            r_FT_inF    = np.array([[0], [0], [ED['TowerBsHt']-zBot]   ]) 
            r_TN_inT    = np.array([[0], [0],[ED['TowerHt']-ED['TowerBsHt']] ])
            if bTiltBeforeNac:
                raise NotImplementedError()
                R_NS0 = np.eye(3)
                R_TN0 = R_y(theta_tilt_y)
            else:
                R_NS0 = R_y(theta_tilt_y)
                R_TN0 = np.eye(3)
                r_NGnac_inN = np.array([[ED['NacCMxn']]             ,[0],[ED['NacCMzn']]                 ])
                r_NS_inN    = np.array([[0]                         ,[0],[ED['Twr2Shft']]                ]) # S on tower axis
            r_SR_inS    = np.array([[ED['OverHang']]            ,[0],[0]]                             ) # S and R
            r_SGhub_inS = np.array([[ED['OverHang']+ED['HubCM']],[0],[0]]                             ) # 

        r_RGhub_inS = - r_SR_inS + r_SGhub_inS


        M_hub   = ED['HubMass']*bHubMass
        M_nac   = ED['NacMass'] *bNacMass
        M_yaw   = ED['YawBrMass']
        IR_hub = np.zeros((3,3))
        I0_nac=np.zeros((3,3)) 

            # TODO, here hub and Gen put together...
        if main_axis=='x':
            IR_hub[2,2] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2
            I0_nac[0,0]= ED['NacYIner']
        elif main_axis=='z':
            IR_hub[0,0] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2
            I0_nac[2,2] = ED['NacYIner']
        IR_hub = IR_hub * bHubMass
        I0_nac = I0_nac * bNacMass

        # Inertias not at COG...
        IG_hub = translateInertiaMatrix(I_A=IR_hub, Mass=M_hub, r_BG=np.array([0,0,0]), r_AG=r_RGhub_inS)
        IG_nac = translateInertiaMatrixToCOG(I0_nac, M_nac, r_NGnac_inN)

        # --------------------------------------------------------------------------------}
        ## --- Creating bodies
        # --------------------------------------------------------------------------------{
        # Bld
        Blds=[]
        Blds.append(YAMSRecFASTBeamBody('blade',ED,bld,Mtop=0,nShapes=nShapes_bld, nSpan=nSpan_bld, main_axis=main_axis, spanFrom0=spanFrom0, massExpected=bladeMassExpected, gravity=gravity, algo=algo)) # NOTE: legacy spanfrom0
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
        # Nacelle Body
        Nac=YAMSRecRigidBody('Nacelle',M_nac,IG_nac,r_NGnac_inN);
        # Yaw Bearing # TODO TODO TODO
        Yaw=YAMSRecRigidBody('YawBearing',M_yaw,(0,0,0),(0,0,0));
        if M_yaw>0:
            print('[WARN] TODO YAW BEARING MASS NOT FULLY IMPLEMENTED IN FTNSB')

        M_rot= sum([B.mass for B in Blds])
        M_RNA= M_rot + Sft.mass + Nac.mass + Yaw.mass
        # Tower Body
        #   None for now
        # Substructure Body
        Fnd = YAMSRecFASTBeamBody('substructure', ED, SD, Mtop=M_RNA, shapes=shapes_sub, nSpan=nSpan_sub, main_axis=main_axis, bStiffening=bStiffening, gravity=gravity, algo=algo)
        #print(Fnd)
        #print('Fnd MM\n',Fnd.MM[6:,6:])
        #print('Fnd KK\n',Fnd.KK[6:,6:])
        # HACK here because doesn't handle this for now
        if SD['GuyanDampMod']==1:
            Fnd.DD[6:,6:] = Fnd.MM[6:,6:]*RayleighCoeff[0] + Fnd.KK[6:,6:]*RayleighCoeff[1] 


        #print('Stiffnening', bStiffening)
        #print('Ttw.KKg   \n', Twr.KKg[6:,6:])
        if DEBUG:
            print('HubMass',Sft.mass)
            print('NacMass',Nac.mass)
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
        r_ET_inE = r_EF_inE
        r_TN_inT = r_FT_inF+r_TN_inT # assume that F and T are in system E here

        if assembly=='manual':
             manual_assembly(Fnd,Yaw,Nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac, fixedShaft=fixedShaft, WT=WT)
        else:
            auto_assembly(Fnd,Yaw,Nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac, fixedShaft=fixedShaft, WT=WT)

        # --- Initial conditions
        omega_init = ED['RotSpeed']*2*np.pi/60 # rad/s
        psi_init   = ED['Azimuth']*np.pi/180   # rad
        FA_init    = ED['TTDspFA']
        Surge_init = ED['PtfmSurge']
        Sway_init  = ED['PtfmSway']
        Heave_init = ED['PtfmHeave']
        Roll_init  = ED['PtfmRoll'] *np.pi/180
        Pitch_init = ED['PtfmPitch']*np.pi/180
        Yaw_init   = ED['PtfmYaw']*np.pi/180

        iPsi     = Struct.iPsi
        nDOFMech = len(Struct.MM)
        q_init   = np.zeros(2*nDOFMech) # x2, state space

        if nShapes_sub>0:
            sub_init = np.array([Surge_init, Sway_init, Heave_init, Roll_init, Pitch_init, Yaw_init])
            for iDOF,iDOFfull in enumerate(shapes_sub):
                q_init[iDOF] = sub_init[iDOFfull]

        if not fixedShaft:
            q_init[iPsi]          = psi_init
            q_init[nDOFMech+iPsi] = omega_init

        Struct.q_init = q_init
        if DEBUG:
            print('Initial conditions:')
            print(q_init)

        # --- Useful data
        WT=self.WT
        WT.DCK = DCK
        WT.FST = FST
        WT.ED  = ED
        WT.SD  = SD
        WT.Hydro     = FST['CompHydro']>0 # FST['CompSeaSt']>0 and 
        WT.HD        = DCK.fst_vt['HydroDyn']
        try:
            WT.WtrDens   = FST['WtrDens']
            WT.WtrDpth   = FST['WtrDpth']
        except:
            WT.WtrDens   = Struct.HD['WtrDens']
            WT.WtrDpth   = Struct.HD['WtrDpth']

        WT.DampMat=DampMat
        WT.RayleighCoeff=RayleighCoeff
        WT.additional_properties +=['DCK', 'FST', 'ED', 'DampMat', 'RayleighCoeff', 'WaterDepth','Hydro']

# --------------------------------------------------------------------------------}
# --- Read Relevant fields from an outb file 
# --------------------------------------------------------------------------------{
def readFASTOut():
    pass




if __name__=='__main__':
    FASTmodel2FNSB()
    pass
