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
from welib.yams.yams_rec import YAMSRecFASTBeamBody, YAMSRecRigidBody, YAMSRecGroundBody
from welib.yams.yams_rec import fB_inB, fB_aug, fBMB, fBMatRecursion, fBMatTranslate
from welib.yams.utils import *
from welib.yams.models.TNSB import TNSBStructure

from welib.weio.fast_input_deck import FASTInputDeck










# --------------------------------------------------------------------------------}
# --- Creating a TNSB model automatically 
# --------------------------------------------------------------------------------{
def auto_assembly(twr, yaw, nac, gen, sft, bld, q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis='x',theta_tilt_y=0,theta_yaw=0,theta_cone_y=0,DEBUG=False,bTiltBeforeNac=False, fixedShaft=False, WT=None):
    # TODO gen

    if main_axis=='x':
        #R_NS     = np.dot(R_y(-tilt_up),R_z(q_psi + np.pi)) # << tilt 
        if bTiltBeforeNac:
            R_cn0 = np.dot(R_x (theta_yaw) , R_y (theta_tilt_y))
            R_cs0 = R_z (np.pi)
        else:
            R_cn0 = R_x (theta_yaw) 
            R_cs0 = np.dot( R_y(theta_tilt_y) , R_z (np.pi)) # Note: OrientBefore
        Shaft_axis='z'
    elif main_axis=='z':
        if bTiltBeforeNac:
            R_cn0 = np.dot(R_z (theta_yaw) , R_y(theta_tilt_y))
            R_cs0 = R_x (np.pi)
        else:
            R_cn0 = R_z (theta_yaw) 
            R_cs0 = np.dot(R_y(theta_tilt_y) , R_x (np.pi) )# Note: OrientBefore
        Shaft_axis='x'
    nB=len(bld)


    # Creating reference frame
    grd = YAMSRecGroundBody()

    # Connections between bodies
    grd.connectTo(twr, Point=r_ET_inE, Type='Rigid')
    twr.connectTo(nac, Point=r_TN_inT, Type='Rigid', RelOrientation = R_cn0 , OrientAfter=True)
    twr.connectTo(yaw, Point=r_TN_inT, Type='Rigid', RelOrientation = R_cn0 , OrientAfter=True)
    if fixedShaft:
        nac.connectTo (sft   , Point=r_NS_inN, Type='Rigid', RelOrientation = R_cs0, OrientAfter=False)
    else:
        nac.connectTo (sft   , Point=r_NS_inN, Type='SphericalJoint',JointRotations=[Shaft_axis],RelOrientation = R_cs0, OrientAfter=False)
    for i,B in enumerate(bld):
        psi_B= -i*2*np.pi/nB # 0 -2pi/2 2pi/3  or 0 pi
        if main_axis=='x':
            R_SB = R_z(0*np.pi + psi_B)
        elif main_axis=='z':
            R_SB = R_x(0*np.pi + psi_B)
        R_SB = np.dot(R_SB, R_y(theta_cone_y))
        sft.connectTo(B, Point=r_SR_inS, Type='Rigid', RelOrientation = R_SB)

    # Setting DOF index for all bodies and connections 
    nq=grd.setupDOFIndex();
    if nq!=len(q):
       print('>>> ',nq,len(q))
       raise Exception('Wrong number of dof')

    grd.updateChildrenKinematicsNonRecursive(q)
    twr.updateChildrenKinematicsNonRecursive(q)
    yaw.updateChildrenKinematicsNonRecursive(q)
    nac.updateChildrenKinematicsNonRecursive(q)
    sft.updateChildrenKinematicsNonRecursive(q)


    # --- Full system
    nq = len(q)
    MM = grd.M
    KK = grd.K
    DD = grd.D

    MM[np.abs(MM)< 1e-09] = 0
    # --- returning everthin in a structure class
    if WT is None:
        WT      = TNSBStructure(main_axis=main_axis,theta_cone=theta_cone_y,theta_tilt=theta_tilt_y,bTiltBeforeNac=bTiltBeforeNac)
    else:
        WT.main_axis      = main_axis
        WT.theta_cone     = theta_cone_y
        WT.theta_tilt     = theta_tilt_y
        WT.bTiltBeforeNac = bTiltBeforeNac
    WT.grd  = grd
    WT.twr  = twr
    WT.yaw  = yaw
    WT.nac  = nac
    WT.hubgen  = sft
    WT.bld = bld
    WT.MM   = MM
    WT.KK   = KK
    WT.DD   = DD
    WT.q    = q
    WT.r_ET_inE=r_ET_inE
    WT.r_TN_inT=r_TN_inT
    WT.r_NS_inN=r_NS_inN
    WT.r_SR_inS=r_SR_inS

    WT.init_trigger()

    return WT




# --------------------------------------------------------------------------------}
# --- Manual assembly of a TNSB model 
# --------------------------------------------------------------------------------{
def manual_assembly(twr, yaw, nac, gen, sft, bld, q,r_ET_inE, r_TN_inT, r_NS_inN, r_SR_inS, main_axis='x',theta_tilt_y=0,theta_cone_y=0,DEBUG=False, bTiltBeforeNac=False, fixedShaft=False, WT=None):

    # Main Parameters
    nDOF = len(q)
#     CyT=- np.array([ twr.PhiV[0][2,-1],  1.5065E-01, 0, 0]) # End value of shapes functions in y direction
    CxT=  np.zeros(twr.nf)
    CyT=  np.zeros(twr.nf)
    CzT=  np.zeros(twr.nf)
    UxT=  np.zeros(twr.nf)
    UyT=  np.zeros(twr.nf)
    UzT=  np.zeros(twr.nf)
    for j,(u,v) in enumerate(zip(twr.PhiU,twr.PhiV)):
        if main_axis=='x':
            CyT[j]=-v[2,-1] # A deflection along z gives a negative angle around y
            CzT[j]= v[1,-1] # A deflection along y gives a positive angle around z # TODO TODO CHECK ME
            UyT[j]= u[1,-1] 
            UzT[j]= u[2,-1] 
            #print('Alpha y - mode {}:'.format(j+1),CyT[j])
        elif main_axis=='z':
            CxT[j]=-v[1,-1] # A deflection along y gives a negative angle around x # TODO TODO CHECK ME
            CyT[j]= v[0,-1] # A deflection along x gives a positive angle around y
            UxT[j]= u[0,-1] 
            UyT[j]= u[1,-1] 
    CyT=CyT[:twr.nf]
    UxT=UxT[:twr.nf]
    # TODO:
    #  Bt_pc=zeros(3,p.nf);
    #  for j=1:p.nf
    #      Bx_pc(:,j)=p.PhiU{j}(:,iNode);
    #      Bt_pc(:,j)=[0; -p.PhiV{j}(3,iNode); p.PhiV{j}(2,iNode)];
    #  end

    # --------------------------------------------------------------------------------}
    ## --- "Manual connection"
    # --------------------------------------------------------------------------------{
    # link E-T
    R_ET     = np.identity(3)
    B_T      = np.array([])
    # B_T      = fBMatRecursion(,np.vstack((Bx_ET,Bt_ET)),R_ET,r_ET_inE)
    B_T_inT  = fB_inB(R_ET, B_T)
    BB_T_inT = fB_aug(B_T_inT, twr.nf)
    MM_T     = fBMB(BB_T_inT,twr.MM)
    KK_T     = fBMB(BB_T_inT,twr.KK)
    DD_T     = fBMB(BB_T_inT,twr.DD)

    twr.r_O    = r_ET_inE
    twr.R_b2g   = R_ET
    twr.B      = B_T    
    twr.B_inB  = B_T_inT
    twr.BB_inB = BB_T_inT

    # ---------------------------------------------
    # Link T-N
    # TODO
    if twr.nf == 0:
        Bx_TN = np.array([])
        Bt_TN = np.array([])
        alpha_y=0
    elif twr.nf == 1:
        if main_axis=='x':
            Bx_TN = np.array([[0],[0],[UzT[0]]])
        elif main_axis=='z':
            Bx_TN = np.array([[UxT[0]],[0],[0]])
        Bt_TN = np.array([[0],[CyT[0]],[0]])
        twr.gzf = q[0,0]
        alpha_y = np.dot(CyT.ravel(), q[0,0].ravel())
    elif twr.nf == 2:
        if main_axis=='x':
            Bx_TN = np.array([[0,0],[0,0],[UzT[0],UzT[1]]])
        elif main_axis=='z':
            Bx_TN = np.array([[UxT[0],UxT[1]],[0,0],[0,0]])
        twr.gzf = q[0:2,0]
        Bt_TN = np.array([[0,0],[CyT[0],CyT[1]],[0,0]])
        alpha_y = np.dot(CyT.ravel() , q[:2,0].ravel())
    else:
        # TODO use CzT
        raise NotImplementedError()
    #print('alpha_y',alpha_y)
    R_TN     = R_y(alpha_y)
    if bTiltBeforeNac:
        R_TN     = np.dot(R_TN, R_y(theta_tilt_y))
    R_EN     = np.dot(R_ET, R_TN)
    B_N      = fBMatRecursion(B_T,Bx_TN,Bt_TN,R_ET,r_TN_inT)
    B_N_inN  = fB_inB(R_EN, B_N)
    BB_N_inN = fB_aug(B_N_inN, nac.nf)
    MM_N     = fBMB(BB_N_inN,nac.MM)
    KK_N     = fBMB(BB_N_inN,nac.KK)

    nac.r_O    = twr.r_O + np.dot(twr.R_b2g, r_TN_inT)
    nac.R_b2g   = R_EN
    nac.B      = B_N    
    nac.B_inB  = B_N_inN
    nac.BB_inB = BB_N_inN
    # TODO YAW
    MM_Y       = fBMB(BB_N_inN, yaw.MM)
    yaw.r_O    = twr.r_O + np.dot(twr.R_b2g, r_TN_inT)
    yaw.R_b2g   = R_EN
    yaw.B      = B_N    
    yaw.B_inB  = B_N_inN
    yaw.BB_inB = BB_N_inN
    # TODO gen

    # ---------------------------------------------
    # Link N-S
    if fixedShaft:
        q_psi = 0 # TODO potential azimuth..
    else:
        iPsi = twr.nf # Index of DOF corresponding to azimuth
        q_psi = q[iPsi,0]
    if main_axis=='x':
        R_NS     = R_z(q_psi + np.pi) 
    elif main_axis=='z':
        R_NS     = R_x(q_psi + np.pi) 
    if not bTiltBeforeNac:
        R_NS     = np.dot(R_y(theta_tilt_y),R_NS)
    R_ES     = np.dot(R_EN, R_NS)
    r_NS     = np.dot(R_EN, r_NS_inN)
    if fixedShaft:
        Bx_NS    = np.array([])
        Bt_NS    = np.array([])
    else:
        Bx_NS    = np.array([[0],[0],[0]])
        if main_axis=='x':
            Bt_NS    = np.array([[0],[0],[1]])
        elif main_axis=='z':
            Bt_NS    = np.array([[1],[0],[0]])
    B_S      = fBMatRecursion(B_N,Bx_NS,Bt_NS,R_EN,r_NS)
    B_S_inS  = fB_inB(R_ES, B_S)
    BB_S_inS = fB_aug(B_S_inS, sft.nf)
    MM_S     = fBMB(BB_S_inS,sft.MM)
    KK_S     = fBMB(BB_S_inS,sft.KK)

    sft.r_O    = nac.r_O + r_NS
    sft.R_b2g   = R_ES
    sft.B      = B_S    
    sft.B_inB  = B_S_inS
    sft.BB_inB = BB_S_inS

    # ---------------------------------------------
    # Link S-B1
    nB   = len(bld)
    # Point R
    r_SR  = np.dot(R_ES, r_SR_inS)
    B_R = fBMatRecursion(B_S,[],[],R_ES,r_SR)
    B_R_bis = fBMatTranslate(B_S, r_SR)
    # Points B1, B2, B3
    MM_B      = np.zeros((nDOF,nDOF))
    KK_B      = np.zeros((nDOF,nDOF))
    DD_B      = np.zeros((nDOF,nDOF))
    nf_done=0
    nf_tot = sum([B.nf for B in bld])
    for i,B in enumerate(bld):
        psi_B= -i*2*np.pi/nB # 0 -2pi/2 2pi/3  or 0 pi
        if main_axis=='x':
            R_SB = R_z(0*np.pi + psi_B)
        elif main_axis=='z':
            R_SB = R_x(0*np.pi + psi_B)
        R_SB = np.dot(R_SB, R_y(theta_cone_y))
        R_EB       = np.dot(R_ES, R_SB)
        B_B_inB    = fB_inB(R_EB, B_R)
        BB_B_inB   = fB_aug(B_B_inB, nf_tot, B.nf, nf_done)

        nf_done   += B.nf
        # Full matrices 
        MM_B +=     fBMB(BB_B_inB,B.MM)
        KK_B +=     fBMB(BB_B_inB,B.KK)
        DD_B +=     fBMB(BB_B_inB,B.DD)

        B.r_O    = sft.r_O + r_SR
        B.B      = B_R
        B.B_inB  = B_B_inB
        B.BB_inB = BB_B_inB
     

    # --- Final assembly
    MM = MM_B.copy()
    if fixedShaft:
        iPsi = MM_S.shape[1]-1
    MM[:iPsi+1,:iPsi+1] += MM_S
    MM[:twr.nf,:twr.nf] += MM_T + MM_N + MM_Y
    KK = KK_B
    KK[:iPsi+1,:iPsi+1] += KK_S
    KK[:twr.nf,:twr.nf] += KK_T + KK_N

    DD = DD_B 
    DD[:twr.nf,:twr.nf] += DD_T
    ## Display to screen
    MM[np.abs(MM)< 1e-09] = 0
    if DEBUG:
        print('--------------------- Geom ---------------------')
        print('r_ET_inE   ',r_ET_inE   .T)
        print('r_TN_inT   ',r_TN_inT   .T)
        print('r_NS_inN   ',r_NS_inN   .T)
        print('r_SR_inS   ',r_SR_inS   .T)
        print('-------------------- Tower ---------------------')
        print('CyT\n',CyT)
        print('alpha_y',alpha_y)
        print('B_T\n',B_T)
        print('B_T_inT\n',B_T_inT)
        print('BB_T_inT\n',BB_T_inT)
#         print('MM_T\n',MM_T)
#         print('KK_T\n',KK_T)
#         print('DD_T\n',DD_T)
        print('------------------- Nacelle --------------------')
        print('B_N\n',B_N)
        print('B_N_inN\n',B_N_inN)
        print('BB_N_inN\n',BB_N_inN)
        print('MM_N\n',MM_N)
#         print('-------------------- Shaft ---------------------')
#         print('R_NS\n',R_NS)
#         print('BB_S_inS\n',BB_S_inS)
#         print('MM_S\n',MM_S)
#         print('------------------- Blades ---------------------')
#         #print('BB_B1_inB1')
#         #print(BB_B1_inB1)
#         print('MM_B\n',MM_B)
#         print('KK_B\n',KK_B)
#         print('DD_B\n',DD_B)
#         print('-------------------- Full ----------------------')
#         print('M ("manually" built)')
#         print(MM)
#         print('K ("manually" build)')
#         print(KK)

    ## Eigenvalue analysis
    #[Q,Lambda]=eig(K,M);
    #Omega2=diag(Lambda);
    #[Omega2,Isort]=sort(Omega2);
    #Q=Q(:,Isort);
    #f_eva= sqrt(Omega2)/(2*pi);
    #for i=1:length(f_eva);
    #    fprintf('f%d = %.3f \n',i,f_eva(i))


    # --- returning everthin in a structure class
    if WT is None:
        WT      = TNSBStructure(main_axis=main_axis,theta_cone=theta_cone_y,theta_tilt=theta_tilt_y,bTiltBeforeNac=bTiltBeforeNac)
    else:
        WT.main_axis      = main_axis
        WT.theta_cone     = theta_cone_y
        WT.theta_tilt     = theta_tilt_y
        WT.bTiltBeforeNac = bTiltBeforeNac
    WT.grd = YAMSRecGroundBody()
    WT.twr  = twr
    WT.yaw  = yaw
    WT.nac  = nac
    WT.hubgen  = sft
    WT.bld  = bld
    WT.MM   = MM
    WT.KK   = KK
    WT.DD   = DD
    WT.q    = q
    WT.r_ET_inE=r_ET_inE
    WT.r_TN_inT=r_TN_inT
    WT.r_NS_inN=r_NS_inN
    WT.r_SR_inS=r_SR_inS
    WT.init_trigger()

    return WT


# --------------------------------------------------------------------------------}
# --- Creating a FNSB model from a FAST model
# --------------------------------------------------------------------------------{
# TODO TODO TODO
# TODO TODO TODO HARMONIZE WITH WINDTURBINE.PY AND TNSB..
# TODO TODO TODO
class FASTmodel2FNSB(FASTWindTurbine):

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
        nShapes_bld=len(shapes_bld)

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
        if SD is None:
            raise Exception('Couldnt read SubDyn file')
        if bld is None:
            raise Exception('Couldnt read blade file')
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
            r_EF_inE    = np.array([zBot                         ,0,0]) 
            r_ET_inE    = np.array([ED['TowerBsHt']              ,0,0]) 
            r_FT_inF    = np.array([ED['TowerBsHt']-zBot         ,0,0]) 
            r_TN_inT    = np.array([ED['TowerHt']-ED['TowerBsHt'],0,0])
            if bTiltBeforeNac:
                raise NotImplementedError()
                R_NS0 = np.eye(3)
                R_TN0 = R_y(theta_tilt_y)
            else:
                R_NS0 = R_y(theta_tilt_y)
                R_TN0 = np.eye(3)
                r_NGnac_inN = np.array([ED['NacCMzn']            ,0,ED['NacCMxn']] )
                r_NS_inN    = np.array([ED['Twr2Shft']           ,0,0]) # S on tower axis
            r_SR_inS    = np.array([0                            ,0,ED['OverHang']] ) # S and R 
            r_SGhub_inS = np.array([0                            ,0,ED['OverHang']+ED['HubCM']]   ) # 
        elif main_axis=='z':
            theta_tilt_y=-ED['ShftTilt']*np.pi/180 # NOTE: tilt has wrong orientation in FAST
            theta_cone_y= ED['Precone(1)']*np.pi/180

            r_EF_inE    = np.array([0, 0,zBot  ]) 
            r_ET_inE    = np.array([0, 0,ED['TowerBsHt']        ]) 
            r_FT_inF    = np.array([0, 0,ED['TowerBsHt']-zBot   ]) 
            r_TN_inT    = np.array([0, 0,ED['TowerHt']-ED['TowerBsHt'] ])
            if bTiltBeforeNac:
                raise NotImplementedError()
                R_NS0 = np.eye(3)
                R_TN0 = R_y(theta_tilt_y)
            else:
                R_NS0 = R_y(theta_tilt_y)
                R_TN0 = np.eye(3)
                r_NGnac_inN = np.array([ED['NacCMxn']             ,0,ED['NacCMzn']                 ])
                r_NS_inN    = np.array([0                         ,0,ED['Twr2Shft']                ]) # S on tower axis
            r_SR_inS        = np.array([ED['OverHang']            ,0,0]                             ) # S and R
            r_SGhub_inS     = np.array([ED['OverHang']+ED['HubCM'],0,0]                             ) # 

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
        Blds.append(YAMSRecFASTBeamBody('blade',ED,bld,Mtop=0,shapes=shapes_bld, nSpan=nSpan_bld, main_axis=main_axis, spanFrom0=spanFrom0, massExpected=bladeMassExpected, gravity=gravity, algo=algo)) # NOTE: legacy spanfrom0
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
            Struct= manual_assembly(Fnd,Yaw,Nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac, fixedShaft=fixedShaft, WT=WT)
        else:
            Struct = auto_assembly(Fnd,Yaw,Nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac, fixedShaft=fixedShaft, WT=WT)

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
