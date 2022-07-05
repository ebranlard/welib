##
import numpy as np
import copy
import matplotlib.pyplot as plt
import os

from welib.yams.yams import FASTBeamBody, RigidBody
from welib.yams.utils import *
from welib.yams.TNSB import manual_assembly, auto_assembly

import welib.weio as weio

# --------------------------------------------------------------------------------}
# --- Creating a FNSB model from a FAST model
# --------------------------------------------------------------------------------{
# TODO TODO TODO
# TODO TODO TODO HARMONIZE WITH WINDTURBINE.PY AND TNSB..
# TODO TODO TODO
def FASTmodel2FNSB(FST_file,shapes_sub=[0,4], nShapes_bld=0,nSpan_sub=101,nSpan_bld=61,bHubMass=1,bNacMass=1,bBldMass=1,DEBUG=False,main_axis ='x',bStiffening=True, assembly='manual', q=None, bTiltBeforeNac=False,
        spanFrom0=True, # TODO for legacy, we keep this for now..
        bladeMassExpected=None
):
    """ 
    Returns the following structure
      Fnd :  BeamBody
      Shft:  RigiBody
      Nac :  RigidBody
      Blds:  List of BeamBodies

      MM, KK, DD : mass, stiffness and damping matrix of full system

      shapes_sub: 0:ux, 1:uy, 2:uz, 3:vx, 4:vy, 5:vz 
              E.g. for surge and pitch (ux, vy): [0,4]


      NOTE/TODO: compare this with "windturbine.py"
    """
    

    # --- Read fst file
    ext=os.path.splitext(FST_file)[1]
    if ext.lower()!='.fst':
        raise Exception('FNSB requires a fst file as input')

    FST=weio.read(FST_file)
    rootdir = os.path.dirname(FST_file)
    EDfile = os.path.join(rootdir,FST['EDFile'].strip('"')).replace('\\','/')
    subfile = os.path.join(rootdir,FST['SubFile'].strip('"')).replace('\\','/')

    # Reading elastodyn file
    ED      = weio.read(EDfile)
    rootdir = os.path.dirname(EDfile)
    bldfile = os.path.join(rootdir,ED['BldFile(1)'].strip('"')).replace('\\','/')
    twrfile = os.path.join(rootdir,ED['TwrFile'].strip('"')).replace('\\','/')
    #twr     = weio.read(twrfile)
    bld     = weio.read(bldfile)

    # Reading SubDyn file
    sub     = weio.read(subfile)
    nShapes_sub = len(shapes_sub)
    graph = sub.toGraph() # NOTE: this is repeated in bodies.py...
    graph.divideElements(sub['NDiv'])
    graph.sortNodesBy('z')
    df = graph.nodalDataFrame()
    zBot = np.min(df['z'])
    zTop = np.max(df['z'])
    RayleighCoeff=None
    DampMat=None
    if sub['GuyanDampMod']==1:
        # Rayleigh Damping
        RayleighCoeff=sub['RayleighDamp']
        #if RayleighCoeff[0]==0:
        #    damp_zeta=omega*RayleighCoeff[1]/2. 
    elif sub['GuyanDampMod']==2:
        # Full matrix
        DampMat = sub['GuyanDampMatrix']
        DampMat=DampMat[np.ix_(shapes,shapes)]



    nB = ED['NumBl']
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
    IR_hub = np.zeros((3,3))
    I0_nac=np.zeros((3,3)) 

    if main_axis=='x':
        IR_hub[2,2] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2
        I0_nac[0,0]= ED['NacYIner']
    elif main_axis=='z':
        IR_hub[0,0] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2
        I0_nac[2,2] = ED['NacYIner']
    IR_hub = IR_hub * bHubMass
    I0_nac = I0_nac * bNacMass

    # Inertias not at COG...
    IG_hub = translateInertiaMatrix(IR_hub, M_hub, np.array([0,0,0]), r_RGhub_inS)
    IG_nac = translateInertiaMatrixToCOG(I0_nac, M_nac, r_NGnac_inN)

    # --------------------------------------------------------------------------------}
    ## --- Creating bodies
    # --------------------------------------------------------------------------------{
    # Bld
    Blds=[]
    Blds.append(FASTBeamBody('blade',ED,bld,Mtop=0,nShapes=nShapes_bld, nSpan=nSpan_bld, main_axis=main_axis, spanFrom0=spanFrom0, massExpected=bladeMassExpected)) # NOTE: legacy spanfrom0
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

    # ShaftHub Body 
    Sft=RigidBody('ShaftHub',M_hub,IG_hub,r_SGhub_inS);
    # Nacelle Body
    Nac=RigidBody('Nacelle',M_nac,IG_nac,r_NGnac_inN);
    M_rot= sum([B.Mass for B in Blds])
    M_RNA= M_rot + Sft.Mass + Nac.Mass;
    # Substructure Body
    Fnd = FASTBeamBody('substructure', ED, sub, Mtop=M_RNA, shapes=shapes_sub, nSpan=nSpan_sub, main_axis=main_axis, bStiffening=bStiffening)
    #print(Fnd)
    #print('Fnd MM\n',Fnd.MM[6:,6:])
    #print('Fnd KK\n',Fnd.KK[6:,6:])
    # HACK here because doesn't handle this for now
    if sub['GuyanDampMod']==1:
        Fnd.DD[6:,6:] = Fnd.MM[6:,6:]*RayleighCoeff[0] + Fnd.KK[6:,6:]*RayleighCoeff[1] 


    #print('Stiffnening', bStiffening)
    #print('Ttw.KKg   \n', Twr.KKg[6:,6:])
    if DEBUG:
        print('HubMass',Sft.Mass)
        print('NacMass',Nac.Mass)
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
        Struct = manual_assembly(Fnd,Nac,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac)
    else:
        Struct = auto_assembly(Fnd,Nac,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis=main_axis,theta_tilt_y=theta_tilt_y,theta_cone_y=theta_cone_y,DEBUG=DEBUG, bTiltBeforeNac=bTiltBeforeNac)

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

    q_init[iPsi]          = psi_init
    q_init[nDOFMech+iPsi] = omega_init

    Struct.q_init = q_init
    if DEBUG:
        print('Initial conditions:')
        print(q_init)

    # --- Useful data
    Struct.ED=ED

    Struct.DampMat=DampMat
    Struct.RayleighCoeff=RayleighCoeff


    return Struct


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
    #StructA= FASTmodel2TNSB('', nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    StructA= FASTmodel2TNSB('examples/_F0T2RNA/Spar_ED_ForED.dat', nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
    assembly='manual'
#     assembly='auto'
#     main_axis='x'
#     #StructM= FASTmodel2TNSB('', nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, assembly=assembly , q=q, main_axis=main_axis, bStiffening=bStiffening)
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
