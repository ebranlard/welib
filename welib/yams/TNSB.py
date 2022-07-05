"""
TNSB refer to : Tower Nacelle, Shaft, Blades

This scripts provides some helper functions to simulates such an assembly of bodies using the Rayleigh-Rizt approximation and joint coordinates. 

The theory is provided in the reference below. The article also contains an example in the its section, which is reproduced in the file test_TNSB.py.


Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

##
import numpy as np
import copy
import os

try:
    from .yams import *
except:
    from yams import *

from welib.yams.windturbine import rigidBlades

class Structure():
    def __init__(self,main_axis='x',theta_tilt=0,theta_yaw=0,theta_cone=0,bTiltBeforeNac=False):
        self.main_axis      = main_axis
        self.theta_tilt     = theta_tilt
        self.theta_cone     = theta_cone
        self.theta_yaw      = theta_yaw
        self.bTiltBeforeNac = bTiltBeforeNac

    def compute_RNA(s):
        s.M_rot= sum([B.mass for B in s.Blds])
        s.M_RNA= s.M_rot + s.Sft.mass + s.Nac.mass;
        s.r_NGnac_inN = s.Nac.s_G_inB
        s.r_NGhub_inN = s.r_NS_inN + np.dot(s.Nac.R_0b.T, np.dot(s.Sft.R_0b, s.Sft.s_G_inB))

        try:
            # --------------------------------------------------------------------------------}
            # --- New method, requires blade to have "toRigidBody"
            # --------------------------------------------------------------------------------{
            # --- Rigid Blades (with origin R, using N as "global" ref)
            blds_rigid = rigidBlades(s.Blds, r_O = [0,0,0]) # TODO blade origins might be wrong in TNSB
            blds_rigid.pos_global = s.r_NR_inN.ravel()
            R_NS = R_y(s.theta_tilt)  # Rotation fromShaft to Nacelle
            blds_rigid.R_b2g      = R_NS


            # --- Creating "hub" with respect to point N (Sft is wrt S)
            #M_hub  = ED['HubMass']
            #JxxHub_atR = ED['HubIner']
            #hub = RigidBody('Hub', M_hub, (JxxHub_atR,0,0), s_OG=r_SGhub_inS, R_b2g=R_NS, s_OP=r_SR_inS, r_O=r_NS_inN) 
            r_NS_inN=s.Sft.r_O.ravel()-s.Nac.r_O.ravel()
            r_SGhub_inS=s.Sft.masscenter # In body coordinates, S, titled
            #hub = RigidBody('Hub', s.Sft.mass, J=s.Sft.masscenter_inertia, s_OG=r_SGhub_inS, R_b2g=R_NS, r_O=r_NS_inN) 

            hub = RigidBody('Hub', s.Sft.mass, J_G=s.Sft.masscenter_inertia, rho_G=s.Sft.masscenter)
            #hub.shiftOrigin(R_NS.T.dot(-r_NS_inN))
            hub.R_b2g      = R_NS
            hub.pos_global = r_NS_inN



            # --- Rotor = Hub + Blades (with origin R, using N as global ref)
            #rot = blds_rigid.combine(s.Sft, R_b2g=R_NS, r_O=blds_rigid.pos_global)
            rot = blds_rigid.combine(hub, R_b2g=R_NS, r_O=blds_rigid.pos_global)
            rot.name='rotor'
            #rotgen = rot.combine(gen, R_b2g=R_NS, r_O=blades.pos_global)

            #RNA = rot.combine(gen).combine(nac,r_O=[0,0,0])
            RNA = rot.combine(s.Nac,r_O=[0,0,0]).combine(s.Yaw, r_O=[0,0,0])
            s.RNA=RNA
            s.r_NGrna_inN=RNA.masscenter # in N
            s.r_NGrot_inN = rot.masscenter

            # Temp storage
            s.blds_rigid = blds_rigid
            s.rot        = rot
            s.hub        = hub


        except:
            print('[WARN] TNSB: Fail to compute RNA with new method, using legacy')
            s.r_NGrot_inN = s.r_NR_inN   # NOTE approximation neglecting cone, putting all rotor mass at R
            s.r_NGrna_inN = 1./s.M_RNA * (s.Nac.mass*s.r_NGnac_inN + s.Sft.mass*s.r_NGhub_inN +  s.M_rot*s.r_NGrot_inN)



    def init_trigger(s):
        s.alpha = s.Twr.alpha_couplings
        s.iPsi  = s.Twr.nf # Index of DOF corresponding to azimuth

        # Useful for load computation
        s.r_NR_inN = s.r_NS_inN + np.dot(s.Nac.R_0b.T, np.dot(s.Sft.R_0b, s.r_SR_inS))
        s.gravity  = s.Twr.gravity
        s.compute_RNA()
        s.nDOF = len(s.q)
        s.nShapes_twr = s.Twr.nf
        s.nShapes_bld = s.Blds[0].nf


    def GF(s,T,x,alpha_y_fact=1):
        """ 
        T is the force along the shaft
        """
        if (s.nShapes_twr!=1):
            raise NotImplementedError('Number of shape function not 1')
        if s.main_axis=='x':
            raise NotImplementedError('Main axis along x')

        # update tower kinematics
        s.Twr.gzf=x 
        alpha = s.Twr.alpha_couplings
        alpha_y=alpha[1]*alpha_y_fact

        rhoN_x = s.r_NGrna_inN[0,0]
        rhoN_z = s.r_NGrna_inN[2,0]
        rNR_x  = s.r_NR_inN[0,0]
        rNR_z  = s.r_NR_inN[2,0]
        g      = s.gravity
        ux1c   = 1
        vy1c   = s.Twr.Bhat_t_bc[1,0]  # Bhat_t_bc[1,j]= self.PhiV[j][0,iNode]
        ky1c   = s.Twr.PhiK[0][0,-1]

        Fz_inE =-T*sin(alpha_y + s.theta_tilt) - s.M_RNA*g # TODO potential softening correction

        Fx_inE = T*cos(alpha_y + s.theta_tilt)
#         Fx_inE = T*cos(s.theta_tilt)


        # --- Softening is already in K
        #         k  = x[0] * s.Twr.PhiK[0][0,:]
        #         vL = x[0] * vy1c
        #         U  =        s.Twr.PhiU[0][0,:]
        #         GF_soft=0
        #         GF_soft+= np.trapz(+U*k*Fz_inE   , s.Twr.s_span)
        #         GF_soft+=          -vL*Fz_inE*ux1c

        My_inE = 0
        My_inE += s.M_RNA*g*( rhoN_x*cos(alpha_y) + rhoN_z*sin(alpha_y))
        My_inE +=T*(rNR_x*sin(s.theta_tilt) + rNR_z*cos(s.theta_tilt) )


        GF =0
        GF += Fx_inE        * ux1c
        GF += vy1c* My_inE
        # GF = GF_soft  

        return GF

    def GF_lin(s,T,x,bFull=True):
        """ 
        T is the force along the shaft
        Fisrt linearization: assumes the sum of alpha_y small
        """
        if (s.nShapes_twr!=1):
            raise NotImplementedError('Number of shape function not 1')
        if s.main_axis=='x':
            raise NotImplementedError('Main axis along x')

        rhoN_x = s.r_NGrna_inN[0,0]
        rhoN_z = s.r_NGrna_inN[2,0]
        rNR_x  = s.r_NR_inN[0,0]
        rNR_z  = s.r_NR_inN[2,0]
        g      = s.gravity
        ux1c   = s.Twr.Bhat_x_bc[1,0]
        vy1c   = s.Twr.Bhat_t_bc[1,0]  # Bhat_t_bc[1,j]= self.PhiV[j][0,iNode]

        GF  =   T*cos(s.theta_tilt) 
        if bFull:
            GF += - T* vy1c * sin(s.theta_tilt) * x[0]
            GF += (vy1c**2 * s.M_RNA*g * rhoN_z) * x[0]
            GF +=  T*vy1c*(rNR_x*sin(s.theta_tilt) + rNR_z*cos(s.theta_tilt) ) 
            GF += vy1c * s.M_RNA*g * rhoN_x
        return GF



    def print_info(s):
        print('----------------------------------------------------------------')
        print('main_axis :', s.main_axis)
        print('gravity   :', s.gravity)
        print('tilt      :', s.theta_tilt*180/np.pi)
        print('cone      :', s.theta_cone*180/np.pi)
        print('yaw       :', s.theta_yaw *180/np.pi)

    def print_origins(s):
        print('Origin T :',s.Twr.r_O.T)
        print('Origin N :',s.Nac.r_O.T)
        print('Origin R :',s.Blds[0].r_O.T)
        print('Origin S :',s.Sft.r_O.T)

    def print_RNA(s):
        print('----------------- RNA ---------------------------------------')
        print('M_RNA      ', s.M_RNA)
        print('r_NGrna_inN',s.r_NGrna_inN.T)
        print('     r_NGnac_inN ',s.r_NGnac_inN.T , 'M_nac',s.Nac.mass)
        print('     r_NGhub_inN ',s.r_NGhub_inN.T , 'M_hub',s.Sft.mass)
        print('     r_NGrot_inN ',s.r_NGrot_inN.T , 'M_rot',s.M_rot)

    def print_couplings(s):
        print('---------------Couplings ---------------------------------------')
        print('Constant: (Bhat_t)')
        print(s.Twr.Bhat_t_bc) # Constant
        print('Time varying:')
        print(s.Twr.alpha_couplings) # Time varying function of Twr.gzf

    def __repr__(self):
        self.print_info()
        self.print_origins()
        self.print_RNA()
        self.print_couplings()
        return ''








# --------------------------------------------------------------------------------}
# --- Creating a TNSB model automatically 
# --------------------------------------------------------------------------------{
def auto_assembly(Twr,Yaw,Nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis='x',theta_tilt_y=0,theta_yaw=0,theta_cone_y=0,DEBUG=False,bTiltBeforeNac=False):
    # TODO Gen

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
    nB=len(Blds)


    # Creating reference frame
    Grd = GroundBody()

    # Connections between bodies
    Grd.connectTo(Twr, Point=r_ET_inE, Type='Rigid')
    Twr.connectTo(Nac, Point=r_TN_inT, Type='Rigid', RelOrientation = R_cn0 , OrientAfter=True)
    Twr.connectTo(Yaw, Point=r_TN_inT, Type='Rigid', RelOrientation = R_cn0 , OrientAfter=True)
    Nac.connectTo (Sft   , Point=r_NS_inN, Type='SphericalJoint',JointRotations=[Shaft_axis],RelOrientation = R_cs0, OrientAfter=False)
    for i,B in enumerate(Blds):
        psi_B= -i*2*np.pi/nB # 0 -2pi/2 2pi/3  or 0 pi
        if main_axis=='x':
            R_SB = R_z(0*np.pi + psi_B)
        elif main_axis=='z':
            R_SB = R_x(0*np.pi + psi_B)
        R_SB = np.dot(R_SB, R_y(theta_cone_y))
        Sft.connectTo(B, Point=r_SR_inS, Type='Rigid', RelOrientation = R_SB)

    # Setting DOF index for all bodies and connections 
    nq=Grd.setupDOFIndex();
    if nq!=len(q):
       print('>>> ',nq,len(q))
       raise Exception('Wrong number of dof')

    Grd.updateChildrenKinematicsNonRecursive(q)
    Twr.updateChildrenKinematicsNonRecursive(q)
    Yaw.updateChildrenKinematicsNonRecursive(q)
    Nac.updateChildrenKinematicsNonRecursive(q)
    Sft.updateChildrenKinematicsNonRecursive(q)


    # --- Full system
    nq = len(q)
    MM = Grd.M
    KK = Grd.K
    DD = Grd.D

    MM[np.abs(MM)< 1e-09] = 0
    # --- returning everthin in a structure class
    Struct      = Structure(main_axis=main_axis,theta_cone=theta_cone_y,theta_tilt=theta_tilt_y,bTiltBeforeNac=bTiltBeforeNac)
    Struct.Grd  = Grd
    Struct.Twr  = Twr
    Struct.Yaw  = Yaw
    Struct.Nac  = Nac
    Struct.Sft  = Sft
    Struct.Blds = Blds
    Struct.MM   = MM
    Struct.KK   = KK
    Struct.DD   = DD
    Struct.q    = q
    Struct.r_ET_inE=r_ET_inE
    Struct.r_TN_inT=r_TN_inT
    Struct.r_NS_inN=r_NS_inN
    Struct.r_SR_inS=r_SR_inS

    Struct.init_trigger()

    return Struct




# --------------------------------------------------------------------------------}
# --- Manual assembly of a TNSB model 
# --------------------------------------------------------------------------------{
def manual_assembly(Twr,Yaw,Nac,Gen,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis='x',theta_tilt_y=0,theta_cone_y=0,DEBUG=False, bTiltBeforeNac=False):

    # Main Parameters
    nDOF = len(q)
    iPsi = Twr.nf # Index of DOF corresponding to azimuth
#     CyT=- np.array([ Twr.PhiV[0][2,-1],  1.5065E-01, 0, 0]) # End value of shapes functions in y direction
    CxT=  np.zeros(Twr.nf)
    CyT=  np.zeros(Twr.nf)
    CzT=  np.zeros(Twr.nf)
    UxT=  np.zeros(Twr.nf)
    UyT=  np.zeros(Twr.nf)
    UzT=  np.zeros(Twr.nf)
    for j,(u,v) in enumerate(zip(Twr.PhiU,Twr.PhiV)):
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
    CyT=CyT[:Twr.nf]
    UxT=UxT[:Twr.nf]
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
    BB_T_inT = fB_aug(B_T_inT, Twr.nf)
    MM_T     = fBMB(BB_T_inT,Twr.MM)
    KK_T     = fBMB(BB_T_inT,Twr.KK)
    DD_T     = fBMB(BB_T_inT,Twr.DD)

    Twr.r_O    = r_ET_inE
    Twr.R_0b   = R_ET
    Twr.B      = B_T    
    Twr.B_inB  = B_T_inT
    Twr.BB_inB = BB_T_inT

    # ---------------------------------------------
    # Link T-N
    # TODO
    if Twr.nf == 0:
        Bx_TN = np.array([])
        Bt_TN = np.array([])
        alpha_y=0
    elif Twr.nf == 1:
        if main_axis=='x':
            Bx_TN = np.array([[0],[0],[UzT[0]]])
        elif main_axis=='z':
            Bx_TN = np.array([[UxT[0]],[0],[0]])
        Bt_TN = np.array([[0],[CyT[0]],[0]])
        Twr.gzf = q[0,0]
        alpha_y = np.dot(CyT.ravel(), q[0,0].ravel())
    elif Twr.nf == 2:
        if main_axis=='x':
            Bx_TN = np.array([[0,0],[0,0],[UzT[0],UzT[1]]])
        elif main_axis=='z':
            Bx_TN = np.array([[UxT[0],UxT[1]],[0,0],[0,0]])
        Twr.gzf = q[0:2,0]
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
    BB_N_inN = fB_aug(B_N_inN, Nac.nf)
    MM_N     = fBMB(BB_N_inN,Nac.MM)
    KK_N     = fBMB(BB_N_inN,Nac.KK)

    Nac.r_O    = Twr.r_O + np.dot(Twr.R_0b, r_TN_inT)
    Nac.R_0b   = R_EN
    Nac.B      = B_N    
    Nac.B_inB  = B_N_inN
    Nac.BB_inB = BB_N_inN
    # TODO YAW
    MM_Y       = fBMB(BB_N_inN,Yaw.MM)
    Yaw.r_O    = Twr.r_O + np.dot(Twr.R_0b, r_TN_inT)
    Yaw.R_0b   = R_EN
    Yaw.B      = B_N    
    Yaw.B_inB  = B_N_inN
    Yaw.BB_inB = BB_N_inN
    # TODO Gen

    # ---------------------------------------------
    # Link N-S
    q_psi = q[iPsi,0]
    if main_axis=='x':
        R_NS     = R_z(q_psi + np.pi) 
    elif main_axis=='z':
        R_NS     = R_x(q_psi + np.pi) 
    if not bTiltBeforeNac:
        R_NS     = np.dot(R_y(theta_tilt_y),R_NS)
    R_ES     = np.dot(R_EN, R_NS)
    r_NS     = np.dot(R_EN, r_NS_inN)
    Bx_NS    = np.array([[0],[0],[0]])
    if main_axis=='x':
        Bt_NS    = np.array([[0],[0],[1]])
    elif main_axis=='z':
        Bt_NS    = np.array([[1],[0],[0]])
    B_S      = fBMatRecursion(B_N,Bx_NS,Bt_NS,R_EN,r_NS)
    B_S_inS  = fB_inB(R_ES, B_S)
    BB_S_inS = fB_aug(B_S_inS, Sft.nf)
    MM_S     = fBMB(BB_S_inS,Sft.MM)
    KK_S     = fBMB(BB_S_inS,Sft.KK)

    Sft.r_O    = Nac.r_O + r_NS
    Sft.R_0b   = R_ES
    Sft.B      = B_S    
    Sft.B_inB  = B_S_inS
    Sft.BB_inB = BB_S_inS

    # ---------------------------------------------
    # Link S-B1
    nB   = len(Blds)
    # Point R
    r_SR  = np.dot(R_ES, r_SR_inS)
    B_R = fBMatRecursion(B_S,[],[],R_ES,r_SR)
    B_R_bis = fBMatTranslate(B_S,r_SR)
    # Points B1, B2, B3
    MM_B      = np.zeros((nDOF,nDOF))
    KK_B      = np.zeros((nDOF,nDOF))
    DD_B      = np.zeros((nDOF,nDOF))
    nf_done=0
    nf_tot = sum([B.nf for B in Blds])
    for i,B in enumerate(Blds):
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

        B.r_O    = Sft.r_O + r_SR
        B.B      = B_R
        B.B_inB  = B_B_inB
        B.BB_inB = BB_B_inB
     

    # --- Final assembly
    MM = MM_B.copy()
    MM[:iPsi+1,:iPsi+1] += MM_S
    MM[:Twr.nf,:Twr.nf] += MM_T + MM_N + MM_Y

    KK = KK_B
    KK[:iPsi+1,:iPsi+1] += KK_S
    KK[:Twr.nf,:Twr.nf] += KK_T + KK_N

    DD = DD_B 
    DD[:Twr.nf,:Twr.nf] += DD_T
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
    Struct      = Structure(main_axis=main_axis,theta_cone=theta_cone_y,theta_tilt=theta_tilt_y,bTiltBeforeNac=bTiltBeforeNac)
    Struct.Grd = GroundBody()
    Struct.Twr  = Twr
    Struct.Yaw  = Yaw
    Struct.Nac  = Nac
    Struct.Sft  = Sft
    Struct.Blds = Blds
    Struct.MM   = MM
    Struct.KK   = KK
    Struct.DD   = DD
    Struct.q    = q
    Struct.r_ET_inE=r_ET_inE
    Struct.r_TN_inT=r_TN_inT
    Struct.r_NS_inN=r_NS_inN
    Struct.r_SR_inS=r_SR_inS
    Struct.init_trigger()

    return Struct



if __name__=='__main__':
    np.set_printoptions(linewidth=500)
