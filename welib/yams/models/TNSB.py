"""

TNSB refer to : Tower Nacelle, Shaft, Blades

Assembles the bodies using the Rayleigh-Ritz approximation and joint coordinates. 

Relies on YAMSRec for the bodies

The theory is provided in the reference below. The article also contains an example in the its section, which is reproduced in the file test_TNSB.py.

Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

##
import numpy as np
import copy
import os

from welib.yams.yams_rec import fB_inB, fB_aug, fBMB, fBMatRecursion, fBMatTranslate
from welib.yams.yams_rec import YAMSRecGroundBody, YAMSRecRigidBody
from welib.yams.rotations import R_x, R_y, R_z
from welib.yams.windturbine import WindTurbineStructure, rigidBlades

class TNSBStructure(WindTurbineStructure):

    def __init__(self, 
                 main_axis='x',
                 bTiltBeforeNac=False,
                 ):

        # --- Calling parent constructor
        WindTurbineStructure.__init__(self)
        # From parent
        #self.grd  = None
        #self.twr  = None
        #self.yaw  = None
        #self.nac  = None
        #self.sft  = None
        #self.bld = None
        #self.MM   = None
        #self.KK   = None
        #self.DD   = None
        #self.q    = None
        #self.r_ET_inE=None
        #self.r_TN_inT=None
        #self.r_NS_inN=None
        #self.r_SR_inS=None

        # --- Additional properties not in Parent class
        self.main_axis      = main_axis
        self.nac_yaw        = 0 # TODO
        self.bTiltBeforeNac = bTiltBeforeNac
        self.additional_properties=[] # for user output in __repr__, so I remember what we have in the object


    def compute_RNA(s, verbose=False):
        s.M_rot= sum([B.mass for B in s.bld])
        s.M_RNA= s.M_rot + s.hubgen.mass + s.nac.mass;
        s.r_NGnac_inN = s.nac.s_G_inB.ravel()
        s.r_NGhub_inN = s.r_NS_inN.ravel() + np.dot(s.nac.R_b2g.T, np.dot(s.hubgen.R_b2g, s.hubgen.s_G_inB))

        try:
            # --------------------------------------------------------------------------------}
            # --- New method, requires blade to have "toRigidBody"
            # --------------------------------------------------------------------------------{
            # --- Rigid Blades (with origin R, using N as "global" ref)
            blds_rigid = rigidBlades(s.bld, r_O = [0,0,0]) # TODO blade origins might be wrong in TNSB
            blds_rigid.pos_global = s.r_NR_inN.ravel()
            R_NS = R_y(s.theta_tilt)  # Rotation fromShaft to Nacelle
            blds_rigid.R_b2g      = R_NS


            # --- Creating "hub" with respect to point N (hubgen is wrt S)
            #M_hub  = ED['HubMass']
            #JxxHub_atR = ED['HubIner']
            #hub = RigidBody('Hub', M_hub, (JxxHub_atR,0,0), s_OG=r_SGhub_inS, R_b2g=R_NS, s_OP=r_SR_inS, r_O=r_NS_inN) 
            r_NS_inN=s.hubgen.pos_global.ravel()-s.nac.pos_global.ravel()
            r_SGhub_inS=s.hubgen.masscenter # In body coordinates, S, titled
            #hub = RigidBody('Hub', s.hubgen.mass, J=s.hubgen.masscenter_inertia, s_OG=r_SGhub_inS, R_b2g=R_NS, r_O=r_NS_inN) 

            hub = YAMSRecRigidBody('Hub', s.hubgen.mass, J=s.hubgen.masscenter_inertia, rho_G=s.hubgen.masscenter)
            #hub.shiftOrigin(R_NS.T.dot(-r_NS_inN))
            hub.R_b2g      = R_NS
            hub.pos_global = r_NS_inN

            # --- Rotor = Hub + Blades (with origin R, using N as global ref)
            #rot = blds_rigid.combine(s.hubgen, R_b2g=R_NS, r_O=blds_rigid.pos_global)
            rot = blds_rigid.combine(hub, R_b2g=R_NS, r_O=blds_rigid.pos_global)
            rot.name='rotor'
            #rotgen = rot.combine(gen, R_b2g=R_NS, r_O=blades.pos_global)

            #RNA = rot.combine(gen).combine(nac,r_O=[0,0,0])
            RNA = rot.combine(s.nac,r_O=[0,0,0]).combine(s.yaw, r_O=[0,0,0])
            s.RNA=RNA
            s.r_NGrna_inN=RNA.masscenter # in N
            s.r_NGrot_inN = rot.masscenter

            # Temp storage
            s.blds_rigid = blds_rigid
            s.rot        = rot
            s.hub        = hub

        except:
            # NOTE: due to test test_TNSB_article which uses a dummy YAMSRecBody for a flexbile body, we do not have access to "toRigidBody" routine
            if verbose:
                print('[WARN] TNSB: Fail to compute RNA with new method, using legacy')
            s.r_NGrot_inN = s.r_NR_inN   # NOTE approximation neglecting cone, putting all rotor mass at R
            s.r_NGrna_inN = 1./s.M_RNA * (s.nac.mass*s.r_NGnac_inN + s.hubgen.mass*s.r_NGhub_inN +  s.M_rot*s.r_NGrot_inN)



    def init_trigger(s):
        s.alpha = s.twr.alpha_couplings
        s.iPsi  = s.twr.nf # Index of DOF corresponding to azimuth

        # Useful for load computation
        s.r_NR_inN = s.r_NS_inN.ravel() + np.dot(s.nac.R_b2g.T, np.dot(s.hubgen.R_b2g, s.r_SR_inS.ravel()))
        s.gravity  = s.twr.gravity
        s.compute_RNA()
        s.nDOF = len(s.q)
        s.nShapes_twr = s.twr.nf
        s.nShapes_bld = s.bld[0].nf


    def GF(s,T,x,alpha_y_fact=1):
        """ 
        T is the force along the shaft
        """
        if (s.nShapes_twr!=1):
            raise NotImplementedError('Number of shape function not 1')
        if s.main_axis=='x':
            raise NotImplementedError('Main axis along x')

        # update tower kinematics
        s.twr.gzf=x 
        alpha = s.twr.alpha_couplings
        alpha_y=alpha[1]*alpha_y_fact

        rhoN_x = s.r_NGrna_inN[0,0]
        rhoN_z = s.r_NGrna_inN[2,0]
        rNR_x  = s.r_NR_inN[0,0]
        rNR_z  = s.r_NR_inN[2,0]
        g      = s.gravity
        ux1c   = 1
        vy1c   = s.twr.Bhat_t_bc[1,0]  # Bhat_t_bc[1,j]= self.PhiV[j][0,iNode]
        ky1c   = s.twr.PhiK[0][0,-1]

        Fz_inE =-T*sin(alpha_y + s.theta_tilt) - s.M_RNA*g # TODO potential softening correction

        Fx_inE = T*cos(alpha_y + s.theta_tilt)
#         Fx_inE = T*cos(s.theta_tilt)


        # --- Softening is already in K
        #         k  = x[0] * s.twr.PhiK[0][0,:]
        #         vL = x[0] * vy1c
        #         U  =        s.twr.PhiU[0][0,:]
        #         GF_soft=0
        #         GF_soft+= trapezoid(+U*k*Fz_inE   , s.twr.s_span)
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

        rhoN_x = s.r_NGrna_inN.flatten()[0]
        rhoN_z = s.r_NGrna_inN.flatten()[2]
        rNR_x  = s.r_NR_inN.flatten()[0]
        rNR_z  = s.r_NR_inN.flatten()[2]
        g      = s.gravity
        ux1c   = s.twr.Bhat_x_bc[1,0]
        vy1c   = s.twr.Bhat_t_bc[1,0]  # Bhat_t_bc[1,j]= self.PhiV[j][0,iNode]

        GF  =   T*cos(s.theta_tilt) 
        if bFull:
            GF += - T* vy1c * sin(s.theta_tilt) * x[0]
            GF += (vy1c**2 * s.M_RNA*g * rhoN_z) * x[0]
            GF +=  T*vy1c*(rNR_x*sin(s.theta_tilt) + rNR_z*cos(s.theta_tilt) ) 
            GF += vy1c * s.M_RNA*g * rhoN_x
        return GF

    def __repr__(self):
        s='<TNSB {} object> with fields:\n'.format(type(self).__name__)
        s+=f' - grd twr yaw nac hubgen bld: RigidBody or FASTBeamBody\n'
        s+=f' - MM KK DD: matrices\n'
        s+=f' - q : {self.q.flatten()}\n'
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
        s+=f' * Origin T  : {self.twr.pos_global.T}\n'
        s+=f' * Origin N  : {self.nac.pos_global.T}\n'
        s+=f' * Origin R  : {self.bld[0].pos_global.T}\n'
        s+=f' * Origin S  : {self.hubgen.pos_global.T}\n'
        s+='----------------- RNA ---------------------------------------\n'
        s+=f'M_RNA       {self.M_RNA:.4f}\n'
        s+=f'r_NGrna_inN {np.asarray(self.r_NGrna_inN).flatten()}\n'
        s+=f'     r_NGnac_inN {np.asarray(self.r_NGnac_inN).flatten().round(4)} M_nac {self.nac.mass:.4f}\n'
        s+=f'     r_NGhub_inN {np.asarray(self.r_NGhub_inN).flatten().round(4)} M_hub {self.hubgen.mass:.4f}\n'
        s+=f'     r_NGrot_inN {np.asarray(self.r_NGrot_inN).flatten().round(4)} M_rot {self.M_rot:.4f}\n'
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
        theta_yaw    =  self.nac_yaw
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
        elif self.main_axis=='z':
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
        grd.connectTo(twr, Point=r_ET_inE, Type='Rigid')
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
#         if WT is None:
#             WT      = TNSBStructure(main_axis=main_axis,theta_cone=self.blade_cone,theta_tilt=self.shaft_tilt,bTiltBeforeNac=bTiltBeforeNac)
#         else:
#         WT.main_axis      = main_axis
#         WT.theta_cone     = self.blade_cone
#         WT.theta_tilt     = self.shaft_tilt
#         self.bTiltBeforeNac = bTiltBeforeNac
        self.fixedShaft = fixedShaft # TODO?

        self.grd  = grd # TODO?
        self.MM   = MM
        self.KK   = KK
        self.DD   = DD
        self.q    = q
        self.init_trigger()


    def manual_assembly(self, q=None, DEBUG=False, fixedShaft=False):
        if q is None:
            raise NotImplementedError('Figure out DOF dimension, zero q vector in windturbine')

        # --- Aliases
        sft = self.hubgen
        gen = self.gen
        nac = self.nac
        yaw = self.yawBr
        bld = self.bld
        twr = self.twr
        r_ET_inE     =  self.r_ET_inE
        r_TN_inT     =  self.r_TN_inT    
        r_NS_inN     =  self.r_NS_inN   
        r_SR_inS     =  self.r_SR_inS    

        # Main Parameters
        nDOF = len(q)
        #   CyT=- np.array([ twr.PhiV[0][2,-1],  1.5065E-01, 0, 0]) # End value of shapes functions in y direction
        CxT=  np.zeros(twr.nf)
        CyT=  np.zeros(twr.nf)
        CzT=  np.zeros(twr.nf)
        UxT=  np.zeros(twr.nf)
        UyT=  np.zeros(twr.nf)
        UzT=  np.zeros(twr.nf)
        for j,(u,v) in enumerate(zip(twr.PhiU,twr.PhiV)):
            if self.main_axis=='x':
                CyT[j]=-v[2,-1] # A deflection along z gives a negative angle around y
                CzT[j]= v[1,-1] # A deflection along y gives a positive angle around z # TODO TODO CHECK ME
                UyT[j]= u[1,-1] 
                UzT[j]= u[2,-1] 
                #print('Alpha y - mode {}:'.format(j+1),CyT[j])
            elif self.main_axis=='z':
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
            if self.main_axis=='x':
                Bx_TN = np.array([[0],[0],[UzT[0]]])
            elif self.main_axis=='z':
                Bx_TN = np.array([[UxT[0]],[0],[0]])
            Bt_TN = np.array([[0],[CyT[0]],[0]])
            twr.gzf = q[0,0]
            alpha_y = np.dot(CyT.ravel(), q[0,0].ravel())
        elif twr.nf == 2:
            if self.main_axis=='x':
                Bx_TN = np.array([[0,0],[0,0],[UzT[0],UzT[1]]])
            elif self.main_axis=='z':
                Bx_TN = np.array([[UxT[0],UxT[1]],[0,0],[0,0]])
            twr.gzf = q[0:2,0]
            Bt_TN = np.array([[0,0],[CyT[0],CyT[1]],[0,0]])
            alpha_y = np.dot(CyT.ravel() , q[:2,0].ravel())
        else:
            # TODO use CzT
            raise NotImplementedError()
        #print('alpha_y',alpha_y)
        R_TN     = R_y(alpha_y)
        if self.bTiltBeforeNac:
            R_TN     = np.dot(R_TN, R_y(self.shaft_tilt))
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
        if self.main_axis=='x':
            R_NS     = R_z(q_psi + np.pi) 
        elif self.main_axis=='z':
            R_NS     = R_x(q_psi + np.pi) 
        if not self.bTiltBeforeNac:
            R_NS     = np.dot(R_y(self.shaft_tilt),R_NS)
        R_ES     = np.dot(R_EN, R_NS)
        r_NS     = np.dot(R_EN, r_NS_inN)
        if fixedShaft:
            Bx_NS    = np.array([])
            Bt_NS    = np.array([])
        else:
            Bx_NS    = np.array([[0],[0],[0]])
            if self.main_axis=='x':
                Bt_NS    = np.array([[0],[0],[1]])
            elif self.main_axis=='z':
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
            if self.main_axis=='x':
                R_SB = R_z(0*np.pi + psi_B)
            elif self.main_axis=='z':
                R_SB = R_x(0*np.pi + psi_B)
            R_SB = np.dot(R_SB, R_y(self.blade_cone))
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
        self.grd = YAMSRecGroundBody() # TODO?
        self.MM   = MM
        self.KK   = KK
        self.DD   = DD
        self.q    = q
        self.init_trigger()


if __name__=='__main__':
    np.set_printoptions(linewidth=500)
