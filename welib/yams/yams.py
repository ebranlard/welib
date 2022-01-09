"""

Yams numerical and recursive formulation.

Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

import numpy as np
from .utils import *
from .bodies import Body         as GenericBody
from .bodies import RigidBody    as GenericRigidBody
from .bodies import FlexibleBody as GenericFlexibleBody
from .bodies import BeamBody     as GenericBeamBody
from .bodies import FASTBeamBody as GenericFASTBeamBody
from .bodies import InertialBody as GenericInertialBody

# --- To ease comparison with sympy version
from numpy import eye, cross, cos ,sin
def Matrix(m):
	return np.asarray(m)
def colvec(v): 
    return np.asarray(v).ravel().reshape(3,1)


# --------------------------------------------------------------------------------}
# --- Connections 
# --------------------------------------------------------------------------------{
class Connection():
    def __init__(self, Type, RelPoint=None, RelOrientation=None, JointRotations=None, OrientAfter=True, parentNode=None, parentBody=None):
        if RelOrientation is None:
            RelOrientation=eye(3)
        if RelPoint is None:
            RelPoint=colvec([0,0,0])
        else:
            RelPoint=colvec(RelPoint)

        self.Type=Type
        
        self.s_C_0_inB = RelPoint
        self.s_C_inB   = self.s_C_0_inB
        self.R_ci_0    = RelOrientation
        self.R_ci      = self.R_ci_0     
        self.OrientAfter= OrientAfter
        self.parentNode = parentNode
        self.parentBody = parentBody
        self.I_DOF= None  # < Index of joints DOF in global DOF vector

        if self.Type=='Rigid':
            self.nj=0
        elif self.Type=='SphericalJoint':
            self.JointRotations=JointRotations;
            self.nj=len(self.JointRotations);
        else:
            raise NotImplementedError()

    def updateKinematics(j, q):
        j.B_ci=Matrix(np.zeros((6,j.nj)))
        if j.Type=='Rigid':
            j.R_ci=j.R_ci_0
        elif j.Type=='SphericalJoint':
            R=eye(3)
            myq    = q   [j.I_DOF,0];
            #myqdot = qdot[j.I_DOF];

            for ir,rot in enumerate(j.JointRotations):
                if rot=='x':
                    I=np.array([1,0,0])
                    Rj=R_x( myq[ir] )
                elif rot=='y':
                    I=np.array([0,1,0])
                    Rj=R_y( myq[ir] )
                elif rot=='z':
                    I=np.array([0,0,1])
                    Rj=R_z( myq[ir] )
                else:
                    raise Exception()
                # Setting Bhat column by column
                j.B_ci[3:,ir] = np.dot(R,I) # NOTE: needs to be done before R updates
                # Updating rotation matrix
                R      = np.dot(R , Rj )
                if j.OrientAfter:
                    j.R_ci = np.dot(R, j.R_ci_0 )
                else:
                    j.R_ci = np.dot(j.R_ci_0, R )

        # TODO this is done twice since it's done when parent.updateKinematics is called. CHOSE!
        if j.parentNode is not None:
            #print('>>>> Joint Kinematics. Updating joint position based on parent node position')
            iNode=j.parentNode
            j.s_C_inB = (j.parentBody.s_P[:,iNode]).reshape(3,1)

    def __repr__(self):
        s ='<Connection object>:\n'
        s+='|Properties:\n'
        s+='| - Type:     {} \n'.format(self.Type)
        s+='| - OrientAfter:  {} \n'.format(self.OrientAfter)
        s+='| - s_C_0_inB:(init pos. of conn. in body)   \n{} \n'.format(self.s_C_0_inB)
        s+='| - s_C_inB:  (current pos. of conn. in body)\n{} \n'.format(self.s_C_inB)
        s+='| - R_ci_0:   (init rot. ro conn. in body)   \n{} \n'.format(self.R_ci_0)
        s+='| - R_ci:     (current rot. ro conn. in body)\n{} \n'.format(self.R_ci)
        s+='|Methods: updateKinematics\n'
        s+='|Usefull getters: None \n'
        return s


# --------------------------------------------------------------------------------}
# --- Bodies 
# --------------------------------------------------------------------------------{
class Body(GenericBody):
    def __init__(B,name=''):
        GenericBody.__init__(B, name=name)
        B.Children    = []
        B.Connections = []
        B.MM     = None
        B.B           = [] # Velocity transformation matrix
        B.updatePosOrientation(colvec([0,0,0]), eye(3))
        B.I_DOF = None
#         B.r_O   = None
        B.gzf   = None
#         B.R_0b  = None
#         B.RR_0b = None

    def __repr__(B):
        s='<YAMSRec Body {} object>:\n'.format(B.name)
        s+='|Properties:\n'
        try: 
            names = [c.name for c in B.Children]
        except:
            names=''
        s+='| - Children: {} {}\n'.format(len(B.Children), names)
        try: 
            types = [c.Type for c in B.Connections]
        except:
            types=''
        s+='| - Connections: {} {}\n'.format(len(B.Connections), types)
        s+='| - I_DOF:  {}\n'.format(B.I_DOF)
        s+='| - r_O  :  {}\n'.format(B.r_O.flatten())
        s+='| - gzf  :  {}\n'.format(B.gzf)
        s+='| - R_0b : \n{}\n'.format(B.R_0b)
        s+='| * nf  : {}\n'.format(B.nf)
        s+='| * R_bc: \n{}\n'.format(B.R_bc)
        s+='| * Bhat_x_bc: \n{}\n'.format(B.Bhat_x_bc)
        s+='| * Bhat_t_bc: \n{}\n'.format(B.Bhat_t_bc)
        s+='| * mass     :   {}\n'.format(B.mass) 
        s+='|Methods: connectTo, updateChildrenKinematicsNonRecursive \n'
        s+='|         updatePosOrientation'
        return s

    def updatePosOrientation(o,x_0,R_0b):
        o.r_O = x_0      # position of body origin in global coordinates
        o.R_0b=R_0b      # transformation matrix from body to global

    def connectTo(self, Child, Point=None, Type=None, BodyPoint=None, RelOrientation=None, JointRotations=None, OrientAfter=True):
        if BodyPoint is not None and Type =='Rigid':
            if BodyPoint is not None:
                if BodyPoint == 'FirstPoint':
                    i_C_inB=0;
                elif BodyPoint == 'LastPoint':
                    i_C_inB=self.nSpan-1
                else:
                    raise NotImplementedError()
                RelPoint = self.s_P0[:,i_C_inB].reshape(3,1)
                c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations, OrientAfter=OrientAfter, parentNode=i_C_inB, parentBody=self)
        elif Type =='Rigid':
            c=Connection(Type, RelPoint=Point, RelOrientation = RelOrientation)
        else: # TODO first node, last node
            c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations, OrientAfter=OrientAfter)

        self.Children.append(Child)
        self.Connections.append(c)

    def _setupDOFIndex(o,n):
        nForMe=o.nf
        # Setting my dof index
        o.I_DOF=n+ np.arange(nForMe) 
        # Update
        n=n+nForMe
        for child,conn in zip(o.Children,o.Connections):
            # Connection first
            nForConn=conn.nj;
            conn.I_DOF=n+np.arange(nForConn)
            # Update
            n=n+nForConn;
            # Then Children
            n=child._setupDOFIndex(n)
        return n

    def updateChildrenKinematicsNonRecursive(p,q, qdot=None, qddot=None):
        if p.I_DOF is None:
            raise Exception('Call setupDOFIndex on the reference body first')
        if qdot is None:
            qdot = q*0
        if qddot is None:
            qddot = q*0
        # At this stage all the kinematics of the body p are known
        # Useful variables
        R_0p =  p.R_0b
        B_p  =  p.B
        r_0p  = p.r_O  # Position of body origin in global coordinates

        nf_all_children=sum([child.nf for child in p.Children])

        for ic,(body_i,conn_pi) in enumerate(zip(p.Children,p.Connections)):
            # Flexible influence to connection point
            R_pc  = p.R_bc
            #print('R_pc')
            #print(R_pc)
            Bx_pc = p.Bhat_x_bc
            Bt_pc = p.Bhat_t_bc
            # Joint influence to next body (R_ci, B_ci)
            conn_pi.updateKinematics(q) # TODO
            #print('R_ci',p.name)
            #print(conn_pi.R_ci)

            # Full connection p and j
            R_pi   = np.dot(R_pc, conn_pi.R_ci )
            if conn_pi.B_ci.shape[1]>0:
                Bx_pi  = np.column_stack((Bx_pc, np.dot(R_pc,conn_pi.B_ci[:3,:])))
                Bt_pi  = np.column_stack((Bt_pc, np.dot(R_pc,conn_pi.B_ci[3:,:])))
            else:
                Bx_pi  = Bx_pc
                Bt_pi  = Bt_pc
              
            # Rotation of body i is rotation due to p and j
            R_0i = np.dot( R_0p , R_pi )
            #print('R_pi',p.name)
            #print(R_pi)
            #print('R_0p',p.name)
            #print(R_0p)
            #print('R_0i',p.name)
            #print(R_0i)

            # Position of connection point in P and 0 system
            r_pi_inP= conn_pi.s_C_inB
            r_pi    = np.dot (R_0p , r_pi_inP )
            #print('Body {} to {}'.format(p.name, body_i.name))
            #print('r_pi:',r_pi_inP.flatten())
            #print('r_pi')
            #print(r_pi)
            #print('Bx_pi')
            #print(Bx_pi)
            #print('Bt_pi')
            #print(Bt_pi)
            B_i      = fBMatRecursion(B_p, Bx_pi, Bt_pi, R_0p, r_pi)
            B_i_inI  = fB_inB(R_0i, B_i)
            BB_i_inI = fB_aug(B_i_inI, body_i.nf)

            body_i.B      = B_i    
            body_i.B_inB  = B_i_inI
            body_i.BB_inB = BB_i_inI

            # --- Updating Position and orientation of child body 
            r_0i = r_0p + r_pi  # in 0 system
            body_i.R_pb = R_pi 
            body_i.updatePosOrientation(r_0i,R_0i)

            # TODO flexible dofs and velocities/acceleration
            body_i.gzf  = q[body_i.I_DOF,0] # TODO use updateKinematics

#            TODO TODO TODO: remaining from matlab?????
            gzf  = q[body_i.I_DOF,0]
#             gz   = q    (i.I_DOF);
#             gzp  = qdot (i.I_DOF);
#             gzpp = qddot(i.I_DOF);
# 
#             % Velocity
#             v_pi_inP  = Bx_pi*qdot(IHat)    ;
#             om_pi_inP = Bt_pi*qdot(IHat)    ;
# 
#             % All velocities
#             %v_0(1:3)      = B_i(1:3,:) * qdot(1:n_to_i);
#             %v_0(4:6)      = B_i(4:6,:) * qdot(1:n_to_i);
#             %v_0(7:6+i.nf) = qdot(n_to_i+1:i.nf);
#             v_i_in0 = BB_i_inI * qdot(I_with_i);
            v_i_in0 = np.zeros((6+body_i.nf,1))
#             %  
#             % V-Accelerations in P and 0
#             a_i_v_in0=zeros(6+i.nf,1); % All v accelerations
#             a_i_v_inP =p.a_O_v_inB + ...
#                      cross(p.om_O_inB,cross(p.om_O_inB, r_pi_inP) )+...
#                      2*cross(p.om_O_inB,  v_pi_inP ) + ...
#                      cross(p.omp_O_v_inB, r_pi_inP);
# 
#             omp_i_v_inP =p.omp_O_v_inB +  cross(p.om_O_inB, om_pi_inP);
#             omp_i_v_inI_bis = R_pi'*omp_i_v_inP; % <<<<<<<< ALTERNATIVE
#             omp_i_v_inI =zeros(3,1);
#             dom=zeros(3,1);
#             om=zeros(3,1);
#             for j=size(B_i_inI,2):-1:1
#                 dom         = B_i_inI(4:6,j)*qdot(j)      ;
#                 omp_i_v_inI = omp_i_v_inI +  cross(dom,om);
#                 om          = om + dom                    ;
#             end
#             om_in0     = B_i(4:6,:) * qdot(1:n_to_i);
#             om_inI_bis = R_0i *om_in0; % <<<<<<<< ALTERNATIVE
# 
#             % Accelerations in 0
            a_i_v_in0=np.zeros((6+body_i.nf,1))
#             a_i_v_in0(1:3)=R_0p*a_i_v_inP;
#             a_i_v_in0(4:6)=R_0p*omp_i_v_inP;
#             a_i_v_in0(4:6)=R_0i*omp_i_v_inI; % <<<<<<<< ALTERNTIVE
#             a_i_v_in0(7:6+i.nf)=gzpp;
# 
# 
#             i.R_pb = R_pi ;
            body_i.updateKinematics(r_0i, R_0i, gzf, v_i_in0, a_i_v_in0)
#             if ~isequal(i.Type,'Rigid')
#                 if i.nf>0
#                     i.computeMassMatrix(); % Triggers 
#                     i.computeInertiaForces();
#                 end
#             end
# %             nf_N=0;
# %             BB_N_inN =[B_N_inN zeros(6,nf_N) ; zeros(nf_N, size(B_N_inN,2)) eye(nf_N)];
# %             % Update of mass matrix
# %             MM_N= BB_N_inN'*Nac.MM*BB_N_inN;
# %             KK_N= BB_N_inN'*Nac.KK*BB_N_inN;


    def updateKinematics(o,x_0,R_0b,gz,v_0,a_v_0):
        # NOTE: this is overriden by BeamBody
        # Updating position of body origin in global coordinates
        o.r_O = x_0[0:3]
        o.gzf = gz
        # Updating Transformation matrix
        o.R_0b=R_0b
        # Updating rigid body velocity and acceleration
        o.v_O_inB     = np.dot(R_0b, v_0[0:3])
        o.om_O_inB    = np.dot(R_0b, v_0[3:6])
        o.a_O_v_inB   = np.dot(R_0b, a_v_0[0:3])
        o.omp_O_v_inB = np.dot(R_0b, a_v_0[3:6])

    @property
    def _positions_global(B): # todo rename
        # NOTE: this is overriden by BeamBody
        return B.r_O # for rigid bodies



    def _getFullM(o,M):
        if isinstance(o,GroundBody):
            raise Exception('Not intended to be called for Ground body')
        MqB      = fBMB(o.BB_inB,o.MM)
        n        = MqB.shape[0]
        M[:n,:n] = M[:n,:n]+MqB     
        for c in o.Children:
            M=c._getFullM(M)
        return M
        
    def _getFullK(o,K):
        if isinstance(o,GroundBody):
            raise Exception('Not intended to be called for Ground body')
        KqB      = fBMB(o.BB_inB,o.KK)
        n        = KqB.shape[0]
        K[:n,:n] = K[:n,:n]+KqB     
        for c in o.Children:
            K=c._getFullK(K)
        return K
        
    def _getFullD(o,D):
        if isinstance(o,GroundBody):
            raise Exception('Not intended to be called for Ground body')
        DqB      = fBMB(o.BB_inB,o.DD)
        n        = DqB.shape[0]
        D[:n,:n] = D[:n,:n]+DqB     
        for c in o.Children:
            D=c._getFullD(D)
        return D

    @property
    def bodies(o):
        """ List of bodies recursively (children of children, etc)"""
        bodies = [o]
        for c in o.Children:
            bodies += c.bodies
        return bodies


    @property
    def R_bc(self):
        return eye(3);
    @property
    def Bhat_x_bc(self):
        return Matrix(np.zeros((3,0)))
    @property
    def Bhat_t_bc(self):
        return Matrix(np.zeros((3,0)))
    @property
    def nf(B):
        if hasattr(B,'PhiU'):
            return len(B.PhiU)
        else:
            return 0
    @property
    def mass(B):
        if B.MM is None:
            return 0
        return B.MM[0,0]


# --------------------------------------------------------------------------------}
# --- Ground Body 
# --------------------------------------------------------------------------------{
class GroundBody(Body, GenericInertialBody):
    def __init__(B):
        Body.__init__(B, 'Grd')
        GenericInertialBody.__init__(B)
        # We'll use the GroundBody object for global "assembly"
        B.nq = 0
        B.q  = []

    def setupDOFIndex(o):
        n=0
        o.nq = o._setupDOFIndex(n)
        return o.nq

    def setDOF(o, q):  
        q   = q.reshape(o.nq,1)
        o.q = q
        if o.I_DOF is None:
            o.setupDOFIndex()
        # Update kinematics of all bodies
        for b in o.bodies:
            b.updateChildrenKinematicsNonRecursive(o.q)


    def __repr__(self):
        s='<YAMSRec GroundBody {} object>:\n'.format(self.name)
        #s+='|Inherits:\n'
        #s+='||'+'\n|'.join(Body.__repr__(self).split('\n'))+'\n'
        #s+='||'+'\n|'.join(GenericInertialBody.__repr__(self).split('\n'))+'\n'
        s+='|Properties:\n'
        s+='|- nq: {}\n'.format(self.nq)
        s+='|- q:  {}\n'.format(self.q)
        try:
            bnames = [b.name for b in self.bodies]
        except:
            bnames=''
        s+='|* bodies: {}\n'.format(len(self.bodies))
        for b in self.bodies:
            s+='|    name : {:8s}, I_DOF : {}\n'.format(b.name, b.I_DOF)
        s+='|Derived properties: M, K, D'
        return s

    def eva(o):
        """ Perform eigenvalue analysis based on system matrices"""
        from welib.system.eva import eigMCK
        MM = o.M
        KK = o.K
        DD = o.D
        freq_d, zeta, Q, freq_0 = eigMCK(MM, DD, KK, method='full_matrix', sort=True)
        return freq_d, zeta, Q, freq_0

    def modes(o, norm='tip_norm'):
        # Backup current q
        q_before = o.q
        # Perform EVA
        freq_d, zeta, Q, freq_0 = o.eva()
        # Apply each mode, to compute full position of structure
        Modes=[]
        for q in Q.T: # loop though columns
            o.setDOF(q)
            mode = o._all_positions_global
            Modes.append(mode)
            # --- Mode scaling
            # TODO figure out main "dimension"
            # Sript below assumes x is main dimension
            # TODO TODO normalization is not bullet proof..
            Uy = mode[1,:]
            Uz = mode[2,:]
            maxAmp  = [np.max(np.abs(Uy)), np.max(np.abs(Uz))]
            iMaxAmp = np.mod(np.argmax(maxAmp)+1,3)
            iOther  = 1 if iMaxAmp==2 else 2
            if norm=='tip_norm':
                tipVal = mode[iMaxAmp, -1]
                fact = tipVal
                mode[iMaxAmp,:] /= fact
                mode[iOther,:] /= fact
            elif norm=='max':
                maxVal = np.max(np.abs(mode[iMaxAmp, :]))
                iMaxVal = np.argmax(np.abs(mode[iMaxAmp,:]))
                fact = 1 / mode[iMaxAmp, iMaxVal]
                mode[iMaxAmp,:] /= fact
                mode[iOther,:] /= fact
            elif norm=='mode_mass':
                raise Exception()
                pass


        # Restore current q
        o.setDOF(q_before)

        return Modes

    @property
    def _all_positions_global(o): # TODO rename 
        """ Return shape of full structure"""
        for ib, b in enumerate(o.bodies):
            pos = b._positions_global
            if ib==0:
                all_pos = pos
            else:
                all_pos = np.column_stack((all_pos,pos))
        return all_pos


    @property
    def M(o):
        M = np.zeros((o.nq, o.nq))
        for c in o.Children:
            M=c._getFullM(M)
        return M
        
    @property
    def K(o):
        K = np.zeros((o.nq, o.nq))
        for c in o.Children:
            K=c._getFullK(K)
        return K
        
    @property
    def D(o):
        D = np.zeros((o.nq, o.nq))
        for c in o.Children:
            D=c._getFullD(D)
        return D


# --------------------------------------------------------------------------------}
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class RigidBody(Body,GenericRigidBody):
    def __init__(B, name, mass, J_G, rho_G):
        """
        Creates a rigid body 
        """
        Body.__init__(B,name)
        GenericRigidBody.__init__(B, name, mass, J_G, rho_G)
        B.s_G_inB = B.masscenter
        B.J_G_inB = B.masscenter_inertia
        B.J_O_inB = translateInertiaMatrixFromCOG(B.J_G_inB, mass, -B.s_G_inB)
        B.MM = rigidBodyMassMatrix(mass, B.J_O_inB, B.s_G_inB) # TODO change interface
        B.DD = np.zeros((6,6))
        B.KK = np.zeros((6,6))



# --------------------------------------------------------------------------------}
# --- Beam Body 
# --------------------------------------------------------------------------------{
class BeamBody(GenericBeamBody, Body):
    def __init__(B, s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=None, s_G0=None, 
            s_min=None, s_max=None,
            bAxialCorr=False, bOrth=False, Mtop=0, bStiffening=True, gravity=None,main_axis='z',
            massExpected=None,
            name='dummy'
            ):
        """ 
          Points P0 - Undeformed mean line of the body
        """
        # --- nherit from BeamBody and Body 
        Body.__init__(B)
        GenericBeamBody.__init__(B,name, s_span, s_P0, m, EI, PhiU, PhiV, PhiK, jxxG=jxxG, s_G0=s_G0, s_min=s_min, s_max=s_max,
                 bAxialCorr=bAxialCorr, bOrth=bOrth, Mtop=Mtop, bStiffening=bStiffening, gravity=gravity, main_axis=main_axis,
                 massExpected=massExpected
                )

        B.gzf   = np.zeros((B.nf,1))
        B.gzpf  = np.zeros((B.nf,1))
        B.gzppf = np.zeros((B.nf,1))

        # TODO
        B.V0         = np.zeros((3,B.nSpan))
        B.K0         = np.zeros((3,B.nSpan))
        B.rho_G0_inS = np.zeros((3,B.nSpan)) # location of COG in each cross section
        #[o.PhiV,o.PhiK] = fBeamSlopeCurvature(o.s_span,o.PhiU,o.PhiV,o.PhiK,1e-2);
        #[o.V0,o.K0]     = fBeamSlopeCurvature(o.s_span,o.s_P0,o.V0,o.K0,1e-2)    ;
        #if isempty(o.s_G0); o.s_G0=o.s_P0; end;
        #if isempty(o.rho_G0_inS); o.rho_G0_inS=np.zeros(3,o.nSpan); end;
        #if isempty(o.rho_G0    ); 
        #    o.rho_G0 =np.zeros(3,o.nSpan);
        #    for i=1:o.nSpan
        #        o.rho_G0(1:3,i) =R_x(o.V0(1,i))*o.rho_G0_inS(:,i);

    @property
    def alpha_couplings(self):
        return  np.dot(self.Bhat_t_bc , self.gzf).ravel()

    @property
    def R_bc(self):
        alpha = self.alpha_couplings

        if self.main_axis=='x':
            return np.dot(R_y(alpha[1]),R_z(alpha[2]))
        elif self.main_axis=='z':
            return np.dot(R_x(alpha[0]),R_y(alpha[1]))
        else:
            raise NotImplementedError()

    def updateKinematics(o,x_0,R_0b,gz,v_0,a_v_0):
        super(BeamBody,o).updateKinematics(x_0,R_0b,gz,v_0,a_v_0)
        # --- Calculation of deformations wrt straight beam axis, curvature (K) and velocities (UP)
        if o.nf>0:
            o.gzpf  = v_0[6:]
            o.gzppf = a_v_0[6:]
            # Deflections shape
            o.U  = np.zeros((3,o.nSpan));
            o.V  = np.zeros((3,o.nSpan));
            o.K  = np.zeros((3,o.nSpan));
            #o.U(1,:) = o.s_span; 
            o.UP = np.zeros((3,o.nSpan));
            for j in range(o.nf):
                o.U [0:3,:] = o.U [0:3,:] + o.gzf[j]  * o.PhiU[j][0:3,:]
                o.UP[0:3,:] = o.UP[0:3,:] + o.gzpf[j] * o.PhiU[j][0:3,:]
                o.V [0:3,:] = o.V [0:3,:] + o.gzf[j]  * o.PhiV[j][0:3,:]
                o.K [0:3,:] = o.K [0:3,:] + o.gzf[j]  * o.PhiK[j][0:3,:]
            o.V_tot=o.V+o.V0;
            o.K_tot=o.K+o.K0;

            # Position of mean line in body coordinates
            o.s_P=o.s_P0+o.U;

            # Position of deflected COG in body coordinates
            # TODO TODO TODO mean_axis not x
            o.rho_G      = np.zeros((3,o.nSpan))
            if o.main_axis=='x':
                o.rho_G[1,:] = o.rho_G0_inS[1,:]*np.cos(o.V_tot[0,:])-o.rho_G0_inS[2,:]*np.sin(o.V_tot[0,:]);
                o.rho_G[2,:] = o.rho_G0_inS[1,:]*np.sin(o.V_tot[0,:])+o.rho_G0_inS[2,:]*np.cos(o.V_tot[0,:]);
            else:
                print('>>>> YAMS: NotImplemented beam along z, wathc out for your results.')
                #raise NotImplementedError()
                #o.rho_G[1,:] = o.rho_G0_inS[1,:]*np.cos(o.V_tot[0,:])-o.rho_G0_inS[2,:]*np.sin(o.V_tot[0,:]);
                #o.rho_G[2,:] = o.rho_G0_inS[1,:]*np.sin(o.V_tot[0,:])+o.rho_G0_inS[2,:]*np.cos(o.V_tot[0,:]);
            o.s_G = o.s_P+o.rho_G; 
            # Alternative:
            #rho_G2     = zeros(3,o.nSpan);
            #rho_G2(2,:) = o.rho_G0(2,:).*cos(o.V(1,:))-o.rho_G0(3,:).*sin(o.V(1,:));
            #rho_G2(3,:) = o.rho_G0(2,:).*sin(o.V(1,:))+o.rho_G0(3,:).*cos(o.V(1,:));
            #compare(o.rho_G,rho_G2,'rho_G');
            # Position of connection point
            for ic, conn in enumerate(o.Connections):
                if conn.parentNode is not None:
                    # TODO: this is done twice see Connection
                    iNode=conn.parentNode;
                    conn.s_C_inB = o.s_P[:,iNode]

    @property
    def _positions_global(B): # TODO rename
        displ_g = B.R_0b.dot(B.s_P) # TODO reference line or COG
        pos_g   = B.r_O + displ_g
        return pos_g

    @property
    def nSpan(B):
        return len(B.s_span)

    def __repr__(self):
        s='<YAMSRec BeamBody {} object>:\n'.format(self.name)
        s+='|Inherits:\n'
        s+='|'+'\n||'.join(Body.__repr__(self).split('\n'))+'\n'
        #s+='|'+'\n||'.join(GenericBeamBody.__repr__(self).split('\n'))+'\n'
        s+='|Properties:\n'
        s+='| * nSpan: {}\n'.format(self.nSpan)
        return s


# --------------------------------------------------------------------------------}
# --- Uniform Beam Body 
# --------------------------------------------------------------------------------{
class UniformBeamBody(BeamBody):
    def __init__(B, name, nShapes, nSpan, L, EI0, m, Mtop=0, jxxG=None, GKt=None, bAxialCorr=True, bCompatibility=False, bStiffnessFromGM=False, bStiffening=True, gravity=None, main_axis='x'):

        import welib.beams.theory as bt
        if jxxG is None:
            jxxG=0
        if GKt is None:
            GKt=0
        if gravity is None and (bStiffening or Mtop>0):
            raise Exception('`gravity` is None, but `bStiffening` is true or `Mtop`>0. Please provide `gravity`.')

        A=1; rho=A*m;
        x=np.linspace(0,L,nSpan);
        # Mode shapes
        freq,s_span,U,V,K = bt.UniformBeamBendingModes('unloaded-topmass-clamped-free',EI0,rho,A,L,x=x,Mtop=Mtop, nModes=nShapes)
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        PhiV = np.zeros((nShapes,3,nSpan)) # Slope
        PhiK = np.zeros((nShapes,3,nSpan)) # Curvature
        if main_axis=='x':
            iModeAxis=2      # Setting modes along z
        elif main_axis=='z':
            iModeAxis=0      # Setting modes along x
        for j in np.arange(nShapes):  
            PhiU[j][iModeAxis,:] = U[j,:] 
            PhiV[j][iModeAxis,:] = V[j,:]
            PhiK[j][iModeAxis,:] = K[j,:]
        m       = m    * np.ones(nSpan)
        jxxG    = jxxG * np.ones(nSpan)
        EI      = np.zeros((3,nSpan))
        if main_axis=='x':
            EI[1,:] = EI0
            EI[2,:] = EI0
        elif main_axis=='z':
            EI[0,:] = EI0
            EI[1,:] = EI0

        GKt     = GKt  * np.ones(nSpan)
        
        # --- Straight undeflected shape (and COG)
        s_P0      = np.zeros((3,nSpan))
        if main_axis=='x':
            s_P0[0,:] = x
        elif main_axis=='z':
            s_P0[2,:] = x

	# Create a beam body
        super(UniformBeamBody,B).__init__(s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=jxxG, bAxialCorr=bAxialCorr, Mtop=Mtop, bStiffening=bStiffening,
                gravity=gravity, main_axis=main_axis, name=name)


    def __repr__(self):
        s='<YAMSRec UniformBeamBody {} object>:\n'.format(self.name)
        s+='|Inherits:\n'
        s+='|'+'\n||'.join(BeamBody.__repr__(self).split('\n'))+'\n'
        return s


# --------------------------------------------------------------------------------}
# --- FAST Beam body 
# --------------------------------------------------------------------------------{
class FASTBeamBody(BeamBody, GenericFASTBeamBody):
    def __init__(B, body_type, ED, inp, Mtop=0, shapes=None, nShapes=None, main_axis='x',nSpan=None,bAxialCorr=False,bStiffening=True, 
            spanFrom0=False, massExpected=None, gravity=None
            ):
        """ 
        """
        if shapes is None:
            if nShapes==2:
                shapes=[0,1]
            elif nShapes==0:
                shapes=[]
            elif nShapes==1:
                shapes=[0]
            else:
                raise NotImplementedError('>> TODO')
        GenericFASTBeamBody.__init__(B, ED, inp, Mtop=Mtop, shapes=shapes, main_axis=main_axis, nSpan=nSpan, bAxialCorr=bAxialCorr, bStiffening=bStiffening, 
                spanFrom0=spanFrom0,
                massExpected=massExpected,
                gravity=gravity
                )
        # We need to inherit from "YAMS" Beam not just generic Beam
        BeamBody.__init__(B, B.s_span, B.s_P0, B.m, B.PhiU, B.PhiV, B.PhiK, B.EI, jxxG=B.jxxG, s_G0=B.s_G0, 
                # NOTE: r_O, r_b2g is lost here
                s_min=B.s_min, s_max=B.s_max,
                bAxialCorr=bAxialCorr, bOrth=B.bOrth, Mtop=Mtop, bStiffening=bStiffening, gravity=B.gravity,main_axis=main_axis,
                massExpected=massExpected
                )

# --------------------------------------------------------------------------------}
# --- B Matrices 
# --------------------------------------------------------------------------------{
def fB_inB(R_EI, B_I):
    """ Transfer a global B_I matrix (body I at point I) into a matrix in it's own coordinate.
    Simply multiply the top part and bottom part of the B matrix by the 3x3 rotation matrix R_EI
    e.g.
         B_N_inN = [R_EN' * B_N(1:3,:);  R_EN' * B_N(4:6,:)];
    """ 
    if len(B_I)==0:
        B_I_inI = Matrix(np.array([]))
    else:
        B_I_inI = Matrix(np.vstack(( np.dot(R_EI.T, B_I[:3,:]), np.dot(R_EI.T , B_I[3:,:]))))
    return B_I_inI

def fB_aug(B_I_inI, nf_I, nf_Curr=None, nf_Prev=None):
    """
    Augments the B_I_inI matrix, to include nf_I flexible degrees of freedom.
    This returns the full B matrix on the left side of Eq.(11) from [1], 
    based on the Bx and Bt matrices on the right side of this equation
    """
    if len(B_I_inI)==0:
        if nf_I>0:
            BB_I_inI = Matrix(np.vstack( (np.zeros((6,nf_I)), np.eye(nf_I))) )
        else:
            BB_I_inI= Matrix(np.zeros((6,0)))
    else:
        if nf_Curr is not None:
            # Case of several flexible bodies connected to one point (i.e. blades)
            nf_After=nf_I-nf_Prev-nf_Curr
            I = np.block( [np.zeros((nf_Curr,nf_Prev)), np.eye(nf_Curr), np.zeros((nf_Curr,nf_After))] )
        else:
            nf_Curr=nf_I
            I=np.eye(nf_I)

        BB_I_inI = np.block([ [B_I_inI, np.zeros((6,nf_I))], [np.zeros((nf_Curr,B_I_inI.shape[1])), I]]);

    return Matrix(BB_I_inI)


def fBMatRecursion(Bp, Bhat_x, Bhat_t, R0p, r_pi):
    """ Recursive formulae for B' and Bhat 
    See discussion after Eq.(12) and (15) from [1]
    """
    # --- Safety checks
    if len(Bp)==0:
        n_p = 0
    elif len(Bp.shape)==2:
        n_p = Bp.shape[1]
    else:
        raise Exception('Bp needs to be empty or a 2d array')
    if len(Bhat_x)==0:
        ni = 0
    elif len(Bhat_x.shape)==2:
        ni = Bhat_x.shape[1]
    else:
        raise Exception('Bi needs to be empty or a 2d array')

    r_pi=r_pi.reshape(3,1)

    # TODO use Translate here
    Bi = Matrix(np.zeros((6,ni+n_p)))
    for j in range(n_p):
        Bi[:3,j] = Bp[:3,j]+cross(Bp[3:,j],r_pi.ravel()) # Recursive formula for Bt mentioned after Eq.(15)
        Bi[3:,j] = Bp[3:,j] # Recursive formula for Bx mentioned after Eq.(12)
    if ni>0:
        Bi[:3,n_p:] = np.dot(R0p, Bhat_x[:,:]) # Recursive formula for Bx mentioned after Eq.(15)
        Bi[3:,n_p:] = np.dot(R0p, Bhat_t[:,:]) # Recursive formula for Bt mentioned after Eq.(12)
    return Bi

def fBMatTranslate(Bp,r_pi):
    """
    Rigid translation of a B matrix to another point, i.e. transfer the velocities from a point to another: 
      - translational velocity:  v@J = v@I + om@I x r@IJ
      - rotational velocity   : om@J = om@I
    """
    Bi=np.zeros(Bp.shape)
    if Bp.ndim==1:
        raise NotImplementedError

    for j in range(Bp.shape[1]):
        Bi[0:3,j] = Bp[0:3,j]+np.cross(Bp[3:6,j],r_pi.ravel());
        Bi[3:6,j] = Bp[3:6,j]
    return Bi


def fBMB(BB_I_inI,MM):
    """ Computes the body generalized matrix: B'^t M' B 
    See Eq.(8) of [1] 
    """
    MM_I = np.dot(np.transpose(BB_I_inI), MM).dot(BB_I_inI)
    return MM_I

if __name__=='__main__':
    pass
