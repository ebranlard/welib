"""

Yams numerical and recursive formulation.

NOTE: this file is intended to be comparable with yams_sympy.py

Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

import numpy as np
from .utils import buildRigidBodyMassMatrix
from .utils import translateInertiaMatrixFromCOG
from .bodies import Body         as GenericBody
from .bodies import RigidBody    as GenericRigidBody
from .bodies import FlexibleBody as GenericFlexibleBody
from .bodies import BeamBody     as GenericBeamBody
from .bodies import FASTBeamBody as GenericFASTBeamBody
from .bodies import InertialBody as GenericInertialBody

# --- To ease comparison with sympy version
# from numpy import eye, cross, cos ,sin

import sympy as sp
from sympy import Matrix, symbols

# --------------------------------------------------------------------------------}
# --- Sympy harmony 
# --------------------------------------------------------------------------------{
def R_x(t):
    if isinstance(t, sp.Basic):
        return Matrix( [[1,0,0], [0,sp.cos(t),-sp.sin(t)], [0,sp.sin(t),sp.cos(t)]])
    else:
        return np.array( [[1,0,0], [0,np.cos(t),-np.sin(t)], [0,np.sin(t),np.cos(t)]])

def R_y(t):
    if isinstance(t, sp.Basic):
        return Matrix( [[sp.cos(t),0,sp.sin(t)], [0,1,0], [-sp.sin(t),0,sp.cos(t)] ])
    else:
        return np.array( [[np.cos(t),0,np.sin(t)], [0,1,0], [-np.sin(t),0,np.cos(t)] ])

def R_z(t):
    if isinstance(t, sp.Basic):
        return Matrix( [[sp.cos(t),-sp.sin(t),0], [sp.sin(t),sp.cos(t),0], [0,0,1]])
    else:
        return np.array( [[np.cos(t),-np.sin(t),0], [np.sin(t),np.cos(t),0], [0,0,1]])

def cross(u, v):
    if hasattr(u, "cross"):
        return u.cross(v)
    else:
        return np.cross(u, v)


# --------------------------------------------------------------------------------}
# --- Connections 
# --------------------------------------------------------------------------------{
class Connection():
    def __init__(self, Type, RelPoint=None, RelOrientation=None, JointRotations=None, OrientAfter=True, parentNode=None, parentBody=None, sympy=False):
        self.sympy = sympy
        if RelOrientation is None:
            RelOrientation=self.eye(3)
        if RelPoint is None:
            RelPoint = [0,0,0]

        self.Type=Type
        
        self.s_C_0_inB = self.vec3(RelPoint)
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

    # --- Generic Tools to work with Sympy and Numpy
    def vec3(self, v):
        if self.sympy:
            return Matrix([[v[0]],[v[1]],[v[2]]])
        else:
            v = np.asarray(v).ravel()
            if len(v)!=3:
                raise Exception('Vector should be of length 3')
            return v

    def Matrix(self, m):
        if self.sympy:
            return Matrix(m)
        else:
            return np.asarray(m)

    def cross(self, V1, V2):
        if self.sympy:
            return [V1[1]*V2[2]-V1[2]*V2[1], V1[2]*V2[0]-V1[0]*V2[2], (V1[0]*V2[1]-V1[1]*V2[0]) ]
        else:
            return np.cross(V1, V2) 

    def eye(self, n): 
        if self.sympy:
            return Matrix( np.eye(n).astype(int) )
        else:
            return np.eye(n)
    # --- End generic tools


    def updateKinematics(j, q):
        j.B_ci = j.Matrix(np.zeros((6,j.nj)))
        if j.Type=='Rigid':
            j.R_ci=j.R_ci_0
        elif j.Type=='SphericalJoint':
            R = np.eye(3)
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
                j.B_ci[3:,ir] = R @ I # NOTE: needs to be done before R updates
                # Updating rotation matrix
                R      = R @ Rj
                if j.OrientAfter:
                    j.R_ci = j.Matrix(R @ j.R_ci_0)
                else:
                    j.R_ci = j.Matrix(j.R_ci_0 @ R)

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
class YAMSRecBody(GenericBody):
    def __init__(B, name='', sympy=False):
        GenericBody.__init__(B, name=name, sympy=sympy)
        B.Children    = []
        B.Connections = []
        B.MM          = None
        B.B           = []     # Velocity transformation matrix
        B.B_inB       = None
        B.BB_inB      = None
#         B.Bhat_x_bc   = None
#         B.Bhat_t_bc   = None
        B.I_DOF       = None
        B.gzf         = None

    def __repr__(B):
        s='<YAMSRecBody {} object>:\n'.format(B.name)
        s+='|Inherits:\n'
        s+='|'+'\n||'.join(GenericBody.__repr__(B).split('\n'))+'--->\n'
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
        s+='| - gzf       :  {}\n'.format(B.gzf)
        s+='| * nf  : {}\n'.format(B.nf)
        s+='| * R_bc: \n{}\n'.format(B.R_bc)
        s+='| * Bhat_x_bc: \n{}\n'.format(B.Bhat_x_bc)
        s+='| * Bhat_t_bc: \n{}\n'.format(B.Bhat_t_bc)
        s+='|Methods: connectTo, updateChildrenKinematicsNonRecursive \n'
        return s

    def connectTo(self, Child, Point=None, Type=None, BodyPoint=None, RelOrientation=None, JointRotations=None, OrientAfter=True):
        if BodyPoint is not None and Type =='Rigid':
            if BodyPoint is not None:
                if BodyPoint == 'FirstPoint':
                    i_C_inB=0;
                elif BodyPoint == 'LastPoint':
                    i_C_inB=self.nSpan-1
                else:
                    raise NotImplementedError()
                RelPoint = self.s_P0[:,i_C_inB]
                c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations, OrientAfter=OrientAfter, parentNode=i_C_inB, parentBody=self, sympy=self.sympy)
        elif Type =='Rigid':
            c=Connection(Type, RelPoint=Point, RelOrientation = RelOrientation, sympy=self.sympy)
        else: # TODO first node, last node
            c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations, OrientAfter=OrientAfter, sympy=self.sympy)

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
        R_0p =  p.R_b2g
        B_p  =  p.B
        r_0p  = p.pos_global  # Position of body origin in global coordinates

        nf_all_children=sum([child.nf for child in p.Children])

        for ic,(body_i,conn_pi) in enumerate(zip(p.Children,p.Connections)):
            # Flexible influence to connection point
            R_pc  = p.R_bc
            Bx_pc = p.Bhat_x_bc
            Bt_pc = p.Bhat_t_bc
            # Joint influence to next body (R_ci, B_ci)
            conn_pi.updateKinematics(q) # TODO

            # Full connection p and j
            R_pi   = R_pc @ conn_pi.R_ci
            if conn_pi.B_ci.shape[1]>0:
                Bx_pi  = p.Matrix(np.column_stack((Bx_pc, R_pc @ conn_pi.B_ci[:3,:])))
                Bt_pi  = p.Matrix(np.column_stack((Bt_pc, R_pc @ conn_pi.B_ci[3:,:])))
            else:
                Bx_pi  = Bx_pc
                Bt_pi  = Bt_pc
              
            # Rotation of body i is rotation due to p and j
            R_0i = R_0p @ R_pi

            # Position of connection point in P and 0 system
            r_pi_inP= conn_pi.s_C_inB
            r_pi    = R_0p @ r_pi_inP 
            B_i      = fBMatRecursion(B_p, Bx_pi, Bt_pi, R_0p, r_pi, sympy=p.sympy)
            B_i_inI  = fB_inB(R_0i, B_i, sympy=p.sympy)
            BB_i_inI = fB_aug(B_i_inI, body_i.nf, sympy=p.sympy)

            body_i.B      = B_i    
            body_i.B_inB  = B_i_inI
            body_i.BB_inB = BB_i_inI

            # --- Updating Position and orientation of child body 
            r_0i = r_0p + r_pi  # in 0 system
            body_i.R_pb = R_pi 
            body_i.pos_global = r_0i
            body_i.R_b2g = R_0i

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


    def updateKinematics(o, x_0=None, R_b2g=None, gz=None, v_0=None, a_v_0=None):
        # NOTE: this is overriden by BeamBody
        # Updating position of body origin in global coordinates
        if x_0 is not None:
            o.pos_global = x_0[0:3]
        if gz is not None:
            o.gzf = gz
        # Updating Transformation matrix
        if R_b2g is not None:
            o.R_b2g=R_b2g
        else:
            R_b2g = o.R_b2g
        # Updating rigid body velocity and acceleration
        if v_0 is not None:
            o.v_O_inB     = R_b2g @ v_0[0:3]
            o.om_O_inB    = R_b2g @ v_0[3:6]
        if a_v_0 is not None:
            o.a_O_v_inB   = R_b2g @ a_v_0[0:3]
            o.omp_O_v_inB = R_b2g @ a_v_0[3:6]

    @property
    def _positions_global(B): # todo rename
        # NOTE: this is overriden by BeamBody
        return B.pos_global # for rigid bodies

    def _getFullM(o, M):
        if isinstance(o, YAMSRecGroundBody):
            raise Exception('Not intended to be called for Ground body')
        MqB      = fBMB(o.BB_inB, o.MM, sympy=o.sympy, name=o.name)
        n        = MqB.shape[0]
        M[:n,:n] = M[:n,:n]+MqB     
        for c in o.Children:
            M=c._getFullM(M)
        return M
        
    def _getFullK(o, K):
        if isinstance(o, YAMSRecGroundBody):
            raise Exception('Not intended to be called for Ground body')
        KqB      = fBMB(o.BB_inB, o.KK, sympy=o.sympy, name=o.name)
        n        = KqB.shape[0]
        K[:n,:n] = K[:n,:n]+KqB     
        for c in o.Children:
            K=c._getFullK(K)
        return K
        
    def _getFullD(o, D):
        if isinstance(o, YAMSRecGroundBody):
            raise Exception('Not intended to be called for Ground body')
        DqB      = fBMB(o.BB_inB,o.DD, sympy=o.sympy, name=o.name)
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
        return self.eye(3);
    @property
    def Bhat_x_bc(self):
        return self.Matrix(np.zeros((3,0)))
    @property
    def Bhat_t_bc(self):
        return self.Matrix(np.zeros((3,0)))
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
# --- Ground/inertial Body 
# --------------------------------------------------------------------------------{
class YAMSRecGroundBody(YAMSRecBody, GenericInertialBody):
    """ 
    Ground body is used to traverse the tree and hold the full mass matrix
    """
    def __init__(B, sympy=False):
        YAMSRecBody.__init__(B, name='Grd', sympy=sympy)
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
        s='<YAMSRecGroundBody {} object>:\n'.format(self.name)
        s+='|Inherits:\n'
        s+='||'+'\n|'.join(YAMSRecBody.__repr__(self).split('\n'))+'--->\n'
        s+='||'+'\n|'.join(GenericInertialBody.__repr__(self).split('\n'))+'--->\n'
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
        M = o.Matrix(np.zeros((o.nq, o.nq)))
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
# --- YAMSRec Rigid Body 
# --------------------------------------------------------------------------------{
class YAMSRecRigidBody(YAMSRecBody,GenericRigidBody):
    def __init__(B, name, mass, J, rho_G=None, s_OP=None, r_O=None, R_b2g=None, sympy=False):
        """
        Creates a rigid body for YAMSRec

        NOTE:
          - Legacy call was always with J_G and rho_G
          - We are adding s_OP

        """
        if s_OP is not None:
            raise Exception('[INFO] You are using the new interface, it should work, but lets debug it')
        YAMSRecBody.__init__(B, name, sympy=sympy)
        #              Interface:(name, mass, J, s_OG, r_O=[0,0,0], R_b2g=np.eye(3), s_OP=None):
        GenericRigidBody.__init__(B, name=name, mass=mass, J=J, s_OG=rho_G, r_O=r_O, R_b2g=R_b2g, s_OP=s_OP)

        B.s_G_inB = B.masscenter
        B.J_G_inB = B.masscenter_inertia
        B.J_O_inB = translateInertiaMatrixFromCOG(B.J_G_inB, mass, -B.s_G_inB)
        B.MM = buildRigidBodyMassMatrix(mass, B.J_O_inB, B.s_G_inB) # TODO change interface
        B.DD = np.zeros((6,6))
        B.KK = np.zeros((6,6))
        # END - YAMSRec RigidBody

    def __repr__(self):
        s='<YAMSRecRigidBody {} object>:\n'.format(self.name)
        s+='|Inherits:\n'
        s+='|'+'\n||'.join(YAMSRecBody.__repr__(self).split('\n'))+'--->\n'
        s+='|'+'\n||'.join(GenericRigidBody.__repr__(self).split('\n'))+'--->\n'
        s+='|Properties:\n'
        s+='| - MM:\n{}\n'.format(self.MM)
        return s




# --------------------------------------------------------------------------------}
# --- YAMS Recursive Beam Body 
# --------------------------------------------------------------------------------{
class YAMSRecBeamBody(GenericBeamBody, YAMSRecBody): 
    def __init__(B, 
                 s_span=None, s_P0=None, m=None, PhiU=None, PhiV=None, PhiK=None, EI=None, jxxG=None, s_G0=None, 
            s_min=None, s_max=None,
            bAxialCorr=False, bOrth=False, Mtop=0, bStiffening=True, gravity=None,main_axis='z',
            massExpected=None,
            damp_zeta=None,
            name='dummyYAMSRecBeamBody',
            algo='',
            directions=None,
            sympy=False
            ):
        """ 
          Points P0 - Undeformed mean line of the body
        """
        int_method    = 'Flex'
        if algo=='OpenFAST': 
            int_method='OpenFAST'
        # --- Inherit from BeamBody and YAMSRecBody 
        YAMSRecBody.__init__(B, sympy=sympy)

        if sympy:
            if directions is None: 
                raise Exception('directions shouldnt be None with sympy')
            nf = len(directions)
            # --- TODO WE CREATE A FAKE INTERFACE
            B.main_axis=main_axis
            B.directions = directions
            nSpan = 2
            B.s_span = [0,symbols('L')]
            B.PhiU = []
            B.PhiV = []
            B.gzf = B.Matrix([0]*nf)
            for j in range(nf):
                PhiU = B.Matrix(np.zeros((3,nSpan)))
                PhiV = B.Matrix(np.zeros((3,nSpan)))
                direction = B.directions[j]
                nD = len(B.directions[j])
                if 'x' in direction:
                    PhiU[0,-1] = symbols('ux{:d}c'.format(j+1))
                    PhiV[0,-1]=symbols('vy{:d}c'.format(j+1))
                if 'y' in direction:
                    PhiU[1,-1] = symbols('uy{:d}c'.format(j+1))
                    PhiV[1,-1] = symbols('ux{:d}c'.format(j+1))
                B.PhiU.append(PhiU)
                B.PhiV.append(PhiV)
        else:
            GenericBeamBody.__init__(B,name, s_span, s_P0, m, EI, PhiU, PhiV, PhiK, jxxG=jxxG, s_G0=s_G0, s_min=s_min, s_max=s_max,
                 bAxialCorr=bAxialCorr, bOrth=bOrth, Mtop=Mtop, bStiffening=bStiffening, gravity=gravity, main_axis=main_axis,
                 damp_zeta=damp_zeta,
                 massExpected=massExpected,
                 int_method=int_method
                )

        B.gzf   = B.Matrix(np.zeros((B.nf,1)))
        B.gzpf  = B.Matrix(np.zeros((B.nf,1)))
        B.gzppf = B.Matrix(np.zeros((B.nf,1)))

        # TODO
        B.V0         = B.Matrix(np.zeros((3,B.nSpan)))
        B.K0         = B.Matrix(np.zeros((3,B.nSpan)))
        B.rho_G0_inS = B.Matrix(np.zeros((3,B.nSpan))) # location of COG in each cross section
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
        gzf = np.atleast_1d(self.gzf)
        if self.sympy:
            return self.Bhat_t_bc @ gzf
        else:
            return (self.Bhat_t_bc @ gzf).ravel()

    @property
    def R_bc(self):
        if self.sympy:
            # We use analytical couplings
            if self.main_axis=='x':
                alpha_y= symbols('alpha_y') #-p.V(3,iNode);
                alpha_z= symbols('alpha_z') # p.V(2,iNode);
                return R_y(alpha_y) @ R_z(alpha_z)

            elif self.main_axis=='z':
                alpha_x= symbols('alpha_x') #-p.V(2,iNode);
                alpha_y= symbols('alpha_y') # p.V(1,iNode);
                return R_x(alpha_x)*R_y(alpha_y)
            else:
                raise NotImplementedError()
        else:
            alpha = self.alpha_couplings

            if self.main_axis=='x':
                return R_y(alpha[1]) @ R_z(alpha[2])
            elif self.main_axis=='z':
                return R_x(alpha[0]) @ R_y(alpha[1])
            else:
                raise NotImplementedError()

    def updateKinematics(o,x_0,R_b2g,gz,v_0,a_v_0, verbose=False):
        super(YAMSRecBeamBody,o).updateKinematics(x_0, R_b2g, gz, v_0, a_v_0)
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
            if not o.sympy:
                # TODO for sympy
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
                    if verbose:
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
        displ_g = B.R_b2g.dot(B.s_P) # TODO reference line or COG
        pos_g   = B.pos_global + displ_g
        return pos_g

    @property
    def nSpan(B):
        return len(B.s_span)

    def __repr__(self):
        s='<YAMSRec BeamBody {} object>:\n'.format(self.name)
        s+='|Inherits:\n'
        s+='|'+'\n||'.join(YAMSRecBody.__repr__(self).split('\n'))+'\n'
        s+='|'+'\n||'.join(GenericBeamBody.__repr__(self).split('\n'))+'\n'
        s+='|Properties:\n'
        s+='| * nSpan: {}\n'.format(self.nSpan)
        return s


# --------------------------------------------------------------------------------}
# --- Uniform Beam Body 
# --------------------------------------------------------------------------------{
class YAMSRecUniformBeamBody(YAMSRecBeamBody):
    def __init__(B, name, nShapes, nSpan, L, EI0, m, Mtop=0, jxxG=None, GKt=None, 
            bAxialCorr=True, bCompatibility=False, bStiffnessFromGM=False, bStiffening=True, 
            gravity=None, main_axis='x',
            shapeFunctions='masslessbeam',
            bottomBC='clamped',
            topBC='free',
            sympy=False
            ):

        import welib.beams.theory as bt
        if jxxG is None:
            jxxG=0
        if GKt is None:
            GKt=0
        if gravity is None and (bStiffening or Mtop>0):
            raise Exception('`gravity` is None, but `bStiffening` is true or `Mtop`>0. Please provide `gravity`.')

        A=1; rho=A*m;
        x=np.linspace(0,L,nSpan);
        BC = '{}-{}'.format(bottomBC, topBC)
        # Mode shapes
        if shapeFunctions=='masslessbeam':
            freq,s_span,U,V,K = bt.UniformBeamBendingModes('unloaded-topmass-{}'.format(BC),EI0,rho,A,L,x=x,Mtop=Mtop, nModes=nShapes)
        elif shapeFunctions=='Guyan':
            if BC=='clamped-free':
                freq,s_span,U,V,K  = bt.UniformBeamGuyanModes(EI0, rho, A, L, x= x, nModes = nShapes)
            else:
                raise NotImplementedError('{} {}'.format(shapeFunctions, BC))

        elif shapeFunctions=='FEM':
            #from welib.FEM.fem_beam import *
            raise NotImplementedError('{} {}'.format(shapeFunctions, BC))
            #nel      = 10             # Number of elements along the beam
            #element  = 'frame3d'      # Type of element used in FEM
            #if TopMass:
            #    Mtop = 50000  # Top mass [kg]
            #    M_tip= rigidBodyMassMatrixAtP(m=Mtop, J_G=None, Ref2COG=None)
            #else:
            #    M_tip=None
            ## --- Structural data for uniform beam
            #E   = 210e9     # Young modulus [Pa] [N/m^2]
            #G   = 79.3e9    # Shear modulus. Steel: 79.3  [Pa] [N/m^2]
            #L   = 100       # Beam Length [m]
            #EIy0= 1.654e+12 # Planar second moment of area [m^4]
            #m0  = 1.026e+04 # Mass per length [kg/m]
            #EIx0= EIy0*2    # Polar second moment of area [m^4]
            #A   = 1.00      # Area [m^2] 
            #Kt  = EIy0/E*10 # Torsion constant [m^4]
            ## --- Compute FEM model and mode shapes
            #FEM=cbeam(L,m=m0,EIx=EIx0,EIy=EIy0,EIz=EIy0,EA=E*A,A=A,E=E,G=G,Kt=Kt,
            #        element=element, nel=nel, BC=BC, M_tip=M_tip)
            #x =FEM['xNodes'][0,:]
            #Q =FEM ['Q']
            #QY1 = [Q [1::6, iMode] for iMode in range(10) if FEM ['modeNames'][iMode].startswith('uy')]

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
        super(YAMSRecUniformBeamBody,B).__init__(s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=jxxG, bAxialCorr=bAxialCorr, Mtop=Mtop, bStiffening=bStiffening,
                gravity=gravity, main_axis=main_axis, name=name, sympy=sympy)


    def __repr__(self):
        s='<YAMSRec UniformBeamBody {} object>:\n'.format(self.name)
        s+='|Inherits:\n'
        s+='|'+'\n||'.join(YAMSRecBeamBody.__repr__(self).split('\n'))+'\n'
        return s


# --------------------------------------------------------------------------------}
# --- FAST Beam body 
# --------------------------------------------------------------------------------{
class YAMSRecFASTBeamBody(YAMSRecBeamBody, GenericFASTBeamBody):
    def __init__(B, body_type, ED, inp, Mtop=0, shapes=None, nShapes=None, main_axis='x',nSpan=None,bAxialCorr=False,bStiffening=True, 
            spanFrom0=False, massExpected=None, gravity=None,
            algo='', # TODO OpenFAST
            sympy=False
            ):
        """ 
        """
        if nShapes is not None:
            raise Exception('nShapes is depreciated use shapes instead')
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
                gravity=gravity,
                algo=algo
                )
        # We need to inherit from "YAMS" Beam not just generic Beam
        # NOTE: TODO TODO TODO: This will result in "YAMSBeamBody to be called twice...)
        YAMSRecBeamBody.__init__(B, B.s_span, B.s_P0, B.m, B.PhiU, B.PhiV, B.PhiK, B.EI, jxxG=B.jxxG, s_G0=B.s_G0, 
                # NOTE: r_O, r_b2g is lost here
                damp_zeta=B.damp_zeta, # Important otherwise lost
                s_min=B.s_min, s_max=B.s_max,
                bAxialCorr=bAxialCorr, bOrth=B.bOrth, Mtop=Mtop, bStiffening=bStiffening, gravity=B.gravity,main_axis=main_axis,
                massExpected=massExpected,
                algo=algo,
                sympy=sympy
                )

# --------------------------------------------------------------------------------}
# --- B Matrices 
# --------------------------------------------------------------------------------{
def fB_inB(R_EI, B_I, sympy=False):
    """ Transfer a global B_I matrix (body I at point I) into a matrix in it's own coordinate.
    Simply multiply the top part and bottom part of the B matrix by the 3x3 rotation matrix R_EI
    e.g.
         B_N_inN = [R_EN' * B_N(1:3,:);  R_EN' * B_N(4:6,:)];
    """ 
    def MatrixLoc(m):
        if sympy:
            return Matrix(m)
        else:
            return np.asarray(m)

    if len(B_I)==0:
        B_I_inI = MatrixLoc(np.array([]))
    else:
        B_I_inI = MatrixLoc(np.vstack((R_EI.T @ B_I[:3,:], R_EI.T @ B_I[3:,:])))
    return B_I_inI

def fB_aug(B_I_inI, nf_I, nf_Curr=None, nf_Prev=None, sympy=False):
    """
    Augments the B_I_inI matrix, to include nf_I flexible degrees of freedom.
    This returns the full B matrix on the left side of Eq.(11) from [1], 
    based on the Bx and Bt matrices on the right side of this equation
    """
    def MatrixLoc(m):
        if sympy:
            return Matrix(m)
        else:
            return np.asarray(m)    
    if len(B_I_inI)==0:
        if nf_I>0:
             BB_I_inI = np.vstack( (np.zeros((6,nf_I)), np.eye(nf_I))) 
        else:
             BB_I_inI = np.zeros((6,0))
        if sympy:
             BB_I_inI = Matrix( BB_I_inI.astype(int) )
    else:
        if nf_Curr is not None:
            # Case of several flexible bodies connected to one point (i.e. blades)
            nf_After=nf_I-nf_Prev-nf_Curr
            I = np.block( [np.zeros((nf_Curr,nf_Prev)), np.eye(nf_Curr), np.zeros((nf_Curr,nf_After))] )
        else:
            nf_Curr=nf_I
            I=np.eye(nf_I)

        BB_I_inI = np.block([ [B_I_inI, np.zeros((6,nf_I))], [np.zeros((nf_Curr,B_I_inI.shape[1])), I]]);

    return MatrixLoc(BB_I_inI)


def fBMatRecursion(Bp, Bhat_x, Bhat_t, R0p, r_pi, sympy=False):
    """ Recursive formulae for B' and Bhat 
    See discussion after Eq.(12) and (15) from [1]
    """
    def MatrixLoc(m):
        if sympy:
            return Matrix(m)
        else:
            return np.asarray(m)
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

    #r_pi=vec3(r_pi)

    # TODO use Translate here
    Bi = MatrixLoc(np.zeros((6,ni+n_p)))
    for j in range(n_p):
        Bi[:3,j] = Bp[:3,j] + cross(Bp[3:,j],r_pi) # Recursive formula for Bt mentioned after Eq.(15)
        Bi[3:,j] = Bp[3:,j] # Recursive formula for Bx mentioned after Eq.(12)
    if ni>0:
        Bi[:3,n_p:] = R0p @ Bhat_x[:,:] # Recursive formula for Bx mentioned after Eq.(15)
        Bi[3:,n_p:] = R0p @ Bhat_t[:,:] # Recursive formula for Bt mentioned after Eq.(12)
    return Bi

def fBMatTranslate(Bp, r_pi, sympy=False):
    """
    Rigid translation of a B matrix to another point, i.e. transfer the velocities from a point to another: 
      - translational velocity:  v@J = v@I + om@I x r@IJ
      - rotational velocity   : om@J = om@I
    """
    Bi=np.zeros(Bp.shape)
    if Bp.ndim==1:
        raise NotImplementedError

    for j in range(Bp.shape[1]):
        Bi[0:3,j] = Bp[0:3,j] + cross(Bp[3:6,j],r_pi)
        Bi[3:6,j] = Bp[3:6,j]
    return Bi


def fBMB(BB_I_inI, MM, sympy=False, name=''):
    """ Computes the body generalized matrix: B'^t M' B 
    See Eq.(8) of [1] 
    """
    if MM is None:
        raise Exception(f'MM is None for body {name}')
    MM_I = (np.transpose(BB_I_inI) @ MM) @ BB_I_inI
    return MM_I

if __name__=='__main__':
    pass
