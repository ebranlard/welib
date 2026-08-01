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
from welib.tools.strings import INFO

from welib.yams.utils import skew

# --- To ease comparison with sympy version
# from numpy import eye, cross, cos ,sin

import sympy as sp
from sympy import Matrix, symbols

from welib.tools.strings import prettyMat

def pm(M, var=None, **kwargs):
    if isinstance(M, sp.Basic):
        return M
    else:
        return prettyMat(M, var, **kwargs, digits=3)

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
    def __init__(self, Type, RelPoint=None, RelOrientation=None, JointRotations=None, JointTranslations=None, OrientAfter=True, parentNode=None, parentBody=None, sympy=False):
        self.sympy = sympy
        if RelOrientation is None:
            RelOrientation=self.eye(3)
        if RelPoint is None:
            RelPoint = [0,0,0]

        self.Type=Type
        
        self.s_C0_inB  = self.vec3(RelPoint) # CONSTANT - s_BC0 @t=0
        self.s_C_inB   = self.s_C0_inB       #             s_BC  @t
        self.R_ci_0    = RelOrientation      #           
        self.R_ci      = self.R_ci_0         # CONSTANT
        self.OrientAfter= OrientAfter
        # Related to flexible parent
        self.parentNode = parentNode
        self.parentBody = parentBody
        self.s_P0C0_inB = None               # CONSTANT - See SKETCHC1 from parent node to connection point @ t=0 (when C_0_inB /= P)
        self.s_P0_inB   = None               # CONSTANT Stored for completness, not necessary
        # Related to rigid body joint rotations (e.g. dynamic shaft, but potentially yaw or tilt, less common)
        self.I_DOF = None  # Index of joints DOF in global DOF vector
        self.q     = None  # Index of joints DOF in global DOF vector
        self.JointRotations = JointRotations
        self.JointTranslations = JointTranslations

        if self.parentNode is not None:
            if self.parentBody is None:
                raise Exception('ParentBody unknown but parent Node specified')
            self.s_P0_inB = self.vec3(self.parentBody.s_P0[:,self.parentNode])
            self.s_P0C0_inB = self.s_C0_inB - self.s_P0_inB # See SKETCHC1
            if np.sum(np.abs(self.s_P0C0_inB))!=0:
                INFO(f'Connection point offset is P0C0={self.s_P0C0_inB} from {self.parentBody.name}.')

        if self.Type=='Rigid':
            self.nj=0
        elif self.Type in ['SphericalJoint', 'Joint']:
            if self.JointRotations is None:
                raise Exception('JointRotations is required for Joint/SphericalJoint')
            self.nj=len(self.JointRotations)
        elif self.Type=='Free':
            if self.JointTranslations is None:
                self.JointTranslations = ['x', 'y', 'z']
            if self.JointRotations is None:
                self.JointRotations = ['x', 'y', 'z']
            self.nj = len(self.JointTranslations) + len(self.JointRotations)
        else:
            raise NotImplementedError()

    def _axis_unit(j, rot):
        if rot=='x':
            return j.vec3([1,0,0])
        elif rot=='y':
            return j.vec3([0,1,0])
        elif rot=='z':
            return j.vec3([0,0,1])
        raise Exception('Unknown axis {}'.format(rot))

    def _axis_rotation(j, rot, qval):
        if rot=='x':
            return R_x(qval)
        elif rot=='y':
            return R_y(qval)
        elif rot=='z':
            return R_z(qval)
        raise Exception('Unknown axis {}'.format(rot))

    def updateConnectionKinematics(j, q):
        """ Connection/joint updateConnectionKinematics

        Update Rotation between bodies (from true joints)

        NOTE: position of connection point (s_C_inB) is updated by parent kinematics!

        """
        j.B_ci = j.Matrix(np.zeros((6,j.nj)))
        if j.Type=='Rigid':
            j.R_ci=j.R_ci_0

        elif j.Type in ['SphericalJoint', 'Joint', 'Free']:
            R = j.eye(3)
            if j.sympy:
                myq = [q[int(i), 0] for i in np.asarray(j.I_DOF).ravel()]
            else:
                myq = q[j.I_DOF, 0]
            j.q = myq
            iq = 0

            if j.Type=='Free':
                s = j.vec3(j.s_C0_inB)
                for itr,tra in enumerate(j.JointTranslations):
                    I = j._axis_unit(tra)
                    qtr = myq[itr]
                    j.B_ci[:3,iq] = I
                    s = s + qtr * I
                    iq += 1
                j.s_C_inB = s

            for rot in j.JointRotations:
                I = j._axis_unit(rot)
                Rj = j._axis_rotation(rot, myq[iq])
                # Setting Bhat column by column (before R updates)
                j.B_ci[3:,iq] = R @ I
                # Updating rotation matrix
                R      = R @ Rj
                iq += 1

            if j.OrientAfter:
                j.R_ci = j.Matrix(R @ j.R_ci_0)
            else:
                j.R_ci = j.Matrix(j.R_ci_0 @ R)

        else:
            raise NotImplementedError('Joint Type' + j.Type)


    def updateConnectionPointPosition(j, R_pc, silent=False):
        r""" 
        SKETCHC1 for Connection not at Parent Node P
                                       uP0P
         C0                        C0--------->.    C
          \ s_P0C0                  \          \ R / 
           \                         \    uP0P  \^/
            P0           ->           P0-------->P
            |                         |     _/
            |                         |  _ /
            |                         | /
            B                         B

        Two main equivalent formulations ("from P" or "from C0"):
          - s_C = s_P(t)         +   R(t)      @ s_P0C0     ("from P")
          - s_C = s_C0 + uP0P(t) + [ R(t) -I ] @ s_P0C0     ("from C0")

        INPUTS:
          - R_pc: The rotation the the flexible parent induce
          -  silent: 
        """
        if j.parentNode is not None:
            #print('>>>> Joint Kinematics. Updating joint position based on parent node position')
            iNode = j.parentNode
            # How much the parentBody node has translated in parent body
            #uP0P = j.vec3(j.parentBody.s_P[:,iNode])-j.vec3(j.parentBody.s_P0[:,iNode])
            s_P = j.vec3(j.parentBody.s_P[:,iNode])
            #uP0P = s_P-j.s_PP_0_inB  # KEEP me Translation of P
            s_PC = R_pc @ j.s_P0C0_inB
            j.s_C_inB = s_P  + s_PC
            #if not silent:
            #    print('>>> s_PC', s_PC)
            #    print(j)

    def addLinVelJacobianContrib(j, Bhat_x, Bhat_t, R_pc):
        """ Add contribution due to connection offset from body extremity point"""
        if j.s_P0C0_inB is not None:
            s_PC = R_pc @ j.s_P0C0_inB 
            #print('>>>> Bhat_x\n', Bhat_x)
            Bhat_x_c = -skew(s_PC, symb = j.sympy) @ Bhat_t
            #print('>>>> Bhat_x_c\n', Bhat_x_c)
            Bhat_x  += Bhat_x_c
            #print('>>>> Bhat_x\n', Bhat_x)
        return Bhat_x

    def __repr__(self):
        s ='<Connection object>:\n'
        s+='|Properties:\n'
        s+='| - Type:     {} \n'.format(self.Type)
        s+='| - OrientAfter:  {} \n'.format(self.OrientAfter)
        s+='| - s_C0_inB: (init pos.    of C    in body) : {} \n'.format(pm(self.s_C0_inB))
        s+='| - s_C_inB:  (current pos. of C    in body) : {} \n'.format(pm(self.s_C_inB))
        s+='| - R_ci_0:   (init rot. ro conn. in body)   :\n{} \n'.format(pm(self.R_ci_0))
        s+='| - R_ci:     (current rot. ro conn. in body):\n{} \n'.format(pm(self.R_ci))
        s+='|Related to joint rotations (independent of flex. body)\n'
        s+='| - I_DOF:    (DOFs involved in joint       ): {} \n'.format(pm(self.I_DOF))
        s+='| - q_c  :    (DOFs involved in joint       ): {} \n'.format(pm(self.q))
        s+='|Related to parent body (flexbled body)\n'
        s+='| - s_P0C0_inB: (offset from body)           : {} \n'.format(pm(self.s_P0C0_inB))
        s+='| - parentNode index "P"                     : {} \n'.format(self.parentNode)
        s+='| - parentBody name:                         : {} \n'.format(self.parentBody.name)
        s+='|Methods: updateKinematics\n'
        s+='|Usefull getters: None \n'
        return s

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
    def eye(self, n): 
        if self.sympy:
            return Matrix( np.eye(n).astype(int) )
        else:
            return np.eye(n)
    # --- End generic tools



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
        #  B.Bhat_x_bc   = None
        #  B.Bhat_t_bc   = None
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

    def kinematics_export(B):
        """Return a canonical kinematics payload for cross-flavor comparisons.

        This intentionally exposes a stable dictionary schema so the numeric
        recursive and sympy-mechanics flavors can be compared without forcing
        implementation unification yet.
        """
        def _safe_get(attr, default=None):
            return getattr(B, attr) if hasattr(B, attr) else default

        return {
            'name': B.name,
            'flavor': 'yams_rec',
            'sympy': bool(B.sympy),
            'nf': int(B.nf) if hasattr(B, 'nf') else 0,
            'I_DOF': _safe_get('I_DOF', None),
            'pos_global': _safe_get('pos_global', None),
            'R_b2g': _safe_get('R_b2g', None),
            'R_g2b': _safe_get('R_g2b', None),
            'R_bc': _safe_get('R_bc', None),
            'Bhat_x_bc': _safe_get('Bhat_x_bc', None),
            'Bhat_t_bc': _safe_get('Bhat_t_bc', None),
            'B': _safe_get('B', None),
            'B_inB': _safe_get('B_inB', None),
            'BB_inB': _safe_get('BB_inB', None),
        }

    def kinematics_export_tree(B):
        """Return canonical kinematics payload for this body and descendants."""
        return [b.kinematics_export() for b in B.bodies]

    def B_matrix(B, in_body=False):
        """Return kinematic B matrix at the body origin.

        Parameters
        ----------
        in_body : bool
            If True, return body-coordinate matrix (B_inB), otherwise global (B).
        """
        return B.B_inB if in_body else B.B

    def BB_matrix(B):
        """Return augmented body-coordinate BB matrix."""
        return B.BB_inB

    def generalized_mass_matrix(B):
        """Return generalized mass contribution of this body.

        This is the body-level contribution $B'^T M' B'$ (or $BB^T M' BB$ for
        flexible bodies via ``BB_inB``), consistent with recursive assembly.
        """
        if isinstance(B, YAMSRecGroundBody):
            raise Exception('Ground body has no standalone generalized mass contribution')
        return fBMB(B.BB_inB, B.MM, sympy=B.sympy, name=B.name)

    def Bhat_matrix(B, kind='x'):
        """Return connection-point Bhat matrix for flexible coupling.

        Parameters
        ----------
        kind : {'x','t'}
            'x' for translational part, 't' for rotational part.
        """
        if kind == 'x':
            return B.Bhat_x_bc
        if kind == 't':
            return B.Bhat_t_bc
        raise ValueError("kind should be 'x' or 't'")

    def connectTo(self, Child, Point=None, Type=None, BodyPoint=None, RelOrientation=None, JointRotations=None, JointTranslations=None, OrientAfter=True):
        """ 
         - BodyPoint: in FirstPoint or LastPoint 
        """
        if Type in ['SphericalJoint', 'Joint', 'Free']:
            c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations, JointTranslations=JointTranslations, OrientAfter=OrientAfter, sympy=self.sympy)

        elif Type == 'Rigid': 
            i_C_inB    = None
            parentBody = None
            if BodyPoint is not None:
                # We find index of closest parent node
                if BodyPoint == 'FirstPoint':
                    i_C_inB=0;
                elif BodyPoint == 'LastPoint':
                    i_C_inB=self.nSpan-1
                else:
                    raise NotImplementedError()
                if Point is None:
                    Point = self.s_P0[:,i_C_inB]

            c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations, JointTranslations=JointTranslations, OrientAfter=OrientAfter, parentNode=i_C_inB, parentBody=self, sympy=self.sympy)

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

        q = p.Matrix( np.asarray(q).reshape((-1, 1)) )

        # At this stage all the kinematics of the body p are known
        # Useful variables
        R_0p =  p.R_b2g
        B_p  =  p.B
        r_0p  = p.pos_global  # Position of body origin in global coordinates

        nf_all_children=sum([child.nf for child in p.Children])

        for ic,(body_i,conn_pi) in enumerate(zip(p.Children,p.Connections)):
            #print(f'Kinematics connections {p.name} > {body_i.name}')
            # Flexible influence to connection point
            R_pc  = p.R_bc        # Contain influence of alpha couplings
            Bx_pc = p.Bhat_x_bc   # NOTE: missing a contributionwhen we have an offset
            Bt_pc = p.Bhat_t_bc   #
            conn_pi.updateConnectionPointPosition(p.R_bc) # Update position of connection point 
            conn_pi.addLinVelJacobianContrib(Bx_pc, Bt_pc, p.R_bc) # Add contribution if rigid body offset

            # Joint influence to next body (R_ci, B_ci) (not a function of flexible body motion!)
            conn_pi.updateConnectionKinematics(q) # TODO

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
            #print('p.name',p.name)
            #print('r_0p',r_0p)
            #print('r_pi',r_pi)
            #print('r_0i',r_0i)
            body_i.R_pb       = R_pi
            body_i.pos_global = r_0i
            body_i.R_b2g      = R_0i

            # TODO flexible dofs and velocities/acceleration
            if len(body_i.I_DOF)>0:
                # NOTE: won't work for sympy if indexing is empty
                if p.sympy:
                    body_i.gzf = Matrix([q[int(i), 0] for i in np.asarray(body_i.I_DOF).ravel()])
                else:
                    body_i.gzf  = q[body_i.I_DOF,0] # TODO use updateKinematics

#            TODO TODO TODO: remaining from matlab?????
            gzf  = body_i.gzf
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
        """ YAMSRecBody updateKinematics"""
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
        """ Ground/Rigid Body - For flexible look in this file below"""
        return self.eye(3);

    @property
    def Bhat_x_bc(self):
        """ Ground/Rigid Body - For flexible go to bodies.py"""
        return self.Matrix(np.zeros((3,0)))

    @property
    def Bhat_t_bc(self):
        """ Ground/Rigid Body - For flexible go to bodies.py """
        return self.Matrix(np.zeros((3,0)))

    @property
    def nf(B):
        """ Generic"""
        if hasattr(B,'PhiU'):
            return len(B.PhiU)
        else:
            return 0
    @property
    def mass(B):
        """ Generic"""
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
        GenericInertialBody.__init__(B, sympy=sympy)
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
        s+='|- q:  {}\n'.format(pm(self.q))
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

    def system_mass_matrix(o):
        """Return full assembled generalized mass matrix."""
        return o.M
        
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
        if rho_G is None:
            rho_G = (0,0,0)
        YAMSRecBody.__init__(B, name, sympy=sympy)
        #              Interface:(name, mass, J, s_OG, r_O=[0,0,0], R_b2g=np.eye(3), s_OP=None):
        GenericRigidBody.__init__(B, name=name, mass=mass, J=J, s_OG=rho_G, r_O=r_O, R_b2g=R_b2g, s_OP=s_OP, sympy=sympy)

        B.s_G_inB = B.masscenter
        B.J_G_inB = B.masscenter_inertia
        B.J_O_inB = translateInertiaMatrixFromCOG(B.J_G_inB, mass, -B.s_G_inB)
        B.MM = buildRigidBodyMassMatrix(mass, B.J_O_inB, B.s_G_inB, symb=sympy) # TODO change interface
        if sympy:
            B.DD = Matrix(np.zeros((6,6)).astype(int))
            B.KK = Matrix(np.zeros((6,6)).astype(int))
        else:
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
                    # yams_sympy convention: x-deflection contributes alpha_y via v_y*
                    PhiV[1,-1]=symbols('vy{:d}c'.format(j+1))
                if 'y' in direction:
                    PhiU[1,-1] = symbols('uy{:d}c'.format(j+1))
                    # yams_sympy convention: y-deflection contributes alpha_x via v_x*
                    PhiV[0,-1] = symbols('vx{:d}c'.format(j+1))
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
        if self.sympy:
            return self.Bhat_t_bc @ self.gzf
        else:
            gzf = np.atleast_1d(self.gzf)
            return (self.Bhat_t_bc @ gzf).ravel()

    @property
    def Bhat_t_bc(self):
        """Flexible connection rotational coupling Jacobian.

        In sympy mode, align with yams_sympy conventions from
        YAMSFlexibleBody.defineExtremity:
          - direction 'x' -> alpha_y term via v_y*
          - direction 'y' -> alpha_x term via v_x*
        """
        if self.sympy:
            Bhat_t_bc = self.Matrix(np.zeros((3, self.nf)))
            for j, direction in enumerate(self.directions):
                if 'x' in direction:
                    Bhat_t_bc[1, j] = symbols('vy{:d}c'.format(j+1))
                if 'y' in direction:
                    Bhat_t_bc[0, j] = symbols('vx{:d}c'.format(j+1))
                if 'z' in direction:
                    Bhat_t_bc[2, j] = symbols('vz{:d}c'.format(j+1))
            return Bhat_t_bc
        return super(YAMSRecBeamBody, self).Bhat_t_bc

    def R_bc_matrix(self, use_symbolic_alpha=None):
        """Return flexible connection rotation matrix.

        Parameters
        ----------
        use_symbolic_alpha : bool or None
            - True:  use compact symbolic alpha placeholders (legacy sympy-rec form)
            - False: use expanded alpha couplings from Bhat_t_bc @ gzf
            - None:  defaults to True in sympy mode, False otherwise
        """
        if use_symbolic_alpha is None:
            use_symbolic_alpha = bool(self.sympy)

        if use_symbolic_alpha:
            if self.main_axis=='x':
                alpha_y = symbols('alpha_y')
                alpha_z = symbols('alpha_z')
                return R_y(alpha_y) @ R_z(alpha_z)
            elif self.main_axis=='z':
                alpha_x = symbols('alpha_x')
                alpha_y = symbols('alpha_y')
                return R_x(alpha_x) * R_y(alpha_y)
            else:
                raise NotImplementedError()

        alpha = self.alpha_couplings
        if self.main_axis=='x':
            return R_y(alpha[1]) @ R_z(alpha[2])
        elif self.main_axis=='z':
            return R_x(alpha[0]) @ R_y(alpha[1])
        else:
            raise NotImplementedError()

    @property
    def R_bc(self):
        """Flexible-body connection rotation matrix.

        Backward-compatible property equivalent to:
            R_bc_matrix(use_symbolic_alpha=self.sympy)
        """
        return self.R_bc_matrix(use_symbolic_alpha=self.sympy)

    def updateKinematics(o,x_0,R_b2g,gz,v_0,a_v_0, verbose=False):
        """ YAMSRec BeamBody updateKinematics"""
        super(YAMSRecBeamBody,o).updateKinematics(x_0, R_b2g, gz, v_0, a_v_0)
        # --- Calculation of deformations wrt straight beam axis, curvature (K) and velocities (UP)
        #print(f'>>>>>>>>>>>>>>>>>> Update Kin flexible body {o.name} nf={o.nf} nc={len(o.Connections)} sympy={o.sympy}')
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

            # --- Update Connection Point location
            # Position of connection point has changed due to flexible bodymotion
            #R_pc  = o.R_bc
            #Bx_pc = o.Bhat_x_bc
            #Bt_pc = o.Bhat_t_bc
            # We do it for consistency, it is not necessary
            for ic, conn in enumerate(o.Connections):
                conn.updateConnectionPointPosition(o.R_bc, silent=True) 

    @property
    def _positions_global(B): # TODO rename
        displ_g = B.R_b2g.dot(B.s_P) # TODO reference line or COG
        r_O = B.pos_global
        pos_g = displ_g
        pos_g[0,:] += r_O[0]
        pos_g[1,:] += r_O[1]
        pos_g[2,:] += r_O[2]
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
            bAxialCorr=True, bStiffnessFromGM=False, bStiffening=True, 
            gravity=None, main_axis='x',
            shapeFunctions='masslessbeam',
            bottomBC='clamped',
            topBC='free',
            MtopInertia=None,
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
        # Optional extra tip-mass inertia contribution handled in BeamBody.computeMassMatrix.
        # Default is no extra contribution so top mass can be represented by an attached rigid body.
        if MtopInertia is None:
            MtopInertia = 0.0
        B.MtopInertia = MtopInertia
        #print('>>> BC', BC)
        # Mode shapes
        if shapeFunctions=='masslessbeam':

            freq,s_span,U,V,K = bt.UniformBeamBendingModes('unloaded-topmass-{}'.format(BC),EI0,rho,A,L,x=x,Mtop=Mtop, nModes=nShapes)

        elif shapeFunctions=='admissible':
            # For internal sub-beam interfaces, enforcing free-end natural BC can be
            # overly restrictive. This Ritz basis enforces only essential BC.
            if bottomBC!='clamped':
                raise NotImplementedError('admissible basis currently supports bottomBC=clamped only')

            if BC=='clamped-free':
                # Keep analytical clamped-free basis for true free-end terminal segment.
                freq,s_span,U,V,K = bt.UniformBeamBendingModes('unloaded-topmass-clamped-free',EI0,rho,A,L,x=x,Mtop=Mtop, nModes=nShapes)
            else:
                freq,s_span,U,V,K = bt.UniformBeamRitzShapeFunctions(BC, L, x=x, nModes=nShapes, norm='tip')

        elif shapeFunctions=='Guyan':
            if BC=='clamped-free':
                if nShapes>2:
                    raise Exception('Guyan only valid for 2 shapes')
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

        GenericFASTBeamBody.__init__(B, ED=ED, inp=inp, Mtop=Mtop, shapes=shapes, main_axis=main_axis, nSpan=nSpan, bAxialCorr=bAxialCorr, bStiffening=bStiffening, 
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
                sympy=sympy, name=body_type
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
        Bi[0:3,j] = Bp[0:3,j] + cross(Bp[3:6,j],r_pi) # Recursive formula for Bt mentioned after Eq.(15)
        Bi[3:6,j] = Bp[3:6,j] # Recursive formula for Bx mentioned after Eq.(12)
    if ni>0:
        Bi[0:3,n_p:] = R0p @ Bhat_x[:,:] # Recursive formula for Bx mentioned after Eq.(15)
        Bi[3:6,n_p:] = R0p @ Bhat_t[:,:] # Recursive formula for Bt mentioned after Eq.(12)
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
