"""

Yams recursive formulation with sympy.
For Kane's formalism (similar), see yams_kane.py

NOTE: this file is intended to be comparable with yams.py

Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""
import numpy as np
import sympy
import sympy as sp
from sympy import Symbol, symbols
from sympy import Matrix, Function, diff
from sympy.printing import lambdarepr
from sympy import init_printing
from sympy import lambdify
#from sympy.abc import *
from sympy import trigsimp
from sympy import cos,sin
from sympy import zeros, transpose

# from sympy.physics.mechanics import Body as SympyBody
from sympy.physics.mechanics import RigidBody as SympyRigidBody
from sympy.physics.mechanics import Point, ReferenceFrame, inertia, dynamicsymbols
from sympy.physics.mechanics.functions import msubs

from sympy.physics.vector import init_vprinting, vlatex

# Local
from welib.yams.yams_sympy_tools import exprHasFunction, skew, colvec, cross #,ete
from collections import OrderedDict 

#init_vprinting(use_latex='mathjax', pretty_print=False)
#
#display=lambda x: sympy.pprint(x, use_unicode=False,wrap_line=False)


__all__ = ['YAMSBody','YAMSInertialBody','YAMSRigidBody','YAMSFlexibleBody'] # New general implementation
__all__+= ['YAMSRecSPBody','RigidBody','GroundBody'] # Old "recursive" implementation. TODO merge the two
__all__+= ['skew', 'rotToDCM', 'DCMtoOmega']

# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def ensureMat(x, nr, nc):
    """ Ensures that the input is a matrix of shape nr, nc"""
    if not isinstance(x,Matrix):
        x=Matrix(x)
    return x.reshape(nr, nc)

def ensureList(x, nr):
    """ Ensures that the input is a list of length nr"""
    x = list(x)
    if len(x)!=nr:
        raise Exception('Wrong dimension, got {}, expected {}'.format(len(x),nr))
    return x
            
def coord2vec(M31, e):
    """ Ugly conversion from a matrix or vector coordinates (implicit frame) to a vector (in a given frame) """
    M31 = ensureList(M31, 3)
    return M31[0] * e.x + M31[1] * e.y + M31[2] * e.z


def rotToDCM(rot_type, rot_amounts, rot_order=None):
    """
    return matrix from a ref_frame to another frame rotated by specified amounts

    INPUTS (see sympy.physics.vector.ReferenceFrame.orient):
     - rot_type    : The method used to generate the direction cosine matrix. Supported methods are
              'SmallRot': small angle rotations (ADDED in YAMS)
              'Axis': simple rotations about a single common axis
              'DCM': for setting the direction cosine matrix directly
              'Body': three successive rotations about new intermediate axes, also called "Euler and Tait-Bryan angles"
              'Space': three successive rotations about the parent frames' unit vectors
              'Quaternion': rotations defined by four parameters which result in a singularity free direction cosine matrix
     - rot_amounts : expressions defining the rotation angles or direction cosine matrix. These must match the rot_type. 
                The input types are:
                'Axis': 2-tuple (expr/sym/func, Vector)
                'DCM': Matrix, shape(3,3)
                'Body': 3-tuple of expressions, symbols, or functions
                'Space': 3-tuple of expressions, symbols, or functions
                'Quaternion': 4-tuple of expressions, symbols, or functions

     - rot_order: string or int. If applicable, the order of the successive of rotations. 
                  The string '123', 'XYZ' and integer 123 are equivalent.
                  Required for 'Body' and 'Space'.
    
    New type added: SmallRot
    
    see sympy.orientnew
        rotToDCM('Axis', (3, N.x)       )
        rotToDCM('Body', (x,y,z), 'XYZ' )
        rotToDCM('DCM' , M )
        rotToDCM('SmallRot' , (x,y,z) )
    """
    ref_frame = ReferenceFrame('dummyref')
    if rot_type =='SmallRot':
            M=-skew(rot_amounts)
            M[0,0]=1
            M[1,1]=1
            M[2,2]=1
            return M
    elif rot_type in ['Body','Space']:
        frame = ref_frame.orientnew('dummy', rot_type, rot_amounts, rot_order)
    else:
        frame = ref_frame.orientnew('dummy', rot_type, rot_amounts)
    return frame.dcm(ref_frame) # from parent to frame
    
def DCMtoOmega(DCM, ref_frame=None):
    """
    Given a DCM matrix, returns the rotational velocity: omega = R' * R^t

        DCM = ref_frame.dcm(body.frame) # from body to inertial
           -> Omega of body wrt ref_frame expressed in ref frame

    If ref_frame is None, only the coordinates of omega are given 
    otherwise, a vector is returned, expressed in ref_frame

    """
    t = dynamicsymbols._t
    OmSkew = (DCM.diff(t) *  DCM.transpose()).simplify()
    if ref_frame is None:
        return (OmSkew[2,1], OmSkew[0,2], OmSkew[1,0])
    else:
        return OmSkew[2,1] * ref_frame.x + OmSkew[0,2]*ref_frame.y + OmSkew[1,0] * ref_frame.z



# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def   SmpMat(bodyname, varname, nr, nc, nq, rname=None, cname=None, noZeroExp=False, singleDOFNumbering=True):
        if rname is None:
            if singleDOFNumbering or nr>1:
                rname=list(np.arange(nr)+1)
            else:
                rname=['']
        if cname is None:
            if singleDOFNumbering or nc>1:
                cname=list(np.arange(nc)+1)
            else:
                cname=['']
        if len(cname)!=nc:
            raise Exception('cname length should match nc for Taylor {} {}'.format(bodyname, varname))
        if len(rname)!=nr:
            raise Exception('rname length should match nr for Taylor {} {}'.format(bodyname, varname))
            
        M0=Matrix(np.zeros((nr,nc)).astype(int))
        for i in np.arange(nr):
            for j in np.arange(nc):
                if noZeroExp:
                    M0[i,j] = symbols('{}_{}_{}{}'.format(varname,bodyname,rname[i],cname[j])) 
                else:
                    M0[i,j] = symbols('{}^0_{}_{}{}'.format(varname,bodyname,rname[i],cname[j])) 
        return M0


class Taylor(object):
    r""" 
    A Taylor object contains a Taylor expansion of a variable as function of q
        M = M^0 + \sum_j=1^nq M^1_j q_j
    where M, M^0, M^1_j are matrices of dimension nr x nc
    See Wallrapp 1993/1994
    """
    def __init__(self, bodyname, varname, nr, nc, nq, rname=None, cname=None, q=None, order=2, noZeroExp=False, singleDOFNumbering=True):
        if rname is None:
            if singleDOFNumbering or nr>1:
                rname=list(np.arange(nr)+1)
            else:
                rname=['']
        if cname is None:
            if singleDOFNumbering or nc>1:
                cname=list(np.arange(nc)+1)
            else:
                cname=['']
        if len(cname)!=nc:
            raise Exception('cname length should match nc for Taylor {} {}'.format(bodyname, varname))
        if len(rname)!=nr:
            raise Exception('rname length should match nr for Taylor {} {}'.format(bodyname, varname))
            
        self.varname=varname
        self.bodyname=bodyname
        self.nq=nq
            
        self.M0=Matrix(np.zeros((nr,nc)).astype(int))
        for i in np.arange(nr):
            for j in np.arange(nc):
                if noZeroExp:
                    self.M0[i,j] = symbols('{}_{}_{}{}'.format(varname,bodyname,rname[i],cname[j])) 
                else:
                    self.M0[i,j] = symbols('{}^0_{}_{}{}'.format(varname,bodyname,rname[i],cname[j])) 
                
        if order==2: 
            self.M1=[]
            for k in np.arange(nq): # DOF number
                self.M1.append(Matrix(np.zeros((nr,nc)).astype(int)))
                for i in np.arange(nr):
                    for j in np.arange(nc):
                        self.M1[k][i,j] = symbols('{}^1_{}_{}_{}{}'.format(varname,k+1,bodyname,rname[i],cname[j])) 
        if order>2:
            raise NotImplementedError('Order 3 not implemented')
    def eval(self, q=None, order=None):
        """ evaluate the taylor series """
        if q is None:
            q=[symbols('q_{}'.format(j+1)) for j in np.arange(self.nq)]

        M = self.M0
        if hasattr(self,'M1'): 
            nq = len(self.M1)
            if len(q)!=nq:
                raise Exception('Inconsistent dimension between q ({}) and M1 ({}) for Taylor {} {}'.format(len(q),nq,self.bodyname, self.varname))
            for k in np.arange(nq):
                M +=self.M1[k]*q[k]
        return M

    def get(self, dof=0, order=0):
        """ return a given order """
        if order ==0:
            return self.M0
        elif order ==1:
            return self.M1[dof]
        else:
            raise NotImplementedError()
    
    def setOrder(self,order):
        if order==1:
            if hasattr(self,'M1'): 
                del self.M1
        else:
            raise Exception('set order for now mainly removes the 2nd order term')

    def __repr__(self):
        print('Note: just call object.eval() to see the terms M0 and M1 with a generic q instead of __repr__')
        M= self.eval()
        return str(M)


#Me = Taylor('T','Me', 3, 3, nq=2, rname='xyz', cname='xyz')
#Me.M1
#Me.eval([x,y])
#Md = Taylor('T','M_d', 3, 1, nq=2, rname='xyz', cname=[''])
#skew(Md.M0)
# --------------------------------------------------------------------------------}
# --- Connections 
# --------------------------------------------------------------------------------{
class Connection():
    def __init__(self, Type, RelPoint=None, RelOrientation=None, JointRotations=None, OrientAfter=True, parentNode=None, parentBody=None, sympy=True):
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


    def updateKinematics(j,q):
        j.B_ci = j.Matrix(np.zeros((6,j.nj)))
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
class YAMSBody(object):
    def __init__(self, name):
        """
           Origin point have no velocities in the body frame! 
        """
        self.frame     = ReferenceFrame('e_'+name)
        self.origin    = Point('O_'+name)
        self.masscenter= Point('G_'+name)
        self.name=name
        self.origin.set_vel(self.frame,0*self.frame.x)
        
        self.parent = None # Parent body, assuming a tree structure
        self.children = [] # children bodies
        self.inertial_frame = None # storing the typical inertial frame use for computation
        self.inertial_origin = None # storing the typical inertial frame use for computation

        self.viz_opts={}

    @property
    def name(self):
        """The name of the body."""
        return self._name

    @name.setter
    def name(self,name):
        """The name of the body."""
        self._name=name

    def __repr__(self):
        s='<{} object "{}" with attributes:>\n'.format(type(self).__name__,self.name)
        s+=' - origin:       {}\n'.format(self.origin)
        s+=' - frame:        {}\n'.format(self.frame)
        s+=' - inertial_frame: {} \n'.format(self.inertial_frame)
        try:
            posI=self.pos_global
        except:
            posI = 'unknown'
        try:
            posP=self.pos_parent
        except:
            posP = 'unknown'
        try:
            omI=self.omega_inertial
        except:
            omI = 'unknown'
        try:
            omP=self.omega_parent
        except:
            omP = 'unknown'
        try:
            velI=self.vel_inertial
        except:
            velI = 'unknown'
        try:
            velP=self.vel_parent
        except:
            velP = 'unknown'
        try:
            accI=self.acc_inertial
        except:
            accI = 'unknown'
        try:
            rotaccI=self.ang_acc_inertial
        except:
            rotaccI = 'unknown'

            # --- Step 3/4: Velocities and accelerations
            #if not isinstance(body, Particle):
            #    alpha = omega.diff(t, N)
            #OmSkew = (R.diff(t) *  R.transpose()).simplify()
            #omega_ident = OmSkew[2,1] * N.x + OmSkew[0,2]*N.y + OmSkew[1,0] * N.z

        s+=' * pos_parent:     {} (origin wrt to parent)\n'.format(posP)
        s+=' * pos_global:     {} (origin)\n'.format(posI)
        s+=' * vel_inertial:   {} (origin)\n'.format(velI)
        s+=' * vel_parent:     {} (origin)\n'.format(velP)
        s+=' * omega_inertial: {} \n'.format(omI)
        s+=' * omega_parent:   {} \n'.format(omP)
        s+=' * acc_inertial:   {} (origin)\n'.format(accI)
        s+=' * ang_acc_inertial:{} (origin)\n'.format(rotaccI)
        s+=' * R_b2g, R_g2b, _alt\n'
        return s

    def __str__(self):
        return self.__repr__()

    # --------------------------------------------------------------------------------}
    # --- Useful getters
    # --------------------------------------------------------------------------------{
    @property
    def pos_global(self):
        """ return position of body origin in compared to inertial origin"""
        return self.origin.pos_from(self.inertial_origin)

    @property
    def pos_parent(self):
        """ return position of body origin with respect to parent origin"""
        return self.origin.pos_from(self.parent.origin)

    @property
    def vel_inertial(self):
        """ return velocity of origin in inertial frame """
        return self.origin.vel(self.inertial_frame)
    @property
    def vel_inertial_alt(self):
        """ return velocity of origin in inertial frame, alternative formulation dr/dt """
        return  self.pos_global.diff(dynamicsymbols._t, self.inertial_frame) # dr/dt

    @property
    def vel_parent(self):
        """ return velocity of origin in parent frame """
        return self.origin.vel(self.parent.frame)

    @property
    def omega_inertial(self):
        """ return rotational velocity of body in inertial frame """
        # NOTE: in kane omega = zero_uaux(body.frame.ang_vel_in(N))
        return self.ang_vel_in(self.inertial_frame)
    @property
    def omega_inertial_alt(self):
        """ return rotational velocity of body in inertial frame, alternative formulation """
        R      = self.R_b2g
        N      = self.inertial_frame
        OmSkew = R.diff(dynamicsymbols._t) *  R.transpose() # dR/dt * R^t
        return OmSkew[2,1] * N.x + OmSkew[0,2]*N.y + OmSkew[1,0] * N.z

    @property
    def omega_parent(self):
        """ return rotational velocity of body wrt parent frame """
        return self.ang_vel_in(self.parent.frame)
    
    @property
    def acc_inertial(self):
        """ return acceleration velocity of body in inertial frame """
        # NOTE: in kane: acc = zero_udot_uaux(P.acc(N))
        return self.origin.acc(self.inertial_frame)
    @property
    def acc_inertial_alt(self):
        """ return acceleration velocity of body in inertial frame, alternative formulation, acc= dr^2/dt^2 """
        return self.pos_global.diff(dynamicsymbols._t, self.inertial_frame).diff(dynamicsymbols._t,self.inertial_frame)

    @property
    def ang_acc_inertial(self):
        """ return acceleration velocity of body in inertial frame """
        return self.ang_acc_in(self.inertial_frame)


    @property
    def R_b2g(self):
        """ Transformation matrix from body to global """
        return Matrix(self.inertial_frame.dcm(self.frame)) # from body to global

    @property
    def R_g2b(self):
        """ Transformation matrix from global to body """
        return self.R_b2g.transpose() 

    def ang_vel_in(self,frame_or_body):
        """ Angular velocity of body wrt to another frame or body
        This is just a wrapper for the ReferenceFrame ang_vel_in function
        """
        if isinstance(frame_or_body,ReferenceFrame):
            return self.frame.ang_vel_in(frame_or_body)
        else:
            if issubclass(type(frame_or_body),YAMSBody):
                return self.frame.ang_vel_in(frame_or_body.frame)
            else:
                raise Exception('Unknown class type, use ReferenceFrame of YAMSBody as argument')

    def ang_acc_in(self,frame_or_body):
        """ Angular acceleration of body wrt to another frame or body
        This is just a wrapper for the ReferenceFrame ang_acc_in function
        """
        if isinstance(frame_or_body,ReferenceFrame):
            return self.frame.ang_acc_in(frame_or_body)
        else:
            if issubclass(type(frame_or_body),YAMSBody):
                return self.frame.ang_acc_in(frame_or_body.frame)
            else:
                raise Exception('Unknown class type, use ReferenceFrame of YAMSBody as argument')
                
    # --------------------------------------------------------------------------------}
    # --- Connection between bodies
    # --------------------------------------------------------------------------------{
    def connectTo(parent, child, type='Rigid', rel_pos=None, rel_pos_b=None, rot_type='Body', rot_amounts=None, rot_order=None, dynamicAllowed=False, ref_frame=None, verbose=False):
        """ 
        Define a connection between parent and child bodies

        INPUTS:
         - type: 
             'Rigid': the two bodies are rigidly connected
             'Free' : the two bodies have free motions (e.g. a floating platform, or two bodies with spring)
             'Joint': the two bodies have a rotational joint between them, rot_amounts needs to be provided

         - rel_pos: array like of shape 1, 2 or 3, for x, y, z component of relative position between 
                    parent and child expressed in the PARENT COORDINATE SYSTEM. Examples:
                        rel_pos = (x_PC(t), y_PC, 0)   where x_PC is a dynamic symbol

         - rel_pos_b: array like of shape 1, 2 or 3, for x, y, z component of relative position between 
                    parent and child expressed in the CHILD COORDINATE  SYSTEM. Examples:
                        rel_pos = (x_PC(t), y_PC, 0)   where x_PC is a dynamic symbol


        ROTATION/ORIENTATION INPUTS: (see sympy.physics.vector.ReferenceFrame.orient):
         - rot_type    : The method used to generate the direction cosine matrix. Supported methods are
                  'SmallRot': small angle rotations (ADDED in YAMS)
                  'Axis': simple rotations about a single common axis
                  'DCM': for setting the direction cosine matrix directly
                  'Body': three successive rotations about new intermediate axes, also called "Euler and Tait-Bryan angles"
                  'Space': three successive rotations about the parent frames' unit vectors
                  'Quaternion': rotations defined by four parameters which result in a singularity free direction cosine matrix
         - rot_amounts : expressions defining the rotation angles or direction cosine matrix. These must match the rot_type. 
                    The input types are:
                    'Axis': 2-tuple (expr/sym/func, Vector)
                    'DCM': Matrix, shape(3,3)
                    'Body': 3-tuple of expressions, symbols, or functions
                    'Space': 3-tuple of expressions, symbols, or functions
                    'Quaternion': 4-tuple of expressions, symbols, or functions
    
         - rot_order: string or int. If applicable, the order of the successive of rotations. 
                      The string '123', 'XYZ' and integer 123 are equivalent.
                      Required for 'Body' and 'Space'.

         - frame_ref: when provided, the rotations are actually about the ref_frame
                      otherwise, rotations are about the parent frame

        """
        # --- Safety checks
        if rel_pos is None and rel_pos_b is None:
            raise Exception('Provide either rel_pos or rel_pos_b')
        if rel_pos is not None and rel_pos_b is not None:
            raise Exception('Provide either rel_pos or rel_pos_b')

        # --- register parent/child relationship
        child.parent = parent
        parent.children.append(child)
        if isinstance(parent, YAMSInertialBody):
            parent.inertial_frame  = parent.frame
            child.inertial_frame   = parent.frame
            parent.inertial_origin = parent.origin
            child.inertial_origin  = parent.origin
        else:
            if parent.inertial_frame is None:
                raise Exception('Parent body was not connected to an inertial frame. Bodies needs to be connected in order, starting from inertial frame.')
            else:
                child.inertial_frame  = parent.inertial_frame # the same frame is used for all connected bodies
                child.inertial_origin = parent.inertial_origin


        # --- Deal with orientation first (creating a path connecting frames together)
        if ref_frame is None:
            ref_frame = parent.frame
        if rot_amounts is None:
            child.frame.orient(ref_frame, 'Axis', (0, ref_frame.x) ) 
        else:
            if rot_type in ['Body','Space']:
                child.frame.orient(ref_frame, rot_type, rot_amounts, rot_order) # <<< 
            else: 
                # rot_type is DCM
                child.frame.orient(ref_frame, rot_type, rot_amounts) # <<< 


        # --- Then deal with position and velocity
        t = dynamicsymbols._t
        if rel_pos is not None:
            pos_frame = parent.frame
            rel_pos_in_child_frame = False
        else:
            rel_pos = rel_pos_b
            pos_frame = child.frame
            rel_pos_in_child_frame = True
        if len(rel_pos)!=3:
            raise Exception('rel_pos needs to be an array of size 3')

        pos = 0 * pos_frame.x
        vel = 0 * pos_frame.x

        if type=='Free':
            # --- "Free", "floating" connection
            if not isinstance(parent, YAMSInertialBody):
                raise Exception('Parent needs to be inertial body for a free connection')
            # Defining relative position and velocity of child wrt parent
            for d,e in zip(rel_pos[0:3], (pos_frame.x, pos_frame.y, pos_frame.z)):
                if d is not None:
                    pos += d * e
                    vel += diff(d,t) * e

        elif type=='Rigid':
            # Defining relative position and velocity of child wrt parent
            for d,e in zip(rel_pos[0:3], (pos_frame.x, pos_frame.y, pos_frame.z)):
                if d is not None:
                    pos += d * e
                    d_ = diff(d,t)
                    #if exprHasFunction(d) and not dynamicAllowed:
                    if d_!=0 and not dynamicAllowed:
                        raise Exception('Position variable cannot be a dynamic variable for a rigid connection: variable {}'.format(d))
                    if dynamicAllowed:
                        vel += d_ * e
        elif type=='Joint':
            # Defining relative position and velocity of child wrt parent
            for d,e in zip(rel_pos[0:3], (pos_frame.x, pos_frame.y, pos_frame.z)):
                if d is not None:
                    pos += d * e
                    if exprHasFunction(d) and not dynamicAllowed:
                        raise Exception('Position variable cannot be a dynamic variable for a joint connection, variable: {}'.format(d))
                    if dynamicAllowed:
                        vel += diff(d,t) * e
            #  Orientation
            if rot_amounts is None:
                raise Exception('rot_amounts needs to be provided with Joint connection')
            if rot_type.lower()=='axis':
                rot_vars = [rot_amounts[0]]
            else:
                rot_vars = rot_amounts
            for d in rot_vars:
                if d!=0 and d!=1 and not exprHasFunction(d):
                    if verbose:
                        print('>>> WARNING: Rotation amount variable is not a dynamic variable for joint connection, variable: {}'.format(d))
                    #raise Exception('Rotation amount variable should be a dynamic variable for a joint connection, variable: {}'.format(d))
        else:
            raise Exception('Unsupported joint type: {}'.format(type))

        if rel_pos_in_child_frame:
            # NOTE: because the position is expressed wrt child frame, diff above does not work
            vel =  diff(pos, t, parent.frame)
        # Position of child origin wrt parent origin
        child.origin.set_pos(parent.origin, pos)
        # Velocity of child origin frame wrt parent frame (0 for rigid or joint)
        child.origin.set_vel(parent.frame, vel);

        parent.update_kinematics_trigger()

    def update_kinematics_trigger(parent):
        """ """
        for child in parent.children:
            # Velocity of child masscenter wrt parent frame, based on origin vel (NOTE: for rigid body only, should be overriden for flexible body)
            child.masscenter.v2pt_theory(child.origin, parent.frame, child.frame);
            # Velocity of child origin wrt inertial frame, using parent origin/frame as intermediate
            child.origin.v1pt_theory(parent.origin, child.inertial_frame, parent.frame)
            # Velocity of child masscenter wrt inertial frame, using parent origin/frame as intermediate
            child.masscenter.v1pt_theory(parent.origin, child.inertial_frame, parent.frame)

        #r_OB = child.origin.pos_from(child.inertial_origin)
        #vel_OB = r_OB.diff(t, child.inertial_frame)

    def update_origin_pos(body, PointRef, pos):
        """ 
        Set position of origin and perform necessary triggers
        TODO
        """
        pass
        #body.origin.set_pos(PointRef, pos)
        #body.origin.v2pt_theory(ref.origin, ref.frame, body.frame)

    def update_ang_vel(body, omega_symbs, ref_frame=None, out_frame=None):
        # Default values
        if ref_frame is None:
            ref_frame = body.inertial_frame
            #ref_frame = body.parent.frame # <<<
        if out_frame is None:
            out_frame = body.frame

        omega = omega_symbs[0]*out_frame.x + omega_symbs[1]*out_frame.y + omega_symbs[2]*out_frame.z

        body.frame.set_ang_vel(ref_frame, omega)
        # Update origin velocity (origin may move compared to parent origin, so use v1pt)
        body.origin.v1pt_theory    (body.parent.origin, ref_frame, body.frame)   # TODO
        body.origin.v1pt_theory    (body.inertial_origin, ref_frame, body.frame) # TODO
        # Update masscenter velocity (COG and origin fixed in bod frame, so use v2pt)
        body.masscenter.v2pt_theory(body.origin, ref_frame,  body.frame)

    def kdeqsSubsOmega(body, omega_symbs, ref_frame=None, out_frame=None, simplify=True):
        """
        return kinematic differential substitutions for "omega"

        omega is the rotational velocity from the ref_frame to the body frame

        The three symbols omega_symbs are expected to be the components of omega
        in the out_frame

        """
        # Default values
        if ref_frame is None:
            ref_frame = body.inertial_frame
        if out_frame is None:
            out_frame = body.frame
        omega = body.frame.ang_vel_in(ref_frame) 
        kdeqsSubs = []
        if simplify:
            kdeqsSubs += [(omega_symbs[0], omega.dot(out_frame.x).simplify())]
            kdeqsSubs += [(omega_symbs[1], omega.dot(out_frame.y).simplify())] 
            kdeqsSubs += [(omega_symbs[2], omega.dot(out_frame.z).simplify())]
        else:
            kdeqsSubs += [(omega_symbs[0], omega.dot(out_frame.x))]
            kdeqsSubs += [(omega_symbs[1], omega.dot(out_frame.y))] 
            kdeqsSubs += [(omega_symbs[2], omega.dot(out_frame.z))]
        return kdeqsSubs


    def printPosVel(body, origin=True, masscenter=True, inertial=True, kdeqsSubs=None,):
        """ print position and velocity """
        from IPython.display import display
        # Origin position and velocity
        if origin:
            print('Position of body origin from parent origin')
            display(body.origin.pos_from(body.parent.origin))
            if inertial:
                print('Position of body origin from inertial origin')
                display(body.origin.pos_from(body.inertial_origin))
            print('Velocity of body origin from parent frame')
            display(body.origin.vel(body.parent.frame))
            if inertial:
                print('Velocity of body origin from inertial frame')
                display(body.origin.vel(body.inertial_frame))
                if kdeqsSubs:
                    display(body.origin.vel(body.inertial_frame).express(body.inertial_frame).subs(kdeqsSubs).simplify())

        # Mass center position and velocity
        if masscenter:
            print('Position of body masscenter from parent origin')
            display(body.masscenter.pos_from(body.parent.origin))
            if inertial:
                print('Position of body masscenter from inertial origin')
                display(body.masscenter.pos_from(body.inertial_origin))
            print('Velocity of body masscenter from parent frame')
            display(body.masscenter.vel(body.parent.frame))
            if inertial:
                print('Velocity of body masscenter from inertial frame')
                display(body.masscenter.vel(body.inertial_frame))
        # Angular velocity
        print('Angular velocity of body frame wrt parent frame')
        display(body.frame.ang_vel_in(body.parent.frame))
        if kdeqsSubs is not None:
            print('Angular velocity of body frame wrt parent frame (with kdeqs substitutions)')
            display(body.frame.ang_vel_in(body.parent.frame).subs(kdeqsSubs))
        if inertial:
            print('Angular velocity of body frame wrt inertial frame')
            display(body.frame.ang_vel_in(body.inertial_frame))
            if kdeqsSubs is not None:
                print('Angular velocity of body frame wrt inertial frame (with kdeqs substitutions)')
                display(body.frame.ang_vel_in(body.inertial_frame).subs(kdeqsSubs))


    # --------------------------------------------------------------------------------}
    # --- Visualization 
    # --------------------------------------------------------------------------------{
    
    def vizOrigin(self, radius=1.0, color='black', format='pydy'):
        if format=='pydy':
            from pydy.viz.shapes import Sphere
            from pydy.viz.visualization_frame import VisualizationFrame
            return VisualizationFrame(self.frame, self.origin, Sphere(color=color, radius=radius))

    def vizCOG(self, radius=1.0, color='red', format='pydy'):
        if format=='pydy':
            from pydy.viz.shapes import Sphere
            from pydy.viz.visualization_frame import VisualizationFrame
            return VisualizationFrame(self.frame, self.masscenter, Sphere(color=color, radius=radius))

    def vizFrame(self, radius=0.1, length=1.0, format='pydy'):
        if format=='pydy':
            from pydy.viz.shapes import Cylinder
            from pydy.viz.visualization_frame import VisualizationFrame
            from sympy.physics.mechanics import Point
            X_frame  = self.frame.orientnew('ffx', 'Axis', (-np.pi/2, self.frame.z) ) # Make y be x
            Z_frame  = self.frame.orientnew('ffz', 'Axis', (+np.pi/2, self.frame.x) ) # Make y be z
            X_shape   = Cylinder(radius=radius, length=length, color='red') # Cylinder are along y
            Y_shape   = Cylinder(radius=radius, length=length, color='green')
            Z_shape   = Cylinder(radius=radius, length=length, color='blue')
            X_center=Point('X'); X_center.set_pos(self.origin, length/2 * X_frame.y)
            Y_center=Point('Y'); Y_center.set_pos(self.origin, length/2 * self.frame.y)
            Z_center=Point('Z'); Z_center.set_pos(self.origin, length/2 * Z_frame.y)
            X_viz_frame = VisualizationFrame(X_frame, X_center, X_shape)
            Y_viz_frame = VisualizationFrame(self.frame, Y_center, Y_shape)
            Z_viz_frame = VisualizationFrame(Z_frame, Z_center, Z_shape)
        return X_viz_frame, Y_viz_frame, Z_viz_frame

    def vizAsCylinder(self, radius, length, axis='z', color='blue', offset=0, format='pydy'):
        """ """
        if format=='pydy':
            # pydy cylinder is along y and centered at the middle of the cylinder
            from pydy.viz.shapes import Cylinder
            from pydy.viz.visualization_frame import VisualizationFrame
            if axis=='y':
                e = self.frame
                a = self.frame.y
            elif axis=='z':
                e = self.frame.orientnew('CF_'+self.name, 'Axis', (np.pi/2, self.frame.x) ) 
                a = self.frame.z
            elif axis=='x':
                e = self.frame.orientnew('CF_'+self.name, 'Axis', (np.pi/2, self.frame.z) ) 
                a = self.frame.x

            shape = Cylinder(radius=radius, length=length, color=color)
            center=Point('CC_'+self.name); center.set_pos(self.origin, (length/2 +offset) * a)
            return VisualizationFrame(e, center, shape)
        else:
            raise NotImplementedError()


    def vizAsRotor(self, radius=0.1, length=1, nB=3,  axis='x', color='white', format='pydy'):
        # --- Bodies visualization
        if format=='pydy':
            from pydy.viz.shapes import Cylinder
            from pydy.viz.visualization_frame import VisualizationFrame
            blade_shape = Cylinder(radius=radius, length=length, color=color)
            viz=[]
            if axis=='x':
                for iB in np.arange(nB):
                    frame  = self.frame.orientnew('b', 'Axis', (-np.pi/2+(iB-1)*2*np.pi/nB , self.frame.x) ) # Y pointing along blade
                    center=Point('RB'); 
                    center.set_pos(self.origin, length/2 * frame.y)
                    viz.append( VisualizationFrame(frame, center, blade_shape) )
                return viz
            else:
                raise NotImplementedError()
            
class YAMSRecSPBody(object):
    def __init__(B, name='', sympy=True):
        B.sympy = True
        B.pos_global  = B.vec3([0,0,0])
        B.R_b2g       = B.eye(3)
        B.name        = name
        B.Children    = []
        B.Connections = []
        B.MM     = None
        B.B           = [] # Velocity transformation matrix
        B.B_inB       = None
        B.BB_inB      = None
#         B.Bhat_x_bc   = None
#         B.Bhat_t_bc   = None
        B.I_DOF       = None
        B.gzf         = None        

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


    def __repr__(B):
        s='<YAMSRecSP Body {} object>:\n'.format(B.name)
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
        s+='| - pos_global:  {}\n'.format(B.pos_global)
        s+='| - gzf       :  {}\n'.format(B.gzf)
        s+='| - R_b2g : \n{}\n'.format(B.R_b2g)
        s+='| * nf  : {}\n'.format(B.nf)
        s+='| * R_bc: \n{}\n'.format(B.R_bc)
        s+='| * Bhat_x_bc: \n{}\n'.format(B.Bhat_x_bc)
        s+='| * Bhat_t_bc: \n{}\n'.format(B.Bhat_t_bc)
        s+='| * mass     :   {}\n'.format(B.mass) 
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
        if isinstance(o, GroundBody):
            raise Exception('Not intended to be called for Ground body')
        MqB      = fBMB(o.BB_inB, o.MM, sympy=o.sympy, name=o.name)
        n        = MqB.shape[0]
        M[:n,:n] = M[:n,:n]+MqB     
        for c in o.Children:
            M=c._getFullM(M)
        return M
        
    def _getFullK(o, K):
        if isinstance(o, GroundBody):
            raise Exception('Not intended to be called for Ground body')
            KqB      = fBMB(o.BB_inB, o.KK, sympy=o.sympy, name=o.name)
            n        = KqB.shape[0]
            K[:n,:n] = K[:n,:n]+KqB     
        for c in o.Children:
            K=c._getFullK(K)
        return K
        
    def _getFullD(o, D):
        if isinstance(o, GroundBody):
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
class YAMSInertialBody(YAMSBody):
    """ Inertial body / ground/ earth 
    Typically only one used
    """
    def __init__(self, name='E'): # "Earth"
        YAMSBody.__init__(self,name)
    

class GroundBody(YAMSRecSPBody):
    """ 
    Ground body is used to traverse the tree and hold the full mass matrix
    """
    def __init__(B):
        super(GroundBody,B).__init__('Grd')
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
        s+='|Inherits:\n'
        s+='||'+'\n|'.join(YAMSRecSPBody.__repr__(self).split('\n'))+'\n'
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
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class YAMSRigidBody(YAMSBody,SympyRigidBody):
    def __init__(self, name, mass=None, J_G=None, rho_G=None, J_form='full', J_at_Origin=False,
             name_for_var = None):
        """
        Define a rigid body and introduce symbols for convenience.
        
           Origin point have no velocities in the body frame! 
            
        
        INPUTS:
            name: string (can be one character), make sure this string is unique between all bodies
        
        OPTIONAL INPUTS:
            mass : scalar, body mass
            J_G  : 3x3 array or 3-array defining the coordinates of the inertia tensor in the body frame at the COG
            rho_G: array-like of length 3 defining the coordinates of the COG in the body frame
            J_form:  form of the inertial tensor
                   'full': all components are included 
                   'diag': only diagonal component  
                   'cross': diagonal and anti-diagonal components
            name_for_var: name to use for variables names. If None, name is used.
        
        
        """
        # YAMS Body creates a default "origin", "masscenter", and "frame"
        # We give the properties to SympyRigidBody later below
        super().__init__(name =name)
        
        if name_for_var is None:
            name_for_var = name
            self.name_for_var = name_for_var

        # --- Mass
        if mass is None:
            mass=Symbol('M_'+name_for_var)
        
        # --- Inertia, creating a dyadic using our frame and G
        if J_G is not None:
            if len(list(J_G))==3:
                ixx=J_G[0]
                iyy=J_G[1]
                izz=J_G[2]
                ixy, iyz, izx =0,0,0
            else:
                J_G = ensureMat(J_G, 3, 3)
                ixx = J_G[0,0]
                iyy = J_G[1,1]
                izz = J_G[2,2]
                izx = J_G[2,0]
                ixy = J_G[0,1]
                iyz = J_G[1,2]
        else:
            ixx = Symbol('J_xx_'+name_for_var)
            iyy = Symbol('J_yy_'+name_for_var)
            izz = Symbol('J_zz_'+name_for_var)
            izx = Symbol('J_zx_'+name_for_var)
            ixy = Symbol('J_xy_'+name_for_var)
            iyz = Symbol('J_yz_'+name_for_var)
        if J_form=='full':
            pass
        elif J_form=='diag':
            ixy, iyz, izx =0,0,0
        elif J_form=='cross':
            ixy, iyz =0,0
        else:
            raise NotImplementedError('J_form ',J_form)
            
        #inertia: dyadic : (inertia(frame, *list), point)
        if J_at_Origin:
            self._inertia_point = self.origin
            _inertia = (inertia(self.frame, ixx, iyy, izz, ixy, iyz, izx), self.origin)
            print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NOTE: inertia at origin')
        else:
            self._inertia_point = self.origin
            _inertia = (inertia(self.frame, ixx, iyy, izz, ixy, iyz, izx), self.masscenter)
            
        # --- Position of COG in body frame
        if rho_G is None: 
            rho_G=symbols('x_G_'+name_for_var+ ', y_G_'+name_for_var+ ', z_G_'+name_for_var)
        self.setGcoord(rho_G)
        self.masscenter.set_vel(self.frame, 0 * self.frame.x)
        # Init Sympy Rigid Body 
        #SympyRigidBody.__init__(self, name, self.masscenter, self.frame, mass, _inertia)
        self.mass = mass # Need to be  before inertia
        #self.masscenter # already defined
        #self.frame # already defined
        self.inertia = _inertia

        # For harmony with flexible bodies
        self.shapeNormSubs= []

    def bodyMassMatrix(self, form='regular', point='origin'):
        """ Body mass matrix in body coordinates M'(q)
        form is ['symbolic' , 'regular', 'TaylorExpanded']
        """
        self.M = zeros(6, 6)
        # Mxx
        self.M[0,0] = self.mass
        self.M[1,1] = self.mass
        self.M[2,2] = self.mass
        print('>>> bodyMassMatrix for rigid bodies is in Beta')

        if form=='TaylorExpanded':
            """ Return the term of the mass matrix at a given order.
            For order 1, terms are different for each dof"""
            # TODO
            self.M[0,0] = 0
            self.M[1,1] = 0
            self.M[2,2] = 0
            # Mxr, Mrx
            self.M[0:3,3:6] = - skew(self.mdCM.get(dof, order)) # =-skew(mdCM)  NOTE:2026 Changed sign
            self.M[3:6,0:3] = self.M[0:3,3:6].transpose()       # = skew(mdCM)
            # Mrr
            self.M[3:6,3:6] = self.J.get(dof,order)

        elif form=='regular':
            from welib.yams.utils import buildRigidBodyMassMatrix 
            self.M[0,0] = self.mass
            self.M[1,1] = self.mass
            self.M[2,2] = self.mass
            self.s_G_inB
            s_OG = self.masscenter_pos_local
            s_OG_ = s_OG.to_matrix(self.frame)
            self.M =  buildRigidBodyMassMatrix(self.mass, self.origin_inertia_matrix, s_OG_, symb=True) 

        elif form=='symbolic':
            # Mxr
            for i in np.arange(0,3):
                for j in np.arange(3,6):
                    self.M[i,j]=Symbol('M_{}{}{}'.format(self.name_for_var,i+1,j+1))
            self.M[0,3]=0
            self.M[1,4]=0
            self.M[2,5]=0
            # Mrr
            char='xyz'
            for i in np.arange(3,6):
                for j in np.arange(3,6):
                    self.M[i,j]=Symbol('J_{}{}{}'.format(self.name_for_var,char[i-3],char[j-3]))
            # Symmetry
            nq=0 # Rigid body
            for i in np.arange(0,6+nq):
                for j in np.arange(0,6+nq):
                    self.M[j,i]=self.M[i,j]
            pass
        else:
            raise NotImplementedError()
        return self.M

            
    def inertiaIsInPrincipalAxes(self):
        """ enforce the fact that the frame is along the principal axes"""
        D=self._inertia.to_matrix(self.frame)
        self.inertia=(inertia(self.frame, D[0,0], D[1,1], D[2,2]), self._inertia_point)
            
    def setGcoord(self, rho_G):
        """
        INPUTS:
            rho_G: array-like of length 3 defining the coordinates of the COG in the body frame
        """
        rho_G = ensureList(rho_G,3)
        self.s_G_inB = rho_G[0]*self.frame.x + rho_G[1]*self.frame.y+ rho_G[2]*self.frame.z # coordinates of 
        
        self.masscenter.set_pos(self.origin, self.s_G_inB)
        
    @property    
    def origin_inertia(self):
        return self.parallel_axis(self.origin)

    @property    
    def origin_inertia_matrix(self):
        return self.origin_inertia.to_matrix(self.frame)
    
    @property    
    def inertia_matrix(self):
        """ Returns inertia matrix in body frame at mass center"""
        J_G= self.parallel_axis(self.masscenter)
        return J_G.to_matrix(self.frame)

    @property
    def masscenter_pos_global(self):
        """ return masscenter position from inertial frame """
        return self.masscenter.pos_from(self.inertial_origin)

    @property
    def masscenter_pos_local(self):
        """ return masscenter position in body coordinate """
        return self.masscenter.pos_from(self.origin)

    @property
    def masscenter_vel_inertial(self):
        """ return velocity of masscenter in inertial frame """
        return self.masscenter.vel(self.inertial_frame)

    @property
    def masscenter_acc_inertial(self):
        """ return acceleration velocity of body COG in inertial frame """
        return self.masscenter.acc(self.inertial_frame)

    @property
    def kinetic_energy_inertial(self):
        return self.kinetic_energy(self.inertial_frame)
    
    def __repr__(self):
        # rigid body
        s=YAMSBody.__repr__(self)
        s+=' - mass:         {}\n'.format(self.mass)
        s+=' * inertia_matrix: {}\n'.format(self.inertia_matrix)
        s+='   (defined at point {})\n'.format(self.inertia[1])
        s+=' - masscenter:   {}\n'.format(self.masscenter)
        try:
            s+=' * masscenter_pos_global:   {}\n'.format(self.masscenter_pos_global)
        except:
            s+=' * masscenter_pos_global:   {}\n'.format('unknown')
        try:
            s+=' * masscenter_pos_local:    {}\n'.format(self.masscenter_pos_local)
        except:
            s+=' * masscenter_pos_local:    {}\n'.format('unknown')
        try:
            s+=' * masscenter_vel_inertial: {}\n'.format(self.masscenter_vel_inertial)
        except:
            s+=' * masscenter_vel_inertial: {}\n'.format('unknown')
        try:
            s+=' * masscenter_acc_inertial: {}\n'.format(self.masscenter_acc_inertial)
        except:
            s+=' * masscenter_acc_inertial: {}\n'.format('unknown')

        s+='Useful getters: origin_inertia, inertia_matrix, origin_inertia_matrix\n'
        s+='                masscenter_inertia\n'
        s+='                kinetic_energy_inertial\n'
        s+='Useful setters: noMass, noInertia, setGcoord\n'
        s+='Useful functions:\n'
        s+='  - kinetic_energy(frame)\n'
        s+='  - linear_momentum(point, frame)\n'
        s+='  - angular_momentum(point, frame)\n'
        s+='  - bodyMassMatrix(q=None, form="TaylorExpanded", order=None, dof=None)\n'
        return s



    def noInertia(self):
        """ set inertia to zero"""
        self.inertia=(inertia(self.frame, 0,0,0), self._inertia_point)

    def noMass(self):
        """ set inertia to zero"""
        self.mass=0

    #def zeroOrigin(self):
    #    """ set origin to zero"""
    #    self.origin.set_pos=(0,0,0)
    def kinetic_energy(self, frame):
        """ Taken from sympy.physics.mechanics.rigidbody.RigidBody.kinetic_energy"""
        from sympy.physics.vector  import dot
        from sympy import S
        rotational_KE = S.Half * dot(
            self.frame.ang_vel_in(frame),
            dot(self.central_inertia, self.frame.ang_vel_in(frame)))
        translational_KE = S.Half * self.mass * dot(self.masscenter.vel(frame), self.masscenter.vel(frame))
        return rotational_KE + translational_KE


        
class RigidBody(YAMSRecSPBody):
    def __init__(B, Name, Mass, J_G, rho_G):
        """
        Creates a rigid body 
        """
        super(RigidBody,B).__init__(Name)
        B.s_G_inB = rho_G
        B.J_G_inB = J_G  
        B.Mass    = Mass 

# --------------------------------------------------------------------------------}
# --- Flexible body/Beam Body 
# --------------------------------------------------------------------------------{
class YAMSFlexibleBody(YAMSBody):
    def __init__(self, name, nq, directions=None, orderMM=2, orderH=2, predefined_kind=None, name_for_var=None, name_for_DOF=None, tip_unit_deflect=False, tip_rotate=True, 
                 noZeroExp=False, singleDOFNumbering=True):
        """ 
        name:  name used for object name, origin
        name_for_var: name/string used for inertial variable names
        name_for_dof: name/string used for degree of freedom and connections

        """
        YAMSBody.__init__(self, name)
        if name_for_var is None:
            name_for_var=name
        if name_for_DOF is None:
            name_for_DOF=name # Important, not name_for_var

        self.name=name
        self.name_for_var = name_for_var
        self.name_for_DOF = name_for_DOF
        self.singleDOFNumbering = singleDOFNumbering
        self.L     = symbols('L_'+name_for_var)
        self.q     = []                         # DOF
        self.qd    = []                         # DOF velocity as "anonymous" variables
        self.qdot  = []                         # DOF velocities
        self.qddot = []                         # DOF accelerations
        t=dynamicsymbols._t
        if nq==1 and not singleDOFNumbering:
            for i in np.arange(nq):
                self.q.append   (dynamicsymbols('q_{}'. format(name_for_DOF)))
                self.qd.append  (dynamicsymbols('qd_{}'.format(name_for_DOF)))
        else:
            for i in np.arange(nq):
                self.q.append   (dynamicsymbols('q_{}{}'. format(name_for_DOF,i+1)))
                self.qd.append  (dynamicsymbols('qd_{}{}'.format(name_for_DOF,i+1)))
        for i in np.arange(nq):
            self.qdot.append(diff(self.q[i],t))
            self.qddot.append(diff(self.qdot[i],t))
        # --- Mass matrix related
        self.mass=symbols('M_{}'.format(name_for_var))
        self.J   = Taylor(name_for_var,'J'  , 3 , 3 , nq=nq, rname='xyz', cname='xyz', order=orderMM, noZeroExp=noZeroExp)
        self.Ct  = Taylor(name_for_var,'C_t', nq, 3 , nq=nq, rname=None , cname='xyz', order=orderMM, noZeroExp=noZeroExp        , singleDOFNumbering=singleDOFNumbering ) # TODO no expansion
        self.Cr  = Taylor(name_for_var,'C_r', nq, 3 , nq=nq, rname=None , cname=['x','y','z'], order=orderMM, noZeroExp=noZeroExp, singleDOFNumbering=singleDOFNumbering)
        self.Me  = Taylor(name_for_var,'M_e', nq, nq, nq=nq, rname=None , cname=None, order=orderMM, noZeroExp=noZeroExp, singleDOFNumbering=singleDOFNumbering)
        self.mdCM= Taylor(name_for_var,'M_d', 3,  1 , nq=nq, rname='xyz', cname=[''], order=orderMM, noZeroExp=noZeroExp)
        # --- h-omega related terms
        self.Gr=[0]*nq
        self.Ge=[0]*nq
        for i in np.arange(nq):
            self.Gr[i] = Taylor(name_for_var, 'G_r_{}'.format(i+1), 3,  3,  nq=nq, rname='xyz', cname='xyz', order=orderH, noZeroExp=noZeroExp)
            self.Ge[i] = SmpMat(name_for_var, 'G_e_{}'.format(i+1), nq, 3,  nq=nq, rname=None, cname='xyz', noZeroExp=noZeroExp, singleDOFNumbering=singleDOFNumbering)
#           self.Ge[i] = Taylor(name_for_var, 'G_e_{}'.format(i+1), nq, 3,  nq=nq, rname=None, cname='xyz', order=orderH, noZeroExp=noZeroExp    , singleDOFNumbering=singleDOFNumbering)
        self.Oe = Taylor(name_for_var, 'O_e', nq, 6,  nq=nq, rname=None, cname=['xx','yy','zz','xy','yz','xz'], order=orderH, noZeroExp=noZeroExp, singleDOFNumbering=singleDOFNumbering)
        # --- Stiffness and damping
        self.Ke  = Taylor(name_for_var,'K_e', nq, nq, nq=nq, rname=None , cname=None, order=1, noZeroExp=noZeroExp, singleDOFNumbering=singleDOFNumbering)
        self.De  = Taylor(name_for_var,'D_e', nq, nq, nq=nq, rname=None , cname=None, order=1, noZeroExp=noZeroExp, singleDOFNumbering=singleDOFNumbering)
        
        if len(directions)<nq:
            raise Exception(f'Number of directions should be at least {nq}. Directions provided are {directions}')
        self.directions=directions[:nq]
        self.defineExtremity(directions, unit_deflect=tip_unit_deflect, rotate=tip_rotate)
        
        self.origin = Point('O_'+self.name)
        # Properties from Rigid body
        # inertia

        self.shapeNormSubs= [(v,1) for v in self.ucList]
        
        # NOTE: masscenter put at origin for flexible bodies for now
        self.masscenter.set_pos(self.origin, 0*self.frame.x)

        self.predefined_kind=predefined_kind
        self.applyKindSimplification()
 
    def __repr__(self):
        # YAMS Flexible body
        s=YAMSBody.__repr__(self)
        s+=' - mass:         {}\n'.format(self.mass)
        #s+=' - inertia:      {}\n'.format(self.inertia[0].to_matrix(self.frame))
        #s+='   (defined at point {})\n'.format(self.inertia[1])
        #s+=' - masscenter:   {}\n'.format(self.masscenter)
        #s+='   (position from origin: {})\n'.format(self.masscenter.pos_from(self.origin))
        s+=' - q:            {}\n'.format(self.q)
        s+=' - qd:           {}\n'.format(self.qd)
        s+=' - ucList:       {}\n'.format(self.ucList)
        s+=' - vcList:       {}\n'.format(self.vcList)
        s+=' - uc    :       {}\n'.format(self.uc)
        s+=' - alpha :       {}\n'.format(self.alpha)
        s+=' - directions:   {}\n'.format(self.directions)
        s+='Useful functions:\n'
        s+='  - bodyMassMatrix(form="symbolic", point="origin")\n'
        return s


    def defineExtremity(self, directions=None, unit_deflect=False, rotate=True): 
        if directions is None:
            directions=['xyz']*len(self.q)
        # Hard coding 1 connection at beam extremity
        self.alpha =[dynamicsymbols('alpha_x{}'.format(self.name_for_DOF)),dynamicsymbols('alpha_y{}'.format(self.name_for_DOF)),dynamicsymbols('alpha_z{}'.format(self.name_for_DOF))]
#         if unit_deflect:
#             self.uc    =[0,0,0]
#         else:
        self.uc    =[dynamicsymbols('u_x{}c'.format(self.name_for_DOF)),dynamicsymbols('u_y{}c'.format(self.name_for_DOF)),0]
        alphax = 0
        alphay = 0
        alphaz = 0
        uxc = 0
        uyc = 0
        vList =[]
        uList =[]
        for i in np.arange(len(self.q)):
            u=0
            v=0
            if 'x' in directions[i]:
                if unit_deflect:
                    u = 1
                else:
                    u = symbols('u_x{}{}c'.format(self.name_for_DOF,i+1))
                if rotate:
                    v = symbols('v_y{}{}c'.format(self.name_for_DOF,i+1)) 
                else:
                    v = 0
                uxc   += u * self.q[i]
                alphay+= v * self.q[i]
            if 'y' in directions[i]:
                if unit_deflect:
                    u = 1
                else:
                    u = symbols('u_y{}{}c'.format(self.name_for_DOF,i+1))
                if rotate:
                    v = symbols('v_x{}{}c'.format(self.name_for_DOF,i+1))
                else:
                    v = 0
                uyc   += u * self.q[i]
                alphax+= v * self.q[i]
            if 'z' in directions[i]:
                v = symbols('v_z{}{}c'.format(self.name_for_DOF,i+1))
                alphaz+= v * self.q[i]
            vList.append(v)
            uList.append(u)
        
        self.alphaSubs=[(self.alpha[0], alphax ), (self.alpha[1],alphay), (self.alpha[2],alphaz)]
        self.ucSubs   =[(self.uc[0], uxc ), (self.uc[1],uyc)]
        self.v2Subs = [(v1*v2,0) for v1 in vList for v2 in vList] # Second order terms put to 0

        self.ucList =uList # "PhiU" values at connection point for each mode
        self.vcList =vList # "PhiV" values at connection point for each mode

    def bodyMassMatrix(self, q=None, form='TaylorExpanded', order=None, dof=None):
        """ Body mass matrix in body coordinates M'(q)
        form is ['symbolic' , 'TaylorExpanded']
        """
        nq=len(self.q)
        if q is None:
            q = self.q
        else:
            if len(q)!=len(self.q):
                raise Exception('Inconsistent dimension between q ({}) and body nq ({}) for body {}'.format(len(q),nq,self.name))
        self.M=zeros(6+nq,6+nq)
        # Mxx
        self.M[0,0] = self.mass
        self.M[1,1] = self.mass
        self.M[2,2] = self.mass

        if order is not None and dof is not None:
            """ Return the term of the mass matrix at a given order.
            For order 1, terms are different for each dof"""
            self.M[0,0] = 0
            self.M[1,1] = 0
            self.M[2,2] = 0
            # Mrx, Mxr
            self.M[0:3,3:6] = - skew(self.mdCM.get(dof, order)) # = -skew(mdCM)
            self.M[3:6,0:3] = self.M[0:3,3:6].transpose()       # =  skew(mdCM)
            # Mrr
            self.M[3:6,3:6] = self.J.get(dof,order)
            # Mgx, Mxg
            self.M[6:6+nq,0:3] = self.Ct.get(dof,order)
            self.M[0:3,6:6+nq] = self.M[6:6+nq,0:3].transpose()
            # Mrg
            self.M[6:6+nq,3:6] = self.Cr.get(dof,order)
            self.M[3:6,6:6+nq] = self.M[6:6+nq,3:6].transpose()
            # Mgg
            self.M[6:6+nq,6:6+nq] = self.Me.get(dof,order)

        elif form=='symbolic':
            # Mxr, Mxg
            for i in np.arange(0,3):
                for j in np.arange(3,6+nq):
                    self.M[i,j]=Symbol('M_{}{}{}'.format(self.name_for_var,i+1,j+1))
            self.M[0,3]=0
            self.M[1,4]=0
            self.M[2,5]=0
            # Mrr
            char='xyz'
            for i in np.arange(3,6):
                for j in np.arange(3,6):
                    self.M[i,j]=Symbol('J_{}{}{}'.format(self.name_for_var,char[i-3],char[j-3]))
            # Mrg
            for i in np.arange(3,6):
                for j in np.arange(6,6+nq):
                    self.M[i,j]=Symbol('M_{}{}{}'.format(self.name_for_var,i+1,j+1))
            # Mgg
            for i in np.arange(6,6+nq):
                for j in np.arange(6,6+nq):
                    self.M[i,j]=Symbol('GM_{}{}{}'.format(self.name_for_var,i-5,j-5))

            for i in np.arange(0,6+nq):
                for j in np.arange(0,6+nq):
                    self.M[j,i]=self.M[i,j]
            pass

        elif form=='TaylorExpanded':
            # We evaluate
            # Mrx, Mxr
            self.M[0:3,3:6] = - skew(self.mdCM.eval(q))    # = - skew(mdCM)
            self.M[3:6,0:3] = self.M[0:3,3:6].transpose()  # =   skew(mdCM)
            # Mrr
            self.M[3:6,3:6] = self.J.eval(q)
            # Mgx, Mxg
            self.M[6:6+nq,0:3] = self.Ct.eval(q)
            self.M[0:3,6:6+nq] = self.M[6:6+nq,0:3].transpose()
            # Mrg
            self.M[6:6+nq,3:6] = self.Cr.eval(q)
            self.M[3:6,6:6+nq] = self.M[6:6+nq,3:6].transpose()
            # Mgg
            self.M[6:6+nq,6:6+nq] = self.Me.eval(q)

        else:
            raise Exception('Unknown mass matrix form option `{}`. Allowed are: `TaylorExpanded` or `symbolic`'.format(form))

        if self.predefined_kind is not None:
            if self.predefined_kind=='twr-z':
                # --- We assume a tower symmetric along z
                # "x-theta"
                self.M[0,3]   = 0
                self.M[0,5]   = 0
                self.M[1,4]   = 0
                self.M[1,5]   = 0
                self.M[2,3]   = 0
                self.M[2,4]   = 0
                self.M[2,5]   = 0
                self.M[1,3] = -self.M[0,4]
                # "theta-theta"
                # "x-theta"
                self.M[3,4]   = 0
                self.M[3,5]   = 0
                self.M[4,3]   = 0
                self.M[4,5]   = 0
                # shape
                for iq in np.arange(nq):
                    xyz = self.directions[iq]
                    if xyz=='y':
                        self.M[0,6+iq]=0
                        self.M[2,6+iq]=0
                        self.M[4,6+iq]=0
                        self.M[5,6+iq]=0
                    if xyz=='x':
                        self.M[1,6+iq]=0
                        self.M[2,6+iq]=0
                        self.M[3,6+iq]=0
                        self.M[5,6+iq]=0
                    for jq in np.arange(nq):
                        xyz2 = self.directions[jq]
                        # Assume no coupling x-y
                        if xyz2!=xyz:
                            self.M[6+iq,6+jq]=0

                # symmetry
                for i in np.arange(self.M.shape[0]):
                    for j in np.arange(self.M.shape[0]):
                        if i<j:
                            self.M[j,i]=self.M[i,j]
            elif self.predefined_kind=='bld-z':
                print('>>> yams_simpy, TODO simplifications of mass matrix for blade')
                # symmetry inertial tensor
                for i in np.arange(3,6):
                    for j in np.arange(3,6):
                        if i<j:
                            self.M[j,i]=self.M[i,j]
            elif self.predefined_kind=='bld-z-straight':
                #print('>>> yams_simpy, predifined kind',self.predefined_kind)
                # inertial tensor diagonal
                for i in np.arange(3,6):
                    for j in np.arange(3,6):
                        if i!=j:
                            self.M[j,i]=0
                            self.M[i,j]=0
                # COG purely on z
                self.M[5,0]=0 # MdBy =0
                self.M[0,5]=0
                self.M[3,2]=0 # MdBy =0
                self.M[2,3]=0
                self.M[5,1]=0 # MdBx =0
                self.M[1,5]=0
                self.M[4,2]=0 # MdBx =0
                self.M[2,4]=0
                # shape
                #print('>>>> DIRECTIONS', self.directions)
                for iq in np.arange(nq):
                    xyz = self.directions[iq]
                    if xyz=='y':
                        self.M[0,6+iq]=0
                        self.M[2,6+iq]=0
                        self.M[4,6+iq]=0
                        self.M[5,6+iq]=0
                    if xyz=='x':
                        self.M[1,6+iq]=0
                        self.M[2,6+iq]=0
                        self.M[3,6+iq]=0
                        self.M[5,6+iq]=0
                    for jq in np.arange(nq):
                        xyz2 = self.directions[jq]
                        # Assume no coupling x-y
                        if xyz2!=xyz:
                            self.M[6+iq,6+jq]=0
                # symmetry
                for i in np.arange(self.M.shape[0]):
                    for j in np.arange(self.M.shape[0]):
                        if i<j:
                            self.M[j,i]=self.M[i,j]
                

            else:
                raise NotImplementedError()
        
        return self.M


    def bodyMassMatrixNonLin(self, form='TaylorExpanded'):
        nq = len(self.q)
        MNL = zeros(6+nq,6+nq)

        MMloc = self.bodyMassMatrix(form=form) # TODO need symbolic form with Wallrapp notation

        # --- NEW CORRECTION
        M_theta_theta_1, M_theta_theta_2 = self.M_theta_theta_expansion()

        for j, q_j in enumerate(self.q):
            # --- First order terms
            # --- M x_theta_1
            C_tj = MMloc[0:3, 6 + j] # Vector
            MNL[0:3,3:6] += - skew(C_tj) * q_j
            # --- M_theta_theta_1
            MNL[3:6,3:6] += M_theta_theta_1[j] * q_j
            # --- M_theta_e_1
            MNL[3:6,6:] += sp.S.Half * self.Ge[j].T * q_j 
            # --- Second order terms
            # --- M_theta_theta_2
            for k, q_k in enumerate(self.q):
                if j==k:
                    # TODO TODO
                    MNL[3:6,3:6] += M_theta_theta_2[j][k] * q_j * q_k
        # --- Make it symmetric
        MNL[3:6, 0:3] = MNL[0:3,3:6].T
        MNL[6: , 3:6] = MNL[3:6, 6:].T
        return MNL

    def M_theta_theta_expansion(self, Mform='TaylorExpanded'):
        """ Return nonlinear terms of M theta theta 
        Note: computed analytically on 20/7/2026 as part of 673 extract notes, to be placed elsewhere.
        """
        MMloc = self.bodyMassMatrix(form=Mform)

        M_theta_theta_1 = [zeros(3,3)] * len(self.q)
        M_theta_theta_2 = [[zeros(3,3)] * len(self.q)] * len(self.q)
        # Loop on body DOFs
        for j, q_j in enumerate(self.q):
            # --- M_theta_theta_1
            C_rj = MMloc[3:6, 6 + j] # Cr row j, a vector of length 3
            M_theta_theta_1[j]= - skew(C_rj)

            # --- M_theta_theta_2
            for k, q_k in enumerate(self.q):
                M_ejk = MMloc[6 + j, 6 + k] # Modal mass
                if 'z' in np.array(self.directions).flatten():
                    raise NotImplementedError('M_theta_theta correction only implemented for z-beam')
                M_theta_theta_2[j][k][2, 2] = M_ejk # Good as long as beam along z
                if len(np.unique(np.array(self.directions).flatten()))!=len(self.q):
                    raise NotImplementedError('Shape function directions needs to be unique for M_theta_tehta correction for now..')
                # NOTE: simplifications for now, assume all modes are in different directions
                if k==j:
                    if self.directions[j]=='x':
                        M_theta_theta_2[j][k][1, 1] = M_ejk
                    elif self.directions[j]=='y':
                        M_theta_theta_2[j][k][0, 0] = M_ejk
                    else:
                        raise NotImplementedError('Only pure x and y directions supported for M_theta_theta correction')
        return M_theta_theta_1, M_theta_theta_2


    
    def bodyQuadraticForce(self, omega, q, qd, form='TaylorExpanded', nonLinCorr=False):
        r""" Body quadratic force  k_\omega (or h_omega)  (centrifugal and gyroscopic)
        inputs:
           omega: angular velocity of the body wrt to the inertial frame, expressed in body coordinates
           q,qd: generalied coordinates and speeds for this body
        """
        # --- Safety
        q     = ensureMat(q , len(q), 1)
        qd    = ensureMat(qd, len(qd), 1)
        omega = ensureMat(omega, 3, 1)
        nq=len(self.q)
        if len(q)!=nq:
            raise Exception('Inconsistent dimension between q ({}) and body nq ({}) for body {}'.format(len(q),nq,self.name))
            
        # --- Init
        k_omega = Matrix(np.zeros((6+nq,1)).astype(int)) 
        om_til = skew(omega)
        ox,oy,oz=omega[0,0],omega[1,0],omega[2,0]
        omega_q = Matrix([ox**2, oy**2,oz**2, ox*oy, oy*oz, ox*oz]).reshape(6,1)
        
        if not nonLinCorr:
            if form=='TaylorExpanded':
                # --- k_omega_t
                k_omega[0:3,0] =  2 *om_til * transpose(self.Ct.eval(q)) * qd # NOTE we use star instead of dot because Ct and qd are sympy Matrix
                # k_omega[0:3,0] += om_til * skew(self.mdCM.eval(q)) * omega # NOTE: Wrong sign convention mdCM sign convention is opposite Wallrap
                k_omega[0:3,0] += om_til * (om_til * self.mdCM.eval(q))  # True expression 
                #k_omega[0:3,0] += - om_til * skew(self.mdCM.eval(q)) * omega # Alternative from true expression by reverting the cross product

                # --- k_omega_r 
                k_omega[3:6,0] = om_til * self.J.eval(q) * omega
                for k in np.arange(nq):
                    k_omega[3:6,0] += self.Gr[k].eval(q) * qd[k] * omega

                # --- k_omega_e
                k_omega[6:6+nq,0] = self.Oe.eval(q) * omega_q
                for k in np.arange(nq):
                    #k_omega[6:6+nq,0] += self.Ge[k].eval(q) * qd[k] * omega
                    k_omega[6:6+nq,0] += self.Ge[k] * qd[k] * omega
            else:
                raise NotImplementedError()
        else:
            M_theta_theta_1, M_theta_theta_2 = self.M_theta_theta_expansion()

            # --- NOTE: only the nonlinCorr here
            # --- k_omega_t
            #k_omega[0:3,0] =  2 *om_til * transpose(self.Ct.eval(q)) * qd # No expansion of Ct
            for j, q_j in enumerate(q):
                mdCM_1_j = self.Ct.M0[j,:].T * q_j
                k_omega[0:3,0] += om_til * (om_til *  mdCM_1_j )  # True expression 

            # --- k_omega_r - Term 1
            M_theta_theta_ = zeros(3, 3)
            for j, q_j in enumerate(q):
                M_theta_theta_ += M_theta_theta_1[j] * q_j
                for k, q_k in enumerate(q):
                    M_theta_theta_ += M_theta_theta_2[j][k] * q_j * q_k
            k_omega[3:6,0] += om_til * ( M_theta_theta_ ) * omega

            # --- k_omega_r - Term 2
            for j in np.arange(nq):
                Gr_ = zeros(3, 3)
                for k, q_k in enumerate(q):
                    Gr_ += 2* M_theta_theta_2[k][j] * q_k
                k_omega[3:6,0] += Gr_ * qd[j] * omega

            # --- k_omega_e
            for j in np.arange(nq):
                Oe_1 = zeros(3,3)
                for k, q_k in enumerate(q):
                    Oe_1 += - M_theta_theta_2[k][j] * q_k
                k_omega[6+j,0] +=  (omega.T * Oe_1 * omega)[0,0]

        return k_omega


    
    def bodyElasticForce(self, q, qd):
        # --- Safety
        q     = ensureMat(q , len(q), 1)
        qd    = ensureMat(qd, len(qd), 1)
        nq=len(self.q)
        if len(q)!=nq:
            raise Exception('Inconsistent dimension between q ({}) and body nq ({}) for body {}'.format(len(q),nq,self.name))
        # --- Init
        ke = Matrix(np.zeros((6+nq,1)).astype(int)) 
        ke[6:6+nq,0] = self.Ke.eval(q)*q + self.De.eval(q)*qd
        return ke

    def bodyGravitationalForce(self, g_vect, q, form='TaylorExpanded'):
        r""" Body gravity force  h_g, acts as an external force on the right hand side
            M a + k_{\omega} + k_e = f_{ext,other} + h_g

        inputs:
           g_vect: gravity vector, expressed in body coordinates
           q: generalized coordinates for this body
        """
        # --- Safety
        q     = ensureMat(q , len(q), 1)
        g_vect = ensureMat(g_vect, 3, 1)
        nq=len(self.q)
        if len(q)!=nq:
            raise Exception('Inconsistent dimension between q ({}) and body nq ({}) for body {}'.format(len(q),nq,self.name))
            
        # --- Init
        M33 = Matrix(np.zeros((3,3)))
        M33[0,0] = self.mass
        M33[1,1] = self.mass
        M33[2,2] = self.mass

        h_g = Matrix(np.zeros((6+nq,1)).astype(int)) 
        if form=='TaylorExpanded':
            # h_g,t
            h_g[0:3,0] = M33.dot(g_vect)                      # f_g    = \int m g_vect
            # h_g,t
            h_g[3:6,0] =  skew(self.mdCM.eval(q)).dot(g_vect) # \tau_g = \int r_{CG} \times (m g_vect) = \int (m r_{CG}) \times g_vect
            # h_g_e
            h_g[6:6+nq,0] =  (self.Ct.eval(q)).dot(g_vect)    #f_{g,e} = \int   m(z) \Phi(z)^T g_vect dz = C_t \, g_vect
        else:
            raise NotImplementedError()

        return h_g


    
    def connectToTip(parent, child, type, rel_pos, rot_type='Body', rot_amounts=None, rot_order=None, rot_order_elastic='XYZ', rot_type_elastic='SmallRot', doSubs=True):
        """
        The connection between a flexible body and another body is similar to the connections between rigid bodies
        The relative position and rotations between bodies are modified in this function to include the elastic 
        displacements and rotations.
        
        For now, we only allow 1 connection for each flexible body.
        
        rel_pos: for flexible bodies, this is the position when the body is undeflected!
        rot_amounts: for flexible bodies, these are the rotations when the body is undeflected!
                     These rotations are assumed to occur AFTER the rotations from the deflection
                     
        s_PC  : vector from parent origin to child origin
        s_PC0 : vector from parent origin to child origin when undeflected
        
        s_PC = s_PC0 + u

        NOTE:
          - by default rot_type_elatic is "SmallRot" 
                e.g. 
                   R_b2g = [ 1     0  nu*q ]
                           [ 0     1    0  ]
                           [-nu*q  0    1  ]
        
          - by default doSubs is True
              alpha_y - > nu y*q 
        """
        rel_pos = [r + u for r,u in zip(rel_pos, parent.uc)]
        # Computing DCM due to elastic motion
        M_B2e = rotToDCM(rot_type_elastic, rot_amounts = parent.alpha, rot_order=rot_order_elastic) # from parent to deformed parent 
        # Insert elastic DOF
        if doSubs:
            rel_pos = [r.subs(parent.ucSubs) for r in rel_pos]
            M_B2e   = M_B2e.subs(parent.alphaSubs)
        # Full DCM
        if rot_amounts is None:
            M_c2B = M_B2e.transpose()
        else:
            M_e2c = rotToDCM(rot_type,         rot_amounts = rot_amounts,  rot_order=rot_order) # from parent to deformed parent 
            M_c2B = (M_e2c.transpose() * M_B2e.transpose() )
        
        #print('Rel pos with Elastic motion:', rel_pos)
        #print('Rel rot with Elastic motion:', rot_amounts)
            
        YAMSBody.connectTo(parent, child, type, rel_pos=rel_pos, rot_type='DCM', rot_amounts=M_c2B, dynamicAllowed=True)
        #YAMSBody.connectTo(parent, child, type, rel_pos=rel_pos, dynamicAllowed=True)
        

    def applyKindSimplification(self):
        """ set inertia to zero"""
        if self.predefined_kind is not None:
            if self.predefined_kind=='twr-z' or self.predefined_kind=='bld-z-straight':
                nq=len(self.q)
                for iq in np.arange(nq):
                    xyz = self.directions[iq]
                    if xyz=='y' or xyz=='x':
                        # Ge
                        for jq in np.arange(nq):
                            self.Ge[iq][jq,0]=0
                            self.Ge[iq][jq,1]=0
                            xyz2 = self.directions[jq]
                            if xyz==xyz2:
                                self.Ge[iq][jq,2]=0
                            if xyz!=xyz2:
                                self.Me.M0[iq,jq]=0
                                self.De.M0[iq,jq]=0
                                self.Ke.M0[iq,jq]=0
                        # Gr
                        for i,xyz2 in enumerate(['x','y','z']):
                            self.Gr[iq].M0[i,0]=0
                            self.Gr[iq].M0[i,1]=0
                            if xyz!=xyz2:
                                self.Gr[iq].M0[i,2]=0
                    else:
                        raise NotImplementedError()
                    # Oe
                    for j in np.arange(4):
                        self.Oe.M0[iq,j]=0
                    for j,xyz2 in [(4,'x'),(5,'y')]:
                        if xyz2==xyz:
                            self.Oe.M0[iq,j]=0
                    # Ct, Cr
                    if xyz=='y':
                        self.Ct.M0[iq,0]=0
                        self.Cr.M0[iq,1]=0
                    if xyz=='x':
                        self.Ct.M0[iq,1]=0
                        self.Cr.M0[iq,0]=0
                    self.Ct.M0[iq,2]=0
                    self.Cr.M0[iq,2]=0
                # J
                self.J.M0[0,1]= 0
                self.J.M0[0,2]= 0
                self.J.M0[1,0]= 0
                self.J.M0[1,2]= 0
                self.J.M0[2,0]= 0
                self.J.M0[2,1]= 0
                # Mdcm
                self.mdCM.M0[0]=0
                self.mdCM.M0[1]=0


    def replaceDict(self, rd=None, form='TaylorExpanded'):
        if rd is None:
            rd = OrderedDict()
        else:
            rd.update(rd)

        # --- Mass matrix
        MM = self.bodyMassMatrix(q = [0]*len(self.q), form = form)

        rd[repr(MM[0,0])] = ('MM_{}'.format(self.name_for_var), [0,0])
        rd[repr(self.mass)] = ('MM_{}'.format(self.name_for_var), [0,0])

        for i in np.arange(0,MM.shape[0]):
            for j in np.arange(3,MM.shape[0]):
                s=repr(MM[i,j])
                if len(s)>1:
                    if s[0]!='-':
                        rd[s] = ('MM_{}'.format(self.name_for_var), [i,j])
        if form=='TaylorExpanded':
            if hasattr(self.Me, 'M1'):
                # We retrieve first order terms
                for iq in np.arange(len(self.q)):
                    MM1 =  self.bodyMassMatrix(dof=iq, order=1)
                    for i in np.arange(0,MM1.shape[0]):
                        for j in np.arange(3,MM1.shape[0]):
                            s=repr(MM1[i,j])
                            if len(s)>1:
                                if s[0]!='-':
                                    rd[s] = ('MM1_{}'.format(self.name_for_var), [iq,i,j])

        for iq in np.arange(len(self.Gr)):
            for i in np.arange(self.Gr[iq].M0.shape[0]):
                for j in np.arange(self.Gr[iq].M0.shape[1]):
                    s=repr(self.Gr[iq].M0[i,j])
                    if len(s)>1:
                        rd[s] = ('Gr_{}'.format(self.name_for_var), [iq,i,j])

        for iq in np.arange(len(self.Ge)):
            for i in np.arange(self.Ge[iq].shape[0]):
                for j in np.arange(self.Ge[iq].shape[1]):
                    s=repr(self.Ge[iq][i,j])
                    if len(s)>1:
                        rd[s] = ('Ge_{}'.format(self.name_for_var), [iq,i,j])

        for i in np.arange(self.Oe.M0.shape[0]):
            for j in np.arange(self.Oe.M0.shape[1]):
                s=repr(self.Oe.M0[i,j])
                if len(s)>1:
                    rd[s] = ('Oe_{}'.format(self.name_for_var), [i,j])

        for i in np.arange(self.Ke.M0.shape[0]):
            for j in np.arange(self.Ke.M0.shape[1]):
                s=repr(self.Ke.M0[i,j])
                if len(s)>1:
                    rd[s] = ('KK_{}'.format(self.name_for_var), [i+6,j+6])

        for i in np.arange(self.De.M0.shape[0]):
            for j in np.arange(self.De.M0.shape[1]):
                s=repr(self.De.M0[i,j])
                if len(s)>1:
                    rd[s] = ('DD_{}'.format(self.name_for_var), [i+6,j+6])
        return rd
    

    @property
    def curvilinear_coord(self):
        s=Symbol('s') # TODO, could use name for var
        return s

    def Phi(self, form='function', var=None, full=False):
        """ Return matrix of shape function displacement field"""
        if var is None:
            var = self.curvilinear_coord
        if full:
            directions=['xyz']*len(self.q)
        else:
            directions = self.directions

        Phis = zeros( 3, len(self.q) )
        if form =='function':
            for iq, axes in enumerate(directions):
                if 'x' in axes:
                    Phis[iq, 0] = Function('Phi_'+self.name_for_var+str(iq+1)+'_x')(var)
                if 'y' in axes:
                    Phis[iq, 1] = Function('Phi_'+self.name_for_var+str(iq+1)+'_y')(var)
        else:
            for iq, axes in enumerate(self.directions):
                if 'x' in axes:
                    Phis[iq, 0] = Symbol('Phi_'+self.name_for_var+str(iq+1)+'_x') 
                if 'y' in axes:
                    Phis[iq, 1] = Symbol('Phi_'+self.name_for_var+str(iq+1)+'_y') 
        return Phis

    def uP(self, form='function', var=None):
        """ Return displacement field at point P"""
        uP = zeros(3,1)
        Phis = self.Phi(form=form, var=var)
        for j, qj in enumerate(self.q):
            uP += qj*Phis[:, j]
        return uP

    def udotP(self, form='function', var=None):
        # udot = sum qdot_j * Phi_j
        Phis = self.Phi(form = form, var=var)
        udotP = sp.zeros(3, 1)
        for j, qdotj in enumerate(self.qdot):
            udotP += qdotj * Phis[:, j]
        return udotP


    def r0(self):
        """ 
        r_0 = (x0, y0, z0)_B is the undeformed position vector, with coords in body frame
        """
        x0,y0,z0 = sp.symbols('x_0, y_0 z_0')
        if self.predefined_kind == 'twr-z':
            x0=0
            y0=0
#         else:
#             raise NotImplementedError('predefined kind {self.predefined_kind}')
        r0 = sp.Matrix([x0, y0, z0])
        return r0

    def rP(self, form='function', var=None):
        """ 
         r_P = r_0 + u = r_0 + sum q_j Phi_j
        """
        r0 = self.r0()
        if self.predefined_kind == 'twr-z':
            var = r0[2,0] # z0
        rP = r0 + self.uP(form=form, var=var)
        return rP
            

    def KE_origin_vel(self):
        # NOTE: those better be unique symbols across the framework and bodies
        vOx, vOy, vOz = sp.symbols('v_Ox, v_Oy, v_Oz') # Body Origin Velocity in body frame
        omx, omy, omz = sp.symbols('omega_x, omega_y, omega_z') # Body angular velocity in body frame
        vO_symb = sp.Matrix([vOx, vOy, vOz])
        om_symb = sp.Matrix([omx, omy, omz])
        if self.inertial_frame is None:
            subs_dict = None
        else:
            vO_coord = self.vel_inertial.to_matrix(self.frame).simplify()
            om_coord = self.omega_inertial.to_matrix(self.frame).simplify()
            subs_dict = {
                vOx: vO_coord[0], vOy: vO_coord[1], vOz: vO_coord[2],
                omx: om_coord[0], omy: om_coord[1], omz: om_coord[2]
            }
        #print('>>> Origin velocity')
        #print(vO_symb)
        #print(om_symb)
        return vO_symb, om_symb, subs_dict

    def kinetic_energy(self, frame=None, subs=False, method='analytical_atoms'):

        # ---  Define arbitrary velocity of body using symbols
        vO_symb, om_symb, origin_vel_subs = self.KE_origin_vel()

        # --- Kinematics of point P
        rP = self.rP()
        r0 = self.r0()
        if self.predefined_kind == 'twr-z':
            var = r0[2, 0] 
        else:
            var= None
        Phis  = self.Phi(var=var)
        udotP = self.udotP(var=var)
        #print(Phis)
        #print(udotP)
        
        # --- Compute individual atoms
        if method =='analytical_atoms':
            T1 = self._KE_atom1_translation(vO_symb)
            T2 = self._KE_atom2_trans_rot(vO_symb, om_symb, rP)
            T3 = self._KE_atom3_trans_elastic(vO_symb, udotP)
            T4 = self._KE_atom4_rotation(om_symb, rP)
            T5 = self._KE_atom5_rot_elastic(om_symb, rP, udotP)
            T6 = self._KE_atom6_pure_elastic(udotP)
            T_symb = T1 + T2 + T3 + T4 + T5 +  T6
        elif method == 'direct':
            T_symb = self._KE_direct_identification(vO_symb, om_symb)

        # --- Final step: replace placeholders with true body-frame kinematics
        if self.inertial_frame is None:
            print('>>>>> KE: inertial_frame is None, keeping things symbolic')
            return T_symb
        else:
            print('subs_dict', origin_vel_subs)
            if subs:
                return T_symb.subs(origin_vel_subs).simplify()
            else:
                return T_symb


    def _KE_atom1_translation(self, vO):
        """Atom 1: 0.5 * m * (vO . vO)"""
        return sp.Rational(1, 2) * self.mass * vO.dot(vO)

    def _KE_atom2_trans_rot(self, vO, om, rP):
        """Atom 2: vO . (om x integral(rP dm))"""
        # integral(rP dm) = S + sum(q_j * Ct_j)
        S_total = self.mdCM.M0
        
        # Add elastic contribution to center of mass moment if Ct exists
        for j, qj in enumerate(self.q):
            S_total += qj * self.Ct.M0[j, :].T
                
        return vO.dot(om.cross(S_total))

    def _KE_atom3_trans_elastic(self, vO, udotP):
        """Atom 3: vO . integral(udotP dm) = vO . sum(qdot_j * Ct_j)"""
        udot_integrated = sp.zeros(3, 1)
        for j, qdotj in enumerate(self.qdot):
            udot_integrated += qdotj * self.Ct.M0[j, :].T
        return vO.dot(udot_integrated)

    def _KE_atom4_rotation(self, om, rP):
        """Atom 4: 0.5 * om^T * J(q) * om"""
        # Retrieve or assemble J(q) = J0 + sum(q_j * Je_j) + ...
        J_q = self.J.M0
        # TODO
        #for j, qj in enumerate(self.q):
        #    J_q += qj * self.Je[j]
        # 1st-order extension: J_q += sum_j (q_j * Je_j)
        #for j, qj in enumerate(self.q):
        #    J_q += qj * self.Je[j]
                
        # 2nd-order extension: J_q += sum_{j,k} (q_j * q_k * Je_jk)
        #for j, qj in enumerate(self.q):
        #    for k, qk in enumerate(self.q):
        #        J_q += qj * qk * self.Je_jk[j][k]
                
        return sp.Rational(1, 2) * (om.T * J_q * om)[0]

    def _KE_atom5_rot_elastic(self, om, rP, udotP):
        """Atom 5: om . integral(rP x udotP dm)"""
        # Evaluates to om . (sum(qdot_k * Cr_k) + sum(q_j * qdot_k * Cr_jk))
        rot_elastic_term = sp.zeros(3, 1)
        
        for k, qdotk in enumerate(self.qdot):
            rot_elastic_term += qdotk * self.Cr.M0[k, :].T
                
        # TODO
        #for j, qj in enumerate(self.q):
        #    for k, qdotk in enumerate(self.qdot):
        #        rot_elastic_term += qj * qdotk * self.Cr_jk[j][k]
                    
        return om.dot(rot_elastic_term)

    def _KE_atom6_pure_elastic(self, udotP):
        """Atom 6: 0.5 * qdot^T * Me * qdot"""
        qdot_vec = sp.Matrix(self.qdot)
        return sp.Rational(1, 2) * (qdot_vec.T * self.Me.M0 * qdot_vec)[0]


    def bodyMassMatrixFromKE(self, T=None, method='analytical_atoms'):
        """
        Computes the generic body mass matrix M(q) from the kinetic energy T(q, v_O, om, qdot)
        by differentiating with respect to the body velocity vector:
            nu = [v_Ox, v_Oy, v_Oz, omega_x, omega_y, omega_z, qdot_1, ..., qdot_N]^T
        
        Returns:
            sp.Matrix of shape (6 + n_modes, 6 + n_modes)
        """
        # 1. Compute kinetic energy with symbolic placeholders if not supplied
        if T is None:
            T = self.kinetic_energy(subs=False, method=method)
            
        # 2. Re-create the exact symbolic placeholders used in kinetic_energy
        vOx, vOy, vOz = sp.symbols('v_Ox, v_Oy, v_Oz')
        omx, omy, omz = sp.symbols('omega_x, omega_y, omega_z')
        
        # 3. Assemble the full velocity vector nu
        nu = sp.Matrix([vOx, vOy, vOz, omx, omy, omz] + list(self.qdot))
        n_dof = len(nu)
        
        # 4. Compute M_ij = d^2(T) / (d nu_i d nu_j)
        M = sp.zeros(n_dof, n_dof)
        for i in range(n_dof):
            # First gradient element dT / d(nu_i)
            dT_dnu_i = sp.diff(T, nu[i])
            
            # Second derivative for lower triangle + diagonal
            for j in range(i, n_dof):
                val = sp.diff(dT_dnu_i, nu[j])
                M[i, j] = val
                if i != j:
                    M[j, i] = val # Enforce symmetry
                    
        return M


    def _KE_direct_identification(self, vO_symb, om_symb):
            """
            Computes 1/2 * v_P . v_P symbolically, expands spatial terms, 
            and matches spatial integrands against shape integral definitions.
            """
            # 1. Point P kinematics
            r0 = self.r0()
            if self.predefined_kind == 'twr-z':
                var = r0[2, 0] 
            else:
                var = None
            rP = self.rP(var=var) # 3x1 Matrix depending on spatial var (e.g. z_0)
            udotP = self.udotP(var=var) # 3x1 Matrix depending on qdot and Phis
            
            # 2. Local point velocity vector v_P
            vP = vO_symb + om_symb.cross(rP) + udotP
            
            # 3. Scalar kinetic energy density = 1/2 * vP . vP
            # Expanding scalar dot product ensures additive terms
            T_density = sp.Rational(1, 2) * vP.dot(vP)
            T_expanded = sp.expand(T_density)
            
            # Identify spatial variable (e.g. z_0)
            z0 = self.r0()[2, 0] if self.predefined_kind == 'twr-z' else sp.Symbol('z_0')
            
            # Identify all time-dependent variables to hold constant
            time_symbols = set(vO_symb) | set(om_symb) | set(self.q) | set(self.qdot)
            
            # Split into additive scalar terms
            terms = sp.Add.make_args(T_expanded)
            
            T_integrated = 0
            for term in terms:
                # Separate time-dependent factors from spatial factors
                spatial_part, time_part = term.as_independent(*time_symbols, as_coeff_prod=True)
                # Pull out pure numerical coefficients (e.g., 1/2, 2, -1) from spatial_part
                num_coeff, pure_spatial = spatial_part.as_coeff_Mul()
                
                # Identify shape integral symbol for the pure spatial integrand
                shape_integral_symb = self._identify_shape_integral(pure_spatial, z0)

                T_integrated += num_coeff * time_part * shape_integral_symb
                
            return T_integrated

    def _identify_shape_integral(self, spatial_expr, s_var):
        """
        Pattern matches spatial expressions inside integral( rho * spatial_expr dz )
        to generate symbolic shape integrals.
        """
        # 1. Constant term (int(1 dz) -> M_T)
        # Pure constant spatial factor
        if spatial_expr == 1:
             return self.mass

        # 2. First mass moment (int(z dz) -> S_z)
        if spatial_expr == s_var:
            return self.mdCM.M0[2]
            #return sp.Symbol('M_dTz')

        # 3. Second mass moment / Inertia terms (int(z^2 dz) -> J_z2)
        if spatial_expr == s_var**2:
            return sp.Symbol('J_z2')

        # 4. Shape function terms
        Phis = self.Phi(var=s_var)
        for j, qj in enumerate(self.q):
            phi_j = Phis[:, j]

            # Linear translation coupling: int(Phi_j_dim dz) -> C_t_j_dim
            for dim_idx, dim_name in enumerate(['x', 'y', 'z']):
                if spatial_expr == phi_j[dim_idx]:
                    return self.Ct.M0[j, dim_idx]
                #return sp.Symbol(f'C_t_{j+1}_{dim_name}')
                if spatial_expr == s_var * phi_j[dim_idx]:
                    if dim_name == 'y':
                        return -self.Cr.M0[j, 0] # Crx
                    elif dim_name == 'x':
                        return self.Cr.M0[j, 1] # Cry
                    else:
#                         if dim_name == 'x':
#                     print('>>>>>>>>>>>> TODO TODO C_t_z')
                        return sp.Symbol(f'C_t_z_{j+1}_{dim_name}')

            # Modal mass coupling: int(Phi_j . Phi_k dz) -> M_e_j_k
            for k, qk in enumerate(self.q):
                phi_k = Phis[:, k]

                # Check vector dot product match
                if spatial_expr == phi_j.dot(phi_k):
                    return self.Me.M0[j,k]
                #return sp.Symbol(f'M_e_{j+1}_{k+1}')

                # Check component-wise products (e.g., Phi_j_x * Phi_k_x)
                for dim_idx in range(3):
                    if spatial_expr == phi_j[dim_idx] * phi_k[dim_idx]:
                        #return sp.Symbol(f'M_e_{j+1}_{k+1}')
                        return self.Me.M0[j,k]
                #return sp.Symbol(f'M_e_{j+1}_{k+1}')

        print(f"Unrecognized spatial integrand structure: {spatial_expr}")
        rho = sp.symbols(r'\rho')
        L   = sp.symbols(r'L')
        return sp.Integral(rho * spatial_expr, (s_var, 0, L))




# --------------------------------------------------------------------------------}
# --- Beam Body 
# --------------------------------------------------------------------------------{
class BeamBody(YAMSRecSPBody):
    def __init__(B,Name,nf,main_axis='z',nD=2):
        super(BeamBody,B).__init__(Name)
        B.PhiU = [None] * nf # B.nf has no setter 
        B.nD  = nD
        B.main_axis = main_axis
    @property
    def alpha_couplings(self):
        return self.Bhat_t_bc @ self.gzf

    @property
    def R_bc(self):
        if self.main_axis=='x':
            alpha_y= symbols('alpha_y') #-p.V(3,iNode);
            alpha_z= symbols('alpha_z') # p.V(2,iNode);
            return R_y(alpha_y)*R_z(alpha_z)

        elif self.main_axis=='z':
            alpha_x= symbols('alpha_x') #-p.V(2,iNode);
            alpha_y= symbols('alpha_y') # p.V(1,iNode);
            return R_x(alpha_x)*R_y(alpha_y)
        else:
            raise NotImplementedError()

    @property
    def Bhat_x_bc(self):
        #      Bx_pc(:,j)=p.PhiU{j}(:,iNode);
        Bhat_x_bc = Matrix(np.zeros((3,self.nf)).astype(int))
        if self.main_axis=='z':
            for j in np.arange(self.nf):
                if j<self.nf/2 or self.nD==1:
                    Bhat_x_bc[0,j]=symbols('ux{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along x
                else:
                    Bhat_x_bc[1,j]=symbols('uy{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along y
        elif self.main_axis=='x':
            for j in np.arange(self.nf):
                if j<self.nf/2 or self.nD==1:
                    Bhat_x_bc[2,j]=symbols('uz{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along z
                else:
                    Bhat_x_bc[1,j]=symbols('uy{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along y
        return Bhat_x_bc

    @property
    def Bhat_t_bc(self):
        #      Bt_pc(:,j)=[0; -p.PhiV{j}(3,iNode); p.PhiV{j}(2,iNode)];
        Bhat_t_bc = Matrix(np.zeros((3,self.nf)).astype(int))
        if self.main_axis=='z':
            for j in np.arange(self.nf):
                if j<self.nf/2  or self.nD==1:
                    Bhat_t_bc[1,j]=symbols('vy{:d}c'.format(j+1))
                else:
                    Bhat_t_bc[0,j]=-symbols('vx{:d}c'.format(j+1))
        elif self.main_axis=='x':
            for j in np.arange(self.nf):
                if j<self.nf/2 or self.nD==1:
                    Bhat_t_bc[1,j]=-symbols('vz{:d}c'.format(j+1))
                else:
                    Bhat_t_bc[2,j]=symbols('vy{:d}c'.format(j+1))
        return Bhat_t_bc




# --------------------------------------------------------------------------------}
# --- Rotation 
# --------------------------------------------------------------------------------{
def R_x(t):
    return Matrix( [[1,0,0], [0,cos(t),-sin(t)], [0,sin(t),cos(t)]])
def R_y(t):
    return Matrix( [[cos(t),0,sin(t)], [0,1,0], [-sin(t),0,cos(t)] ])
def R_z(t): 
    return Matrix( [[cos(t),-sin(t),0], [sin(t),cos(t),0], [0,0,1]])
# --------------------------------------------------------------------------------}
# --- B Matrices 
# --------------------------------------------------------------------------------{
def fB_inB(R_EI, B_I, sympy=True):
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

def fB_aug(B_I_inI, nf_I, nf_Curr=None, nf_Prev=None, sympy=True):
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


def fBMatRecursion(Bp, Bhat_x, Bhat_t, R0p, r_pi, sympy=True):
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

    #r_pi=colvec(r_pi)

    # TODO use Translate here
    Bi = MatrixLoc(np.zeros((6,ni+n_p)))
    for j in range(n_p):
        Bi[:3,j] = Bp[:3,j]+cross(Bp[3:,j],r_pi) # Recursive formula for Bt mentioned after Eq.(15)
        Bi[3:,j] = Bp[3:,j] # Recursive formula for Bx mentioned after Eq.(12)
    if ni>0:
        Bi[:3,n_p:] = R0p @ Bhat_x[:,:] # Recursive formula for Bx mentioned after Eq.(15)
        Bi[3:,n_p:] = R0p @ Bhat_t[:,:] # Recursive formula for Bt mentioned after Eq.(12)
    return Bi

def fBMatTranslate(Bp, r_pi, sympy=True):
    """
    Rigid translation of a B matrix to another point, i.e. transfer the velocities from a point to another: 
      - translational velocity:  v@J = v@I + om@I x r@IJ
      - rotational velocity   : om@J = om@I
    """
    Bi=np.zeros(Bp.shape)
    if Bp.ndim==1:
        raise NotImplementedError

    for j in range(Bp.shape[1]):
        Bi[0:3,j] = Bp[0:3,j]+np.cross(Bp[3:6,j],r_pi)
        Bi[3:6,j] = Bp[3:6,j]
    return Bi


def fBMB(BB_I_inI, MM, sympy=True, name=''):
    """ Computes the body generalized matrix: B'^t M' B 
    See Eq.(8) of [1] 
    """
    if MM is None:
        raise Exception(f'MM is None for body {name}')
    MM_I = (np.transpose(BB_I_inI) @ MM) @ BB_I_inI
    return MM_I



if __name__ == "__main__":
    x, y, z = dynamicsymbols('x, y, z')
    phi_x, phi_y, phi_z = dynamicsymbols('phi_x, phi_y, phi_z')
    ref = YAMSInertialBody('E')

    rot = YAMSRigidBody('R')
    twr = YAMSFlexibleBody('T', 1, directions=['x'], predefined_kind='twr-z')
#     twr = YAMSFlexibleBody('T', nDOF_twr, directions=opts['twrDOFDir'][:nDOF_twr], orderMM=opts['orderMM'], orderH=opts['orderH'], 
#                            predefined_kind='twr-z', tip_unit_deflect=opts['twr_tip_unit_deflect'], tip_rotate=opts['twr_tip_rotate'],
#                            noZeroExp=not opts['singleExpNumbering'], singleDOFNumbering=opts['singleDOFNumbering'])

#     ref.connectTo(twr, type='Free' , rel_pos=[x,y,z], rot_type='Body', rot_amounts=[phi_x,phi_y,phi_z], rot_order='XYZ')


    ref.connectTo(rot, type='Free', rel_pos=[x,y,z], rot_type='Body', rot_amounts=[phi_x,phi_y,phi_z], rot_order='XYZ')

    print(rot.kinetic_energy(ref.frame))
