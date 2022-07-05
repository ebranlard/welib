"""
Generic bodies classes
These classes will be used for more advanced classes:
    - new and old YAMS body classes for Sympy
    - YAMS body for numerical yams
"""
from welib.yams.utils import translateInertiaMatrixToCOG, translateInertiaMatrixFromCOG
from welib.yams.utils import rigidBodyMassMatrix 
from welib.yams.utils import R_x, R_y, R_z
from welib.yams.flexibility import GMBeam, GKBeam, GKBeamStiffnening, GeneralizedMCK_PolyBeam
from welib.yams.flexibility import checkRegularNode
# from welib.yams.utils import skew


__all__ = ['Body','InertialBody','RigidBody','FlexibleBody']

# --- For harmony with sympy
import numpy as np
from numpy import eye, cross, cos ,sin
def Matrix(m):
    return np.asarray(m)
def zeros(m,n):
    return np.zeros((m,n))

# --------------------------------------------------------------------------------}
# --- Generic Body 
# --------------------------------------------------------------------------------{
class Body(object):
    """
    Base class for rigid bodies and flexible bodies
    """
    def __init__(self, name='', r_O=[0,0,0], R_b2g=np.eye(3)):
        self.name = name
        self._r_O   = np.asarray(r_O).ravel()
        self._R_b2g = np.asarray(R_b2g)

        self._mass=None
        self.MM  = None # To be defined by children

    def __repr__(self):
        s='<Generic Body object>:\n'
        return s

    @property
    def Mass(self):
        raise Exception('`Mass` is an old interface, use `mass` instead')

    @property
    def mass(self):
        return self._mass

    @property    
    def pos_global(self):
        """ Position of origin in global coordinates """
        return self._r_O

    @pos_global.setter
    def pos_global(self, r_O):
        self._r_O = np.asarray(r_O).ravel()

    @property
    def R_b2g(self):
        """ Transformation matrix from body to global """
        return self._R_b2g

    @R_b2g.setter
    def R_b2g(self, R_b2g):
        self._R_b2g = R_b2g

    @property
    def R_g2b(self):
        """ Transformation matrix from global to body """
        return self._R_b2g.transpose() 

    def pos_local(self, x_gl):
        """ return position vector from origin of body, in body coordinates, of a point in global """
        return self.R_g2b.dot(x_gl - self._r_O)


# --------------------------------------------------------------------------------}
# --- Ground Body 
# --------------------------------------------------------------------------------{
class InertialBody(Body):
    def __init__(self, name='Grd'):
        Body.__init__(self, name=name)


# --------------------------------------------------------------------------------}
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class RigidBody(Body):
    def __init__(self, name, mass, J, s_OG, r_O=[0,0,0], R_b2g=np.eye(3), s_OP=None):
        """
        Creates a rigid body 

        INPUTS:
         - name: name of body (string)
         - mass: body mass (float)
         - J: inertia tensor (array-like) in body frame, defined either at:
               - center of mass G, located at s_OG from the body origin
               - OR point P, located at s_OP from the body origin
             J may be defined as:
               - a 3x3 matrix
               - a 3-vector (Jxx, Jyy, Jzz) representing the diagonal values
               - a 6-vector (Jxx, Jyy, Jzz, Jxy, Jyz, Jzx) representing the diagonal values
         - s_OG: vector from body origin to body COG in body frame 

         - s_OP: vector from body origin to point where inertia is defined,
                 in body frame
                 (only if inertia is not defined at COG).
         - r_O: vector from global origin to body origin, in global coordinates
         - R_b2g : transformation matrix from body to gobal coordinates

        """
        Body.__init__(self, name, r_O=r_O, R_b2g=R_b2g)
        self._mass  = mass
        self._s_OG = np.asarray(s_OG).ravel()

        # Ensuring a 3x3 inertia matrix
        J = np.asarray(J)
        Jflat=J.ravel()
        if len(Jflat)==3:
            J = np.diag(Jflat)
        elif len(Jflat)==6:
            J = np.diag(Jflat[:3])
            J[0,1]=J[1,0]=Jflat[3]
            J[1,2]=J[2,1]=Jflat[4]
            J[1,3]=J[3,1]=Jflat[5]
            
        # inertia at COG
        if s_OP is not None:
            s_PG=  self._s_OG - s_OP
            self._J_G = translateInertiaMatrixToCOG(J, mass, s_PG)
        else:
            self._J_G = J

    def shiftOrigin(self, s_OOnew):
        """ change body origin
        s_OOnew: vector from old origin to new origin
        """
        s_OnewG    = -np.asarray(s_OOnew) + self._s_OG
        self._s_OG = s_OnewG

    # --------------------------------------------------------------------------------
    # --- Inertia
    # --------------------------------------------------------------------------------
    @property    
    def masscenter(self):
        """ Position of mass center in body frame"""
        return self._s_OG

    @property
    def masscenter_pos_global(self):
        """ return masscenter position from inertial frame """
        try:
            return self._r_O + self.R_b2g.dot(self._s_OG)
        except:
            import pdb; pdb.set_trace()

    @property    
    def inertia(self):
        """ Inertia at body origin in body coordinates"""
        return self.inertia_at([0,0,0])

    @property    
    def masscenter_inertia(self):
        """ Returns inertia matrix at COG in body frame"""
        return self._J_G

    def inertia_at(self, s_OP, R_f2g=None):
        """ returns body inertia at a given point, and given frame (default body frame)
        INPUTS:
         - s_OP: point coordinates from body origin in body coordinates
         - R_f2g: transformation matrix from a given frame when inertia is wanted to global
        """
        # 
        s_GP =   np.asarray(s_OP) - self._s_OG
        J = translateInertiaMatrixFromCOG(self._J_G, self.mass, s_GP)
        if R_f2g is not None:
            R_b2f = np.dot(R_f2g.T, self.R_b2g)
            J = R_b2f.dot(J).dot(R_b2f.T)
        return J

    @property
    def mass_matrix(self):
        """ Body mass matrix at origin"""
        return rigidBodyMassMatrix(self.mass, self.inertia, self._s_OG) # TODO change interface

    def mass_matrix_at(self, s_OP):
        """ Body mass matrix at a given point"""
        J = self.inertia_at(s_OP)
        s_PG = -np.asarray(s_OP)+ self._s_OG
        return rigidBodyMassMatrix(self.mass, J, s_PG) # TODO change interface

    def __repr__(self):
        s='<RigidBody object>:\n'.format(self.name)
        s+=' * pos_global:            {} (origin)\n'.format(np.around(self.pos_global,6))
        s+=' * masscenter:            {} (body frame)\n'.format(np.around(self.masscenter,6))
        s+=' * masscenter_pos_global: {} \n'.format(np.around(self.masscenter_pos_global,6))
        s+=' - mass:         {}\n'.format(self.mass)
        s+=' * R_b2g: \n {}\n'.format(self.R_b2g)
        s+=' * masscenter_inertia: \n{}\n'.format(np.around(self.masscenter_inertia,6))
        s+=' * inertia: (at origin)\n{}\n'.format(np.around(self.inertia,6))
        s+='Useful getters: inertia_at, mass_matrix\n'
        return s

    def combine(self, other, name=None, R_b2g=np.eye(3), r_O=None):
        """ Combine two rigid bodies and form a new rigid body
        
        """
        M = self.mass + other.mass
        x_G = (self.mass * self.masscenter_pos_global + other.mass * other.masscenter_pos_global)/M

        if name is None:
            name=self.name + other.name

        # Inertias in new body frame and at new COG
        s_O1_G = self.pos_local(x_G)
        s_O2_G = other.pos_local(x_G)
        J1 = self.inertia_at(s_O1_G, R_b2g)
        J2 = other.inertia_at(s_O2_G, R_b2g)
        #print('s_O1_G ',s_O1_G)
        #print('s_O2_G ',s_O2_G)
        #print('J1\n ',J1)
        #print('J2\n ',J2)
        #print('J12\n ',J1+J2)

        if r_O is None:
            # Putting origin of new body at COG of common body
            r_O  = x_G
            s_OG = [0,0,0]
        else:
            s_OG = (R_b2g.T).dot(x_G-r_O)
        return RigidBody(name, M, J1+J2, s_OG, r_O=r_O, R_b2g=R_b2g)


# --------------------------------------------------------------------------------}
# --- Flexible Body 
# --------------------------------------------------------------------------------{
class FlexibleBody(Body):
    def __init__(self, name, 
            r_O=[0,0,0], R_b2g=np.eye(3) # Position and orientation in global
            ):
        """
        Creates a Flexible body 
        """
        Body.__init__(self, name, r_O=r_O, R_b2g=R_b2g)

# --------------------------------------------------------------------------------}
# --- Beam Body 
# --------------------------------------------------------------------------------{
class BeamBody(FlexibleBody):
    def __init__(self, name, s_span, s_P0, m, EI, PhiU, PhiV, PhiK, jxxG=None, s_G0=None, 
            s_min=None, s_max=None,
            r_O=[0,0,0], R_b2g=np.eye(3), # Position and orientation in global
            damp_zeta=None, RayleighCoeff=None, DampMat=None,
            bAxialCorr=False, bOrth=False, Mtop=0, Omega=0, bStiffening=True, gravity=None, main_axis='z', massExpected=None,
            int_method='Flex'
            ):
        """
        Creates a Flexible Beam body 
          Points P0 - Undeformed mean line of the body
        """
        FlexibleBody.__init__(self, name, r_O=r_O, R_b2g=R_b2g)
        self.main_axis = main_axis
        self.s_span = s_span
        if s_min is None:
            self.s_min = np.min(s_span)
        else:
            self.s_min = s_min
        if s_max is None:
            self.s_max = np.max(s_span)
        else:
            self.s_max = s_max
        self.m      = m
        self.s_G0   = s_G0
        self.PhiU   = PhiU
        self.PhiV   = PhiV
        self.PhiK   = PhiK
        self.jxxG   = jxxG
        self.s_P0   = s_P0
        self.EI     = EI
        if jxxG is None:
            self.jxxG   = 0*m
        if self.s_G0 is None:
            self.s_G0=self.s_P0
    
        self.s_G    = self.s_G0
        self.bAxialCorr = bAxialCorr
        self.bOrth      = bOrth
        self.bStiffening= bStiffening
        self.Mtop       = Mtop
        self.Omega      = Omega # rad/s
        self.gravity    = gravity
        self.int_method = int_method

        self.damp_zeta  = damp_zeta
        self.RayleighCoeff  = RayleighCoeff
        self.DampMat        = DampMat

        if massExpected is not None:
            self.computeMassMatrix()
            Mass = self.MM[0,0]
            factor = Mass/massExpected
            if np.abs(factor-1)>1e-5:
                print('>>>BeamBody: Scaling mass distribution with factor {:.4f} in order to get a desired mass of {}'.format(factor,massExpected))
                self.m /= factor

        self.computeMassMatrixTaylor()
        self.computeMassMatrix(s_G = self.s_G0, inPlace=True)
        self.computeStiffnessMatrix()
        self.computeDampingMatrix(damp_zeta)

        # TODO
        #self.V0         = np.zeros((3,self.nSpan))
        #self.K0         = np.zeros((3,self.nSpan))
        #self.rho_G0_inS = np.zeros((3,self.nSpan)) # location of COG in each cross section
        #[o.PhiV,o.PhiK] = fBeamSlopeCurvature(o.s_span,o.PhiU,o.PhiV,o.PhiK,1e-2);
        #[o.V0,o.K0]     = fBeamSlopeCurvature(o.s_span,o.s_P0,o.V0,o.K0,1e-2)    ;

    def toRigidBody(self):
        """ Create a rigid body from a flexible body """
        return RigidBody(self.name+'_rigid', self.mass, self.masscenter_inertia, self.masscenter, r_O=self.pos_global, R_b2g=self.R_b2g)

    def computeStiffnessMatrix(B, Mtop=None, Omega=None):
        B.KK0 = GKBeam(B.s_span, B.EI, B.PhiK, bOrth=B.bOrth, method=B.int_method)

        if Mtop is not None:
            B.Mtop=Mtop
        if Omega is not None:
            B.Omega=Omega

        if B.bStiffening:
            B.KKg_self = GKBeamStiffnening(B.s_span, B.PhiV, B.gravity, B.m, B.Mtop, B.Omega, main_axis=B.main_axis, bSelfWeight=True,  bMtop=False, bRot=False, method=B.int_method)
            B.KKg_Mtop = GKBeamStiffnening(B.s_span, B.PhiV, B.gravity, B.m, B.Mtop, B.Omega, main_axis=B.main_axis, bSelfWeight=False, bMtop=True , bRot=False, method=B.int_method)
            B.KKg_rot  = GKBeamStiffnening(B.s_span, B.PhiV, B.gravity, B.m, B.Mtop, B.Omega, main_axis=B.main_axis, bSelfWeight=False, bMtop=False, bRot=True , method=B.int_method)
            B.KKg = B.KKg_self+B.KKg_Mtop+B.KKg_rot
        else:
            B.KKg      = B.KK0*0
            B.KKg_self = B.KK0*0
            B.KKg_Mtop = B.KK0*0
            B.KKg_rot  = B.KK0*0

        B.KK=B.KK0+B.KKg
        if len(np.isnan(B.KK))>0:
            #print('>>> WARNING, some stiffness matrix values are nan, replacing with 0')
            B.KK[np.isnan(B.KK)]=0

    def computeDampingMatrix(self, damp_zeta=None):
        self.DD = np.zeros((6+self.nf,6+self.nf))
        if damp_zeta is None:
            return
        for j,zeta in enumerate(damp_zeta):
            gm = self.MM[6+j,6+j]
            gk = self.KK[6+j,6+j]
            if gk<0:
                print('[WARN] Bodies: Shape function has negative stiffness (likely due to geometrical stiffening)')
                om = 0
            else:
                om            = np.sqrt(gk/gm)
            xi = zeta*2*np.pi
            c  = xi * gm * om / np.pi
            self.DD[6+j,6+j] = c

        if len(np.isnan(self.DD))>0:
            #print('>>> WARNING, some damping matrix values are nan, replacing with 0')
            self.DD[np.isnan(self.DD)]=0


    @property    
    def start_pos(self):
        """ start of body wrt origin """
        return self.s_P0[:,0]
    @property    
    def end_pos(self):
        """ end of body wrt origin """
        return self.s_P0[:,-1]
    # --------------------------------------------------------------------------------}
    # --- Inertia 
    # --------------------------------------------------------------------------------{
    @property    
    def mass(self):
        """ Body mass"""
        return self.MM[0,0]

    @property    
    def masscenter(self):
        """ Position of mass center in body frame"""
        if self.mass>0:
            if self.int_method =='OpenFAST':
                S = self.first_moment_inertia
                #FirstMom = sum(p['BElmntMass']*p['RNodes'])
                return np.array(S/self.mass)
            else:
                return  np.trapz(self.m*self.s_G0,self.s_span)/self.mass
        else:
            return np.array([0,0,0])

    @property
    def masscenter_pos_global(self):
        """ return masscenter position from inertial frame """
        return self._r_O + self.R_b2g.dot(self.masscenter)

    @property    
    def inertia(self):
        """ Returns inertia matrix at Origin in body frame"""
        return self.MM[3:6,3:6]

    @property    
    def masscenter_inertia(self):
        """ Returns inertia matrix at COG in body frame. 
        NOTE: this is approximate for flexible bodies
        """
        return translateInertiaMatrixToCOG(self.inertia, self.mass, self.masscenter)

    def inertia_at(self, s_OP, R_f2g=None):
        """ returns body inertia at a given point, and given frame (default body frame)
        NOTE: this is approximate for flexible bodies
        INPUTS:
         - s_OP: point coordinates from body origin in body coordinates
         - R_f2g: transformation matrix from a given frame when inertia is wanted to global
        """
        # 
        s_GP =   np.asarray(s_OP) - self.masscenter
        J = translateInertiaMatrixFromCOG(self.masscenter_inertia, self.mass, s_GP)
        if R_f2g is not None:
            R_b2f = np.dot(R_f2g.T, self.R_b2g)
            J = R_b2f.dot(J).dot(R_b2f.T)
        return J

    def inertia_global_at(self, r_P_gl):
        """ returns body inertia in global frame at a point given in global coordinates
        INPUTS:
         - s_P_gl: point coordinates in global coordinates
        """
        r_O  = self.pos_global
        r_OP = np.array(r_P_gl)-r_O
        s_OP = self.R_g2b.dot(r_OP) 
        return self.inertia_at(s_OP, R_f2g=np.eye(3))


    @property    
    def first_moment_inertia(self):
        """ Returns first moment of inertia from body origin"""
        if self.int_method =='OpenFAST':
            dr = checkRegularNode(self.s_span)
            s_span = self.s_span[1:-1] # NOTE: temporary, m shouldn't me used with this method
            m      = self.m[1:-1] *dr   # Important HACK 
            s_G    = self.s_G[:,1:-1]
            #np.sum(yy) 
            #p['FirstMom']  = sum(p['BElmntMass']*p['RNodes'])    + p['TipMass']*p['BldFlexL']               # wrt blade root    
            S1x = np.sum(s_G[0,:]*m)
            S1y = np.sum(s_G[1,:]*m)
            S1z = np.sum(s_G[2,:]*m)
            return S1x, S1y, S1z

        else:
            raise NotImplementedError()
    @property
    def first_moment_inertia_from_start(self):
        """ Returns first moment of inertia from start position of body (not origin)"""
        if self.int_method =='OpenFAST':
            dr = checkRegularNode(self.s_span)
            s_span = self.s_span[1:-1] # NOTE: temporary, m shouldn't me used with this method
            m      = self.m[1:-1] *dr   # Important HACK 
            s_G    = self.s_G[:,1:-1]
            #np.sum(yy) 
            #p['FirstMom']  = sum(p['BElmntMass']*p['RNodes'])    + p['TipMass']*p['BldFlexL']               # wrt blade root    
            S0 = self.start_pos
            S1x = np.sum((s_G[0,:]-S0[0])*m)
            S1y = np.sum((s_G[1,:]-S0[1])*m)
            S1z = np.sum((s_G[2,:]-S0[2])*m)
            return S1x, S1y, S1z

        else:
            raise NotImplementedError()

    @property
    def mass_matrix(self):
        """ Body mass matrix at origin"""
        return self.MM

    def mass_matrix_at(self, s_OP):
        """ Body mass matrix at a ginve point"""
        J = self.inertia_at(s_OP)
        s_PG = -np.asarray(s_OP)+ self._s_OG
        return rigidBodyMassMatrix(self.mass, J, s_PG) # TODO change interface


    def updateFlexibleKinematics(B, qe, qep):
        """ see yams.py updateKinematics"""
        # Deflections shape
        B.U  = np.zeros((3,B.nSpan));
        B.V  = np.zeros((3,B.nSpan));
        B.K  = np.zeros((3,B.nSpan));
        B.UP = np.zeros((3,B.nSpan));
        for j in range(B.nf):
            B.UP[0:3,:] +=  qep[j] * B.PhiU[j][0:3,:]
            B.U [0:3,:] +=  qe[j]  * B.PhiU[j][0:3,:]
            B.V [0:3,:] +=  qe[j]  * B.PhiV[j][0:3,:]
            B.K [0:3,:] +=  qe[j]  * B.PhiK[j][0:3,:]
        #B.V_tot=B.V+B.V0;
        #B.K_tot=B.K+B.K0;
        # Position of mean line
        #B.s_P=B.s_P0+B.U;
        # Position of COG TODO see yams.py for a better treatment
        B.s_G = B.s_G0 + B.U;

    def computeMassMatrixTaylor(B):
        """ Compute Taylor expansion of the mass matrix """
        s_G_bkp = B.s_G

        delta = 0.1
        qep   = np.zeros(B.nf)
        B.MM1 =[]
        B.Gr1 =[]
        B.Ge1 =[]
        B.Oe1 =[]
        B.Oe61=[]
        for j in np.arange(B.nf):
            # Positive perturbation
            qe  = np.zeros(B.nf); qe[j]=+delta; B.updateFlexibleKinematics(qe,qep)
            MMp, ITp = B.computeMassMatrix(s_G=B.s_G, inPlace=False)
            # Negative perturbation
            qe  = np.zeros(B.nf); qe[j]=-delta; B.updateFlexibleKinematics(qe,qep)
            MMm, ITm = B.computeMassMatrix(s_G=B.s_G, inPlace=False)
            # Derivative
            ITd = dict()
            for k in ['Gr','Ge','Oe','Oe6']:
                ITd[k] = (ITp[k]-ITm[k])/(2*delta)
            MM = (MMp-MMm)/(2*delta)
            # Store in object
            B.MM1.append(MM)
            B.Gr1.append (ITd['Gr'])
            B.Ge1.append (ITd['Ge'])
            B.Oe1.append (ITd['Oe'])
            B.Oe61.append(ITd['Oe6'])
        # Revert back 
        B.s_G = s_G_bkp


    def computeMassMatrix(B, s_G = None, inPlace=True):
        if s_G is None:
            s_G = B.s_G
        MM, IT = GMBeam(s_G, B.s_span, B.m, B.PhiU, jxxG=B.jxxG, method=B.int_method, main_axis=B.main_axis, bAxialCorr=B.bAxialCorr, bOrth=B.bOrth, rot_terms=True)
        if len(np.isnan(MM))>0:
            #print('>>> WARNING, some mass matrix values are nan, replacing with 0')
            MM[np.isnan(MM)]=0
        if inPlace:
            B.MM, B.Gr, B.Ge, B.Oe, B.Oe6 = MM, IT['Gr'], IT['Ge'], IT['Oe'], IT['Oe6']
        return MM, IT


    @property
    def length(B):
        return B.s_max-B.s_min

    @property
    def nSpan(B):
        return len(B.s_span)

    @property
    def nf(B):
        return len(B.PhiU)

    @property
    def Bhat_x_bc(self,iNode=-1):
        Bhat_x_bc = Matrix(np.zeros((3,self.nf)))
        for j in np.arange(self.nf):
            Bhat_x_bc[:,j]=self.PhiU[j][:,iNode] #  along x
        return Bhat_x_bc

    @property
    def Bhat_t_bc(self,iNode=-1):
        """ unit "alpha" couplings """
        Bhat_t_bc = Matrix(np.zeros((3,self.nf)))
        for j in np.arange(self.nf):
            if self.main_axis=='x':
                Bhat_t_bc[0,j]=0                      # torsion
                Bhat_t_bc[1,j]=-self.PhiV[j][2,iNode]
                Bhat_t_bc[2,j]= self.PhiV[j][1,iNode]
            elif self.main_axis=='z':
                Bhat_t_bc[0,j]=-self.PhiV[j][1,iNode]
                Bhat_t_bc[1,j]= self.PhiV[j][0,iNode]
                Bhat_t_bc[2,j]= 0                     # torsion
        return Bhat_t_bc

    def __repr__(self):
        s='<BeamBody {} object>:\n'.format(self.name)
        s+=' * pos_global:            {} (origin)\n'.format(np.around(self.pos_global,6))
        s+=' * masscenter:            {} (body frame)\n'.format(np.around(self.masscenter,6))
        s+=' * masscenter_pos_global: {} \n'.format(np.around(self.masscenter_pos_global,6))
        s+=' - mass:         {}\n'.format(self.mass)
        s+=' * length:      {}\n'.format(self.length)
        s+=' * R_b2g: \n {}\n'.format(self.R_b2g)
        s+=' * masscenter_inertia: \n{}\n'.format(np.around(self.masscenter_inertia,6))
        s+=' * inertia: (at origin)\n{}\n'.format(np.around(self.inertia,6))
        s+=' - Properties: s_span, m, EI, Mtop, PhiU, PhiV, PhiW\n'
        s+='               MM, KK, KK0, KKg, KKg_Mtop, KKg_self\n'
        s+='Usefull getters: inertia_at, mass_matrix_at, toRigidBody \n'
        return s

# --------------------------------------------------------------------------------}
# --- FAST Beam body 
# --------------------------------------------------------------------------------{
class FASTBeamBody(BeamBody):
    def __init__(self, ED, inp, Mtop=0, shapes=None, main_axis='z', nSpan=None, bAxialCorr=False, bStiffening=True, jxxG=None, Omega=0,
            spanFrom0=False,
            bldStartAtRotorCenter=True,
            massExpected=None,
            gravity=None,
            algo=''):
        """ 
        INPUTS:
           ED: ElastoDyn inputs as read from weio
           inp: blade or tower file, as read by weio
           Mtop: top mass if any
           nSpan: number of spanwise station used (interpolated from input)
                  Use -1 or None to use number of stations from input file
        """
        damp_zeta     = None
        RayleighCoeff = None
        DampMat       = None
        s_start       = 0
        int_method    = 'Flex'
        if algo=='OpenFAST': 
            int_method='OpenFAST'

        # --- Reading properties, coefficients
        exp = np.arange(2,7)
        if 'BldProp' in inp.keys():
            # --- Blade
            name      = 'bld'
            shapeBase = ['BldFl1','BldFl2','BldEdg']
            if shapes is None:
                shapes=[0,1,2]
            coeff = np.zeros((len(exp), len(shapes)))
            for iishape, ishape in enumerate(shapes):
                base=shapeBase[ishape]
                coeff[0, iishape] = inp[base+'Sh(2)']
                coeff[1, iishape] = inp[base+'Sh(3)']
                coeff[2, iishape] = inp[base+'Sh(4)']
                coeff[3, iishape] = inp[base+'Sh(5)']
                coeff[4, iishape] = inp[base+'Sh(6)']
            damp_zeta = np.array([ inp['BldFlDmp(1)'], inp['BldFlDmp(2)'], inp['BldEdDmp(1)']])/100
            damp_zeta=damp_zeta[shapes]
            mass_fact = inp['AdjBlMs']   # Factor to adjust blade mass density (-)
            prop      = inp['BldProp']  
            s_bar, m, EIFlp, EIEdg  =prop[:,0], prop[:,3], prop[:,4], prop[:,5]
            """
            BldBodyStartAtRoot    BldBodyAndSpanStartAtR     BldBodyAndSpanStartAtR
                  /                          /                       /
                 /                          /                       /
                BS0                        /                       S!=0
                |                         /                       /
                R                        RBS0                    RB

            """

            if algo=='OpenFAST': 
                if (nSpan is None or nSpan==-1):
                    nSpan = ED['BldNodes']
                    #print('Using nSpan = BldNodes = ',nSpan)
                s_span  = s_bar*(ED['TipRad']-ED['HubRad'])
                if not bldStartAtRotorCenter:
                    raise Exception('For OpenFAST algo, bldStartAtRotorCenter should be true')
                s_start = ED['HubRad'] # Will be added later
                r_O = [0,0,0] # NOTE: blade defined wrt point R

            else:
#                 if not bldStartAtRotorCenter:
#                     # TODO TODO TODO, do we want the offset here, or when connecting bld to hub/rot?
#                     r_O = [0,0,ED['HubRad']] # NOTE: blade defined wrt point BldRoot
                # TODO we need two or three options with better naming
                if spanFrom0:
                    s_span=s_bar*(ED['TipRad']-ED['HubRad']) + ED['HubRad'] # NOTE: span starting at HubRad
                    if np.abs(s_span[0])<1e-6:
                        pass    
                    else:
                        # We add two positions with zero before
                        s_span = np.concatenate(([0,s_span[0]*0.99],s_span))
                        m      = np.concatenate(([0,0],m))
                        EIFlp  = np.concatenate(([0,0],EIFlp))
                        EIEdg  = np.concatenate(([0,0],EIEdg))
                    #s_span=s_bar*ED['TipRad'] # NOTE: this is a wrong scaling
                else:
                    s_span=s_bar*(ED['TipRad']-ED['HubRad']) + ED['HubRad'] # NOTE: span starting at HubRad
                r_O = [0,0,0] # NOTE: blade defined wrt point R for now
            #print(s_span)

            psi_B= 0
            if main_axis=='x':
                R_SB = R_z(0*np.pi + psi_B)
            elif main_axis=='z':
                R_SB = R_x(0*np.pi + psi_B)
            R_SB = np.dot(R_SB, R_y(ED['PreCone(1)']*np.pi/180))  # Blade 2 shaft
            R_b2g= R_SB

        elif 'TowProp' in inp.keys():
            # --- Tower
            name      = 'twr'
            shapeBase = ['TwFAM1','TwFAM2','TwSSM1','TwSSM2']
            if shapes is None:
                shapes=[0,1,2,3]
            coeff = np.zeros((len(exp), len(shapes)))
            for iishape, ishape in enumerate(shapes):
                base=shapeBase[ishape]
                coeff[0, iishape] = inp[base+'Sh(2)']
                coeff[1, iishape] = inp[base+'Sh(3)']
                coeff[2, iishape] = inp[base+'Sh(4)']
                coeff[3, iishape] = inp[base+'Sh(5)']
                coeff[4, iishape] = inp[base+'Sh(6)']
            damp_zeta = np.array([inp['TwrFADmp(1)'], inp['TwrFADmp(2)'], inp['TwrSSDmp(1)'], inp['TwrSSDmp(2)']])/100 # structural damping ratio 
            damp_zeta=damp_zeta[shapes]

            mass_fact = inp['AdjTwMa']                                              # Factor to adjust tower mass density (-)
            prop     = inp['TowProp']
            span_max = ED['TowerHt']-ED['TowerBsHt']
            s_bar, m, EIFlp, EIEdg  = prop[:,0], prop[:,1], prop[:,2], prop[:,3]
            r_O = [0,0,ED['TowerBsHt']]
            R_b2g=np.eye(3)
            s_span=s_bar*span_max

            if algo=='OpenFAST': 
                if (nSpan is None or nSpan==-1):
                    nSpan = ED['TwrNodes']
                    #print('Using nSpan = TwrNodes = ',nSpan)

        elif 'SttcSolve' in inp.keys():
            # --- Substructure / fnd
            from welib.fast.subdyn import SubDyn   
            name = 'fnd'
            sd = SubDyn(sdData = inp)
            p, damp_zeta, RayleighCoeff, DampMat = sd.toYAMSData(shapes)
            r_O   = p['r_O']
            R_b2g = p['R_b2g']

        else:
            print(inp.keys())
            raise Exception('Body type not supported, key `BldProp`, `TowProp`, or `SttcSolve` not found in file')

        try:
            gravity=ED['Gravity']
        except:
            if gravity is None:
                print('[WARN] yams/bodies.py: gravity is no longer present in elastodyn file, provide it as input')

        if name in ['twr','bld']:
            m *= mass_fact
            p = GeneralizedMCK_PolyBeam(s_span, m, EIFlp, EIEdg, coeff, exp, damp_zeta, jxxG=jxxG, 
                    gravity=gravity, Mtop=Mtop, Omega=Omega, nSpan=nSpan, bAxialCorr=bAxialCorr, bStiffening=bStiffening, main_axis=main_axis, shapes=shapes, algo=algo, s_start=s_start)
#             from welib.fast.elastodyn import bladeParameters
#             pGM=p
#             pED = bladeParameters(ED.filename)
#             import pdb; pdb.set_trace()


        elif name in ['fnd']:
            pass

        else:
            raise NotImplementedError()

        # TODO TODO sort out span for Blades and HubRad 

        BeamBody.__init__(self, name, p['s_span'], p['s_P0'], p['m'], p['EI'], p['PhiU'], p['PhiV'], p['PhiK'], jxxG=p['jxxG'], 
                s_min=p['s_min'], s_max=p['s_max'],
                r_O = r_O, R_b2g=R_b2g,  # NOTE: this is lost in YAMS
                damp_zeta=damp_zeta, RayleighCoeff=RayleighCoeff, DampMat=DampMat,
                bAxialCorr=bAxialCorr, bOrth=name=='bld', gravity=gravity, Mtop=Mtop, Omega=Omega, bStiffening=bStiffening, main_axis=main_axis,
                massExpected=massExpected,
                int_method=int_method
                )
