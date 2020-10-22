"""
Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

import numpy as np
from .flexibility import GMBeam, GKBeam, GKBeamStiffnening, polymode
from .utils import *

# --- To ease comparison with sympy version
from numpy import eye, cross, cos ,sin
def Matrix(m):
	return np.asarray(m)
def colvec(v): 
    v=np.asarray(v).ravel()
    return np.array([[v[0]],[v[1]],[v[2]]])


# --------------------------------------------------------------------------------}
# --- Connections 
# --------------------------------------------------------------------------------{
class Connection():
    def __init__(self,Type,RelPoint=None,RelOrientation=None,JointRotations=None, OrientAfter=True):
        if RelOrientation is None:
            RelOrientation=eye(3)
        if RelPoint is None:
            RelPoint=colvec([0,0,0])

        self.Type=Type
        
        self.s_C_0_inB = RelPoint
        self.s_C_inB   = self.s_C_0_inB
        self.R_ci_0    = RelOrientation
        self.R_ci      = self.R_ci_0     
        self.OrientAfter= OrientAfter

        if self.Type=='Rigid':
            self.nj=0
        elif self.Type=='SphericalJoint':
            self.JointRotations=JointRotations;
            self.nj=len(self.JointRotations);
        else:
            raise NotImplementedError()

    def updateKinematics(j,q):
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


# --------------------------------------------------------------------------------}
# --- Bodies 
# --------------------------------------------------------------------------------{
class Body(object):
    def __init__(B,Name=''):
        B.Children    = []
        B.Connections = []
        B.Name        = Name
        B.MM     = None
        B.B           = [] # Velocity transformation matrix
        B.updatePosOrientation(colvec([0,0,0]), eye(3))

    def updatePosOrientation(o,x_0,R_0b):
        o.r_O = x_0      # position of body origin in global coordinates
        o.R_0b=R_0b      # transformation matrix from body to global

    def connectTo(self, Child, Point=None, Type=None, RelOrientation=None, JointRotations=None, OrientAfter=True):
        if Type =='Rigid':
            c=Connection(Type, RelPoint=Point, RelOrientation = RelOrientation)
        else: # TODO first node, last node
            c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations, OrientAfter=OrientAfter)
        self.Children.append(Child)
        self.Connections.append(c)

    def setupDOFIndex(o,n):
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
            n=child.setupDOFIndex(n)
        return n

    def __repr__(B):
        pass

    @property
    def R_bc(self):
        return eye(3);
    @property
    def Bhat_x_bc(self):
        return Matrix(np.zeros((3,0)))
    @property
    def Bhat_t_bc(self):
        return Matrix(np.zeros((3,0)))

    def updateChildrenKinematicsNonRecursive(p,q):
        # At this stage all the kinematics of the body p are known
        # Useful variables
        R_0p =  p.R_0b
        B_p  =  p.B
        r_0p  = p.r_O

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
            #print('R_ci',p.Name)
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
            #print('R_pi',p.Name)
            #print(R_pi)
            #print('R_0p',p.Name)
            #print(R_0p)
            #print('R_0i',p.Name)
            #print(R_0i)

            # Position of connection point in P and 0 system
            r_pi_inP= conn_pi.s_C_inB
            r_pi    = np.dot (R_0p , r_pi_inP )
            #print('r_pi')
            #print(r_pi_inP)
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

    def getFullM(o,M):
        if not isinstance(o,GroundBody):
            MqB      = fBMB(o.BB_inB,o.MM)
            n        = MqB.shape[0]
            M[:n,:n] = M[:n,:n]+MqB     
        for c in o.Children:
            M=c.getFullM(M)
        return M
        
    def getFullK(o,K):
        if not isinstance(o,GroundBody):
            KqB      = fBMB(o.BB_inB,o.KK)
            n        = KqB.shape[0]
            K[:n,:n] = K[:n,:n]+KqB     
        for c in o.Children:
            K=c.getFullK(K)
        return K
        
    def getFullD(o,D):
        if not isinstance(o,GroundBody):
            DqB      = fBMB(o.BB_inB,o.DD)
            n        = DqB.shape[0]
            D[:n,:n] = D[:n,:n]+DqB     
        for c in o.Children:
            D=c.getFullD(D)
        return D

    @property
    def nf(B):
        if hasattr(B,'PhiU'):
            return len(B.PhiU)
        else:
            return 0

    @property
    def Mass(B):
        if B.MM is None:
            return 0
        return B.MM[0,0]



    def updateKinematics(o,x_0,R_0b,gz,v_0,a_v_0):
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

# --------------------------------------------------------------------------------}
# --- Ground Body 
# --------------------------------------------------------------------------------{
class GroundBody(Body):
    def __init__(B):
        super(GroundBody,B).__init__('Grd')

# --------------------------------------------------------------------------------}
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class RigidBody(Body):
    def __init__(B, Name, Mass, J_G, rho_G):
        """
        Creates a rigid body 
        """
        super(RigidBody,B).__init__()
        B.s_G_inB = rho_G
        B.J_G_inB = J_G
        B.J_O_inB = translateInertiaMatrixFromCOG(B.J_G_inB, Mass, B.s_G_inB)
        B.MM = fGMRigidBody(Mass,B.J_O_inB,B.s_G_inB)
        B.DD = np.zeros((6,6))
        B.KK = np.zeros((6,6))


def fGMRigidBody(Mass,J,rho):
    """ Generalized mass matrix for a rigid body (i.e. mass matrix) Eq.(15) of [1] """
    S=Mass*skew(rho)
    MM=np.zeros((6,6))
    MM[0:3,0:3] = Mass*np.eye(3);
    MM[0:3,3:6] = -S;
    MM[3:6,0:3] = S ; # transpose(S)=-S;
    MM[3:6,3:6] = J ;
    return MM


# --------------------------------------------------------------------------------}
# --- Beam Body 
# --------------------------------------------------------------------------------{
class BeamBody(Body):
    def __init__(B, s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=None, s_G0=None, bAxialCorr=False, bOrth=False, Mtop=0, bStiffening=True, gravity=None,main_axis='z'):
        """ 
          Points P0 - Undeformed mean line of the body
        """
        super(BeamBody,B).__init__()
        B.main_axis = main_axis

        B.s_span = s_span
        B.m      = m
        B.s_G0   = s_G0
        B.PhiU   = PhiU
        B.PhiV   = PhiV
        B.PhiK   = PhiK
        B.jxxG   = jxxG
        B.s_P0   = s_P0
        B.EI     = EI
        if jxxG is None:
            B.jxxG   = 0*m
        if B.s_G0 is None:
            B.s_G0=B.s_P0
    
        B.s_G    = B.s_G0
        B.bAxialCorr = bAxialCorr
        B.bOrth      = bOrth
        B.bStiffening= bStiffening
        B.Mtop       = Mtop
        B.gravity    = gravity

        B.computeMassMatrix()
        B.computeStiffnessMatrix()
        B.DD = np.zeros((6+B.nf,6+B.nf))

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

    @property
    def Bhat_x_bc(self,iNode=-1):
        Bhat_x_bc = Matrix(np.zeros((3,self.nf)))
        for j in np.arange(self.nf):
            Bhat_x_bc[:,j]=self.PhiU[j][:,iNode] #  along x
        return Bhat_x_bc

    @property
    def Bhat_t_bc(self,iNode=-1):
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




    def computeStiffnessMatrix(B):
        B.KK0 = GKBeam(B.s_span, B.EI, B.PhiK, bOrth=B.bOrth)
        if B.bStiffening:
            B.KKg = GKBeamStiffnening(B.s_span, B.PhiV, B.gravity, B.m, B.Mtop, main_axis=B.main_axis)
        else:
            B.KKg=B.KK0*0

        B.KK=B.KK0+B.KKg

    def computeMassMatrix(B):
        B.MM = GMBeam(B.s_G, B.s_span, B.m, B.PhiU, jxxG=B.jxxG, bUseIW=True, main_axis=B.main_axis, bAxialCorr=B.bAxialCorr, bOrth=B.bOrth)

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

            # Position of mean line
            o.s_P=o.s_P0+o.U;

            # Position of deflected COG
            # TODO TODO TODO mean_axis not x
            o.rho_G      = np.zeros((3,o.nSpan))
            if o.main_axis=='x':
                o.rho_G[1,:] = o.rho_G0_inS[1,:]*np.cos(o.V_tot[0,:])-o.rho_G0_inS[2,:]*np.sin(o.V_tot[0,:]);
                o.rho_G[2,:] = o.rho_G0_inS[1,:]*np.sin(o.V_tot[0,:])+o.rho_G0_inS[2,:]*np.cos(o.V_tot[0,:]);
            else:
                raise NotImplementedError()
                o.rho_G[1,:] = o.rho_G0_inS[1,:]*np.cos(o.V_tot[0,:])-o.rho_G0_inS[2,:]*np.sin(o.V_tot[0,:]);
                o.rho_G[2,:] = o.rho_G0_inS[1,:]*np.sin(o.V_tot[0,:])+o.rho_G0_inS[2,:]*np.cos(o.V_tot[0,:]);
            o.s_G = o.s_P+o.rho_G;
            # Alternative:
            #rho_G2     = zeros(3,o.nSpan);
            #rho_G2(2,:) = o.rho_G0(2,:).*cos(o.V(1,:))-o.rho_G0(3,:).*sin(o.V(1,:));
            #rho_G2(3,:) = o.rho_G0(2,:).*sin(o.V(1,:))+o.rho_G0(3,:).*cos(o.V(1,:));
            #compare(o.rho_G,rho_G2,'rho_G');
            # Position of connection point
            print('TODO connection points')
            #for ic=1:length(o.Connections)
            #    iNode=o.Connections{ic}.ParentNode;
            #    %o.Connections{ic}.s_C_inB = o.U(1:3,iNode);
            #    o.Connections{ic}.s_C_inB = o.s_P(1:3,iNode);

    @property
    def nSpan(B):
        return len(B.s_span)


# --------------------------------------------------------------------------------}
# --- Uniform Beam Body 
# --------------------------------------------------------------------------------{
class UniformBeamBody(BeamBody):
    def __init__(B, Name, nShapes, nSpan, L, EI0, m, Mtop=0, jxxG=None, GKt=None, bAxialCorr=True, bCompatibility=False, bStiffnessFromGM=False, bStiffening=True, gravity=None, main_axis='x'):

        import welib.beams.theory as bt
        if jxxG is None:
            jxxG=0
        if GKt is None:
            GKt=0

        A=1; rho=A*m;
        x=np.linspace(0,L,nSpan);
        # Mode shapes
        freq,s_span,U,V,K = bt.UniformBeamBendingModes('unloaded-topmass-clamped-free',EI0,rho,A,L,x=x,Mtop=Mtop)
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
        super(UniformBeamBody,B).__init__(s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=jxxG, bAxialCorr=bAxialCorr, Mtop=Mtop, bStiffening=bStiffening, gravity=gravity, main_axis=main_axis)


# --------------------------------------------------------------------------------}
# --- FAST Beam body 
# --------------------------------------------------------------------------------{
class FASTBeamBody(BeamBody):
    def __init__(B,body_type,ED,inp,Mtop,nShapes=2,main_axis='x',nSpan=40,bAxialCorr=False,bStiffening=True):
        """ 
        INPUTS:
           nSpan: number of spanwise station used (interpolated from input)
                  Use -1 or None to use number of stations from input file
        """
        # --- Reading coefficients
        exp = np.arange(2,7)
        if body_type.lower()=='blade':
            coeff = np.array([[ inp['BldFl1Sh(2)'], inp['BldFl2Sh(2)'], inp['BldEdgSh(2)']],
                              [ inp['BldFl1Sh(3)'], inp['BldFl2Sh(3)'], inp['BldEdgSh(3)']],
                              [ inp['BldFl1Sh(4)'], inp['BldFl2Sh(4)'], inp['BldEdgSh(4)']],
                              [ inp['BldFl1Sh(5)'], inp['BldFl2Sh(5)'], inp['BldEdgSh(5)']],
                              [ inp['BldFl1Sh(6)'], inp['BldFl2Sh(6)'], inp['BldEdgSh(6)']]])

            damp_zeta  = np.array([ inp['BldFlDmp(1)'], inp['BldFlDmp(2)'], inp['BldEdDmp(1)']])/100
            mass_fact = inp['AdjBlMs']                                              # Factor to adjust blade mass density (-)
          
            prop     = inp['BldProp']
            span_max = ED['TipRad']   # TODO TODO do somthing about hub rad
            s_bar, m, EIFlp, EIEdg  =prop[:,0], prop[:,3], prop[:,4], prop[:,5]

        elif body_type.lower()=='tower':
            coeff = np.array([[ inp['TwFAM1Sh(2)'], inp['TwFAM2Sh(2)'], inp['TwSSM1Sh(2)'], inp['TwSSM2Sh(2)']],
                              [ inp['TwFAM1Sh(3)'], inp['TwFAM2Sh(3)'], inp['TwSSM1Sh(3)'], inp['TwSSM2Sh(3)']],
                              [ inp['TwFAM1Sh(4)'], inp['TwFAM2Sh(4)'], inp['TwSSM1Sh(4)'], inp['TwSSM2Sh(4)']],
                              [ inp['TwFAM1Sh(5)'], inp['TwFAM2Sh(5)'], inp['TwSSM1Sh(5)'], inp['TwSSM2Sh(5)']],
                              [ inp['TwFAM1Sh(6)'], inp['TwFAM2Sh(6)'], inp['TwSSM1Sh(6)'], inp['TwSSM2Sh(6)']]])
            damp_zeta = np.array([inp['TwrFADmp(1)'], inp['TwrFADmp(2)'], inp['TwrSSDmp(1)'], inp['TwrSSDmp(2)']])/100 # structural damping ratio 
            mass_fact = inp['AdjTwMa']                                              # Factor to adjust tower mass density (-)

            prop     = inp['TowProp']
            span_max = ED['TowerHt']-ED['TowerBsHt']
            s_bar, m, EIFlp, EIEdg  = prop[:,0], prop[:,1], prop[:,2], prop[:,3]


        else:
            raise Exception('Body type not supported {}'.format(body_type))
        nShpMax=coeff.shape[1]
        if nShapes>nShpMax:
            raise Exception('A maximum of {} shapes function possible with FAST {} body'.format(nShpMax,body_type))

        gravity=ED['Gravity']

        # --- Interpolating structural properties
        m *= mass_fact
        if nSpan is None or nSpan<0:
            nSpan  = len(prop[:,0])
        s_span     = np.linspace(0,span_max,nSpan)
        m          = np.interp(s_span/span_max,s_bar,m)
        EIFlp      = np.interp(s_span/span_max,s_bar,EIFlp)
        EIEdg      = np.interp(s_span/span_max,s_bar,EIEdg)

        # --- Definition of main directions
        ShapeDir=np.zeros(nShpMax).astype(int)
        DirNames=['x','y','z']
        EI =np.zeros((3,nSpan))
        if main_axis=='x':
            iMain = 0  # longitudinal axis along x
            ShapeDir[0:2] = 2 # First two shapes are along z (flapwise/Fore-Aft)
            ShapeDir[2:]  = 1 # Third shape along along y (edgewise/Side-Side) Sign...
            EI[2,:] = EIFlp
            EI[1,:] = EIEdg
        elif main_axis=='z':
            iMain = 2  # longitudinal axis along z
            ShapeDir[0:2] = 0 # First two shapes are along x (flapwise/Fore-Aft)
            ShapeDir[2:]  = 1 # Third shape along along y (edgewise/Side-Side) Sign...
            EI[0,:] = EIFlp
            EI[1,:] = EIEdg
        else:
            raise NotImplementedError()

        # --- Undeflected shape
        s_P0          = np.zeros((3,nSpan))
        s_P0[iMain,:] = s_span 
        # TODO blade COG
        jxxG=m*0

        # --- Shape functions
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        PhiV = np.zeros((nShapes,3,nSpan)) # Slope
        PhiK = np.zeros((nShapes,3,nSpan)) # Curvature
        for j in np.arange(nShapes):
            iAxis = ShapeDir[j]
            PhiU[j][iAxis,:], PhiV[j][iAxis,:], PhiK[j][iAxis,:] = polymode(s_span,coeff[:,j],exp)


        super(FASTBeamBody,B).__init__(s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=jxxG, bAxialCorr=bAxialCorr, bOrth=body_type=='blade', gravity=gravity,Mtop=Mtop, bStiffening=bStiffening, main_axis=main_axis)

        # Damping
        B.DD=np.zeros((6+nShapes,6+nShapes))

        # Using diagonal damping
        for j in range(nShapes):
            m             = B.MM[6+j,6+j]
            k             = B.KK[6+j,6+j]
            om            = np.sqrt(k/m)
            xi            = damp_zeta[j]*2*np.pi
            c             = xi * m * om / np.pi
            B.DD[6+j,6+j] = c

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
